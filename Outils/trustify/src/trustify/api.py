"""Public programmatic API for trustify.

All CLI commands dispatch through this module. Callers writing Python code
should import from here:

    from trustify import generate_schema, check, batch_check, load_dataset, init_config

(The top-level `trustify/__init__.py` re-exports these names.)
"""

import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from trustify.cache import (
    atomic_write_dir,
    cache_lock,
    compute_content_hash,
    default_cache_root,
    is_cached,
    pyd_module_filename,
)
from trustify.core.misc_utilities import import_parser_module
from trustify.core.trad2_pydantic import generate_pyd_and_pars
from trustify.core.trad2_utilities import TRAD2Content
from trustify.projects import Project, effective_projects, effective_trust_root, resolve_projects, validate_projects

# Sentinel `__version__` value emitted by trustify/__init__.py when
# `importlib.metadata.version("trustify")` raises (broken dist info,
# partial reinstall). Detected at cache-write time and used to refuse
# the shared cache — see `generate_schema`. Audit 6.11.
_UNKNOWN_VERSION = "0.0.0+unknown"


def _version() -> str:
    # Single source of truth lives in trustify/__init__.py — itself derived
    # from `importlib.metadata.version("trustify")`. Late-import here because
    # api.py is loaded by trustify/__init__.py itself and the attribute is
    # set above the import block so it's already in place by the time
    # this function actually runs.
    from trustify import __version__

    return __version__


def generate_schema(
    projects: list[str] | None = None, trust_root: str | None = None, out: str | Path | None = None
) -> Path:
    """Generate (or reuse-cached) the TRUST schema for the given project set.

    Returns the directory containing `trustify_gen.py`,
    `trustify_gen_pyd_<digest>.py`, `TRAD2_trustify`,
    `source_locations.json`, and `provenance.json`.

    - `projects` and `trust_root` are normalized via `projects.resolve_projects`.
      When `trust_root` is None it falls back to `$TRUST_ROOT`, matching
      `load_dataset` / `batch_check` / `generate_keywords`.
    - `out` overrides the default `~/.cache/trustify/<hash>/` location.
    - If the target directory already contains the generated files, nothing
      is regenerated.
    - `provenance.json` records the at-generation absolute paths +
      timestamp + trustify version. It is display-only metadata —
      consumed by the `trustify cache` subcommand but excluded from
      the content hash so two clones of the same project set at
      different absolute paths share one cache entry.
    """
    import datetime
    import json
    import time

    from trustify.core.source_location import build_projects_dict

    # Refuse the shared cache when the trustify version cannot be
    # determined (audit 6.11). The sentinel `__version__` would
    # otherwise become part of the cache hash, collapsing every
    # broken install into one shared bucket and risking cross-tree
    # schema reuse. `out=<dir>` bypasses the cache root entirely, so
    # the user-controlled out path is unaffected.
    if out is None and _version() == _UNKNOWN_VERSION:
        raise RuntimeError(
            "trustify cannot determine its installed version "
            f"(__version__={_UNKNOWN_VERSION!r}) — refusing to use the shared "
            "schema cache so two broken installs cannot serve each other's "
            "generated schema. Fix the install (re-run `pip install` to "
            "restore dist info), or pass --out=<dir> to write the schema "
            "to a user-controlled directory that bypasses the cache."
        )

    step_seconds: dict = {}

    def _step(name):
        class _Ctx:
            def __enter__(self_):
                self_._t0 = time.perf_counter()
                return self_

            def __exit__(self_, *exc):
                step_seconds[name] = round(time.perf_counter() - self_._t0, 3)

        return _Ctx()

    projects = effective_projects(projects)
    trust_root = effective_trust_root(trust_root)
    resolved: list[Project] = resolve_projects(projects, trust_root)
    # Every non-trust_root entry must be a baltik (project.cfg present).
    # Raises ConfigError with an actionable message otherwise.
    validate_projects(resolved)
    # `resolve_projects` prepends trust_root to the project list when
    # the caller passed one, so `resolved[0]` is the trust_root iff
    # `trust_root` is truthy.
    if trust_root:
        proj_dict = build_projects_dict(
            trust_root=str(resolved[0].path),
            projects=[str(p.path) for p in resolved[1:]],
        )
    else:
        proj_dict = build_projects_dict(
            trust_root=None,
            projects=[str(p.path) for p in resolved],
        )
    # Every project (trust_root included) contributes its `src/` subtree.
    # `scanSourceFiles` walks `src_dirs` in order and lets later entries
    # override earlier ones on basename collision — keeping the BALTIK
    # overlay precedence the previous `build/src` layout used to provide.
    src_dirs = [str(p.path / "src") for p in resolved]
    tr_for_scan = str(resolved[0].path) if trust_root else None
    with _step("scan_sources"):
        content = TRAD2Content.BuildFromOrgAndSources(
            None,
            src_dirs,
            projects=proj_dict,
            trust_root=tr_for_scan,
        )

    digest = compute_content_hash(content, _version())
    final_dir = Path(out) if out is not None else default_cache_root() / digest
    if is_cached(final_dir):
        # Cache hit: original timings remain in the existing
        # provenance.json from the first generation.
        return final_dir

    def _write_and_trim():
        # Write everything into a sibling temp dir; atomic_write_dir renames
        # it onto final_dir only after provenance.json (the completion
        # sentinel) is on disk. A SIGKILL inside the with-block leaves the
        # temp dir cleaned up and final_dir untouched (audit 2.1).
        with atomic_write_dir(final_dir) as out_dir:
            trad2_path = out_dir / "TRAD2_trustify"
            # The pyd module name carries the content digest so distinct
            # schemas loaded in the same process land in distinct
            # `sys.modules` keys. Using the content hash (not an abs-path
            # hash) keeps the generated parser file byte-identical across
            # clones — its `from ... import *` line references this same
            # digest-suffixed module name.
            pyd_path = out_dir / pyd_module_filename(digest)
            pars_path = out_dir / "trustify_gen.py"
            # `TRAD2Content.toTRAD2(f)` writes `f` (TRAD2 text) and a
            # sibling `source_locations.json` in the same directory.
            # Passing `None` as the locations arg tells
            # `generate_pyd_and_pars` to auto-detect that sibling JSON.
            with _step("write_trad2"):
                content.toTRAD2(str(trad2_path))
            with _step("generate_pyd_and_pars"):
                generate_pyd_and_pars(trad2_path, None, pyd_path, pars_path)

            total = round(sum(step_seconds.values()), 3)
            provenance = {
                "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
                "trustify_version": _version(),
                "projects": proj_dict,
                "generation_time_seconds": total,
                "generation_steps_seconds": step_seconds,
            }
            # provenance.json MUST be written last — `is_cached` uses
            # its presence as the completion sentinel.
            (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

        # Auto-trim: in the default-cache-root case, delete any sibling
        # entries that came from the exact same project set (typical
        # when source line numbers shift during active TRUST / baltik
        # development). Skipped when `out=` overrides the cache root
        # (user-controlled territory) or when `.trustify.json`
        # `cache.auto_trim: false` / `TRUSTIFY_NO_AUTO_TRIM` opts out.
        if out is None:
            from trustify.cache import should_auto_trim, trim_same_project_set

            if should_auto_trim():
                evicted = trim_same_project_set(final_dir)
                if evicted:
                    print(
                        f"trustify cache auto-trim: evicted {len(evicted)} "
                        f"superseded entr{'y' if len(evicted) == 1 else 'ies'} "
                        f"from the same project set ({', '.join(e[:8] for e in evicted)})",
                        file=sys.stderr,
                    )

    if out is None:
        # Cache mode: serialise writers + the trim that follows under
        # an exclusive cache lock so a concurrent generate_schema
        # cannot rmtree a sibling we're about to use (audit 2.2). The
        # re-check inside the lock catches the common case where
        # another writer completed our entry while we were waiting.
        with cache_lock("exclusive"):
            if is_cached(final_dir):
                return final_dir
            _write_and_trim()
    else:
        # Override mode: user-controlled destination, no cross-process
        # coordination needed.
        _write_and_trim()
    return final_dir


def _resolve_schema_module(schema, projects, trust_root):
    """Pick the schema module to parse against: either the pre-built one
    pointed at by `schema=<dir>`, or the result of `generate_schema(...)`.
    Centralised so `load_dataset` and `check` apply the exact same rule.

    The actual import is wrapped in a SHARED cache lock so a concurrent
    generate_schema (which takes EXCLUSIVE for write+trim) cannot
    rmtree the directory we're importing from. A FileNotFoundError
    inside the SHARED block means a writer trimmed the dir in the
    narrow window between generate_schema's release and our SHARED
    acquire — we re-call generate_schema once to rebuild and retry.
    """
    if schema is not None:
        with cache_lock("shared"):
            return import_parser_module(Path(schema) / "trustify_gen.py")
    # $TRUST_ROOT fallback is centralised in generate_schema → projects.effective_trust_root.
    last_exc: FileNotFoundError | None = None
    for _attempt in range(2):
        gen_dir = generate_schema(projects=projects, trust_root=trust_root)
        try:
            with cache_lock("shared"):
                return import_parser_module(gen_dir / "trustify_gen.py")
        except FileNotFoundError as exc:
            last_exc = exc
    raise RuntimeError(
        "trustify: persistent race with concurrent cache trim — "
        "two regenerate-then-import attempts both lost the dir before "
        "the SHARED lock acquired. Re-run; if it persists, set "
        "TRUSTIFY_NO_AUTO_TRIM to disable auto-trim."
    ) from last_exc


def _parse_text(text: str, schema_mod, filename: str = "") -> tuple[Any, Any]:
    """Tokenize `text` and parse it into a dataset pydantic model.

    Returns `(dataset, parser)` so the caller can reuse the parser
    (e.g. to skip a second tokenization in `prune_after_end`).
    """
    from trustify.core import trust_parser

    parser = trust_parser.TRUSTParser()
    parser.tokenize(text)
    stream = trust_parser.TRUSTStream(parser=parser, file_nam=filename)
    cls = schema_mod.get_parser_class("Dataset")
    return cls.ReadFromTokens(stream), parser


def load_dataset(
    filename: str | Path,
    projects: list[str] | None = None,
    trust_root: str | None = None,
    schema: str | Path | None = None,
) -> Any:
    """Load a TRUST dataset and return it as a pydantic model.

    Two input modes (resolved in order):
      1. `schema=` — directory produced by a previous `generate()`. Reused as-is.
      2. `projects=` / `trust_root=` — generate (or reuse-cached) the schema, then load.

    If neither is supplied, falls back to `$TRUST_ROOT` from the environment.
    """
    schema_mod = _resolve_schema_module(schema, projects, trust_root)
    with open(filename, encoding="utf-8") as f:
        text = f.read().rstrip()
    ds, _parser = _parse_text(text, schema_mod, filename=str(filename))
    return ds


@dataclass
class CheckResult:
    """Single-dataset check outcome."""

    path: str
    status: str  # "PASSED" | "FAILED" | "SKIPPED"
    message: str = ""  # error or skip reason; empty when PASSED
    obsolete_marker: bool = False
    """True when the dataset carries a `# TRUSTIFY NOT #` opt-out marker
    AND `no_skip=True` was passed AND the parse PASSED. Signals to the
    caller that the marker has become stale and can be removed."""
    marker_justification: str = ""
    """Free-form text from a `# TRUSTIFY NOT: <text> #` opt-out tag.
    Empty when no tag is present, when the bare `# TRUSTIFY NOT #` form
    is used, or when the inline text is blank after whitespace
    collapsing. Populated on every return path that saw the tag —
    SKIPPED, PASSED-with-obsolete-marker, AND FAILED-with-no-skip — so
    every renderer can read it uniformly."""
    diff: str = ""
    """Unified diff between the serialized round-trip output and the
    expected (pruned) source, populated only when the round-trip
    equality check fails. Empty for every other FAILED reason (parse
    error, malformed marker, OSError) and for PASSED / SKIPPED. Captured
    so the batch-check summary (and `--summary-out` file) can surface
    the actual mismatch without users having to re-run on each file."""


@dataclass
class BatchCheckResult:
    """Aggregated outcome over a batch of datasets."""

    total: int
    passed: int
    failed: int
    skipped: int
    results: list[CheckResult]


_MARKER_RE = None  # compiled lazily — well-formed TRUSTIFY NOT tag
_MARKER_PREFIX_RE = None  # compiled lazily — bare sentinel; matches even malformed tags


def _iter_hash_comments(text: str):
    """Yield `(start, substring)` for each standalone `#...#` comment.

    `#` characters that appear inside a `/* ... */` block or a
    `"..."` quoted string are NOT comment openers — they're skipped.
    An unterminated `#` (no closing `#` before EOF) is dropped
    silently. Used by `_extract_marker` to apply the TRUSTIFY-NOT
    regex only inside genuine hash-comment contexts (audit 1.7).
    """
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]
        if ch == "#":
            start = i
            i += 1
            while i < n and text[i] != "#":
                i += 1
            if i >= n:
                return  # unterminated comment
            yield (start, text[start : i + 1])
            i += 1
        elif ch == "/" and i + 1 < n and text[i + 1] == "*":
            i += 2
            while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                i += 1
            i += 2  # past `*/`
        elif ch == '"':
            i += 1
            while i < n and text[i] != '"':
                i += 1
            i += 1  # past closing `"`
        else:
            i += 1


def _extract_marker(text: str) -> tuple[bool, str, str]:
    """Return `(has_marker, justification, error)`.

    - `has_marker` is True if a well-formed `# TRUSTIFY NOT[: ...] #`
      comment is present anywhere in `text`.
    - `justification` is the captured text after the optional `:`,
      whitespace-collapsed (`' '.join(captured.split())`) — empty when
      no `:`-suffix, when the suffix is blank, or when no marker is
      present.
    - `error` is a non-empty message when a TRUSTIFY NOT sentinel is
      found but the surrounding tag is malformed (e.g. `# TRUSTIFY NOT
      toto #` with no colon). In that case `has_marker` is False and
      the caller MUST fail the check — silently ignoring the malformed
      form would be a surprising correctness hole.

    The regex is applied ONLY to standalone `#...#` comments;
    `# TRUSTIFY NOT #` embedded in a `/* ... */` block or in a
    `"..."` quoted string does NOT count (audit 1.7).
    """
    global _MARKER_RE, _MARKER_PREFIX_RE
    if _MARKER_RE is None:
        import re

        _MARKER_RE = re.compile(r"#\s*TRUSTIFY\s+NOT(?:\s*:\s*([^#]*))?\s*#")
        _MARKER_PREFIX_RE = re.compile(r"#\s*TRUSTIFY\s+NOT\b")
    first_match: re.Match | None = None
    for cmt_start, cmt_text in _iter_hash_comments(text):
        # `_MARKER_RE` includes the `#` markers, so it can be applied
        # directly to the full comment substring.
        m = _MARKER_RE.search(cmt_text)
        prefix_m = _MARKER_PREFIX_RE.search(cmt_text)
        if prefix_m and m is None:
            # Sentinel-without-well-formed-tag inside a real comment:
            # malformed marker.
            line_no = text.count("\n", 0, cmt_start + prefix_m.start()) + 1
            return (
                False,
                "",
                f"ill-formed TRUSTIFY NOT tag at line {line_no}: "
                f"expected '# TRUSTIFY NOT #' or '# TRUSTIFY NOT: <text> #'",
            )
        if m is not None and first_match is None:
            first_match = m
    if first_match is None:
        return (False, "", "")
    captured = first_match.group(1) or ""
    return (True, " ".join(captured.split()), "")


def check(
    data_file: str | Path,
    projects: list[str] | None = None,
    trust_root: str | None = None,
    schema: str | Path | None = None,
    no_skip: bool = False,
) -> CheckResult:
    """Parse `data_file` against the generated/cached schema and verify a
    round-trip write produces the same dataset. Returns a `CheckResult` —
    `FAILED` results carry the exception type+message rather than raising,
    so `batch_check` can keep iterating after a single failure.

    `no_skip=True` disables the `TRUSTIFY NOT` opt-out marker; the
    dataset is checked regardless. Useful when verifying a fix to a
    previously-skipped file."""
    from trustify.core import misc_utilities

    path = str(data_file)
    # Read + marker extract sit OUTSIDE the parse try/except so the
    # FAILED return path under --no-skip can also carry the
    # justification (consumers render it on the obsolete WARNING line).
    try:
        with open(path, encoding="utf-8") as f:
            text = f.read()
    except OSError as e:
        return CheckResult(path=path, status="FAILED", message=f"{type(e).__name__}: {e}")
    has_marker, justification, marker_error = _extract_marker(text)
    if marker_error:
        # Malformed tag — fail fast with a clear message. The user
        # very likely intended a TRUSTIFY NOT opt-out and would not
        # want it silently treated as if no marker were present.
        return CheckResult(path=path, status="FAILED", message=marker_error)
    if not no_skip and has_marker:
        return CheckResult(
            path=path,
            status="SKIPPED",
            message="dataset opted out via 'TRUSTIFY NOT' marker",
            marker_justification=justification,
        )
    try:
        # Parse once and reuse the resulting parser for prune_after_end below.
        # Going through load_dataset would re-read the file *and* prune_after_end
        # would tokenize the same content a second time (~2 ms x N files on a
        # batch_check run).
        schema_mod = _resolve_schema_module(schema, projects, trust_root)
        ds, parser = _parse_text(text.rstrip(), schema_mod, filename=path)
        written = "".join(ds.toDatasetTokens())
        # Compare against the original with everything after `end`/`fin` pruned —
        # the parser stops at the end-of-dataset keyword, so trailing whitespace
        # or comments after it must not cause spurious round-trip mismatches.
        expected = misc_utilities.prune_after_end(text, parser=parser)
        # `print_on_diff=False` suppresses the inline `logger.error` so we
        # don't double-emit: the diff is now captured into CheckResult.diff
        # and rendered explicitly by the CLI (single-file `check` output
        # and `batch-check` summary block / `--summary-out` file).
        eq = misc_utilities.check_str_equality(written, expected, print_on_diff=False)
        if not eq.ok:
            return CheckResult(
                path=path,
                status="FAILED",
                message="dataset is not bit-identical after round-trip",
                marker_justification=justification,
                diff=eq.why,
            )
        # Marker was present, we bypassed it with no_skip, and the parse
        # passed → the marker is stale and can be removed.
        return CheckResult(
            path=path,
            status="PASSED",
            obsolete_marker=(no_skip and has_marker),
            marker_justification=justification,
        )
    except Exception as e:
        return CheckResult(
            path=path,
            status="FAILED",
            message=f"{type(e).__name__}: {e}",
            marker_justification=justification,
        )


def _check_worker(args: tuple) -> CheckResult:
    """Pool worker entry point. Module-level so it's picklable across forks."""
    file, schema_dir, no_skip = args
    return check(file, schema=schema_dir, no_skip=no_skip)


def _pick_mp_start_method() -> str:
    """Choose a `multiprocessing` start method that's safe given that
    trustify's parent process has already imported pydantic and the
    schema parser module (both pull in C extensions). Audit 2.7.

    Order:
      1. `$TRUSTIFY_MP_START_METHOD` (any non-empty value, NO_COLOR
         convention) — escape hatch for users who need to force
         `fork` for perf, or `spawn` to debug forkserver issues.
      2. `forkserver` if the platform advertises it (Linux + macOS):
         workers fork from a clean child process, so the parent's
         threads and C-extension internal locks never get into the
         children. Slightly slower than `fork` (each worker re-imports
         the schema) but no deadlock risk.
      3. `spawn` as the portable fallback (Windows).

    `fork` is deliberately NEVER picked silently: even on Linux it
    can deadlock when a pydantic-held C-extension lock is in an
    inopportune state at fork time, and Python 3.14 moves macOS off
    the fork default for the same reason.
    """
    import multiprocessing

    env = os.environ.get("TRUSTIFY_MP_START_METHOD")
    if env:
        return env
    if "forkserver" in multiprocessing.get_all_start_methods():
        return "forkserver"
    return "spawn"


def batch_check(
    data_files: list[str | Path],
    projects: list[str] | None = None,
    trust_root: str | None = None,
    schema: str | Path | None = None,
    no_skip: bool = False,
    verbose: bool = False,
    jobs: int = 1,
) -> BatchCheckResult:
    """Run `check` on every entry. The schema is generated (or reused from
    cache) ONCE before iteration; every per-file `check` reuses that schema
    via `schema=`, so we don't pay the source-scan cost N times.

    `no_skip=True` forwards to each per-file `check`, disabling the
    `TRUSTIFY NOT` opt-out marker.

    `verbose=True` prints one `[i/N] STATUS path[: message]` progress line to
    stdout per finished dataset (flushed). Off by default so library callers
    stay silent; the CLI flips it on (with `--quiet` to opt out).

    `jobs` (>=0) controls process-level parallelism. `jobs=1` (the
    default) runs every check inline in the current process. `jobs>=2`
    farms the per-file checks out to a `multiprocessing.Pool` of that
    many workers, using the start method chosen by
    `_pick_mp_start_method` (forkserver where available, else spawn —
    never plain `fork` silently, which can deadlock on pydantic's
    C-extension locks; see audit 2.7). Under forkserver/spawn the
    workers do NOT inherit the parent via fork+CoW, so each one
    re-imports the schema on first use — but the parent primes it
    first, so that re-import is a cache hit (a single filesystem stat
    + module load per worker, under the SHARED cache lock), not a
    fresh source scan. `jobs=0` is the conventional "use all available
    CPUs" shorthand: it resolves to `os.cpu_count()` (or 1 when that
    returns None). Results stream back in submission order so the
    `[i/N]` progress lines stay correctly numbered. Ctrl-C aborts the
    pool cleanly.
    """
    if jobs < 0:
        raise ValueError(f"batch_check: jobs must be >= 0, got {jobs}")
    if jobs == 0:
        import os

        jobs = os.cpu_count() or 1
    if schema is None:
        # $TRUST_ROOT fallback is centralised in generate_schema.
        schema = generate_schema(projects=projects, trust_root=trust_root)
    total = len(data_files)

    def _emit(i: int, r: CheckResult) -> None:
        if verbose:
            from trustify._color import colorize

            suffix = f": {r.message}" if r.message else ""
            if r.marker_justification:
                suffix += f" — {r.marker_justification}"
            print(colorize(f"[{i}/{total}] {r.status} {r.path}{suffix}"), flush=True)

    results: list[CheckResult] = []
    if jobs == 1:
        for i, f in enumerate(data_files, 1):
            r = check(f, schema=schema, no_skip=no_skip)
            results.append(r)
            _emit(i, r)
    else:
        import multiprocessing

        # Prime the schema in the parent so cache locks + import are
        # warm — workers will re-import the schema on first use under
        # `forkserver` / `spawn` (no fork+CoW inheritance), but the
        # cache hit + SHARED-lock path makes the cost a single
        # file-system stat + module load per worker.
        _resolve_schema_module(schema, None, None)
        work = [(str(f), str(schema), no_skip) for f in data_files]
        # Amortise IPC overhead across many tiny tasks while keeping
        # head-of-line latency small enough that verbose progress stays
        # responsive.
        chunksize = max(1, min(16, total // (jobs * 4) or 1))
        # `forkserver` (POSIX) over plain `fork` to avoid inheriting
        # C-extension internal locks from the parent — pydantic +
        # parser-module imports can deadlock under fork-then-imap.
        # See `_pick_mp_start_method` and audit 2.7.
        ctx = multiprocessing.get_context(_pick_mp_start_method())
        pool = ctx.Pool(jobs)
        try:
            for i, r in enumerate(pool.imap(_check_worker, work, chunksize=chunksize), 1):
                results.append(r)
                _emit(i, r)
        except KeyboardInterrupt:
            pool.terminate()
            pool.join()
            raise
        else:
            pool.close()
            pool.join()

    passed = sum(1 for r in results if r.status == "PASSED")
    skipped = sum(1 for r in results if r.status == "SKIPPED")
    failed = sum(1 for r in results if r.status == "FAILED")
    return BatchCheckResult(total=len(results), passed=passed, failed=failed, skipped=skipped, results=results)


def init_config(
    directory: str | Path = ".",
    projects: list[str] | None = None,
    trust_root: str | None = None,
    overwrite: bool = False,
) -> Path:
    """Write a `.trustify.json` template into `directory`. Raises
    `FileExistsError` if the file already exists and `overwrite` is False.
    Returns the path to the written file.

    Emits the full scaffold so users can discover every available knob:
    `workspace.projects` (list, defaults to []), `workspace.trust_root`
    (string, defaults to ""), `workspace.auto_resolve_dependencies`
    (bool, defaults to true), and `lsp.enum_dedup_threshold` (int,
    defaults to 4). Provided values override the empty placeholders.
    """
    out = Path(directory) / ".trustify.json"
    if out.exists() and not overwrite:
        raise FileExistsError(f"{out} already exists. Pass overwrite=True to replace it.")
    payload = {
        "workspace": {
            "projects": list(projects) if projects else [],
            "trust_root": trust_root or "",
            "auto_resolve_dependencies": True,
        },
        "lsp": {"enum_dedup_threshold": 4},
    }
    out.write_text(json.dumps(payload, indent=2) + "\n")
    return out


# Enum-attribute payload pattern: `chaine(into=["a","b",...])`. Capture
# the comma-separated quoted-string list so we can pull each value out.
_ENUM_INTO_RE = re.compile(r"into=\[(.*?)\]")
_QUOTED_STRING_RE = re.compile(r'"([^"]*)"')

# Tokens that break the vim/gedit syntax files generated downstream — see
# the legacy `create_Keywords.sh` filter (`sed -i "/{/d; /*/d; /#/d"`).
_UNSAFE_KEYWORD_CHARS = set("{}#*")


def _collect_keywords(content) -> list[str]:
    """Walk a `TRAD2Content` and return a sorted, deduplicated keyword list.

    Includes: every block name and its synonyms, every attribute name and
    its synonyms, every enum value found in `chaine(into=[...])` attribute
    types. Drops abstract `_base` / `_deriv` classes (they cannot appear
    in a dataset) and tokens containing characters that would break the
    downstream vim/gedit syntax files.
    """

    def _accept(kw: str) -> bool:
        if not kw:
            return False
        if kw.endswith(("_base", "_deriv")):
            return False
        return not any(c in _UNSAFE_KEYWORD_CHARS for c in kw)

    kws = set()
    for block in content.data:
        if _accept(block.name):
            kws.add(block.name)
        for s in block.synos:
            if _accept(s):
                kws.add(s)
        for attr in block.attrs:
            if _accept(attr.name):
                kws.add(attr.name)
            for s in attr.synos:
                if _accept(s):
                    kws.add(s)
            m = _ENUM_INTO_RE.search(attr.type)
            if m:
                for v in _QUOTED_STRING_RE.findall(m.group(1)):
                    if _accept(v):
                        kws.add(v)
    return sorted(kws)


def generate_keywords(
    out: str | Path = ".",
    projects: list[str] | None = None,
    trust_root: str | None = None,
    schema: str | Path | None = None,
) -> Path:
    """Generate `Keywords.txt` and `Keywords.Vim` into `out`.

    Source of truth is the `TRAD2_trustify` file in the schema directory
    (either `schema=` directly, or one freshly generated from
    `projects` / `trust_root`). The output format mirrors
    `doc/TRUST/Keywords.{txt,Vim}`:

    - `Keywords.txt`: one keyword per line, prefixed with `|`.
    - `Keywords.Vim`: a single `syntax keyword TRUSTLanguageKeywords  ...`
      line, keywords double-space-separated.

    Returns the output directory.
    """
    # $TRUST_ROOT fallback is centralised in generate_schema.
    schema_dir = Path(schema) if schema is not None else generate_schema(projects=projects, trust_root=trust_root)

    trad2_path = schema_dir / "TRAD2_trustify"
    if not trad2_path.exists():
        raise FileNotFoundError(f"TRAD2_trustify not found in schema dir {schema_dir}")
    content = TRAD2Content.BuildContentFromTRAD2(str(trad2_path))

    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    kws = _collect_keywords(content)
    (out_dir / "Keywords.txt").write_text("".join(f"|{k}\n" for k in kws))
    (out_dir / "Keywords.Vim").write_text("syntax keyword TRUSTLanguageKeywords  " + "  ".join(kws) + "\n")
    return out_dir


class _PreconditionError(ValueError):
    """Raised by ``api.generate_markdown`` for precondition violations (--out shape).

    Subclass of ``ValueError`` for back-compatibility; the CLI catches
    this specific type to route precondition errors to argparse's
    ``p.error`` (exit 2). Content errors (missing extras, duplicates,
    subdirs in extras) raise plain ``ValueError`` so they bubble
    through the generic exception handler instead → exit 1.
    """


def generate_markdown(
    out: str | Path,
    projects: list[str] | None = None,
    trust_root: str | None = None,
    schema: str | Path | None = None,
    kw_ref_path: str | Path | None = None,
) -> Path:
    """Generate the keyword reference manual into `out`.

    Steps:
      1. Validate `out` (must not be a file; must be missing or empty).
      2. Resolve `schema_dir` — `schema` if given, else `generate(...)`.
      3. Build the `OriginClassifier` from `projects` + `trust_root`.
      4. Build the `ExtrasIndex` from `projects` + `trust_root`
         (scans each ``<origin>/docs/trustify/extras/``).
      5. Run `DocGenerator(...).generate(out, kw_ref_path)`. The
         renderer resolves ``\\input`` + ``\\includeimage`` directives
         through the extras index.
      6. Copy every used image into `<out>/figures/`.
      7. Emit a stderr warning if any extras files went unused.

    Precondition violations (``--out`` shape) raise
    ``_PreconditionError`` which the CLI maps to exit code 2. All
    other failures (missing referenced file, duplicate, subdir,
    image-copy failure, schema import failure) bubble as their native
    exception type → exit 1.
    """
    import shutil
    import sys

    from trustify.doc.extras import ExtrasIndex
    from trustify.doc.generator import DocGenerator
    from trustify.doc.origin import OriginClassifier

    out_path = Path(out)
    if out_path.exists():
        if not out_path.is_dir():
            raise _PreconditionError(f"--out path is not a directory: {out_path}")
        if any(out_path.iterdir()):
            raise _PreconditionError(f"--out directory is not empty: {out_path}")
    else:
        out_path.mkdir(parents=True)

    # Resolve once and reuse for every callee — generate_schema does the
    # env-var fallback internally, but OriginClassifier / ExtrasIndex don't,
    # and they previously silently skipped the TRUST origin when only
    # $TRUST_ROOT (no --trust-root flag) was set. Same situation for
    # $project_directory: a baltik invocation must reach all three callees.
    projects = effective_projects(projects)
    trust_root = effective_trust_root(trust_root)
    schema_dir = Path(schema) if schema is not None else generate_schema(projects=projects, trust_root=trust_root)

    classifier = OriginClassifier.from_projects(projects, trust_root)
    extras = ExtrasIndex.from_projects(projects, trust_root)
    DocGenerator(
        schema_dir,
        classifier,
        extras,
    ).generate(out_path, kw_ref_path=Path(kw_ref_path) if kw_ref_path else None)

    # Copy every used image into <out>/figures/. Create the dir on
    # demand so it doesn't exist when nothing was resolved.
    used = extras.used_images()
    if used:
        figures_dst = out_path / "figures"
        figures_dst.mkdir(exist_ok=True)
        for src in used:
            shutil.copy2(src, figures_dst / src.name)

    # Warn about unused extras files. Each name is tagged with its
    # project of origin — any baltik may ship extras, so the origin is
    # what makes the warning actionable.
    unused = extras.unused_files()
    if unused:
        names = ", ".join(f"{p.name} ({extras.origin_of(p) or 'unknown project'})" for p in unused)
        print(
            f"trustify generate_markdown: warning: {len(unused)} unused extras file(s): {names}",
            file=sys.stderr,
        )

    return out_path


def generate_pdf(
    out: str | Path,
    from_markdown: str | Path | None = None,
    projects: list[str] | None = None,
    trust_root: str | None = None,
    schema: str | Path | None = None,
    doxygen: str | Path | None = None,
    keep_build: bool = False,
) -> Path:
    """Build a single PDF of the keyword reference manual.

    End-to-end by default: resolve the schema, run `generate_markdown`
    into a throwaway temp dir, then doxygen -> LaTeX -> PDF. Pass
    `from_markdown` to reuse an existing `generate_markdown` output and
    skip schema resolution + generation.

    `out` is the PDF file path (always required; its parent dir must
    already exist). Precondition violations (missing out-parent, an
    invalid `from_markdown` dir) raise `_PreconditionError` (CLI exit
    2); toolchain / build failures raise RuntimeError (CLI exit 1).

    Requires doxygen and a LaTeX toolchain (pdflatex, make) on PATH;
    both are probed up-front so a missing install fails fast.
    """
    import shutil
    import sys
    import tempfile

    from trustify import pdf as _pdf

    def _log(msg: str) -> None:
        # Progress goes to stderr so stdout stays just the final PDF
        # path (the CLI prints it there for scripting). flush so each
        # phase shows up immediately during the multi-minute passes.
        print(f"trustify generate_pdf: {msg}", file=sys.stderr, flush=True)

    out_path = Path(out).resolve()
    if not out_path.parent.is_dir():
        raise _PreconditionError(f"parent of --out does not exist: {out_path.parent}")

    # Probe the toolchain before the (slow) markdown + doxygen passes.
    doxygen_bin = _pdf.resolve_doxygen(str(doxygen) if doxygen else None)
    _pdf.probe_latex_toolchain()

    cleanup_md = False
    if from_markdown is not None:
        try:
            md_dir = _pdf.validate_markdown_dir(from_markdown)
        except ValueError as exc:
            raise _PreconditionError(str(exc)) from exc
        _log(f"[1/3] reusing markdown tree {md_dir}")
    else:
        projects = effective_projects(projects)
        trust_root = effective_trust_root(trust_root)
        md_dir = Path(tempfile.mkdtemp(prefix="trustify-pdf-md-"))
        cleanup_md = True
        _log("[1/3] generating markdown reference (schema + manual)...")
        generate_markdown(out=md_dir, projects=projects, trust_root=trust_root, schema=schema)

    build_dir = Path(tempfile.mkdtemp(prefix="trustify-pdf-build-"))
    try:
        _log("[2/3] running doxygen (this can take a few minutes)...")
        _pdf.run_doxygen(doxygen_bin, md_dir, build_dir)
        _log("[3/3] running LaTeX (pdflatex via make)...")
        _pdf.run_latex(build_dir)
        _pdf.finalize_pdf(build_dir, out_path)
        _log(f"wrote {out_path} ({out_path.stat().st_size // 1024} KiB)")
    finally:
        if keep_build:
            _log(f"kept build dir {build_dir}")
            if cleanup_md:
                _log(f"kept markdown dir {md_dir}")
        else:
            shutil.rmtree(build_dir, ignore_errors=True)
            if cleanup_md:
                shutil.rmtree(md_dir, ignore_errors=True)

    return out_path


def format_dataset_file(filename: str | Path) -> str:
    """Read `filename` and return its formatted content as a string.

    Convenience wrapper for the CLI. Programmatic callers operating on
    text should use `trustify.format_dataset(text)` directly.
    """
    from trustify.formatter import format_dataset as _fmt

    with open(filename, encoding="utf-8") as f:
        return _fmt(f.read())


@dataclass
class FormatFileResult:
    """Per-file outcome from `batch_format`.

    `diff` is the unified diff between the source content and the
    formatted content; empty when the file is already canonical.
    """

    path: str
    diff: str = ""


@dataclass
class BatchFormatResult:
    files: list[FormatFileResult]

    @property
    def changed_files(self) -> list[FormatFileResult]:
        return [f for f in self.files if f.diff]


def batch_format(data_files: list[str | Path], apply: bool = False) -> BatchFormatResult:
    """Format every `.data` file in `data_files`.

    Default is **dry-run**: read each file, compute the formatted
    output, and return a report carrying the unified diff per file
    that would change. Nothing is written to disk.

    `apply=True` switches to in-place rewrite. Each changed file is
    written atomically (temp-file + os.replace, matching modernize)
    so a SIGINT or disk-full mid-write leaves the source content
    untouched (audit 2.3 pattern). Audit finding 1.10 — destructive
    is opt-in, never the default.
    """
    import difflib

    from trustify.formatter import format_dataset as _fmt
    from trustify.modernize import _atomic_write_text

    results: list[FormatFileResult] = []
    for f in data_files:
        path = str(f)
        with open(path, encoding="utf-8") as fh:
            src = fh.read()
        new = _fmt(src)
        diff = ""
        if new != src:
            diff = "".join(
                difflib.unified_diff(
                    src.splitlines(keepends=True),
                    new.splitlines(keepends=True),
                    fromfile=path,
                    tofile=path,
                )
            )
            if apply:
                _atomic_write_text(path, new)
        results.append(FormatFileResult(path=path, diff=diff))
    return BatchFormatResult(files=results)


def modernize(
    projects: list[str] | None = None,
    trust_root: str | None = None,
    apply: bool = False,
):
    """Apply XD-tag modernization rules to the resolved workspace.

    Thin pass-through to `trustify.modernize.modernize` so the public
    API surface mirrors the other top-level entry points (`check`,
    `batch_check`, `generate_schema`, ...).
    """
    from trustify.modernize import modernize as _impl

    return _impl(projects=projects, trust_root=trust_root, apply=apply)
