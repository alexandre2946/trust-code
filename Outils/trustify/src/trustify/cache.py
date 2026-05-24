"""Content-hash cache for generated trustify parser/pydantic modules.

The cache key is computed from:
  1. The installed trustify package version. Without this, `pip install -U
     trustify` would not invalidate stale caches.
  2. The canonical text serialization of the in-memory TRAD2 content
     (i.e. what `TRAD2Content.toTRAD2Text()` would produce). This captures
     both the source content and the project ordering — both materially
     affect downstream pydantic/parser generation.
  3. The canonical JSON serialization of the per-keyword source-location
     map (`TRAD2Content.serialize_source_locations()`). Including this
     ensures the cache key is sensitive to intra-project source moves
     (same TRAD2 text, different file/line) but stays portable: two
     clones of the same project set at different absolute paths produce
     byte-identical JSON and therefore share one cache entry.

Default cache root: `$TRUSTIFY_CACHE_DIR` if set, else
`~/.cache/trustify/`. Per-invocation override via the `out=`
argument of `api.generate_schema` (wired to `--out` on the CLI).
"""

import fcntl
import hashlib
import json
import os
import shutil
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal


def default_cache_root() -> Path:
    """Where the cache lives by default.

    Reads `$TRUSTIFY_CACHE_DIR` lazily on every call so a test (or a
    `make` recipe) can override it without re-importing the module.
    Falls back to `~/.cache/trustify/`. Empty / whitespace-only env
    values fall back to the default too — silently honouring `""` would
    cause the cache to land in the process's CWD, which is rarely what
    a user setting the var to nothing actually meant.

    HPC sites with a quota'd / read-only `$HOME`, ephemeral CI
    runners, and multi-user shared TRUST installs use this to
    redirect the cache.
    """
    raw = os.environ.get("TRUSTIFY_CACHE_DIR", "")
    if raw.strip():
        return Path(os.path.expanduser(raw))
    return Path(os.path.expanduser("~")) / ".cache" / "trustify"


def compute_content_hash(content: Any, version: str) -> str:
    """Return a 16-char hex digest of:
      1. trustify version,
      2. canonical TRAD2 text (`content.toTRAD2Text()`),
      3. canonical source-locations JSON (`content.serialize_source_locations()`),
         or empty string when the method is absent (legacy stubs / tests).

    `content` must expose `toTRAD2Text()`. The
    `serialize_source_locations()` method is the new contract introduced
    by the source-locations refactor; falling back to an empty string
    keeps pre-existing test fakes that only mock `toTRAD2Text` working
    (their hash bucket then differs from real callers, which is fine).
    """
    h = hashlib.sha256()
    h.update(version.encode("utf-8"))
    h.update(b"\0")
    h.update(content.toTRAD2Text().encode("utf-8"))
    h.update(b"\0")
    locations_text = ""
    if hasattr(content, "serialize_source_locations"):
        locations_text = content.serialize_source_locations()
    h.update(locations_text.encode("utf-8"))
    return h.hexdigest()[:16]


def cache_dir_for(content: Any, version: str, override: Path | None = None, root: Path | None = None) -> Path:
    """Return the cache directory for this content+version pair, creating it
    if missing. If `override` is supplied (corresponds to CLI `--out`), use it
    instead of the hash-derived location. `root` lets tests redirect the
    default cache base away from `~/.cache/trustify/`.
    """
    if override is not None:
        out = Path(override)
    else:
        base = Path(root) if root is not None else default_cache_root()
        out = base / compute_content_hash(content, version)
    out.mkdir(parents=True, exist_ok=True)
    return out


@contextmanager
def cache_lock(mode: Literal["shared", "exclusive"], root: Path | None = None) -> Iterator[None]:
    """Hold an advisory POSIX `flock` on `<root>/.lock` for the duration
    of the with block. SHARED for readers (`is_cached` + the schema
    import in `core/misc_utilities.import_parser_module`), EXCLUSIVE
    for writers (`generate_schema` write + `trim_same_project_set`).

    Multiple SHARED holders coexist; an EXCLUSIVE holder excludes both
    other EXCLUSIVE holders and SHARED holders. Kernel-managed, so a
    process dying with the lock held releases it automatically — no
    stale-lock cleanup heuristics.

    The `.lock` sentinel file is created on first use and intentionally
    left behind on exit; removing it would race with another process
    trying to acquire.
    """
    base = Path(root) if root is not None else default_cache_root()
    base.mkdir(parents=True, exist_ok=True)
    lockfile = base / ".lock"
    op = fcntl.LOCK_SH if mode == "shared" else fcntl.LOCK_EX
    fd = os.open(str(lockfile), os.O_CREAT | os.O_RDWR, 0o644)
    try:
        fcntl.flock(fd, op)
        try:
            yield
        finally:
            fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)


@contextmanager
def atomic_write_dir(final_dir: Path) -> Iterator[Path]:
    """Write a cache entry into a sibling temp dir, then atomically
    rename it onto `final_dir`.

    Guarantee: the contents of `final_dir` go from "old or absent" to
    "complete new entry" in a single atomic `os.replace` — readers
    racing the write never see a half-written tree. A SIGKILL between
    the schema's first and last write leaves a `.<name>.tmp.<pid>/`
    sibling but never a half-written `final_dir`.

    On exception inside the `with` block the temp dir is removed and
    `final_dir` is left untouched.
    """
    final_dir = Path(final_dir)
    temp = final_dir.parent / f".{final_dir.name}.tmp.{os.getpid()}"
    if temp.exists():
        shutil.rmtree(temp)
    temp.mkdir(parents=True)
    try:
        yield temp
    except BaseException:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    # POSIX `rename(2)` (and Python's `os.replace`) refuses to overwrite
    # a non-empty dir. The is_cached gate upstream ensures we only reach
    # here when no COMPLETE entry exists at `final_dir`; if a partial
    # leftover (from a previous crash) is present, swap it for the fresh
    # entry. This brief rmtree/replace gap is closed by the exclusive
    # cache lock in `generate_schema` (see audit finding 2.2).
    if final_dir.exists():
        shutil.rmtree(final_dir)
    os.replace(temp, final_dir)


# Pyd module name pattern. Three call sites used to construct this
# string literally — api.generate_schema (write), is_cached (read via
# glob), and the test/fixture path in core.misc_utilities — and CLAUDE.md
# called the lockstep out as a hazard. Centralised here so the
# write-side filename and the read-side glob can't drift.
_PYD_MODULE_PREFIX = "trustify_gen_pyd_"
PYD_MODULE_GLOB = f"{_PYD_MODULE_PREFIX}*.py"


def pyd_module_filename(digest: str) -> str:
    """Filename for the digest-suffixed pyd module: `trustify_gen_pyd_<digest>.py`.

    The digest is what makes two concurrent schemas occupy distinct
    `sys.modules` slots (each generated parser module's
    `from ... import *` references this name). Match-via-glob on the
    read side uses `PYD_MODULE_GLOB`.
    """
    return f"{_PYD_MODULE_PREFIX}{digest}.py"


def is_cached(directory: Path) -> bool:
    """Return True when `directory` holds a COMPLETE cache entry.

    Three sentinels must all be present:
      - `trustify_gen.py` (the parser module),
      - `trustify_gen_pyd_<digest>.py` (the pydantic module — glob-based
        so concurrent schemas in distinct `sys.modules` slots all match),
      - `provenance.json` (written last by `generate_schema` — its
        presence proves the write was not interrupted).

    The sentinel rule is what makes interrupted writes safe: a process
    killed between writing `trustify_gen.py` and `provenance.json`
    leaves a directory that future runs correctly treat as not-cached.
    """
    d = Path(directory)
    return any(d.glob(PYD_MODULE_GLOB)) and (d / "trustify_gen.py").exists() and (d / "provenance.json").exists()


# ---------------------------------------------------------------------------
# Cache inspection / manipulation API (consumed by `trustify cache`)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CacheEntry:
    """One generated-schema directory under the cache root.

    `created_at`, `trustify_version`, `projects`,
    `generation_time_seconds`, and `generation_steps_seconds` all come
    from `provenance.json` (written last by `api.generate_schema`, and
    required by `is_cached`). When fields are absent from that file —
    e.g. legacy provenance written by an older trustify — the values
    are `None` / `{}` and the CLI renders them as `(unknown)`.
    """

    id: str  # 16-char content hash (dir basename)
    path: Path  # absolute path to the entry dir
    created_at: str | None  # ISO timestamp or None
    trustify_version: str | None
    projects: dict[str, str] = field(default_factory=dict)
    size_bytes: int = 0
    generation_time_seconds: float | None = None
    generation_steps_seconds: dict[str, float] = field(default_factory=dict)


def _entry_size(path: Path) -> int:
    """Sum of file sizes (bytes) under `path`. Symlinks are not followed —
    `lstat` reports the link slot itself, never the target (audit 5.5).
    `os.path.getsize` would dereference, so a broken file symlink raised
    OSError (silently caught → understated total) and a live symlink
    double-counted (target's bytes once for the regular file, once for
    the link).
    """
    total = 0
    for root, _dirs, files in os.walk(path, followlinks=False):
        for name in files:
            # Per-file tolerance: one bad entry must not abort the whole sum.
            try:
                total += os.lstat(os.path.join(root, name)).st_size
            except OSError:  # noqa: PERF203
                # Concurrent delete between os.walk's readdir and lstat —
                # ignore and continue.
                continue
    return total


def _load_entry(entry_dir: Path) -> CacheEntry:
    """Build a `CacheEntry` from a cache-dir path. Reads `provenance.json`
    when present; falls back to `None` fields when absent (older entries).
    """
    prov_path = entry_dir / "provenance.json"
    created_at: str | None = None
    trustify_version: str | None = None
    projects: dict[str, str] = {}
    generation_time_seconds: float | None = None
    generation_steps_seconds: dict[str, float] = {}
    if prov_path.is_file():
        try:
            data = json.loads(prov_path.read_text())
        except (OSError, json.JSONDecodeError):
            data = {}
        ca = data.get("generated_at")
        if isinstance(ca, str):
            created_at = ca
        tv = data.get("trustify_version")
        if isinstance(tv, str):
            trustify_version = tv
        proj = data.get("projects")
        if isinstance(proj, dict):
            projects = {str(k): str(v) for k, v in proj.items()}
        gt = data.get("generation_time_seconds")
        if isinstance(gt, (int, float)):
            generation_time_seconds = float(gt)
        steps = data.get("generation_steps_seconds")
        if isinstance(steps, dict):
            generation_steps_seconds = {str(k): float(v) for k, v in steps.items() if isinstance(v, (int, float))}
    return CacheEntry(
        id=entry_dir.name,
        path=entry_dir.resolve(),
        created_at=created_at,
        trustify_version=trustify_version,
        projects=projects,
        size_bytes=_entry_size(entry_dir),
        generation_time_seconds=generation_time_seconds,
        generation_steps_seconds=generation_steps_seconds,
    )


def iter_entries(root: Path | None = None) -> list[CacheEntry]:
    """Return every cache entry under `root` (default
    `~/.cache/trustify/`), sorted newest first by `created_at`.

    Only directories containing both `trustify_gen.py` and a
    digest-suffixed `trustify_gen_pyd_*.py` count as entries.
    Pre-refactor caches that lack `provenance.json` still appear; their
    timestamps sort last.
    """
    base = Path(root) if root is not None else default_cache_root()
    if not base.is_dir():
        return []
    out: list[CacheEntry] = []
    for child in sorted(base.iterdir()):
        if not child.is_dir():
            continue
        if not is_cached(child):
            continue
        out.append(_load_entry(child))
    # Descending sort by ISO timestamp. None / missing timestamps map to
    # "" via `or ""` and sort last (any real ISO date > "" lexically),
    # giving the documented newest-first ordering. Audit 6.2.
    out.sort(key=lambda e: e.created_at or "", reverse=True)
    return out


class CacheEntryNotFound(LookupError):
    """No cache entry matches the requested ID prefix."""


class AmbiguousCacheEntry(LookupError):
    """Multiple cache entries match the requested ID prefix."""

    def __init__(self, prefix: str, matches: list[str]):
        self.prefix = prefix
        self.matches = list(matches)
        super().__init__(f"ambiguous cache id prefix {prefix!r}: matches {', '.join(sorted(matches))}")


def resolve_id(prefix: str, root: Path | None = None) -> Path:
    """Resolve `prefix` to a single cache entry directory.

    Accepts any unique prefix of an entry's 16-char hash (git
    short-SHA style). Empty prefix never matches. Raises
    `CacheEntryNotFound` when no entry matches, `AmbiguousCacheEntry`
    when more than one does.
    """
    if not prefix:
        raise CacheEntryNotFound("empty cache id prefix")
    base = Path(root) if root is not None else default_cache_root()
    if not base.is_dir():
        raise CacheEntryNotFound(f"no cache entry matches {prefix!r}")
    matches = [d for d in base.iterdir() if d.is_dir() and is_cached(d) and d.name.startswith(prefix)]
    if not matches:
        raise CacheEntryNotFound(f"no cache entry matches {prefix!r}")
    if len(matches) > 1:
        raise AmbiguousCacheEntry(prefix, [m.name for m in matches])
    return matches[0].resolve()


def _projects_key(projects: dict[str, str]) -> str:
    """Canonical hash-stable representation of a `projects` mapping
    (`{name: abs_path}`). Used to group cache entries by the project
    set they were generated from.
    """
    return json.dumps(projects, sort_keys=True, separators=(",", ":"))


def trim_same_project_set(keep_path: Path, root: Path | None = None) -> list[str]:
    """Delete every cache entry under `root` other than `keep_path`
    whose `provenance.projects` mapping exactly matches `keep_path`'s.

    Used as the auto-trim hook in `api.generate_schema`: each freshly
    written entry supersedes prior entries built from the same
    project set (typical during active development where only source
    line numbers shift between regenerations).

    Returns the IDs of the deleted entries. Entries with missing or
    unreadable `provenance.json` (pre-refactor caches, or the new
    entry's siblings that haven't recorded provenance) are left
    alone — no provenance, no match.
    """
    import shutil

    keep_path = Path(keep_path).resolve()
    keep_entry = _load_entry(keep_path)
    if not keep_entry.projects:
        return []
    target_key = _projects_key(keep_entry.projects)
    deleted: list[str] = []
    base = Path(root) if root is not None else default_cache_root()
    if not base.is_dir():
        return []
    for child in base.iterdir():
        if not child.is_dir() or not is_cached(child):
            continue
        if child.resolve() == keep_path:
            continue
        entry = _load_entry(child)
        if not entry.projects:
            continue
        if _projects_key(entry.projects) != target_key:
            continue
        shutil.rmtree(child)
        deleted.append(entry.id)
    return deleted


def trim_all_duplicates(root: Path | None = None) -> dict[str, list[str]]:
    """Walk every cache entry, group by `provenance.projects`, keep the
    newest per group (by `generated_at`; ties broken by entry id), and
    delete the rest. Returns `{kept_id: [deleted_ids]}` for groups that
    actually had duplicates.

    Entries without a `provenance.projects` mapping (legacy / corrupt)
    are excluded from grouping and left untouched.
    """
    import shutil

    groups: dict[str, list[CacheEntry]] = {}
    for entry in iter_entries(root):
        if not entry.projects:
            continue
        groups.setdefault(_projects_key(entry.projects), []).append(entry)
    result: dict[str, list[str]] = {}
    for entries in groups.values():
        if len(entries) <= 1:
            continue
        # iter_entries already sorts newest-first; the first element is
        # the keeper.
        keeper, *to_delete = entries
        deleted: list[str] = []
        for entry in to_delete:
            shutil.rmtree(entry.path)
            deleted.append(entry.id)
        if deleted:
            result[keeper.id] = deleted
    return result


def should_auto_trim(start_dir: Path | None = None) -> bool:
    """Return whether `api.generate_schema` should auto-trim same-project-set
    cache entries after writing a new one.

    Defaults to True. Disabled when:
      - environment variable ``TRUSTIFY_NO_AUTO_TRIM`` is set to any
        non-empty value (matching NO_COLOR's convention used in
        `_color.py` — presence-of-a-value, no deny-list parsing, so
        `TRUSTIFY_NO_AUTO_TRIM=0` disables, not the surprising reverse);
      - the nearest `.trustify.json` (walking up from `start_dir`) has
        `cache.auto_trim: false`.
    """
    if os.environ.get("TRUSTIFY_NO_AUTO_TRIM"):
        return False
    try:
        from trustify.projects import ConfigError, discover_config
    except ImportError:
        return True
    try:
        cfg = discover_config(start_dir)
    except (ConfigError, ValueError, OSError):
        return True
    return True if cfg is None else cfg.cache.auto_trim


def copy_entry(src: Path, dest_dir: Path) -> Path:
    """Copy the generated files from `src` (a cache entry dir) into
    `dest_dir`. `__pycache__` directories and `*.pyc` files are
    skipped. Raises `FileExistsError` when `dest_dir` exists and is
    non-empty (the existing dest is left untouched).

    Returns the path of the populated destination.
    """
    src = Path(src)
    dest_dir = Path(dest_dir)
    if dest_dir.exists():
        if not dest_dir.is_dir():
            raise FileExistsError(f"{dest_dir} exists and is not a directory")
        if any(dest_dir.iterdir()):
            raise FileExistsError(f"{dest_dir} exists and is not empty")
    else:
        dest_dir.mkdir(parents=True)
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    for item in src.iterdir():
        if item.name == "__pycache__" or item.name.endswith(".pyc"):
            continue
        target = dest_dir / item.name
        if item.is_dir():
            shutil.copytree(item, target, ignore=ignore, symlinks=False)
        else:
            shutil.copy2(item, target)
    return dest_dir
