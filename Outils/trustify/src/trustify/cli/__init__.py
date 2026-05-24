# PYTHON_ARGCOMPLETE_OK
"""trustify CLI — thin argparse dispatcher over `trustify.api`.

Entry point registered in `pyproject.toml` as console_scripts `trustify`.
Also runnable as `python -m trustify.cli ...`.

Project-resolution flags (`--projects`, `--trust-root`, `--schema`)
are GLOBAL: they live on the top-level parser and must come BEFORE the
subcommand. `--schema` is mutually exclusive with both `--projects` and
`--trust-root` (manual check; argparse mutex groups handle only one
set).
"""

import argparse
import fnmatch
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path

from trustify import api
from trustify.core.misc_utilities import TrustifyInternalError
from trustify.projects import (
    ConfigError,
    TrustifyConfig,
    detect_baltik_root,
    discover_config,
    effective_projects,
    effective_trust_root,
    expand_dependencies,
)

# Subcommands that actually consume a pre-generated schema directory.
# `--schema` is rejected on any other command (e.g. generate_schema would
# overwrite what `--schema` points at; format is pure text).
_SCHEMA_CONSUMERS = {
    "check",
    "batch-check",
    "batch-check-projects",
    "generate_markdown",
    "generate_keywords",
    "generate_pdf",
}

# Subcommands that need the full four-signal workspace resolution
# (config discovery + baltik auto-detection + env fallback + dependency
# expansion). init-config is included because its job is to *freeze* the
# resolved tree into a file: same resolution as the schema-consuming
# commands, just written to disk instead of fed to the parser. `projects`
# is included because it dumps the resolved list for downstream scripts.
_NEEDS_PROJECT_RESOLUTION = {
    "generate_schema",
    "check",
    "batch-check",
    "batch-check-projects",
    "generate_keywords",
    "generate_markdown",
    "init-config",
    "projects",
}
# `modernize` deliberately opts OUT of the four-signal resolver: it
# rewrites .cpp / .xd sources in place, so its scope is narrower than
# schema-consuming commands (see `_scope_modernize` for the rules).


@dataclass(frozen=True)
class _CliResolution:
    """Result of `_resolve_project_args`.

    `projects` is the post-expansion list passed to `api.*` calls.
    `trust_root` is the resolved trust source tree (or None).
    `primaries` is the **pre-expansion** project list — the explicit
    projects the user asked for via CLI flag / `.trustify.json` /
    baltik auto-detection / env var, before transitive `[dependencies]`
    were folded in. The `projects` subcommand uses it to tag each entry
    of `projects` as either `primary` (was in this set) or `dependency`
    (added by `expand_dependencies`).

    `api_kwargs` returns the `(projects, trust_root)` pair callers
    forward to `api.*` — it deliberately omits `primaries`, which is a
    CLI-internal detail.
    """

    projects: list[str] | None
    trust_root: str | None
    primaries: list[str] | None

    @property
    def api_kwargs(self) -> dict:
        return {"projects": self.projects, "trust_root": self.trust_root}


def _get_version() -> str:
    from trustify import __version__

    return __version__


_PROJECT_CTX_NOTE = (
    "Project context — these GLOBAL flags MUST come BEFORE the subcommand:\n"
    "  -p / --projects DIR    (repeat for multiple, defaults to $project_directory)\n"
    "  --trust-root DIR       defaults to $TRUST_ROOT\n"
    "  --schema DIR           reuse an already-generated schema\n"
    "If neither flag is given, the nearest .trustify.json (walking up from CWD)\n"
    "is consulted; failing that, $project_directory and $TRUST_ROOT from the\n"
    "environment are used as default --projects / --trust-root values."
)

_FORMAT_NOTE = (
    "Pure text manipulation — no schema/project context needed.\n"
    "Accepts a single file. Default: print to stdout. Use --out to\n"
    "write to a chosen path, or --in-place to rewrite the file."
)

_BATCH_FORMAT_NOTE = (
    "Pure text manipulation — no schema/project context needed.\n"
    "Reformats every listed file in place. For single-file formatting\n"
    "with stdout / --out / --in-place flags, use `trustify format`."
)

_PATH_EXPANSION_NOTE = (
    "Each DATA_PATH is either a .data file or a directory; directories\n"
    "are recursed for *.data files. Use --exclude PATTERN (repeatable) to\n"
    "skip files or subdirectories: a pattern without '/' matches basenames\n"
    "anywhere in the tree, with '/' matches a path relative to the listed\n"
    "root, and a trailing '/' restricts the match to directories. Patterns\n"
    "use fnmatch glob syntax (*, ?, [abc]). Explicit files on the command\n"
    "line bypass --exclude."
)


_ENV_HELP = """\
Environment variables that affect trustify's runtime behavior:

  TRUST_ROOT               Path to a checked-out TRUST source tree. Used
                           as the schema root when --trust-root is absent.
  project_directory        Legacy lowest-priority `--projects` fallback.
                           Overridden by --projects, .trustify.json, and
                           baltik auto-detection (project.cfg walk-up).
  TRUSTIFY_CACHE_DIR       Override the schema cache root
                           (default: ~/.cache/trustify/). Use on HPC sites
                           with a quota'd $HOME or for shared installs.
  TRUSTIFY_NO_AUTO_TRIM    Any non-empty value disables the dedup-by-
                           project-set trim that `generate_schema` runs
                           after each write. Same as .trustify.json
                           `cache.auto_trim: false`.
  TRUSTIFY_FORCE_FORMAT    Any non-empty value makes `batch-format` apply
                           in place by default (equivalent to --apply).
  TRUSTIFY_MP_START_METHOD Override the multiprocessing start method for
                           `batch-check` workers (default: forkserver on
                           Linux, spawn on Windows). `fork` risks deadlock.
  NO_COLOR                 Any non-empty value disables ANSI colors
                           (https://no-color.org/ convention).
"""


class _HelpEnvAction(argparse.Action):
    """Print _ENV_HELP and exit, mirroring how argparse's own `--version` works."""

    def __init__(self, option_strings, dest=argparse.SUPPRESS, default=argparse.SUPPRESS, help=None):
        super().__init__(option_strings=option_strings, dest=dest, default=default, nargs=0, help=help)

    def __call__(self, parser, _namespace, _values, _option_string=None):
        parser.exit(message=_ENV_HELP)


def _make_common() -> argparse.ArgumentParser:
    """Parent parser holding the global project-resolution flags.

    Attached only to the main parser (never to subparsers): the flags
    must come BEFORE the subcommand name and cannot be silently
    overridden later in the command line.

    `action="append"` for `--projects` avoids the `nargs="*"` greedy-eat
    bug where `-p A check` would consume `check` as a value. Syntax is
    `-p A -p B` (one value per occurrence).
    """
    common = argparse.ArgumentParser(add_help=False)
    g = common.add_argument_group("project resolution")
    g.add_argument(
        "--projects",
        "-p",
        action="append",
        metavar="DIR",
        help="Project directory. Repeat for multiple (`-p A -p B`). Combined with --trust-root.",
    )
    g.add_argument("--trust-root", metavar="DIR", help="TRUST source tree. Defaults to $TRUST_ROOT.")
    g.add_argument(
        "--schema",
        metavar="DIR",
        help="Reuse an already-generated schema. Mutually exclusive with --projects and --trust-root.",
    )
    g.add_argument(
        "--no-auto-deps",
        dest="no_auto_deps",
        action="store_true",
        default=None,
        help=(
            "Do not auto-resolve transitive [dependencies] from each project's project.cfg. "
            "By default trustify mirrors baltik_configure and pulls every declared dependency in. "
            "Equivalent to setting `dependencies.auto_resolve: false` in .trustify.json."
        ),
    )
    return common


def _add_batch_check_opts(sp: argparse.ArgumentParser) -> None:
    """Add the options shared by `batch-check` and `batch-check-deps`
    (everything except how the input paths are obtained)."""
    sp.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="PATTERN",
        help="Glob pattern to skip while recursing directories (repeatable).",
    )
    sp.add_argument(
        "--no-skip", action="store_true", help="Check every dataset even if it carries a 'TRUSTIFY NOT' opt-out marker."
    )
    sp.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress the per-dataset progress line (printed by default).",
    )
    sp.add_argument(
        "--summary-out",
        metavar="PATH",
        help="Also write the final summary (Summary line + SKIPPED/FAILED list + obsolete-marker warnings) to PATH.",
    )
    sp.add_argument(
        "--jobs",
        "-j",
        type=int,
        default=1,
        metavar="N",
        help="Run per-file checks across N worker processes (default: 1, sequential). "
        "N>=2 uses a fork-based multiprocessing.Pool; workers inherit the "
        "already-loaded schema, so there is no per-worker re-import cost. "
        "N=0 means 'use all available CPUs' (resolved via os.cpu_count()).",
    )
    sp.add_argument(
        "--follow-symlinks",
        action="store_true",
        help=(
            "Follow symbolic links during directory recursion (both dir "
            "and file symlinks). Default skips them — a symlinked file "
            "at a leaf emits a one-line skip warning to stderr so "
            "out-of-scope links don't silently sneak into the check."
        ),
    )


def _build_parser() -> argparse.ArgumentParser:
    common = _make_common()

    p = argparse.ArgumentParser(
        prog="trustify",
        description="TRUST dataset tooling: schema generation, validation, format conversion.",
        epilog="Run `trustify <command> --help` for command-specific help, or "
        "`trustify --help-env` for environment variables. "
        "Project-resolution flags (--projects / --trust-root / "
        "--schema) must be specified BEFORE the subcommand.",
        parents=[common],
    )
    p.add_argument("--version", action="version", version=_get_version())
    p.add_argument(
        "--help-env",
        action=_HelpEnvAction,
        help="Show environment variables that affect trustify and exit.",
    )

    sub = p.add_subparsers(dest="command", required=True, metavar="<command>")
    fmt = argparse.RawDescriptionHelpFormatter

    # generate_schema
    sp = sub.add_parser(
        "generate_schema",
        formatter_class=fmt,
        description=_PROJECT_CTX_NOTE,
        help="Generate the parser/pydantic schema.",
    )
    sp.set_defaults(cmd="generate_schema")
    sp.add_argument("--out", "-o", metavar="DIR", help="Write dir. Defaults to ~/.cache/trustify/<hash>.")

    # check
    sp = sub.add_parser(
        "check",
        formatter_class=fmt,
        description=_PROJECT_CTX_NOTE,
        help="Validate one TRUST dataset against the schema.",
    )
    sp.set_defaults(cmd="check")
    sp.add_argument("data_file", metavar="DATA_FILE")
    sp.add_argument(
        "--no-skip", action="store_true", help="Check the dataset even if it carries a 'TRUSTIFY NOT' opt-out marker."
    )

    # batch-check
    sp = sub.add_parser(
        "batch-check",
        formatter_class=fmt,
        description=_PROJECT_CTX_NOTE + "\n\n" + _PATH_EXPANSION_NOTE,
        help="Validate multiple datasets and report an aggregate.",
    )
    sp.set_defaults(cmd="batch-check")
    sp.add_argument("data_paths", nargs="+", metavar="DATA_PATH")
    _add_batch_check_opts(sp)

    # batch-check-projects
    sp = sub.add_parser(
        "batch-check-projects",
        formatter_class=fmt,
        description=(
            _PROJECT_CTX_NOTE + "\n\n"
            "Run batch-check on a subdirectory of each resolved project. The set of "
            "projects is the same one `trustify projects` dumps, and `--only` filters "
            "it by role exactly like that command (trust_root, primary, dependency; "
            "no `--only` means every resolved project). Each selected project root has "
            "`--subdir` appended; projects lacking that directory are skipped with a "
            "stderr note.\n\n"
            "Note: `--exclude` here is the batch-check file-glob filter (below), NOT a "
            "role filter — use `--only` for roles.\n\n" + _PATH_EXPANSION_NOTE
        ),
        help="Validate the datasets shipped by resolved projects (filter with --only).",
    )
    sp.set_defaults(cmd="batch-check-projects")
    sp.add_argument(
        "--only",
        metavar="ROLE[,ROLE...]",
        help="Keep only projects with one of these roles (trust_root, primary, "
        "dependency), like `trustify projects --only`. Default: all resolved projects. "
        "E.g. --only dependency.",
    )
    sp.add_argument(
        "--subdir",
        metavar="REL_PATH",
        default="tests/Reference",
        help="Directory appended to each selected project root to locate datasets "
        "(default: tests/Reference). E.g. --subdir tests/Validation.",
    )
    _add_batch_check_opts(sp)

    # format
    sp = sub.add_parser("format", formatter_class=fmt, description=_FORMAT_NOTE, help="Format a single .data file.")
    sp.set_defaults(cmd="format")
    sp.add_argument("data_file", metavar="DATA_FILE")
    out_or_inplace = sp.add_mutually_exclusive_group()
    out_or_inplace.add_argument("--out", "-o", metavar="PATH", help="Write formatted output to PATH.")
    out_or_inplace.add_argument("--in-place", "-i", action="store_true", help="Edit the file in place.")

    # batch-format
    sp = sub.add_parser(
        "batch-format",
        formatter_class=fmt,
        description=_BATCH_FORMAT_NOTE + "\n\n" + _PATH_EXPANSION_NOTE,
        help="Format multiple .data files in place.",
    )
    sp.set_defaults(cmd="batch-format")
    sp.add_argument("data_paths", nargs="+", metavar="DATA_PATH")
    sp.add_argument(
        "--exclude",
        action="append",
        default=[],
        metavar="PATTERN",
        help="Glob pattern to skip while recursing directories (repeatable).",
    )
    sp.add_argument(
        "--apply",
        action="store_true",
        help=(
            "Write changes in place (default is dry-run: print the diff "
            "and exit 1 when changes are pending). "
            "TRUSTIFY_FORCE_FORMAT in the env (any non-empty value) has "
            "the same effect."
        ),
    )
    sp.add_argument(
        "--follow-symlinks",
        action="store_true",
        help=(
            "Follow symbolic links during directory recursion (both dir "
            "and file symlinks). Default skips them — a symlinked file "
            "at a leaf emits a one-line skip warning to stderr."
        ),
    )

    # modernize
    sp = sub.add_parser(
        "modernize",
        formatter_class=fmt,
        description=(
            "Rewrites legacy XD-tag forms in the current top-level project's "
            "`<project>/src/` tree: numeric brace flags (-3/-2/-1/0/1) -> "
            "BRACE/NO_BRACE/INHERITS_BRACE, numeric opt flags (0/1) -> REQ/OPT, "
            "and lines > 120 chars split via XD_CONT. Dry-run by default "
            "(exits 1 if changes pending, like `ruff format --check`); `--apply` "
            "writes the rewrite.\n\n"
            "Scope (narrower than schema-consuming commands — modernize rewrites "
            "sources in place):\n"
            "  --projects A [B ...]    those exact projects, no transitive deps, "
            "no trust_root overlay; a warning is emitted at the end.\n"
            "  (default)               the auto-detected baltik — project.cfg "
            "walk-up from CWD, or $project_directory env. Ignores --trust-root.\n"
            "  --trust-root T (alone)  modernize the TRUST tree itself, only when "
            "no baltik resolves.\n"
            "`.trustify.json` workspace.projects / workspace.trust_root is "
            "rejected — modernize requires a developer environment, not a "
            "user-environment config."
        ),
        help="Rewrite legacy XD-tag forms in the current top-level project's `src/` tree.",
    )
    sp.set_defaults(cmd="modernize")
    sp.add_argument(
        "--apply",
        action="store_true",
        help="Write the rewrite (default: dry-run with diff on stdout).",
    )

    # init-config
    sp = sub.add_parser(
        "init-config",
        formatter_class=fmt,
        description=_PROJECT_CTX_NOTE + "\n\nFreezes the currently-resolved workspace "
        "(--projects / --trust-root flags, $TRUST_ROOT / $project_directory env vars, "
        "baltik auto-detection via project.cfg, transitive [dependencies]) into the "
        "target directory's .trustify.json. Future invocations from there will use the "
        "frozen list without re-resolving — handy for decoupling a dataset folder from "
        "whatever env happens to be active.",
        help="Freeze the currently-resolved workspace into a .trustify.json file.",
    )
    sp.set_defaults(cmd="init-config")
    sp.add_argument("directory", nargs="?", default=".", metavar="DIR", help="Target directory (default: cwd).")
    sp.add_argument("--force", action="store_true", help="Overwrite an existing .trustify.json.")

    # generate_markdown
    sp = sub.add_parser(
        "generate_markdown",
        formatter_class=fmt,
        description=_PROJECT_CTX_NOTE,
        help="Generate the keyword reference manual.",
    )
    sp.set_defaults(cmd="generate_markdown")
    sp.add_argument(
        "--out",
        "-o",
        required=True,
        metavar="DIR",
        help="Output directory. Must not exist, or exist as an empty directory.",
    )
    sp.add_argument(
        "--kw-ref-path",
        metavar="PATH",
        help="Optional override for keyword_reference.md location. Default: <out>/keyword_reference.md.",
    )

    # generate_pdf
    sp = sub.add_parser(
        "generate_pdf",
        formatter_class=fmt,
        description=(
            _PROJECT_CTX_NOTE + "\n\n"
            "Build a single hyperlinked PDF of the keyword reference manual. "
            "By default runs the whole pipeline (schema -> markdown -> doxygen "
            "-> LaTeX -> PDF) using throwaway temp directories. Pass "
            "--from-markdown to reuse an existing `generate_markdown` output and "
            "skip schema resolution + generation.\n\n"
            "Requires doxygen and a LaTeX toolchain (pdflatex, make) on PATH; "
            "both are probed up-front so a missing install fails fast."
        ),
        help="Build a PDF of the keyword reference manual.",
    )
    sp.set_defaults(cmd="generate_pdf")
    sp.add_argument(
        "--out",
        "-o",
        required=True,
        metavar="PATH.pdf",
        help="PDF output path (a file). The parent directory must already exist.",
    )
    sp.add_argument(
        "--from-markdown",
        metavar="DIR",
        default=None,
        help="Reuse an existing `trustify generate_markdown --out` directory; "
        "skips schema resolution and markdown generation.",
    )
    sp.add_argument(
        "--keep-build",
        action="store_true",
        help="Keep the scratch build dir (and the temp markdown dir in "
        "end-to-end mode) for inspection; their paths are printed to stderr.",
    )
    sp.add_argument(
        "--doxygen",
        default=None,
        metavar="PATH",
        help="Path to the doxygen binary (overrides $TRUST_DOXYGEN_BINARY and PATH).",
    )

    # generate_keywords
    sp = sub.add_parser(
        "generate_keywords",
        formatter_class=fmt,
        description=_PROJECT_CTX_NOTE,
        help="Generate Keywords.txt / Keywords.Vim.",
    )
    sp.set_defaults(cmd="generate_keywords")
    sp.add_argument("--out", "-o", metavar="DIR", default=".", help="Output directory (default: cwd).")

    # projects
    sp = sub.add_parser(
        "projects",
        formatter_class=fmt,
        description=(
            _PROJECT_CTX_NOTE + "\n\n"
            "Dump the resolved project list, in overlay order (trust_root first if any, "
            "then transitive [dependencies], then primary baltiks last). Output is one path "
            "per line by default, intended for shell scripting:\n"
            "  e.g. `BALTIK_DEPENDENCIES := $(shell trustify projects --only=dependency)`.\n\n"
            "Roles: `trust_root` (the TRUST source tree), `primary` (explicitly listed by "
            "the user — CLI flag, .trustify.json, or auto-detection), `dependency` (pulled "
            "in transitively via project.cfg [dependencies])."
        ),
        help="List resolved project paths (for downstream scripts).",
    )
    sp.set_defaults(cmd="projects")
    sp.add_argument("--json", dest="as_json", action="store_true", help="Emit JSON instead of one-path-per-line.")
    sp.add_argument(
        "--only",
        metavar="ROLE[,ROLE...]",
        help="Keep only entries with one of these roles (trust_root, primary, dependency). "
        "Mutually exclusive with --exclude.",
    )
    sp.add_argument(
        "--exclude",
        metavar="ROLE[,ROLE...]",
        help="Drop entries with one of these roles. Mutually exclusive with --only.",
    )

    # cache
    sp = sub.add_parser(
        "cache",
        formatter_class=fmt,
        description="Inspect and manage the trustify schema cache (~/.cache/trustify/).",
        help="Inspect and manage the trustify schema cache.",
    )
    sp.set_defaults(cmd="cache")
    cache_sub = sp.add_subparsers(dest="cache_action", required=True, metavar="<action>")
    cl = cache_sub.add_parser(
        "list",
        formatter_class=fmt,
        description="List all cached schemas with creation time, version, projects, and on-disk size.",
        help="List cache entries.",
    )
    cl.add_argument(
        "--json", dest="as_json", action="store_true", help="Emit machine-readable JSON instead of the table."
    )
    cc = cache_sub.add_parser(
        "clean",
        formatter_class=fmt,
        description="Delete cache entries. With no IDs, "
        "prompts before wiping every entry "
        "(use --force to skip the prompt). "
        "With explicit IDs, deletes only "
        "those entries without prompting.",
        help="Delete cache entries.",
    )
    cc.add_argument(
        "ids", nargs="*", metavar="ID", help="Cache entry IDs (any unique prefix). Omit to target all entries."
    )
    cc.add_argument("--force", "-f", action="store_true", help="Skip the confirmation prompt when wiping all entries.")
    cp = cache_sub.add_parser(
        "copy",
        formatter_class=fmt,
        description="Copy a cache entry's generated "
        "files into DEST_DIR. DEST_DIR is "
        "created if absent and must be "
        "empty if it already exists. "
        "`__pycache__` and `*.pyc` are "
        "filtered out.",
        help="Copy a cache entry to a directory.",
    )
    cp.add_argument("id", metavar="ID", help="Cache entry ID (any unique prefix).")
    cp.add_argument("dest_dir", metavar="DEST_DIR", help="Destination directory.")
    cph = cache_sub.add_parser(
        "path",
        formatter_class=fmt,
        description="Print the absolute path of a "
        "cache entry. Useful in shell "
        "scripting: `--schema $(trustify "
        "cache path <id>)`.",
        help="Print absolute path of a cache entry.",
    )
    cph.add_argument("id", metavar="ID", help="Cache entry ID (any unique prefix).")
    ctr = cache_sub.add_parser(
        "trim",
        formatter_class=fmt,
        description="Deduplicate cache entries: group "
        "by the project set they were "
        "generated from and delete every "
        "entry but the newest in each "
        "group. The same logic runs "
        "automatically after each "
        "`generate_schema` against the "
        "default cache root; this command "
        "lets you trigger it manually.",
        help="Deduplicate cache entries by project set.",
    )
    _ = ctr  # subparser already registered

    # install-completion
    sp = sub.add_parser(
        "install-completion",
        formatter_class=fmt,
        description="Write a shell completion script so `trustify <TAB>` completes\n"
        "subcommands and flags. Defaults to bash; the file is dropped at\n"
        "a per-shell lazy-load path picked up automatically by\n"
        "bash-completion v2 (or by fish). For zsh, the file goes to a\n"
        "directory you must add to $fpath before `compinit`; the command\n"
        "prints the exact incantation after writing.\n\n"
        "Use --print to dump the script to stdout instead — handy when\n"
        "wiring it into a custom rc file or shipping it elsewhere.",
        help="Install shell completion for the trustify CLI.",
    )
    sp.set_defaults(cmd="install-completion")
    sp.add_argument(
        "--shell",
        choices=["bash", "zsh", "fish"],
        default=None,
        help="Target shell. Default: detected from $SHELL (fallback bash).",
    )
    sp.add_argument("--dest", metavar="PATH", help="Override the destination file path.")
    sp.add_argument(
        "--print",
        action="store_true",
        dest="do_print",
        help="Print the completion script to stdout instead of writing it.",
    )
    sp.add_argument("--force", "-f", action="store_true", help="Overwrite an existing destination file.")

    return p


def _die_if_missing(files: list[Path], parser: argparse.ArgumentParser, cmd: str) -> None:
    """Abort via `parser.error` if any path in `files` is not an existing
    regular file. All missing paths are reported in one message.
    """
    missing = [str(f) for f in files if not f.is_file()]
    if missing:
        parser.error(f"{cmd}: no such file: {', '.join(missing)}")


def _match_exclude(patterns: list[str], *, name: str, rel_path: str, is_dir: bool) -> bool:
    """Return True if `name` / `rel_path` is excluded by any glob in
    `patterns`.

    - trailing `/` → directory-only; strip it then match
    - pattern contains `/` → matched against `rel_path` (relative to the
      listed root, using forward slashes)
    - otherwise → matched against `name` (the basename)
    """
    for pat in patterns:
        if pat.endswith("/"):
            if not is_dir:
                continue
            pat = pat[:-1]
        if "/" in pat:
            if fnmatch.fnmatch(rel_path, pat):
                return True
        else:
            if fnmatch.fnmatch(name, pat):
                return True
    return False


def _walk_data_files(root: Path, patterns: list[str], follow_symlinks: bool = False, cmd: str = "") -> list[Path]:
    """Recursively yield `*.data` files under `root` in sorted order,
    with `patterns` pruning matched directories and skipping matched
    files.

    `follow_symlinks=False` (the default) drops both symlinked
    directories (via `os.walk(followlinks=False)`) AND symlinked
    files at the leaves. Skipped file symlinks emit a one-line
    warning to stderr so the user notices their tree contains
    out-of-scope links — silently dropping them was the audit-2.6
    footgun. `follow_symlinks=True` opts back in: dir symlinks are
    recursed and file symlinks are processed without warning.
    `cmd` only affects the warning prefix.
    """
    out: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=follow_symlinks):
        rel_dir = os.path.relpath(dirpath, root).replace(os.sep, "/")
        # Prune matched subdirectories in-place so os.walk skips them.
        # Sort so iteration order is deterministic across platforms.
        dirnames.sort()
        dirnames[:] = [
            d
            for d in dirnames
            if not _match_exclude(
                patterns,
                name=d,
                rel_path=d if rel_dir == "." else f"{rel_dir}/{d}",
                is_dir=True,
            )
        ]
        for fn in sorted(filenames):
            if not fn.endswith(".data"):
                continue
            rel = fn if rel_dir == "." else f"{rel_dir}/{fn}"
            if _match_exclude(patterns, name=fn, rel_path=rel, is_dir=False):
                continue
            full = Path(dirpath) / fn
            if not follow_symlinks and full.is_symlink():
                prefix = f"{cmd}: " if cmd else ""
                print(
                    f"{prefix}skipping symlinked file {full} (pass --follow-symlinks to include)",
                    file=sys.stderr,
                )
                continue
            out.append(full)
    return out


def _expand_data_paths(
    paths: list[Path],
    excludes: list[str],
    parser: argparse.ArgumentParser,
    cmd: str,
    follow_symlinks: bool = False,
) -> list[Path]:
    """Expand a mix of file and directory paths into `.data` files.

    `parser.error` is raised if any path doesn't exist on disk. Existing
    files are kept verbatim — explicit CLI files bypass `excludes` AND
    the symlink filter (the audit-2.6 concern is about UNINTENDED links
    surfaced by directory recursion, not about files the user explicitly
    named). Directories are recursed; `excludes` prune the walk;
    `follow_symlinks` controls whether symlinks discovered during
    recursion are followed (off by default — symlinked files emit a
    skip warning to stderr).
    """
    missing = [str(p) for p in paths if not p.exists()]
    if missing:
        parser.error(f"{cmd}: no such file or directory: {', '.join(missing)}")
    out: list[Path] = []
    for p in paths:
        if p.is_dir():
            out.extend(_walk_data_files(p, excludes, follow_symlinks=follow_symlinks, cmd=cmd))
        else:
            out.append(p)
    return out


def _format_size(num_bytes: int) -> str:
    """Render `num_bytes` as a human-readable size (B / KB / MB / GB)."""
    n = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024.0 or unit == "TB":
            if unit == "B":
                return f"{int(n)} {unit}"
            return f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{num_bytes} B"


def _format_created_at(iso_ts: str | None) -> str:
    """Trim an ISO-8601 timestamp to `YYYY-MM-DD HH:MM:SS`; `(unknown)` for None."""
    if not iso_ts:
        return "(unknown)"
    s = iso_ts.rstrip("Z").replace("T", " ")
    # Drop any sub-second fraction.
    dot = s.find(".")
    if dot != -1:
        s = s[:dot]
    return s


def _format_duration(seconds: float | None) -> str:
    """Render a duration as `12.3s` (or `1m23s` past one minute).
    `None` → `(unknown)`.
    """
    if seconds is None:
        return "(unknown)"
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    return f"{m}m{s:02d}s"


def _format_projects(projects: dict, max_width: int = 60) -> str:
    """Render the projects mapping as a comma-separated list of project
    names — paths are omitted (the table view would have to truncate
    them past the first one or two anyway). Empty dict → `(unknown)`.
    Use `cache list --json` to see the at-generation paths.
    """
    if not projects:
        return "(unknown)"
    s = ", ".join(projects.keys())
    if len(s) > max_width:
        s = s[: max_width - 1] + "..."
    return s


def _print_cache_table(entries: list) -> None:
    """Render the `cache list` table to stdout."""
    if not entries:
        print("No cache entries.")
        return
    rows = [
        (
            e.id[:8],
            _format_created_at(e.created_at),
            e.trustify_version or "(unknown)",
            _format_duration(e.generation_time_seconds),
            _format_size(e.size_bytes),
            _format_projects(e.projects),
        )
        for e in entries
    ]
    headers = ("ID", "CREATED", "VERSION", "TIME", "SIZE", "PROJECTS")
    cols = list(zip(headers, *rows, strict=False))
    widths = [max(len(str(cell)) for cell in col) for col in cols]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    for row in rows:
        print(fmt.format(*row))


def _path_args(args) -> list[Path]:
    """Return the positional path args of the subcommand (or `[]` if it
    doesn't take any). Used as anchors for the four-signal resolution —
    `.trustify.json` and `project.cfg` are looked up relative to the
    dataset being acted on, not relative to CWD.
    """
    paths: list[Path] = []
    for attr in ("data_file", "data_paths"):
        val = getattr(args, attr, None)
        if isinstance(val, str):
            paths.append(Path(val))
        elif isinstance(val, list):
            paths.extend(Path(p) for p in val)
    return paths


def _is_inside(child: Path, parent: Path) -> bool:
    """True when `child` resolves under `parent` (or equals it).

    Both sides are `.resolve()`-d first so symlink + relative-path noise
    can't fool the comparison. Used by `_resolve_project_args` to detect
    "anchor lives inside $TRUST_ROOT" and suppress stale baltik signals.
    """
    try:
        return child.resolve().is_relative_to(parent.resolve())
    except (OSError, ValueError):
        return False


def _discover_signals_for(anchor: Path) -> tuple[TrustifyConfig | None, Path | None]:
    """Look up the two upward-walking signals for a single anchor.
    Returns `(config, baltik_root)`. Caller decides what to do when
    multiple anchors disagree.
    """
    return (discover_config(anchor), detect_baltik_root(anchor))


def _resolve_project_args(args, parser: argparse.ArgumentParser) -> _CliResolution:
    """Apply the four-signal precedence and return kwargs for the
    programmatic API.

    Per field (projects, trust_root) the precedence is:

      CLI flag  >  .trustify.json  >  baltik auto-detection  >  env var

    For commands that take dataset paths, each path is used as a
    separate discovery anchor (the closest `.trustify.json` /
    `project.cfg` to that path wins). If two paths discover different
    configs (or different baltik roots), this errors out via
    `parser.error` (Q2 answer in the design discussion: be strict, not
    surprising). If the subcommand has no path args, CWD is the
    anchor.

    LSP / cache settings from a discovered config aren't returned here —
    only the schema-relevant fields. The full `TrustifyConfig` instance
    isn't exposed by this helper; consumers that need the LSP block
    call `discover_config` directly via the public `trustify` import.
    """
    anchors = [_start_anchor(p) for p in _path_args(args)] or [Path.cwd()]
    # The ConfigError raised by TrustifyConfig.from_file already carries
    # the offending file path in its message — pull the try outside the
    # loop so ruff PERF203 stays happy.
    try:
        signals = [_discover_signals_for(a) for a in anchors]
    except (ConfigError, ValueError) as e:
        parser.error(f"invalid .trustify.json: {e}")

    # Require all anchors to agree on both signals; mismatch means the user
    # is mixing workspaces and the schema choice would be ambiguous.
    distinct_configs = {(s[0].source_path if s[0] else None) for s in signals}
    distinct_baltiks = {s[1] for s in signals}
    if len(distinct_configs) > 1:
        seen = ", ".join(sorted(str(p) for p in distinct_configs if p)) or "(none)"
        parser.error(
            f"path arguments resolve to different .trustify.json files ({seen}); "
            "run trustify once per workspace, or pass --projects / --trust-root explicitly."
        )
    if len(distinct_baltiks) > 1:
        seen = ", ".join(sorted(str(p) for p in distinct_baltiks if p)) or "(none)"
        parser.error(
            f"path arguments resolve to different baltik roots ({seen}); "
            "run trustify once per baltik, or pass --projects explicitly."
        )

    config, baltik_root = signals[0]

    # Resolve trust_root first — projects resolution depends on knowing
    # whether the anchors live inside the TRUST source tree.
    trust_root = args.trust_root
    if trust_root is None and config and config.workspace.trust_root is not None:
        trust_root = config.workspace.trust_root
    trust_root = effective_trust_root(trust_root)

    # When every anchor sits inside the TRUST source tree, an implicit
    # baltik overlay (auto-detected project.cfg, or a `$project_directory`
    # left over from a previous `source <baltik>/env.sh`) is suspicious
    # but no longer dropped: overlaying a baltik whose sources live
    # inside TRUST — e.g. ICoCo — onto TRUST datasets is a legitimate
    # "do my baltik's XD tags break TRUST?" check. We keep the overlay
    # and warn (see below), leaving the call to the user. Explicit
    # signals (`--projects` flag, `.trustify.json workspace.projects`)
    # never warn — the user clearly asked for them.
    inside_trust_root = False
    if trust_root:
        tr = Path(trust_root).resolve()
        inside_trust_root = all(_is_inside(a, tr) for a in anchors)

    projects = args.projects
    explicit_projects = projects is not None
    if projects is None and config and config.workspace.projects is not None:
        projects = config.workspace.projects
        explicit_projects = True
    if projects is None and baltik_root is not None:
        projects = [str(baltik_root)]
    projects = effective_projects(projects)

    # Snapshot the pre-expansion list — the `projects` subcommand needs it
    # to distinguish primary (user-asked) entries from those pulled in by
    # transitive dependency expansion below.
    primaries = list(projects) if projects else None

    # An implicit baltik overlay active while every anchor sits inside
    # $TRUST_ROOT is the likely-stale-env case — warn, but proceed.
    if inside_trust_root and not explicit_projects and primaries:
        _warn_overlay_inside_trust_root(trust_root, baltik_root, primaries)

    # Auto-resolve transitive [dependencies] from each project's
    # project.cfg, unless the user opts out — `--no-auto-deps` wins,
    # else `.trustify.json` `workspace.auto_resolve_dependencies`,
    # else True.
    if args.no_auto_deps:
        auto_resolve = False
    elif config is not None:
        auto_resolve = config.workspace.auto_resolve_dependencies
    else:
        auto_resolve = True
    if auto_resolve:
        try:
            projects = expand_dependencies(projects)
        except ConfigError as e:
            parser.error(f"dependency resolution failed: {e}")

    _emit_resolution_summary(args, config, baltik_root, trust_root, primaries, projects)
    return _CliResolution(projects=projects, trust_root=trust_root, primaries=primaries)


def _warn_overlay_inside_trust_root(trust_root, baltik_root, primaries) -> None:
    """Warn that an implicit baltik overlay is active even though every
    dataset anchor sits inside $TRUST_ROOT.

    This used to be a hard suppression — the overlay was silently
    dropped. It is now only a warning: overlaying a baltik whose sources
    live inside the TRUST tree onto TRUST datasets is a legitimate way to
    test whether the baltik's XD tags break those datasets, so the choice
    is left to the user. An accidental stale `$project_directory` from a
    previous `source <baltik>/env.sh` stays visible rather than silently
    steering the result.
    """
    via = "auto-detected project.cfg" if baltik_root is not None else "$project_directory env"
    overlaid = ", ".join(primaries)
    print(
        f"trustify: WARNING: dataset anchor is inside $TRUST_ROOT ({trust_root}) but a baltik "
        f"overlay is active ({overlaid}, via {via}). Proceeding with the overlay. If this is a "
        "stale environment, pass --trust-root alone (or unset $project_directory) to check "
        "against plain TRUST.",
        file=sys.stderr,
    )


def _warn_modernize_inside_trust_root(trust_root, target, via) -> None:
    """Warn that modernize will rewrite a baltik's sources in place even
    though CWD sits inside $TRUST_ROOT.

    Like the schema resolver's inside-$TRUST_ROOT warning, this replaces
    a former hard suppression — modernizing a baltik whose sources live
    inside TRUST (e.g. ICoCo) is now allowed. The warning is sterner
    because modernize edits files in place: a stale `$project_directory`
    could otherwise rewrite the wrong tree.
    """
    print(
        f"trustify: WARNING: CWD is inside $TRUST_ROOT ({trust_root}) but modernize will rewrite "
        f"sources in {target} (via {via}). Proceeding. If this is a stale environment, pass "
        "--projects explicitly (or unset $project_directory) to choose the target, or run from "
        "inside the TRUST tree with no baltik env to modernize TRUST itself.",
        file=sys.stderr,
    )


def _emit_resolution_summary(args, config, baltik_root, trust_root, primaries, projects) -> None:
    """Print a one-shot stderr summary of how the workspace was resolved
    whenever any signal was implicit (auto-discovered config, baltik
    auto-detection, env var, transitive deps). Skipped entirely when
    every entry came from the CLI itself, when nothing resolved at all,
    or for the `projects` subcommand (which is itself a stdout dump of
    the same data).
    """
    if getattr(args, "cmd", None) == "projects":
        return

    rows: list[tuple[str, str, str]] = []  # (kind, path, source_label)

    if trust_root:
        if args.trust_root:
            src = "--trust-root flag"
        elif config and config.workspace.trust_root is not None:
            src = f"config {config.source_path}"
        else:
            src = "$TRUST_ROOT env"
        rows.append(("trust_root", str(trust_root), src))

    primaries_set = set(primaries or [])
    for p in primaries or []:
        if args.projects:
            src = "--projects flag"
        elif config and config.workspace.projects is not None:
            src = f"config {config.source_path}"
        elif baltik_root and str(baltik_root) == p:
            src = "auto-detected (project.cfg walk-up)"
        else:
            src = "$project_directory env"
        rows.append(("project", p, src))

    rows.extend(("dependency", p, "transitive [dependencies]") for p in projects or [] if p not in primaries_set)

    # Silent when nothing implicit kicked in — explicit flags are already
    # visible in the command line, no need to echo them back.
    explicit_sources = {"--trust-root flag", "--projects flag"}
    if not any(src not in explicit_sources for (_, _, src) in rows):
        return

    width = max(len(kind) for kind, _, _ in rows) if rows else 0
    print("trustify: resolved workspace:", file=sys.stderr)
    for kind, path, src in rows:
        print(f"  {kind:<{width}}  {path}  ({src})", file=sys.stderr)


def _scope_modernize(args, parser: argparse.ArgumentParser) -> tuple[list[str] | None, str | None, bool]:
    """Return `(projects, trust_root, warn_explicit)` for `trustify modernize`.

    Modernize rewrites `.cpp` / `.xd` sources in place, so it deliberately
    bypasses the four-signal precedence used by schema-consuming commands.
    The rules below pin the scope to the user's "current top-level project"
    by default and forbid implicit cross-checkout rewrites.

    Precedence (first match wins):

      1. ``--projects A [B ...]`` → exactly those projects, NO transitive
         ``[dependencies]`` expansion, NO `--trust-root` overlay. Returns
         ``warn_explicit=True`` so the handler emits a multi-project
         banner — schema-style scope makes sense for *reading* sources
         but is surprising for *writing* them.
      2. Auto-detected baltik — `project.cfg` walk-up from CWD or the
         ``$project_directory`` env var. Modernize that baltik only.
         ``--trust-root`` and ``$TRUST_ROOT`` are honoured for schema
         purposes elsewhere but ignored here.
      3. Only ``--trust-root`` / ``$TRUST_ROOT`` resolves (no baltik, no
         ``$project_directory``) → modernize the TRUST tree itself.

    Two failure modes are surfaced via ``parser.error``:

    - ``.trustify.json`` contributed ``workspace.projects`` or
      ``workspace.trust_root`` to the resolution. ``.trustify.json``
      captures a *user environment* (decoupling dataset folders from
      sourced env vars); driving in-place source rewrites from it is
      the wrong tool — the developer should be inside a baltik (or
      use ``--projects`` to opt in explicitly).
    - Nothing resolves at all.

    Anchor-inside-trust_root behaviour: when CWD lives under
    ``$TRUST_ROOT`` and an implicit baltik signal (auto-detected
    ``project.cfg`` or the ``$project_directory`` env) is active,
    modernize targets that baltik and prints a WARNING — it no longer
    suppresses the signal in favour of the whole TRUST tree. This mirrors
    the schema-resolving path (`_resolve_project_args`): a baltik whose
    sources live inside TRUST (e.g. ICoCo) is a valid modernize target.
    The warning is sterner here because modernize rewrites sources in
    place, so a stale env from a previous ``source <baltik>/env.sh``
    could rewrite the wrong tree; the choice is left to the user.
    """
    if args.projects:
        return list(args.projects), None, True

    try:
        config, baltik_root = _discover_signals_for(Path.cwd().resolve())
    except (ConfigError, ValueError) as e:
        parser.error(f"invalid .trustify.json: {e}")

    if config and (config.workspace.projects is not None or config.workspace.trust_root is not None):
        parser.error(
            f"modernize: workspace was resolved from {config.source_path}. "
            "This command rewrites source files and requires a developer "
            "environment: run from inside a baltik (project.cfg ancestry), "
            "set $project_directory, or pass --projects explicitly. Use "
            "--trust-root / $TRUST_ROOT alone to modernize the TRUST tree. "
            ".trustify.json is for user-environment management and must "
            "not drive in-place source rewrites."
        )

    trust_root = effective_trust_root(args.trust_root)
    inside_trust_root = bool(trust_root) and _is_inside(Path.cwd().resolve(), Path(trust_root).resolve())

    # Implicit baltik signals (auto-detected project.cfg first, then the
    # $project_directory env fallback) are honoured even when CWD sits
    # inside $TRUST_ROOT: a baltik whose sources live inside TRUST (e.g.
    # ICoCo) is a valid modernize target. We warn first so a stale env
    # from a previous `source <baltik>/env.sh` stays visible — but,
    # unlike a read-only schema build, modernize then rewrites that
    # baltik's sources in place.
    if baltik_root is not None:
        if inside_trust_root:
            _warn_modernize_inside_trust_root(trust_root, str(baltik_root), "auto-detected project.cfg")
        return [str(baltik_root)], None, False
    pd_env = os.environ.get("project_directory")  # noqa: SIM112
    if pd_env:
        if inside_trust_root:
            _warn_modernize_inside_trust_root(trust_root, pd_env, "$project_directory env")
        return [pd_env], None, False

    if trust_root:
        return [], trust_root, False

    parser.error(
        "modernize: nothing resolved to modernize. Pass --projects to list "
        "baltiks explicitly, run from inside a baltik (project.cfg ancestry), "
        "set $project_directory, or set $TRUST_ROOT to modernize the TRUST "
        "tree."
    )


def _start_anchor(path: Path) -> Path:
    """Map a CLI path arg to a directory the upward-walk helpers can
    take. A directory is returned as-is; a file becomes its parent.
    Non-existent paths are treated as files (parent), matching the
    "user is asking about THIS file" intent.
    """
    if path.is_dir():
        return path
    return path.parent


def main(argv: list[str] | None = None) -> int:
    p = _build_parser()
    # Shell completion hook: in $_ARGCOMPLETE mode this writes completions
    # to FD 8 and exits before parse_args. Soft import so stale TRUST envs
    # that haven't refreshed deps still work.
    try:
        import argcomplete

        argcomplete.autocomplete(p)
    except ImportError:
        pass
    args = p.parse_args(argv)

    if args.schema is not None and (args.projects or args.trust_root):
        p.error("--schema is mutually exclusive with --projects / --trust-root")

    if args.schema is not None and args.cmd not in _SCHEMA_CONSUMERS:
        p.error(f"--schema cannot be used with {args.cmd}")

    # Project resolution is expensive (filesystem walks, dependency
    # expansion, multi-path conflict detection) and only meaningful for
    # commands that actually consume the resolved list to build a
    # schema. Skip it for text-only / config-writing / cache-management
    # commands so spurious workspace mismatches don't block legitimate
    # invocations.
    resolved = _resolve_project_args(args, p) if args.cmd in _NEEDS_PROJECT_RESOLUTION else None
    kw = resolved.api_kwargs if resolved is not None else None

    try:
        if args.cmd == "generate_schema":
            d = api.generate_schema(out=args.out, **kw)
            print(d)
        elif args.cmd == "check":
            from trustify._color import colorize

            _die_if_missing([Path(args.data_file)], p, "check")
            r = api.check(args.data_file, schema=args.schema, no_skip=args.no_skip, **kw)
            line = f"{r.status}: {r.path}"
            if r.message:
                line += f" — {r.message}"
            if r.marker_justification:
                line += f" — {r.marker_justification}"
            print(colorize(line))
            if r.diff:
                # Round-trip mismatch: surface the captured diff so the
                # user can see what changed without re-running.
                print(r.diff)
            if r.obsolete_marker:
                warn = (
                    f"WARNING: obsolete '# TRUSTIFY NOT #' marker — "
                    f"{r.path} passes the check and the marker can be removed"
                )
                if r.marker_justification:
                    warn += f" — {r.marker_justification}"
                print(colorize(warn + "."), file=sys.stderr)
            # SKIPPED is a deliberate user opt-out (TRUSTIFY NOT marker)
            # — not an error. Only FAILED yields a non-zero exit.
            return 1 if r.status == "FAILED" else 0
        elif args.cmd in ("batch-check", "batch-check-projects"):
            if args.cmd == "batch-check-projects":
                data_path_objs = _resolve_check_project_paths(args, p, resolved)
                if data_path_objs:
                    n_dirs = len(data_path_objs)
                    print(f"batch-check-projects: checking {n_dirs} {'directory' if n_dirs == 1 else 'directories'}:")
                    for d_ in data_path_objs:
                        print(f"  {d_}")
            else:
                data_path_objs = [Path(p_) for p_ in args.data_paths]
            files = _expand_data_paths(
                data_path_objs,
                args.exclude,
                p,
                args.cmd,
                follow_symlinks=args.follow_symlinks,
            )
            if not files:
                print(f"{args.cmd}: no .data files found")
                return 0
            verbose = not args.quiet
            r = api.batch_check(files, schema=args.schema, no_skip=args.no_skip, verbose=verbose, jobs=args.jobs, **kw)
            obsolete = [c for c in r.results if c.obsolete_marker]
            # Collect every summary line in a buffer; emit to stdout/stderr
            # below and (optionally) write the same content to --summary-out.
            summary_lines: list[str] = []
            summary_lines.append(
                f"Summary: {r.passed}/{r.total} test(s) passed successfully "
                f"({r.skipped} skipped, {r.failed} failed)" + (f", {len(obsolete)} obsolete-marker" if obsolete else "")
            )
            # Always reprint the SKIPPED/FAILED list at the end — in verbose
            # mode they were printed during the run but get lost in the
            # progress output, and in quiet mode they weren't printed at all.
            # Sort non-passed entries so FAILED ones land last — that
            # keeps the actionable items in the terminal viewport when
            # the summary scrolls. Stable secondary sort by path keeps
            # the output deterministic across runs. SKIPPED first
            # (status alphabetical "F" < "S" → invert with a tuple
            # key), FAILED last.
            non_passed = sorted(
                (c for c in r.results if c.status != "PASSED"),
                key=lambda c: (c.status == "FAILED", c.path),
            )

            def _with_just(line: str, c) -> str:
                # Append the marker justification with ` — <text>` when
                # present; uniform across SKIPPED summary lines and the
                # obsolete-marker WARNING block.
                return f"{line} — {c.marker_justification}" if c.marker_justification else line

            for c in non_passed:
                summary_lines.append(_with_just(f"  [{c.status}] {c.path}: {c.message}", c))
                if c.diff:
                    # Captured round-trip diff — indent under the FAILED
                    # line so summary readers (and `--summary-out`
                    # consumers) can see the actual mismatch without
                    # re-running per-file. Empty lines stay empty (don't
                    # add stray indentation).
                    summary_lines.extend(f"    {ln}" if ln else "" for ln in c.diff.splitlines())
            # Snapshot the boundary BEFORE appending the obsolete-marker
            # block — stdout gets everything up to here, stderr gets the
            # rest. Snapshotting (rather than recomputing as 1 + N) keeps
            # this correct now that each FAILED entry may carry a
            # variable-length diff trailer.
            n_summary = len(summary_lines)
            if obsolete:
                summary_lines.append("")
                summary_lines.append(
                    f"WARNING: {len(obsolete)} dataset(s) carry an obsolete "
                    "'# TRUSTIFY NOT #' marker — they pass the check and "
                    "the marker can be removed:"
                )
                summary_lines.extend(_with_just(f"  [obsolete] {c.path}", c) for c in obsolete)

            # Stdout gets the summary + per-dataset list (with any diffs);
            # the obsolete-marker warning goes to stderr (preserving
            # existing test contracts). `colorize` is applied at print
            # time only — `summary_lines` itself stays ANSI-free so the
            # `--summary-out` file below is clean.
            from trustify._color import colorize

            for line in summary_lines[:n_summary]:
                print(colorize(line))
            for line in summary_lines[n_summary:]:
                print(colorize(line), file=sys.stderr)

            if args.summary_out:
                Path(args.summary_out).write_text("\n".join(summary_lines) + "\n")
            return 0 if r.failed == 0 else 1
        elif args.cmd == "format":
            f = Path(args.data_file)
            _die_if_missing([f], p, "format")
            formatted = api.format_dataset_file(f)
            if args.in_place:
                f.write_text(formatted)
            elif args.out is not None:
                Path(args.out).write_text(formatted)
            else:
                print(formatted, end="")
        elif args.cmd == "batch-format":
            files = _expand_data_paths(
                [Path(p_) for p_ in args.data_paths],
                args.exclude,
                p,
                "batch-format",
                follow_symlinks=args.follow_symlinks,
            )
            if not files:
                print("batch-format: no .data files found")
                return 0
            # Destructive operation is opt-in (audit 1.10). `--apply` or
            # TRUSTIFY_FORCE_FORMAT (any non-empty value, NO_COLOR
            # convention) writes in place; the default is a dry-run
            # that prints diffs and exits non-zero when any file would
            # change — CI-friendly, mirrors `modernize`.
            apply = args.apply or bool(os.environ.get("TRUSTIFY_FORCE_FORMAT"))
            result = api.batch_format(files, apply=apply)
            changed = result.changed_files
            if not apply:
                for fr in changed:
                    sys.stdout.write(fr.diff)
            verb = "formatted" if apply else "would change"
            if changed:
                print(f"summary: {len(changed)} files {verb}", file=sys.stderr)
            else:
                print(f"summary: 0 files {verb}", file=sys.stderr)
            if not apply and changed:
                return 1
            return 0
        elif args.cmd == "modernize":
            mp_projects, mp_trust_root, warn_explicit = _scope_modernize(args, p)
            report = api.modernize(projects=mp_projects, trust_root=mp_trust_root, apply=args.apply)
            changed = report.changed_files
            if not args.apply:
                for fr in changed:
                    sys.stdout.write(fr.diff)
            verb = "modernized" if args.apply else "would change"
            if changed:
                detail = (
                    f"{report.summary['brace']} brace flags, "
                    f"{report.summary['opt']} opt flags, "
                    f"{report.summary['splits']} splits, "
                    f"{report.summary['unsplittable']} unsplittable"
                )
                print(
                    f"summary: {len(changed)} files {verb} — {detail}",
                    file=sys.stderr,
                )
            else:
                print(f"summary: 0 files {verb}", file=sys.stderr)
            if warn_explicit:
                scope = "\n  ".join(mp_projects or [])
                n = len(mp_projects or [])
                print(
                    f"warning: --projects given; modernize scope is {n} explicitly-listed "
                    f"project(s):\n  {scope}\n"
                    "default scope is the current top-level project (project.cfg ancestry "
                    "or $project_directory). verify each listed path is intended — "
                    "modernize rewrites <project>/src/ in place.",
                    file=sys.stderr,
                )
            if not args.apply and changed:
                return 1
            return 0
        elif args.cmd == "init-config":
            out = api.init_config(
                directory=args.directory, projects=kw["projects"], trust_root=kw["trust_root"], overwrite=args.force
            )
            print(out)
        elif args.cmd == "generate_keywords":
            d = api.generate_keywords(out=args.out, schema=args.schema, **kw)
            print(d)
        elif args.cmd == "generate_markdown":
            try:
                d = api.generate_markdown(out=args.out, schema=args.schema, kw_ref_path=args.kw_ref_path, **kw)
            except api._PreconditionError as e:
                p.error(f"generate_markdown: {e}")
            print(d)
        elif args.cmd == "generate_pdf":
            # Resolve the workspace only when we must generate markdown;
            # --from-markdown reuses an existing tree and needs no schema.
            if args.from_markdown is not None:
                gp_kw = {"projects": None, "trust_root": None}
            else:
                gp_kw = _resolve_project_args(args, p).api_kwargs
            try:
                d = api.generate_pdf(
                    out=args.out,
                    from_markdown=args.from_markdown,
                    schema=args.schema,
                    doxygen=args.doxygen,
                    keep_build=args.keep_build,
                    **gp_kw,
                )
            except api._PreconditionError as e:
                p.error(f"generate_pdf: {e}")
            print(d)
        elif args.cmd == "projects":
            return _run_projects(args, p, resolved)
        elif args.cmd == "cache":
            return _run_cache(args, p)
        elif args.cmd == "install-completion":
            return _run_install_completion(args, p)
    except TrustifyInternalError:
        # Programmer-bug / schema-configuration error inside trustify
        # itself — re-raise unconditionally so the user gets the full
        # stack. The one-line CLI message would point at the trustify
        # call site, not the actual bug; the stack is the only
        # actionable signal. Audit 5.4.
        raise
    except Exception as e:
        # End-user mode: one-line stderr, no traceback. Under
        # TRUSTIFY_DEBUG (which already enables verbose parser/generator
        # logging — its documented purpose is "debug trustify itself")
        # re-raise so Python prints the full stack — contributors hitting
        # an internal crash should not have to fall back to `python -c`
        # to recover the traceback. Audit 5.1.
        if os.environ.get("TRUSTIFY_DEBUG"):
            raise
        print(f"trustify: {type(e).__name__}: {e}", file=sys.stderr)
        return 1
    return 0


def _detect_shell() -> str:
    """Pick a default shell for `install-completion` from `$SHELL`.

    Falls back to bash if the env var is missing or names something we
    don't support — bash is the original env_TRUST.sh case and the most
    common interactive shell.
    """
    name = Path(os.environ.get("SHELL", "")).name.lower()
    return name if name in ("bash", "zsh", "fish") else "bash"


def _default_completion_dest(shell: str) -> Path:
    """Per-shell lazy-load path for completion files.

    Bash and fish auto-discover these locations (bash-completion v2 for
    bash). Zsh has no equivalent lazy-load convention — we write a
    `_trustify` autoload function and rely on the user wiring the
    directory into `$fpath`; the post-install hint spells out how.
    """
    home = Path.home()
    if shell == "bash":
        base = Path(os.environ.get("XDG_DATA_HOME") or home / ".local" / "share")
        return base / "bash-completion" / "completions" / "trustify"
    if shell == "fish":
        base = Path(os.environ.get("XDG_CONFIG_HOME") or home / ".config")
        return base / "fish" / "completions" / "trustify.fish"
    if shell == "zsh":
        base = Path(os.environ.get("XDG_DATA_HOME") or home / ".local" / "share")
        return base / "zsh" / "site-functions" / "_trustify"
    raise ValueError(f"unknown shell {shell!r}")


def _post_install_hint(shell: str, dest: Path) -> str:
    if shell == "bash":
        return (
            "Open a new shell (or `exec bash`) to activate. Requires "
            "bash-completion v2+,\nwhich lazy-loads files from this path on first tab."
        )
    if shell == "fish":
        return "Open a new shell to activate; fish auto-discovers this path."
    if shell == "zsh":
        return (
            f"Add to your ~/.zshrc, BEFORE `compinit` runs:\n"
            f"  fpath=({dest.parent} $fpath)\n"
            f"Then open a new shell (or run `autoload -U compinit && compinit`)."
        )
    return ""


def _run_install_completion(args, parser: argparse.ArgumentParser) -> int:
    try:
        import argcomplete
    except ImportError:
        parser.error("install-completion: argcomplete is not installed (pip install argcomplete)")

    shell = args.shell or _detect_shell()
    code = argcomplete.shellcode(["trustify"], shell=shell)

    if args.do_print:
        sys.stdout.write(code)
        if not code.endswith("\n"):
            sys.stdout.write("\n")
        return 0

    dest = Path(args.dest).expanduser() if args.dest else _default_completion_dest(shell)
    if dest.exists() and not args.force:
        parser.error(f"install-completion: {dest} already exists; rerun with --force to overwrite.")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(code)
    print(f"Installed {shell} completion to {dest}")
    hint = _post_install_hint(shell, dest)
    if hint:
        print(hint)
    return 0


_PROJECT_ROLES = ("trust_root", "primary", "dependency")


def _resolved_project_roles(resolved: _CliResolution) -> list[tuple[str, str]]:
    """Tag each resolved entry with its role, in overlay order: trust_root
    first (if any), then `resolved.projects` (pre-expansion primaries vs
    transitively-added dependencies). Shared by `projects` and
    `batch-check-projects`."""
    primaries_set = {str(Path(p).resolve()) for p in (resolved.primaries or [])}
    entries: list[tuple[str, str]] = []
    if resolved.trust_root:
        entries.append((str(Path(resolved.trust_root).resolve()), "trust_root"))
    entries.extend((path, "primary" if path in primaries_set else "dependency") for path in (resolved.projects or []))
    return entries


def _parse_project_roles(raw: str, flag: str, parser: argparse.ArgumentParser) -> set[str]:
    """Parse a comma-separated ROLE list (for --only / --exclude), erroring
    on any unknown role."""
    roles = {r.strip() for r in raw.split(",") if r.strip()}
    invalid = roles - set(_PROJECT_ROLES)
    if invalid:
        parser.error(f"invalid role(s) in {flag}: {sorted(invalid)}; valid roles are {list(_PROJECT_ROLES)}")
    return roles


def _resolve_check_project_paths(args, parser: argparse.ArgumentParser, resolved: _CliResolution) -> list[Path]:
    """Map each resolved project to its `args.subdir` directory for
    `batch-check-projects`, filtered by `--only` role exactly like
    `trustify projects --only=...` (no `--only` → every resolved project).
    Projects lacking the subdir are skipped with a stderr note; an empty
    result short-circuits cleanly in the caller. Errors only when the
    workspace cannot be resolved at all."""
    if resolved is None or (not resolved.projects and not resolved.trust_root):
        parser.error(
            "batch-check-projects: no projects to resolve — run inside a baltik, or pass "
            "--projects / --trust-root (see `trustify projects`)."
        )
    entries = _resolved_project_roles(resolved)
    if args.only is not None:
        allowed = _parse_project_roles(args.only, "--only", parser)
        entries = [(pth, role) for pth, role in entries if role in allowed]
    out: list[Path] = []
    for pth, _role in entries:
        cand = Path(pth) / args.subdir
        if cand.is_dir():
            out.append(cand)
        else:
            print(f"batch-check-projects: skipping '{pth}' — no '{args.subdir}' directory.", file=sys.stderr)
    return out


def _run_projects(args, parser: argparse.ArgumentParser, resolved: _CliResolution) -> int:
    """Dump the resolved project list for shell scripting.

    Output is one path per line (overlay order: trust_root first if any,
    then deps in post-order, then primaries last). With `--json`, the
    output is a JSON array of `{path, role}` objects.

    `--only`/`--exclude` filter by role. Errors (exit 2 via parser.error)
    when nothing can be resolved at all — empty output after role
    filtering is NOT an error (it's just "no entries match that role").
    """
    if args.only is not None and args.exclude is not None:
        parser.error("--only and --exclude are mutually exclusive")

    if not resolved.projects and not resolved.trust_root:
        parser.error(
            "no projects to resolve: pass --projects, --trust-root, or set $TRUST_ROOT / "
            "$project_directory in the environment."
        )

    entries = _resolved_project_roles(resolved)
    if args.only is not None:
        allowed = _parse_project_roles(args.only, "--only", parser)
        entries = [(p, r) for p, r in entries if r in allowed]
    elif args.exclude is not None:
        excluded = _parse_project_roles(args.exclude, "--exclude", parser)
        entries = [(p, r) for p, r in entries if r not in excluded]

    if args.as_json:
        print(json.dumps([{"path": p, "role": r} for p, r in entries], indent=2))
    else:
        for path, _ in entries:
            print(path)
    return 0


def _run_cache(args, parser: argparse.ArgumentParser) -> int:
    """Dispatch the `cache <action>` sub-subcommand."""
    import json as _json
    import shutil

    from trustify import cache as _cache

    action = args.cache_action

    if action == "list":
        entries = _cache.iter_entries()
        if args.as_json:
            payload = [
                {
                    "id": e.id,
                    "path": str(e.path),
                    "created_at": e.created_at,
                    "trustify_version": e.trustify_version,
                    "projects": e.projects,
                    "size_bytes": e.size_bytes,
                    "generation_time_seconds": e.generation_time_seconds,
                    "generation_steps_seconds": e.generation_steps_seconds,
                }
                for e in entries
            ]
            print(_json.dumps(payload, indent=2, sort_keys=True))
        else:
            _print_cache_table(entries)
        return 0

    if action == "clean":
        try:
            paths = _resolve_clean_targets(args.ids)
        except (_cache.CacheEntryNotFound, _cache.AmbiguousCacheEntry) as e:
            parser.error(f"cache clean: {e}")
        if not paths:
            print("No cache entries to clean.")
            return 0
        if not args.ids and not args.force:
            total = sum(_cache._entry_size(p) for p in paths)
            try:
                reply = input(
                    f"Delete {len(paths)} cache entr{'y' if len(paths) == 1 else 'ies'} ({_format_size(total)})? [y/N] "
                )
            except EOFError:
                # Closed stdin (e.g. piped from /dev/null in a script) is
                # semantically a refusal — same outcome as the user
                # answering anything other than 'y'. Audit 5.2.
                reply = ""
            if reply.strip().lower() not in ("y", "yes"):
                print("Aborted.")
                return 1
        for path in paths:
            shutil.rmtree(path)
            print(f"Deleted {path.name}")
        return 0

    if action == "copy":
        try:
            src = _cache.resolve_id(args.id)
        except (_cache.CacheEntryNotFound, _cache.AmbiguousCacheEntry) as e:
            parser.error(f"cache copy: {e}")
        try:
            dest = _cache.copy_entry(src, Path(args.dest_dir))
        except FileExistsError as e:
            parser.error(f"cache copy: {e}")
        print(dest)
        return 0

    if action == "path":
        try:
            path = _cache.resolve_id(args.id)
        except (_cache.CacheEntryNotFound, _cache.AmbiguousCacheEntry) as e:
            parser.error(f"cache path: {e}")
        print(path)
        return 0

    if action == "trim":
        groups = _cache.trim_all_duplicates()
        if not groups:
            print("No duplicate cache entries.")
            return 0
        total_deleted = sum(len(v) for v in groups.values())
        for keeper, deleted in groups.items():
            print(f"Kept {keeper[:8]}; deleted {', '.join(d[:8] for d in deleted)}")
        print(
            f"Trimmed {total_deleted} duplicate "
            f"entr{'y' if total_deleted == 1 else 'ies'} across "
            f"{len(groups)} project set"
            f"{'' if len(groups) == 1 else 's'}."
        )
        return 0

    parser.error(f"cache: unknown action {action!r}")  # unreachable


def _resolve_clean_targets(ids: list[str]) -> list[Path]:
    """Return the set of cache directories to delete for `cache clean`.

    With no IDs: every entry under the cache root. With explicit IDs:
    one path per ID (each resolved via `resolve_id`).
    """
    from trustify import cache as _cache

    if not ids:
        return [e.path for e in _cache.iter_entries()]
    return [_cache.resolve_id(prefix) for prefix in ids]


if __name__ == "__main__":
    sys.exit(main())
