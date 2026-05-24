"""Project resolution + `.trustify.json` discovery.

A "project", in trustify's model, is a directory whose `src/` subtree
contains C++ sources (`*.cpp`) and/or `*.xd` files with `// XD` comments.
trustify recursively scans `<project>/src` and aggregates the XD tags
(missing `src/` is a hard error).

This module is pure logic: it does not import from `trustify.core` and
does not touch the schema. See `Outils/trustify/docs/project-resolution.md`
for the user-facing description of the four-signal precedence
(CLI flag > .trustify.json > baltik auto-detection > env var) and the
CLI-vs-API split.
"""

import configparser
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

CONFIG_FILENAME = ".trustify.json"
BALTIK_MARKER_FILENAME = "project.cfg"


class ConfigError(Exception):
    """Raised when a .trustify.json file is malformed or contains an invalid value."""


@dataclass(frozen=True)
class Project:
    """A trustify input project. `path` is always absolute and `~`/$VAR-expanded.

    `is_trust_root` distinguishes the (at most one) Project that came
    from the `trust_root` channel: TRUST itself has no `project.cfg`
    at its root, so the baltik-shape validation in `validate_projects`
    skips it. Every other entry IS expected to be a baltik (i.e. to
    contain a `project.cfg`).
    """

    path: Path
    is_trust_root: bool = False


@dataclass(frozen=True)
class CacheConfig:
    """Cache-related knobs from a `.trustify.json`."""

    auto_trim: bool = True


@dataclass(frozen=True)
class WorkspaceConfig:
    """Project-resolution knobs from a `.trustify.json` `workspace`
    section. Groups the three options that together describe *which*
    sources trustify should scan: the dependency baltiks (`projects`),
    the TRUST source tree (`trust_root`), and whether transitive
    `[dependencies]` from each project's `project.cfg` should be
    auto-resolved (`auto_resolve_dependencies`).

    `auto_resolve_dependencies` (default True): when False,
    `expand_dependencies` becomes a no-op and trustify uses only the
    projects explicitly listed in --projects / .trustify.json /
    auto-detection — without transitively pulling in their
    `[dependencies]`. Equivalent to passing the `--no-auto-deps` CLI
    flag.
    """

    projects: list[str] | None = None
    trust_root: str | None = None
    auto_resolve_dependencies: bool = True


@dataclass(frozen=True)
class LspConfig:
    """LSP-related knobs from a `.trustify.json`. Forwarded as-is to the
    trustify-lsp server — this module doesn't consume them, just parses
    them so the LSP doesn't need to duplicate the discovery logic.
    """

    schema_dir: str | None = None
    enum_dedup_threshold: int | None = None


@dataclass(frozen=True)
class TrustifyConfig:
    """Parsed `.trustify.json`.

    Built via `TrustifyConfig.from_file(path)` or, more usually, found
    automatically via the top-level `discover_config` helper. The
    canonical four-signal precedence (see docs) layers this on top of
    auto-detection and env vars; CLI flags still win when present.

    `source_path` is the absolute path of the file the config was loaded
    from, useful for error messages and debugging.
    """

    source_path: Path
    workspace: WorkspaceConfig = field(default_factory=WorkspaceConfig)
    cache: CacheConfig = field(default_factory=CacheConfig)
    lsp: LspConfig = field(default_factory=LspConfig)
    # Original parsed payload — exposed so the LSP (or any future
    # consumer) can read keys we don't model here without going through
    # another loader.
    raw: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_file(cls, path: Path) -> "TrustifyConfig":
        """Parse + validate a `.trustify.json` file. Raises `ConfigError`
        when any recognized key carries the wrong type, or when the file
        sits next to a `project.cfg` AND declares
        `workspace.projects` / `workspace.trust_root` (see
        `_check_no_project_overlap_with_baltik` for why that combination
        is rejected). Unparseable JSON raises the underlying
        `json.JSONDecodeError` (caller may want to wrap that into a
        friendlier message).
        """
        resolved = Path(path).resolve()
        data = _load_and_validate(resolved)
        _check_no_project_overlap_with_baltik(resolved, data)
        ws_d = data.get("workspace") or {}
        cache_d = data.get("cache") or {}
        lsp_d = data.get("lsp") or {}
        return cls(
            source_path=resolved,
            workspace=WorkspaceConfig(
                projects=list(ws_d["projects"]) if "projects" in ws_d else None,
                trust_root=ws_d.get("trust_root"),
                auto_resolve_dependencies=ws_d.get("auto_resolve_dependencies", True),
            ),
            cache=CacheConfig(auto_trim=cache_d.get("auto_trim", True)),
            lsp=LspConfig(
                schema_dir=lsp_d.get("schema_dir"),
                enum_dedup_threshold=lsp_d.get("enum_dedup_threshold"),
            ),
            raw=data,
        )

    def to_api_kwargs(self) -> dict:
        """Return `(projects, trust_root)` kwargs suitable for
        `api.generate_schema` / `api.load_dataset` / etc.

        Missing fields are returned as None so the caller can decide
        whether to fill them in (typically by chaining with env-var
        fallbacks via `effective_projects` / `effective_trust_root`).
        """
        return {"projects": self.workspace.projects, "trust_root": self.workspace.trust_root}


def _expand(s: str) -> Path:
    """Expand `~`, `$VAR`, and `${VAR}` references, then return an absolute Path."""
    return Path(os.path.expandvars(os.path.expanduser(s))).resolve()


# Sentinel name for `ConfigParser(default_section=...)`. Python's
# `configparser` injects every option from `[DEFAULT]` into every other
# section, which baltik_configure does NOT do (it treats `[DEFAULT]` as
# a regular, unused section). Redirecting the "default section" slot to
# an obviously-unused name turns `[DEFAULT]` into a normal section and
# kills the inheritance magic. See the `baltik_configure parity` appendix
# in docs/project-resolution.md.
_UNUSED_DEFAULT_SECTION = "_trustify_no_default_section_"


def _read_baltik_config(cfg_path: Path) -> configparser.ConfigParser | None:
    """Parse a baltik `project.cfg` with the same lenience as
    `baltik_configure`. Returns `None` if the file is malformed.

    The configparser instance is preconfigured to match baltik:

    - `inline_comment_prefixes=("#",)` mirrors `sed 's/#.*//'` — a
      trailing `# comment` after a value is dropped, on any line.
    - `default_section=_UNUSED_DEFAULT_SECTION` makes `[DEFAULT]` a
      regular (ignored) section instead of a magic per-section default
      provider.
    - `optionxform = str` keeps option keys case-sensitive — baltik
      matches keys with case-sensitive `sed`, so `Name = Foo` does NOT
      match a lookup of `name`. Symmetric with mixed-case dependency
      keys like `TrioCFD : ...`.
    """
    cp = configparser.ConfigParser(
        inline_comment_prefixes=("#",),
        default_section=_UNUSED_DEFAULT_SECTION,
    )
    cp.optionxform = str
    try:
        cp.read(cfg_path)
    except configparser.Error as e:
        # Don't silently degrade — a typo in [description].name used to
        # fall back to the directory basename with no diagnostic, so the
        # baltik silently shipped a wrong project name into the cache
        # hash. Warn at least once per call site (audit 5.6).
        from trustify.core.misc_utilities import logger

        logger.warning("malformed baltik config %s: %s — falling back to defaults", cfg_path, e)
        return None
    return cp


def project_name_for(path) -> str:
    """Return the canonical project name for the directory at ``path``.

    Reads ``[description].name`` from ``<path>/project.cfg`` (the baltik
    config file that ``bin/trust baltik`` and ``baltik_configure`` lay
    down). Falls back to the directory basename when no project.cfg
    exists, the file is malformed, or the ``name`` field is missing or
    empty — that keeps the helper safe to call on every project
    regardless of whether it is a baltik. (TRUST itself has no
    project.cfg at its root; the reserved ``"trust"`` name is assigned
    separately at the caller.)
    """
    p = Path(path)
    cfg = p / BALTIK_MARKER_FILENAME
    if cfg.is_file():
        cp = _read_baltik_config(cfg)
        if cp is not None:
            name = cp.get("description", "name", fallback="").strip()
            if name:
                return name
    return p.name or "project"


def effective_trust_root(trust_root: str | None) -> str | None:
    """Return `trust_root` if truthy, otherwise the value of `$TRUST_ROOT`
    in the environment, otherwise None.

    Lowest-priority `trust_root` source — the CLI layers config-file and
    baltik auto-detection above this. Programmatic API callers use this
    helper directly to get env-var fallback only (no implicit
    discovery).
    """
    return trust_root or os.environ.get("TRUST_ROOT")


def effective_projects(projects: list[str] | None) -> list[str] | None:
    """Return `projects` when the caller supplied a value (including the
    empty list), otherwise fall back to `[$project_directory]` if that
    env var is set, otherwise `None`.

    Lowest-priority `projects` source — the CLI layers config-file and
    baltik auto-detection (via `detect_baltik_root`) above this.
    `$project_directory` and `detect_baltik_root` are intentionally
    redundant: detection is the more reliable path (no env hygiene
    required) but the env var still works if a user `cd`'d out of the
    baltik before invoking trustify.

    `None` vs `[]` semantics are deliberate:
      - `None` means "caller didn't specify anything" — env-var
        fallback applies.
      - `[]` means "caller explicitly resolved to no baltik overlay"
        — env-var fallback is suppressed. The CLI uses this to honour
        the anchor-inside-trust_root rule even when it delegates to
        `api.generate_schema(projects=...)`, which would otherwise
        re-read `$project_directory` and overlay the baltik back in.
    """
    if projects is not None:
        return projects
    pd = os.environ.get("project_directory")  # noqa: SIM112
    return [pd] if pd else None


def resolve_projects(projects: list[str] | None, trust_root: str | None) -> list[Project]:
    """Return the final ordered project list.

    `trust_root` (if provided) is prepended to `projects` — this matches the
    overlay precedence used by the legacy `install.sh` flow, where TRUST
    sources are scanned first and baltik sources override them.

    Raises `ValueError` if both `projects` and `trust_root` are empty/None.
    """
    out: list[Project] = []
    if trust_root:
        out.append(Project(path=_expand(trust_root), is_trust_root=True))
    out.extend(Project(path=_expand(p), is_trust_root=False) for p in projects or [])
    if not out:
        raise ValueError(
            "No projects to resolve: pass --projects, --trust-root, or set $TRUST_ROOT in the environment."
        )
    return out


def validate_projects(resolved: list[Project]) -> None:
    """Validate that every non-trust_root project is a proper baltik.

    A baltik is required to carry:

    - a top-level `project.cfg` declaring its name in `[description].name`
      (`bin/trust baltik` and `baltik_configure` lay this down on
      creation), and
    - a top-level `src/` directory — that's where trustify's source
      scanner looks for `*.cpp` / `*.xd` files. An empty `src/` is fine
      (the baltik simply contributes no XD tags and the schema gets
      everything from trust_root + other overlays); a *missing* one
      would surface mid-pipeline as a bare scanner exception, so we
      reject it up-front with a clearer message.

    Treating a bare directory as a project (silently falling back to its
    basename for the project name) used to work but masked
    misconfigurations — most notably `.trustify.json` `projects`
    entries pointing at the wrong folder, or developers accidentally
    passing a build directory as a `--projects` value.

    Raises `ConfigError` with an actionable message at the first
    offending project. The trust_root entry is always skipped — TRUST
    itself has its own layout.
    """
    for p in resolved:
        if p.is_trust_root:
            continue
        marker = p.path / BALTIK_MARKER_FILENAME
        if not marker.is_file():
            raise ConfigError(
                f"trustify project {p.path} is not a valid baltik: missing "
                f"`{BALTIK_MARKER_FILENAME}` at its root. Every directory passed "
                "as --projects (or listed under `projects` in .trustify.json) must "
                "be a baltik with `[description].name` set in its project.cfg — "
                "run `baltik_configure` from inside the baltik to create one. "
                "If you intended to point at the TRUST source tree, pass it via "
                "--trust-root (or $TRUST_ROOT) instead of --projects."
            )
        src_dir = p.path / "src"
        if not src_dir.is_dir():
            raise ConfigError(
                f"trustify project {p.path}: missing `src/` directory. Every "
                "baltik must have a top-level `src/` subtree — trustify's "
                "source scanner looks there for `*.cpp` / `*.xd` files. An "
                "empty `src/` is fine (the baltik then contributes no XD "
                "tags), but a missing one is a misconfiguration. Run "
                "`baltik_configure` from inside the baltik to scaffold the "
                "layout, or `mkdir src/` if the baltik genuinely carries no "
                "sources."
            )


def _start_dir_for(start: Path | None) -> Path:
    """Pick the directory to anchor an upward walk on.

    `None` → CWD. A directory → that directory. A file → its parent.
    Always returned as a resolved absolute path.
    """
    if start is None:
        return Path.cwd().resolve()
    start = Path(start)
    if start.is_dir():
        return start.resolve()
    if start.is_file():
        return start.parent.resolve()
    # Path that doesn't exist on disk — treat its `.parent` as the
    # anchor, the caller may be probing for a yet-to-be-created file.
    return start.parent.resolve()


def _discover_config_path(start_dir: Path) -> Path | None:
    """Walk from `start_dir` up to filesystem root, returning the first
    `.trustify.json` found, or None when no such file exists in any
    ancestor.

    Private — public callers should use `discover_config` (returns the
    parsed config, not the path).
    """
    d = Path(start_dir).resolve()
    while True:
        candidate = d / CONFIG_FILENAME
        if candidate.is_file():
            return candidate
        if d.parent == d:
            return None
        d = d.parent


def discover_config(start: Path | str | None = None) -> TrustifyConfig | None:
    """Walk up from `start` (or CWD if None) looking for a
    `.trustify.json` and return its parsed `TrustifyConfig`. Returns
    None when no config file is found in any ancestor.

    `start` may be a directory or a file path; in the file case the walk
    starts from its parent directory. Useful when resolving the workspace
    for a specific dataset file rather than the current shell.

    Programmatic API users who want CLI-style auto-discovery should call
    this helper explicitly and pass the result through `to_api_kwargs()`
    into `generate_schema` / `load_dataset` / etc. The API itself never
    walks the filesystem for configs.

    Raises `ConfigError` on malformed JSON content (validated keys);
    the underlying `json.JSONDecodeError` leaks through for unparseable
    files.
    """
    start_dir = _start_dir_for(Path(start) if start is not None else None)
    path = _discover_config_path(start_dir)
    return TrustifyConfig.from_file(path) if path else None


def detect_baltik_root(start: Path | str | None = None) -> Path | None:
    """Walk up from `start` (or CWD if None) looking for a directory
    containing `project.cfg`. Returns the absolute path of that
    directory, or None when no such directory is found.

    Identifies the "current baltik" without depending on the build
    environment having set `$project_directory`. Used by the CLI as the
    third-priority `projects` source (after explicit CLI flags and a
    `.trustify.json`); see the docs for the full precedence ladder.
    """
    start_dir = _start_dir_for(Path(start) if start is not None else None)
    d = start_dir.resolve()
    while True:
        if (d / BALTIK_MARKER_FILENAME).is_file():
            return d
        if d.parent == d:
            return None
        d = d.parent


def _load_and_validate(path: Path) -> dict:
    """Parse a `.trustify.json` file and validate every recognized key's
    type. Raises `ConfigError` with the file path + key name when
    something is off.

    Validated keys (top-level):
      - `workspace`: dict, with:
          - `workspace.projects`: list[str]
          - `workspace.trust_root`: str
          - `workspace.auto_resolve_dependencies`: bool
      - `cache`: dict, with `cache.auto_trim`: bool
      - `lsp`: dict, with `lsp.schema_dir`: str, `lsp.enum_dedup_threshold`: int

    Unrecognized keys are silently kept in the returned dict (forward
    compatibility with future fields, particularly LSP-side ones).
    """
    with open(path) as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: top-level value must be a JSON object, got {type(data).__name__}")
    if "workspace" in data:
        ws = data["workspace"]
        if not isinstance(ws, dict):
            raise ConfigError(f"{path}: 'workspace' must be a JSON object")
        if "projects" in ws:
            if not isinstance(ws["projects"], list):
                raise ConfigError(f"{path}: 'workspace.projects' must be a list of strings")
            if not all(isinstance(x, str) for x in ws["projects"]):
                raise ConfigError(f"{path}: 'workspace.projects' entries must be strings")
        if "trust_root" in ws and not isinstance(ws["trust_root"], str):
            raise ConfigError(f"{path}: 'workspace.trust_root' must be a string")
        if "auto_resolve_dependencies" in ws and not isinstance(ws["auto_resolve_dependencies"], bool):
            raise ConfigError(f"{path}: 'workspace.auto_resolve_dependencies' must be a boolean")
    if "cache" in data:
        if not isinstance(data["cache"], dict):
            raise ConfigError(f"{path}: 'cache' must be a JSON object")
        if "auto_trim" in data["cache"] and not isinstance(data["cache"]["auto_trim"], bool):
            raise ConfigError(f"{path}: 'cache.auto_trim' must be a boolean")
    if "lsp" in data:
        if not isinstance(data["lsp"], dict):
            raise ConfigError(f"{path}: 'lsp' must be a JSON object")
        if "schema_dir" in data["lsp"] and not isinstance(data["lsp"]["schema_dir"], str):
            raise ConfigError(f"{path}: 'lsp.schema_dir' must be a string")
        if "enum_dedup_threshold" in data["lsp"] and not isinstance(data["lsp"]["enum_dedup_threshold"], int):
            raise ConfigError(f"{path}: 'lsp.enum_dedup_threshold' must be an integer")
    return data


def _check_no_project_overlap_with_baltik(config_path: Path, data: dict) -> None:
    """Reject a `.trustify.json` that declares `workspace.projects` or
    `workspace.trust_root` while sitting next to a `project.cfg`.

    Inside a baltik (where the marker file lives at the baltik root),
    the project IS auto-detected from `project.cfg`. A coexisting
    `.trustify.json` should carry only LSP / cache settings; if it
    redeclares the project list or `trust_root`, the intent is
    ambiguous and was previously silently shadowing the auto-detected
    baltik. Surface the misconfiguration loudly so the author can fix
    it.

    `config_path` must be an absolute resolved path.
    """
    if not (config_path.parent / BALTIK_MARKER_FILENAME).is_file():
        return
    ws = data.get("workspace") or {}
    declared = [k for k in ("projects", "trust_root") if k in ws]
    if not declared:
        return
    keys = " and ".join(f"`workspace.{k}`" for k in declared)
    raise ConfigError(
        f"{config_path}: declares {keys} but sits next to a `{BALTIK_MARKER_FILENAME}`. "
        "When trustify is run from inside a baltik, the project is auto-detected from "
        "project.cfg — a co-located .trustify.json must carry only LSP / cache settings "
        "(no `workspace.projects`, no `workspace.trust_root`). Either remove the "
        "redundant keys from this file, or move the .trustify.json outside the baltik "
        "tree if you want it to override the auto-detected project."
    )


def load_config(path: Path) -> dict:
    """Back-compat shim returning the raw dict, for callers that don't
    want a `TrustifyConfig` instance. New code should use
    `TrustifyConfig.from_file(path)` directly.
    """
    return _load_and_validate(path)


def _parse_dependencies(project_cfg_path: Path) -> dict[str, str]:
    """Return the `[dependencies]` section of a `project.cfg` as a dict
    of ``{dep_name: raw_path_string}``. Raw path strings are returned
    verbatim — the caller is responsible for substituting shell
    constructs (``\\`pwd\\```, env vars) and resolving against the
    baltik root via ``_resolve_dep_path``.

    Returns an empty dict when the `[dependencies]` section is missing,
    or when the project.cfg cannot be read.
    """
    if not project_cfg_path.is_file():
        return {}
    cp = _read_baltik_config(project_cfg_path)
    if cp is None or not cp.has_section("dependencies"):
        return {}
    return dict(cp.items("dependencies"))


def _resolve_dep_path(raw_path: str, baltik_root: Path) -> Path:
    """Resolve a dependency path declaration as written in a baltik's
    `project.cfg`. baltik_configure evaluates these as shell strings;
    trustify supports the subset that occurs in practice:

    - ``\\`pwd\\``` — substituted with `baltik_root` (the directory
      containing the declaring `project.cfg`). baltik_configure runs
      with cwd == baltik_root, so this matches the runtime behaviour.
    - ``$VAR``, ``${VAR}``, ``~`` — standard env-var / home expansion.
    - relative paths — resolved against `baltik_root`.
    - absolute paths — used as-is.

    Any other shell construct (`` `cmd` `` for cmd != pwd, ``$(...)``,
    pipes, etc.) is left intact and will probably resolve to a
    nonexistent path — the caller should surface that as a clean
    ConfigError rather than crashing.

    Matching surrounding `"`/`'` quotes are stripped — `baltik_configure`'s
    `eval "echo $entry"` strips them naturally, so `dep = "/some path"`
    must resolve to `/some path` and not to a literal `"/some path"`.
    """
    expanded = raw_path.strip()
    if (
        len(expanded) >= 2
        and expanded[0] == expanded[-1]
        and expanded[0] in ('"', "'")
        # Refuse to strip when the inner content itself contains the
        # same quote — `"foo"bar"` (three `"`) would yield `foo"bar`,
        # silently mangling already-malformed input into something
        # different and still malformed. Audit 6.9.
        and expanded[0] not in expanded[1:-1]
    ):
        expanded = expanded[1:-1]
    expanded = expanded.replace("`pwd`", str(baltik_root))
    expanded = os.path.expandvars(os.path.expanduser(expanded))
    p = Path(expanded)
    if not p.is_absolute():
        p = baltik_root / p
    return p.resolve()


def expand_dependencies(projects: list[str] | None) -> list[str] | None:
    """Recursively expand each project's `[dependencies]`. Returns the
    full overlay-correct list of project paths (deduplicated, in
    post-order: every dependency appears before the baltik that
    declares it).

    For each path in `projects`:
      - Read its `project.cfg` `[dependencies]` section.
      - Resolve each `dep_name : raw_path` pair via `_resolve_dep_path`.
      - Recurse into each resolved dep.
      - Append the declaring baltik AFTER its deps (post-order — so the
        baltik overrides its dependencies in `scanSourceFiles`).

    Validates that `dep_name` matches the dep's own
    `[description].name` (mirrors `baltik_configure`'s
    `invalid_dependency_name_error`), and that each declared dep
    points at an existing directory. Also rejects two distinct paths
    that both declare the same `[description].name` — mirrors
    `baltik_configure`'s `too_many_paths_error` (see
    `baltik_dependencies_management::check_dependencies`). Raises
    `ConfigError` with an actionable message in either case. Trust_root
    paths (which never appear in this list) are handled by the caller,
    not here.

    Empty / None input returns unchanged — keeps callers from having
    to special-case the no-projects case.
    """
    if not projects:
        return projects
    seen: set[Path] = set()
    seen_by_name: dict[str, Path] = {}
    result: list[str] = []

    def _visit(raw_path_str: str, declared_by: Path | None = None, declared_as: str | None = None) -> None:
        path = Path(raw_path_str).resolve()
        if path in seen:
            return
        if not path.is_dir():
            if declared_by is None:
                # Top-level project — error coming from --projects /
                # .trustify.json[projects] / auto-detect / env-var.
                raise ConfigError(f"project list: {path} does not exist or is not a directory.")
            raise ConfigError(
                f"{declared_by / BALTIK_MARKER_FILENAME}: dependency {declared_as!r} points "
                f"at {path}, which does not exist or is not a directory."
            )
        cfg_path = path / BALTIK_MARKER_FILENAME
        actual_name = project_name_for(path)
        # Cross-check the project name when this entry was declared as a
        # dependency — the parent's `[dependencies]` key must match the
        # dep's own `[description].name`, exactly like baltik_configure
        # enforces.
        if declared_as is not None and actual_name != declared_as:
            raise ConfigError(
                f"{declared_by / BALTIK_MARKER_FILENAME}: dependency declared as "
                f"`{declared_as}` but {cfg_path} sets `[description].name = {actual_name}`. "
                "These must match (same rule as baltik_configure)."
            )
        # Two different baltik directories must not both claim the same
        # `[description].name`. baltik_configure's `check_dependencies`
        # bails out on `too_many_paths_error` for this case; we mirror it.
        previous = seen_by_name.get(actual_name)
        if previous is not None and previous != path:
            raise ConfigError(
                f"two distinct baltiks both declare `[description].name = {actual_name}`: "
                f"{previous} and {path}. Same rule as baltik_configure's "
                "`too_many_paths_error` — pick one or rename."
            )
        seen.add(path)
        seen_by_name[actual_name] = path
        for dep_name, raw in _parse_dependencies(cfg_path).items():
            dep_path = _resolve_dep_path(raw, path)
            _visit(str(dep_path), declared_by=path, declared_as=dep_name)
        result.append(str(path))

    for p in projects:
        _visit(p)
    return result
