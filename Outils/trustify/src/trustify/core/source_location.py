"""Source-origin metadata for TRAD2 declarations.

Each XD declaration (block or attribute) records the file/line it was
parsed from. We store this as a `SourceLocation` triple
`(project, path, line)` where:

- `project` is a project name (the reserved literal `"trust"` for the
  TRUST source tree, or for a baltik the value read from
  `<baltik>/project.cfg` `[description].name`, falling back to the
  directory basename);
- `path` is the file path relative to that project's root;
- `line` is 1-indexed.

The absolute filesystem path is recovered at runtime by
`SourceLocation.resolve()` using a `projects` dict mapping project
names to absolute roots. The dict is provided by the consumer — at LSP
runtime it comes from the LSP's own configuration, so two clones of
the same project set produce byte-identical cache entries.
"""

import os
from dataclasses import dataclass
from pathlib import Path

TRUST_PROJECT_NAME = "trust"


class ProjectCollisionError(Exception):
    """Two projects resolve to the same name — pick distinct names or set
    a unique `[description].name` in one of their `project.cfg` files."""


@dataclass(frozen=True)
class SourceLocation:
    project: str
    path: str
    line: int

    def to_dict(self) -> dict:
        return {"project": self.project, "path": self.path, "line": self.line}

    @classmethod
    def from_dict(cls, d: dict) -> "SourceLocation":
        return cls(project=d["project"], path=d["path"], line=int(d["line"]))

    def resolve(self, projects: dict[str, str]) -> tuple[str, int] | None:
        """Return (absolute_path, line) using `projects[self.project]` as
        the root. Returns ``None`` when the project is not configured.
        """
        root = projects.get(self.project)
        if root is None:
            return None
        return str(Path(root) / self.path), self.line


def build_projects_dict(trust_root: str | None, projects: list[str]) -> dict[str, str]:
    """Return `{project_name: absolute_path}`.

    - `trust_root` (if provided) is assigned the reserved name `"trust"`.
    - Each entry in `projects` is named by
      `trustify.projects.project_name_for(path)` — reads
      `[description].name` from `<path>/project.cfg` when present, falls
      back to the directory basename otherwise.
    - Raises `ProjectCollisionError` on name collision (any two
      entries — including trust_root — that resolve to the same name).

    `trust_root` may be `None` or an empty string to indicate no trust
    project — both are treated the same way (no ``"trust"`` entry is added).
    This matches the LSP default setting where ``trust_root`` defaults to
    ``""``.
    """
    from trustify.projects import project_name_for

    out: dict[str, str] = {}

    def _add(name: str, path: str) -> None:
        if name in out and out[name] != path:
            raise ProjectCollisionError(
                f"project name collision on {name!r}: {out[name]} vs {path}. Rename or reorganise."
            )
        out[name] = path

    if trust_root:
        _add(TRUST_PROJECT_NAME, str(trust_root))
    for p in projects:
        _add(project_name_for(p), str(p))
    return out


def relativize(absolute_path: str, projects: dict[str, str]) -> tuple[str, str]:
    """Find which project owns `absolute_path` and return
    `(project_name, relative_path)`.

    Resolves longest matching project root first (so a nested project
    wins over its parent). Raises `ValueError` if no project owns the
    path.

    Both `absolute_path` and each project root are canonicalized via
    `os.path.realpath` before the prefix-match test (audit 4.3). Two
    clones with different symlink layouts (e.g. one passing the
    canonical trust root and one passing it via a symlink) must
    therefore produce identical output, so the resulting
    `source_locations.json` is byte-identical across clones and the
    cache hash stays stable.
    """
    path_norm = os.path.realpath(absolute_path).rstrip("/")
    candidates = sorted(
        ((name, os.path.realpath(root).rstrip("/")) for name, root in projects.items()),
        key=lambda kv: len(kv[1]),
        reverse=True,
    )
    for name, root_norm in candidates:
        if path_norm == root_norm or path_norm.startswith(root_norm + "/"):
            rel = path_norm[len(root_norm) + 1 :] if len(path_norm) > len(root_norm) else ""
            return name, rel
    raise ValueError(f"{path_norm!r} is not under any configured project root ({list(projects.values())})")
