"""Origin classifier — maps a `SourceLocation` (or its bare `project`
field) to a rendered origin label.

After the source-locations refactor, `_infoMain` is a
`(project, relative_path, line)` triple; the classifier just needs to
turn the `project` name into a human-readable label. The reserved
project name `"trust"` renders as the legacy uppercase label `"TRUST"`
so existing generated docs in baltiks keep working; every other
project renders as the name read from its ``project.cfg``
(``[description].name``), falling back to the directory basename.
"""

import os
import re

_TRUST_PROJECT_NAME = "trust"
_TRUST_LABEL = "TRUST"


class OriginClassifier:
    """Map a project name to its rendered origin label.

    Build via `from_projects(projects, trust_root)`; lookup via
    `classify_project(project_name)` or `classify_location(loc)`.
    """

    def __init__(self, project_to_label: dict[str, str]):
        # Preserve insertion order — `labels` is order-sensitive (it
        # drives the ordering of per-origin pages in the generated
        # markdown manual).
        self._project_to_label: dict[str, str] = dict(project_to_label)

    @property
    def labels(self) -> list[str]:
        """Distinct labels in insertion order."""
        seen: list[str] = []
        for label in self._project_to_label.values():
            if label not in seen:
                seen.append(label)
        return seen

    @property
    def is_multi(self) -> bool:
        return len(self.labels) > 1

    def classify_project(self, project_name: str) -> str | None:
        """Return the rendered label for `project_name`, or `None` when
        the project is not configured.
        """
        if not project_name:
            return None
        return self._project_to_label.get(project_name)

    def classify_location(self, loc) -> str | None:
        """Return the rendered label for `loc.project`, or `None` when
        unconfigured. `loc` is a `SourceLocation`-shaped object."""
        name = getattr(loc, "project", None)
        if not name:
            return None
        return self.classify_project(name)

    @classmethod
    def from_projects(cls, projects: list[str] | None, trust_root: str | None) -> "OriginClassifier":
        """Build the canonical layout: each `projects` entry first (so a
        baltik nested in `$TRUST_ROOT` appears before the TRUST entry in
        `labels` ordering), then the `trust_root` entry as `"TRUST"`.

        Each project's name is taken from
        `trustify.projects.project_name_for(path)` — reads
        `[description].name` from `<path>/project.cfg` when present,
        falls back to the directory basename. The label equals the
        name. `trust_root` is special-cased to the reserved project
        name `"trust"` and the legacy uppercase label `"TRUST"`.

        NOTE: this insertion order (projects → trust_root) is the
        opposite of `trustify.projects.resolve_projects` and
        `trustify.doc.extras.ExtrasIndex.from_projects`, which both
        put trust_root first. That is deliberate, not a bug — those
        two implement **overlay precedence** (TRUST is the base, later
        entries shadow earlier ones on basename collision), while this
        classifier drives **presentation order** in the generated
        manual (the user's own baltik pages appear before the TRUST
        reference pages). Keep the two conventions distinct.

        Collisions are not handled here —
        `trustify.core.source_location.build_projects_dict` upstream
        errors before we reach this point.
        """
        from trustify.projects import project_name_for

        mapping: dict[str, str] = {}
        for raw in projects or []:
            name = project_name_for(os.path.realpath(raw))
            mapping[name] = name
        if trust_root:
            mapping[_TRUST_PROJECT_NAME] = _TRUST_LABEL
        return cls(mapping)


def origin_slug(label: str) -> str:
    """Lowercase, page-id-safe slug derived from an origin label.

    Used to build per-origin page IDs (``KW_<slug>_<family>``) and
    filenames (``kw_<slug>_<family>.md``) in multi-origin mode.
    """
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")
