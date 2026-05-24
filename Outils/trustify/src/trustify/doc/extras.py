"""Per-project extras-dir scanner + lookup index.

Each origin (the trust_root path + each --projects entry) may ship a
``docs/trustify/extras/`` subdirectory. ``ExtrasIndex.from_projects``
walks every present dir, partitions files into a markdown map (keyed by
stem, ``.md`` only) and an image map (keyed by full filename,
everything else), and detects cross-origin duplicates eagerly.

During rendering, ``resolve_input(name)`` and ``resolve_image(filename)``
return file content / source path and mark the entry as used.
``api.generate_markdown`` then copies only ``used_images()`` to ``<out>/figures/``
and warns about anything left over (``unused_files()``).
"""

import os
from pathlib import Path

# Origin label for the trust_root tree. Must match the label
# OriginClassifier renders for the reserved "trust" project
# (trustify.doc.origin._TRUST_LABEL) so the two agree on how TRUST is
# named in user-facing output.
_TRUST_LABEL = "TRUST"


class ExtrasIndex:
    """Catalogue of files scanned from every origin's ``docs/trustify/extras/``
    directory, with per-entry usage tracking.

    Use ``from_projects`` to build; the explicit constructor is for tests.
    """

    def __init__(
        self,
        markdown_files: dict[str, Path],
        image_files: dict[str, Path],
        file_origins: dict[Path, str] | None = None,
    ) -> None:
        self._md = dict(markdown_files)
        self._img = dict(image_files)
        # Scanned-file path -> project-of-origin label. Empty when built
        # via the bare constructor (tests); populated by from_projects.
        self._file_origins: dict[Path, str] = dict(file_origins or {})
        self._used_md: set = set()
        self._used_img: set = set()

    def origin_of(self, path: Path) -> str | None:
        """Project-of-origin label for a scanned file — the same label
        OriginClassifier renders (``"TRUST"`` for the trust_root tree, the
        project name otherwise). ``None`` for files not tracked (e.g. an
        index built through the bare constructor)."""
        return self._file_origins.get(path)

    def used_images(self) -> list[Path]:
        """Source paths of images that have been resolved at least once."""
        return [self._img[name] for name in sorted(self._used_img)]

    def unused_files(self) -> list[Path]:
        """Paths of scanned files that were never resolved.

        Markdown files first (sorted by stem), then images (sorted by
        filename).
        """
        unused_md = [self._md[name] for name in sorted(self._md) if name not in self._used_md]
        unused_img = [self._img[name] for name in sorted(self._img) if name not in self._used_img]
        return unused_md + unused_img

    def resolve_input(self, name: str) -> str:
        """Return the markdown content of ``name.md`` and mark it used.

        Raises ``ValueError`` if no ``name.md`` was scanned from any
        origin.
        """
        path = self._md.get(name)
        if path is None:
            raise ValueError(
                f"trustify generate_markdown: \\input{{{{{name}}}}}: no extras/{name}.md found in any scanned project"
            )
        self._used_md.add(name)
        return path.read_text()

    def resolve_image(self, filename: str) -> Path:
        """Return the source path of ``filename`` and mark it used.

        Raises ``ValueError`` if no ``filename`` was scanned from any
        origin.
        """
        path = self._img.get(filename)
        if path is None:
            raise ValueError(
                f"trustify generate_markdown: \\includeimage{{{{{filename}}}}}: no extras/"
                f"{filename} found in any scanned project"
            )
        self._used_img.add(filename)
        return path

    @classmethod
    def from_projects(cls, projects: list[str] | None, trust_root: str | None) -> "ExtrasIndex":
        """Scan ``docs/trustify/extras/`` under each origin (trust_root
        first, then each --projects entry) and assemble the index.

        Raises ``ValueError`` if any scanned dir contains a subdirectory
        or if the same filename appears in two origins' extras dirs.
        Missing extras dirs are skipped silently.
        """
        from trustify.projects import project_name_for

        markdown_files: dict[str, Path] = {}
        image_files: dict[str, Path] = {}
        # (origin dir, project-of-origin label). The label matches what
        # OriginClassifier renders: "TRUST" for trust_root, the project
        # name (project.cfg [description].name, else dir basename) for
        # each --projects entry.
        origins: list[tuple[Path, str]] = []
        if trust_root:
            origins.append((Path(trust_root), _TRUST_LABEL))
        origins.extend((Path(raw), project_name_for(os.path.realpath(raw))) for raw in projects or [])
        # Track which origin contributed each name, to produce a
        # readable error on cross-origin collision.
        md_origin: dict[str, Path] = {}
        img_origin: dict[str, Path] = {}
        file_origins: dict[Path, str] = {}
        for origin, label in origins:
            extras = origin / "docs" / "trustify" / "extras"
            if not extras.is_dir():
                continue
            for entry in sorted(extras.iterdir()):
                if entry.is_dir():
                    raise ValueError(f"trustify generate_markdown: extras dir cannot contain subdirectories: {entry}")
                if not entry.is_file():
                    # Symlinks to non-files, sockets, etc — skip silently.
                    continue
                if entry.suffix == ".md":
                    key = entry.stem
                    if key in markdown_files:
                        raise ValueError(
                            f"trustify generate_markdown: duplicate extras file "
                            f"'{entry.name}' in "
                            f"'{md_origin[key]}' and '{extras}'"
                        )
                    markdown_files[key] = entry
                    md_origin[key] = extras
                else:
                    key = entry.name
                    if key in image_files:
                        raise ValueError(
                            f"trustify generate_markdown: duplicate extras file "
                            f"'{entry.name}' in "
                            f"'{img_origin[key]}' and '{extras}'"
                        )
                    image_files[key] = entry
                    img_origin[key] = extras
                file_origins[entry] = label
        return cls(markdown_files, image_files, file_origins)
