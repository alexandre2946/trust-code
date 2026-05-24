"""Test fixtures for building mock baltiks on disk.

This module is intentionally tiny and stays under ``tests/`` — it is
NOT public API. The aim is to remove the boilerplate that was
duplicating across ~6 test modules:

    baltik = parent / name
    baltik.mkdir(parents=True)
    (baltik / "project.cfg").write_text("[description]\\nname = ...\\n[dependencies]\\n...")
    (baltik / "src" / "marker.cpp").write_text("// XD ...")

...and to make multi-baltik / nested-baltik / dep-graph fixtures readable.

A mock baltik produced by ``make_baltik`` is **shape-identical** to a
real baltik: same ``project.cfg`` layout, same ``src/`` subtree
convention, same ``// XD ...`` comment syntax. Tests built on it
exercise the real ``trustify.projects`` / scanner code paths — no
mocking of internal helpers.

Usage
-----

Basic baltik with no dependencies::

    from tests._baltik_fixtures import make_baltik

    with tempfile.TemporaryDirectory() as d:
        baltik = make_baltik(Path(d), "my_baltik")
        # → d/my_baltik/project.cfg, name = my_baltik

With dependencies::

    dep = make_baltik(root, "dep")
    main = make_baltik(root, "main", deps={"dep": dep})

Nested baltiks (sub inside parent)::

    parent = make_baltik(root, "parent")
    sub = make_baltik(parent, "sub")
    # → root/parent/project.cfg + root/parent/sub/project.cfg

With XD content (each block becomes a separate ``.xd`` file in
``<baltik>/src/``)::

    make_baltik(
        root,
        "kw_provider",
        xd_blocks=[
            "// XD foo objet_lecture foo NO_BRACE Foo keyword.\\n",
            "// XD bar foo bar NO_BRACE Bar refines foo.\\n",
        ],
    )

The ``mock_workspace`` context manager bundles the tempdir + a small
``.baltik(...)`` builder for tests that need several baltiks::

    with mock_workspace() as ws:
        dep = ws.baltik("dep")
        main = ws.baltik("main", deps={"dep": dep})
        ws.run_under(main, ...)
"""

from __future__ import annotations

import contextlib
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping


def make_baltik(
    root: Path,
    name: str,
    *,
    deps: Mapping[str, Path | str] | None = None,
    xd_blocks: list[str] | None = None,
    description_extras: Mapping[str, str] | None = None,
) -> Path:
    """Create a mock baltik at ``root/name`` with a real ``project.cfg``.

    :param root: containing directory. May itself be a baltik directory
        (returned by an earlier ``make_baltik`` call) — that nests the
        new baltik inside the parent on disk. Created if missing.
    :param name: baltik name. Written into ``project.cfg``
        ``[description].name`` *and* used as the directory basename.
    :param deps: optional ``{dep_name: path_or_string}`` mapping. Each
        entry is written verbatim into ``project.cfg`` ``[dependencies]``,
        so paths may be ``Path`` objects, absolute strings, relative
        strings like ``"../sibling_baltik"``, or even baltik-style shell
        expressions like ``"`pwd`/sub"``. The dep dir must exist *if*
        the test that uses this fixture exercises validation —
        ``expand_dependencies`` will raise if it doesn't.
    :param xd_blocks: optional list of XD declaration strings. Each
        string is dropped into ``<baltik>/src/markers_<i>.xd`` verbatim.
        Use this when the test depends on the schema actually carrying
        baltik-specific keywords (e.g. multi-origin doc generation).
    :param description_extras: optional ``{key: value}`` pairs added to
        the ``[description]`` section beyond ``name`` (typical extras:
        ``author``, ``executable``, ``cpp_flags``, ``ld_flags``,
        ``kernel``). Skip when defaults are fine.
    :returns: the absolute path of the created baltik directory.
    """
    baltik = (root / name).resolve()
    baltik.mkdir(parents=True, exist_ok=True)

    description_lines = [f"name : {name}"]
    if description_extras:
        description_lines.extend(f"{k} : {v}" for k, v in description_extras.items())

    cfg_parts = ["[description]", *description_lines, ""]
    if deps:
        cfg_parts.append("[dependencies]")
        cfg_parts.extend(f"{dep_name} : {dep_path}" for dep_name, dep_path in deps.items())
        cfg_parts.append("")
    (baltik / "project.cfg").write_text("\n".join(cfg_parts))

    # `src/` is part of the canonical baltik layout — `validate_projects`
    # rejects baltiks that lack it. Always create it so the fixture
    # matches what `baltik_configure` produces in practice; tests that
    # need keyword content add `.xd` files via `xd_blocks`.
    src = baltik / "src"
    src.mkdir(exist_ok=True)
    if xd_blocks:
        for i, block in enumerate(xd_blocks):
            content = block if block.endswith("\n") else block + "\n"
            (src / f"markers_{i}.xd").write_text(content)

    return baltik


@dataclass
class Workspace:
    """Tempdir-backed container yielded by ``mock_workspace``.

    ``root`` is the temp directory (an absolute, resolved ``Path``).
    Call ``.baltik(name, ...)`` to create baltiks under it; signature
    matches ``make_baltik`` minus the ``root`` argument.
    """

    root: Path
    _baltiks: dict[str, Path] = field(default_factory=dict)

    def baltik(
        self,
        name: str,
        *,
        parent: Path | str | None = None,
        deps: Mapping[str, Path | str] | None = None,
        xd_blocks: list[str] | None = None,
        description_extras: Mapping[str, str] | None = None,
    ) -> Path:
        """Create a baltik under this workspace.

        :param parent: optional containing directory (must be inside the
            workspace). Defaults to ``workspace.root``; pass an existing
            baltik path to nest the new one inside it.
        """
        anchor = self.root if parent is None else Path(parent)
        path = make_baltik(
            anchor,
            name,
            deps=deps,
            xd_blocks=xd_blocks,
            description_extras=description_extras,
        )
        self._baltiks[name] = path
        return path


@contextlib.contextmanager
def mock_workspace() -> Iterator[Workspace]:
    """Yield a ``Workspace`` backed by a fresh ``tempfile.TemporaryDirectory``.

    The directory is cleaned up automatically on exit. Most tests in
    this codebase open a ``TemporaryDirectory`` and immediately build
    one or more baltiks under it — this context manager wraps both
    halves so tests indent one level shallower and the baltik shape
    stays consistent.
    """
    with tempfile.TemporaryDirectory() as d:
        yield Workspace(root=Path(d).resolve())
