"""Unit tests for trustify.core.source_location."""

import unittest
from dataclasses import FrozenInstanceError

from trustify.core.source_location import (
    ProjectCollisionError,
    SourceLocation,
    build_projects_dict,
    relativize,
)


class TestSourceLocation(unittest.TestCase):
    def test_round_trip_to_dict(self):
        loc = SourceLocation(project="trust", path="src/Foo.cpp", line=42)
        self.assertEqual(SourceLocation.from_dict(loc.to_dict()), loc)

    def test_resolve_returns_absolute_path_and_line(self):
        loc = SourceLocation(project="trust", path="src/Foo.cpp", line=42)
        abs_root, line = loc.resolve({"trust": "/abs/root"})
        self.assertEqual(abs_root, "/abs/root/src/Foo.cpp")
        self.assertEqual(line, 42)

    def test_resolve_returns_none_for_unconfigured_project(self):
        loc = SourceLocation(project="baltik_x", path="src/Foo.cpp", line=1)
        self.assertIsNone(loc.resolve({"trust": "/abs/root"}))

    def test_frozen_dataclass(self):
        loc = SourceLocation(project="trust", path="x", line=1)
        with self.assertRaises(FrozenInstanceError):
            loc.project = "other"  # frozen


class TestBuildProjectsDict(unittest.TestCase):
    def test_trust_root_only(self):
        d = build_projects_dict(trust_root="/home/u/trust-code", projects=[])
        self.assertEqual(d, {"trust": "/home/u/trust-code"})

    def test_trust_root_plus_baltik(self):
        d = build_projects_dict(
            trust_root="/home/u/trust-code",
            projects=["/home/u/baltiks/baltik_a"],
        )
        self.assertEqual(
            d,
            {
                "trust": "/home/u/trust-code",
                "baltik_a": "/home/u/baltiks/baltik_a",
            },
        )

    def test_basename_collision_raises(self):
        with self.assertRaises(ProjectCollisionError) as cm:
            build_projects_dict(
                trust_root="/home/u/trust",
                projects=["/elsewhere/trust"],
            )
        # message names both colliding paths
        self.assertIn("/home/u/trust", str(cm.exception))
        self.assertIn("/elsewhere/trust", str(cm.exception))

    def test_no_trust_root_uses_basenames_only(self):
        d = build_projects_dict(trust_root=None, projects=["/abs/baltik_a"])
        self.assertEqual(d, {"baltik_a": "/abs/baltik_a"})

    def test_baltik_name_comes_from_project_cfg_when_present(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            baltik = Path(tmp) / "lower-case-dir"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname : MyBaltik\n")
            d = build_projects_dict(trust_root=None, projects=[str(baltik)])
            self.assertEqual(d, {"MyBaltik": str(baltik)})

    def test_baltik_falls_back_to_basename_without_project_cfg(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            baltik = Path(tmp) / "no_cfg_here"
            baltik.mkdir()
            d = build_projects_dict(trust_root=None, projects=[str(baltik)])
            self.assertEqual(d, {"no_cfg_here": str(baltik)})

    def test_name_collision_via_project_cfg_raises(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            a = Path(tmp) / "dir_a"
            b = Path(tmp) / "dir_b"
            a.mkdir()
            b.mkdir()
            # Both baltiks declare the same name in their project.cfg —
            # build_projects_dict must surface the collision the same way
            # it does for basename collisions.
            (a / "project.cfg").write_text("[description]\nname : Same\n")
            (b / "project.cfg").write_text("[description]\nname : Same\n")
            with self.assertRaises(ProjectCollisionError):
                build_projects_dict(trust_root=None, projects=[str(a), str(b)])


class TestRelativize(unittest.TestCase):
    def test_finds_owning_project(self):
        projects = {"trust": "/home/u/trust-code", "baltik_a": "/home/u/baltik_a"}
        proj, rel = relativize("/home/u/trust-code/src/Foo/Bar.cpp", projects)
        self.assertEqual(proj, "trust")
        self.assertEqual(rel, "src/Foo/Bar.cpp")

    def test_baltik_path_resolved(self):
        projects = {"trust": "/home/u/trust-code", "baltik_a": "/home/u/baltik_a"}
        proj, rel = relativize("/home/u/baltik_a/src/Foo.cpp", projects)
        self.assertEqual(proj, "baltik_a")
        self.assertEqual(rel, "src/Foo.cpp")

    def test_path_outside_all_projects_raises(self):
        projects = {"trust": "/home/u/trust-code"}
        with self.assertRaises(ValueError) as cm:
            relativize("/somewhere/else/x.cpp", projects)
        self.assertIn("/somewhere/else/x.cpp", str(cm.exception))

    def test_longest_prefix_wins(self):
        # If one project root is nested inside another, the nested one wins.
        projects = {
            "trust": "/home/u/trust",
            "embedded": "/home/u/trust/baltiks/embedded",
        }
        proj, rel = relativize("/home/u/trust/baltiks/embedded/src/X.cpp", projects)
        self.assertEqual(proj, "embedded")
        self.assertEqual(rel, "src/X.cpp")

    def test_path_equals_project_root(self):
        # A path exactly equal to the project root produces an empty relative path.
        proj, rel = relativize("/home/u/trust", {"trust": "/home/u/trust"})
        self.assertEqual(proj, "trust")
        self.assertEqual(rel, "")

    def test_path_equals_project_root_trailing_slash(self):
        # A trailing slash on the input is normalized before comparison.
        proj, rel = relativize("/home/u/trust/", {"trust": "/home/u/trust"})
        self.assertEqual(proj, "trust")
        self.assertEqual(rel, "")

    def test_relativize_resolves_symlinks_on_both_sides(self):
        """Audit 4.3: `relativize` must canonicalize the input path
        and each project root via `os.path.realpath` before the
        prefix-match test. Otherwise a clone whose `trust_root` was
        passed in symlink form (or whose file paths came in via a
        symlink) produces a different (project, relative) string than
        a clone that uses the canonical form, which would leak into
        `source_locations.json` and break cache-hash stability across
        clones — the CLAUDE.md "byte-identical schema files across
        clones" guarantee.
        """
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            canonical = Path(tmp) / "trust-canonical"
            (canonical / "src").mkdir(parents=True)
            (canonical / "src" / "Foo.cpp").write_text("")

            link = Path(tmp) / "trust-link"
            link.symlink_to(canonical)

            # Scenario 1: projects dict has the canonical root; the file
            # path came in via the symlink. Result must match scenario 2.
            proj_a, rel_a = relativize(
                str(link / "src" / "Foo.cpp"),
                {"trust": str(canonical)},
            )
            # Scenario 2: projects dict has the symlinked root; the file
            # path is canonical.
            proj_b, rel_b = relativize(
                str(canonical / "src" / "Foo.cpp"),
                {"trust": str(link)},
            )
            self.assertEqual((proj_a, rel_a), (proj_b, rel_b))
            self.assertEqual(proj_a, "trust")
            # The relative part is the only thing that should appear in
            # the eventual JSON — confirm it's the on-disk relative form,
            # not something containing the symlinked component name.
            self.assertEqual(rel_a, "src/Foo.cpp")

            # And the trivial no-symlink case must still work — both
            # sides canonical, no surprises.
            proj_c, rel_c = relativize(
                str(canonical / "src" / "Foo.cpp"),
                {"trust": str(canonical)},
            )
            self.assertEqual(proj_c, "trust")
            self.assertEqual(rel_c, "src/Foo.cpp")


if __name__ == "__main__":
    unittest.main()
