"""Unit tests for trustify.doc.extras — pure logic, synthetic filesystem fixtures."""

import tempfile
import unittest
from pathlib import Path


class TestExtrasIndexCore(unittest.TestCase):
    """Constructor + empty-state behaviour."""

    def test_empty_index_has_no_used_or_unused(self):
        from trustify.doc.extras import ExtrasIndex

        idx = ExtrasIndex({}, {})
        self.assertEqual(idx.used_images(), [])
        self.assertEqual(idx.unused_files(), [])

    def test_constructor_stores_maps(self):
        from trustify.doc.extras import ExtrasIndex

        md = {"foo": Path("/x/foo.md")}
        img = {"bar.png": Path("/x/bar.png")}
        idx = ExtrasIndex(md, img)
        # Internal access is fine — sanity check the wiring.
        self.assertEqual(idx._md, md)
        self.assertEqual(idx._img, img)

    def test_origin_of_defaults_to_none(self):
        # Bare constructor (no file_origins) — origin lookup is None, so
        # callers degrade gracefully when the index wasn't built via
        # from_projects.
        from trustify.doc.extras import ExtrasIndex

        idx = ExtrasIndex({"foo": Path("/x/foo.md")}, {})
        self.assertIsNone(idx.origin_of(Path("/x/foo.md")))

    def test_constructor_stores_file_origins(self):
        from trustify.doc.extras import ExtrasIndex

        p = Path("/x/foo.md")
        idx = ExtrasIndex({"foo": p}, {}, {p: "MyProject"})
        self.assertEqual(idx.origin_of(p), "MyProject")


class TestExtrasIndexFromProjects(unittest.TestCase):
    """Scanning behaviour: partitioning, error policy, edge cases."""

    @staticmethod
    def _make_extras(root: Path, *names: str) -> Path:
        """Create ``<root>/docs/trustify/extras/`` and touch each file
        (with content ``# <name>`` for .md files, empty otherwise).
        Returns the extras dir path.
        """
        extras = root / "docs" / "trustify" / "extras"
        extras.mkdir(parents=True)
        for name in names:
            p = extras / name
            if p.suffix == ".md":
                p.write_text(f"# {p.stem}\n")
            else:
                p.write_bytes(b"")
        return extras

    def test_no_origins_returns_empty_index(self):
        from trustify.doc.extras import ExtrasIndex

        idx = ExtrasIndex.from_projects([], None)
        self.assertEqual(idx.unused_files(), [])
        self.assertEqual(idx.used_images(), [])

    def test_trust_root_extras_scanned(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            self._make_extras(tr, "foo.md", "bar.png")
            idx = ExtrasIndex.from_projects(None, str(tr))
            # Stems map for .md; filenames map for images.
            self.assertIn("foo", idx._md)
            self.assertIn("bar.png", idx._img)

    def test_projects_extras_scanned(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            proj = Path(d) / "baltik"
            proj.mkdir()
            self._make_extras(proj, "guide.md")
            idx = ExtrasIndex.from_projects([str(proj)], None)
            self.assertIn("guide", idx._md)

    def test_from_projects_records_origin_per_file(self):
        # Each scanned file remembers which project it came from, so the
        # unused-extras warning can name the origin (any baltik may now
        # ship extras). trust_root renders as 'TRUST'; a baltik renders
        # as its project name (here the dir basename, no project.cfg).
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            proj = Path(d) / "MyBaltik"
            proj.mkdir()
            self._make_extras(tr, "fromtrust.md")
            self._make_extras(proj, "frombaltik.png")
            idx = ExtrasIndex.from_projects([str(proj)], str(tr))
            labels = {p.name: idx.origin_of(p) for p in idx.unused_files()}
            self.assertEqual(labels["fromtrust.md"], "TRUST")
            self.assertEqual(labels["frombaltik.png"], "MyBaltik")

    def test_md_files_indexed_by_stem(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            self._make_extras(tr, "abc.md")
            idx = ExtrasIndex.from_projects(None, str(tr))
            # Keyed by stem, NOT by full filename.
            self.assertIn("abc", idx._md)
            self.assertNotIn("abc.md", idx._md)

    def test_image_files_indexed_by_full_filename(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            self._make_extras(tr, "diagram.svg")
            idx = ExtrasIndex.from_projects(None, str(tr))
            # Keyed by filename including extension.
            self.assertIn("diagram.svg", idx._img)
            self.assertNotIn("diagram", idx._img)

    def test_subdir_in_extras_raises(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            extras = self._make_extras(tr, "foo.md")
            (extras / "subdir").mkdir()
            with self.assertRaises(ValueError) as ctx:
                ExtrasIndex.from_projects(None, str(tr))
            self.assertIn("subdir", str(ctx.exception))

    def test_duplicate_md_across_origins_raises(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            proj = Path(d) / "baltik"
            proj.mkdir()
            self._make_extras(tr, "shared.md")
            self._make_extras(proj, "shared.md")
            with self.assertRaises(ValueError) as ctx:
                ExtrasIndex.from_projects([str(proj)], str(tr))
            msg = str(ctx.exception)
            self.assertIn("shared.md", msg)
            self.assertIn(str(tr), msg)
            self.assertIn(str(proj), msg)

    def test_duplicate_image_across_origins_raises(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            proj = Path(d) / "baltik"
            proj.mkdir()
            self._make_extras(tr, "logo.png")
            self._make_extras(proj, "logo.png")
            with self.assertRaises(ValueError):
                ExtrasIndex.from_projects([str(proj)], str(tr))

    def test_md_and_image_share_basename_is_ok(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            proj = Path(d) / "baltik"
            proj.mkdir()
            self._make_extras(tr, "thing.md")
            self._make_extras(proj, "thing.png")
            idx = ExtrasIndex.from_projects([str(proj)], str(tr))
            # md["thing"] is the .md file; img["thing.png"] is the image.
            # No collision: they're in separate namespaces.
            self.assertIn("thing", idx._md)
            self.assertIn("thing.png", idx._img)

    def test_missing_extras_dir_skipped(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            tr = Path(d) / "trust"
            tr.mkdir()
            # No docs/trustify/extras/ created.
            idx = ExtrasIndex.from_projects(None, str(tr))
            self.assertEqual(idx._md, {})
            self.assertEqual(idx._img, {})


class TestExtrasIndexResolve(unittest.TestCase):
    """Lookup methods: hit semantics, miss errors, used-marking."""

    def test_resolve_input_returns_file_content(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            md = Path(d) / "foo.md"
            md.write_text("# Hello\n")
            idx = ExtrasIndex({"foo": md}, {})
            self.assertEqual(idx.resolve_input("foo"), "# Hello\n")

    def test_resolve_input_marks_used(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            md = Path(d) / "foo.md"
            md.write_text("X")
            idx = ExtrasIndex({"foo": md}, {})
            idx.resolve_input("foo")
            self.assertEqual(idx.unused_files(), [])

    def test_resolve_input_missing_raises(self):
        from trustify.doc.extras import ExtrasIndex

        idx = ExtrasIndex({}, {})
        with self.assertRaises(ValueError) as ctx:
            idx.resolve_input("nope")
        self.assertIn("nope", str(ctx.exception))
        self.assertIn("input", str(ctx.exception))

    def test_resolve_image_returns_source_path(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            img = Path(d) / "foo.png"
            img.write_bytes(b"")
            idx = ExtrasIndex({}, {"foo.png": img})
            self.assertEqual(idx.resolve_image("foo.png"), img)

    def test_resolve_image_marks_used(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            img = Path(d) / "foo.png"
            img.write_bytes(b"")
            idx = ExtrasIndex({}, {"foo.png": img})
            idx.resolve_image("foo.png")
            self.assertEqual(idx.used_images(), [img])
            self.assertEqual(idx.unused_files(), [])

    def test_resolve_image_missing_raises(self):
        from trustify.doc.extras import ExtrasIndex

        idx = ExtrasIndex({}, {})
        with self.assertRaises(ValueError) as ctx:
            idx.resolve_image("missing.png")
        self.assertIn("missing.png", str(ctx.exception))
        self.assertIn("includeimage", str(ctx.exception))

    def test_resolve_idempotent(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            img = Path(d) / "foo.png"
            img.write_bytes(b"")
            idx = ExtrasIndex({}, {"foo.png": img})
            idx.resolve_image("foo.png")
            idx.resolve_image("foo.png")
            # Used once, but still tracked once (set semantics).
            self.assertEqual(idx.used_images(), [img])

    def test_unused_lists_md_before_images_sorted(self):
        from trustify.doc.extras import ExtrasIndex

        with tempfile.TemporaryDirectory() as d:
            md_b = Path(d) / "b.md"
            md_b.write_text("")
            md_a = Path(d) / "a.md"
            md_a.write_text("")
            img_z = Path(d) / "z.png"
            img_z.write_bytes(b"")
            img_y = Path(d) / "y.png"
            img_y.write_bytes(b"")
            idx = ExtrasIndex(
                {"b": md_b, "a": md_a},
                {"z.png": img_z, "y.png": img_y},
            )
            # md first (sorted by stem), then images (sorted by filename).
            self.assertEqual(idx.unused_files(), [md_a, md_b, img_y, img_z])


if __name__ == "__main__":
    unittest.main()
