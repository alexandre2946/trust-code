"""Tests for trustify.doc.origin.OriginClassifier.

After the source-locations refactor, _infoMain triples carry a project
name directly, so the classifier reduces to a project_name -> label dict
lookup. These tests cover the new API surface only — the legacy
path-prefix matcher is gone.
"""

import unittest

from trustify.core.source_location import SourceLocation
from trustify.doc.origin import OriginClassifier, origin_slug


class TestOriginClassifierCore(unittest.TestCase):
    def test_constructor_stores_mapping(self):
        c = OriginClassifier({"trust": "TRUST", "baltik_a": "baltik_a"})
        self.assertEqual(c.classify_project("trust"), "TRUST")
        self.assertEqual(c.classify_project("baltik_a"), "baltik_a")

    def test_constructor_copies_mapping(self):
        m = {"trust": "TRUST"}
        c = OriginClassifier(m)
        m["trust"] = "OTHER"
        self.assertEqual(c.classify_project("trust"), "TRUST")

    def test_labels_preserves_insertion_order(self):
        c = OriginClassifier({"baltik_a": "baltik_a", "baltik_b": "baltik_b", "trust": "TRUST"})
        self.assertEqual(c.labels, ["baltik_a", "baltik_b", "TRUST"])

    def test_is_multi_false_for_single_label(self):
        c = OriginClassifier({"trust": "TRUST"})
        self.assertFalse(c.is_multi)

    def test_is_multi_true_for_multiple_labels(self):
        c = OriginClassifier({"trust": "TRUST", "baltik_a": "baltik_a"})
        self.assertTrue(c.is_multi)


class TestOriginClassifierClassifyProject(unittest.TestCase):
    def test_returns_label_for_known_project(self):
        c = OriginClassifier({"trust": "TRUST"})
        self.assertEqual(c.classify_project("trust"), "TRUST")

    def test_returns_none_for_unknown_project(self):
        c = OriginClassifier({"trust": "TRUST"})
        self.assertIsNone(c.classify_project("baltik_a"))

    def test_returns_none_for_empty_project(self):
        c = OriginClassifier({"trust": "TRUST"})
        self.assertIsNone(c.classify_project(""))


class TestOriginClassifierClassifyLocation(unittest.TestCase):
    def test_dispatches_to_classify_project(self):
        c = OriginClassifier({"trust": "TRUST", "baltik_a": "baltik_a"})
        loc = SourceLocation(project="baltik_a", path="src/x.cpp", line=1)
        self.assertEqual(c.classify_location(loc), "baltik_a")

    def test_returns_none_for_unconfigured_project(self):
        c = OriginClassifier({"trust": "TRUST"})
        loc = SourceLocation(project="baltik_x", path="src/x.cpp", line=1)
        self.assertIsNone(c.classify_location(loc))

    def test_returns_none_for_loc_without_project_attr(self):
        c = OriginClassifier({"trust": "TRUST"})
        # Object without a `project` attribute returns None gracefully.
        self.assertIsNone(c.classify_location(object()))


class TestOriginClassifierFromProjects(unittest.TestCase):
    def test_no_projects_with_trust_root(self):
        c = OriginClassifier.from_projects(projects=None, trust_root="/abs/trust-code")
        self.assertEqual(c.classify_project("trust"), "TRUST")
        self.assertFalse(c.is_multi)

    def test_baltik_label_is_basename(self):
        c = OriginClassifier.from_projects(projects=["/abs/baltik_a"], trust_root="/abs/trust-code")
        self.assertEqual(c.classify_project("baltik_a"), "baltik_a")
        self.assertEqual(c.classify_project("trust"), "TRUST")

    def test_projects_come_before_trust_in_label_order(self):
        c = OriginClassifier.from_projects(
            projects=["/abs/baltik_a", "/abs/baltik_b"],
            trust_root="/abs/trust-code",
        )
        self.assertEqual(c.labels, ["baltik_a", "baltik_b", "TRUST"])

    def test_no_trust_root_yields_baltiks_only(self):
        c = OriginClassifier.from_projects(projects=["/abs/baltik_a"], trust_root=None)
        self.assertEqual(c.classify_project("baltik_a"), "baltik_a")
        self.assertIsNone(c.classify_project("trust"))

    def test_no_inputs_yields_empty(self):
        c = OriginClassifier.from_projects(projects=None, trust_root=None)
        self.assertEqual(c.labels, [])
        self.assertFalse(c.is_multi)

    def test_baltik_label_comes_from_project_cfg_when_present(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            baltik = Path(tmp) / "lower-case-dir"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname : TrioCFD\n")
            c = OriginClassifier.from_projects(projects=[str(baltik)], trust_root="/abs/trust-code")
            # The doc-side classifier and the schema-side build_projects_dict
            # MUST agree on the same name, otherwise the doc generator's
            # _origin_of() can't find the keyword's project label.
            self.assertEqual(c.classify_project("TrioCFD"), "TrioCFD")
            self.assertIsNone(c.classify_project("lower-case-dir"))
            self.assertEqual(c.labels, ["TrioCFD", "TRUST"])

    def test_baltik_label_falls_back_to_basename_without_project_cfg(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            baltik = Path(tmp) / "no_cfg_baltik"
            baltik.mkdir()
            c = OriginClassifier.from_projects(projects=[str(baltik)], trust_root="/abs/trust-code")
            self.assertEqual(c.classify_project("no_cfg_baltik"), "no_cfg_baltik")


class TestOriginSlug(unittest.TestCase):
    def test_uppercase_label_lowercased(self):
        self.assertEqual(origin_slug("TRUST"), "trust")

    def test_special_chars_collapsed(self):
        self.assertEqual(origin_slug("baltik a/b-c"), "baltik_a_b_c")

    def test_trailing_underscores_stripped(self):
        self.assertEqual(origin_slug("__foo--"), "foo")


if __name__ == "__main__":
    unittest.main()
