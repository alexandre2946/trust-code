"""Focused unit tests for trustify.doc.render helpers (no $TRUST_ROOT)."""

import unittest
from types import SimpleNamespace


class TestDefaultAnnotation(unittest.TestCase):
    """`_default_annotation` surfaces the documented default recorded by the
    code generator in `json_schema_extra['trust_default']`."""

    def test_renders_trust_default(self):
        from trustify.doc.render import _default_annotation

        fld = SimpleNamespace(json_schema_extra={"trust_default": "yy"})
        self.assertEqual(_default_annotation(fld), ", *default:* `yy`")

    def test_renders_int_trust_default(self):
        from trustify.doc.render import _default_annotation

        # entier(into=...) / entier(max=...) record an int default.
        fld = SimpleNamespace(json_schema_extra={"trust_default": 2})
        self.assertEqual(_default_annotation(fld), ", *default:* `2`")

    def test_empty_without_metadata(self):
        from trustify.doc.render import _default_annotation

        self.assertEqual(_default_annotation(SimpleNamespace(json_schema_extra=None)), "")
        self.assertEqual(_default_annotation(SimpleNamespace(json_schema_extra={})), "")
        # An unrelated extra key must not produce an annotation.
        self.assertEqual(_default_annotation(SimpleNamespace(json_schema_extra={"other": 1})), "")


if __name__ == "__main__":
    unittest.main()
