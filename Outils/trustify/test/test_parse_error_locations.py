"""Unit tests for TrustifyParseError (structured location for parser errors).

Each test crafts a minimal `.data` snippet that triggers exactly one
`kind` of parser error and asserts the structured fields. Run via:

    cd Outils/trustify && make test
"""

import typing
import unittest

from trustify.misc_utilities import ClassFactory, TrustifyException, TrustifyParseError

from test.test_rw_elementary import TestCase as ElementaryTestCase
from trustify import base
from trustify.base import Dataset_Parser


class TestStructuredFields(ElementaryTestCase):
    """One test per `kind`. Each parses a hand-crafted snippet via the
    same `import_and_gen_stream` helper test_rw_elementary uses, then
    asserts the captured `TrustifyParseError`'s structured fields.
    """

    def _parse_and_capture(self, snippet):
        """Run Dataset_Parser.ReadFromTokens; return the raised error."""
        stream = self.import_and_gen_stream(snippet, simplify=False)
        try:
            Dataset_Parser.ReadFromTokens(stream)
        except TrustifyParseError as exc:
            return exc
        self.fail("expected TrustifyParseError but no exception was raised")

    def _parse_keyword_and_capture(self, snippet):
        """Mirror of test_rw_elementary.generic_test: dispatch to the
        first keyword's parser class instead of Dataset_Parser. Needed
        for snippets whose top-level keyword is not an Interprete (e.g.
        `coucou`, `read_med_bidon`) and would otherwise be swallowed by
        Declaration_Parser.
        """
        stream = self.import_and_gen_stream(snippet, simplify=False)
        try:
            cls_nam = stream.probeNextLow()
            ze_cls = ClassFactory.GetParserClassFromName(cls_nam)
            ze_cls.ReadFromTokens(stream)
        except TrustifyParseError as exc:
            return exc
        self.fail("expected TrustifyParseError but no exception was raised")

    def test_invalid_keyword(self):
        snippet = "dimension 2\nzzznotakeyword\n"
        err = self._parse_and_capture(snippet)
        self.assertEqual(err.kind, "invalid-keyword")
        self.assertEqual(err.token, "zzznotakeyword")
        self.assertEqual(err.line, 1)
        self.assertIsNotNone(err.col)

    def test_unexpected_attribute(self):
        # `coucou` is the keyword used at test_rw_elementary:654 for
        # the unexpected-attribute regex assertion.
        snippet = "coucou { not_an_attr 0 }\n"
        err = self._parse_keyword_and_capture(snippet)
        self.assertEqual(err.kind, "unexpected-attribute")
        self.assertEqual(err.token, "not_an_attr")
        self.assertEqual(err.attr_name, "not_an_attr")

    def test_mandatory_missing(self):
        # `read_med_bidon` is the keyword used at test_rw_elementary:474
        # for the mandatory-attribute regex assertion.
        snippet = "read_med_bidon { }\n"
        err = self._parse_keyword_and_capture(snippet)
        self.assertEqual(err.kind, "mandatory-missing")
        self.assertIsNotNone(err.attr_name)

    def test_expected_brace(self):
        # Same shape as test_rw_elementary:311.
        snippet = "read_med_bidon not_a_brace\n"
        err = self._parse_keyword_and_capture(snippet)
        self.assertEqual(err.kind, "expected-brace")
        self.assertEqual(err.token, "not_a_brace")

    def test_invalid_value_int(self):
        # Builtin int parser, mirroring test_rw_elementary:176.
        with self.assertRaises(TrustifyParseError) as cm:
            self.builtin_test(int, "toto")
        self.assertEqual(cm.exception.kind, "invalid-value")
        self.assertEqual(cm.exception.token, "toto")

    def test_dimension_undefined(self):
        # Mirrors test_rw_elementary.test_builtin_list(fixed=True): the
        # `size_is_dim` annotation routes through AbstractSizeIsDim.readListSize,
        # which raises TrustifyParseError(kind="dimension-undefined") when no
        # dimension was set first.
        base.Dimension_Parser._DIMENSION = -1
        with self.assertRaises(TrustifyParseError) as cm:
            self.builtin_test(
                typing.Annotated[typing.List[float], "size_is_dim"],
                "1.0 2.0 3.0",
                simplify=False,
            )
        self.assertEqual(cm.exception.kind, "dimension-undefined")
        self.assertEqual(cm.exception.line, 0)


class TestTrustifyParseErrorClass(unittest.TestCase):
    """Smoke tests for the dataclass-like structured exception itself."""

    def test_subclass_of_trustify_exception(self):
        self.assertTrue(issubclass(TrustifyParseError, TrustifyException))

    def test_constructs_with_all_fields(self):
        err = TrustifyParseError(
            "boom",
            file_name="x.data",
            line=3,
            col=12,
            end_line=3,
            end_col=20,
            token="bogus",
            attr_name=None,
            kind="invalid-keyword",
        )
        self.assertEqual(str(err), "boom")
        self.assertEqual(err.file_name, "x.data")
        self.assertEqual(err.line, 3)
        self.assertEqual(err.col, 12)
        self.assertEqual(err.end_line, 3)
        self.assertEqual(err.end_col, 20)
        self.assertEqual(err.token, "bogus")
        self.assertIsNone(err.attr_name)
        self.assertEqual(err.kind, "invalid-keyword")

    def test_existing_handlers_still_match(self):
        """Backwards-compat: `except TrustifyException` keeps catching."""
        err = TrustifyParseError(
            "boom",
            file_name="x.data",
            line=0,
            col=None,
            end_line=None,
            end_col=None,
            token=None,
            attr_name=None,
            kind="invalid-keyword",
        )
        try:
            raise err
        except TrustifyException as caught:
            self.assertIs(caught, err)
        else:
            self.fail("TrustifyException did not catch TrustifyParseError")
