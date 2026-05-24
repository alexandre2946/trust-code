"""Unit tests for trustify.core.naming — pure string-mangling helpers."""

import unittest


class TestNaming(unittest.TestCase):
    def test_to_pyd_name_lowercase_capitalize(self):
        from trustify.core.naming import ToPydName

        self.assertEqual(ToPydName("Read_MED"), "Read_med")
        self.assertEqual(ToPydName("interprete"), "Interprete")

    def test_to_parser_name_appends_suffix(self):
        from trustify.core.naming import ToParserName

        self.assertEqual(ToParserName("interprete"), "Interprete_Parser")

    def test_parser_suffix_constant(self):
        from trustify.core.naming import _PARSER_SUFFIX

        self.assertEqual(_PARSER_SUFFIX, "_Parser")


class TestValidVariableName(unittest.TestCase):
    """Audit 4.8: `valid_variable_name` used to silently strip every
    character outside `[A-Za-z0-9_]`, so `étude` became `tude` —
    potentially colliding with a sibling attribute called `tude`,
    with no diagnostic at generation time. TRUST keywords are
    documented as `[A-Za-z0-9_]+` only; non-conforming names should
    therefore fail loudly, pointing at the offending characters so
    the user can fix the XD tag in the source."""

    def test_pass_through_for_valid_identifier(self):
        from trustify.core.trad2_pydantic import valid_variable_name

        self.assertEqual(valid_variable_name("foo_bar_42"), "foo_bar_42")

    def test_leading_digit_is_prefixed(self):
        # Pre-existing behaviour preserved: leading digits get an `i`
        # prefix to produce a valid Python identifier.
        from trustify.core.trad2_pydantic import valid_variable_name

        self.assertEqual(valid_variable_name("9th_attr"), "i9th_attr")

    def test_python_keyword_is_suffixed(self):
        # Pre-existing behaviour preserved: Python keywords get a
        # trailing underscore so the generated field name is legal.
        from trustify.core.trad2_pydantic import valid_variable_name

        self.assertEqual(valid_variable_name("lambda"), "lambda_")

    def test_accented_char_raises(self):
        from trustify.core.trad2_pydantic import valid_variable_name

        with self.assertRaises(ValueError) as cm:
            valid_variable_name("étude")
        msg = str(cm.exception)
        # The message must name the offending characters AND the
        # original name so the XD tag is easy to locate in the source.
        self.assertIn("étude", msg)
        self.assertIn("é", msg)

    def test_punctuation_raises(self):
        from trustify.core.trad2_pydantic import valid_variable_name

        with self.assertRaises(ValueError) as cm:
            valid_variable_name("foo-bar")
        self.assertIn("foo-bar", str(cm.exception))
        self.assertIn("-", str(cm.exception))

    def test_whitespace_raises(self):
        # A space inside the name is also illegal — TRUST keywords
        # never contain whitespace, so a stray space points at an
        # XD-parse error upstream.
        from trustify.core.trad2_pydantic import valid_variable_name

        with self.assertRaises(ValueError):
            valid_variable_name("foo bar")


if __name__ == "__main__":
    unittest.main()
