"""Tests for source-range plumbing in trustify (line/char positions on
parsed slices). See docs/superpowers/specs/2026-05-08-trustify-source-ranges-design.md
in the trustify-lsp repo.
"""

import unittest

from trustify.trust_parser import SourceRange


class TestCase(unittest.TestCase):

    # ----- SourceRange dataclass -----------------------------------------

    def test_source_range_fields(self):
        r = SourceRange(0, 1, 2, 3)
        self.assertEqual(r.start_line, 0)
        self.assertEqual(r.start_char, 1)
        self.assertEqual(r.end_line, 2)
        self.assertEqual(r.end_char, 3)

    def test_source_range_frozen(self):
        r = SourceRange(0, 0, 0, 1)
        with self.assertRaises(Exception):
            r.start_line = 5  # frozen → FrozenInstanceError

    def test_source_range_equality(self):
        self.assertEqual(SourceRange(0, 0, 0, 1), SourceRange(0, 0, 0, 1))
        self.assertNotEqual(SourceRange(0, 0, 0, 1), SourceRange(0, 0, 0, 2))


    # ----- Tokenizer content positions -----------------------------------

    def test_tokenize_simple_positions(self):
        from trustify.trust_parser import TRUSTParser
        p = TRUSTParser()
        p.tokenize("a b\n  c")
        # Confirm the existing tokenization behavior is preserved:
        self.assertEqual(p.tabToken, ['a', ' b', '\n', ' ', ' c'])
        # contentLine / contentCol track where each token's content
        # (after leading whitespace) starts, 0-based.
        self.assertEqual(p.contentLine, [0, 0, 1, 1, 1])
        self.assertEqual(p.contentCol,  [0, 2, 0, 1, 2])

    def test_tokenize_internal_newline_in_leading_ws(self):
        # Tokens whose leading whitespace contains '\n' report the line
        # *after* the newline.  The tokenizer splits on space before newline,
        # so "foo\n  bar" produces 4 tokens.
        from trustify.trust_parser import TRUSTParser
        p = TRUSTParser()
        p.tokenize("foo\n  bar")
        # tabToken: ['foo', '\n', ' ', ' bar']
        self.assertEqual(p.tabToken, ['foo', '\n', ' ', ' bar'])
        self.assertEqual(p.contentLine, [0, 1, 1, 1])
        self.assertEqual(p.contentCol,  [0, 0, 1, 2])

    def test_tokenize_tab_handling(self):
        from trustify.trust_parser import TRUSTParser
        p = TRUSTParser()
        p.tokenize("\tfoo")
        # tabToken: ['', '\tfoo']
        self.assertEqual(p.tabToken, ['', '\tfoo'])
        self.assertEqual(p.contentLine, [0, 0])
        self.assertEqual(p.contentCol,  [0, 1])


    # ----- TRUSTTokens.range + lastReadTokens ----------------------------

    def test_trust_tokens_default_range_is_none(self):
        from trustify.trust_parser import TRUSTTokens
        self.assertIsNone(TRUSTTokens().range)
        self.assertIsNone(TRUSTTokens(low=["x"], orig=["x"]).range)

    def test_trust_tokens_explicit_range(self):
        from trustify.trust_parser import TRUSTTokens
        r = SourceRange(0, 0, 0, 3)
        t = TRUSTTokens(low=["foo"], orig=["foo"], range=r)
        self.assertEqual(t.range, r)

    def test_last_read_tokens_populates_range(self):
        from trustify.trust_parser import TRUSTParser, TRUSTStream
        p = TRUSTParser()
        p.tokenize("foo bar")
        s = TRUSTStream(p)
        s.nextLow()                              # consume 'foo'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 0, 0, 3))
        s.nextLow()                              # consume 'bar'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 4, 0, 7))

    def test_last_read_tokens_multiline(self):
        from trustify.trust_parser import TRUSTParser, TRUSTStream
        p = TRUSTParser()
        p.tokenize("foo\n  bar")
        s = TRUSTStream(p)
        s.nextLow()                              # consume 'foo'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 0, 0, 3))
        s.nextLow()                              # consume 'bar'
        self.assertEqual(s.lastReadTokens().range, SourceRange(1, 2, 1, 5))

    def test_last_read_tokens_empty_stream(self):
        from trustify.trust_parser import TRUSTParser, TRUSTStream
        p = TRUSTParser()
        p.tokenize("foo")
        s = TRUSTStream(p)
        # Nothing consumed yet → empty TRUSTTokens with range=None
        t = s.lastReadTokens()
        self.assertIsNone(t.range)
        self.assertEqual(t.low(), [])
        self.assertEqual(t.orig(), [])

    # ----- TRUSTTokens.Join ----------------------------------------------

    def test_join_with_ranges_spans(self):
        from trustify.trust_parser import TRUSTTokens
        a = TRUSTTokens(low=["a"], orig=["a"], range=SourceRange(0, 0, 0, 1))
        b = TRUSTTokens(low=["b"], orig=["b"], range=SourceRange(0, 2, 0, 3))
        j = TRUSTTokens.Join([a, b])
        self.assertEqual(j.range, SourceRange(0, 0, 0, 3))
        self.assertEqual(j.low(), ["a", "b"])
        self.assertEqual(j.orig(), ["a", "b"])

    def test_join_with_none_yields_none(self):
        from trustify.trust_parser import TRUSTTokens
        a = TRUSTTokens(low=["a"], orig=["a"], range=SourceRange(0, 0, 0, 1))
        b = TRUSTTokens(low=["b"], orig=["b"])  # range=None
        self.assertIsNone(TRUSTTokens.Join([a, b]).range)
        self.assertIsNone(TRUSTTokens.Join([b, a]).range)

    def test_join_empty_yields_none(self):
        from trustify.trust_parser import TRUSTTokens
        self.assertIsNone(TRUSTTokens.Join([]).range)

    def test_join_multiline(self):
        from trustify.trust_parser import TRUSTTokens
        a = TRUSTTokens(low=["a"], orig=["a"], range=SourceRange(0, 0, 0, 1))
        b = TRUSTTokens(low=["b"], orig=["b"], range=SourceRange(2, 4, 2, 5))
        j = TRUSTTokens.Join([a, b])
        self.assertEqual(j.range, SourceRange(0, 0, 2, 5))

    # ----- Abstract_Parser.get_token_range -------------------------------

    def test_get_token_range_returns_range_or_none(self):
        # Use a fake subclass — no TRAD2 / generated module needed for
        # this unit test. Real-parse coverage lives in trustify-lsp's
        # integration tests where the generated module is already wired up.
        from trustify.base import Abstract_Parser
        from trustify.trust_parser import TRUSTTokens

        class FakeParser(Abstract_Parser):
            pass

        p = FakeParser()
        p._tokens["foo"] = TRUSTTokens(low=["foo"], orig=["foo"],
                                       range=SourceRange(0, 0, 0, 3))
        p._tokens["synthetic"] = TRUSTTokens(low=["{"], orig=[" {\n"])

        self.assertEqual(p.get_token_range("foo"), SourceRange(0, 0, 0, 3))
        self.assertIsNone(p.get_token_range("synthetic"))
        self.assertIsNone(p.get_token_range("missing"))

    # ----- Save / restore correctness ------------------------------------

    def test_range_after_restore_then_consume(self):
        # After save → probe → restore → consume, the next lastReadTokens
        # range must reflect the *consumed* slice, not the probed one.
        from trustify.trust_parser import TRUSTParser, TRUSTStream
        p = TRUSTParser()
        p.tokenize("foo bar baz")
        s = TRUSTStream(p)
        s.nextLow()                              # consume 'foo'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 0, 0, 3))
        s.save("X")
        s.nextLow()                              # consume 'bar'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 4, 0, 7))
        s.restore("X")
        s.nextLow()                              # re-consume 'bar' from the saved point
        # After re-consumption, _prevIdx is updated and the range matches
        # the freshly consumed slice (which is the same 'bar' content).
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 4, 0, 7))


if __name__ == "__main__":
    unittest.main()
