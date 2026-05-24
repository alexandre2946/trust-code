"""Tests for source-range plumbing in trustify (line/char positions on
parsed slices).
"""

import unittest
from dataclasses import FrozenInstanceError

from trustify import SourceRange


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
        with self.assertRaises(FrozenInstanceError):
            r.start_line = 5

    def test_source_range_equality(self):
        self.assertEqual(SourceRange(0, 0, 0, 1), SourceRange(0, 0, 0, 1))
        self.assertNotEqual(SourceRange(0, 0, 0, 1), SourceRange(0, 0, 0, 2))

    # ----- Tokenizer content positions -----------------------------------

    def test_tokenize_simple_positions(self):
        from trustify import TRUSTParser

        p = TRUSTParser()
        p.tokenize("a b\n  c")
        # Confirm the existing tokenization behavior is preserved:
        self.assertEqual(p.tabToken, ["a", " b", "\n", " ", " c"])
        # contentLine / contentCol track where each token's content
        # (after leading whitespace) starts, 0-based.
        self.assertEqual(p.contentLine, [0, 0, 1, 1, 1])
        self.assertEqual(p.contentCol, [0, 2, 0, 1, 2])

    def test_tokenize_internal_newline_in_leading_ws(self):
        # Tokens whose leading whitespace contains '\n' report the line
        # *after* the newline.  The tokenizer splits on space before newline,
        # so "foo\n  bar" produces 4 tokens.
        from trustify import TRUSTParser

        p = TRUSTParser()
        p.tokenize("foo\n  bar")
        # tabToken: ['foo', '\n', ' ', ' bar']
        self.assertEqual(p.tabToken, ["foo", "\n", " ", " bar"])
        self.assertEqual(p.contentLine, [0, 1, 1, 1])
        self.assertEqual(p.contentCol, [0, 0, 1, 2])

    def test_tokenize_tab_handling(self):
        from trustify import TRUSTParser

        p = TRUSTParser()
        p.tokenize("\tfoo")
        # tabToken: ['', '\tfoo']
        self.assertEqual(p.tabToken, ["", "\tfoo"])
        self.assertEqual(p.contentLine, [0, 0])
        self.assertEqual(p.contentCol, [0, 1])

    # ----- TRUSTTokens.range + lastReadTokens ----------------------------

    def test_trust_tokens_default_range_is_none(self):
        from trustify.core.trust_parser import TRUSTTokens

        self.assertIsNone(TRUSTTokens().range)
        self.assertIsNone(TRUSTTokens(low=["x"], orig=["x"]).range)

    def test_trust_tokens_explicit_range(self):
        from trustify.core.trust_parser import TRUSTTokens

        r = SourceRange(0, 0, 0, 3)
        t = TRUSTTokens(low=["foo"], orig=["foo"], range=r)
        self.assertEqual(t.range, r)

    def test_last_read_tokens_populates_range(self):
        from trustify import TRUSTParser, TRUSTStream

        p = TRUSTParser()
        p.tokenize("foo bar")
        s = TRUSTStream(p)
        s.nextLow()  # consume 'foo'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 0, 0, 3))
        s.nextLow()  # consume 'bar'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 4, 0, 7))

    def test_last_read_tokens_multiline(self):
        from trustify import TRUSTParser, TRUSTStream

        p = TRUSTParser()
        p.tokenize("foo\n  bar")
        s = TRUSTStream(p)
        s.nextLow()  # consume 'foo'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 0, 0, 3))
        s.nextLow()  # consume 'bar'
        self.assertEqual(s.lastReadTokens().range, SourceRange(1, 2, 1, 5))

    def test_last_read_tokens_multiline_quoted(self):
        # _mergeQuoted rejoins a quoted string that spans newlines into ONE
        # token whose body contains '\n'. The reported end must land on the
        # line where the token actually ends (audit #2) — not on the start
        # line at an out-of-bounds column.
        from trustify import TRUSTParser, TRUSTStream

        p = TRUSTParser()
        p.tokenize('system "echo hello\ngoodbye" 42')
        s = TRUSTStream(p)
        s.nextLow()  # consume 'system'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 0, 0, 6))
        s.nextLow()  # consume the multi-line quoted token
        # body spans line 0 (`"echo hello`) into line 1 (`goodbye"`, 8 chars)
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 7, 1, 8))
        s.nextLow()  # consume '42' — on the second line, after `goodbye" `
        self.assertEqual(s.lastReadTokens().range, SourceRange(1, 9, 1, 11))

    def test_stream_pos_multiline_quoted(self):
        # base._stream_pos shares the same end-of-token computation and
        # feeds TrustifyParseError.end_line/end_col. It must also place the
        # end on the line the token ends on (audit #2). Returns a 4-tuple
        # (line, col, end_line, end_col).
        from trustify import TRUSTParser, TRUSTStream
        from trustify.core.base import _stream_pos

        p = TRUSTParser()
        p.tokenize('system "echo hello\ngoodbye" 42')
        s = TRUSTStream(p)
        s.nextLow()  # consume 'system'
        s.probeNextLow()  # probe the quoted token → it is the "current" token
        self.assertEqual(_stream_pos(s), (0, 7, 1, 8))

    def test_last_read_tokens_empty_stream(self):
        from trustify import TRUSTParser, TRUSTStream

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
        from trustify.core.trust_parser import TRUSTTokens

        a = TRUSTTokens(low=["a"], orig=["a"], range=SourceRange(0, 0, 0, 1))
        b = TRUSTTokens(low=["b"], orig=["b"], range=SourceRange(0, 2, 0, 3))
        j = TRUSTTokens.Join([a, b])
        self.assertEqual(j.range, SourceRange(0, 0, 0, 3))
        self.assertEqual(j.low(), ["a", "b"])
        self.assertEqual(j.orig(), ["a", "b"])

    def test_join_with_none_yields_none(self):
        from trustify.core.trust_parser import TRUSTTokens

        a = TRUSTTokens(low=["a"], orig=["a"], range=SourceRange(0, 0, 0, 1))
        b = TRUSTTokens(low=["b"], orig=["b"])  # range=None
        self.assertIsNone(TRUSTTokens.Join([a, b]).range)
        self.assertIsNone(TRUSTTokens.Join([b, a]).range)

    def test_join_empty_yields_none(self):
        from trustify.core.trust_parser import TRUSTTokens

        self.assertIsNone(TRUSTTokens.Join([]).range)

    def test_join_multiline(self):
        from trustify.core.trust_parser import TRUSTTokens

        a = TRUSTTokens(low=["a"], orig=["a"], range=SourceRange(0, 0, 0, 1))
        b = TRUSTTokens(low=["b"], orig=["b"], range=SourceRange(2, 4, 2, 5))
        j = TRUSTTokens.Join([a, b])
        self.assertEqual(j.range, SourceRange(0, 0, 2, 5))

    # ----- Abstract_Parser.get_token_range -------------------------------

    def test_get_token_range_returns_range_or_none(self):
        # Use a fake subclass — no TRAD2 / generated module needed for
        # this unit test. Real-parse coverage lives in trustify-lsp's
        # integration tests where the generated module is already wired up.
        from trustify.core.base import Abstract_Parser
        from trustify.core.trust_parser import TRUSTTokens

        class FakeParser(Abstract_Parser):
            pass

        p = FakeParser()
        p._tokens["foo"] = TRUSTTokens(low=["foo"], orig=["foo"], range=SourceRange(0, 0, 0, 3))
        p._tokens["synthetic"] = TRUSTTokens(low=["{"], orig=[" {\n"])

        self.assertEqual(p.get_token_range("foo"), SourceRange(0, 0, 0, 3))
        self.assertIsNone(p.get_token_range("synthetic"))
        self.assertIsNone(p.get_token_range("missing"))

    # ----- Save / restore correctness ------------------------------------

    def test_range_after_restore_then_consume(self):
        # After save → probe → restore → consume, the next lastReadTokens
        # range must reflect the *consumed* slice, not the probed one.
        from trustify import TRUSTParser, TRUSTStream

        p = TRUSTParser()
        p.tokenize("foo bar baz")
        s = TRUSTStream(p)
        s.nextLow()  # consume 'foo'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 0, 0, 3))
        s.save("X")
        s.nextLow()  # consume 'bar'
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 4, 0, 7))
        s.restore("X")
        s.nextLow()  # re-consume 'bar' from the saved point
        # After re-consumption, _prevIdx is updated and the range matches
        # the freshly consumed slice (which is the same 'bar' content).
        self.assertEqual(s.lastReadTokens().range, SourceRange(0, 4, 0, 7))


class TestTokenizeQuotedAndCommented(unittest.TestCase):
    """Audit 1.4: `_mergeQuoted` used to count `"` on every token,
    including those inside `# ... #` and `/* ... */` comments — a
    standalone `"` token within a comment would toggle `in_quot` and
    slurp subsequent (real) tokens into a phantom quoted span. The
    position walker also assumed token bodies never contain `\\n`,
    which broke for merged multi-line quoted strings.

    TRUST C++ supports both: System::interpreter reads tokens until
    one ends with `"` (joining with single spaces, so `"line1\\nline2"`
    legally turns into the shell command `line1 line2`), and the C++
    tokenizer has zero quote-awareness so quotes inside `#...#`
    comments are not special. trustify must match both behaviours.
    """

    def test_quote_inside_hash_comment_does_not_merge_across(self):
        # `# foo " bar #` would, pre-fix, toggle in_quot on the
        # standalone `"` token and slurp `bar`, `#`, and everything
        # up to the next stray quote into a phantom merged token —
        # corrupting tokenization.
        from trustify import TRUSTParser

        p = TRUSTParser()
        p.tokenize('# foo " bar #\nreal_token')
        # `real_token` must survive as its own token at line 1 — the
        # merge mustn't reach across the closing `#`.
        self.assertIn("real_token", [t.strip() for t in p.tabToken])
        idx = next(i for i, t in enumerate(p.tabToken) if t.strip() == "real_token")
        self.assertEqual(p.contentLine[idx], 1)

    def test_quote_inside_slash_star_comment_does_not_merge_across(self):
        # Same hazard for `/* ... " ... */` comments.
        from trustify import TRUSTParser

        p = TRUSTParser()
        p.tokenize('/* a " b */\nreal_token')
        self.assertIn("real_token", [t.strip() for t in p.tabToken])
        idx = next(i for i, t in enumerate(p.tabToken) if t.strip() == "real_token")
        self.assertEqual(p.contentLine[idx], 1)

    def test_quoted_string_outside_comment_still_merges(self):
        # Regression guard for the `system "rm -rf foo"` use case the
        # merge was designed to handle.
        from trustify import TRUSTParser

        p = TRUSTParser()
        p.tokenize('system "rm -rf foo"')
        # The merged tabToken stream should contain one token for the
        # whole `"rm -rf foo"` span.
        merged = [t.strip() for t in p.tabToken if t.strip()]
        self.assertIn('"rm -rf foo"', merged)

    def test_multi_line_quoted_string_positions_track_internal_newlines(self):
        # TRUST C++ accepts `System "line1\nline2"` (newlines collapse
        # to single spaces in the shell command). trustify must keep
        # contentLine/contentCol accurate for tokens AFTER the merged
        # multi-line span — the position walker has to handle `\n`
        # inside merged token bodies, not just leading whitespace.
        from trustify import TRUSTParser

        p = TRUSTParser()
        p.tokenize('system "line1\nline2"\nafter')
        # Find the token containing `after` — it lives on line 2 of
        # the input (`system`, then a 2-line merged quoted string,
        # then `after`).
        idx = next(i for i, t in enumerate(p.tabToken) if t.strip() == "after")
        self.assertEqual(p.contentLine[idx], 2, f"after-token line wrong; tabToken={p.tabToken!r}")
        self.assertEqual(p.contentCol[idx], 0)


if __name__ == "__main__":
    unittest.main()
