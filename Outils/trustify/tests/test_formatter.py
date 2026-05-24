"""Unit tests for the TRUST .data formatter — ported from the TypeScript extension."""

from __future__ import annotations

import unittest

from trustify import format_dataset
from trustify.formatter import (
    FoldingRegion,
    FormatContext,
    analyze_line,
    compute_format_context,
    find_folding_regions,
    format_on_type,
)


class TestAnalyzeLine(unittest.TestCase):
    """Ported from suite `analyzeLine` in clients/vscode/src/test/extension.test.ts."""

    def test_plain_text_has_no_braces(self):
        assert analyze_line("dimension 3") == (0, 0, 0, False)

    def test_standalone_open_brace(self):
        assert analyze_line("{") == (1, 0, 0, False)

    def test_standalone_close_brace(self):
        assert analyze_line("}") == (0, 1, 1, False)

    def test_keyword_followed_by_open_brace(self):
        assert analyze_line("conditions_initiales {") == (1, 0, 0, False)

    def test_inline_balanced_block_no_net_change(self):
        assert analyze_line("convection { negligeable }") == (1, 1, 0, False)

    def test_nested_inline_balanced_blocks(self):
        result = analyze_line("solveur_pression GCP { precond ssor { omega -1.5 } seuil 1.e+20 impr }")
        assert result == (2, 2, 0, False)

    def test_single_line_comment_braces_inside_ignored(self):
        assert analyze_line("# tmax  10. { } #") == (0, 0, 0, False)

    def test_single_line_comment_then_code_with_brace(self):
        assert analyze_line("# ignored # keyword {") == (1, 0, 0, False)

    def test_line_starting_with_hash_begins_multiline_comment(self):
        assert analyze_line("#") == (0, 0, 0, True)

    def test_multiline_comment_start_with_text_and_braces_all_ignored(self):
        assert analyze_line("# BEGIN PARTITION") == (0, 0, 0, True)

    def test_leading_close_brace_not_counted_when_preceded_by_text(self):
        assert analyze_line("foo }") == (0, 1, 0, False)

    def test_multiple_leading_close_braces(self):
        assert analyze_line("} }") == (0, 2, 2, False)

    # /* ... */ comments — same brace-suppression contract as # ... #.

    def test_single_line_slash_star_comment_braces_inside_ignored(self):
        assert analyze_line("/* { } */") == (0, 0, 0, False)

    def test_slash_star_comment_then_code_with_brace(self):
        assert analyze_line("/* ignored */ keyword {") == (1, 0, 0, False)

    def test_line_starting_with_open_slash_star_begins_multiline_comment(self):
        assert analyze_line("/*") == (0, 0, 0, True)

    def test_resume_in_slash_star_state(self):
        # When in_comment="slash", braces stay suppressed until `*/`.
        assert analyze_line("{ still } */ then {", in_comment="slash") == (
            1,
            0,
            0,
            False,
        )


class TestComputeContextBefore(unittest.TestCase):
    """Ported from suite `computeContextBefore` in extension.test.ts."""

    def test_empty_text_gives_zero_level(self):
        assert compute_format_context("") == FormatContext(level=0, in_block_comment=False)

    def test_flat_lines_give_zero_level(self):
        assert compute_format_context("a\nb\n") == FormatContext(level=0, in_block_comment=False)

    def test_open_brace_increases_level(self):
        assert compute_format_context("Lire pb\n{\n") == FormatContext(level=1, in_block_comment=False)

    def test_nested_open_braces_accumulate_level(self):
        assert compute_format_context("a\n{\nb\n{\n") == FormatContext(level=2, in_block_comment=False)

    def test_closed_block_returns_to_previous_level(self):
        assert compute_format_context("a\n{\nb\n}\n") == FormatContext(level=0, in_block_comment=False)

    def test_detects_open_block_comment(self):
        assert compute_format_context("a\n#\n") == FormatContext(level=0, in_block_comment=True, comment_mode="hash")

    def test_detects_open_slash_star_comment(self):
        assert compute_format_context("a\n/*\n") == FormatContext(level=0, in_block_comment=True, comment_mode="slash")

    def test_closed_slash_star_comment_is_not_open(self):
        assert compute_format_context("a\n/*\ncomment\n*/\n") == FormatContext(
            level=0, in_block_comment=False, comment_mode=None
        )

    def test_braces_inside_slash_star_comment_do_not_count(self):
        assert compute_format_context("/*\n{\n*/\n") == FormatContext(
            level=0, in_block_comment=False, comment_mode=None
        )

    def test_closed_block_comment_is_not_open(self):
        assert compute_format_context("a\n#\ncomment\n#\n") == FormatContext(level=0, in_block_comment=False)

    def test_braces_inside_block_comment_do_not_count(self):
        assert compute_format_context("#\n{\n#\n") == FormatContext(level=0, in_block_comment=False)


def _fmt_expanded(text: str) -> str:
    # Drop the trailing `\n` that `format_dataset` now adds on every
    # non-empty output (audit 6.4 — guaranteed-terminator invariant
    # covered by TestFormatEolNormalization). The TestFormat*
    # assertions below were written before that invariant existed and
    # check indentation / comment handling only; the rstrip here keeps
    # those assertions tight on what they actually want to verify.
    return format_dataset(text).removesuffix("\n")


class TestFormatBasicIndentation(unittest.TestCase):
    """Ported from suite `formatTrustData — basic indentation`."""

    def test_flat_content_is_unchanged(self):
        text = "dimension 3\nProbleme_FT pb\nDomaine DOM"
        assert _fmt_expanded(text) == text

    def test_standalone_braces_indent_their_content(self):
        text = "Lire sch\n{\ntinit 0.\n}"
        assert _fmt_expanded(text) == "Lire sch\n{\n    tinit 0.\n}"

    def test_open_brace_at_end_of_keyword_line_indents_content(self):
        text = "conditions_initiales {\nfoo\n}"
        assert _fmt_expanded(text) == "conditions_initiales {\n    foo\n}"

    def test_nested_blocks_produce_increasing_indentation(self):
        text = "outer\n{\ninner\n{\ndeep\n}\n}"
        assert _fmt_expanded(text) == "outer\n{\n    inner\n    {\n        deep\n    }\n}"

    def test_inline_balanced_block_does_not_change_indent_level(self):
        text = "Lire sch\n{\nconvection { negligeable }\n}"
        assert _fmt_expanded(text) == "Lire sch\n{\n    convection { negligeable }\n}"

    def test_leading_whitespace_stripped_and_replaced(self):
        text = "   Lire sch\n   {\n      tinit 0.\n   }"
        assert _fmt_expanded(text) == "Lire sch\n{\n    tinit 0.\n}"

    def test_empty_lines_adjacent_to_braces_are_removed(self):
        text = "Lire sch\n{\n\ntinit 0.\n\n}"
        assert _fmt_expanded(text) == "Lire sch\n{\n    tinit 0.\n}"

    def test_indent_level_never_goes_negative(self):
        assert _fmt_expanded("}\nfoo") == "}\nfoo"


def _fmt(text: str) -> str:
    # Drop the trailing `\n` that `format_dataset` now adds on every
    # non-empty output (audit 6.4). Same rationale as `_fmt_expanded`
    # above: existing assertions focus on whitespace / comment /
    # blank-line behaviour, not on the new EOL invariant.
    return format_dataset(text).removesuffix("\n")


class TestFormatAlreadyInlineBlocks(unittest.TestCase):
    """Already-inline single-line blocks pass through unchanged."""

    def test_already_inline_block_unchanged(self):
        assert _fmt("convection { negligeable }") == "convection { negligeable }"


class TestFormatCommentHandling(unittest.TestCase):
    """Ported from suite `formatTrustData — comment handling`."""

    def test_single_line_comment_prevents_block_collapse(self):
        text = "Lire sch\n{\n# tmax 10. #\ntinit 0.\n}"
        expected = "Lire sch\n{\n    # tmax 10. #\n    tinit 0.\n}"
        assert _fmt(text) == expected

    def test_braces_inside_single_line_comment_do_not_affect_indent(self):
        text = "Lire sch\n{\n# { fake block } #\ntinit 0.\n}"
        expected = "Lire sch\n{\n    # { fake block } #\n    tinit 0.\n}"
        assert _fmt(text) == expected

    def test_multiline_comment_braces_inside_ignored_block_not_collapsed(self):
        text = "foo\n#\n{ ignored brace }\n#\nbar"
        assert _fmt(text) == "foo\n#\n{ ignored brace }\n#\nbar"

    def test_commented_partition_section(self):
        text = "# BEGIN PARTITION\nPartition DOM\n{\nNb_parts 2\n}\nFin\nEND PARTITION #\nfoo"
        # Entire section between # markers is a block comment — unchanged
        assert _fmt(text) == text

    def test_indented_content_inside_hash_comment_is_preserved(self):
        # User-authored indentation inside a `# ... #` comment must survive
        # formatting unchanged: comments often contain commented-out
        # valid syntax that the user expects to be able to uncomment.
        text = "# BEGIN PARTITION\nPartition DOM\n{\n    Nb_parts 2\n    zones_name DOM\n}\nEND PARTITION #\nfoo"
        assert _fmt(text) == text

    def test_indented_content_inside_slash_star_comment_is_preserved(self):
        text = "/*\nPartition DOM\n{\n    Nb_parts 2\n}\n*/\nfoo"
        assert _fmt(text) == text

    def test_tabbed_content_inside_comment_is_preserved_verbatim(self):
        # The formatter does not rewrite indentation inside comments at
        # all: tabs, spaces, mixed indents — all left alone.
        text = "# BEGIN\n\tindented_with_tab\n\t  mixed_tab_and_spaces\nEND #"
        assert _fmt(text) == text

    def test_comment_inside_real_block_keeps_inner_indent(self):
        # The outer block re-indents normally, but content inside the
        # comment retains the original column the user typed it at.
        text = (
            "outer\n"
            "{\n"
            "# BEGIN INNER\n"
            "        deeply_indented_in_comment\n"
            "    less_indented_in_comment\n"
            "END INNER #\n"
            "    real_kw value\n"
            "}"
        )
        expected = (
            "outer\n"
            "{\n"
            "    # BEGIN INNER\n"
            "        deeply_indented_in_comment\n"
            "    less_indented_in_comment\n"
            "END INNER #\n"
            "    real_kw value\n"
            "}"
        )
        assert _fmt(text) == expected

    def test_commented_block_not_collapsed_even_when_shape_looks_collapsible(self):
        # `Postraiter_domaine\n{\n    x\n}` is exactly the shape Pass 2
        # collapses — but since the whole section sits inside a
        # `# BEGIN ... END #` comment the collapse must be refused.
        text = "# BEGIN PARTITION\nPostraiter_domaine\n{\n    fichier mesh.lata\n}\nEND PARTITION #"
        assert _fmt(text) == text

    def test_indent_preservation_inside_comment_is_idempotent(self):
        text = "# BEGIN PARTITION\nPartition DOM\n{\n    Nb_parts 2\n}\nEND PARTITION #"
        once = _fmt(text)
        twice = _fmt(once)
        assert once == twice == text

    def test_range_starting_in_comment_preserves_inner_indent(self):
        # When formatting a range whose start sits inside an open
        # multi-line comment, the same preservation rule applies.
        text = "    still in comment\n    closing line #"
        ctx = FormatContext(level=0, in_block_comment=True, comment_mode="hash")
        assert format_dataset(text, ctx) == text


class TestFormatIdempotency(unittest.TestCase):
    """Ported from suite `formatTrustData — idempotency`."""

    def test_formatting_twice_gives_same_result(self):
        text = "Lire pb\n{\n    solved_equations\n    {\n        Navier_Stokes ns\n    }\n}"
        once = _fmt(text)
        twice = _fmt(once)
        assert once == twice

    def test_realistic_excerpt_is_idempotent(self):
        text = "\n".join(
            [
                "Lire pb",
                "{",
                "solved_equations",
                "{",
                "Navier_Stokes_FT_Disc              hydraulique",
                "Transport_Interfaces_FT_Disc       balle",
                "}",
                "hydraulique",
                "{",
                "modele_turbulence sous_maille_wale",
                "{",
                "Cw 0.",
                "}",
                "solveur_pression GCP { precond ssor { omega -1.5 } seuil 1.e+20 impr }",
                "convection           { negligeable }",
                "}",
                "}",
            ]
        )
        once = _fmt(text)
        twice = _fmt(once)
        assert once == twice


class TestFormatEolNormalization(unittest.TestCase):
    """Audit 6.4: every formatted file ends with exactly one `\\n` and
    contains no `\\r` (CRLF inputs are normalised to LF). git and POSIX
    cat treat newline-less files as malformed; without normalisation,
    formatting a CRLF file under Windows-CRLF habits left mixed-EOL
    output that round-trips badly.

    These tests call `format_dataset` directly (bypassing the `_fmt`
    helper above) because the helper rstrips the trailing newline to
    keep the older indentation-focused assertions tight — that would
    defeat the very invariant this class is testing.
    """

    def test_input_without_trailing_newline_gets_one(self):
        # Common shape from editor "no newline at end of file" warnings.
        out = format_dataset("Lire pb {\nfoo bar\n}")
        self.assertTrue(out.endswith("\n"), f"expected trailing newline, got: {out!r}")

    def test_input_with_trailing_newline_keeps_exactly_one(self):
        # No doubling on already-newline-terminated input.
        out = format_dataset("Lire pb {\nfoo bar\n}\n")
        self.assertTrue(out.endswith("\n"))
        self.assertFalse(out.endswith("\n\n"), f"trailing newline doubled: {out!r}")

    def test_empty_input_stays_empty(self):
        # No spurious newline added to an empty file.
        self.assertEqual(format_dataset(""), "")

    def test_crlf_input_normalised_to_lf(self):
        out = format_dataset("Lire pb {\r\nfoo bar\r\n}\r\n")
        self.assertNotIn("\r", out)
        self.assertTrue(out.endswith("\n"))

    def test_trailing_newline_added_is_idempotent(self):
        # format(format(x)) == format(x) even when the first pass had
        # to ADD the trailing newline.
        text = "Lire pb {\nfoo bar\n}"  # no trailing newline
        once = format_dataset(text)
        twice = format_dataset(once)
        self.assertEqual(once, twice)


class TestFormatBlankLineNormalization(unittest.TestCase):
    """Pass 1 normalization of blank lines."""

    def test_single_blank_line_between_content_is_preserved(self):
        assert _fmt("a\n\nb") == "a\n\nb"

    def test_multiple_consecutive_blank_lines_collapsed_to_one(self):
        assert _fmt("a\n\n\nb") == "a\n\nb"
        assert _fmt("a\n\n\n\nb") == "a\n\nb"

    def test_blanks_between_content_inside_block_normalized(self):
        text = "Lire sch\n{\nfoo\n\n\nbar\n}"
        assert _fmt(text) == "Lire sch\n{\n    foo\n\n    bar\n}"

    def test_block_with_only_blank_lines_becomes_empty_braces(self):
        text = "Lire sch\n{\n\n\n}"
        assert _fmt(text) == "Lire sch\n{\n}"


class TestFormatRangeFormatting(unittest.TestCase):
    """Ported from suite `formatTrustData — range formatting via startContext`."""

    def test_content_at_level_1_indented_correctly(self):
        text = "foo\nbar"
        ctx = FormatContext(level=1, in_block_comment=False)
        assert format_dataset(text, ctx) == "    foo\n    bar"

    def test_closing_brace_within_range_dedents_from_start_context_level(self):
        text = "foo\n}"
        ctx = FormatContext(level=1, in_block_comment=False)
        assert format_dataset(text, ctx) == "    foo\n}"

    def test_block_inside_range_indented_correctly(self):
        text = "keyword\n{\nfoo\n}"
        ctx = FormatContext(level=1, in_block_comment=False)
        assert format_dataset(text, ctx) == "    keyword\n    {\n        foo\n    }"

    def test_content_starting_inside_block_comment_emitted_verbatim(self):
        text = "still in comment\n#\nnormal line"
        ctx = FormatContext(level=0, in_block_comment=True)
        assert format_dataset(text, ctx) == text


class TestFindFoldingRegions(unittest.TestCase):
    """Folding ranges derived from brace depth and multi-line `#...#` comments."""

    def test_empty_text_has_no_regions(self):
        assert find_folding_regions("") == []

    def test_flat_content_has_no_regions(self):
        assert find_folding_regions("foo\nbar\nbaz") == []

    def test_inline_block_on_single_line_has_no_region(self):
        assert find_folding_regions("foo { bar }") == []

    def test_multiline_block_with_standalone_open_brace(self):
        text = "foo\n{\nbar\n}"
        assert find_folding_regions(text) == [FoldingRegion(1, 3, False)]

    def test_multiline_block_with_trailing_open_brace(self):
        text = "foo {\nbar\n}"
        assert find_folding_regions(text) == [FoldingRegion(0, 2, False)]

    def test_nested_blocks_produce_both_outer_and_inner_regions(self):
        text = "outer {\n    inner {\n        deep\n    }\n}"
        regions = find_folding_regions(text)
        assert FoldingRegion(0, 4, False) in regions
        assert FoldingRegion(1, 3, False) in regions
        assert len(regions) == 2

    def test_multiline_block_comment_is_tagged_as_comment(self):
        text = "foo\n# start\nmiddle\nend #\nbar"
        assert find_folding_regions(text) == [FoldingRegion(1, 3, True)]

    def test_single_line_hash_pair_has_no_region(self):
        assert find_folding_regions("foo # comment # bar") == []

    def test_braces_inside_block_comment_are_ignored(self):
        text = "#\n{ not real }\n#"
        assert find_folding_regions(text) == [FoldingRegion(0, 2, True)]

    def test_unmatched_open_brace_has_no_region(self):
        assert find_folding_regions("foo {\nbar") == []

    def test_unmatched_close_brace_has_no_region(self):
        assert find_folding_regions("foo\n}") == []

    def test_comment_and_block_coexist_at_correct_lines(self):
        text = "outer {\n#\ninside\n#\n}"
        regions = find_folding_regions(text)
        assert FoldingRegion(0, 4, False) in regions
        assert FoldingRegion(1, 3, True) in regions
        assert len(regions) == 2

    def test_sibling_blocks_do_not_nest(self):
        text = "a {\nfoo\n}\nb {\nbar\n}"
        regions = find_folding_regions(text)
        assert FoldingRegion(0, 2, False) in regions
        assert FoldingRegion(3, 5, False) in regions
        assert len(regions) == 2

    def test_multiline_slash_star_comment_is_tagged_as_comment(self):
        text = "foo\n/* start\nmiddle\nend */\nbar"
        assert find_folding_regions(text) == [FoldingRegion(1, 3, True)]

    def test_single_line_slash_star_has_no_region(self):
        assert find_folding_regions("foo /* comment */ bar") == []

    def test_braces_inside_slash_star_block_comment_are_ignored(self):
        text = "/*\n{ not real }\n*/"
        assert find_folding_regions(text) == [FoldingRegion(0, 2, True)]


class TestFormatOnType(unittest.TestCase):
    """Single-line re-indent after typing a character that changes indent level."""

    def test_closing_brace_outdented_to_outer_level(self):
        text = "foo\n{\n    }"
        assert format_on_type(text, 2) == ("}", len("    }"))

    def test_closing_brace_already_at_correct_indent_returns_none(self):
        text = "foo\n{\n    bar\n}"
        assert format_on_type(text, 3) is None

    def test_nested_closing_brace_indented_to_middle_level(self):
        text = "outer\n{\n    inner {\n}"
        assert format_on_type(text, 3) == ("    }", len("}"))

    def test_line_without_leading_close_brace_returns_none(self):
        text = "foo\n{\n    bar"
        assert format_on_type(text, 2) is None

    def test_line_inside_block_comment_returns_none(self):
        text = "#\n}\nfoo"
        assert format_on_type(text, 1) is None

    def test_line_out_of_range_returns_none(self):
        assert format_on_type("foo", 5) is None

    def test_empty_line_returns_none(self):
        text = "foo\n{\n"
        assert format_on_type(text, 2) is None

    def test_text_after_leading_close_brace_preserved(self):
        text = "foo\n{\n    } extra"
        assert format_on_type(text, 2) == ("} extra", len("    } extra"))


class TestStripBlanksAdjacentToBraces(unittest.TestCase):
    """Blank lines adjacent to structural braces are removed."""

    def test_blank_right_after_open_brace_is_removed(self):
        text = "keyword\n{\n\nfoo\n}"
        assert _fmt_expanded(text) == "keyword\n{\n    foo\n}"

    def test_blank_right_before_close_brace_is_removed(self):
        text = "keyword\n{\nfoo\n\n}"
        assert _fmt_expanded(text) == "keyword\n{\n    foo\n}"

    def test_multiple_blanks_after_open_brace_all_removed(self):
        text = "keyword\n{\n\n\n\nfoo\n}"
        assert _fmt_expanded(text) == "keyword\n{\n    foo\n}"

    def test_multiple_blanks_before_close_brace_all_removed(self):
        text = "keyword\n{\nfoo\n\n\n\n}"
        assert _fmt_expanded(text) == "keyword\n{\n    foo\n}"

    def test_trailing_blank_after_content_is_preserved_as_final_newline(self):
        # `foo\n\n` becomes `foo\n` (trailing newline preserved, run of blanks
        # collapsed to one). The expected `foo\n` is the single trailing
        # newline both the old blank-collapsing logic and the new
        # always-terminate-the-file invariant produce; bypass _fmt_expanded
        # (which would strip that newline) and assert against format_dataset
        # directly so the equality lands exactly as written.
        text = "foo\n\n"
        assert format_dataset(text) == "foo\n"

    def test_trailing_blanks_after_open_brace_at_eof_are_dropped(self):
        # File ends with `{\n\n` — every trailing blank is "right after `{`",
        # so they all go and the file loses its trailing newline.
        text = "keyword\n{\n\n"
        assert _fmt_expanded(text) == "keyword\n{"

    def test_blank_after_inline_open_brace_keyword_is_removed(self):
        text = "keyword {\n\nfoo\n}"
        assert _fmt_expanded(text) == "keyword {\n    foo\n}"

    def test_blank_before_close_brace_with_content_after_is_removed(self):
        # `} something` still starts with `}` so the blank before it is dropped.
        text = "keyword\n{\nfoo\n\n} bar"
        assert _fmt_expanded(text) == "keyword\n{\n    foo\n} bar"

    def test_blank_between_sibling_blocks_is_preserved(self):
        # The blank sits after a `}` and before a non-`}` line; neither
        # side fires the strip rule, so the blank survives.
        text = "keyword {\nfoo\n}\n\nother {\nbar\n}"
        assert _fmt_expanded(text) == "keyword {\n    foo\n}\n\nother {\n    bar\n}"

    def test_blank_inside_multiline_comment_is_preserved(self):
        # The `{` is inside the comment, so the blank is not adjacent to a
        # structural brace and must be kept.
        text = "# comment opens here {\n\nstill in comment #"
        assert _fmt_expanded(text) == text

    def test_blank_in_range_starting_in_block_comment_is_preserved(self):
        # Range formatting where the slice starts inside a multi-line comment.
        # The first line ends with `{` but it's inside the comment, so the
        # blank that follows must not be stripped.
        text = "still in comment {\n\nmore comment #"
        ctx = FormatContext(level=0, in_block_comment=True)
        assert format_dataset(text, ctx) == text

    def test_brace_inside_single_line_comment_does_not_trigger_strip(self):
        # `# ... #` self-closes on the same line; the trailing `{` is real.
        # Sanity check that the existing behavior on the structural brace
        # still kicks in.
        text = "a # ignored # b {\n\nfoo\n}"
        assert _fmt_expanded(text) == "a # ignored # b {\n    foo\n}"


if __name__ == "__main__":
    unittest.main()
