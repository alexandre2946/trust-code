"""Unit tests for the shared comment scanner."""

from __future__ import annotations

import unittest

from trustify.core.comment_scanner import (
    comment_state_at_line,
    is_position_in_comment,
    scan_line,
)


class TestScanLineHashStyle(unittest.TestCase):
    def test_plain_text_no_ranges(self):
        ranges, end = scan_line("convection negligeable")
        assert ranges == []
        assert end is None

    def test_single_hash_comment(self):
        ranges, end = scan_line("# hello #")
        assert ranges == [(0, 9)]
        assert end is None

    def test_open_hash_no_close_state_carries(self):
        ranges, end = scan_line("foo # open")
        assert ranges == [(4, 10)]
        assert end == "hash"

    def test_resume_in_hash_state_until_close(self):
        ranges, end = scan_line("body close # tail", in_comment="hash")
        assert ranges == [(0, 12)]
        assert end is None

    def test_two_hash_comments_on_one_line(self):
        ranges, end = scan_line("# a # mid # b #")
        # `# a #` covers cols 0..4 inclusive → half-open (0, 5);
        # `# b #` covers cols 10..14 inclusive → half-open (10, 15).
        assert ranges == [(0, 5), (10, 15)]
        assert end is None


class TestScanLineSlashStyle(unittest.TestCase):
    def test_single_slash_star_comment(self):
        ranges, end = scan_line("/* hi */")
        assert ranges == [(0, 8)]
        assert end is None

    def test_empty_slash_star_comment(self):
        ranges, end = scan_line("/**/")
        assert ranges == [(0, 4)]
        assert end is None

    def test_open_slash_star_no_close_state_carries(self):
        ranges, end = scan_line("foo /* open")
        assert ranges == [(4, 11)]
        assert end == "slash"

    def test_resume_in_slash_state_until_close(self):
        ranges, end = scan_line("body close */ tail", in_comment="slash")
        assert ranges == [(0, 13)]
        assert end is None

    def test_lone_star_not_close_in_slash_state(self):
        # A single `*` followed by anything but `/` does NOT close.
        ranges, end = scan_line("body * still open", in_comment="slash")
        assert ranges == [(0, 17)]
        assert end == "slash"


class TestScanLineMixed(unittest.TestCase):
    def test_hash_inside_open_slash_is_plain_text(self):
        # `#` inside `/* */` must not toggle anything.
        ranges, end = scan_line("/* hash # inside */")
        assert ranges == [(0, 19)]
        assert end is None

    def test_slash_open_inside_hash_is_plain_text(self):
        # `/*` inside `# ... #` must not open a slash comment.
        ranges, end = scan_line("# /* fake */ #")
        assert ranges == [(0, 14)]
        assert end is None

    def test_adjacent_mixed_comments(self):
        # `# h # /* s */` — two comments, different styles.
        ranges, end = scan_line("# h # /* s */")
        assert ranges == [(0, 5), (6, 13)]
        assert end is None


class TestScanLineUpTo(unittest.TestCase):
    def test_up_to_truncates_open_range_to_boundary(self):
        ranges, end = scan_line("# open and far past", up_to=5)
        assert ranges == [(0, 5)]
        assert end == "hash"

    def test_up_to_at_zero_returns_empty(self):
        ranges, end = scan_line("# anything", up_to=0)
        assert ranges == []
        assert end is None

    def test_up_to_past_end_treated_as_eol(self):
        ranges, end = scan_line("# x #", up_to=100)
        assert ranges == [(0, 5)]
        assert end is None

    def test_up_to_inside_slash_marker(self):
        # `up_to=1` falls inside the 2-char `/*` opener — must not enter
        # slash state since the `*` isn't reached.
        ranges, end = scan_line("/* foo */", up_to=1)
        assert ranges == []
        assert end is None


class TestCommentStateAtLine(unittest.TestCase):
    def test_line_zero_is_always_none(self):
        assert comment_state_at_line(["# foo", "bar"], 0) is None

    def test_carries_hash_open_across_line(self):
        lines = ["foo # open", "still", "close # tail"]
        assert comment_state_at_line(lines, 1) == "hash"
        assert comment_state_at_line(lines, 2) == "hash"
        # Past the closing `#` on line 2 we'd be back to None; the
        # function reports the state at the *start* of the target line.

    def test_carries_slash_open_across_line(self):
        lines = ["/* open", "middle", "close */"]
        assert comment_state_at_line(lines, 1) == "slash"
        assert comment_state_at_line(lines, 2) == "slash"

    def test_closed_before_target_is_none(self):
        lines = ["# closed #", "clean", ""]
        assert comment_state_at_line(lines, 1) is None
        assert comment_state_at_line(lines, 2) is None


class TestIsPositionInComment(unittest.TestCase):
    def test_outside_any_comment(self):
        assert is_position_in_comment(["foo bar"], 0, 0) is False
        assert is_position_in_comment(["foo bar"], 0, 4) is False

    def test_inside_hash_comment(self):
        assert is_position_in_comment(["# hi #"], 0, 3) is True

    def test_marker_char_is_inside(self):
        # Opening `#` is at col 0 — the marker counts as inside.
        assert is_position_in_comment(["# hi #"], 0, 0) is True

    def test_inside_slash_comment(self):
        assert is_position_in_comment(["/* hi */"], 0, 4) is True

    def test_inside_multi_line_slash_comment(self):
        lines = ["/* open", "middle", "close */"]
        assert is_position_in_comment(lines, 1, 3) is True
        assert is_position_in_comment(lines, 2, 0) is True


if __name__ == "__main__":
    unittest.main()
