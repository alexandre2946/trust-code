"""Pure-text formatter for TRUST `.data` files.

Single-pass indentation: strips existing indentation, re-indents based
on `{` / `}` depth, preserves multi-line comments verbatim, collapses
consecutive blank lines, drops blanks immediately after `{` or before
`}`. Already-inline blocks (`convection { negligeable }`) stay on one
line. Multi-line blocks are never collapsed.

Format settings are HARD-CODED — every `.data` file gets 4-space
indentation, no tabs — so multi-project collaborations agree on
whitespace by construction. Callers cannot override.

Public surface:
- `format_dataset(text, start_context=None) -> str` — full format
- `compute_format_context(text) -> FormatContext` — get the indent
  context at the end of `text`; used by range-formatting callers.
- `analyze_line(line) -> LineAnalysis` — single-line brace/comment scan.
- `find_folding_regions(text) -> list[FoldingRegion]` — fold ranges
  derived from `{...}` blocks.
- `format_on_type(text, line_num) -> tuple[str, int] | None` — re-indent
  a single line after the user types `}`.
"""

from __future__ import annotations

from typing import NamedTuple

from trustify.core.comment_scanner import CommentMode, scan_line

# Locked formatting parameters — see module docstring for rationale.
_TAB_SIZE = 4
_INSERT_SPACES = True


class LineAnalysis(NamedTuple):
    opens: int
    closes: int
    leading_closes: int
    ends_in_comment: bool


def _analyze_line_with_state(line: str, in_comment: CommentMode = None) -> tuple[LineAnalysis, CommentMode]:
    """Internal version of :func:`analyze_line` that also returns the
    precise end comment mode. Used by callers that need to thread
    ``"hash"``/``"slash"`` state across line boundaries (in particular
    to resume a multi-line ``/* ... */`` comment correctly during range
    formatting).
    """
    comment_ranges, end_state = scan_line(line, in_comment)
    masked = set()
    for s, e in comment_ranges:
        masked.update(range(s, e))

    opens = 0
    closes = 0
    leading_closes = 0
    leading_done = False

    for i, c in enumerate(line):
        if i in masked:
            # Markers and comment bodies count as "non-whitespace content
            # already seen" so a trailing close-brace stops being treated
            # as a leading close.
            leading_done = True
            continue
        if c == "}":
            closes += 1
            if not leading_done:
                leading_closes += 1
        elif c == "{":
            opens += 1
            leading_done = True
        elif c not in (" ", "\t"):
            leading_done = True

    return (
        LineAnalysis(opens, closes, leading_closes, end_state is not None),
        end_state,
    )


def analyze_line(line: str, in_comment: CommentMode = None) -> LineAnalysis:
    """Count brace and comment state for one line.

    ``in_comment`` is the comment mode entering the line (``None`` if
    not in a comment). Recognises both ``# ... #`` and ``/* ... */``.
    """
    info, _ = _analyze_line_with_state(line, in_comment)
    return info


class FormatContext(NamedTuple):
    """Indent level and comment state at a point in the document.

    ``in_block_comment`` is kept as a ``bool`` for backwards
    compatibility with callers that only need a yes/no answer.
    ``comment_mode`` carries the precise mode (``"hash"`` / ``"slash"``
    / ``None``) so range formatting can resume a multi-line comment of
    the right style.
    """

    level: int
    in_block_comment: bool
    comment_mode: CommentMode = None


def _resolve_mode(ctx: FormatContext) -> CommentMode:
    """Best-effort extraction of the comment mode from a context.

    A context constructed with ``in_block_comment=True`` but no explicit
    ``comment_mode`` defaults to ``"hash"`` — that matches the pre-/* */
    behaviour where the formatter only knew about hash comments.
    """
    if ctx.comment_mode is not None:
        return ctx.comment_mode
    if ctx.in_block_comment:
        return "hash"
    return None


def compute_format_context(text: str) -> FormatContext:
    """Return the FormatContext at the end of ``text`` — used by callers
    that need to format a range of a larger document with knowledge of
    its surrounding indent level / open comment state.
    """
    level = 0
    state: CommentMode = None

    for raw_line in text.split("\n"):
        line = raw_line.strip()
        if line == "":
            continue

        if state is not None:
            # Inside a multi-line comment — only advance the comment
            # state, no brace accounting.
            _, state = scan_line(line, state)
            continue

        info, state = _analyze_line_with_state(line, state)
        level = max(0, level - info.leading_closes)
        level = max(0, level + info.opens - info.closes + info.leading_closes)

    return FormatContext(level=level, in_block_comment=state is not None, comment_mode=state)


def format_dataset(
    text: str,
    start_context: FormatContext | None = None,
) -> str:
    """Format TRUST .data text.

    Normalizes indentation based on brace depth, collapses consecutive
    blank lines to at most one, and drops blank lines immediately after
    `{` or before `}`. Tab size / use-spaces are locked; see module
    docstring.

    Trailing whitespace is dropped on every line — including lines that
    open a multi-line comment and lines BETWEEN the comment delimiters.
    This is intentional (audit 6.3): trailing whitespace carries no
    syntactic meaning in TRUST, every standard editor / linter strips
    it, and applying the rule uniformly avoids surprising the user with
    different rules inside vs. outside comment blocks. Leading
    whitespace inside `/* ... */` blocks IS preserved verbatim so
    commented-out source the user might uncomment stays readable
    (line 188).

    EOL normalisation (audit 6.4): in full-file mode (no
    `start_context` argument) the output always ends with exactly one
    `\\n` when the input has any content (none added to a truly empty
    input). CRLF inputs collapse to LF because `.strip()`/`.rstrip()`
    on every line removes the trailing `\\r` left by `split("\\n")`.
    git and POSIX cat treat newline-less files as malformed, so
    enforcing the trailing newline keeps formatted output well-formed
    everywhere and idempotent under repeated formatting. Range mode
    (any `start_context`) SKIPS the trailing-newline rule because the
    snippet may be embedded inside a larger document — appending `\\n`
    would corrupt the surrounding text.
    """
    full_file_mode = start_context is None
    if start_context is None:
        start_context = FormatContext(level=0, in_block_comment=False, comment_mode=None)

    indent_str = " " * _TAB_SIZE if _INSERT_SPACES else "\t"
    lines = text.split("\n")
    expanded: list[str] = []
    level = start_context.level
    state: CommentMode = _resolve_mode(start_context)
    pending_blank = False
    prev_opened_block = False

    for raw_line in lines:
        line = raw_line.strip()

        if line == "":
            pending_blank = True
            continue

        line_entered_in_comment = state is not None

        starts_with_close = (not line_entered_in_comment) and line.startswith("}")
        if pending_blank and not prev_opened_block and not starts_with_close:
            expanded.append("")
        pending_blank = False

        if line_entered_in_comment:
            # Inside a multi-line comment we keep the user's indentation
            # verbatim — commented-out source may contain valid syntax the
            # user wants to be able to read or uncomment unchanged (e.g.
            # the `# BEGIN PARTITION ... END PARTITION #` mechanism).
            expanded.append(raw_line.rstrip())
            _, state = scan_line(line, state)
            prev_opened_block = False
            continue

        info, state = _analyze_line_with_state(line, state)

        level = max(0, level - info.leading_closes)
        expanded.append(indent_str * level + line)
        level = max(0, level + info.opens - info.closes + info.leading_closes)

        prev_opened_block = not info.ends_in_comment and line.endswith("{")

    if pending_blank and not prev_opened_block:
        expanded.append("")

    result = "\n".join(expanded)
    if full_file_mode and result and not result.endswith("\n"):
        result += "\n"
    return result


class FoldingRegion(NamedTuple):
    start_line: int
    end_line: int
    is_comment: bool


def find_folding_regions(text: str) -> list[FoldingRegion]:
    """Return foldable regions for TRUST .data text.

    Emits one region per matching `{`/`}` pair that spans more than one
    line, one per multi-line `# ... #` comment, and one per multi-line
    `/* ... */` comment. Single-line pairs are skipped; an unmatched brace
    is silently dropped. Regions are returned sorted by start line.
    """
    regions: list[FoldingRegion] = []
    brace_stack: list[int] = []
    state: CommentMode = None
    comment_start_line = 0

    lines = text.split("\n")
    for line_num, raw_line in enumerate(lines):
        entered_in_comment = state is not None
        comment_ranges, state_after = scan_line(raw_line, state)

        # Detect transitions of the comment state on this line. Multiple
        # ranges may close and re-open on the same line — we only care
        # about the *first* open if not already inside one (start of
        # range) and the *last* close if we end outside any comment.
        if not entered_in_comment and comment_ranges:
            # Comment opened on this line (no prior open).
            comment_start_line = line_num
        # The active multi-line comment closed somewhere on this line.
        if entered_in_comment and state_after is None and comment_start_line < line_num:
            regions.append(FoldingRegion(comment_start_line, line_num, True))

        masked = set()
        for s, e in comment_ranges:
            masked.update(range(s, e))

        for col, ch in enumerate(raw_line):
            if col in masked:
                continue
            if ch == "{":
                brace_stack.append(line_num)
            elif ch == "}" and brace_stack:
                start = brace_stack.pop()
                if start < line_num:
                    regions.append(FoldingRegion(start, line_num, False))

        # If a fresh comment opened *and* did not close on the same line,
        # remember its start for the eventual close detection on a later
        # line.
        if not entered_in_comment and state_after is not None:
            comment_start_line = line_num

        state = state_after

    regions.sort(key=lambda r: (r.start_line, r.end_line))
    return regions


def format_on_type(
    text: str,
    line_num: int,
) -> tuple[str, int] | None:
    """Re-indent a single line after a character is typed.

    Tab size / use-spaces are locked; see module docstring.

    Only acts when the target line starts with one or more `}` — otherwise
    returns None. Preserves any content after the leading `}` sequence.
    Returns None when no re-indent is needed (line already correct, line
    inside a block comment, or line out of range); otherwise returns
    `(new_line, original_line_length)`.
    """
    lines = text.split("\n")
    if line_num < 0 or line_num >= len(lines):
        return None

    target = lines[line_num]
    trimmed = target.strip()
    if not trimmed:
        return None

    info = analyze_line(trimmed)
    if info.leading_closes == 0:
        return None

    before = "\n".join(lines[:line_num]) + "\n" if line_num > 0 else ""
    ctx = compute_format_context(before)
    if ctx.in_block_comment:
        return None

    new_level = max(0, ctx.level - info.leading_closes)
    indent_str = " " * _TAB_SIZE if _INSERT_SPACES else "\t"
    new_line = indent_str * new_level + trimmed
    if new_line == target:
        return None
    return (new_line, len(target))
