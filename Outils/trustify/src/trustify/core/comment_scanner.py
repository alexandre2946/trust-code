"""Comment-state scanner for TRUST `.data` files.

Comment recogniser for the **formatter and the LSP** (trustify-lsp:
diagnostics, semantic-token highlighting, document-symbol extraction,
etc.) — anywhere a callsite needs to know whether a character position
sits inside a `# ... #` or `/* ... */` comment.

Divergence from the runtime parser (read this before assuming this
module is THE comment authority):

  This scanner is CHARACTER-based — any `#` opens a hash comment and any
  `/*` opens a block comment, regardless of surrounding whitespace, so
  `val#x#` is seen as containing a comment.

  The actual dataset parser (`core/trust_parser.py::TRUSTParser.tokenize`)
  is TOKEN-based: it whitespace-splits first, then only a token that is
  *exactly* `#` / `/*` / `*/` toggles comment state. So today a glued
  `#end#` is a single DATA token to the parser, not a comment — the two
  recognisers disagree on glued delimiters.

  This is a known, deliberately one-sided divergence, and it resolves
  from the TRUST side: TRUST's own reader is moving to accept glued forms
  like `#end#` as valid comments, and the trustify parser will follow
  once it does. This scanner already implements that forward-looking
  char-based behaviour; when the parser converges on it the two will
  agree and this note can be deleted. No fix is wanted on the trustify
  side in the meantime.

  Not a problem in practice: the LSP — this scanner's main consumer —
  already emits a diagnostic when comment delimiters are not
  space-separated, steering authors to the space-delimited form that
  BOTH recognisers agree on today.

Primitive: `scan_line(line, in_comment=None, up_to=None)` returns the
comment ranges on the line plus the comment mode at the scan boundary.
Callers thread the returned state into the next call to scan across
lines.

Helpers: `comment_state_at_line(lines, target_line)`,
`is_position_in_comment(lines, line, col)`. Both take a PRE-SPLIT
`list[str]` rather than the raw text (audit 6.5). The intended
consumer is the LSP, which probes "is line N in a comment?" for
many lines per document version; accepting raw text would invite an
O(file²) cost as every call re-split. Callers split once at the
document/text boundary (typically alongside the LSP's per-version
document cache) and reuse the list across helpers.
"""

from __future__ import annotations

from typing import Literal

CommentMode = Literal["hash", "slash"] | None


def scan_line(
    line: str,
    in_comment: CommentMode = None,
    up_to: int | None = None,
) -> tuple[list[tuple[int, int]], CommentMode]:
    """Scan ``line[:up_to]`` and return ``(comment_ranges, end_state)``.

    ``in_comment`` is the comment mode entering the line: ``None`` if not
    inside a comment, ``"hash"`` or ``"slash"`` otherwise. The returned
    ``end_state`` is the mode at the scan boundary (column ``up_to`` or
    end of line); thread it into the next ``scan_line`` call to scan the
    next line.

    ``up_to`` is an exclusive column index. ``None`` means scan the whole
    line. A column past end-of-line is treated as end-of-line.

    ``comment_ranges`` is a list of ``(start, end_exclusive)`` half-open
    intervals. Each interval includes the marker characters (``#`` or
    ``/*`` / ``*/``).
    """
    end = min(up_to, len(line)) if up_to is not None else len(line)
    if end <= 0:
        # ``up_to=0`` means scan nothing — even an entering ``in_comment``
        # state collapses to an empty range because the range would be
        # ``(0, 0)``.
        return [], None if up_to == 0 else in_comment
    ranges: list[tuple[int, int]] = []
    comment_start = 0 if in_comment is not None else -1
    col = 0
    while col < end:
        ch = line[col]
        if in_comment == "hash":
            if ch == "#":
                ranges.append((comment_start, col + 1))
                in_comment = None
                comment_start = -1
            col += 1
            continue
        if in_comment == "slash":
            if ch == "*" and col + 1 < end and line[col + 1] == "/":
                ranges.append((comment_start, col + 2))
                in_comment = None
                comment_start = -1
                col += 2
            else:
                col += 1
            continue
        if ch == "#":
            comment_start = col
            in_comment = "hash"
            col += 1
            continue
        if ch == "/" and col + 1 < end and line[col + 1] == "*":
            comment_start = col
            in_comment = "slash"
            col += 2
            continue
        col += 1
    if in_comment is not None:
        ranges.append((comment_start, end))
    return ranges, in_comment


def comment_state_at_line(lines: list[str], target_line: int) -> CommentMode:
    """Return the comment state at the start of ``lines[target_line]``
    (0-indexed). Reconstructed by scanning every prior line; ``None``
    if ``target_line`` is 0 or negative.

    ``lines`` is the pre-split document — typically
    ``text.split("\\n")`` cached once per document version. See module
    docstring (audit 6.5) for why the helper does NOT accept raw text.
    """
    if target_line <= 0:
        return None
    state: CommentMode = None
    upper = min(target_line, len(lines))
    for i in range(upper):
        _, state = scan_line(lines[i], state)
    return state


def is_position_in_comment(lines: list[str], line: int, col: int) -> bool:
    """Whether ``lines[line][col]`` is strictly inside a comment.

    ``col`` is half-open: ``col == 0`` on a line not entering in a
    comment is not in any comment; the opening marker character itself
    is considered inside the comment range it opens.

    ``lines`` is the pre-split document — see module docstring
    (audit 6.5).
    """
    state = comment_state_at_line(lines, line)
    if line < 0 or line >= len(lines):
        return state is not None
    ranges, _ = scan_line(lines[line], state)
    return any(s <= col < e for s, e in ranges)
