#!/usr/bin/env python3
"""
Refactor legacy managed-pointer checks to operator bool:

    ptr.est_nul()   ->  !ptr
    ptr.non_nul()   ->   ptr

Handles mixed cases like !ptr.est_nul() (becomes `ptr`, collapsing the
double negation) and !ptr.non_nul() (becomes `!ptr`).

The "receiver" of the method call is recovered by walking backwards from
the `.est_nul()`/`.non_nul()` token through a postfix-expression chain:
identifiers, member accesses (`.`, `->`, `::`), function calls `(...)`,
and subscripts `[...]`, respecting nesting of `()`, `[]`, `{}` and
skipping over string/char literals and comments.

String literals, character literals, and comments (// and /* */) are
never modified.

Usage:
    python3 refactor_ptr.py [--dry-run] [--write] [--ext .cpp,.h,...] PATH [PATH ...]

By default the script prints a unified diff for each file. Pass --write
to actually modify files in place.
"""

from __future__ import annotations

import argparse
import difflib
import os
import re
import sys
from pathlib import Path
from typing import Iterable

# --------------------------------------------------------------------------- #
# Tokenizer-ish scan                                                           #
# --------------------------------------------------------------------------- #

# We don't need a full C++ lexer. We need to know, at any index in the source,
# whether we're inside a string/char literal or a comment, so we skip those.
# We also need to find `.est_nul()` / `.non_nul()` outside those regions.

METHOD_RE = re.compile(r'\.\s*(est_nul|non_nul)\s*\(\s*\)')


def _skippable_regions(src: str) -> list[tuple[int, int]]:
    """Return sorted list of (start, end) half-open intervals covering every
    string literal, char literal, and comment in `src`. Matches inside these
    regions must be ignored."""
    regions: list[tuple[int, int]] = []
    i = 0
    n = len(src)
    while i < n:
        c = src[i]
        nxt = src[i + 1] if i + 1 < n else ''

        # Line comment
        if c == '/' and nxt == '/':
            j = src.find('\n', i + 2)
            if j == -1:
                j = n
            regions.append((i, j))
            i = j
            continue

        # Block comment
        if c == '/' and nxt == '*':
            j = src.find('*/', i + 2)
            j = n if j == -1 else j + 2
            regions.append((i, j))
            i = j
            continue

        # Raw string literal: R"delim(...)delim"
        # Detect R" possibly prefixed by u8/u/U/L.
        if c == 'R' and nxt == '"' and _is_raw_string_start(src, i):
            end = _scan_raw_string(src, i)
            regions.append((i, end))
            i = end
            continue
        if c in 'uUL' and _is_raw_string_start_with_prefix(src, i):
            # find the R"
            r_idx = src.index('R', i, i + 4)
            end = _scan_raw_string(src, r_idx)
            regions.append((i, end))
            i = end
            continue

        # Regular string literal
        if c == '"':
            j = _scan_quoted(src, i, '"')
            regions.append((i, j))
            i = j
            continue

        # Char literal
        if c == "'":
            j = _scan_quoted(src, i, "'")
            regions.append((i, j))
            i = j
            continue

        i += 1

    return regions


def _is_raw_string_start(src: str, i: int) -> bool:
    # src[i] == 'R', src[i+1] == '"'. Ensure 'R' isn't part of an identifier.
    if i > 0 and (src[i - 1].isalnum() or src[i - 1] == '_'):
        return False
    return True


def _is_raw_string_start_with_prefix(src: str, i: int) -> bool:
    # Check for patterns like u8R"...", uR"...", UR"...", LR"..."
    # src[i] is in 'uUL'.
    rest = src[i:i + 4]
    if rest.startswith('u8R"'):
        return i == 0 or not (src[i - 1].isalnum() or src[i - 1] == '_')
    if len(rest) >= 3 and rest[0] in 'uUL' and rest[1] == 'R' and rest[2] == '"':
        return i == 0 or not (src[i - 1].isalnum() or src[i - 1] == '_')
    return False


def _scan_raw_string(src: str, i: int) -> int:
    # src[i] == 'R', src[i+1] == '"'. Parse delimiter up to '(', then find )delim"
    assert src[i] == 'R' and src[i + 1] == '"'
    delim_start = i + 2
    paren = src.find('(', delim_start)
    if paren == -1:
        return len(src)
    delim = src[delim_start:paren]
    terminator = ')' + delim + '"'
    end = src.find(terminator, paren + 1)
    if end == -1:
        return len(src)
    return end + len(terminator)


def _scan_quoted(src: str, i: int, quote: str) -> int:
    """Scan a regular (non-raw) quoted literal starting at `i`, return index
    just past the closing quote. Handles backslash escapes."""
    j = i + 1
    n = len(src)
    while j < n:
        c = src[j]
        if c == '\\':
            j += 2
            continue
        if c == quote:
            return j + 1
        if c == '\n':
            # unterminated; bail out
            return j
        j += 1
    return n


def _in_regions(pos: int, regions: list[tuple[int, int]]) -> bool:
    # Linear scan is fine; regions are typically few per file and sorted.
    # For very large files we could bisect, but this is not a hot path.
    for start, end in regions:
        if start <= pos < end:
            return True
        if start > pos:
            return False
    return False


# --------------------------------------------------------------------------- #
# Receiver extraction                                                          #
# --------------------------------------------------------------------------- #

# The receiver ends just before the '.' of `.est_nul()`. We parse it as a
# postfix-expression chain walked backwards:
#
#     atom  ( connector  atom )*
#
# where:
#   - atom       ::= identifier | bracketed_group | template_arg_group
#   - connector  ::= '.' | '->' | '::'
#   - bracketed_group matches `(...)`, `[...]`, `{...}`
#
# The backward walker alternates: starting in state "expecting atom", it
# consumes one atom, then looks left for a connector; if found, consume it
# and loop; otherwise stop. Whitespace between tokens is fine, but we
# never cross a connector/atom boundary onto an unrelated token.

_IDENT_CHAR = re.compile(r'[A-Za-z0-9_]')
_IDENT_START = re.compile(r'[A-Za-z_]')


def _skip_ws_back(src: str, i: int) -> int:
    """Return the index of the first non-whitespace char at or before i."""
    while i >= 0 and src[i] in ' \t\r\n':
        i -= 1
    return i


def _consume_atom_back(src: str, i: int,
                       regions: list[tuple[int, int]]) -> int | None:
    """Try to consume one atom ending at index i (inclusive). Return the
    index just before the atom on success, or None if no atom found."""
    if i < 0:
        return None
    if _in_regions(i, regions):
        return None

    c = src[i]

    # Bracketed group: ')', ']', '}'
    if c in ')]}':
        opener = _match_bracket_backwards(src, i, regions)
        if opener < 0:
            return None
        return opener - 1

    # Template arg group: '>' not preceded by '-' (that would be '->').
    if c == '>' and not (i >= 1 and src[i - 1] == '-'):
        opener = _match_angle_backwards(src, i, regions)
        if opener is not None:
            return opener - 1
        return None

    # Identifier.
    if _IDENT_CHAR.match(c):
        j = i
        while j >= 0 and _IDENT_CHAR.match(src[j]):
            j -= 1
        # Must start with a letter or underscore, not a digit, otherwise
        # it's a numeric literal which isn't a valid receiver.
        first = j + 1
        if first <= i and _IDENT_START.match(src[first]):
            return j
        return None

    return None


def _consume_connector_back(src: str, i: int) -> int | None:
    """Try to consume a member/scope connector ending at index i.
    Returns the index just before the connector on success, else None."""
    if i < 0:
        return None
    c = src[i]
    # '->'
    if c == '>' and i >= 1 and src[i - 1] == '-':
        return i - 2
    # '::'
    if c == ':' and i >= 1 and src[i - 1] == ':':
        return i - 2
    # '.'  (but NOT '..' which would be weird, and NOT '->' — already
    #  handled — and the '.' is not part of a number because numeric
    #  literals can't be receivers of a method call anyway.)
    if c == '.':
        return i - 1
    return None


def _match_angle_backwards(src: str, close_pos: int,
                           regions: list[tuple[int, int]]) -> int | None:
    """Heuristically match '<...>' as a template-argument group, scanning
    backwards from a '>' at `close_pos`. Returns the index of the matching
    '<', or None if this doesn't look like template args.

    We scan backwards tracking nested '<>', '()', '[]', '{}'. We bail out
    (returning None) if we see top-level tokens that shouldn't appear in
    template args: ';', '{' at depth 0, or the end of file.
    """
    depth_angle = 1
    depth_paren = 0
    depth_brack = 0
    depth_brace = 0
    i = close_pos - 1
    # Bound the search: template args rarely exceed a few hundred chars.
    limit = max(0, close_pos - 500)
    while i >= limit:
        if _in_regions(i, regions):
            region = next(r for r in regions if r[0] <= i < r[1])
            i = region[0] - 1
            continue
        c = src[i]
        if c == '>' and not (i >= 1 and src[i - 1] == '-'):
            if depth_paren == depth_brack == depth_brace == 0:
                depth_angle += 1
        elif c == '<':
            if depth_paren == depth_brack == depth_brace == 0:
                depth_angle -= 1
                if depth_angle == 0:
                    # Verify the '<' is preceded (after whitespace) by an
                    # identifier char — otherwise it's a less-than, not a
                    # template bracket.
                    k = _skip_ws_back(src, i - 1)
                    if k >= 0 and _IDENT_CHAR.match(src[k]):
                        return i
                    return None
        elif c == ')':
            depth_paren += 1
        elif c == '(':
            if depth_paren == 0:
                return None
            depth_paren -= 1
        elif c == ']':
            depth_brack += 1
        elif c == '[':
            if depth_brack == 0:
                return None
            depth_brack -= 1
        elif c == '}':
            depth_brace += 1
        elif c == '{':
            if depth_brace == 0:
                return None
            depth_brace -= 1
        elif c == ';':
            if depth_paren == depth_brack == depth_brace == 0:
                return None
        i -= 1
    return None


def _find_receiver_start(src: str, dot_pos: int,
                         regions: list[tuple[int, int]]) -> int:
    """Given the index of the '.' in `.est_nul()`, return the index where
    the receiver expression starts.

    Walks backwards, alternating atoms and connectors:

        atom := identifier | '(' ... ')' | '[' ... ']' | '{' ... '}' | '<' ... '>'
        connector := '.' | '->' | '::'

    Bracket/template atoms attach directly to the atom on their left
    (no connector needed — they're call/subscript/template applications).
    Identifier atoms must be joined to the previous atom (if any) by a
    connector.
    """
    i = _skip_ws_back(src, dot_pos - 1)
    if i < 0:
        return dot_pos

    receiver_first = dot_pos  # no atom consumed yet

    def consume_atom(idx: int) -> int | None:
        """Try to consume ONE atom ending at `idx` (inclusive). On success,
        update receiver_first to point to the atom's start and return the
        index just before the atom. On failure, return None."""
        nonlocal receiver_first
        if idx < 0 or _in_regions(idx, regions):
            return None
        ch = src[idx]

        # Bracket group.
        if ch in ')]}':
            opener = _match_bracket_backwards(src, idx, regions)
            if opener < 0:
                return None
            receiver_first = opener
            return opener - 1

        # Template args.
        if ch == '>' and not (idx >= 1 and src[idx - 1] == '-'):
            opener = _match_angle_backwards(src, idx, regions)
            if opener is None:
                return None
            receiver_first = opener
            return opener - 1

        # Identifier.
        if _IDENT_CHAR.match(ch):
            j = idx
            while j >= 0 and _IDENT_CHAR.match(src[j]):
                j -= 1
            ident_start = j + 1
            if not _IDENT_START.match(src[ident_start]):
                return None  # numeric literal
            receiver_first = ident_start
            return j

        return None

    # Consume the first atom (required).
    after = consume_atom(i)
    if after is None:
        return dot_pos

    # Keep consuming: either another bracket/template atom attached
    # directly, or a connector-then-identifier pair, or an identifier
    # that serves as the callee/indexee of the just-consumed bracket.
    last_was_bracket = src[i] in ')]}' or (src[i] == '>' and not (i >= 1 and src[i - 1] == '-'))
    while True:
        k = _skip_ws_back(src, after)
        if k < 0:
            break

        ch = src[k]

        # Directly-attached bracket/template atom.
        if ch in ')]}' or (ch == '>' and not (k >= 1 and src[k - 1] == '-')):
            nxt = consume_atom(k)
            if nxt is None:
                break
            after = nxt
            last_was_bracket = True
            continue

        # If the last atom consumed was a bracket/template, an identifier
        # directly to its left (no connector) is its callee/indexee.
        if last_was_bracket and _IDENT_CHAR.match(ch) and not _in_regions(k, regions):
            nxt = consume_atom(k)
            if nxt is None:
                break
            after = nxt
            last_was_bracket = False
            continue

        # Connector + atom (a mandatory pair). The atom after a connector
        # can be an identifier or a bracket atom (e.g. `foo->bar[i].est_nul()`
        # — walking back from the dot we cross `->` then hit `]`).
        before_conn = _consume_connector_back(src, k)
        if before_conn is None:
            break
        a = _skip_ws_back(src, before_conn)
        if a < 0:
            break
        if _in_regions(a, regions):
            break
        nxt = consume_atom(a)
        if nxt is None:
            break
        after = nxt
        last_was_bracket = src[a] in ')]}' or (src[a] == '>' and not (a >= 1 and src[a - 1] == '-'))

    return receiver_first


def _match_bracket_backwards(src: str, close_pos: int,
                             regions: list[tuple[int, int]]) -> int:
    """Given a closing bracket position, return the index of its matching
    opener, or -1 if unbalanced."""
    pairs = {')': '(', ']': '[', '}': '{'}
    closer = src[close_pos]
    opener = pairs[closer]
    depth = 1
    i = close_pos - 1
    while i >= 0:
        if _in_regions(i, regions):
            region = next(r for r in regions if r[0] <= i < r[1])
            i = region[0] - 1
            continue
        c = src[i]
        if c == closer:
            depth += 1
        elif c == opener:
            depth -= 1
            if depth == 0:
                return i
        elif c in ')]}':
            # nested different bracket type closing; recurse
            i = _match_bracket_backwards(src, i, regions)
            if i < 0:
                return -1
            i -= 1
            continue
        i -= 1
    return -1


_BOOL_CTX_KEYWORDS = {'if', 'while', 'for', 'assert', 'static_assert'}


def _is_in_boolean_context(src: str, span_start: int, span_end: int,
                           regions: list[tuple[int, int]]) -> bool:
    """Return True if the replaced expression at [span_start, span_end)
    sits in a position where contextual boolean conversion is guaranteed
    to happen, so an `explicit operator bool` on the receiver will fire
    without an explicit cast.

    Guaranteed-boolean positions we detect:
      - directly followed by `&&`, `||`, or `?` (ternary condition)
      - directly preceded by `&&`, `||`, or `!`
      - enclosed in parentheses that immediately follow `if`, `while`,
        `for` (as the loop condition, after the first `;`), `!`, or one
        of a few other boolean-forcing keywords
    """
    # --- right context ---
    j = span_end
    n = len(src)
    # skip whitespace and comments
    while j < n:
        if _in_regions(j, regions):
            region = next(r for r in regions if r[0] <= j < r[1])
            j = region[1]
            continue
        if src[j] in ' \t\r\n':
            j += 1
            continue
        break
    if j < n:
        c = src[j]
        nxt = src[j + 1] if j + 1 < n else ''
        if (c == '&' and nxt == '&') or (c == '|' and nxt == '|'):
            return True
        if c == '?':
            # Distinguish ternary-`?` from trigraphs or other oddities.
            # A lone `?` after an expression is practically always ternary.
            return True

    # --- left context ---
    i = span_start - 1
    while i >= 0:
        if _in_regions(i, regions):
            region = next(r for r in regions if r[0] <= i < r[1])
            i = region[0] - 1
            continue
        if src[i] in ' \t\r\n':
            i -= 1
            continue
        break
    if i >= 0:
        c = src[i]
        prev = src[i - 1] if i >= 1 else ''
        if (c == '&' and prev == '&') or (c == '|' and prev == '|'):
            return True
        if c == '!' and prev != '!' and prev != '=':
            # A single `!` immediately on our left — note: the absorption
            # logic already ate this in most cases, but an un-absorbed `!`
            # (e.g. following `!!`) still forces boolean context.
            return True

    # --- enclosing parens of `if`/`while`/`for`/`!`/... ---
    # Find the innermost enclosing '(' (at depth 0 of everything else).
    paren_open = _find_enclosing_paren_open(src, span_start, regions)
    if paren_open is not None:
        # Look at what precedes the '(' (skipping whitespace/comments).
        k = paren_open - 1
        while k >= 0:
            if _in_regions(k, regions):
                region = next(r for r in regions if r[0] <= k < r[1])
                k = region[0] - 1
                continue
            if src[k] in ' \t\r\n':
                k -= 1
                continue
            break
        if k >= 0:
            # `!(...)`
            if src[k] == '!' and (k == 0 or src[k - 1] != '!'):
                return True
            # Keyword before `(`.
            if _IDENT_CHAR.match(src[k]):
                end = k + 1
                kw_start = k
                while kw_start >= 0 and _IDENT_CHAR.match(src[kw_start]):
                    kw_start -= 1
                kw = src[kw_start + 1:end]
                if kw in _BOOL_CTX_KEYWORDS:
                    # For `for (init; cond; step)`, only the middle slot is
                    # boolean. Detect by checking whether our position is
                    # after the first `;` inside the parens.
                    if kw == 'for':
                        return _is_in_for_cond_slot(src, paren_open,
                                                    span_start, regions)
                    return True

    return False


def _find_enclosing_paren_open(src: str, pos: int,
                               regions: list[tuple[int, int]]) -> int | None:
    """Return the index of the '(' that directly encloses `pos`, or None
    if `pos` is not inside parens within the current block.

    Note: `;` is NOT treated as a hard boundary because `for (init; cond;
    step)` uses `;` as an intra-paren separator. We rely on `{`/`}` and
    start-of-file to bound the search.
    """
    depth_paren = 0
    depth_brack = 0
    depth_brace = 0
    i = pos - 1
    while i >= 0:
        if _in_regions(i, regions):
            region = next(r for r in regions if r[0] <= i < r[1])
            i = region[0] - 1
            continue
        c = src[i]
        if c == ')':
            depth_paren += 1
        elif c == '(':
            if depth_paren == 0 and depth_brack == 0 and depth_brace == 0:
                return i
            depth_paren -= 1
        elif c == ']':
            depth_brack += 1
        elif c == '[':
            depth_brack -= 1
        elif c == '}':
            depth_brace += 1
        elif c == '{':
            if depth_brace == 0:
                return None
            depth_brace -= 1
        i -= 1
    return None


def _is_in_for_cond_slot(src: str, paren_open: int, pos: int,
                         regions: list[tuple[int, int]]) -> bool:
    """Return True if `pos` is inside the middle slot (the condition) of a
    C-style `for (init; cond; step)` header."""
    # Count top-level `;` between paren_open and pos.
    depth_paren = 0
    depth_brack = 0
    depth_brace = 0
    semis = 0
    i = paren_open + 1
    while i < pos:
        if _in_regions(i, regions):
            region = next(r for r in regions if r[0] <= i < r[1])
            i = region[1]
            continue
        c = src[i]
        if c == '(':
            depth_paren += 1
        elif c == ')':
            depth_paren -= 1
        elif c == '[':
            depth_brack += 1
        elif c == ']':
            depth_brack -= 1
        elif c == '{':
            depth_brace += 1
        elif c == '}':
            depth_brace -= 1
        elif c == ';' and depth_paren == 0 and depth_brack == 0 and depth_brace == 0:
            semis += 1
        i += 1
    return semis == 1


# --------------------------------------------------------------------------- #
# Transformation                                                               #
# --------------------------------------------------------------------------- #

def refactor(src: str, *, wrap_bool: bool = True) -> tuple[str, int]:
    """Return (new_src, total_num_replacements). Runs to a fixed point so
    that nested occurrences — e.g. ``foo(x.non_nul()).est_nul()`` — are
    all rewritten.

    When `wrap_bool` is True (the default), a `non_nul()` rewrite that
    would produce a bare receiver in a non-boolean context is wrapped in
    ``bool(...)`` to support `explicit operator bool`. When False, the
    rewrite always emits the bare receiver (correct only if the pointer
    type has an implicit `operator bool`)."""
    total = 0
    for _ in range(100):
        src, n = _refactor_once(src, wrap_bool=wrap_bool)
        if n == 0:
            break
        total += n
    return src, total


def _refactor_once(src: str, *, wrap_bool: bool) -> tuple[str, int]:
    """Single pass: rewrite only occurrences whose receivers contain no
    further `.est_nul()` / `.non_nul()` call. Returns (new_src, count)."""
    regions = _skippable_regions(src)

    replacements: list[tuple[int, int, str]] = []

    for m in METHOD_RE.finditer(src):
        dot_pos = m.start()
        if _in_regions(dot_pos, regions):
            continue

        method = m.group(1)
        call_end = m.end()

        recv_start = _find_receiver_start(src, dot_pos, regions)
        if recv_start >= dot_pos:
            continue
        receiver = src[recv_start:dot_pos].rstrip()

        if METHOD_RE.search(receiver):
            continue

        neg_start, has_neg = _find_absorbable_negation(src, recv_start, regions)

        if has_neg:
            span_start = neg_start
            if method == 'est_nul':
                new_text = receiver            # !!x -> x  (still boolean: the outer ! of !!x is gone, but this was
                                               # `!receiver.est_nul()` i.e. `!(!receiver)` so result `receiver`
                                               # is in the SAME context the original expression was in.)
                needs_wrap = (wrap_bool
                              and not _is_in_boolean_context(
                                  src, span_start, call_end, regions))
                if needs_wrap:
                    new_text = f'bool({receiver})'
            else:
                new_text = '!' + receiver      # `!` forces conversion; no wrap needed
            span_end = call_end
        else:
            span_start = recv_start
            if method == 'est_nul':
                new_text = '!' + receiver      # `!` forces conversion
            else:
                needs_wrap = (wrap_bool
                              and not _is_in_boolean_context(
                                  src, span_start, call_end, regions))
                new_text = f'bool({receiver})' if needs_wrap else receiver
            span_end = call_end

        replacements.append((span_start, span_end, new_text))

    if not replacements:
        return src, 0

    # Apply from the end so earlier indices stay valid.
    out_parts: list[str] = []
    last = len(src)
    for start, end, text in reversed(replacements):
        out_parts.append(src[end:last])
        out_parts.append(text)
        last = start
    out_parts.append(src[:last])
    return ''.join(reversed(out_parts)), len(replacements)


def _find_absorbable_negation(src: str, recv_start: int,
                              regions: list[tuple[int, int]]
                              ) -> tuple[int, bool]:
    """Look left of `recv_start` for a leading '!' that can be folded into
    the transformation. Return (position_of_bang, True) on success,
    else (recv_start, False).

    We require the '!' to be a plain logical-NOT: not part of '!=', and not
    itself preceded by another '!' (we don't want to eat just one of '!!').
    We only absorb one '!'; an odd number of '!'s still needs one flip.
    """
    i = recv_start - 1
    while i >= 0 and src[i] in ' \t\r\n':
        i -= 1
    if i < 0 or src[i] != '!':
        return recv_start, False
    if _in_regions(i, regions):
        return recv_start, False
    # Guard: '!=' — shouldn't happen since '=' would be to the right, but check.
    if i + 1 < len(src) and src[i + 1] == '=':
        return recv_start, False
    # Guard: don't absorb a '!' that is itself the right half of '!!', because
    # doing so would leave behind a stray '!'. Simplest: only absorb if the
    # char immediately before (skipping ws) is not '!'.
    k = i - 1
    while k >= 0 and src[k] in ' \t\r\n':
        k -= 1
    if k >= 0 and src[k] == '!':
        return recv_start, False
    return i, True


# --------------------------------------------------------------------------- #
# CLI                                                                          #
# --------------------------------------------------------------------------- #

DEFAULT_EXTS = ('.c', '.cc', '.cpp', '.cxx', '.h', '.hh', '.hpp', '.hxx',
                '.inl', '.ipp', '.tcc')


def iter_files(paths: Iterable[str], exts: tuple[str, ...]) -> Iterable[Path]:
    for p in paths:
        path = Path(p)
        if path.is_file():
            yield path
        elif path.is_dir():
            for f in path.rglob('*'):
                if f.is_file() and f.suffix in exts:
                    yield f


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='+', help='Files or directories to process.')
    ap.add_argument('--write', action='store_true',
                    help='Modify files in place. Without this flag, prints '
                         'a unified diff and leaves files untouched.')
    ap.add_argument('--ext', default=','.join(DEFAULT_EXTS),
                    help='Comma-separated list of file extensions to process '
                         'when a directory is given. Default: %(default)s')
    ap.add_argument('--no-wrap-bool', dest='wrap_bool', action='store_false',
                    help='Do not wrap `non_nul()` rewrites in `bool(...)` '
                         'when the context is not contextually boolean. '
                         'Use this only if your target pointer type has an '
                         'implicit (non-explicit) `operator bool`.')
    ap.set_defaults(wrap_bool=True)
    args = ap.parse_args(argv)

    exts = tuple(e if e.startswith('.') else '.' + e
                 for e in args.ext.split(',') if e.strip())

    total_files = 0
    changed_files = 0
    total_replacements = 0

    for path in iter_files(args.paths, exts):
        total_files += 1
        try:
            original = path.read_text(encoding='utf-8')
        except UnicodeDecodeError:
            print(f'[skip] {path}: not UTF-8', file=sys.stderr)
            continue

        new, n = refactor(original, wrap_bool=args.wrap_bool)
        if n == 0:
            continue

        changed_files += 1
        total_replacements += n

        if args.write:
            path.write_text(new, encoding='utf-8')
            print(f'[edit] {path}: {n} replacement(s)')
        else:
            diff = difflib.unified_diff(
                original.splitlines(keepends=True),
                new.splitlines(keepends=True),
                fromfile=str(path), tofile=str(path),
            )
            sys.stdout.writelines(diff)

    mode = 'wrote' if args.write else 'would change'
    print(f'\n{mode} {changed_files}/{total_files} file(s), '
          f'{total_replacements} replacement(s)', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
