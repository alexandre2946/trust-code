#!/usr/bin/env python

"""TRUST dataset token stream.

`TRUSTParser` whitespace-splits a `.data` text into the
position-preserving `tabToken` list and the lowercase, comment-blanked
`tabTokenLow` list. `TRUSTStream` wraps the result and exposes the
probe / validate / save / restore token-cursor API consumed by the
generated parser classes (see `core/base.py`).
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class SourceRange:
    """LSP-style range. (start_line, start_char) is the first position
    inside the slice; (end_line, end_char) is one past the last position —
    matching lsprotocol.types.Range semantics. All values 0-based.
    """

    start_line: int
    start_char: int
    end_line: int
    end_char: int


def _advance_position(line: int, col: int, body: str) -> tuple[int, int]:
    """Return the 0-based ``(line, col)`` one past the last character of
    ``body``, given that ``body`` (already stripped of leading whitespace)
    starts at ``(line, col)``.

    Walks embedded newlines: ``_mergeQuoted`` can splice a quoted string
    across line breaks into a single token whose body contains ``\\n``, so
    the end may sit on a later line than the start. Mirrors the per-char
    position walk in ``TRUSTParser.tokenize``. With no ``\\n`` this is just
    ``(line, col + len(body))``.
    """
    for ch in body:
        if ch == "\n":
            line += 1
            col = 0
        else:
            col += 1
    return line, col


class TRUSTEndOfStreamException(Exception):
    """Raised when the parser reaches the end of the stream, and ones tries to consume
    more tokens."""

    def __init__(self, msg="Unexpected end of file (or stream)!!"):
        Exception.__init__(self, msg)


class TRUSTTokens:
    """Handy class representing a slice of a TRUSTStream (see below).
    Allows to synchronously deal with lower case and original tokens.
    """

    def __init__(self, low=None, orig=None, range=None):
        if orig is None:
            orig = []
        if low is None:
            low = []
        self._orig = orig[:]
        self._low = low[:]
        self.range = range  # SourceRange | None — None for synthetic / empty slices

    def orig(self):
        return self._orig

    def low(self):
        return self._low

    @classmethod
    def Join(cls, lst):
        """Merge a list of TRUSTTokens objects into a single TRUSTTokens object.
        The merged range spans from the first input's start to the last input's end.
        If any input has range=None or the input list is empty, the merged range
        is None — partial information is not synthesized."""
        l_low, l_orig, ranges = [], [], []
        for ttk in lst:
            l_low.extend(ttk._low)
            l_orig.extend(ttk._orig)
            ranges.append(ttk.range)
        if not ranges or any(r is None for r in ranges):
            joined = None
        else:
            joined = SourceRange(ranges[0].start_line, ranges[0].start_char, ranges[-1].end_line, ranges[-1].end_char)
        return TRUSTTokens(l_low, l_orig, range=joined)


class TRUSTStream:
    """Handy class for scanning the tokens once the parsing has been done.
    Not needed internally but more user friendly for external clients of the parser.
    The parser produces a list of 'original' tokens (as were found in the original dataset)
    and a list of 'lower' tokens, which are stripped, lower-case version of the original tokens, and where the comments
    are compelty removed.
    The two lists are always synchronized (at the same position) and have equal lengths.
    """

    def __init__(self, parser=None, file_nam="??"):
        self.file_nam = file_nam  # For tracking purposes only
        self.tok = []  # Original version, as found in the dataset
        self.tokLow = []  # Lower-case, stripped version, with comments blanked out
        self.lineNum = []  # Line number of the token (for error display)
        self.contentLine = []  # 0-based line where each token's content starts
        self.contentCol = []  # 0-based column where each token's content starts
        if parser is not None:
            self.tok = parser.tabToken[:]  # Better take copies ...
            self.tokLow = parser.tabTokenLow[:]
            self.lineNum = parser.lineNum[:]
            self.contentLine = parser.contentLine[:]
            self.contentCol = parser.contentCol[:]
        self.idx = 0  # Current position in the stream
        # Internal stuff:
        self._probing = False  # whether we are trying to move ahead in the stream, but have not validated the move yet.
        self._probeIdx = 0  # probing position
        self._prevIdx = -1
        self._states = {}  # see  save() / restore() methods

    def clone(self):
        ret = TRUSTStream()
        ret.file_nam = self.file_nam
        ret.tok = self.tok[:]
        ret.tokLow = self.tokLow[:]
        ret.lineNum = self.lineNum[:]
        ret.contentLine = self.contentLine[:]
        ret.contentCol = self.contentCol[:]
        ret.idx = self.idx
        ret._probing = self._probing
        ret._probeIdx = self._probeIdx
        ret._prevIdx = self._prevIdx
        ret._states = dict(self._states)
        return ret

    def __len__(self):
        return len(self.tok)

    def pos(self):
        """Where are we in the stream?"""
        return self.idx

    def setPos(self, i):
        self.idx = i

    def currentTokenIdx(self):
        """Index of the token most relevant for error reporting.

        After ``probeNextLow()`` but before ``validateNext()``, ``self.idx``
        still points at the *previous* validated slot — the slot the probe
        skipped past (over blanks/comments) to find the next real token.
        Error sites raise on the probed token, so they need its slot
        (``_probeIdx - 1``), not the stale ``self.idx``.
        """
        if self._probing:
            return self._probeIdx - 1
        return self.idx

    def currentLine(self):
        """Return the line number of the current token in the data file"""
        idx = self.currentTokenIdx()
        if idx >= len(self.lineNum):
            return self.lineNum[-1]
        return self.lineNum[idx]

    def fileName(self):
        return self.file_nam

    def eof(self):
        """Have we reached end of stream?"""
        # If only blank tokens left, this is EOF:
        try:
            self.probeNextLow()
            self._probing = False
            return False
        except TRUSTEndOfStreamException:
            return True

    def nextLow(self):
        """Returns the next non-blank token (lower case) and advance the stream just **after** it."""
        self._probing = False
        s = self.probeNextLow()
        self.validateNext()
        return s

    def probeNextLow(self):
        """Same as nextLow() but without moving the stream forward. This allows the user to 'peek' into the stream.
        If happy with the result, one can finalize the move with validateNext() below."""
        self._probing = True
        for i, val in enumerate(self.tokLow[self.idx :]):
            if val != "":
                self._probeIdx = self.idx + i + 1
                return val
        self._probing = False
        raise TRUSTEndOfStreamException()

    def validateNext(self):
        """Validate the forward move initiated with probeNextLow()"""
        assert self._probing
        self._prevIdx = self.idx
        self.idx = self._probeIdx
        self._probing = False
        return self.idx

    def lastReadTokens(self):
        """Return the last tokens read tokens when invoking nextLow() or validateNext().
        This returns a TRUSTTokens carrying both the lower-case and original
        slice, plus a SourceRange covering the consumed region (or None if
        nothing has been consumed yet).
        """
        if self._prevIdx < 0:
            return TRUSTTokens()
        # Find the first non-blank token in the consumed slice — that is the
        # content start. (Leading whitespace-only tokens are skipped by
        # probeNextLow, so the actual content starts at the probe hit.)
        first_content_idx = self._prevIdx
        for i in range(self._prevIdx, self.idx):
            if self.tokLow[i] != "":
                first_content_idx = i
                break
        sl, sc = self._zeroBasedStart(first_content_idx)
        last_idx = self.idx - 1
        el, ec = self._zeroBasedEnd(last_idx)
        rng = SourceRange(sl, sc, el, ec)
        return TRUSTTokens(
            low=self.tokLow[self._prevIdx : self.idx], orig=self.tok[self._prevIdx : self.idx], range=rng
        )

    def _zeroBasedStart(self, idx):
        """Return (line, col) of the start of tabToken[idx]'s content
        (after its leading whitespace), 0-based."""
        return self.contentLine[idx], self.contentCol[idx]

    def _zeroBasedEnd(self, idx):
        """Return (line, col) one past the last character of tabToken[idx]'s
        content body, 0-based. A merged multi-line quoted string (see
        _mergeQuoted) can carry embedded '\\n', so walk the body via
        _advance_position rather than assuming the end stays on the start
        line."""
        t = self.tok[idx]
        # Length of the leading whitespace prefix.
        i = 0
        while i < len(t) and t[i] in (" ", "\t", "\n"):
            i += 1
        return _advance_position(self.contentLine[idx], self.contentCol[idx], t[i:])

    def save(self, tag):
        """Save the current state of the stream in a internal dict under key 'tag'. This state can be restored
        using the restore() method. This allows moving forward 'just to try' and then coming back.
        To inspect just one token forward, see 'probeNextLow'
        """
        state = (self.idx, self._probing, self._probeIdx)
        self._states[tag] = state

    def restore(self, tag):
        """See save()."""
        if tag not in self._states:
            raise Exception(f"Can not restore() state, invalid tag '{tag}'!!")
        self.idx, self._probing, self._probeIdx = self._states[tag]

    def dropTail(self, idx):
        """Remove end of stream after index 'idx' (which is itself preserved). Useful for pruning blank stuff at the end of a file for example."""
        self.tok = self.tok[: idx + 1]
        self.tokLow = self.tokLow[: idx + 1]
        self.lineNum = self.lineNum[: idx + 1]
        self.contentLine = self.contentLine[: idx + 1]
        self.contentCol = self.contentCol[: idx + 1]


class TRUSTParser:
    """Main class allowing the parsing of a TRUST dataset."""

    def __init__(self):
        self.tabToken = []  # The list of tokens, with case, tabs and line returns preserved
        self.tabTokenLow = []  # Same as above, but all lowercase, and no tabs, no LR. Comment tokens are also empty.
        self.lineNum = []  # A list with exactly the same size as self.tabToken indicating on which line in the original dataset the token is found
        self.contentLine = []  # 0-based line where each token's content starts (after leading whitespace)
        self.contentCol = []  # 0-based column where each token's content starts (after leading whitespace)

    def _comment_mask(self):
        """Return a bool list — True at index i iff `tabTokenLow[i]` falls
        INSIDE a `# ... #` or `/* ... */` comment. Quotes inside such
        tokens must NOT toggle `_mergeQuoted`'s in-quote state (audit
        1.4 bug 1): TRUST's own tokenizer has zero quote-awareness, so
        comments containing standalone or unbalanced `"` are well-formed
        TRUST data — and `tests/Reference/Decouper_multi_3D/*.data` etc.
        actually use `# BEGIN PARTITION\\nsystem "rm -rf ..."\\n# END
        PARTITION #` idioms that pre-fix would have triggered a phantom
        merge.

        Same state-machine semantics as the post-merge blanking loop
        below (which still runs on `tabTokenLow` to mask comment slots
        for parsing): `#` toggles, `/*` opens, `*/` closes, the two
        styles do NOT nest into each other.
        """
        in_comment: list[bool] = []
        in_hash = False
        star_slash = 0
        for t in self.tabTokenLow:
            # `t` is the closing marker — flag THIS token as inside the
            # comment too (so the closing `#` / `*/` itself is part of
            # the no-quote-toggle zone).
            entering_or_leaving = False
            if t == "#" and star_slash == 0:
                in_hash = not in_hash
                entering_or_leaving = True
            elif t == "/*" and not in_hash:
                star_slash += 1
                entering_or_leaving = True
            elif t == "*/" and not in_hash:
                star_slash = max(0, star_slash - 1)
                entering_or_leaving = True
            in_comment.append(in_hash or star_slash > 0 or entering_or_leaving)
        return in_comment

    def _mergeQuoted(self):
        """Merge back tokens containing quotes so that
        ['ab', 'cd', '"some', ' text', ' in', ' quotes"']
        becomes
        ['ab', 'cd', '"some text in quotes"']

        Tokens inside `# ... #` / `/* ... */` comments are passed
        through unchanged — they do NOT toggle the in-quote state
        (audit 1.4 bug 1). TRUST C++ has no quote-awareness in its
        main tokenizer, so quotes inside comments are well-formed
        TRUST input that trustify must not mis-merge.
        """
        in_comment = self._comment_mask()
        # Merge tokens containing quotes " (used for example in things like 'system "rm -rf toto"')
        in_quot = False
        tt, ttl, lst, lstl = [], [], [], []
        for i, (t, tl) in enumerate(zip(self.tabToken, self.tabTokenLow, strict=False)):
            # Skip quote-toggle bookkeeping on comment-internal tokens.
            if in_comment[i]:
                # Whether in_quot is True or False, a comment-internal
                # token just passes through the current accumulator —
                # but in practice a real .data file never opens a quote
                # outside a comment that stays open into a comment, so
                # in_quot is almost always False here. Preserve the
                # invariant in both cases.
                if in_quot:
                    lst.append(t)
                    lstl.append(tl)
                else:
                    tt.append(t)
                    ttl.append(tl)
                continue
            # If there is an even number of quotes in the token we are already fine, like in '"dom"'
            if tl.count('"') % 2 != 0:
                in_quot = not in_quot
                if not in_quot:  # branch when we just left the quoted mode:
                    tt.append("".join([*lst, t]))
                    ttl.append("".join([*lstl, tl]))
                    lst, lstl = [], []
                    continue
            if in_quot:
                lst.append(t)
                lstl.append(tl)
            else:
                tt.append(t)
                ttl.append(tl)
        self.tabToken, self.tabTokenLow = tt, ttl

    def tokenize(self, txtIn):
        """Split dataset on spaces, tabs and line returns, but keep the spaces, line returns and tabs, in the original tokens
        to be able to reproduce a dataset with them at the end. The two following members are set:
        self.tabToken is the initial stream of tokens, with all spaces, tabs and line returns preserved (if you join it on '' you retrieve
        the initial string)
        self.tabTokenLow is the lower-case, stripped version of self.tabToken
        Once the tokenisation has been performed, a TRUSTStream object can be used to work with the tokens.
        """

        # The idea: split on <space>, then '<tab>', then '<ret>' and put back the separator at the **begining** of the token.
        # (will help in tools like packagepy where we might want to substitute tokens)
        def innerSplit(tab, char):
            ret = []
            for t in tab:
                ts = t.split(char)
                ts2 = [ts[0]] + [char + a for a in ts[1:]] if len(ts) > 1 else ts
                ret.extend(ts2)
            return ret

        a = innerSplit([txtIn], " ")
        a = innerSplit(a, "\t")
        a = innerSplit(a, "\n")
        self.tabToken = a
        self.tabTokenLow = [t.lower().strip() for t in self.tabToken]
        self._mergeQuoted()

        # Extract line numbers (based on '\n'):
        lin = 1
        for t in self.tabToken:
            self.lineNum.append(lin)
            lin += t.count("\n")

        # Extract per-token content positions (0-based line/col where actual
        # content starts after leading whitespace). Most token bodies contain
        # no internal '\n' — only their leading whitespace does (separators
        # were prepended during tokenize). BUT `_mergeQuoted` can produce a
        # body that spans lines for a multi-line quoted string like
        # `"line1\nline2"`, which TRUST C++ accepts in `System "..."` and
        # similar interpreters (audit 1.4 bug 2). So walk the body too,
        # handling embedded `\n` rather than assuming it's absent.
        cur_line, cur_col = 0, 0
        for t in self.tabToken:
            i = 0
            while i < len(t) and t[i] in (" ", "\t", "\n"):
                if t[i] == "\n":
                    cur_line += 1
                    cur_col = 0
                else:
                    cur_col += 1
                i += 1
            self.contentLine.append(cur_line)
            self.contentCol.append(cur_col)
            for ch in t[i:]:
                if ch == "\n":
                    cur_line += 1
                    cur_col = 0
                else:
                    cur_col += 1

        # Now validate the tokens, and empty slots of tabTokenLow where we have comments
        # Rules: - inside a pair of '#' everything is ignored (even '/*' or '*/')
        #        - if comment starts with '/*' then '#' is ignored and opening and closing '*/' must match
        inHash, starSlashCnt = False, 0
        for i, t in enumerate(self.tabTokenLow):
            emptyOnce = False
            if t == "#" and starSlashCnt == 0:
                inHash, emptyOnce = not inHash, True
            if t == "/*" and not inHash:
                starSlashCnt += 1
            if t == "*/" and not inHash:
                if starSlashCnt == 0:
                    raise ValueError(f"Invalid closing '*/' at line #{self.lineNum[i]}")
                starSlashCnt -= 1
                emptyOnce = True
            if starSlashCnt > 0 or inHash or emptyOnce:
                self.tabTokenLow[i] = ""
