#!/usr/bin/env python

"""Take a look at docs/README_for_trustify_devs.md

Module to handle parsing, loading and writing of
    - the XD tags found in the C++ code
    - the TRAD2.org file
    - the loading/writing of the TRAD2 file

Classes TRAD2Attr and TRAD2Block represent pieces of the TRAD2 file.
Class TRAD2Content is the central piece.

Authors: A Bruneton, C Van Wambeke
"""

import contextlib
import functools
import os
import re

from trustify.core.misc_utilities import logger, pretty_error
from trustify.modernize import is_unsplittable

# Single regex covering every XD tag variant — used as a cheap "is there any
# XD tag at all on this line?" probe. Matches:
#   //  XD                          (level 1)
#   //  2XD / 3XD                   (nested levels)
#   //  XD_CONT / XD_ADD_P / XD_ADD_DICO   (with optional 2/3 prefix)
# The tag must follow '//' with at least one space, and must end on whitespace
# or end-of-line. Compiled once at import — the previous per-call re.search
# pattern caused a measurable schema-extraction slowdown.
_ANY_XD_TAG_RE = re.compile(r"//\s+[23]?XD(?:_(?:CONT|ADD_P|ADD_DICO))?(?:\s|$)")

# Dispatch regex used by scanOneCppFile — same shape as _ANY_XD_TAG_RE but
# excludes XD_CONT (already merged into the parent line by
# extractAndGroupXDLines) and captures the level prefix in a named group so
# we can dispatch to scanOneCppLine without iterating all three levels.
_DISPATCH_TAG_RE = re.compile(r"//\s+(?P<lvl>[23]?)XD(?:_(?:ADD_P|ADD_DICO))?(?=\s|$)")
_LVL_INDEX = {"": 0, "2": 1, "3": 2}


@functools.cache
def _compiled_tag_re(tag):
    """Cached per-tag compiled regex. At most ~12 unique tags are ever
    requested (3 levels x {XD, XD_CONT, XD_ADD_P, XD_ADD_DICO}), so the
    cache size is bounded by the call sites.
    """
    return re.compile(r"//\s+" + re.escape(tag) + r"(?:\s|$)")


def _line_carries_tag(line, tag):
    """Return True when ``line`` opens a C++ comment whose first token is
    exactly ``tag``.

    Matches ``// <tag>`` (with one or more spaces between ``//`` and the
    tag) and rejects substring occurrences inside ordinary comments — so a
    sentence like ``// the XD parent is field_base`` no longer trips the
    scanner. The ``//`` may sit anywhere in the line so that inline
    XD_ADD_P / XD_ADD_DICO annotations on a ``param.ajouter*`` call still
    work.

    ``tag`` is the full tag string (e.g. ``"XD"``, ``"2XD"``,
    ``"XD_ADD_P"``, ``"XD_CONT"``). The tag must be followed by
    whitespace or end-of-line — never by another identifier character.

    Fast path: the vast majority of source lines don't contain ``XD`` as
    a substring, so a cheap ``"XD" in line`` test short-circuits before
    the regex search.
    """
    if "XD" not in line:
        return False
    return _compiled_tag_re(tag).search(line) is not None


def convertTyp(typ):
    """Types declared in XD tags of cpp files might differ from what is officially supported in TRAD_2, so convert:"""
    mp = {"int": "entier", "double": "floattant", "flag": "rien", "dico": "chaine(into=[])"}
    return mp.get(typ, typ)


# String / paren-aware C++ tokenization helpers (audit 6.1). The previous
# naive `line.split(",")` + `find("(")` / `find(")")` silently misparsed
# any `ajouter*` call with a nested paren arg like `(this)` or a quoted
# comma, masking Param::REQUIRED in lambda_ortho-style sources and
# truncating synonym strings on `,`.
def _split_cpp_args(s: str) -> list[str]:
    """Split `s` on TOP-LEVEL commas only — commas inside string literals
    or inside `(...)` / `[...]` / `{...}` groups are ignored. Each
    returned arg is stripped of outer whitespace.

    Angle brackets are intentionally NOT treated as nesting: C++ uses
    `<` and `>` for both templates and comparisons, and the contexts
    we tokenize never need template-aware splitting.
    """
    args: list[str] = []
    depth = 0
    in_string = False
    string_char: str | None = None
    escaped = False
    current: list[str] = []
    for ch in s:
        if in_string:
            current.append(ch)
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == string_char:
                in_string = False
            continue
        if ch in ('"', "'"):
            in_string = True
            string_char = ch
            current.append(ch)
            continue
        if ch in "([{":
            depth += 1
            current.append(ch)
            continue
        if ch in ")]}":
            depth -= 1
            current.append(ch)
            continue
        if ch == "," and depth == 0:
            args.append("".join(current).strip())
            current = []
            continue
        current.append(ch)
    if current or args:
        args.append("".join(current).strip())
    return args


def _extract_call_args(s: str, method: str) -> tuple[str, str, str]:
    """Locate `<method>[_suffix](...)` in `s`, returning `(prefix, args, suffix)`.

    `method` is treated as a prefix on the C++ method name (so
    `method="ajouter"` matches `ajouter`, `ajouter_non_std`,
    `ajouter_flag`, ...). The opening `(` is the first one after the
    matched method-name region (allowing whitespace between). The
    closing `)` is found by matched-paren walking that respects string
    literals and nested groups, so a preceding unrelated `(... )` —
    e.g. `if (cond) param.ajouter(...)` — does NOT confuse the
    extraction, unlike the previous `find("(")` / `find(")")` approach.

    Raises ValueError if no `<method>` call is found, or if the parens
    aren't balanced before end-of-string.
    """
    m = re.search(rf"\b{re.escape(method)}\w*", s)
    if m is None:
        raise ValueError(f"no {method!r} call found in {s!r}")
    # First `(` after the method name (allow whitespace).
    i = m.end()
    while i < len(s) and s[i].isspace():
        i += 1
    if i >= len(s) or s[i] != "(":
        raise ValueError(f"no '(' after {method!r} in {s!r}")
    open_idx = i
    depth = 0
    in_string = False
    string_char: str | None = None
    escaped = False
    for j in range(open_idx, len(s)):
        ch = s[j]
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == string_char:
                in_string = False
            continue
        if ch in ('"', "'"):
            in_string = True
            string_char = ch
            continue
        if ch == "(":
            depth += 1
            continue
        if ch == ")":
            depth -= 1
            if depth == 0:
                return s[: m.start()], s[open_idx + 1 : j], s[j + 1 :]
    raise ValueError(f"unbalanced parens in {method!r} call: {s!r}")


def _call_parens_balanced(s: str) -> bool:
    """Return True once `s` holds a complete parenthesised call — at least
    one '(' seen and every '(' matched by a ')'. String literals are
    respected so a paren inside a quoted argument is not counted.

    Used by `extractAndGroupXDLines` to tell whether a macro invocation
    (`Implemente_instanciable(...)` / `Add_synonym(...)`) closes on the
    current line or spills onto continuation lines (clang-format splits
    long class names across several lines).
    """
    depth = 0
    seen = False
    in_string = False
    string_char: str | None = None
    escaped = False
    for ch in s:
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == string_char:
                in_string = False
            continue
        if ch in ('"', "'"):
            in_string = True
            string_char = ch
        elif ch == "(":
            depth += 1
            seen = True
        elif ch == ")":
            depth -= 1
    return seen and depth <= 0


def _net_paren_depth(s: str) -> int:
    """Net '(' minus ')' count of `s`, ignoring parens inside string
    literals. Zero for a balanced or paren-free line; positive when an
    opening paren is left unclosed (a call opens here and spills below);
    negative when a closing paren has no opener on this line (the line
    continues a call opened on an earlier line). Used to recognise a
    `Param::ajouter*()` call split across several lines.
    """
    depth = 0
    in_string = False
    string_char: str | None = None
    escaped = False
    for ch in s:
        if in_string:
            if escaped:
                escaped = False
            elif ch == "\\":
                escaped = True
            elif ch == string_char:
                in_string = False
            continue
        if ch in ('"', "'"):
            in_string = True
            string_char = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
    return depth


# Strict syntax for the `type` field of an XD attr declaration. Anything
# matching here (and nothing else) is a legal type string — the pydantic
# generator dispatches on this shape via .startswith / [13:].split("]")[0],
# so any garbage trailing the closing `)` of a parameterised form would
# be silently dropped (e.g. `chaine(into=["a","b"])toto` produces the
# same schema as `chaine(into=["a","b"])`, losing `toto`). Audit 4.9.
#
# Forms accepted:
#   * plain identifier  — covers typ_map members (entier, floattant,
#     chaine, rien, list*, listchaine*, listentier*), block names from
#     all_blocks, the suppress_param sentinel, and ref_<class> refs.
#   * chaine(into=[ "a", "b", ... ])  — optional trailing comma allowed
#     because the XD_ADD_DICO post-processing in scanOneCppLine appends
#     items as `,...,]` (see trad2_utilities.py:645). Empty list also
#     allowed because convertTyp maps the `dico` shorthand to
#     `chaine(into=[])`. An optional `,default="X"` clause is accepted
#     after the closing `]`; `X` must be one of the listed values. It is
#     recorded by the pydantic generator as display-only Field metadata
#     (`json_schema_extra["trust_default"]`, see
#     `trad2_pydantic._chaine_declared_default`) and surfaced in the
#     generated manual — it does NOT change the runtime default (first
#     choice / None), so parsing and round-trip are unaffected. Membership
#     of `X` is validated at code-generation time, not here, because
#     `dico` enums are still empty when this validator runs. The
#     `entier(into=...)` and `entier(max=N)` forms below take the same
#     `,default=N` clause (int instead of quoted string).
#   * entier(into=[ 1, -2, ... ])     — bare ints OR quoted ints; the
#     pydantic generator strips outer quotes before int(), so both are
#     legal. Negative allowed. Trailing comma allowed for symmetry.
#   * entier(max=N)                   — single int upper bound (negative
#     allowed). `entier(min=M)` is a lower bound, and `entier(min=M,max=N)`
#     a closed range (min always precedes max).
#
# `entier(into=...)`, `entier(max=N)` and `entier(min=M[,max=N])` accept
# the same optional `,default=N` documentation clause as `chaine(into=...)`
# (bare or quoted int). It is recorded as display-only metadata, never the
# runtime default; membership / range is validated at code-generation time.
_VALID_TYPE_RE = re.compile(
    r"\A(?:"
    r"[a-zA-Z_][a-zA-Z0-9_]*"
    r'|chaine\(into=\[(?:"[^"\]]*"(?:,"[^"\]]*")*,?)?\](?:,default="[^"\]]*")?\)'
    r'|entier\(into=\[(?:(?:-?\d+|"-?\d+")(?:,(?:-?\d+|"-?\d+"))*,?)?\](?:,default=(?:-?\d+|"-?\d+"))?\)'
    r'|entier\(max=-?\d+(?:,default=(?:-?\d+|"-?\d+"))?\)'
    r'|entier\(min=-?\d+(?:,max=-?\d+)?(?:,default=(?:-?\d+|"-?\d+"))?\)'
    r")\Z"
)


def _validate_type_string(typ: str, fname: str, lineno: int) -> None:
    """Raise with a clear message if `typ` is not one of the documented
    type forms. See `_VALID_TYPE_RE` above for the accepted alphabet.
    """
    if _VALID_TYPE_RE.match(typ) is None:
        raise Exception(
            pretty_error(
                fname,
                lineno,
                f"illegal type string {typ!r} in 'XD attr' instruction. "
                "Expected a plain identifier (a TRUST class name, or one of "
                "entier|floattant|chaine|rien|list*|listchaine*|listentier*|"
                "suppress_param|ref_<class>), or one of the parameterised "
                'forms chaine(into=["a","b",...]) / entier(into=[1,2,...]) / '
                "entier(max=N) with NOTHING trailing the closing ')'.",
            )
        )


def _validate_xd_name(name: str, fname: str, lineno: int, kind: str) -> None:
    """Reject a keyword/attribute `name` carrying a character outside
    [A-Za-z0-9_]. TRUST names are documented as that class only.

    Validated here at extraction time so the diagnostic carries the real
    source location (the C++/.xd file and line). Otherwise the illegal
    name surfaces only much later — as a bare ValueError in the pydantic
    generator (`valid_variable_name`, for attribute names) or an invalid
    Python identifier in the generated module (for keyword names) —
    neither of which points back at the XD tag that needs fixing.
    """
    bad = sorted(set(re.findall(r"[^0-9A-Za-z_]", name)))
    if bad:
        raise Exception(
            pretty_error(
                fname,
                lineno,
                f"illegal {kind} name {name!r}: contains character(s) {''.join(bad)!r} "
                "outside [A-Za-z0-9_]. Fix the XD tag in the source.",
            )
        )


class TRAD2Attr:
    """An attribute of a block in the TRAD2 logic"""

    def __init__(self):
        self.name = ""  #: Name of the attribute
        self.type = ""  #: Type of the attribute
        self.synos = []  #: List of synonyms for the attribute
        self.is_opt = True  #: Is the attribute optional
        self.desc = ""  #: Description
        self.info = ["", -1]  #: Filename / Lineno where the attribute was declared

    @classmethod
    def BuildFromTab(cls, tab, fname, lineno, convert=False, trust_root=None, synthetic_opt=False):  # noqa: ARG003
        """Build a TRAD2Attr from a single array representing a (splitted) line of the TRAD2, or a C++ XD comment
        @param convert: whether to convert types because we come from C++
        @param synthetic_opt: when True, the opt token in `tab` was not
            written by the user — the caller derived it programmatically
            (typically from a `Param::REQUIRED` C++ marker on an
            ``XD_ADD_P`` / ``XD_ADD_DICO`` line). Suppresses the
            legacy-opt-flag warning; the user has nothing to modernize.

        ``trust_root`` is part of a back-compat keyword chain
        (``BuildContentFromTRAD2`` -> ``_ParseXD`` -> here) — accepted for
        signature consistency, ignored.
        """
        if len(tab) < 4:
            raise Exception(
                pretty_error(fname, lineno, "incomplete 'XD attr' (attribute line) instruction!!")
            ) from None
        a = TRAD2Attr()
        names, _, syno, opt = [t.lower() for t in tab[:4]]

        # For attributes, we accept '|' in both the first slot and the third one:
        a.name = names.split("|")[0]
        _validate_xd_name(a.name, fname, lineno, "attribute")
        syno2 = set((names + "|" + syno).split("|"))
        with contextlib.suppress(KeyError):
            syno2.remove(a.name)
        a.synos = list(syno2)
        a.synos.sort()

        # Type must sometimes preserve upper case ... (chaine=into["X","Y"])
        if tab[1].startswith("chaine="):
            a.type = tab[1]
        else:
            a.type = tab[1].lower()

        # Types are different if coming from TRAD2 of from C++ ...
        if convert:
            a.type = convertTyp(a.type)
        _validate_type_string(a.type, fname, lineno)

        a.desc = " ".join(tab[4:])
        # Flag for the optionality of an attribute : legacy codes (0 for 'REQ' and 1 for 'OPT') are still supported:
        opt = opt.upper()
        if opt not in ["0", "1", "REQ", "OPT"]:
            raise Exception(
                pretty_error(
                    fname,
                    lineno,
                    f"invalid optional flag in 'XD attr' (attribute line) instruction ('{opt}'). It should be 'OPT' or 'REQ'!!",
                )
            ) from None
        if opt in ("0", "1") and not synthetic_opt:
            named = "REQ" if opt == "0" else "OPT"
            logger.warning(
                "%s:%d: legacy numeric opt flag '%s' — run `trustify modernize` to upgrade to '%s'.",
                fname,
                lineno + 1,
                opt,
                named,
            )
        a.is_opt = opt == "OPT" or opt == "1"
        # `info` is left as the absolute path here; the post-scan pass
        # `TRAD2Content._convert_infos_to_source_locations` (run from
        # `BuildFromOrgAndSources`) rewrites it into a SourceLocation
        # using the projects dict.
        a.info = [fname, lineno + 1]
        return a

    def toTRAD2Text(self) -> str:
        """Return the canonical TRAD2 single-line text for this attribute."""
        opt = "OPT" if self.is_opt else "REQ"
        synos = self.synos if self.synos else [self.name]
        synos_s = "|".join(synos)
        return f"  attr {self.name} {self.type} {synos_s} {opt} {self.desc}\n"


class TRAD2Block:
    """A block describing a keyword in the TRAD2 file."""

    def __init__(self):
        self.name = ""  #: Keyword main name
        self.name_base = ""  #: Parent class for the keyword
        self.synos = []  #: List of synonyms
        self.mode = -123  #: Mode of the keyword (with braces, etc... see docs/README_for_trust_devs.md)
        self.desc = ""  #: Description
        self.info = ["", -1]  #: Filename / Lineno where the keyword was defined
        self.attrs = []  #: A list of TRAD2Attr = the attributes expected for the keyword
        self.attr_synos = None  #: Dictionary of (all) attribute synonyms - filled in by the trad2_pydantic module

    @classmethod
    def BuildFromTab(cls, tab, fname, lineno, trust_root=None):  # noqa: ARG003
        """Build from a single array representing a (splitted) line of the TRAD2 or a C++ XD comment.

        ``trust_root`` is part of a back-compat keyword chain
        (``BuildContentFromTRAD2`` -> ``_ParseXD`` -> here) — accepted for
        signature consistency, ignored.
        """
        base_nam = tab[1].lower()
        b = TRAD2BlockList() if base_nam == "listobj" else TRAD2Block()
        b.name, b.name_base, nam2, acco_s = [t.lower() for t in tab[:4]]
        _validate_xd_name(b.name, fname, lineno, "keyword")
        if b.name == b.name_base:
            raise Exception(
                pretty_error(
                    fname,
                    lineno,
                    f"Keyword/class '{b.name}' inherits from itself!! You should put a parent class as second parameter in the XD line.",
                )
            ) from None
        a = None
        # Tag for usage of braces : legacy mode is supported: "INHERITS_BRACE" (-1), "NO_BRACE" (0), "BRACE" (1)
        acco_s = acco_s.upper()
        if acco_s in ["INHERITS_BRACE", "NO_BRACE", "BRACE"]:
            a = acco_s
        else:
            try:
                a_int = int(acco_s)
                a = {
                    -3: "BRACE",  # -3 (resp) -2:  Historically like 1 (resp 0) but also needs to be put after Discretisation (e.g. for all Problems, discretisation must be done before reading them)
                    -2: "NO_BRACE",
                    -1: "INHERITS_BRACE",
                    0: "NO_BRACE",
                    1: "BRACE",
                }.get(a_int)
                if a is not None:
                    logger.warning(
                        "%s:%d: legacy numeric brace flag '%s' — run `trustify modernize` to upgrade to '%s'.",
                        fname,
                        lineno + 1,
                        acco_s,
                        a,
                    )
            except ValueError:
                pass
        if a is None:
            raise Exception(
                pretty_error(
                    fname,
                    lineno,
                    f"option for curly braces should either 'INHERITS_BRACE', or 'NO_BRACE', or 'BRACE', not '{acco_s}'!!",
                )
            ) from None
        b.mode = a
        b.synos = nam2.split("|")
        # See note in TRAD2Attr.BuildFromTab — the post-scan pass
        # rewrites this into a SourceLocation using the projects dict.
        b.info = [fname, lineno + 1]
        b._finishBuild(tab[4:])
        return b

    def _finishBuild(self, tab):
        self.desc = " ".join(tab)

    def toTRAD2Text(self) -> str:
        """Return the canonical TRAD2 text for this block + its attrs."""
        synos = "|".join(self.synos)
        s = f"{self.name} {self.name_base} {synos} {self.mode} {self.desc}\n"
        for a in self.attrs:
            s += a.toTRAD2Text()
        return s


class TRAD2BlockList(TRAD2Block):
    """Specific class for 'listobj' entries in the TRAD2 (describing a list of objects)"""

    def __init__(self):
        TRAD2Block.__init__(self)
        self.itemtype = None  #: the list item type
        self.comma = -123  #: Whether the list takes commas or not

    def _finishBuild(self, tab):
        """Override to extract list-relevant data"""
        self.itemtype = tab[0].lower()
        comma_dict = {"-1": "INHERITS_COMMA", "0": "NO_COMMA", "1": "COMMA"}
        v = tab[1].upper()
        if v not in (list(comma_dict.keys()) + list(comma_dict.values())):
            raise Exception(
                pretty_error(
                    self.info[0],
                    self.info[1],
                    f"option for comma (for a list) should be in {comma_dict.values()}, not '{tab[1]}'!!",
                )
            ) from None
        # Flag for comma in lists: either INHERITS_COMMA (formerly -1), 'NO_COMMA' (formerly 0), or 'COMMA' (formerly 1)
        self.comma = comma_dict.get(v, v)
        self.desc = " ".join(tab[2:])

    def toTRAD2Text(self) -> str:
        """Return the canonical TRAD2 text for this listobj entry.
        ``listobj`` blocks must not carry attributes.
        """
        from trustify.core.source_location import SourceLocation

        if len(self.attrs):
            if isinstance(self.info, SourceLocation):
                _path, _line = self.info.path, self.info.line
            else:
                _path, _line = self.info[0], self.info[1]
            raise Exception(
                pretty_error(_path, _line, "list object description should not have any attribute!!")
            ) from None
        synos = "|".join(self.synos)
        return f"{self.name} {self.name_base} {synos} {self.mode} {self.itemtype} {self.comma} {self.desc}\n"


class TRAD2Content:
    """Main class - this is a memory version of a TRAD2 file.
    A TRAD2Content is made of a list of TRAD2Block, each containing a list of TRAD2Attr.
    """

    def __init__(self, projects=None):
        self.data = []  #: A list of TRAD2Block
        self.projects = dict(projects) if projects else {}
        # Errors found while scanning C++/.xd sources for XD tags. Accumulated
        # rather than raised on first hit (see scanOneCppFile / scanSourceFiles)
        # so one pass reports every malformed tag, not just the first.
        self._scan_errors = []
        self.re_comment = re.compile(r"^\s*//")  #: RegExp to match simple one line comment
        # Tag-anchored: the XD_CONT keyword must directly follow '//', with
        # at least one space. '//XD_CONT' (no space) and 'foo // XD_CONT'
        # (XD_CONT not the first token) are both rejected. Levels 2/3
        # match too so 2XD_CONT / 3XD_CONT fold into their respective
        # 2XD_ADD_P / 3XD_ADD_P openers.
        self.re_cont = re.compile(r"^\s*//\s+[23]?XD_CONT(?:\s|$)")

    @classmethod
    def BuildContentFromTRAD2(cls, trad2, source_locations=None, trust_root=None):
        """Build content from a TRAD2 file, optionally augmenting per-block
        / per-attr `info` with `SourceLocation` triples loaded from a
        sibling `source_locations.json`.

        - `source_locations=None` (default): look for
          `<dir>/source_locations.json` next to `trad2`. Load if
          present; silently fall through to no source locations if
          absent. Hand-written TRAD2 files without a sibling JSON are
          officially supported — `info` stays empty, and consumers
          that display C++-source backrefs simply omit them.
        - `source_locations="<path>"`: load the explicit path. Raises
          if the path does not exist (callers passing an explicit
          path expect the file to be there).
        - `trust_root` is accepted for back-compat but ignored — the
          JSON carries project names directly.
        """
        import json

        from trustify.core.source_location import SourceLocation

        trad2 = str(trad2)
        res = TRAD2Content()
        with open(trad2, encoding="utf-8") as f:
            lines = f.readlines()

        locs = None  # name -> SourceLocation; attrs nested under blk.attrs
        explicit_path = source_locations is not None
        json_path = source_locations if explicit_path else os.path.join(os.path.dirname(trad2), "source_locations.json")
        if os.path.exists(json_path):
            with open(json_path, encoding="utf-8") as f:
                payload = json.load(f)
            if payload.get("version", 1) != 1:
                raise Exception(f"unsupported source_locations.json version: {payload.get('version')!r}")
            locs = payload.get("entries", {})
        elif explicit_path:
            raise FileNotFoundError(f"source_locations file does not exist: {json_path}")

        curr_obj = None
        for lin_n, line in enumerate(lines):
            if line.strip() == "":
                continue
            line = line.strip()
            tab = [t for t in line.split(" ") if t.strip() != ""]
            if len(tab) >= 2 and tab[0] == "ref":
                continue
            obj = cls._ParseXD(trad2, lin_n, tab, convert=False, trust_root=trust_root)
            if isinstance(obj, TRAD2Block):
                res.data.append(obj)
                curr_obj = obj
                if locs is not None and obj.name in locs:
                    e = locs[obj.name]
                    obj.info = SourceLocation(project=e["project"], path=e["path"], line=int(e["line"]))
            else:
                if curr_obj is None:
                    raise Exception(
                        pretty_error(trad2, lin_n, "'XD attr' (attribute line) read before main 'XD' line!")
                    )
                curr_obj.attrs.append(obj)
                if locs is not None and curr_obj.name in locs:
                    a = locs[curr_obj.name].get("attrs", {}).get(obj.name)
                    if a is not None:
                        obj.info = SourceLocation(project=a["project"], path=a["path"], line=int(a["line"]))
        return res

    @classmethod
    def BuildFromOrgAndSources(cls, trad2_org, src_dirs, projects=None, trust_root=None):
        """Build content from TRAD2.org file (optional) and from a set of
        C++/xd source directories.

        `projects` is a `{name: absolute_path}` dict (see
        `trustify.core.source_location.build_projects_dict`). When
        provided, it populates `self.projects` and is threaded into the
        per-declaration `SourceLocation` resolution. When omitted but
        `trust_root` is given, we synthesize `{"trust": trust_root}` so
        legacy callers stop the conversion pass from being skipped
        (which would leave `obj.info` as `[path, line]` lists and break
        `serialize_source_locations`). The legacy `trust_root` keyword
        is kept temporarily; Task 12 will remove it once all callers
        pass `projects=` directly.
        """
        if trad2_org is not None and not os.path.exists(trad2_org):
            raise FileNotFoundError(f"TRAD2.org file not found: {trad2_org!r}")
        # Defensive: `trustify.api.generate_schema` runs `validate_projects`
        # upstream which already rejects baltiks missing `src/` with a
        # clearer `ConfigError`. This check catches callers that build a
        # `TRAD2Content` directly (tests, ad-hoc scripts) and surfaces
        # the bad path explicitly.
        for d in src_dirs:
            if not os.path.exists(d):
                raise FileNotFoundError(
                    f"Source directory not found: {d!r}. "
                    "Every project passed to `BuildFromOrgAndSources` must "
                    "have an existing `src/` subdirectory. If you're calling "
                    "trustify's public API, this should have been caught earlier "
                    "by `validate_projects`."
                )
        if projects is None and trust_root is not None:
            from trustify.core.source_location import TRUST_PROJECT_NAME

            projects = {TRUST_PROJECT_NAME: str(trust_root)}
        ret = TRAD2Content(projects=projects)
        ret.synos = {}
        ret.scanSourceFiles(src_dirs, trust_root=trust_root)
        ret.assemble(trad2_org, trust_root=trust_root)
        if ret.projects:
            ret._convert_infos_to_source_locations()
        return ret

    def serialize_source_locations(self) -> str:
        """Return the canonical JSON text for `source_locations.json`.

        This is the byte-stable form fed into the cache hash, so:
        - `sort_keys=True` for deterministic key order;
        - compact separators (`","`, `":"`) for size + reproducibility;
        - integer `line`, never float.
        """
        import json

        from trustify.core.source_location import SourceLocation

        entries = {}
        for blk in self.data:
            if not isinstance(blk.info, SourceLocation):
                # Should never happen post-conversion; surface loudly.
                raise TypeError(
                    f"serialize_source_locations: block {blk.name!r} has info={blk.info!r} (expected SourceLocation)"
                )
            attrs = {}
            for a in blk.attrs:
                if not isinstance(a.info, SourceLocation):
                    raise TypeError(
                        f"serialize_source_locations: attr {a.name!r} of block {blk.name!r} has info={a.info!r}"
                    )
                attrs[a.name] = a.info.to_dict()
            d = blk.info.to_dict()
            d["attrs"] = attrs
            entries[blk.name] = d
        payload = {"version": 1, "entries": entries}
        return json.dumps(payload, sort_keys=True, separators=(",", ":"))

    def toTRAD2(self, f_nam_out):
        """Write the TRAD2 text and the sibling `source_locations.json`
        into `<dir>/source_locations.json` (a single canonical name per
        cache directory).

        The two files are written sequentially (not atomically) — both
        feed into the cache hash, so a crash between them would leave
        the directory inconsistent. Callers that need atomicity (the
        only production caller is `api.generate_schema`) must wrap the
        invocation in `cache.atomic_write_dir`; the temp-dir-then-rename
        there makes the WHOLE cache entry (both files + the generated
        `.py` modules + `provenance.json`) appear or disappear as one.
        See audit findings 2.1 and 2.5.
        """
        tot_s = [blk.toTRAD2Text() for blk in self.data]
        with open(f_nam_out, "w", encoding="utf-8") as f:
            f.write("".join(tot_s))
        logger.info(f"==> Written file '{f_nam_out}'")
        out_dir = os.path.dirname(os.path.abspath(f_nam_out)) or "."
        json_path = os.path.join(out_dir, "source_locations.json")
        with open(json_path, "w", encoding="utf-8") as f:
            f.write(self.serialize_source_locations())
        logger.info(f"==> Written file '{json_path}'")

    def toTRAD2Text(self) -> str:
        """Return the canonical TRAD2 text serialization as a string.
        Equivalent to what `toTRAD2(fname)` writes to `fname` (the
        sibling `source_locations.json` is NOT included). Used by
        `trustify.cache.compute_content_hash` to derive a stable cache
        key without disk I/O.
        """
        return "".join(blk.toTRAD2Text() for blk in self.data)

    def toJSon(self):
        """This is the future ... :-)"""

    def _parseMacro(self, tag, line):
        """Generic method for parsing a 'implemente_instanciable' or 'add_synonym()'
        All trailing '_64' or '_32_64' are removed. 64b is handled when building the synonyms.

        Uses `_extract_call_args` / `_split_cpp_args` so a synonym
        string containing a comma (`Add_synonym(Foo, "a,b")`) is no
        longer truncated at the embedded `,` (audit 6.1).
        """
        if tag not in line:
            return None
        _, args_text, _ = _extract_call_args(line, method=tag)
        args = _split_cpp_args(args_text)
        cls_nam = args[0].strip()
        if cls_nam.endswith("_32_64"):
            cls_nam = cls_nam[: -len("_32_64")]
        s = args[1].strip().replace('"', "")
        if s.endswith("_64"):
            # Slice off ONLY the trailing `_64` — str.replace would
            # also drop any mid-string `_64` (e.g. `foo_64bit_64`
            # would collapse to `foobit`). Mirrors the `_32_64`
            # stripping above. Audit 1.3.
            s = s[: -len("_64")]
        return cls_nam, s

    @classmethod
    def _ParseXD(cls, f_nam, lin_n, tab, convert=False, trust_root=None):
        """Parse a TRAD2 line, or a '// XD' line in the C++ file
        @param tab The line is already split on whitespace and stored in 'tab'
        @param f_nam, lin_n: the origin of the declaration
        """
        if len(tab) < 5:
            raise Exception(pretty_error(f_nam, lin_n, "incomplete 'XD' instruction (header line)!!"))

        if tab[0] == "attr":
            # XD param tag: "// XD attr sutherland bloc_sutherland sutherland 1 Sutherland law for viscosity ..."
            attr = TRAD2Attr.BuildFromTab(tab[1:], f_nam, lin_n, convert=convert, trust_root=trust_root)
            return attr
        else:
            # XD main tag:  "// XD fluide_quasi_compressible fluide_dilatable_base fluide_quasi_compressible -1 Quasi-compressible flow with a low ...."
            blk = TRAD2Block.BuildFromTab(tab, f_nam, lin_n, trust_root=trust_root)
            return blk

    def _parseXD_ADD_something(
        self, tag="XD_ADD_P", expected_num_args=2, f_nam="", cpp_meth="ajouter", lin_n=-1, tab=None
    ):
        """Parse a line with XD_ADD_P or XD_ADD_DICO -> this implies extracting some
        information from the C++ call: param.ajouter(...).
        Only used when parsing C++ stuff.
        @param tag: XD tag that is expected on this line
        @param expected_num_args: number of args expected after the XD tag (not by tha 'ajouter..' method!!)
        """
        #  Example
        #   "param.ajouter_flag("P0", &alphaE_); // XD_ADD_P rien Pressure nodes"
        if tab is None:
            tab = []
        pos = tab.index(tag)
        assert pos >= 0
        if len(tab[pos + 1 :]) < expected_num_args:
            raise Exception(pretty_error(f_nam, lin_n, f"incomplete {tag} instruction!!"))
        typ = tab[pos + 1].lower()
        desc = " ".join(tab[pos + 2 :]) if expected_num_args > 1 else typ
        # Parse C++ parameters passed to the ajouter_ method
        lin_start = " ".join(tab[:pos])
        # A multi-line C++ call is not supported: the tag is a `//` comment
        # embedded mid-statement, and the call may even open on a prior
        # line that the line-grouping already dropped. An unbalanced paren
        # count on the tag line is the tell-tale (depth > 0: call opens
        # here and spills below; depth < 0: tag sits on the closing line
        # of a call opened earlier).
        depth = _net_paren_depth(lin_start)
        multiline_msg = (
            f"multi-line 'Param::{cpp_meth}*()' call carrying a '{tag}' tag is not supported "
            "- put the call and its tag on a single line"
        )
        if f".{cpp_meth}" not in lin_start:
            if depth != 0:
                raise Exception(pretty_error(f_nam, lin_n, multiline_msg))
            raise Exception(
                pretty_error(
                    f_nam, lin_n, f"{tag} instruction put on a line where 'Param::{cpp_meth}*' method is not called!!"
                )
            )
        # Use balanced-paren + string-aware splitting (audit 6.1): the
        # previous naive `find("(")` / `find(")")` + `split(",")`
        # mis-located the call args whenever the C++ used a nested arg
        # such as `(this)` — e.g. `param.ajouter_non_std("lambda_ortho",
        # (this), Param::REQUIRED)` had `find(")")` stop at the (this)
        # paren, dropping Param::REQUIRED silently and marking the attr
        # OPT when it should be REQ.
        try:
            _, args_text, _ = _extract_call_args(lin_start, method=cpp_meth)
        except ValueError as e:
            # Call opens on this line but never closes -> spilled below.
            if depth > 0:
                raise Exception(pretty_error(f_nam, lin_n, multiline_msg)) from None
            raise Exception(pretty_error(f_nam, lin_n, f"could not parse Param::{cpp_meth}*() call: {e}")) from None
        tb = _split_cpp_args(args_text)
        if len(tb) not in [1, 2, 3]:
            raise Exception(pretty_error(f_nam, lin_n, "wrong number of arguments to Param::ajouter*() method??"))
        # Param name
        param_nam = tb[0].lower().strip().replace('"', "")
        # Required param?
        opt = "1"
        if len(tb) == 3 and "required" in tb[2].lower():
            opt = "0"
        return param_nam, typ, opt, desc

    def _parseXD_ADD_P(self, lvl_s, f_nam, lin_n, tab):
        """Parse a line containing a "XD_ADD_P" tag"""
        tag = f"{lvl_s}XD_ADD_P"
        return self._parseXD_ADD_something(
            tag=tag, expected_num_args=2, cpp_meth="ajouter", f_nam=f_nam, lin_n=lin_n, tab=tab
        )

    def _parseXD_ADD_DICO(self, lvl_s, f_nam, lin_n, tab):
        """Parse a line containing a "XD_ADD_DICO" tag"""
        tag = f"{lvl_s}XD_ADD_DICO"
        return self._parseXD_ADD_something(
            tag=tag, expected_num_args=1, cpp_meth="dictionnaire", f_nam=f_nam, lin_n=lin_n, tab=tab
        )

    def scanOneCppLine(self, lvl, lin, curr_obj, res, f_name, lin_n, trust_root=None):
        """Scan a single line of a C++ file"""
        lvl_s = ["", "2", "3"][lvl]
        XD = f"{lvl_s}XD"  # "XD", or "2XD", or "3XD"

        # Make sure this is a valid XD comment
        tab0 = lin.split(" ")
        tab = [t for t in tab0 if t.strip() != ""]
        if f"//{XD}" in tab or f"//{XD}_ADD_P" in tab:
            raise Exception(
                pretty_error(f_name, lin_n, f"misformatted {XD} comment (should have space before {XD} ...)")
            )

        #
        # This part of the treatment is exactly the same as what we do with raw TRAD_2 data:
        # The _line_carries_tag check requires the tag to follow '//' directly,
        # so a stray 'XD' token inside an English comment can no longer match.
        #
        if _line_carries_tag(lin, XD) and XD in tab:
            pos = tab.index(XD)
            # Skip "XD ref ..." for now
            if len(tab[pos:]) >= 3 and tab[pos + 1].lower() == "ref":
                return
            obj = self._ParseXD(f_name, lin_n, tab[pos + 1 :], convert=True, trust_root=trust_root)
            if isinstance(obj, TRAD2Block):
                res.append(obj)
                curr_obj[lvl] = obj
            else:  # attribute
                if curr_obj[lvl] is None:
                    raise Exception(pretty_error(f_name, lin_n, "'XD attr' read before main 'XD' line!"))
                curr_obj[lvl].attrs.append(obj)
        #
        # Those last two are only found in C++ files: XD_ADD_P and XD_ADD_DICO
        #
        if _line_carries_tag(lin, f"{XD}_ADD_P") and f"{XD}_ADD_P" in tab:
            if curr_obj[lvl] is None:
                raise Exception(pretty_error(f_name, lin_n, "'XD attr' read before main 'XD' line!"))
            nam1, typ, opt, desc = self._parseXD_ADD_P(lvl_s, f_name, lin_n, tab)
            # synthetic_opt=True: the `opt` value here was derived
            # from the C++ `Param::REQUIRED` marker, not written by
            # the user. Suppresses the legacy-opt nudge — there is
            # no `XD_ADD_P` syntax for the user to modernize toward.
            a = TRAD2Attr.BuildFromTab(
                [nam1, typ, nam1, opt, desc],
                f_name,
                lin_n,
                convert=True,
                trust_root=trust_root,
                synthetic_opt=True,
            )
            curr_obj[lvl].attrs.append(a)
        if _line_carries_tag(lin, f"{XD}_ADD_DICO") and f"{XD}_ADD_DICO" in tab:
            if curr_obj[lvl] is None:
                raise Exception(pretty_error(f_name, lin_n, "'XD attr' read before main 'XD' line!"))
            nam1, _, _, _ = self._parseXD_ADD_DICO(lvl_s, f_name, lin_n, tab)
            # Ensure the last attr added in the block is a 'dico':
            last_attr = curr_obj[lvl].attrs[-1]
            typ = last_attr.type
            if not typ.startswith("chaine(into=["):  # bof bof ...
                raise Exception(
                    pretty_error(
                        f_name, lin_n, "'XD_ADD_DICO' read, but no preceding 'XD_ADD_P dico ...' instruction found!!"
                    )
                )
            last_attr.type = typ.replace("]", f'"{nam1}",]')

    def extractAndGroupXDLines(self, f_name, lins):
        """Group multi-lines (XD comments expanding on several lines).
        @return a list of sub-lists, each sub-list having exactly 3 items: line number / (joined) original lines / (joined) stripped, lower-case lines

        Hot path on the schema-extraction critical path: every C++/.xd
        file under \\$TRUST_ROOT is fed through this. A pre-screen skips
        the ~99% of source lines that can't possibly affect the schema —
        no `XD`, no macro substring — before paying for `.lower()` or
        regex matches. Macro names are written in canonical case in
        TRUST (only `Implemente_instanciable` and `Add_synonym(` ever
        appear), so case-sensitive substring checks are sufficient for
        the screen.
        """
        grouped_lines = []  # triplet: line number / original line / stripped, lower-case line
        curr_block = [-1, "", ""]
        in_block = False
        in_macro = False  # mid multi-line Implemente_instanciable / Add_synonym call
        re_cont = self.re_cont  # local alias avoids attribute lookup per iteration
        for lin_n, lin in enumerate(lins):
            # Continuation of a multi-line macro call: these lines carry
            # no XD/macro substring, so the pre-screen below would drop
            # them. Accumulate into the open macro block until its parens
            # balance, then let the normal flow consume the block.
            if in_macro:
                s = lin.strip()
                curr_block[1] = " ".join([curr_block[1], s])
                curr_block[2] = " ".join([curr_block[2], s.lower()])
                if _call_parens_balanced(curr_block[2]):
                    in_macro = False
                continue
            # Cheap pre-screen on the raw (unstripped) line: lines that
            # carry no XD tag and no macro invocation cannot start,
            # continue, or close an XD block in a meaningful way — drop
            # them straight away, BEFORE paying for `.strip()`. The strip
            # call would dominate per-line cost otherwise (~500k calls).
            if "XD" not in lin and "Implemente_instanciable" not in lin and "Add_synonym(" not in lin:
                # A non-relevant line still marks the end of any open XD block;
                # the closed block is recorded at the next opener (or at EOF).
                in_block = False
                continue
            lin = lin.strip()

            is_cont = re_cont.match(lin) if "XD_CONT" in lin else None
            if is_cont:
                # Remove '// XD_CONT' span and re-strip; the line may then
                # become empty or be the actual continuation content.
                start, end = is_cont.span()
                lin = (lin[:start] + lin[end:]).strip()
            l_low = lin.lower()
            is_comment = lin.startswith("//")  # stripped above, so this is exact
            if is_cont and not in_block:
                raise Exception(
                    pretty_error(
                        f_name,
                        lin_n,
                        "'XD_CONT' improperly used! It should be placed just after a line with XD or XD_ADD_P.",
                    )
                )
            # These macros are usually one line, but a long class name may
            # push clang-format to split the call across several lines.
            cond1 = not is_comment and "implemente_instanciable" in l_low
            cond2 = not is_comment and "add_synonym(" in l_low
            if cond1 or cond2:
                in_block = False
                if curr_block[0] >= 0:
                    grouped_lines.append(curr_block)
                curr_block = [lin_n, lin, l_low]
                # Parens don't close on this line -> accumulate the
                # continuation lines (handled at the top of the loop).
                in_macro = not _call_parens_balanced(l_low)
                continue
            # Multi-line XD blocks: only treat the line as one if a real
            # '// <TAG>' opener is present — not just any occurrence of the
            # word XD inside an English comment. The pre-screen already
            # guarantees "XD" appears in `lin`, so the single combined
            # regex is invoked at most once per candidate line.
            cond3 = not is_cont and _ANY_XD_TAG_RE.search(lin) is not None
            if cond3 and is_comment and len(lin) > 120 and not is_unsplittable(lin):
                logger.warning(
                    "%s: long XD line — use tag 'XD_CONT' or run `trustify modernize` to fix.",
                    f_name,
                )
            if cond3:
                # Save previously accumulated block:
                if curr_block[0] >= 0:
                    grouped_lines.append(curr_block)
                curr_block = [lin_n, lin, l_low]  # Start new block
                in_block = True
            else:
                # Ignore the irrelevant lines:
                if not (in_block and is_cont):
                    in_block = (
                        False  # Mark the end of the block (it will be registered at the next new block or at the end)
                    )
                    continue
                # Otherwise accumulate in the current block (without the XD_CONT tag!)
                curr_block[1] = " ".join([curr_block[1], lin])
                curr_block[2] = " ".join([curr_block[2], l_low])
        # Save last block (if any)
        if curr_block[0] >= 0:
            grouped_lines.append(curr_block)
        return grouped_lines

    def scanOneCppFile(self, f_name, trust_root=None):
        """Scan one C++ file for the XD tags. Also extract the synonyms given by the
        'add_synonym' macro, and extend 'self.synos'
        """
        impl = {}
        # C++ file may contain XD, but also 2XD, 3XD tags. Syntax is exactly the same.
        # See docs/README_for_trust_devs.md for details. Max 3 levels.
        curr_obj_per_level = [None, None, None]
        res_per_level = [[], [], []]
        #
        # Frist part, scan the file to extract XD stuff
        #
        logger.debug(f"File: {f_name}")
        with open(f_name, encoding="utf-8") as f:
            content = f.read()

        # File-level pre-screen: a file with no XD tag and no synonym macro
        # cannot contribute anything to the schema. About a third of the
        # TRUST source tree falls in this bucket — skip them entirely instead
        # of paying for splitlines + the per-line pre-screen.
        if "XD" not in content and "Implemente_instanciable" not in content and "Add_synonym(" not in content:
            return []

        # Line grouping can itself raise (broken XD_CONT structure). Record
        # and skip the whole file rather than aborting the entire scan.
        try:
            grouped_lines = self.extractAndGroupXDLines(f_name, content.splitlines())
        except Exception as e:
            self._record_scan_error(e)
            return []

        for lin_n, lin, l_low in grouped_lines:
            # Each grouped line is independent: _scanOneGroupedLine records
            # any error and the scan moves on, so one malformed tag doesn't
            # hide every later one.
            self._scanOneGroupedLine(impl, lin_n, lin, l_low, f_name, curr_obj_per_level, res_per_level, trust_root)

        #
        # Aggregate final result, by puting levels in order:
        #
        res = []
        for r in res_per_level:
            res.extend(r)

        return res

    def _scanOneGroupedLine(self, impl, lin_n, lin, l_low, f_name, curr_obj_per_level, res_per_level, trust_root):
        """Process one grouped source line, recording any error instead of
        raising so the scan continues to the next line. All accumulated
        errors are surfaced together at the end of the scan."""
        try:
            self._scanGroupedLineImpl(impl, lin_n, lin, l_low, f_name, curr_obj_per_level, res_per_level, trust_root)
        except Exception as e:
            self._record_scan_error(e)

    def _scanGroupedLineImpl(self, impl, lin_n, lin, l_low, f_name, curr_obj_per_level, res_per_level, trust_root):
        """Synonym-macro + XD-tag dispatch for a single grouped line. An
        early `return` here is the equivalent of `continue` in the original
        per-line loop in scanOneCppFile."""
        # For each possible instruction, an example is provided:

        #
        # The first two cases handle synonyms
        #
        if "implemente_instanciable" in l_low:
            # "Implemente_instanciable(Terme_Boussinesq_VEF_Face,"Boussinesq_VEF_P1NC",Terme_Boussinesq_base);"
            # A multi-line call is folded into one logical line by
            # extractAndGroupXDLines; if its parens still don't balance
            # here the source is genuinely malformed (call never
            # closed) -> turn the bare ValueError into a located error.
            try:
                cls_nam, kw = self._parseMacro("implemente_instanciable", l_low)
            except (ValueError, IndexError) as e:
                raise Exception(
                    pretty_error(f_name, lin_n, f"could not parse 'Implemente_instanciable' macro: {e}")
                ) from None
            kwt = kw.split("|")  # we might have several synos there already...
            impl[cls_nam] = kwt[0]
            for k in kwt[1:]:
                self.synos.setdefault(kwt[0], []).append(k)
        if "add_synonym(" in l_low:
            #  "Add_synonym(Terme_Boussinesq_VEF_Face,"Boussinesq_temperature_VEF_Face");"
            try:
                cls_nam, s = self._parseMacro("add_synonym", l_low)
            except (ValueError, IndexError) as e:
                raise Exception(pretty_error(f_name, lin_n, f"could not parse 'Add_synonym' macro: {e}")) from None
            if cls_nam.endswith("_64"):
                cls_nam = cls_nam[: -len("_64")]
            if cls_nam not in impl:
                raise Exception(
                    pretty_error(f_name, lin_n, "'Add_synonym' macro used before 'Implemente_instanciable'")
                )
            kw = impl[cls_nam]
            self.synos.setdefault(kw, [])
            if s not in self.synos[kw]:
                self.synos[kw].append(s)

        #
        # Now handle XD tags. One regex search with a level-capture group
        # replaces the previous 9-call `_line_carries_tag` matrix
        # (3 levels x {XD, XD_ADD_P, XD_ADD_DICO}). XD_CONT can never
        # appear here — extractAndGroupXDLines already stripped it.
        #
        if "XD" not in lin:
            return
        m = _DISPATCH_TAG_RE.search(lin)
        if m is None:
            return
        lvl = _LVL_INDEX[m.group("lvl")]
        self.scanOneCppLine(lvl, lin, curr_obj_per_level, res_per_level[lvl], f_name, lin_n, trust_root=trust_root)

    def _record_scan_error(self, exc):
        """Append a scan-time error message to the accumulator."""
        self._scan_errors.append(str(exc))

    def _raise_accumulated_scan_errors(self):
        """If any XD-extraction errors were accumulated, raise a single
        exception listing every one of them. No-op when none were found."""
        if not self._scan_errors:
            return
        n = len(self._scan_errors)
        parts = [f"{n} error(s) while extracting XD tags:"]
        for i, msg in enumerate(self._scan_errors, 1):
            parts.append(f"\n=== error {i}/{n} ==={msg}")
        raise Exception("".join(parts))

    def scanSourceFiles(self, src_dirs, trust_root=None):
        """Scan all C++ source files (.cpp) and standalone XD declaration files (.xd) from a
        list of root directories (typically $TRUST_ROOT/src). If several directories are provided,
        the last ones have priority if the same file is found twice (for BALTIKs overrides, this
        makes sure the BALTIK source file has precedence).
        .xd files use the same `// XD ...` grammar as in-cpp comments; only the file extension differs.
        """
        import glob

        g_all = {}
        for d in src_dirs:
            for ext in ("cpp", "xd"):
                pattern = os.path.join(d, "**", f"*.{ext}")
                # glob.glob returns arbitrary order — sort so the resulting
                # TRAD2 content (and thus the cache hash) is deterministic
                # across runs.
                g = sorted(glob.glob(pattern, recursive=True))
                fnames = [os.path.split(nam)[-1] for nam in g]
                gd = dict(zip(fnames, g, strict=False))
                g_all.update(gd)  # last project wins on basename collision
        self.synos = {}
        self._scan_errors = []
        for f_nam in g_all.values():
            res = self.scanOneCppFile(f_nam, trust_root=trust_root)
            for d in res:
                logger.debug(d.toTRAD2Text())
            self.data.extend(res)
        # Every malformed tag across every file is reported at once here,
        # rather than aborting the scan on the first one.
        self._raise_accumulated_scan_errors()

    def assemble(self, trad2org, trust_root=None):
        """
        Concatenate the loading of TRAD2.org (if provided) and the results of the scan of the C++/xd files.
        At the same time, complete the second name (after the type) with the potential synonyms (TODO: review this?)
        ``trad2org`` may be None when there is no legacy .org file to fold in (all content lives in .cpp/.xd).
        """
        # Complete synonym lookup so that it works in all directions:
        lkp = {}
        for k, v in self.synos.items():
            lkp[k] = k
            for vv in v:
                lkp[vv] = k
        logger.debug("Synonyms are " + str(self.synos))

        # Load TRAD2.org first if provided. `source_locations=None`
        # auto-detects; the .org's sibling JSON (if any) provides the
        # backrefs, otherwise this is graceful no-source-locations.
        if trad2org is not None:
            torg = TRAD2Content.BuildContentFromTRAD2(trad2org, source_locations=None, trust_root=trust_root)
            # ... prepend its entries to the source-scanned ones ...
            self.data = torg.data + self.data

        # complete the synonyms for all entries (TRAD2.org-loaded + source-scanned)
        # (typically 'scheme_euler_explicit')
        for d in self.data:
            nam1, nam2 = d.name, d.synos[0]
            if nam1 in lkp or nam2 in lkp:
                # there are synonyms, put them with a '|' separator in the second slot after type:
                key = lkp[nam1] if nam1 in lkp else lkp[nam2]
                syn = {key, nam2}
                syn.update(self.synos[key])
                with contextlib.suppress(KeyError):
                    syn.remove(nam1)
                # Sort so the synonym ordering is deterministic across
                # runs — Python's set iteration depends on hash
                # randomization, which would otherwise leak into the
                # generated TRAD2 text and break cache-hash stability.
                d.synos = sorted(syn)

    def _convert_infos_to_source_locations(self) -> None:
        """Walk every block + attribute and replace `obj.info` (a
        `[absolute_path, lineno]` list produced by `BuildFromTab`)
        with a `SourceLocation` resolved against `self.projects`.

        Objects whose `info` is already a `SourceLocation` are left
        unchanged. Raises `ValueError` (via `relativize` or the
        ${TRUST_ROOT} expansion guard) when a path falls outside every
        configured project root.
        """
        for blk in self.data:
            self._convert_one_info(blk)
            for attr in blk.attrs:
                self._convert_one_info(attr)

    def _convert_one_info(self, obj) -> None:
        from trustify.core.source_location import SourceLocation, relativize

        if isinstance(obj.info, SourceLocation):
            return
        path_or_template, line = obj.info[0], int(obj.info[1])
        # During the transition the writer may still hand us a
        # "${TRUST_ROOT}/..."-style template (legacy convert_path_to_relative
        # output). Expand it against the trust project root first.
        if path_or_template.startswith("${TRUST_ROOT}"):
            root = self.projects.get("trust")
            if root is None:
                raise ValueError(
                    f"cannot expand '${{TRUST_ROOT}}' in {path_or_template!r}: "
                    f"no 'trust' project configured "
                    f"(projects={list(self.projects.keys())})"
                )
            path_or_template = root + path_or_template[len("${TRUST_ROOT}") :]
        project, rel = relativize(path_or_template, self.projects)
        obj.info = SourceLocation(project=project, path=rel, line=line)
