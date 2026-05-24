#!/usr/bin/env python3
"""
Post-process Doxygen navtree files to flatten the root project-name entry.

Problem
-------
Doxygen wraps the entire navigation tree under a single root node labelled
with the project name (here "TRUST").  In the browser this means every
section of the documentation sits one level too deep — the user has to expand
the "TRUST" node before seeing Quick Start, User Guide, etc.

What this script does
---------------------
1. **Flatten navtreedata.js** — keep the root node as a clickable leaf
   pointing to the main page (so the project name remains in the tree), but
   promote all of its former children to top-level entries.

2. **Fix navtreeindex*.js breadcrumbs** — Doxygen encodes navigation paths as
   integer arrays (e.g. ``[0,3,1]`` means "root → child 3 → child 1").  After
   hoisting all children by one level every path that previously started at
   root's children now lives at the top level, so the first index of each
   breadcrumb must be incremented by 1.

3. **Patch navtree.js** — Doxygen's navtree JavaScript hard-codes
   ``breadcrumbs.unshift(0)`` to always route navigation through the root
   node first.  After flattening this extra step is wrong and must be removed.

Usage
-----
    python3 flatten_navtree.py <doxygen-html-dir>
"""

import glob
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rewrite(path: str, transform) -> None:
    """Read *path*, apply *transform* to its text content, write it back."""
    p = Path(path)
    p.write_text(transform(p.read_text(encoding="utf-8")), encoding="utf-8")


def _find_matching_close_bracket(content: str, open_pos: int) -> int:
    """Return the index of the ``]`` matching the ``[`` at *open_pos*.

    Walks forward counting bracket depth.  Skips characters inside JavaScript
    string literals (both ``'``-quoted and ``"``-quoted, with backslash
    escapes) so a stray ``[`` or ``]`` inside a page title or symbol name
    (e.g. ``'operator[]'``) does not desynchronise the depth counter.
    Raises ``ValueError`` if no matching ``]`` is found before EOF.
    """
    assert content[open_pos] == "[", "open_pos must point at a '[' character"
    depth = 1
    i = open_pos + 1
    n = len(content)
    while i < n:
        c = content[i]
        if c == '"' or c == "'":
            # Skip the string literal entirely (handle backslash-escapes).
            quote = c
            i += 1
            while i < n and content[i] != quote:
                if content[i] == "\\" and i + 1 < n:
                    i += 2
                    continue
                i += 1
            i += 1  # step past the closing quote
            continue
        if c == "[":
            depth += 1
        elif c == "]":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError(f"no matching ']' found starting at offset {open_pos}")


def flatten_navtreedata(content: str) -> str:
    """Restructure navtreedata.js: keep root as a leaf, hoist its children.

    navtreedata.js has the form::

        var NAVTREE = [
          [ "ProjectName", "index.html", [
            [ "Section A", "pageA.html", null ],
            [ "Section B", "pageB.html", [...] ],
            ...
          ]]
        ];

    After this transformation it becomes::

        var NAVTREE = [
          [ "ProjectName", "index.html", null ],   ← leaf, no children
          [ "Section A", "pageA.html", null ],      ← promoted to top level
          [ "Section B", "pageB.html", [...] ],
          ...
        ];
    """
    start = content.index("var NAVTREE")

    # Locate the opening of the root entry: ["Name", "url", [
    m = re.search(r'\[\s*"([^"]+)"\s*,\s*"([^"]+)"\s*,\s*\[', content[start:])
    name, url = m.group(1), m.group(2)

    # Walk forward from the opening bracket of the children array to find
    # where the children array ends.  The walker is string-aware so a stray
    # ``[`` or ``]`` inside a title doesn't desynchronise the depth counter.
    open_bracket = start + m.end() - 1  # position of the '[' captured by the regex
    close_bracket = _find_matching_close_bracket(content, open_bracket)
    children_start = open_bracket + 1
    children_content = content[children_start:close_bracket].strip()
    i = close_bracket + 1  # one past the closing ']'

    # Find the end of the full NAVTREE statement (closing "];").
    navtree_end = content.index("];", i) + 2

    new_navtree = (
        f'var NAVTREE =\n[\n  [ "{name}", "{url}", null ],\n  {children_content}\n];\n'
    )
    return content[:start] + new_navtree + content[navtree_end:].lstrip()


def update_breadcrumbs(content: str) -> str:
    """Increment the first index of every breadcrumb array in navtreeindex*.js.

    Doxygen stores navigation breadcrumbs as ``"<href>":[<int-list>]``
    entries, e.g.::

        "pageA.html":[0,2]    →  root (index 0) → child 2

    After flattening, what was formerly root's child at position N is now
    at top-level position N+1 (because the root node itself still occupies
    position 0 as a leaf).  Incrementing the first index of every breadcrumb
    by 1 corrects all paths without changing the relative structure.

    The regex is anchored to the ``"<href>":[...]`` shape (where ``<href>``
    is a Doxygen-generated ``*.html`` URL, optionally with an ``#anchor``)
    so unrelated ``:[N,...]`` patterns elsewhere in the file are not
    touched — audit H4.
    """
    def inc_first(m: re.Match) -> str:
        href, nums = m.group(1), m.group(2).split(",")
        nums[0] = str(int(nums[0]) + 1)
        return f'"{href}":[{",".join(nums)}]'

    return re.sub(
        r'"([^"]+\.html(?:#[^"]*)?)"\s*:\s*\[(\d+(?:,\d+)*)\]',
        inc_first,
        content,
    )


# Legacy Doxygen versions (≤ 1.9.x) prepended 0 to every breadcrumb via the
# line below, assuming sidebar navigation always started at NAVTREE[0].
# Doxygen 1.16+ removed this and stores absolute breadcrumb paths directly
# in navTreeSubIndices.  Our flatten_navtreedata() hoists root's children
# up one level, so the legacy unshift logic would push the sidebar walk
# one level too deep and silently break highlighting.
_LEGACY_UNSHIFT_RE = re.compile(
    r"\bo\.breadcrumbs\.unshift\s*\(\s*0\s*\)\s*;[^\n]*\n",
)


def patch_navtree_js(content: str) -> str:
    """Strip Doxygen's legacy ``breadcrumbs.unshift(0)`` line if present.

    The line was emitted by Doxygen ≤ 1.9.x and removed in Doxygen 1.16+.
    The pattern is a regex (not a literal-string ``.replace``) so trivial
    upstream changes — different surrounding comment, different whitespace,
    LF vs CRLF — still match.  Audit H3.

    Idempotent: on Doxygen 1.16+ the line is absent and this is a no-op.
    """
    return _LEGACY_UNSHIFT_RE.sub("", content)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def reroute_bibliography(content: str) -> str:
    """Redirect the Bibliography navtree slot from References.html to citelist.html.

    Doxygen generates the bibliography as citelist.html (a special internal page)
    that it never includes in the navtree automatically.  The @page References
    wrapper is used only to obtain a navtree slot; this function replaces its URL
    with citelist.html so the sidebar entry lands directly on the bibliography
    content rather than an intermediate description page.

    Safe to apply to any file: str.replace() is a no-op when the target string
    is absent, so calling this on a file that contains no References.html
    reference leaves the file unchanged.
    """
    return content.replace('"References.html"', '"citelist.html"')


def main() -> None:
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <doxygen-html-dir>", file=sys.stderr)
        sys.exit(1)

    html_dir = sys.argv[1]

    # 1. Flatten the NAVTREE structure in navtreedata.js.
    _rewrite(f"{html_dir}/navtreedata.js", flatten_navtreedata)

    # 2. Fix breadcrumb indices in every navtreeindex*.js shard.
    #    Doxygen splits the index into multiple files for large projects.
    for path in sorted(glob.glob(f"{html_dir}/navtreeindex*.js")):
        _rewrite(path, update_breadcrumbs)

    # 3. Remove the unshift(0) call from the navtree JavaScript engine.
    _rewrite(f"{html_dir}/navtree.js", patch_navtree_js)

    # 4. Redirect the Bibliography navtree entry to citelist.html.
    #
    #    Doxygen ≤ 1.9.x stores the full navtree inline in navtreedata.js.
    #    Doxygen ≥ 1.16.x lazy-loads children into a separate index.js.
    #    Both files are rewritten here; reroute_bibliography() is a no-op
    #    when "References.html" is absent, so rewriting a file that does
    #    not contain it is always safe.
    _rewrite(f"{html_dir}/navtreedata.js", reroute_bibliography)
    idx_js = Path(html_dir) / "index.js"
    if idx_js.exists():
        _rewrite(str(idx_js), reroute_bibliography)
    for path in sorted(glob.glob(f"{html_dir}/navtreeindex*.js")):
        _rewrite(path, reroute_bibliography)


if __name__ == "__main__":
    main()
