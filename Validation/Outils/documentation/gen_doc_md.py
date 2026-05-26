#!/usr/bin/env python3
"""
Generate Doxygen-compatible Markdown API documentation for the trustutils package.

Pipeline position
-----------------
This script is called by ``make md`` in ``Validation/Outils/documentation/``,
which is itself called by ``make trustutils_doc`` in ``docs/``.
It sits entirely outside the trustify / keyword-reference pipeline and
documents the *Python* helpers that users write validation forms with.

How it works
------------
The script uses Python's ``ast`` module to parse each trustutils source file
*statically* — it never imports the package at runtime.  This means it has
zero runtime dependencies (no TRUST environment needed) and works regardless
of which optional packages (numpy, matplotlib, ...) are installed.

For each module listed in ``MODULES``:
1.  Parse the source file with ``ast``.
2.  Extract module, class, method, and function docstrings.
3.  Convert NumPy-style docstring sections (Parameters, Returns, ...) to
    Markdown tables and fenced code blocks.
4.  Emit a Doxygen ``@page`` Markdown file.

Additionally, an index page (``validation_form_api.md``) is generated listing
every module page via ``@subpage``.

Output
------
One ``<page_id>.md`` file per documented module, plus ``validation_form_api.md``.
Files are written to *out_dir* and later copied into the Doxygen source tree
by ``make copy_trustutils`` in ``docs/Makefile``.

Usage
-----
    gen_doc_md.py <trustutils_dir> <out_dir>

Arguments
---------
trustutils_dir
    Path to the ``trustutils/`` package root
    (typically ``$TRUST_ROOT/Validation/Outils/trustutils``).
out_dir
    Directory where the generated ``.md`` files are written.
"""

import sys
import os
import re
import ast
from pathlib import Path

# ---------------------------------------------------------------------------
# Module registry
# ---------------------------------------------------------------------------
# Each entry: (path_relative_to_pkg_root, doxygen_page_id, display_name)
#
# Page IDs MUST start with "TRUSTUTILS_" so that patch_filter.py routes
# them to the "Trustutils" search tab rather than the generic Documentation
# tab.  See the RULES table in docs/scripts/patch_filter.py.

MODULES = [
    ("files.py",                       "TRUSTUTILS_files",            "trustutils.files"),
    ("jupyter/run.py",                 "TRUSTUTILS_jupyter_run",      "trustutils.jupyter.run"),
    ("jupyter/plot.py",                "TRUSTUTILS_jupyter_plot",     "trustutils.jupyter.plot"),
    ("jupyter/filelist.py",            "TRUSTUTILS_jupyter_filelist", "trustutils.jupyter.filelist"),
    ("jupyter/widget.py",              "TRUSTUTILS_jupyter_widget",   "trustutils.jupyter.widget"),
    ("visitutils/tools_for_visit.py",  "TRUSTUTILS_visitutils",       "trustutils.visitutils.tools_for_visit"),
    ("visitutils/message.py",          "TRUSTUTILS_visitutils_msg",   "trustutils.visitutils.message"),
]

# ---------------------------------------------------------------------------
# AST helpers
# ---------------------------------------------------------------------------

def _get_docstring(node):
    """Return the docstring of an AST node, or None if there is none.

    A docstring is the first statement of a module, class, or function body
    when that statement is a bare string expression.  Python 3.8+ represents
    string literals as ``ast.Constant``; older versions used ``ast.Str``.
    ``ast.Str`` is removed in 3.14. 
    We choose to drop pre-3.8 compat (pydantic + trustify require 3.9+ anyway)
    """
    if (node.body
            and isinstance(node.body[0], ast.Expr)
            and isinstance(node.body[0].value, (ast.Constant))):
        val = node.body[0].value
        return val.value if isinstance(val, ast.Constant) else val.s
    return None


def _get_signature(func_node):
    """Reconstruct a human-readable argument list from an ``ast.FunctionDef``.

    Type annotations are intentionally omitted — they tend to be verbose and
    add little value in this documentation context.  Default values are
    represented as ``...`` rather than being evaluated (the AST does not
    preserve the original source text of defaults).

    Handles positional-only arguments (``/`` syntax, Python 3.8+, via the
    ``posonlyargs`` attribute), regular arguments, ``*args``, keyword-only
    arguments, and ``**kwargs``.
    """
    args = func_node.args
    parts = []

    # posonlyargs was added in Python 3.8; fall back to empty list on older ASTs.
    all_pos = list(getattr(args, 'posonlyargs', [])) + list(args.args)
    n_defaults = len(args.defaults)
    # defaults align to the *end* of all_pos, so the first argument with a
    # default is at index (len(all_pos) - n_defaults).
    offset = len(all_pos) - n_defaults

    for i, arg in enumerate(all_pos):
        parts.append(f"{arg.arg}=..." if i >= offset else arg.arg)

    if args.vararg:
        parts.append(f"*{args.vararg.arg}")
    for kw in args.kwonlyargs:
        parts.append(kw.arg)
    if args.kwarg:
        parts.append(f"**{args.kwarg.arg}")

    return f"({', '.join(parts)})"


def _anchor(text):
    """Derive a safe Doxygen anchor identifier from an arbitrary string.

    Doxygen anchors (used in ``{#anchor}`` syntax) must consist of
    alphanumeric characters and underscores only.  This function lowercases
    the input and replaces every other character with ``_``.
    """
    return re.sub(r'[^a-zA-Z0-9]', '_', text.lower())


# ---------------------------------------------------------------------------
# NumPy-style docstring parser
# ---------------------------------------------------------------------------
# NumPy docstring format uses section headers followed by a line of dashes:
#
#   Parameters
#   ----------
#   name : type
#       Description.
#
# _SECTION_RE matches known section names; _UNDERLINE_RE matches the
# separator line.  Together they detect section boundaries.

_SECTION_RE = re.compile(
    r'^(Parameters|Returns?|Attributes|Raises|Notes?|Examples?|'
    r'See Also|References|Yields?|Warns?)\s*$'
)
_UNDERLINE_RE = re.compile(r'^[-=]{2,}\s*$')


def _parse_numpy_docstring(doc):
    """Split a NumPy-style docstring into a dict of named sections.

    The key ``'summary'`` holds everything before the first section header.
    All other keys are section names as they appear in the docstring
    (e.g. ``'Parameters'``, ``'Returns'``).

    Plain docstrings with no section headers are returned as
    ``{'summary': <full text>}``.
    """
    if not doc:
        return {"summary": ""}

    # Normalise indentation: expand tabs, then strip the common leading indent
    # so that section headers appear at column 0 regardless of how the
    # docstring was indented in the source.
    raw_lines = doc.expandtabs().splitlines()
    stripped  = [l for l in raw_lines if l.strip()]
    indent    = min((len(l) - len(l.lstrip()) for l in stripped), default=0)
    lines     = [l[indent:] for l in raw_lines]

    sections      = {}
    current       = "summary"
    current_lines = []
    i = 0
    while i < len(lines):
        line      = lines[i]
        next_line = lines[i + 1] if i + 1 < len(lines) else ""
        # A section boundary is: a recognised section name followed immediately
        # by a line of dashes (or equals signs) of length >= 2.
        if _SECTION_RE.match(line.strip()) and _UNDERLINE_RE.match(next_line.strip()):
            sections[current] = "\n".join(current_lines).strip()
            current       = line.strip()
            current_lines = []
            i += 2          # skip both the header line and the underline
        else:
            current_lines.append(line)
            i += 1
    sections[current] = "\n".join(current_lines).strip()
    return sections


def _format_param_block(text):
    """Convert a NumPy ``Parameters`` (or ``Attributes`` / ``Raises``) block to a Markdown table.

    Input format (one parameter per entry)::

        name : type
            Multi-line description
            that may span several lines.
        other_name : type
            Another description.

    Output: a GFM table with columns Name | Type | Description.

    If parsing fails (unrecognised format), the raw text is returned unchanged
    as a fallback so no information is lost.
    """
    if not text.strip():
        return ""
    lines     = text.splitlines()
    rows      = []
    cur_name  = cur_type = ""
    cur_desc: list = []
    # A parameter definition line: starts at column 0, contains " : ".
    param_re  = re.compile(r'^(\w[\w., *\[\]]*?)\s*:\s*(.*?)\s*$')

    def _flush():
        """Append the current parameter to rows if one is being accumulated."""
        if cur_name:
            rows.append((cur_name, cur_type, " ".join(cur_desc).strip()))

    for line in lines:
        if not line.startswith((" ", "\t")):
            # Unindented line — try to parse as a new parameter definition.
            m = param_re.match(line)
            if m:
                _flush()
                cur_name, cur_type = m.group(1).strip(), m.group(2).strip()
                cur_desc = []
                continue
        if cur_name:
            # Indented continuation of the current parameter's description.
            cur_desc.append(line.strip())
    _flush()

    if not rows:
        return text   # fallback: return raw text unmodified

    out = "| Name | Type | Description |\n|---|---|---|\n"
    for name, typ, desc in rows:
        out += f"| `{name}` | {typ or '—'} | {desc} |\n"
    return out


def _format_returns_block(text):
    """Convert a NumPy ``Returns`` or ``Yields`` block to inline Markdown.

    Input format::

        type_name : optional_type
            Description.

    Each return value is emitted as ``*type_name* — description`` on its own
    line.  If parsing yields no results the raw text is returned unchanged.
    """
    if not text.strip():
        return text
    lines   = text.splitlines()
    results = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip() and not line.startswith((" ", "\t")):
            # Unindented line — start of a return value entry.
            typ, _, rest = line.partition(":")
            desc_parts = [rest.strip()] if rest.strip() else []
            i += 1
            # Gather indented continuation lines.
            while i < len(lines) and lines[i].startswith((" ", "\t")):
                desc_parts.append(lines[i].strip())
                i += 1
            full_desc = " ".join(desc_parts).strip()
            results.append(f"*{typ.strip()}*" + (f" — {full_desc}" if full_desc else ""))
        else:
            i += 1
    return "\n".join(results) if results else text


def _docstring_to_md(doc):
    """Convert a (possibly NumPy-style) docstring to Markdown.

    Handles the following NumPy sections in a fixed display order:
      - Summary (everything before the first section header)
      - Parameters / Attributes  → Markdown table
      - Returns / Yields         → inline type–description pairs
      - Raises / Warns           → Markdown table
      - Note / Notes             → block-quote
      - Example / Examples       → fenced Python code block

    Any section not listed above (e.g. ``See Also``, ``References``) is
    silently dropped.  Plain docstrings with no sections pass through as-is.
    """
    if not doc:
        return ""
    sections = _parse_numpy_docstring(doc)
    parts    = []

    summary = sections.get("summary", "").strip()
    if summary:
        parts.append(summary)

    for key in ("Parameters", "Attributes"):
        val = sections.get(key, "").strip()
        if val:
            parts.append(f"\n**{key}:**\n\n{_format_param_block(val)}")

    for key in list(sections):
        if re.match(r'^Returns?$|^Yields?$', key):
            val = sections[key].strip()
            if val:
                parts.append(f"\n**{key}:** {_format_returns_block(val)}")

    for key in ("Raises", "Warns"):
        val = sections.get(key, "").strip()
        if val:
            parts.append(f"\n**{key}:**\n\n{_format_param_block(val)}")

    for key in ("Note", "Notes"):
        val = sections.get(key, "").strip()
        if val:
            parts.append(f"\n> **Note:** {val}")

    for key in list(sections):
        if re.match(r'^Examples?$', key):
            val = sections[key].strip()
            if val:
                parts.append(f"\n**Example:**\n\n```python\n{val}\n```")

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Markdown document builders
# ---------------------------------------------------------------------------

def _doc_function_md(func_node, heading="###", parent_anchor=""):
    """Generate a Markdown section for a single function or method AST node.

    The heading contains only the bare name (with anchor); the full call
    signature is displayed immediately below as a fenced Python code block so
    that it gets proper syntax highlighting without crowding the heading.

    Parameters
    ----------
    func_node     : ast.FunctionDef or ast.AsyncFunctionDef
    heading       : str
        Markdown heading level (``###`` for both methods and module-level
        functions).
    parent_anchor : str
        Anchor prefix of the containing class or module, used to build a
        unique nested anchor (``<parent>_<funcname>``).

    Returns
    -------
    str
        Markdown text for the function, terminated with a newline.
    """
    name = func_node.name
    sig  = _get_signature(func_node)
    # Nested anchors prevent collisions between same-named methods in different
    # classes within the same @page.
    anc  = f"{parent_anchor}_{_anchor(name)}" if parent_anchor else _anchor(name)
    doc  = _get_docstring(func_node) or ""
    s    = f"{heading} {name} {{#{anc}}}\n\n"
    s   += f"```python\n{name}{sig}\n```\n\n"
    s   += _docstring_to_md(doc) or "*No documentation.*"
    return s + "\n"


def _doc_class_md(cls_node, module_anchor=""):
    """Generate a Markdown section for a class AST node and its public methods.

    The section structure is::

        ## ClassName {#anchor}

        <class docstring>

        **Constructor:** `ClassName(args)` (only if __init__ has a docstring)

        ### method_name {#anchor_method_name}
        ```python
        method_name(args)
        ```
        <docstring>

        ---

        ### next_method {#anchor_next_method}
        ...

    The ``#### Methods`` intermediate header is intentionally omitted — the
    ``###`` method headings are visually sufficient to signal the transition
    from the class description to its members.

    Private methods (name starting with ``_``) and dunder methods other than
    ``__init__`` are omitted.  Methods are sorted alphabetically.

    Parameters
    ----------
    cls_node      : ast.ClassDef
    module_anchor : str
        Anchor prefix of the containing module, used to build a unique
        class anchor (``<module>_<classname>``).
    """
    name = cls_node.name
    anc  = f"{module_anchor}_{_anchor(name)}" if module_anchor else _anchor(name)
    s    = f"## {name} {{#{anc}}}\n\n"

    doc = _get_docstring(cls_node)
    if doc:
        s += _docstring_to_md(doc) + "\n\n"

    # Include __init__ only when it carries its own docstring (the class-level
    # docstring already covers construction in the common case).
    init = next(
        (n for n in cls_node.body
         if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
         and n.name == "__init__"),
        None,
    )
    if init and _get_docstring(init):
        sig  = _get_signature(init)
        s   += f"**Constructor:** `{name}{sig}`\n\n"
        s   += _docstring_to_md(_get_docstring(init)) + "\n\n"

    # Public methods: no leading underscore (excludes dunder and private).
    methods = [
        n for n in cls_node.body
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        and not n.name.startswith("_")
    ]
    for i, m in enumerate(sorted(methods, key=lambda n: n.name)):
        if i > 0:
            s += "\n---\n\n"
        s += _doc_function_md(m, heading="###", parent_anchor=anc) + "\n"

    return s


def _doc_module_md(rel_path, page_id, module_name, source_text):
    """Generate a complete Doxygen ``@page`` Markdown document for one source file.

    The page structure is::

        @page <page_id> <module_name>

        <module docstring>

        ## ClassName {#anchor}
        <class docstring>
        ### method_name {#anchor}
        ```python
        method_name(args)
        ```
        <docstring>
        ---
        ## NextClass {#anchor}
        ...

        ---

        ## Functions {#<anchor>_functions}
        ### function_name {#anchor}
        ```python
        function_name(args)
        ```
        <docstring>

    Classes and functions are each sorted alphabetically.  Only public symbols
    (names not starting with ``_``) are included.

    Parameters
    ----------
    rel_path    : str
        Path of the source file relative to the trustutils package root
        (used only for the error message on parse failure).
    page_id     : str
        Doxygen page identifier, e.g. ``TRUSTUTILS_files``.
    module_name : str
        Human-readable module name shown as the page title.
    source_text : str
        Full Python source text to parse.
    """
    try:
        tree = ast.parse(source_text)
    except SyntaxError as e:
        # Emit a minimal error page rather than aborting the whole build.
        return f"@page {page_id} {module_name}\n\n*Parse error: {e}*\n"

    s = f"@page {page_id} {module_name}\n\n"

    mod_doc = _get_docstring(tree)
    if mod_doc:
        s += _docstring_to_md(mod_doc) + "\n\n"

    # Use the lowercased page_id as the anchor prefix for all sub-anchors on
    # this page to avoid collisions between same-named symbols across pages.
    mod_anchor = _anchor(page_id)

    # Top-level public classes, sorted alphabetically.
    classes = [
        n for n in tree.body
        if isinstance(n, ast.ClassDef) and not n.name.startswith("_")
    ]
    for i, cls_node in enumerate(sorted(classes, key=lambda n: n.name)):
        if i != 0:
            s += "\n---\n\n"
        s += _doc_class_md(cls_node, module_anchor=mod_anchor)

    # Top-level public functions (not nested inside classes), sorted alphabetically.
    funcs = [
        n for n in tree.body
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        and not n.name.startswith("_")
    ]
    if funcs:
        if classes:
            s += "\n---\n\n"
        s += f"## Functions {{#{mod_anchor}_functions}}\n\n"
        for i, func_node in enumerate(sorted(funcs, key=lambda n: n.name)):
            if i > 0:
                s += "\n---\n\n"
            s += _doc_function_md(func_node, heading="###", parent_anchor=mod_anchor) + "\n"

    return s


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def generate_docs(trustutils_dir, out_dir):
    """Generate all documentation pages for the modules listed in MODULES.

    For each module: parse → convert → write ``<page_id>.md``.
    Then write the ``validation_form_api.md`` index page.

    Missing source files produce a warning on stderr and are skipped (the
    index page will simply not list them), so a partial trustutils installation
    does not abort the whole doc build.
    """
    trustutils_dir = Path(trustutils_dir).resolve()
    out_dir        = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    page_ids = []
    for rel_path, page_id, module_name in MODULES:
        src_file = trustutils_dir / rel_path
        if not src_file.exists():
            print(f"WARNING: {src_file} not found, skipping.", file=sys.stderr)
            continue
        source_text = src_file.read_text(encoding="utf-8")
        content     = _doc_module_md(rel_path, page_id, module_name, source_text)
        out_file    = out_dir / f"{page_id.lower()}.md"
        out_file.write_text(content, encoding="utf-8")
        print(f"  Written: {out_file}")
        page_ids.append((page_id, module_name))

    # Index page — this file is copied to
    # docs/content/dev_corner/validation_form_api.md and declares the
    # @page that DevCorner links to via @subpage.
    index  = "@page Dev_ValidationFormAPI Validation form API guide\n\n"
    index += (
        "The `trustutils` Python package provides utilities for writing TRUST validation "
        "forms in Jupyter notebooks. It is located in `Validation/Outils/trustutils/`.\n\n"
    )
    for pid, _mname in page_ids:
        index += f"- @subpage {pid}\n"

    index_file = out_dir / "validation_form_api.md"
    index_file.write_text(index, encoding="utf-8")
    print(f"  Written: {index_file}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: gen_doc_md.py <trustutils_dir> <out_dir>")
        sys.exit(1)
    generate_docs(sys.argv[1], sys.argv[2])
