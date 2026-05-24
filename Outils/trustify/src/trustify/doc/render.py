"""Per-keyword Markdown card renderer.

`render_keyword_card(cls, direct_parent, extras)` returns a
Doxygen-compatible Markdown block (heading + synonyms + inherits-from
+ description + parameter list). Pure function; ``extras`` is a
``trustify.doc.extras.ExtrasIndex`` used to resolve ``\\input{{name}}``
and ``\\includeimage{{name.ext}}`` directives that appear in
keyword docstrings.

Directive resolution order: ``\\input`` first (so the inlined markdown
can itself carry image directives), then ``\\includeimage``.
"""

import re

# ``\\input{{name}}`` — captures the markdown stem.
_RE_INPUT = re.compile(r"\\input{{([a-z0-9]+)}}")

# ``\\includeimage{{name.ext}}`` — captures basename + extension.
_RE_IMG = re.compile(r"\\includeimage{{([a-z0-9]+)\.([a-z0-9]+)}}")


def _resolve_input_directives(text, extras):
    """Replace every ``\\input{{name}}`` in ``text`` with the content of
    ``name.md`` resolved through ``extras``. Raises ``ValueError`` on
    miss (bubbles up to abort the build).
    """

    def sub(match):
        return extras.resolve_input(match.group(1))

    return _RE_INPUT.sub(sub, text)


def _resolve_image_directives(text, extras):
    """Replace every ``\\includeimage{{name.ext}}`` in ``text``
    with ``![name](figures/name.ext)``. Calls ``extras.resolve_image``
    for verification + used-marking; raises ``ValueError`` on miss.
    """

    def sub(match):
        base, ext = match.group(1), match.group(2)
        filename = f"{base}.{ext}"
        # Side effect: verify presence and mark used. Ignore return.
        extras.resolve_image(filename)
        return f"\n![{base}](figures/{filename})\n"

    return _RE_IMG.sub(sub, text)


def render_keyword_card(cls, direct_parent, extras) -> str:
    """Render a single keyword class as a Doxygen ``###`` subsection.

    See module docstring for output shape. ``extras`` is an
    ``ExtrasIndex`` used to resolve in-prose directives.
    """
    from trustify.core.misc_utilities import break_type

    name = cls.__name__.lower()
    s = f"### {name} {{#kw_{name}}}\n\n"

    # Synonyms — skip the legacy 'nul' placeholder.
    synos = cls._synonyms.get(None) or []
    if synos and synos != ["nul"]:
        s += "**Synonyms:** {}\n\n".format(", ".join(synos))

    # Inheritance link — omitted for root keywords and for the invisible
    # root classes Objet_u / Objet_lecture.
    if cls is not direct_parent and direct_parent.__name__ not in ("Objet_u", "Objet_lecture"):
        s += f"**Inherits from:** @ref kw_{direct_parent.__name__.lower()} \n\n"

    # Prose description. Resolve \input first so inlined markdown can
    # carry \includeimage directives that get resolved next.
    doc = cls.__doc__ or ""
    doc = _resolve_input_directives(doc, extras)
    rendered_lines = [_resolve_image_directives(line.lstrip(), extras) for line in doc.split("\n")]
    s += "\n".join(rendered_lines)

    # `List*` classes have no fields — stop here.
    if cls.__name__.startswith("List"):
        return s
    if not hasattr(cls, "model_fields"):
        return s

    # Parameter list.
    if len(cls._synonyms.keys()) > 1:
        s += "\n\n**Parameters:**\n\n"
    for attr, attr_synos in cls._synonyms.items():
        if attr is None:
            continue
        fld_nfo = cls.model_fields.get(attr)
        if fld_nfo is None:
            continue
        all_names = " | ".join([attr] + (attr_synos or []))
        opt, _, _ = break_type(fld_nfo.annotation)
        prefix = f"- **[{all_names}]** " if opt else f"- **{all_names}** "
        typ, desc = _extract_type_and_desc(fld_nfo)
        typ = re.sub(r":ref:`([^`]+)`", r"@ref kw_\1", typ)
        s += f"{prefix} (*type:* {typ}{_default_annotation(fld_nfo)}) {desc}\n\n"
    return s


def _default_annotation(fld_nfo) -> str:
    """Return a `` , *default:* `X` `` fragment when the field declares a
    documented default (``json_schema_extra['trust_default']``, set by the
    code generator for ``chaine(into=[...],default="X")`` attributes), else
    an empty string. Display-only — it does not reflect the pydantic runtime
    default.
    """
    extra = getattr(fld_nfo, "json_schema_extra", None)
    if isinstance(extra, dict) and "trust_default" in extra:
        return f", *default:* `{extra['trust_default']}`"
    return ""


def _extract_type_and_desc(fld_nfo):
    """Convert a Pydantic field annotation into a human-readable type
    string + return the field's description.
    """
    from trustify.core.base import Builtin_Parser
    from trustify.core.misc_utilities import break_type

    ann = fld_nfo.rebuild_annotation()
    if Builtin_Parser.IsBuiltin(ann):
        pars_inst = Builtin_Parser.InstanciateFromBuiltin(ann)
        nice_typ = pars_inst.getFormattedType()
    else:
        _, typ, _ = break_type(ann)
        t = typ[0].__name__.lower()
        nice_typ = "objet_u" if t == "objet_u" else f":ref:`{typ[0].__name__.lower()}`"
    return nice_typ, fld_nfo.description
