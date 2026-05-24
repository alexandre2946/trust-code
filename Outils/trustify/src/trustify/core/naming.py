"""Pure string-mangling helpers for trustify's class naming convention.

These are the only pieces of the legacy ``ClassFactory`` that survive
its removal — they have no state, just convert between the three name
forms used in the codebase:

- TRUST keyword:  ``read_med``        (free-form, lowercase)
- Pydantic class: ``Read_med``        (lowercase + capitalize)
- Parser class:   ``Read_med_Parser`` (pydantic name + suffix)
"""

_PARSER_SUFFIX = "_Parser"


def ToPydName(name):
    """Convert a TRUST keyword to its pydantic class name."""
    return name.lower().capitalize()


def ToParserName(name):
    """Convert a TRUST keyword to its parser class name."""
    return ToPydName(name) + _PARSER_SUFFIX
