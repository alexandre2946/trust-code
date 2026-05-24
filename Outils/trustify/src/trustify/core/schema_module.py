"""Per-schema Python module wrapper.

Replaces the global ``ClassFactory`` singleton. Each loaded schema is
a Python module (loaded via ``trustify.core.misc_utilities.import_parser_module``);
this class wraps it with the lookup helpers consumers need.

One ``SchemaModule`` instance per loaded schema. Multiple schemas
coexist in a single process via ``sys.modules`` keyed by a hash of
the resolved parser-module path.
"""

from typing import Any


class SchemaModule:
    """Thin wrapper over a generated parser module.

    Holds the module object plus convenience accessors. The module's
    own attributes (``_SYNO_ORIG_NAME``, ``_all_constrain_base_pyd``,
    ``_all_constrain_base_parser``) are populated by the footer of the
    generated parser module at import time.
    """

    def __init__(self, module):
        self._module = module

    def __getattr__(self, name):
        """Delegate unknown attribute lookups to the underlying module.

        This allows ``schema.Coucou`` to work as ``schema.module.Coucou``,
        keeping call-sites that set ``self.mod = self._TRUG[slot]`` working
        after the migration from raw module to :class:`SchemaModule`.
        """
        # Avoid infinite recursion during __init__ before _module is set
        if name == "_module":
            raise AttributeError(name)
        return getattr(self._module, name)

    @property
    def module(self):
        """Direct access to the underlying Python module.

        Most consumers should prefer the named accessors below, but
        the doc generator etc. occasionally need attributes the
        wrapper does not expose (e.g. ``Objet_u``).
        """
        return self._module

    def get_pyd_class(self, name):
        """Get a pydantic class by name (lowercase).

        Args:
            name: The class name in lowercase (e.g. "foo" -> "Foo").

        Returns:
            The pydantic class.

        Raises:
            AttributeError: If the class does not exist in the module.
        """
        from trustify.core.naming import ToPydName

        return getattr(self._module, ToPydName(name))

    def get_parser_class(self, name):
        """Get a parser class by name (lowercase).

        Args:
            name: The class name in lowercase (e.g. "foo" -> "Foo_Parser").

        Returns:
            The parser class.

        Raises:
            AttributeError: If the class does not exist in the module.
        """
        from trustify.core.naming import ToParserName

        return getattr(self._module, ToParserName(name))

    def get_pyd_from_parser(self, parser_cls):
        """Given a parser class, return the corresponding pydantic class.

        Args:
            parser_cls: A parser class (e.g. Foo_Parser).

        Returns:
            The corresponding pydantic class (e.g. Foo).

        Raises:
            AttributeError: If the pydantic class does not exist.
        """
        from trustify.core.naming import _PARSER_SUFFIX

        pyd_name = parser_cls.__name__[: -len(_PARSER_SUFFIX)]
        return getattr(self._module, pyd_name)

    def get_parser_from_pyd(self, pyd_cls):
        """Given a pydantic class, return the corresponding parser class.

        Args:
            pyd_cls: A pydantic class (e.g. Foo).

        Returns:
            The corresponding parser class (e.g. Foo_Parser).

        Raises:
            AttributeError: If the parser class does not exist.
        """
        from trustify.core.naming import _PARSER_SUFFIX

        return getattr(self._module, pyd_cls.__name__ + _PARSER_SUFFIX)

    def exists(self, name):
        """Check whether a pydantic or parser class exists.

        Args:
            name: The class name in lowercase.

        Returns:
            True if either the pydantic or parser class exists.
        """
        from trustify.core.naming import ToParserName, ToPydName

        return hasattr(self._module, ToPydName(name)) or hasattr(self._module, ToParserName(name))

    def all_constrain_base_pyd(self):
        """Return the module's precomputed list of pydantic base classes.

        Returns:
            List of pydantic classes descending from Objet_u.
        """
        return self._module._all_constrain_base_pyd

    def all_constrain_base_parser(self):
        """Return the module's precomputed list of parser base classes.

        Returns:
            List of parser classes descending from ConstrainBase_Parser.
        """
        return self._module._all_constrain_base_parser

    def synonyms_for(self, syno):
        """Return the list of pydantic classes for a synonym.

        Args:
            syno: A synonym string.

        Returns:
            List of pydantic classes that have this synonym, or [] if unknown.
        """
        return self._module._SYNO_ORIG_NAME.get(syno, [])


def _build_synonyms(module):
    """Walk the generated module; return ``{syno: [orig_cls, ...]}``.

    Called from the generated module's footer at import time. Mirrors
    the legacy ``ClassFactory.BuildSynonymMap`` logic.
    """
    out: dict[str, list[Any]] = {}
    objet_u = getattr(module, "Objet_u", None)
    if objet_u is None:
        return out
    for v in vars(module).values():
        if isinstance(v, type) and hasattr(v, "_synonyms") and objet_u in v.__mro__:
            for s in v._synonyms.get(None, []):
                out.setdefault(s, []).append(v)
    return out


def _build_all_constrain_base_pyd(module):
    """List of pydantic classes that descend from ``Objet_u``.

    Pre-computed at module-init to avoid linear scans at parse time.
    """
    objet_u = getattr(module, "Objet_u", None)
    if objet_u is None:
        return []
    return [v for v in vars(module).values() if isinstance(v, type) and objet_u in v.__mro__]


def _build_all_constrain_base_parser(module):
    """List of parser classes that descend from ``ConstrainBase_Parser``."""
    from trustify.core.base import ConstrainBase_Parser

    return [v for v in vars(module).values() if isinstance(v, type) and issubclass(v, ConstrainBase_Parser)]
