"""From a TRAD2 file, generates:
 - the Pydantic schema (as a generated Python file containing one Pydantic class per TRUST keyword)
 - the parsing classes (as a generated Python file, containing one Python class per TRUST keyword)

Authors: G Sutra, A Bruneton
"""

import re
import textwrap

import trustify.core.trad2_utilities as tu
from trustify.core.misc_utilities import logger, pretty_error
from trustify.core.naming import ToParserName, ToPydName

# Captures the (quoted) value of the optional `,default="X"` clause on a
# `chaine(into=[...],default="X")` type. The group keeps the surrounding
# quotes so the token can be matched directly against the (also-quoted)
# enum choices and dropped verbatim into the generated metadata.
_CHAINE_DEFAULT_RE = re.compile(r',default=("[^"\]]*")\)\Z')

# Same for the integer `entier(into=[...],default=N)` / `entier(max=N,default=M)`
# forms — the default is a bare or quoted int.
_INT_DEFAULT_RE = re.compile(r',default=(-?\d+|"-?\d+")\)\Z')

# Extract the `min=` / `max=` bounds from an `entier(min=...,max=...)` type.
_ENTIER_MIN_RE = re.compile(r"[(,]min=(-?\d+)[,)]")
_ENTIER_MAX_RE = re.compile(r"[(,]max=(-?\d+)[,)]")

################################################################


def valid_variable_name(s):
    """Make a valid variable name from any str. Useful to avoid variable name colinding with Python
    keywords ('lambda' for example ...).

    Raises `ValueError` when `s` contains any character outside
    `[A-Za-z0-9_]` (audit 4.8). TRUST keyword names are documented as
    that character class only — non-conforming names used to be
    silently stripped (`étude` → `tude`), which could collide with a
    sibling attribute and lose the offending name without any
    diagnostic. Failing loudly points the user back at the XD tag in
    the source that needs fixing.
    """
    import keyword
    import re

    bad = sorted(set(re.findall(r"[^0-9a-zA-Z_]", s)))
    if bad:
        raise ValueError(
            f"TRUST keyword/attribute name {s!r} contains illegal character(s) "
            f"{''.join(bad)!r}; names must match [A-Za-z0-9_]+. Fix the XD tag "
            f"in the source."
        )
    # prefix leading number with 'i'
    s = re.sub(r"^([0-9])", r"i\1", s)
    # suffix reserved keywords with '_'
    if keyword.iskeyword(s):
        s += "_"
    return s


def substitute_nl(s, repl):
    s = s.replace("NL1", repl)
    s = s.replace("NL2", repl * 2)
    return s


def format_docstring(description):
    import itertools

    docstring = f'r"""\n{description}\n"""'
    docstring = substitute_nl(docstring, "\n")
    docstring = docstring.splitlines()
    docstring = [textwrap.wrap(line, width=90) if line else [""] for line in docstring]
    docstring = list(itertools.chain.from_iterable(docstring))
    docstring = [f"    {line.lstrip()}" for line in docstring]
    return docstring


def is_base_or_deriv(cls_nam):
    """The TRUST keywords ending with _base or _deriv usually needs to have their type explicitely read in the dataset"""
    return cls_nam.endswith(("_deriv", "_base")) or cls_nam == "class_generic"


def _info_to_triple(info):
    """Convert a TRAD2Block/Attr ``info`` value to the
    ``(project, path, line)`` triple form embedded in the generated
    parser's ``_infoMain`` / ``_infoAttr`` class variables.

    - ``SourceLocation(project, path, line)`` → ``(project, path, line)``.
    - Legacy ``[abs_path, line]`` (from a no-projects ``BuildFromTab``
      path, e.g. the "simple" hand-written test fixture which has no
      ``source_locations.json``) → ``("", abs_path, line)``. Empty
      project is the signal for "no project resolution available";
      consumers that need an absolute path fall back to using `path`
      as-is.
    - Anything else → ``("", "", -1)``.
    """
    from trustify.core.source_location import SourceLocation

    if isinstance(info, SourceLocation):
        return (info.project, info.path, info.line)
    if isinstance(info, (list, tuple)) and len(info) == 2:
        return ("", info[0], info[1])
    return ("", "", -1)


def _check_nonempty_enum(choices_str: str, attr) -> None:
    """Raise a clear diagnostic when an ``into=[...]`` enum is still empty
    at code-generation time.

    ``_VALID_TYPE_RE`` deliberately accepts ``chaine(into=[])`` because the
    ``dico`` shorthand starts empty and is filled in by ``XD_ADD_DICO``
    lines. If it reaches codegen still empty, the keyword has no usable
    values — emit an actionable error rather than the invalid ``Literal[]``
    / ``default=`` (a SyntaxError on import) or the cryptic ``int('')``
    ValueError the naive path would otherwise produce.
    """
    if choices_str.strip():
        return
    _, _at_path, _at_line = _info_to_triple(attr.info)
    raise Exception(
        pretty_error(
            _at_path,
            _at_line,
            f"attribute {attr.name!r} declares an empty enumeration "
            f"({attr.type!r}); an enum (or a 'dico') must list at least one "
            "value. For a 'dico', add an 'XD_ADD_DICO' line; otherwise put "
            "the allowed values inside into=[...].",
        )
    )


def _chaine_declared_default(attr_type, choices_str, attr):
    """Return the (quoted) value of a ``chaine(into=[...],default="X")``
    clause, or ``None`` when the clause is absent.

    This value is recorded as display-only Field metadata, not as the
    runtime default — the latter stays the first choice (REQ) / None
    (OPT) so parsing and round-trip are unaffected. Membership of ``X``
    in the listed choices is validated here (not at scan time), because
    ``dico`` enums are still empty when the scanner runs
    `_validate_type_string` and only get their values from later
    ``XD_ADD_DICO`` lines.

    Assumes ``_check_nonempty_enum`` already ran, so ``choices_str`` lists
    at least one value.
    """
    m = _CHAINE_DEFAULT_RE.search(attr_type)
    if m is None:
        return None
    default_tok = m.group(1)
    choices = [c.strip() for c in choices_str.split(",") if c.strip()]
    if default_tok not in choices:
        _, _at_path, _at_line = _info_to_triple(attr.info)
        raise Exception(
            pretty_error(
                _at_path,
                _at_line,
                f"attribute {attr.name!r}: default {default_tok} in type "
                f"{attr.type!r} is not one of the declared values {choices}.",
            )
        )
    return default_tok


def _entier_declared_default(attr_type):
    """Return the int value of an ``entier(...,default=N)`` clause, or
    ``None`` when absent. Caller validates the value against the choices
    / bound. Like the chaine variant, this is display-only metadata."""
    m = _INT_DEFAULT_RE.search(attr_type)
    if m is None:
        return None
    return int(m.group(1).strip('"'))


def generate_attribute_synos(block, all_blocks):
    """Generate a dictionary containing:
    - under the 'None' key : the list of synonyms for the keyword itself ;
    - and for **all** the attributes of a class (even the inherited ones) their synonyms.
    Important: the keys (i.e. the attribute names) are ordered in a specific logic:
        - attribute of the child class first
        - and then inherited attributes from the mother class ...
    """
    import sys
    from copy import deepcopy

    assert sys.hexversion >= 0x3060000, "Need Python > 3.6 to guarantee key ordering in dictionary!"

    if block.attr_synos is not None:
        return deepcopy(block.attr_synos)  # Important: should deep copy because we pop from this dict later on
    # Get parent class
    base_syn = {}
    if block.name_base in all_blocks:
        base_syn = generate_attribute_synos(all_blocks[block.name_base], all_blocks)

    # Build final return value
    # 1. First the synos of the keyword itself:
    ret = {None: [s for s in block.synos if s != block.name]}
    # 2. Then the attributes of the current class
    for attr in block.attrs:
        if attr.type != "suppress_param":
            attr_nam = valid_variable_name(attr.name)
            ret[attr_nam] = attr.synos[:]
            # Append original attribute name if it was modified by valid_variable_name():
            if attr_nam != attr.name:
                ret[attr_nam].append(attr.name)
        else:
            # Remove 'suppressed' attributes from inherited part:
            if attr.name not in base_syn:
                _, _path, _line = _info_to_triple(block.info)
                raise Exception(
                    pretty_error(
                        _path,
                        _line,
                        f"invalid 'suppress_param' directive: attribute {attr.name!r} is not present in "
                        f"the parent class '{block.name_base}' and so cannot be suppressed.",
                    )
                )
            else:
                base_syn.pop(attr.name)

    # 3. Then the attributes of the inherited class (after suppressed have been removed, and discarding None key
    # which corresponds to the syno of the keyword itself!)
    for k in base_syn:
        if k is not None:
            ret[k] = base_syn[k]
    block.attr_synos = ret

    return ret


def get_list_type_from_block(block, all_blocks):
    """From a TRAD2 block get the relevant type that should be used in the Pydantic class for lists.
    This also handles nested lists.
    Assumption: terminal type is always a complex type (=another TRUST keyword) in case of nested lists.
    This is checked with an assert.
    """
    assert isinstance(block, tu.TRAD2BlockList)
    it_typ = block.itemtype
    assert it_typ in all_blocks, "List item type is not another TRUST class! Not implemented."
    it_blk = all_blocks[it_typ]
    if isinstance(it_blk, tu.TRAD2BlockList):
        s_typ, most_nested_type = get_list_type_from_block(it_blk, all_blocks)  # Recurse
    else:
        s_typ = ToPydName(it_blk.name)
        most_nested_type = it_blk
    return f'Annotated[List[{s_typ}], "{ToPydName(block.name)}"]', most_nested_type


def write_pyd_block(block, pyd_file, all_blocks, recur_level=0):
    """Write a TRAD2Block as a pydantic class, in the pyd_file"""
    assert isinstance(block, tu.TRAD2Block)
    if block.pyd_written:
        return

    # detect and nicely diagnose infinite recursion
    _, _blk_path, _blk_line = _info_to_triple(block.info)
    if recur_level > 80:
        print(
            pretty_error(
                _blk_path,
                _blk_line,
                f" Infinite recursion! Something wrong with block named '{block.name}'!! Is an attribute of the class, the class itself?",
            )
        )
    if recur_level > 100:
        raise Exception(
            pretty_error(
                _blk_path,
                _blk_line,
                f" Infinite recursion! Something wrong with block named '{block.name}'!! Is an attribute of the class, the class itself?",
            )
        )

    # check base class actually exists!
    if (
        block.name != "objet_u"
        and block.name_base not in ["objet_u"]
        and not block.name.startswith("listobj")
        and block.name_base not in all_blocks
    ):
        raise Exception(
            pretty_error(
                _blk_path,
                _blk_line,
                f" keyword '{block.name}' is declared having parent '{block.name_base}', but this parent does not exist!",
            )
        )

    # dependencies must be written before self (see down below for list items)
    dependencies = [block.name_base] + [a.type for a in block.attrs]
    for dpy in dependencies:
        dpy_block = all_blocks.get(dpy, None)
        if dpy_block:
            write_pyd_block(dpy_block, pyd_file, all_blocks, recur_level=recur_level + 1)

    # Get base class name. If void (like for Objet_U), inherit from TRUSTBaseModel:
    base_cls_n = ToPydName(block.name_base) or "TRUSTBaseModel"
    lines = [
        "#" * 64,
        "",
        f"class {ToPydName(block.name)}({base_cls_n}):",
    ]
    lines += format_docstring(block.desc)

    # Prepare keyword and arg synonyms:
    synonyms = generate_attribute_synos(block, all_blocks)

    typ_map = {
        "entier": ("int", "default=0"),
        "floattant": ("float", "default=0.0"),
        "chaine": ("str", "default=''"),
        "rien": ("bool", "default=False"),
        "list": ("List[float]", "default_factory=list"),
        "listf": ("Annotated[List[float], 'size_is_dim']", "default_factory=list"),
        "listentier": ("List[int]", "default_factory=list"),
        "listentierf": ("Annotated[List[int], 'size_is_dim']", "default_factory=list"),
        "listchaine": ("List[str]", "default_factory=list"),
        "listchainef": ("Annotated[List[str], 'size_is_dim']", "default_factory=list"),
    }

    for attr in block.attrs:
        assert isinstance(attr, tu.TRAD2Attr)
        # `default_factory=<class>` not `lambda: eval("<class>()")`. The
        # `eval` form bought late name binding, but it isn't needed:
        # write_pyd_block above already emits every dependency before
        # self (line 188-193), and any cycle is caught by the recursion
        # guard. Audit 4.7.
        args = f"default_factory={ToPydName(attr.type)}"
        attr_typ = attr.type
        attr_desc = substitute_nl(attr.desc, "")
        # Display-only Field metadata appended after `args`; survives the
        # optional-attribute rewrite below (where `args` -> default=None).
        extra_meta = ""

        if attr.type in all_blocks:  # Complex attribute (=another TRUST class)
            cls = all_blocks[attr.type]
            if isinstance(cls, tu.TRAD2BlockList):
                attr_typ, most_nested_typ = get_list_type_from_block(cls, all_blocks)
                # Make sure list item type is written before the list itself:
                write_pyd_block(most_nested_typ, pyd_file, all_blocks)
                attr_desc = cls.desc
                args = "default_factory=list"
            else:
                attr_typ = ToPydName(attr.type)
        elif attr.type in typ_map:
            attr_typ, args = typ_map[attr.type]
        elif attr.type.startswith("chaine(into="):
            choices = attr.type[13:].split("]")[0]
            _check_nonempty_enum(choices, attr)
            attr_typ = f"Literal[{choices}]"
            args = f"default={choices.split(',')[0]}"
            decl_default = _chaine_declared_default(attr.type, choices, attr)
            if decl_default is not None:
                extra_meta = f', json_schema_extra={{"trust_default": {decl_default}}}'
        elif attr.type.startswith("entier(into="):
            choices = attr.type[13:].split("]")[0]
            _check_nonempty_enum(choices, attr)
            choices = [int(x.strip('"').strip("'")) for x in choices.split(",")]
            attr_typ = f"Literal{choices}"
            args = f"default={choices[0]}"
            decl_default = _entier_declared_default(attr.type)
            if decl_default is not None:
                if decl_default not in choices:
                    _, _at_path, _at_line = _info_to_triple(attr.info)
                    raise Exception(
                        pretty_error(
                            _at_path,
                            _at_line,
                            f"attribute {attr.name!r}: default {decl_default} in type "
                            f"{attr.type!r} is not one of the declared values {choices}.",
                        )
                    )
                extra_meta = f', json_schema_extra={{"trust_default": {decl_default}}}'
        elif attr.type.startswith("entier(min=") or attr.type.startswith("entier(max="):
            attr_typ = "int"
            _m_min = _ENTIER_MIN_RE.search(attr.type)
            _m_max = _ENTIER_MAX_RE.search(attr.type)
            min_val = int(_m_min.group(1)) if _m_min else None
            max_val = int(_m_max.group(1)) if _m_max else None
            if min_val is not None and max_val is not None and min_val > max_val:
                _, _at_path, _at_line = _info_to_triple(attr.info)
                raise Exception(
                    pretty_error(
                        _at_path,
                        _at_line,
                        f"attribute {attr.name!r}: min {min_val} exceeds max {max_val} "
                        f"in type {attr.type!r} (empty range).",
                    )
                )
            # Runtime default: 0 clamped into the [min, max] range so it
            # never violates the emitted ge=/le= constraints.
            default_val = 0
            if min_val is not None:
                default_val = max(default_val, min_val)
            if max_val is not None:
                default_val = min(default_val, max_val)
            parts = [f"default={default_val}"]
            if min_val is not None:
                parts.append(f"ge={min_val}")
            if max_val is not None:
                parts.append(f"le={max_val}")
            args = ", ".join(parts)
            decl_default = _entier_declared_default(attr.type)
            if decl_default is not None:
                if (min_val is not None and decl_default < min_val) or (max_val is not None and decl_default > max_val):
                    _, _at_path, _at_line = _info_to_triple(attr.info)
                    raise Exception(
                        pretty_error(
                            _at_path,
                            _at_line,
                            f"attribute {attr.name!r}: default {decl_default} in type "
                            f"{attr.type!r} is outside the allowed range "
                            f"[{min_val if min_val is not None else '-inf'}, "
                            f"{max_val if max_val is not None else '+inf'}].",
                        )
                    )
                extra_meta = f', json_schema_extra={{"trust_default": {decl_default}}}'
        elif attr.type == "suppress_param":
            # NOTE tricky
            # hide the inherited *instance* attribute by a *class* attribute with same name
            attr_typ = "ClassVar[str]"
            attr_desc = "suppress_param"
            args = 'default="suppress_param"'
        elif attr.type.startswith("ref_"):
            attr_typ = "str"
            args = 'default=""'
        else:
            _, _at_path, _at_line = _info_to_triple(attr.info)
            message = f"unresolved type in XD tags : '{attr.type}' - this was found line {_at_line} of file '{_at_path}' - did you forget to define it?"
            logger.error(message)
            raise NotImplementedError(message)

        # Is the attribute optional?
        if attr.is_opt and attr_desc != "suppress_param":
            attr_typ = f"Optional[{attr_typ}]"
            args = "default=None"  # attribute is not set by default
        # Fix attribute name:
        attr_nam = valid_variable_name(attr.name)

        lines.append(f'    {attr_nam}: {attr_typ} = Field(description=r"""{attr_desc}""", {args}{extra_meta})')

    lines += [
        f"    _synonyms: ClassVar[dict] = {synonyms}",
    ]
    lines.append("\n")
    # actual file writing
    pyd_file.write("\n".join(lines))
    # make sure we won't write the same block again:
    block.pyd_written = True


def write_pars_block(block, pars_file, all_blocks):
    """Write a TRAD2Block as a Parser class, in the pars_file"""
    import trustify.core.base as tb

    assert isinstance(block, tu.TRAD2Block)
    if block.pars_written:
        return

    # dependencies must be written before self
    dependencies = [block.name_base] + [a.type for a in block.attrs]
    for dependency in dependencies:
        dependency = all_blocks.get(dependency, None)
        if dependency and dependency.name not in ["listobj", "listobj_impl"]:
            write_pars_block(dependency, pars_file, all_blocks)

    # Get base class name. If void (like for Objet_U), inherit from base.ConstrainBase_Parser:
    if block.name_base != "":
        if block.name_base == "listobj":
            assert isinstance(block, tu.TRAD2BlockList)
            # yes, list of lists might exist!!
            base_cls_n = "base.ListOfBuiltin_Parser" if block.itemtype.startswith("list") else "base.ListOfBase_Parser"
        else:
            base_cls_n = ToParserName(block.name_base)
    else:
        base_cls_n = "base.ConstrainBase_Parser"
    cls_nam = ToParserName(block.name)
    # If the class is already defined in base.py, skip it:
    if cls_nam in tb.__dict__:
        logger.debug(f"Class {cls_nam} already found in 'base' module - not generated.")
        block.pars_written = True
        return
    lines = [
        "#" * 64,
        "",
        f"class {cls_nam}({base_cls_n}):",
    ]

    info_attr = {}
    for attr in block.attrs:
        info_attr[valid_variable_name(attr.name)] = _info_to_triple(attr.info)

    lines += [f'    _braces: str = "{block.mode}"']

    # The XXX_base class (and similar) will need to read their type directly in the dataset, so need this:
    if is_base_or_deriv(block.name):
        lines += ["    _read_type: bool = True"]

    # Lists need extra information:
    if block.name_base == "listobj":
        assert isinstance(block, tu.TRAD2BlockList)
        lines += [f'    _comma: str = "{block.comma}"', f"    _itemType: Objet_u = {ToPydName(block.itemtype)}"]

    lines += [
        f"    _infoMain: tuple = {_info_to_triple(block.info)}",
        f"    _infoAttr: dict = {info_attr}",
        "    _attributeList: Optional[list] = None",  # Needs to be defined in every class.
        "    _attributeSynos: Optional[dict] = None",  # Idem - see methods _GetAtributeList() and _InvertSynos() in base.py
    ]
    lines.append("\n")
    # actual file writing
    pars_file.write("\n".join(lines))
    block.pars_written = True


def generate_pyd_and_pars(trad2_filename, source_locations_filename, out_pyd_filename, out_pars_filename):
    """Generate both modules (pydantic and parsing ones).

    `source_locations_filename` is passed straight to
    `TRAD2Content.BuildContentFromTRAD2(..., source_locations=...)` —
    pass ``None`` to auto-detect the sibling next to `trad2_filename`
    (silent if absent — schemas without bundled source locations are
    officially supported), or a path to load an explicit
    `source_locations.json`.
    """
    all_blocks = tu.TRAD2Content.BuildContentFromTRAD2(trad2_filename, source_locations=source_locations_filename).data

    # add two base classes that are required and not declared in TRAD2 file.
    # `mode` uses the named token rather than the legacy numeric `1` so that
    # if a future refactor lifts the skip-gates in `write_pars_block` (line
    # 313) or the explicit `listobj_impl` filter (line 597-598), the
    # emitted `_braces` literal stays one of the canonical {BRACE,
    # NO_BRACE, INHERITS_BRACE} keys — `WithBraces()` would otherwise miss
    # the lookup at `base.py:301` and raise "Internal error". Audit 1.5.
    objet_u = tu.TRAD2Block()
    objet_u.name = "objet_u"
    objet_u.mode = "BRACE"
    listobj_impl = tu.TRAD2Block()
    listobj_impl.name = "listobj_impl"
    listobj_impl.mode = "NO_BRACE"
    all_blocks += [objet_u, listobj_impl]

    # make a dict to easily find block by name
    all_blocks = {block.name: block for block in all_blocks}

    # add two properties used during writing of blocks
    for block in all_blocks.values():
        block.pyd_written = False
        block.pars_written = False

    header_com = f"""
        ################################################################
        # This file was generated automatically from :
        # {trad2_filename!s}
        ################################################################

    """
    header_pyd = (
        header_com
        + '''
        from typing import Annotated, ClassVar, List, Literal, Optional, Any, Dict
        import pydantic
        from pydantic import ConfigDict, Field, create_model

        class TRUSTBaseModel(pydantic.BaseModel):
            model_config = ConfigDict(validate_assignment=True, protected_namespaces=())

            # Back-reference to the pyd module this class hierarchy belongs to.
            # Populated by the pyd-module footer (one assignment per generated
            # module). Each pyd module has its own TRUSTBaseModel — schema A's
            # subclasses inherit the schema-A value via MRO, schema B's inherit
            # the schema-B value. This lets parser-side code resolve "which
            # schema does this pyd object belong to" without consulting the
            # shared `Abstract_Parser._schema_module` slot (which only tracks
            # the most-recently activated schema).
            _schema_module: ClassVar[Any] = None

            @classmethod
            def with_fields(cls, **field_definitions):
                """ Overriding this classmethod allows to dynamically create new pydantic models
                with fields added dynamically. This is used in hacks.py for FT problems ...
                See https://github.com/pydantic/pydantic/issues/1937
                """
                return create_model(cls.__name__, __base__=cls, **field_definitions)

            def __init__(self, *args, **kwargs):
                pydantic.BaseModel.__init__(self, *args, **kwargs)
                self._parser = None   #: An instance of AbstractCommon_Parser that is used to build the current pydantic object

            def self_validate(self):
                """ Validate an instance to see if it complies with the pydantic schema """
                dmp = self.__class__.model_dump(self, warnings=False, serialize_as_any=True)
                self.__class__.model_validate(dmp)

            def toDatasetTokens(self):
                """ Convert a pydantic object (self) back to a stream of tokens that can be output in a file to reproduce
                a TRUST dataset. """
                if self._parser is None:
                    from trustify.core.naming import _PARSER_SUFFIX
                    _schema_mod = self.__class__._schema_module
                    _parser_cls = getattr(_schema_mod, self.__class__.__name__ + _PARSER_SUFFIX)
                    self._parser = _parser_cls()
                self._parser._pyd_value = self
                return self._parser.toDatasetTokens()

            def __getattribute__(self, nam):
                """ Override to allow the (scripting) user to use attribute synonyms, for
                example 'pb.post_processing' instead of 'pb.postraitement' ...

                Hot path: most accesses are for actual fields, dunders, or pydantic-
                internal names that can never be user-facing synonyms. Skip the synonym
                walk entirely in those cases. Also access `model_fields` via the class
                (type()) rather than via the instance — instance-level access triggers a
                PydanticDeprecatedSince211 warning on every single call, which used to
                burn ~1 second per 60 datasets in the warnings.warn machinery.
                """
                # Dunders / private names are never synonyms.
                if nam.startswith("_"):
                    return super().__getattribute__(nam)
                cls = super().__getattribute__("__class__")
                if nam in cls.model_fields:
                    return super().__getattribute__(nam)
                # Synonym walk — only reached for non-field, non-private names.
                while cls is not TRUSTBaseModel:
                    for attr_nam, lst_syno in cls._synonyms.items():
                        if attr_nam is not None and nam in lst_syno:
                            return super().__getattribute__(attr_nam)
                    cls = cls.__base__
                return super().__getattribute__(nam)

            def __setattr__(self, nam, val):
                """ Override to allow the (scripting) user to use attribute synonyms, for
                example 'pb.post_processing' instead of 'pb.postraitement' ...

                Same fast-path + class-level model_fields tricks as __getattribute__.
                """
                # Dunders / private names are never synonyms.
                if nam.startswith("_"):
                    return super().__setattr__(nam, val)
                cls = self.__class__
                if nam in cls.model_fields:
                    return super().__setattr__(nam, val)
                # Synonym walk — only reached for non-field, non-private names.
                while cls is not TRUSTBaseModel:
                    for attr_nam, lst_syno in cls._synonyms.items():
                        if attr_nam is not None and nam in lst_syno:
                            return super().__setattr__(nam, val)
                    cls = cls.__base__
                return super().__setattr__(nam, val)

    '''
    )

    footer_pyd = '''
    ##################################################
    ### Classes Declaration and Dataset
    ##################################################
    class Declaration(Objet_u):
        """ Added Pydantic class to handle forward declaration in the TRUST dataset """
        ze_type    : type = Field(description='Class type being read in the forward declaration', default=None)
        identifier : str  = Field(description='Name assigned to the object in the dataset', default='??')

    class Read(Interprete):
        """ The 'read' instruction in a TRUST dataset. Overriden from the automatic generation to make the second argument a Objet_u.
            See also Read_Parser class in base.py module.
        """
        identifier: str = Field(description='Identifier of the class being read. Must match a previous forward Declaration.', default="??")
        obj       : Objet_u = Field(description='The object being read.', default=None)
        _synonyms: ClassVar[dict] = {None: ['lire'], 'identifier': [], 'obj': []}

    class Dataset(Objet_u):
        """ A full TRUST dataset! It is an ordered list of objects. """
        _declarations: Dict[str,Any] = { }  # Private - key: declaration name (like 'pb' in 'pb_conduction pb'), value: a couple (cls, index)
                                            # where 'cls' is a Declaration object, and index is its position in the 'entries' member
        entries   : List[Objet_u]   = Field(description='The objects making up the dataset', default=[])

        def get(self, identifier):
            """ User method - Returns the object associated with a name in the data set """
            from trustify.core.misc_utilities import TrustifyException
            if not identifier in self._declarations:
                raise TrustifyException(f"Invalid identifer '{identifier}'")
            it_num = self._declarations[identifier][1]
            if it_num < 0:
                raise TrustifyException(f"Identifer '{identifier}' has been declared, but has not been read in the dataset (no 'read {identifier} ...' instruction)")
            # Return the object attached into the 'read' instance:
            return self.entries[it_num].obj
    '''

    header_pars = (
        header_com
        + f"""
        import trustify.core.base as base
        from trustify.core.base import *
        # Import all the Pydantic generated classes:
        from {out_pyd_filename.stem} import *
    """
    )

    footer_pars = """
        ################################################################
        # Per-schema-module mechanism: set _schema_module back-refs
        # + pre-compute module-level caches consumed by SchemaModule
        # accessors.
        #
        # `_activate_schema()` is re-callable: when multiple schemas
        # coexist in one process (typical in tests), each schema's
        # `_schema_module` class-level attribute on the shared base
        # parser classes (Dataset_Parser etc.) gets reassigned. A
        # caller can re-activate this module's schema before using it
        # to make sure the shared base classes route lookups through
        # THIS module's namespace.
        ################################################################
        import sys as _sys
        from trustify.core.schema_module import (
            _build_synonyms, _build_all_constrain_base_pyd,
            _build_all_constrain_base_parser)

        _this_module = _sys.modules[__name__]

        def _activate_schema():
            for _val in list(globals().values()):
                if isinstance(_val, type) and issubclass(_val, Abstract_Parser):
                    _val._schema_module = _this_module

        _activate_schema()

        # Pyd-side schema affinity: every TRUSTBaseModel subclass in
        # this module's namespace inherits `_schema_module = this
        # parser module` via MRO. Set once at parser-module-import
        # time on the pyd module's own TRUSTBaseModel class object (so
        # two coexisting schemas keep distinct affinity — each pyd
        # module has its own TRUSTBaseModel). Consumed by
        # `Abstract_Parser._parser_cls_for` to route class-resolution
        # through the pyd value's own schema rather than through the
        # shared parser-base `_schema_module` slot.
        TRUSTBaseModel._schema_module = _this_module

        _SYNO_ORIG_NAME = _build_synonyms(_this_module)
        _all_constrain_base_pyd = _build_all_constrain_base_pyd(_this_module)
        _all_constrain_base_parser = _build_all_constrain_base_parser(_this_module)

        ################################################################
        # Apply hacks.
        ################################################################
    """

    out_pyd_filename = out_pyd_filename or trad2_filename.name + "_pyd.py"
    out_pars_filename = out_pars_filename or trad2_filename.name + "_pars.py"

    # Generate output dirs if needed:
    out_pyd_filename.parents[0].mkdir(parents=True, exist_ok=True)
    out_pars_filename.parents[0].mkdir(parents=True, exist_ok=True)

    # Writing Pydantic schema:
    with open(out_pyd_filename, "w", encoding="utf-8") as pyd_file:
        pyd_file.write(textwrap.dedent(header_pyd).strip())
        pyd_file.write("\n" * 3)
        for block in all_blocks.values():
            write_pyd_block(block, pyd_file, all_blocks)
        # Put implementation of Declaration and Dataset:
        pyd_file.write(textwrap.dedent(footer_pyd).strip())
        pyd_file.write("\n")

    # Writing parsing classes:
    with open(out_pars_filename, "w", encoding="utf-8") as pars_file:
        pars_file.write(textwrap.dedent(header_pars).strip())
        pars_file.write("\n" * 3)
        for block in all_blocks.values():
            if block.name in ["listobj", "listobj_impl"]:  # Those were directly redirected to ListOfBase_Parser
                continue
            write_pars_block(block, pars_file, all_blocks)
        # Register all classes in factory:
        pars_file.write(textwrap.dedent(footer_pars).strip())
        pars_file.write("\n")
        # Apply Phase-1-extracted hacks via the importable hook (see trustify/hacks.py).
        # The body of hacks.py is no longer baked into each generated module.
        pars_file.write(
            textwrap.dedent("""
            import sys as _sys
            from trustify.core.hacks import apply as _trustify_apply_hacks
            _trustify_apply_hacks(_sys.modules[__name__])
        """).strip()
            + "\n"
        )

    return all_blocks
