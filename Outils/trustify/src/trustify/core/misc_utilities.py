"""
Various utilities.
"""

import contextlib
import functools
import logging
import sys
import threading
from typing import ClassVar

from trustify._color import BLUE, GREEN, RED, RESET, YELLOW, colors_enabled, paint


################################################################
# Logger
################################################################
class CustomFormatter(logging.Formatter):
    """Logger formatter that colors records by level when the handler's
    stream (stderr, by default) is a TTY and `NO_COLOR` is not set.
    Colorization is resolved per `format()` call so a runtime change
    to `sys.stderr` or the env is honoured immediately — important
    for tests that monkeypatch stderr to a non-TTY buffer.
    """

    _FMT: ClassVar[str] = "%(levelname)s - [%(filename)s:%(lineno)d] -- %(message)s"
    _COLOR_BY_LEVEL: ClassVar[dict[int, str]] = {
        logging.DEBUG: BLUE,
        logging.INFO: GREEN,
        logging.WARNING: YELLOW,
        logging.ERROR: RED,
        logging.CRITICAL: RED,
    }

    def format(self, record):
        color = self._COLOR_BY_LEVEL.get(record.levelno, "")
        log_fmt = f"{color}{self._FMT}{RESET}" if color and colors_enabled(sys.stderr) else self._FMT
        return logging.Formatter(log_fmt).format(record)


def init_logger():
    import os

    logger = logging.getLogger("trustify")
    if os.getenv("TRUSTIFY_DEBUG"):
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(CustomFormatter())
    logger.addHandler(handler)
    return logger


logger = init_logger()


def pretty_error(fname, lineno, msg):
    """Small handy function for error message generation. Colors the
    filename red and the line number blue when stdout is a TTY (and
    `NO_COLOR` unset); plain text otherwise — `paint()` handles the
    gating, so callers can build error strings without worrying about
    where the result will be printed.
    """
    line_marker = paint(f"{lineno + 1}", BLUE) + "\n"
    return f"\n{paint(fname, RED)}:{line_marker}{msg}"


################################################################
# Exception
################################################################
class TrustifyException(Exception):
    """Raised whenever we encounter an error in the parsing process more tokens."""

    def __init__(self, msg="Error!"):
        Exception.__init__(self, msg)


class TrustifyInternalError(TrustifyException):
    """Programmer-bug / schema-configuration error inside trustify itself.

    Distinct from TrustifyParseError (user-facing parse failure with
    structured source location). The CLI surfaces the full traceback
    for these regardless of TRUSTIFY_DEBUG — these are bugs in
    trustify or in the schema generator, and the stack is the only
    actionable signal for diagnosis. Audit 5.4.
    """


class TrustifyParseError(TrustifyException):
    """Schema/grammar error with structured source location.

    Subclasses TrustifyException so existing `except TrustifyException`
    handlers keep working. The first positional arg is the same human-
    readable message ``GenErr`` produces today; structured fields are
    populated alongside it for tools (LSP) that need exact ranges.

    All location fields use 0-based dataset coordinates. ``col`` /
    ``end_col`` may be ``None`` when the parser is past EOF or when no
    token can be pointed at; consumers must degrade gracefully.
    """

    def __init__(
        self,
        message,
        *,
        file_name,
        line,
        col,
        end_line,
        end_col,
        token,
        attr_name,
        kind,
    ):
        TrustifyException.__init__(self, message)
        self.file_name = file_name
        self.line = line
        self.col = col
        self.end_line = end_line
        self.end_col = end_col
        self.token = token
        self.attr_name = attr_name
        self.kind = kind


################################################################
# Module dynamic import
################################################################


_IMPORT_LOCKS_GUARD = threading.Lock()
_IMPORT_LOCKS: dict[str, threading.Lock] = {}


def _import_lock_for(mod_name: str) -> threading.Lock:
    """Return a per-mod_name lock, creating it on first request under the
    global guard. Concurrent imports of the SAME schema serialize on
    this lock; concurrent imports of DIFFERENT schemas proceed in
    parallel. Audit 6.6.
    """
    with _IMPORT_LOCKS_GUARD:
        lock = _IMPORT_LOCKS.get(mod_name)
        if lock is None:
            lock = threading.Lock()
            _IMPORT_LOCKS[mod_name] = lock
        return lock


def import_parser_module(path):
    """Import a generated parser module from ``path``; return a
    :class:`SchemaModule` wrapper.

    The underlying Python module is registered in ``sys.modules`` under
    a stable hash-based name derived from the resolved path. Distinct
    paths get distinct entries even when their basename is the same
    (every generated parser is named ``trustify_gen.py``).

    Thread-safe (audit 6.6): a per-mod_name lock serializes concurrent
    imports of the same schema so two threads can't both create their
    own module instance and race on the `sys.modules` assignment —
    that race used to leak a half-built module (whichever thread lost
    the assignment race still returned the wrapper around its own
    abandoned mod) and could scramble the cross-module `_schema_module`
    back-refs the footer sets up. Different-schema imports remain
    fully concurrent.
    """
    import hashlib
    import importlib.util
    import sys
    from pathlib import Path

    from trustify.core.schema_module import SchemaModule

    abs_path = Path(path).resolve()
    pkg_dir = str(abs_path.parent)

    digest = hashlib.sha256(str(abs_path).encode()).hexdigest()[:12]
    mod_name = f"trustify_schema_{digest}"

    with _import_lock_for(mod_name):
        if mod_name in sys.modules:
            mod = sys.modules[mod_name]
            # Multiple schemas can coexist in one Python process (e.g. the
            # test suite loads both the "simple" hand-written fixture and
            # the live "full" schema). The base parser classes
            # (Dataset_Parser etc., defined in `trustify.core.base`) carry
            # a class-level `_schema_module` back-reference that whichever
            # schema's footer ran last overrides. Re-activate here so this
            # module's classes route lookups through THIS module again.
            activate = getattr(mod, "_activate_schema", None)
            if activate is not None:
                activate()
            return SchemaModule(mod)

        inserted = pkg_dir not in sys.path
        if inserted:
            sys.path.insert(0, pkg_dir)
        try:
            spec = importlib.util.spec_from_file_location(mod_name, str(abs_path))
            if spec is None or spec.loader is None:
                raise ImportError(f"could not load module from {abs_path}")
            mod = importlib.util.module_from_spec(spec)
            sys.modules[mod_name] = mod
            try:
                spec.loader.exec_module(mod)
            except BaseException:
                sys.modules.pop(mod_name, None)
                raise
            return SchemaModule(mod)
        finally:
            if inserted:
                with contextlib.suppress(ValueError):
                    sys.path.remove(pkg_dir)


##################################################################
## Various datasets utils
##################################################################


def simplify_successive_blanks(aStr):
    """
    reduce all succesive blanks and new lines to one blank
    remove leadind/trailing blanks
    simplify_successive_blanks('  aa \n  bb   ') -> 'aa bb'
    """
    res = aStr.split()
    return " ".join(res)


def prune_after_end(data_ex, parser=None):
    """Produce a new exact same string where all the data after the 'end' keyword is discarded.

    `parser` (optional) is a `TRUSTParser` that has already tokenized
    `data_ex`. When set, the second tokenization that this function
    would normally do is skipped — `TRUSTStream` makes shallow copies
    of the parser's token lists, so independence from the caller's
    stream is preserved. `api.check` passes the parser it created
    during the round-trip parse to avoid tokenizing the same dataset
    twice (~2 ms x N files on a batch_check run).
    """
    from trustify.core.trust_parser import TRUSTParser, TRUSTStream

    if parser is None:
        parser = TRUSTParser()
        parser.tokenize(data_ex)
    stream = TRUSTStream(parser)
    j = -1
    for i, t in enumerate(stream.tokLow):
        if t in ["end", "fin"]:
            j = i
            break
    if j > 0:
        stream.dropTail(j)
    else:
        # Many datasets do not have 'end' so just remove last blanks
        for i, t in enumerate(stream.tokLow[::-1]):
            if t != "":
                j = i
                break
        if j > 0:
            stream.dropTail(len(stream) - j - 1)
    return "".join(stream.tok)


class BoolWithMsg:
    def __init__(self, ok, why=""):
        self.ok = ok
        self.why = why


def check_str_equality(s1, s2, print_on_diff=True):
    """
    Check equality between two (long) strings (with line returns) and produce a diff (as Linux 'diff') if not equal.
    @return a CheckBaseXyz object which has 'check.ok' attribute (boolean), and 'check.why' (string) as diff string.
    """
    import difflib as DIFF

    s1, s2 = str(s1), str(s2)
    if s1 != s2:
        why = ""
        str1 = s1.split("\n")
        str2 = s2.split("\n")
        res = DIFF.unified_diff(str1, str2)
        res = list(res)  # unified_diff 'generator' object is not subscriptable
        res = "\n  " + "\n  ".join(res[2:])  # avoid firsts lines
        why = f"Diff between the two strings:\n{res}"
        if print_on_diff:
            print("")
            logger.error(why)
        return BoolWithMsg(False, why)
    return BoolWithMsg(True, "Strings are equal")


def get_single_parser_base(cls):
    """Return parent class in the Parser hierarchy - there should be only one parent"""
    from trustify.core.naming import _PARSER_SUFFIX

    base_cls = [c for c in cls.__bases__ if c.__name__.endswith(_PARSER_SUFFIX)]
    assert len(base_cls) <= 1, "Too many base classes?!"
    return base_cls[0]


##################################################################
## Various unit tests utils
##################################################################


class UnitUtils:
    """
    Useful stuff for unit tests.
    """

    _TRUG: ClassVar[dict] = {}  # The list of generated modules, as Python module objects
    _TRAD2: ClassVar[dict] = {"simple": "TRAD_2_adr_simple"}  # hand-written XD fixture.
    # The "full" slot is generated live from $TRUST_ROOT via
    # trustify.api.generate_schema — no on-disk fixture.

    _NO_REGENERATE = False  # useful for debugging - avoid regenerating each time (this is a bit long with full TRAD2)

    def generate_python_and_import(self, slot):
        """Generate Python files (Pydantic and parsing) for `slot` and
        store the resulting SchemaModule in ``self._TRUG[slot]``.

        - The ``"full"`` slot is generated live from ``$TRUST_ROOT`` via
          ``trustify.api.generate_schema``. This requires ``TRUST_ROOT``
          to be set; the test is skipped otherwise. No stale snapshot
          is kept in the git tree.
        - All other slots use a hand-written ``tests/trad2/<name>``
          fixture and are generated via ``generate_pyd_and_pars`` into
          ``tests/generated/<slot>/``.
        """
        import os
        import pathlib

        if self._TRUG.get(slot) is not None:
            # Re-activate this slot's schema in case another slot was
            # used since (multi-schema processes share class-level
            # `_schema_module` back-refs on the base parser classes).
            activate = getattr(self._TRUG[slot]._module, "_activate_schema", None)
            if activate is not None:
                activate()
            return

        if slot == "full":
            if "TRUST_ROOT" not in os.environ:
                self.skipTest("'full' slot requires $TRUST_ROOT")
            from trustify.api import generate_schema

            cache_dir = generate_schema(trust_root=os.environ["TRUST_ROOT"])
            file_pars = cache_dir / "trustify_gen.py"
        else:
            tstdir = pathlib.Path(self._test_dir)
            file_in = tstdir / "trad2" / self._TRAD2[slot]
            out_dir = tstdir / "generated" / slot
            out_dir.mkdir(parents=True, exist_ok=True)
            from trustify.cache import pyd_module_filename

            file_pyd = out_dir / pyd_module_filename(slot)
            file_pars = out_dir / f"trustify_gen_{slot}.py"
            if not self._NO_REGENERATE or not file_pars.is_file():
                from trustify.core.trad2_pydantic import generate_pyd_and_pars

                # Hand-written fixtures ship without a sibling
                # source_locations.json; `None` (the default) does the
                # right thing — auto-detect, silent if absent — now
                # that no-source-locations is officially supported.
                generate_pyd_and_pars(file_in, None, file_pyd, file_pars)

        self.__class__._TRUG[slot] = import_parser_module(file_pars)
        self.assertIsNotNone(self._TRUG[slot])


#####################################
# Typing related helper methods
#####################################


@functools.cache
def break_type(typ):
    """For a type like 'Optional[Annotated[List[toto], "fixed_size"]]' :
    @return true if Optional, false otherwise
    @return a list of types, e.g. [list, toto] corresponding to the nesting of types provided and always as Python std types (not typing module)
    @return anything extra put as annotation if any (e.g. "fixed_size"), None otherwise

    Cached: callers (Builtin_Parser, _parseAndSetAttribute, ...) pass
    type objects that are themselves cached pydantic field
    declarations — same typ → same result. Hot in batch-check; the
    typing.get_origin / get_args calls underneath are expensive
    relative to a dict lookup. Callers MUST treat the returned list
    as read-only (no mutation), which they currently do.
    """
    from typing import Annotated, Literal, Union, get_args, get_origin

    orig = get_origin(typ)
    # Non-'typing' type
    if orig is None:
        return False, [typ], None
    # Handle 'Optional[toto]' which is actually an alias for 'Union[toto, None]':
    opt = orig is Union
    if opt:
        args = get_args(typ)
        typ = args[0]
    orig = get_origin(typ)

    # Handle any annotation
    ann = None
    if orig is Annotated:
        typ, ann = get_args(typ)
        orig = get_origin(typ)

    # Break down nested type into a list:
    typs = []
    while orig is not None:
        typs.append(orig)
        args = get_args(typ)
        typ = args[0]
        orig = get_origin(typ)
    # Add most nested type (a builtin) or Literal arguments:
    if len(typs) and typs[-1] is Literal:
        typs.extend(args)
    else:
        typs.append(typ)

    return opt, typs, ann


@functools.cache
def strip_optional(typ):
    """Returns the typ but without the 'Optional' part.

    Cached: same hot-path argument-stability rationale as `break_type`.
    """
    from typing import Union, get_args, get_origin

    orig = get_origin(typ)
    # Non-'typing' type
    if orig is None:
        return typ

    if orig is Union:
        return get_args(typ)[0]
    return typ
