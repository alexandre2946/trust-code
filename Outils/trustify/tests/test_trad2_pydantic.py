"""Focused unit tests for trustify.core.trad2_pydantic (code generator)."""

import tempfile
import unittest
from pathlib import Path


def _generate_from_trad2_text(text: str):
    """Write `text` as a TRAD2 file in a temp dir and run the code
    generator on it. Returns nothing — raises whatever
    `generate_pyd_and_pars` raises. The generated files are discarded
    with the temp dir.
    """
    from trustify.core.trad2_pydantic import generate_pyd_and_pars

    with tempfile.TemporaryDirectory() as d:
        trad2_path = Path(d) / "TRAD2_trustify"
        trad2_path.write_text(text, encoding="utf-8")
        pyd_path = Path(d) / "trustify_gen_pyd_test.py"
        pars_path = Path(d) / "trustify_gen.py"
        generate_pyd_and_pars(trad2_path, None, pyd_path, pars_path)


def _generate_pyd_text(text: str) -> str:
    """Like `_generate_from_trad2_text` but returns the generated pyd
    module source so callers can assert on the emitted field shape."""
    from trustify.core.trad2_pydantic import generate_pyd_and_pars

    with tempfile.TemporaryDirectory() as d:
        trad2_path = Path(d) / "TRAD2_trustify"
        trad2_path.write_text(text, encoding="utf-8")
        pyd_path = Path(d) / "trustify_gen_pyd_test.py"
        pars_path = Path(d) / "trustify_gen.py"
        generate_pyd_and_pars(trad2_path, None, pyd_path, pars_path)
        return pyd_path.read_text(encoding="utf-8")


class TestChaineIntoDefault(unittest.TestCase):
    """`chaine(into=[...],default="X")` records X as display-only Field
    metadata (``json_schema_extra['trust_default']``); X must be one of
    the listed choices or codegen raises. It does NOT change the runtime
    default — the first choice stays the pydantic default (REQ) / None
    (OPT) — so parse and round-trip behaviour are untouched.
    """

    def test_declared_default_recorded_as_metadata_req(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            '  attr myenum chaine(into=["xxx","xx","yy"],default="yy") myenum REQ An enum with a default.\n'
        )
        pyd = _generate_pyd_text(text)
        self.assertIn('Literal["xxx","xx","yy"]', pyd)
        # Runtime default is unchanged: still the FIRST choice.
        self.assertIn('default="xxx"', pyd)
        # The declared default is recorded as display-only metadata.
        self.assertIn('json_schema_extra={"trust_default": "yy"}', pyd)

    def test_declared_default_survives_optional_override(self):
        # The OPT branch rewrites the default to None for round-trip
        # safety; the metadata must survive that rewrite — this is where
        # a documented default matters most.
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            '  attr myenum chaine(into=["xxx","xx","yy"],default="yy") myenum OPT An optional enum.\n'
        )
        pyd = _generate_pyd_text(text)
        self.assertIn("default=None", pyd)
        self.assertIn('json_schema_extra={"trust_default": "yy"}', pyd)

    def test_no_clause_no_metadata(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            '  attr myenum chaine(into=["xxx","xx","yy"]) myenum REQ A plain enum.\n'
        )
        pyd = _generate_pyd_text(text)
        self.assertIn('default="xxx"', pyd)
        self.assertNotIn("trust_default", pyd)

    def test_default_not_among_choices_raises(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            '  attr myenum chaine(into=["xxx","xx"],default="zz") myenum REQ A bad default.\n'
        )
        with self.assertRaises(Exception) as ctx:
            _generate_pyd_text(text)
        msg = str(ctx.exception)
        # The diagnostic must name the attribute and the offending default.
        self.assertIn("myenum", msg)
        self.assertIn("zz", msg)


class TestEntierDefault(unittest.TestCase):
    """`entier(into=[...],default=N)` and `entier(max=M,default=N)` record N
    as display-only Field metadata (``json_schema_extra['trust_default']``,
    as an int). N must be one of the listed choices / within the bound, or
    codegen raises. The runtime default (first choice / 0) is unchanged.
    """

    def test_into_default_recorded_as_metadata(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            "  attr mynum entier(into=[1,2,3],default=2) mynum REQ An int enum with a default.\n"
        )
        pyd = _generate_pyd_text(text)
        self.assertIn("Literal[1, 2, 3]", pyd)
        # Runtime default unchanged: first choice.
        self.assertIn("default=1", pyd)
        # Declared default recorded as an int metadata value.
        self.assertIn('json_schema_extra={"trust_default": 2}', pyd)

    def test_max_default_recorded_as_metadata(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            "  attr mynum entier(max=10,default=3) mynum REQ A bounded int with a default.\n"
        )
        pyd = _generate_pyd_text(text)
        # Runtime constraint + default unchanged.
        self.assertIn("default=0, le=10", pyd)
        self.assertIn('json_schema_extra={"trust_default": 3}', pyd)

    def test_into_default_not_among_choices_raises(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            "  attr mynum entier(into=[1,2,3],default=9) mynum REQ A bad default.\n"
        )
        with self.assertRaises(Exception) as ctx:
            _generate_pyd_text(text)
        msg = str(ctx.exception)
        self.assertIn("mynum", msg)
        self.assertIn("9", msg)

    def test_max_default_exceeds_bound_raises(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n"
            "  attr mynum entier(max=10,default=42) mynum REQ A default over the max.\n"
        )
        with self.assertRaises(Exception) as ctx:
            _generate_pyd_text(text)
        msg = str(ctx.exception)
        self.assertIn("mynum", msg)
        self.assertIn("42", msg)


class TestEntierMinBound(unittest.TestCase):
    """`entier(min=M[,max=N][,default=K])` emits `ge=M` (and `le=N` when a
    max is given). The runtime default is 0 clamped into the range, so it
    never violates the bounds. A declared `default=K` is recorded as
    metadata and must lie within the range.
    """

    def _gen(self, typ, opt="REQ"):
        text = f"mykw objet_u mykw BRACE A test keyword.\n  attr mynum {typ} mynum {opt} A bounded int.\n"
        return _generate_pyd_text(text)

    def test_min_max_emits_ge_le_and_clamps_default(self):
        # 0 is below the min, so the runtime default clamps up to the min.
        pyd = self._gen("entier(min=2,max=10)")
        self.assertIn("default=2, ge=2, le=10)", pyd)

    def test_min_only_emits_ge_no_le(self):
        pyd = self._gen("entier(min=2)")
        self.assertIn("default=2, ge=2)", pyd)

    def test_zero_in_range_kept_as_default(self):
        pyd = self._gen("entier(min=-5,max=5)")
        self.assertIn("default=0, ge=-5, le=5)", pyd)

    def test_min_exceeds_max_raises(self):
        with self.assertRaises(Exception) as ctx:
            self._gen("entier(min=10,max=2)")
        self.assertIn("mynum", str(ctx.exception))

    def test_declared_default_below_min_raises(self):
        with self.assertRaises(Exception) as ctx:
            self._gen("entier(min=2,max=10,default=1)")
        self.assertIn("mynum", str(ctx.exception))

    def test_declared_default_within_range_is_metadata(self):
        pyd = self._gen("entier(min=2,max=10,default=5)")
        # Runtime default is the clamped 0 -> 2; the declared 5 is metadata.
        self.assertIn("default=2, ge=2, le=10", pyd)
        self.assertIn('json_schema_extra={"trust_default": 5}', pyd)


class TestEmptyEnumDiagnostic(unittest.TestCase):
    """An empty `into=[]` enum list passes the type-string validator (the
    `dico` shorthand legitimately starts as `chaine(into=[])` before
    `XD_ADD_DICO` appends values). If it survives to code generation
    still empty, the generator must raise a CLEAR diagnostic naming the
    attribute — not emit invalid Python (`Literal[]` / `default=`) nor
    crash with a cryptic `int('')` ValueError.
    """

    def test_empty_chaine_into_raises_clear_error(self):
        text = (
            "mykw objet_u mykw BRACE A test keyword.\n  attr myenum chaine(into=[]) myenum OPT An empty enumeration.\n"
        )
        with self.assertRaises(Exception) as ctx:
            _generate_from_trad2_text(text)
        msg = str(ctx.exception)
        self.assertIn("myenum", msg)
        self.assertIn("empty", msg.lower())

    def test_empty_entier_into_raises_clear_error(self):
        text = "mykw objet_u mykw BRACE A test keyword.\n  attr mynum entier(into=[]) mynum OPT An empty enumeration.\n"
        with self.assertRaises(Exception) as ctx:
            _generate_from_trad2_text(text)
        msg = str(ctx.exception)
        self.assertIn("mynum", msg)
        self.assertIn("empty", msg.lower())


if __name__ == "__main__":
    unittest.main()
