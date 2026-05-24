"""Tests for the post-Phase-1 hacks.apply(module) hook. The hack effects
themselves are exercised end-to-end by the rw test suites; here we lock in
the public shape and the no-op-on-missing-names behavior.
"""

import types
import unittest


class TestHacksApply(unittest.TestCase):
    def test_apply_is_a_callable_taking_one_module_argument(self):
        from trustify.core import hacks

        self.assertTrue(callable(hacks.apply))

    def test_apply_is_a_noop_on_an_empty_stub_module(self):
        """Missing names must be tolerated — hacks.apply runs against test
        fixtures (TRAD_2_adr_simple) that don't define Milieu_base_Parser
        etc., and must not raise.
        """
        from trustify.core import hacks

        stub = types.ModuleType("stub_generated")
        # No Milieu_base_Parser / Eqn_base_Parser / Objet_u etc. in stub.__dict__
        hacks.apply(stub)  # must not raise

    def test_apply_tolerates_partial_module_missing_Objet_u(self):
        """`_apply_dynamic_liste_equations` guards on Objet_u being present in
        the module. A stub module that provides Eqn_base / Problem_read_generic
        but has no Objet_u must still complete cleanly (early-return guard fires
        before any constrain-base walk).
        """
        from trustify.core import hacks

        stub = types.ModuleType("stub_partial")
        stub.Eqn_base = object  # sentinel
        stub.Problem_read_generic = object  # sentinel
        hacks.apply(stub)  # must not raise — guard fires before any registry walk

    def test_skip_type_reads_restores_polymorphic_base_parsers(self):
        """`_apply_skip_type_reads` zeroes `_read_type` on every subclass of
        Milieu_base_Parser (so concrete classes like Fluide_Incompressible
        don't re-read their type when they're the only valid value of a
        non-polymorphic attribute). It then restores `_read_type=True` on
        the polymorphic bases that ARE used as polymorphic attribute types
        — `Eqn_base_Parser` (Listeqn items in Scalaires_passifs) and
        `Fluide_base_Parser` (Fluide_Diphasique.fluide0 / fluide1).

        Lock that in: build a synthetic Milieu sweep target, apply the hack,
        then assert the polymorphic bases stay True while a concrete leaf
        goes False.
        """
        from trustify.core import hacks

        class Milieu_base_Parser:
            _read_type = True  # would be False after sweep — sweep target

        class Mor_eqn_Parser:
            _read_type = True

        class Eqn_base_Parser(Mor_eqn_Parser):
            _read_type = True

        class Fluide_base_Parser(Milieu_base_Parser):
            _read_type = True

        class Fluide_Incompressible_Parser(Fluide_base_Parser):
            _read_type = True

        stub = types.ModuleType("stub_milieu_sweep")
        stub.Milieu_base_Parser = Milieu_base_Parser
        stub.Mor_eqn_Parser = Mor_eqn_Parser
        stub.Eqn_base_Parser = Eqn_base_Parser
        stub.Fluide_base_Parser = Fluide_base_Parser
        stub._all_constrain_base_parser = [
            Milieu_base_Parser,
            Mor_eqn_Parser,
            Eqn_base_Parser,
            Fluide_base_Parser,
            Fluide_Incompressible_Parser,
        ]
        hacks._apply_skip_type_reads(stub.__dict__, stub)

        # Polymorphic bases are restored to True so attribute reads on them
        # actually consume a discriminator token.
        self.assertTrue(Eqn_base_Parser._read_type)
        self.assertTrue(Fluide_base_Parser._read_type)
        # Concrete leaf (and the Milieu/Mor_eqn root the sweep targets) are
        # False so they don't double-read the type when they're the resolved
        # outcome of polymorphic dispatch.
        self.assertFalse(Fluide_Incompressible_Parser._read_type)
        self.assertFalse(Milieu_base_Parser._read_type)
        self.assertFalse(Mor_eqn_Parser._read_type)

    def test_problem_generic_handler_installs_brace_slurp_branch(self):
        """`_apply_problem_generic_attribute_handler` overrides
        `Problem_read_generic_Parser.handleUnexpectedAttribute` so that

          (a) tokens recognised in `solved_equations` are dispatched to
              the corresponding equation parser (`liste_equations` hack);
          (b) a literal `{` triggers the FTD_IJK inner-param-block slurp
              — Probleme_FTD_IJK_base::readOn calls
              `param.lire_avec_accolades(is)` after the equations to read
              an extra `{ key val ... }` block; trustify reads it as raw
              tokens and splices them back on write-back.

        Lock the override in (function name) and confirm the `{` branch
        is present in the installed body so a regression that removes
        either branch fails this test.
        """
        import inspect

        from trustify.core import hacks

        class Problem_read_generic_Parser:
            def handleUnexpectedAttribute(self, *_a, **_kw):
                raise AssertionError("default handler was kept — hack didn't install")

            # Real Problem_read_generic_Parser inherits toDatasetTokens
            # from the parser base; the FTD_IJK write-back hack wraps it
            # at class load (audit 6.10).
            def toDatasetTokens(self):
                return []

        class TRUSTTokens:
            pass

        class Builtin_Parser:
            @classmethod
            def InstanciateFromBuiltin(cls, _):
                return None

        stub = types.ModuleType("stub_problem_generic")
        stub.Problem_read_generic_Parser = Problem_read_generic_Parser
        stub.Builtin_Parser = Builtin_Parser
        stub.TRUSTTokens = TRUSTTokens
        stub.Eqn_base = object
        hacks._apply_problem_generic_attribute_handler(stub.__dict__, stub)

        # Override took: not the original method anymore.
        handler = Problem_read_generic_Parser.handleUnexpectedAttribute
        self.assertEqual(handler.__name__, "handleUnexpectedAttribute_pb_generic")
        # The installed body carries the FTD_IJK brace-slurp branch and the
        # equation-alias branch (locks in both root-cause fixes from the
        # TrioCFD trustify_check work).
        src = inspect.getsource(handler)
        self.assertIn('if tok == "{"', src)
        self.assertIn("_ftd_ijk_param_block", src)
        self.assertIn("solved_equations", src)

    def test_problem_generic_brace_slurp_does_not_shadow_instance_toDS(self):
        """Audit 6.10: the FTD_IJK write-back used to patch
        `self.toDatasetTokens = patched_toDS` from inside the handler.
        After the patch fired, every `inst.__dict__` carried a closure
        that captured `self`, `_ftd_ijk_param_block`, etc. — brittle
        across instance reuse (test harness, LSP), and inverted the
        normal class-level method override path. The fix installs the
        write-back logic on the CLASS once and dispatches on per-instance
        state.

        Drive the handler with a fake stream that supplies `{ }` so the
        slurp loop terminates immediately. After the handler returns,
        the instance must NOT have a per-instance `toDatasetTokens`
        attribute — the patch must live on the class.
        """
        from trustify.core import hacks

        class _FakeTok:
            def __init__(self, low=()):
                self._low = list(low)
                self._orig = list(low)

            def low(self):
                return list(self._low)

            def orig(self):
                return list(self._orig)

        class _FakeStream:
            """Token stream emitting `{` (already consumed) then `}` —
            enough for the slurp loop's depth counter to balance.
            """

            def __init__(self):
                self._yet_to_read = ["}"]
                self._last = _FakeTok(["{"])

            def lastReadTokens(self):
                return self._last

            def nextLow(self):
                t = self._yet_to_read.pop(0)
                self._last = _FakeTok([t])
                return t

        class TRUSTTokens:
            def __init__(self, low=None, orig=None):
                self._low = low or []
                self._orig = orig or []

            def low(self):
                return list(self._low)

            def orig(self):
                return list(self._orig)

        class Problem_read_generic_Parser:
            def handleUnexpectedAttribute(self, *_a, **_kw):
                raise AssertionError("default handler kept — hack didn't install")

            def toDatasetTokens(self):
                return ["ORIGINAL"]

        class Builtin_Parser:
            @classmethod
            def InstanciateFromBuiltin(cls, _):
                return None

        stub = types.ModuleType("stub_problem_generic_brace_slurp")
        stub.Problem_read_generic_Parser = Problem_read_generic_Parser
        stub.Builtin_Parser = Builtin_Parser
        stub.TRUSTTokens = TRUSTTokens
        stub.Eqn_base = object
        hacks._apply_problem_generic_attribute_handler(stub.__dict__, stub)

        inst = Problem_read_generic_Parser()
        inst._attrInOrder = []
        # Trigger the `tok == "{"` slurp path.
        Problem_read_generic_Parser.handleUnexpectedAttribute(inst, _FakeStream(), "{", "pb")

        # Instance state for the block was captured (block + position).
        self.assertTrue(hasattr(inst, "_ftd_ijk_param_block"))
        self.assertEqual(inst._ftd_ijk_inner_pos, 0)
        # BUT the toDatasetTokens method must NOT have been shadowed at
        # the instance level — the patch belongs on the class.
        self.assertNotIn(
            "toDatasetTokens",
            inst.__dict__,
            "FTD_IJK write-back patched the INSTANCE — should patch the class",
        )


if __name__ == "__main__":
    unittest.main()
