"""Unit tests for trustify.core.schema_module."""

import types
import unittest


class TestSchemaModuleCore(unittest.TestCase):
    def test_holds_module_reference(self):
        from trustify.core.schema_module import SchemaModule

        mod = types.ModuleType("fake")
        sm = SchemaModule(mod)
        self.assertIs(sm.module, mod)


class TestSchemaModuleLookups(unittest.TestCase):
    def _make_module_with(self, **attrs):
        """Build a synthetic module with the given attributes."""
        mod = types.ModuleType("fake")
        for k, v in attrs.items():
            setattr(mod, k, v)
        return mod

    def test_get_pyd_class(self):
        from trustify.core.schema_module import SchemaModule

        FakeCls = type("Foo", (), {})
        sm = SchemaModule(self._make_module_with(Foo=FakeCls))
        self.assertIs(sm.get_pyd_class("foo"), FakeCls)

    def test_get_parser_class(self):
        from trustify.core.schema_module import SchemaModule

        FakePars = type("Foo_Parser", (), {})
        sm = SchemaModule(self._make_module_with(Foo_Parser=FakePars))
        self.assertIs(sm.get_parser_class("foo"), FakePars)

    def test_get_pyd_class_unknown_raises(self):
        from trustify.core.schema_module import SchemaModule

        sm = SchemaModule(self._make_module_with())
        with self.assertRaises(AttributeError):
            sm.get_pyd_class("nope")

    def test_get_pyd_from_parser(self):
        from trustify.core.schema_module import SchemaModule

        FakePyd = type("Foo", (), {})
        FakePars = type("Foo_Parser", (), {})
        sm = SchemaModule(self._make_module_with(Foo=FakePyd, Foo_Parser=FakePars))
        self.assertIs(sm.get_pyd_from_parser(FakePars), FakePyd)

    def test_get_parser_from_pyd(self):
        from trustify.core.schema_module import SchemaModule

        FakePyd = type("Foo", (), {})
        FakePars = type("Foo_Parser", (), {})
        sm = SchemaModule(self._make_module_with(Foo=FakePyd, Foo_Parser=FakePars))
        self.assertIs(sm.get_parser_from_pyd(FakePyd), FakePars)

    def test_exists_true_for_pyd(self):
        from trustify.core.schema_module import SchemaModule

        sm = SchemaModule(self._make_module_with(Foo=object))
        self.assertTrue(sm.exists("foo"))

    def test_exists_true_for_parser(self):
        from trustify.core.schema_module import SchemaModule

        sm = SchemaModule(self._make_module_with(Foo_Parser=object))
        self.assertTrue(sm.exists("foo"))

    def test_exists_false_for_unknown(self):
        from trustify.core.schema_module import SchemaModule

        sm = SchemaModule(self._make_module_with())
        self.assertFalse(sm.exists("nope"))

    def test_all_constrain_base_pyd_passthrough(self):
        from trustify.core.schema_module import SchemaModule

        marker = [object()]
        sm = SchemaModule(self._make_module_with(_all_constrain_base_pyd=marker))
        self.assertIs(sm.all_constrain_base_pyd(), marker)

    def test_synonyms_for_returns_list_or_empty(self):
        from trustify.core.schema_module import SchemaModule

        FakeCls = type("Foo", (), {})
        sm = SchemaModule(self._make_module_with(_SYNO_ORIG_NAME={"alias": [FakeCls]}))
        self.assertEqual(sm.synonyms_for("alias"), [FakeCls])
        self.assertEqual(sm.synonyms_for("missing"), [])


class TestImportParserModuleReturnsSchemaModule(unittest.TestCase):
    """`import_parser_module` should now return a SchemaModule (not the
    raw module), and use a hash-based sys.modules name so distinct paths
    never collide even when their basename is the same."""

    def test_returns_schema_module(self):
        import tempfile
        from pathlib import Path

        from trustify.core.misc_utilities import import_parser_module
        from trustify.core.schema_module import SchemaModule

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "trustify_gen.py"
            p.write_text("# empty parser stub\nObjet_u = type('Objet_u', (), {})\n")
            sm = import_parser_module(p)
            self.assertIsInstance(sm, SchemaModule)

    def test_two_distinct_paths_get_distinct_sys_modules_entries(self):
        import sys
        import tempfile
        from pathlib import Path

        from trustify.core.misc_utilities import import_parser_module

        with tempfile.TemporaryDirectory() as d:
            p1 = Path(d) / "a" / "trustify_gen.py"
            p2 = Path(d) / "b" / "trustify_gen.py"
            p1.parent.mkdir()
            p2.parent.mkdir()
            p1.write_text("MARKER = 'a'\nObjet_u = type('Objet_u', (), {})\n")
            p2.write_text("MARKER = 'b'\nObjet_u = type('Objet_u', (), {})\n")
            sm1 = import_parser_module(p1)
            sm2 = import_parser_module(p2)
            self.assertEqual(sm1.module.MARKER, "a")
            self.assertEqual(sm2.module.MARKER, "b")
            # Sys.modules names must differ.
            self.assertNotEqual(sm1.module.__name__, sm2.module.__name__)
            # Both still in sys.modules.
            self.assertIn(sm1.module.__name__, sys.modules)
            self.assertIn(sm2.module.__name__, sys.modules)

    def test_reimport_same_path_returns_cached_module(self):
        import tempfile
        from pathlib import Path

        from trustify.core.misc_utilities import import_parser_module

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "trustify_gen.py"
            p.write_text("MARKER = 'x'\nObjet_u = type('Objet_u', (), {})\n")
            sm1 = import_parser_module(p)
            sm2 = import_parser_module(p)
            self.assertIs(sm1.module, sm2.module)

    def test_concurrent_imports_of_same_path_share_module(self):
        """Audit 6.6: without serialization, two threads both pass the
        `if mod_name in sys.modules` check, both call `module_from_spec`
        to make THEIR OWN mod, both race on the `sys.modules[mod_name] =
        mod` assignment, and both `exec_module` their separate copy.
        Whichever loses the sys.modules race still holds a reference
        to a half-built leak. The footer's `_this_module =
        _sys.modules[__name__]` makes things worse: it may grab the
        OTHER thread's mod mid-init, scrambling the cross-module back-
        refs.

        The fix serializes per-module-name with a lock. To make the
        race observable here we widen the window by sleeping inside
        `exec_module` via a wrapper around `spec_from_file_location` —
        without the injection, the GIL serializes Python steps tightly
        enough that the test would usually pass even pre-fix.
        """
        import importlib.util
        import sys
        import tempfile
        import threading
        import time
        from concurrent.futures import ThreadPoolExecutor
        from pathlib import Path
        from unittest.mock import patch

        from trustify.core.misc_utilities import import_parser_module

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "trustify_gen.py"
            p.write_text("MARKER = 'concurrent'\nObjet_u = type('Objet_u', (), {})\n")

            # The race-widening shim: sleep INSIDE
            # `spec_from_file_location` the first time it's called.
            # That call happens between the `sys.modules` cache check
            # and the `sys.modules[mod_name] = mod` assignment — exactly
            # the window where a second thread can slip past the cache
            # check, create its own mod, and overwrite/be-overwritten
            # by the first thread's assignment. Without serialization,
            # the two threads end up holding distinct module objects.
            real_spec_from_file_location = importlib.util.spec_from_file_location
            called = threading.Event()

            def shimmed(*args, **kwargs):
                if not called.is_set():
                    called.set()
                    time.sleep(0.05)
                return real_spec_from_file_location(*args, **kwargs)

            # Make sure no prior import is cached for this path.
            import hashlib

            mod_name = "trustify_schema_" + hashlib.sha256(str(p.resolve()).encode()).hexdigest()[:12]
            sys.modules.pop(mod_name, None)

            with (
                patch.object(importlib.util, "spec_from_file_location", shimmed),
                ThreadPoolExecutor(max_workers=2) as ex,
            ):
                futures = [ex.submit(import_parser_module, p) for _ in range(2)]
                results = [f.result() for f in futures]

            self.assertIs(
                results[0].module,
                results[1].module,
                "Two concurrent imports of the same path returned distinct module objects — "
                "one of them references a leaked half-built module that isn't in sys.modules.",
            )
            # Sys.modules holds exactly that shared instance.
            self.assertIs(sys.modules.get(mod_name), results[0].module)


if __name__ == "__main__":
    unittest.main()
