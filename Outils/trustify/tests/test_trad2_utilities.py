"""Focused unit tests for trustify.trad2_utilities."""

import os
import subprocess
import sys
import unittest


class TestModuleImportable(unittest.TestCase):
    # Resolved at import-time relative to this test file; independent of any
    # trustify installation that may be ahead on sys.path.
    _REPO_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))

    def test_module_imports_without_TRUST_ROOT(self):
        """trad2_utilities must be importable when TRUST_ROOT is unset — the
        package needs to work as a normal Python module, not only inside a
        TRUST environment.
        """
        env = {k: v for k, v in os.environ.items() if k != "TRUST_ROOT"}
        env["PYTHONPATH"] = self._REPO_SRC
        result = subprocess.run(
            [sys.executable, "-c", "import trustify.core.trad2_utilities"],
            env=env,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            result.returncode,
            0,
            f"import failed:\nstdout={result.stdout}\nstderr={result.stderr}",
        )


class TestTRAD2ContentProjects(unittest.TestCase):
    """Phase B of the source-locations refactor: TRAD2Content carries a
    projects dict (name -> absolute path)."""

    def test_default_projects_is_empty_dict(self):
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content()
        self.assertEqual(c.projects, {})

    def test_constructor_stores_provided_projects(self):
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content(projects={"trust": "/abs/trust", "baltik_a": "/abs/b"})
        self.assertEqual(c.projects, {"trust": "/abs/trust", "baltik_a": "/abs/b"})

    def test_constructor_copies_projects_dict(self):
        from trustify.core.trad2_utilities import TRAD2Content

        original = {"trust": "/abs/trust"}
        c = TRAD2Content(projects=original)
        original["trust"] = "/changed"
        # The content's dict must not have been aliased.
        self.assertEqual(c.projects, {"trust": "/abs/trust"})

    def test_build_from_org_and_sources_populates_projects(self):
        # Real schema build needs $TRUST_ROOT - this test only checks the
        # API surface (projects param flows through).
        if "TRUST_ROOT" not in os.environ:
            self.skipTest("requires TRUST_ROOT")
        from trustify.core.source_location import TRUST_PROJECT_NAME
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content.BuildFromOrgAndSources(
            None,
            [os.path.join(os.environ["TRUST_ROOT"], "src")],
            projects={TRUST_PROJECT_NAME: os.environ["TRUST_ROOT"]},
        )
        self.assertEqual(c.projects[TRUST_PROJECT_NAME], os.environ["TRUST_ROOT"])

    def test_legacy_trust_root_only_synthesizes_projects_dict(self):
        # Some callers (e.g. trustify.api.generate_schema during the
        # transition) still pass trust_root= without projects=. The
        # factory must synthesize {"trust": trust_root} so the
        # conversion pass runs and serialize_source_locations works
        # downstream.
        import os
        import tempfile

        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "src")
            os.makedirs(src)
            # Empty src dir is fine — we just want to confirm the
            # projects dict was synthesized, no real scan needed.
            from trustify.core.source_location import TRUST_PROJECT_NAME
            from trustify.core.trad2_utilities import TRAD2Content

            c = TRAD2Content.BuildFromOrgAndSources(None, [src], trust_root=d)
            self.assertEqual(c.projects, {TRUST_PROJECT_NAME: d})


class TestSourceLocationConversion(unittest.TestCase):
    def test_convert_infos_to_source_locations_replaces_lists(self):
        from trustify.core.source_location import SourceLocation
        from trustify.core.trad2_utilities import TRAD2Attr, TRAD2Block, TRAD2Content

        c = TRAD2Content(projects={"trust": "/home/u/trust-code"})
        blk = TRAD2Block()
        blk.name = "pb_test"
        blk.info = ["/home/u/trust-code/src/Cas/Pb.cpp", 27]
        attr = TRAD2Attr()
        attr.name = "tinit"
        attr.info = ["/home/u/trust-code/src/Sch/Scheme.cpp", 15]
        blk.attrs = [attr]
        c.data = [blk]

        c._convert_infos_to_source_locations()

        self.assertEqual(
            blk.info,
            SourceLocation(project="trust", path="src/Cas/Pb.cpp", line=27),
        )
        self.assertEqual(
            attr.info,
            SourceLocation(project="trust", path="src/Sch/Scheme.cpp", line=15),
        )

    def test_convert_skips_objects_already_converted(self):
        from trustify.core.source_location import SourceLocation
        from trustify.core.trad2_utilities import TRAD2Block, TRAD2Content

        loc = SourceLocation(project="trust", path="x", line=1)
        c = TRAD2Content(projects={"trust": "/r"})
        blk = TRAD2Block()
        blk.info = loc
        c.data = [blk]
        c._convert_infos_to_source_locations()  # must not crash
        self.assertIs(blk.info, loc)

    def test_convert_raises_when_path_outside_projects(self):
        from trustify.core.trad2_utilities import TRAD2Block, TRAD2Content

        c = TRAD2Content(projects={"trust": "/r"})
        blk = TRAD2Block()
        blk.info = ["/somewhere/else/X.cpp", 1]
        c.data = [blk]
        with self.assertRaises(ValueError):
            c._convert_infos_to_source_locations()

    def test_convert_expands_legacy_trust_root_template(self):
        from trustify.core.source_location import SourceLocation
        from trustify.core.trad2_utilities import TRAD2Block, TRAD2Content

        c = TRAD2Content(projects={"trust": "/home/u/trust-code"})
        blk = TRAD2Block()
        blk.info = ["${TRUST_ROOT}/src/Foo.cpp", 10]
        c.data = [blk]
        c._convert_infos_to_source_locations()
        self.assertEqual(
            blk.info,
            SourceLocation(project="trust", path="src/Foo.cpp", line=10),
        )

    def test_convert_raises_on_unexpandable_trust_root_template(self):
        from trustify.core.trad2_utilities import TRAD2Block, TRAD2Content

        # No "trust" project configured, but a ${TRUST_ROOT} template
        # arrives anyway. Must raise with a clear message.
        c = TRAD2Content(projects={"baltik_a": "/abs/b"})
        blk = TRAD2Block()
        blk.info = ["${TRUST_ROOT}/src/X.cpp", 1]
        c.data = [blk]
        with self.assertRaises(ValueError) as cm:
            c._convert_infos_to_source_locations()
        self.assertIn("TRUST_ROOT", str(cm.exception))
        self.assertIn("trust", str(cm.exception).lower())


class TestSourceLocationsJSON(unittest.TestCase):
    def _make_content(self):
        from trustify.core.source_location import SourceLocation
        from trustify.core.trad2_utilities import TRAD2Attr, TRAD2Block, TRAD2Content

        c = TRAD2Content(projects={"trust": "/r"})
        blk = TRAD2Block()
        blk.name = "pb_test"
        blk.info = SourceLocation(project="trust", path="src/Cas/Pb.cpp", line=27)
        a1 = TRAD2Attr()
        a1.name = "tinit"
        a1.info = SourceLocation(project="trust", path="src/Sch/Scheme.cpp", line=15)
        a2 = TRAD2Attr()
        a2.name = "tmax"
        a2.info = SourceLocation(project="trust", path="src/Sch/Scheme.cpp", line=18)
        blk.attrs = [a1, a2]
        c.data = [blk]
        return c

    def test_serialize_source_locations_returns_canonical_json(self):
        c = self._make_content()
        s = c.serialize_source_locations()
        import json

        data = json.loads(s)
        self.assertEqual(data["version"], 1)
        self.assertIn("pb_test", data["entries"])
        e = data["entries"]["pb_test"]
        self.assertEqual(e["project"], "trust")
        self.assertEqual(e["path"], "src/Cas/Pb.cpp")
        self.assertEqual(e["line"], 27)
        self.assertEqual(set(e["attrs"]), {"tinit", "tmax"})
        self.assertEqual(e["attrs"]["tinit"]["line"], 15)

    def test_serialize_is_stable_byte_for_byte(self):
        c1 = self._make_content()
        c2 = self._make_content()
        self.assertEqual(c1.serialize_source_locations(), c2.serialize_source_locations())

    def test_serialize_uses_sorted_keys(self):
        # `pb_test`'s attrs are inserted in [tinit, tmax] order; the
        # serialized form must be alphabetically sorted so the byte
        # output is stable regardless of insertion order.
        c = self._make_content()
        s = c.serialize_source_locations()
        # `tinit` precedes `tmax` alphabetically, and sort_keys=True
        # guarantees this even if attrs were inserted in reverse.
        self.assertLess(s.find('"tinit"'), s.find('"tmax"'))

    def test_toTRAD2_writes_json_file_instead_of_nfo(self):
        import os
        import tempfile

        c = self._make_content()
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "TRAD2_trustify")
            c.toTRAD2(out)
            self.assertTrue(os.path.exists(os.path.join(d, "source_locations.json")))
            self.assertFalse(os.path.exists(out + ".nfo"))


class TestSourceLocationsJSONReader(unittest.TestCase):
    def test_build_content_from_trad2_reads_json_sibling(self):
        import os
        import tempfile

        from trustify.core.source_location import SourceLocation
        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            trad2_path = os.path.join(d, "TRAD2_trustify")
            with open(trad2_path, "w") as f:
                f.write("pb_test objet_lecture pb_test BRACE A test.\n  attr tinit floattant t_init OPT initial time\n")
            with open(os.path.join(d, "source_locations.json"), "w") as f:
                f.write(
                    '{"version":1,"entries":{"pb_test":'
                    '{"attrs":{"tinit":{"line":15,"path":"src/Sch.cpp","project":"trust"}},'
                    '"line":27,"path":"src/Cas/Pb.cpp","project":"trust"}}}'
                )
            c = TRAD2Content.BuildContentFromTRAD2(trad2_path)
        self.assertEqual(len(c.data), 1)
        blk = c.data[0]
        self.assertEqual(blk.info, SourceLocation(project="trust", path="src/Cas/Pb.cpp", line=27))
        self.assertEqual(blk.attrs[0].info, SourceLocation(project="trust", path="src/Sch.cpp", line=15))

    def test_missing_json_warns_but_loads(self):
        # When only the TRAD2 text is available (e.g. a hand-written
        # fixture in the test suite), BuildContentFromTRAD2 should still
        # produce content — but obj.info will be the parser's own
        # `[abs_or_template, lineno]` form, not SourceLocation. The
        # generator does not use `info` in that path, so this is fine.
        import os
        import tempfile

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            trad2_path = os.path.join(d, "TRAD2_trustify")
            with open(trad2_path, "w") as f:
                f.write("pb_test objet_lecture pb_test BRACE A test.\n")
            c = TRAD2Content.BuildContentFromTRAD2(trad2_path)
        self.assertEqual(len(c.data), 1)


_HASH_DETERMINISM_SCRIPT = """
import os, sys
# Silence the scanner's "really long XD line" warnings so the only
# stdout line is the hash. Anything else (warnings, info) goes to
# stderr — see capture in the test below.
import logging
logging.disable(logging.CRITICAL)
class _Sink:
    def write(self_, *args, **kwargs): pass
    def flush(self_): pass
_real_stdout = sys.stdout
sys.stdout = _Sink()
from trustify.cache import compute_content_hash
from trustify.core.source_location import build_projects_dict
from trustify.core.trad2_utilities import TRAD2Content
trust_root = os.environ["TRUST_ROOT"]
proj_dict = build_projects_dict(trust_root=trust_root, projects=[])
content = TRAD2Content.BuildFromOrgAndSources(
    None, [os.path.join(trust_root, "src")],
    projects=proj_dict, trust_root=trust_root,
)
sys.stdout = _real_stdout
print(compute_content_hash(content, "v-test"))
"""


class TestSchemaGenerationIsDeterministic(unittest.TestCase):
    """Schema generation MUST produce byte-identical canonical text on
    every run against the same source — otherwise the content-hash
    cache misses every time and the cache fills up with orphaned
    "duplicate" entries.

    Past regression: `assemble()` built synonym lists from a Python
    `set()`, whose iteration order depends on hash randomization. The
    fix is `sorted(set(...))`. This test exists to keep that property
    locked: any future change that introduces nondeterministic
    iteration (set / unordered dict / glob without sort / etc.) into
    a step contributing to the hash will surface here as a loud
    failure.

    To actually exercise hash randomization the test must spawn TWO
    Python processes with different ``PYTHONHASHSEED`` values — set
    iteration is deterministic within a single process, so a same-
    process double-build wouldn't catch this regression.

    Skipped when `$TRUST_ROOT` is not available.
    """

    def test_hash_is_stable_across_python_hash_seeds(self):
        if "TRUST_ROOT" not in os.environ:
            self.skipTest("requires TRUST_ROOT")
        env_base = os.environ.copy()
        # Two different seeds → different `hash()` behavior for str/bytes
        # objects. If anything in the canonical-text pipeline leaks
        # set/dict iteration order, these two runs will diverge.
        env_a = {**env_base, "PYTHONHASHSEED": "1"}
        env_b = {**env_base, "PYTHONHASHSEED": "2"}
        run_a = subprocess.run(
            [sys.executable, "-c", _HASH_DETERMINISM_SCRIPT],
            env=env_a,
            capture_output=True,
            text=True,
        )
        run_b = subprocess.run(
            [sys.executable, "-c", _HASH_DETERMINISM_SCRIPT],
            env=env_b,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            run_a.returncode,
            0,
            f"seed=1 child failed:\n{run_a.stderr}",
        )
        self.assertEqual(
            run_b.returncode,
            0,
            f"seed=2 child failed:\n{run_b.stderr}",
        )
        hash_a = run_a.stdout.strip()
        hash_b = run_b.stdout.strip()
        self.assertEqual(
            hash_a,
            hash_b,
            f"compute_content_hash differs across Python processes:\n"
            f"  PYTHONHASHSEED=1 → {hash_a}\n"
            f"  PYTHONHASHSEED=2 → {hash_b}\n"
            f"Something in the canonical-text pipeline depends on hash\n"
            f"randomization. Likely culprit: a set() or unordered dict\n"
            f"whose iteration order ends up in the TRAD2 text or in\n"
            f"`serialize_source_locations()`. Check `assemble()`'s\n"
            f"synonym fixup first.",
        )


class TestXDTagAnchored(unittest.TestCase):
    """XD tag detection must trigger ONLY when the tag opens its
    comment (`// XD ...`) — never when 'XD' merely appears as a word
    inside English prose. Regression for the case where 'NB: the XD
    parent below is field_base' was being parsed as a malformed
    declaration.
    """

    def test_xd_inside_english_comment_is_not_a_tag(self):
        from trustify.core.trad2_utilities import _line_carries_tag

        # All of these mention XD inside ordinary prose — none should match.
        for line in (
            "// NB: the XD parent below is field_base",
            "// see XD declaration above for details",
            "// some text that ends with XD",
            "  // this is XD-like but really isn't",
        ):
            with self.subTest(line=line):
                self.assertFalse(
                    _line_carries_tag(line, "XD"),
                    f"line should NOT be flagged as an XD declaration: {line!r}",
                )

    def test_real_xd_opener_is_detected(self):
        from trustify.core.trad2_utilities import _line_carries_tag

        cases = [
            ("// XD foo bar BRACE description", "XD"),
            ("  // XD attr nb_comp entier nb_comp REQ description", "XD"),
            ("param.ajouter(...); // XD_ADD_P rien description", "XD_ADD_P"),
            ("// 2XD foo bar BRACE description", "2XD"),
            ("// 3XD_ADD_DICO option_name", "3XD_ADD_DICO"),
            ("// XD_CONT continuation text", "XD_CONT"),
        ]
        for line, tag in cases:
            with self.subTest(line=line, tag=tag):
                self.assertTrue(
                    _line_carries_tag(line, tag),
                    f"line should be flagged as a {tag} declaration: {line!r}",
                )

    def test_no_space_between_slashes_and_tag_does_not_match(self):
        """'//XD ...' (no space) must not be detected as a valid opener —
        the legacy scanner errored on this; the new check just refuses to
        match, leaving the existing downstream error in place.
        """
        from trustify.core.trad2_utilities import _line_carries_tag

        self.assertFalse(_line_carries_tag("//XD foo bar BRACE", "XD"))


class TestModernizeWarnings(unittest.TestCase):
    """The scanner emits a WARNING for every modernizable pattern, each
    pointing at `trustify modernize`."""

    def _scan(self, source: str) -> str:
        """Write `source` to a temp .cpp file, run scanSourceFiles, and
        capture the logger output."""
        import io
        import logging
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "fixture.cpp"
            f.write_text(source)
            buf = io.StringIO()
            handler = logging.StreamHandler(buf)
            logger = logging.getLogger("trustify")
            logger.addHandler(handler)
            try:
                content = TRAD2Content()
                content.scanSourceFiles([d])
            finally:
                logger.removeHandler(handler)
            return buf.getvalue()

    def test_long_line_warning_fires_at_120(self):
        word = "x" * 5
        long_desc = " ".join([word] * 25)
        source = f"// XD a objet_u a BRACE {long_desc}\n"
        self.assertGreater(len(source.rstrip()), 120)
        out = self._scan(source)
        self.assertIn("long xd line", out.lower())
        self.assertIn("trustify modernize", out)

    def test_legacy_brace_flag_warning(self):
        source = "// XD a objet_u a 1 a short description\n"
        out = self._scan(source)
        self.assertIn("legacy", out.lower())
        self.assertIn("brace", out.lower())
        self.assertIn("trustify modernize", out)

    def test_legacy_opt_flag_warning(self):
        source = "// XD a objet_u a BRACE block\n// XD attr nam typ syn 0 description\n"
        out = self._scan(source)
        self.assertIn("legacy", out.lower())
        self.assertIn("opt", out.lower())
        self.assertIn("trustify modernize", out)

    def test_named_forms_emit_no_warning(self):
        source = "// XD a objet_u a BRACE short description\n// XD attr nam typ syn REQ short description\n"
        out = self._scan(source)
        self.assertNotIn("modernize", out.lower())
        self.assertNotIn("legacy", out.lower())
        self.assertNotIn("long xd line", out.lower())

    def test_xd_add_p_does_not_trigger_opt_warning(self):
        # `// XD_ADD_P` lines carry no user-written opt flag — the
        # scanner derives it from the C++ `Param::REQUIRED` marker.
        # The legacy-opt nudge must NOT fire on these lines; users
        # have nothing to modernize on the XD_ADD_P side.
        source = '// XD foo objet_u foo BRACE block opener\nparam.ajouter("k", &x); // XD_ADD_P int description text\n'
        out = self._scan(source)
        self.assertNotIn("legacy numeric opt flag", out)

    def test_xd_add_dico_does_not_trigger_opt_warning(self):
        # Same rationale for XD_ADD_DICO — opt is auto-derived. The
        # XD_ADD_DICO line itself doesn't even produce a fresh attr;
        # it extends the preceding XD_ADD_P's allowed-values list.
        # The preceding XD_ADD_P must be typed `dico` for the chain
        # to parse.
        source = (
            "// XD foo objet_u foo BRACE block opener\n"
            'param.ajouter("k", &x); // XD_ADD_P dico description\n'
            'param.dictionnaire("alpha"); // XD_ADD_DICO alpha\n'
        )
        out = self._scan(source)
        self.assertNotIn("legacy numeric opt flag", out)

    def test_long_unsplittable_xd_attr_does_not_warn(self):
        # XD attr with a huge `chaine(into=[...])` type — structured
        # prefix alone > 120, so modernize cannot fix it. Warning must
        # be suppressed to avoid unactionable noise.
        values = ",".join(f'"v_{i}"' for i in range(20))
        source = f"// XD block objet_u block BRACE opener\n// XD attr nam chaine(into=[{values}]) syn OPT desc\n"
        out = self._scan(source)
        self.assertNotIn("long xd line", out.lower())

    def test_long_unsplittable_block_header_does_not_warn(self):
        # Block header with a very long name -> structured prefix > 120.
        long_name = "a" * 100
        source = f"// XD {long_name} objet_u {long_name} BRACE desc\n"
        out = self._scan(source)
        self.assertNotIn("long xd line", out.lower())

    def _scan_to_content(self, source: str):
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "fixture.cpp"
            f.write_text(source)
            content = TRAD2Content()
            # Must NOT raise — that's the regression the test guards.
            content.scanSourceFiles([d])
        return content

    def test_level_2_add_p_with_2xd_cont_recovers_description(self):
        # Reproduces the post-Rule-4 canonical form for 2XD_ADD_P:
        # opener on the C++ line, description on the 2XD_CONT line.
        # The scanner must fold the 2XD_CONT into the opener so
        # _parseXD_ADD_P sees both type AND description.
        source = (
            "// XD foo objet_u foo BRACE block opener\n"
            'Param& p = param.ajouter_param("S"); // XD_ADD_P sub\n'
            "// XD_CONT subblock\n"
            "// 2XD sub objet_lecture nul INHERITS_BRACE not_set\n"
            'p.ajouter("C", &C_); // 2XD_ADD_P double\n'
            "// 2XD_CONT description text\n"
        )
        content = self._scan_to_content(source)
        subs = [b for b in content.data if b.name == "sub"]
        self.assertEqual(len(subs), 1)
        attrs = [a for a in subs[0].attrs if a.name == "c"]
        self.assertEqual(len(attrs), 1)
        self.assertEqual(attrs[0].desc, "description text")

    def test_level_3_add_p_with_3xd_cont_recovers_description(self):
        # Same as the level-2 case but one level deeper. 3XD_CONT must
        # also fold into 3XD_ADD_P.
        source = (
            "// XD outer objet_u outer BRACE outer opener\n"
            'Param& p2 = param.ajouter_param("S2"); // XD_ADD_P sub2\n'
            "// XD_CONT level-1 subblock\n"
            "// 2XD sub2 objet_lecture nul INHERITS_BRACE level-2 desc\n"
            'Param& p3 = p2.ajouter_param("S3"); // 2XD_ADD_P sub3\n'
            "// 2XD_CONT level-2 subsubblock\n"
            "// 3XD sub3 objet_lecture nul INHERITS_BRACE level-3 desc\n"
            'p3.ajouter("D", &D_); // 3XD_ADD_P double\n'
            "// 3XD_CONT deepest description\n"
        )
        content = self._scan_to_content(source)
        deepest = [b for b in content.data if b.name == "sub3"]
        self.assertEqual(len(deepest), 1)
        attrs = [a for a in deepest[0].attrs if a.name == "d"]
        self.assertEqual(len(attrs), 1)
        self.assertEqual(attrs[0].desc, "deepest description")

    def test_long_splittable_description_still_warns(self):
        # Structured prefix is short; the DESCRIPTION is what overflows.
        # Modernize CAN fix this -> warning must still fire so the
        # user knows to run it.
        long_desc = " ".join(["word"] * 40)
        source = f"// XD a objet_u a BRACE {long_desc}\n"
        out = self._scan(source)
        self.assertIn("long xd line", out.lower())
        self.assertIn("trustify modernize", out)


class TestParseMacroSuffixStripping(unittest.TestCase):
    """Audit 1.3: `_parseMacro` strips a trailing `_64` from synonym
    names (TRUST's 64-bit shim uses `<class>_64` aliases that should
    fold back into the unsuffixed name). The original implementation
    used `s.replace('_64', '')` which substitutes EVERY occurrence —
    a synonym with a mid-string `_64` (e.g. `foo_64bit_64`) silently
    lost both copies, mangling the synonym table.

    Latent in current TRUST sources (no synonym has an inner `_64`
    today), but the sibling stripping of `_32_64` at the same call
    site already uses slicing — the divergence was the smell."""

    def test_trailing_64_stripped(self):
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content()
        cls_nam, syno = c._parseMacro("Implemente_instanciable", 'Implemente_instanciable(MyClass, "foo_64")')
        self.assertEqual(cls_nam, "MyClass")
        self.assertEqual(syno, "foo")

    def test_mid_string_64_preserved(self):
        # Bug repro: `foo_64bit_64` ends with `_64` so the trailing
        # strip fires — but the inner `_64bit_` substring must NOT be
        # touched.
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content()
        cls_nam, syno = c._parseMacro("Implemente_instanciable", 'Implemente_instanciable(MyClass, "foo_64bit_64")')
        self.assertEqual(cls_nam, "MyClass")
        self.assertEqual(syno, "foo_64bit", "_parseMacro must strip only the TRAILING `_64`, not every occurrence")

    def test_no_64_suffix_left_alone(self):
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content()
        _cls_nam, syno = c._parseMacro("Implemente_instanciable", 'Implemente_instanciable(MyClass, "plain_name")')
        self.assertEqual(syno, "plain_name")


class TestTypeStringValidation(unittest.TestCase):
    """Audit 4.9: every `attr.type` value must conform to one of the
    documented forms. The pydantic generator dispatches on the type
    string via `attr.type.startswith("chaine(into=")` etc. and then
    extracts the payload with `[13:].split("]")[0]` — so any garbage
    AFTER the closing `])` of the legal form is silently dropped, and
    the produced schema does not reflect the XD tag in the source.

    Example: an XD attr line containing `chaine(into=["a","b"])toto`
    used to be accepted verbatim and produced the same schema as the
    canonical `chaine(into=["a","b"])`, losing the `toto` suffix
    without any diagnostic. The fix tightens the type-string syntax to
    `[A-Za-z_][A-Za-z0-9_]*` OR one of the parameterised forms with
    nothing trailing the closing `)`.
    """

    def _valid_types(self):
        return [
            # plain identifiers (typ_map members, block names, suppress_param, ref_*)
            "entier",
            "floattant",
            "chaine",
            "rien",
            "list",
            "listf",
            "listentier",
            "listentierf",
            "listchaine",
            "listchainef",
            "bloc_lecture",
            "Conduction",
            "Pb_thermohydraulique_QC",
            "suppress_param",
            "ref_pb_thermo",
            # chaine(into=...) forms — including the empty form produced by
            # `XD_ADD_P dico` (via convertTyp), and the trailing-comma form
            # produced by the XD_ADD_DICO append in trad2_utilities.py:645
            "chaine(into=[])",
            'chaine(into=["a"])',
            'chaine(into=["a","b"])',
            'chaine(into=["a","b","c"])',
            'chaine(into=["a",])',
            'chaine(into=["a","b",])',
            'chaine(into=["<="])',
            'chaine(into=["="])',
            # chaine(into=...) with a trailing `,default="..."` clause. Only
            # the SYNTAX is checked here; that the value is one of the listed
            # choices is validated later, at code-generation time.
            'chaine(into=["a"],default="a")',
            'chaine(into=["a","b"],default="b")',
            'chaine(into=["xxx","xx","xxx","xx"],default="xx")',
            'chaine(into=["a","b",],default="a")',  # trailing comma in list still OK
            # entier(into=...) forms — bare or quoted ints, negative allowed
            "entier(into=[1])",
            "entier(into=[1,2,3])",
            "entier(into=[-1,0,1])",
            'entier(into=["0","1"])',
            # entier(into=...) / entier(max=...) with a trailing `,default=N`
            # clause (bare or quoted int). Membership / bound is validated
            # later, at code-generation time.
            "entier(into=[1,2,3],default=2)",
            "entier(into=[-1,0,1],default=-1)",
            'entier(into=[1,2],default="2")',
            # entier(max=...) — single int, negative allowed
            "entier(max=2)",
            "entier(max=-5)",
            "entier(max=10,default=3)",
            "entier(max=10,default=-3)",
            # entier(min=...) — lower bound, optionally combined with max
            # (min always precedes max) and/or a default clause.
            "entier(min=2)",
            "entier(min=-5)",
            "entier(min=2,max=10)",
            "entier(min=-5,max=5)",
            "entier(min=2,max=10,default=5)",
            "entier(min=2,default=4)",
        ]

    def _invalid_types(self):
        return [
            # garbage trailing the closing `)`
            'chaine(into=["a","b"])toto',
            'chaine(into=["a"])X',
            "entier(into=[1,2])bar",
            "entier(max=2)baz",
            # malformed brackets / punctuation
            'chaine(into=["a","b")',
            'chaine(into=["a","b"]',
            "entier(max=)",
            # malformed / unsupported `default=` clauses
            'chaine(into=["a"],default="a")X',  # garbage after the closing `)`
            'chaine(into=["a"],default=)',  # default value not quoted
            'chaine(into=["a"],default="a"',  # missing closing paren
            'chaine(into=["a"]default="a")',  # missing comma before default
            "entier(into=[1],default=1)X",  # garbage after the closing `)`
            "entier(into=[1],default=a)",  # default not an int
            "entier(max=5,default=)",  # missing default value
            'chaine(into=["a"],default=1)',  # int default on a chaine enum
            "entier(min=)",  # missing min value
            "entier(max=10,min=2)",  # wrong order — min must precede max
            "entier(min=2,max=10)Z",  # garbage after the closing `)`
            # punctuation in identifier
            "foo-bar",
            "foo.bar",
            "foo bar",  # embedded space (whitespace-split would block this earlier, but still illegal)
            # leading digit (not a valid identifier shape)
            "9foo",
            # unbalanced/mismatched quotes
            'chaine(into=["a)',
            # non-ascii
            "étude",
            # empty
            "",
        ]

    def test_valid_type_strings_pass(self):
        from trustify.core.trad2_utilities import _validate_type_string

        for typ in self._valid_types():
            _validate_type_string(typ, "fake.cpp", 1)  # must not raise

    def test_invalid_type_strings_raise(self):
        from trustify.core.trad2_utilities import _validate_type_string

        for typ in self._invalid_types():
            with self.assertRaises(Exception, msg=f"expected {typ!r} to be rejected") as cm:
                # `pretty_error` reports lineno+1 (TRAD2 convention), so
                # passing 41 here makes the formatted message say `:42`.
                _validate_type_string(typ, "fake.cpp", 41)
            # The error message must name the offending type AND the
            # source location so the user can find the XD tag.
            msg = str(cm.exception)
            self.assertIn(typ, msg, f"error message for {typ!r} should quote the value")
            self.assertIn("fake.cpp", msg)
            self.assertIn("42", msg)

    def test_build_from_tab_rejects_trailing_garbage_in_type(self):
        # End-to-end: a real XD attr line with trailing garbage on the
        # type field must fail at parse time, not silently produce a
        # truncated schema.
        from trustify.core.trad2_utilities import TRAD2Attr

        with self.assertRaises(Exception) as cm:
            TRAD2Attr.BuildFromTab(
                ["weird", 'chaine(into=["a","b"])toto', "weird", "OPT", "some desc"],
                "fixture.cpp",
                10,
            )
        msg = str(cm.exception)
        self.assertIn("chaine(into", msg)
        self.assertIn("toto", msg)


class TestCppArgSplitting(unittest.TestCase):
    """Audit 6.1: `_parseMacro` and `_parseXD_ADD_something` used to
    tear C++ macro arguments apart with `line.split(",")` and the
    FIRST `(` / `)` found in the string. Any nested call (`(this)`,
    `make_name("foo")`) or quoted comma silently misparsed — the
    real lambda_ortho example below produced opt='1' (OPTIONAL)
    while the C++ enforced `Param::REQUIRED`. The fix swaps both
    sites to balanced-paren + string-aware splitting helpers.
    """

    def test_split_cpp_args_top_level_commas_only(self):
        from trustify.core.trad2_utilities import _split_cpp_args

        # Bare commas at the top level: ordinary split.
        self.assertEqual(_split_cpp_args("a, b, c"), ["a", "b", "c"])
        # Nested parens: comma inside (...) does NOT split.
        self.assertEqual(_split_cpp_args('"name", (this), Param::REQUIRED'), ['"name"', "(this)", "Param::REQUIRED"])
        # Comma inside a string literal does NOT split.
        self.assertEqual(_split_cpp_args('"a,b", &v'), ['"a,b"', "&v"])
        # Square / curly brackets behave like parens.
        self.assertEqual(_split_cpp_args("foo[1,2], bar{3,4}, baz"), ["foo[1,2]", "bar{3,4}", "baz"])

    def test_extract_call_args_balanced_parens(self):
        from trustify.core.trad2_utilities import _extract_call_args

        # Pick the OUTER paren pair starting at the named call.
        prefix, args, suffix = _extract_call_args('foo.bar("x", (this))', method="bar")
        self.assertEqual(args, '"x", (this)')
        self.assertEqual(prefix, "foo.")
        self.assertEqual(suffix, "")
        # Skip preceding parens (e.g. `if (cond) foo.bar(...)`).
        _, args, _ = _extract_call_args("if (cond) foo.bar(x, y)", method="bar")
        self.assertEqual(args, "x, y")
        # String literal with parens inside is not confused.
        _, args, _ = _extract_call_args('foo.bar("(", "x")', method="bar")
        self.assertEqual(args, '"(", "x"')

    def test_parse_xd_add_p_with_nested_paren_and_required(self):
        # The real bug — lambda_ortho in
        # src/VEF/Sources/PDC/Perte_Charge_Circulaire_VEF_P1NC.cpp:56.
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content()
        line = 'param.ajouter_non_std("lambda_ortho",(this),Param::REQUIRED); // XD_ADD_P chaine'
        tab = line.split()
        nam, typ, opt, _desc = c._parseXD_ADD_something(
            tag="XD_ADD_P", expected_num_args=1, cpp_meth="ajouter", f_nam="x.cpp", lin_n=0, tab=tab
        )
        self.assertEqual(nam, "lambda_ortho")
        self.assertEqual(typ, "chaine")
        self.assertEqual(opt, "0", "Param::REQUIRED was masked by the (this) nested-paren split")

    def test_parse_xd_add_p_common_case_still_works(self):
        # Regression guard for the most common XD_ADD_P shape.
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content()
        line = 'param.ajouter("longueur_maille", &methode_); // XD_ADD_P chaine'
        tab = line.split()
        nam, typ, opt, _desc = c._parseXD_ADD_something(
            tag="XD_ADD_P", expected_num_args=1, cpp_meth="ajouter", f_nam="x.cpp", lin_n=0, tab=tab
        )
        self.assertEqual(nam, "longueur_maille")
        self.assertEqual(typ, "chaine")
        self.assertEqual(opt, "1")

    def test_parse_macro_handles_synonym_with_comma(self):
        # `add_synonym(foo, "a,b")` used to truncate to `a` silently
        # because line.split(",") cut inside the quoted string. The
        # production scanner always lower-cases the line before
        # invoking `_parseMacro` (see scanOneCppFile), so the test
        # mirrors that.
        from trustify.core.trad2_utilities import TRAD2Content

        c = TRAD2Content()
        cls_nam, syno = c._parseMacro("add_synonym", 'add_synonym(foo,"a,b");')
        self.assertEqual(cls_nam, "foo")
        self.assertEqual(syno, "a,b")


class TestMultiLineMacro(unittest.TestCase):
    """A macro call (`Implemente_instanciable` / `Add_synonym`) may be
    split across several lines — clang-format does this for long class
    names. The scanner's per-line pre-screen used to keep only the
    truncated first line, so `_parseMacro` then hit an unbalanced '('
    and raised a bare `ValueError` with no file/line context. The
    grouping step must instead accumulate continuation lines until the
    call's parens balance.
    """

    def _scan_to_content(self, source: str):
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "fixture.cpp"
            f.write_text(source)
            content = TRAD2Content()
            content.scanSourceFiles([d])
        return content

    def test_call_parens_balanced_helper(self):
        from trustify.core.trad2_utilities import _call_parens_balanced

        self.assertTrue(_call_parens_balanced("Implemente_instanciable(A, B, C);"))
        self.assertFalse(_call_parens_balanced("Implemente_instanciable(A,"))
        self.assertFalse(_call_parens_balanced('Add_synonym(Foo, "b'))  # opened, no close
        # A ')' inside a string literal must not count as a real close.
        self.assertFalse(_call_parens_balanced('Add_synonym(Foo, "a)b"'))
        # No paren at all -> not a (complete) call.
        self.assertFalse(_call_parens_balanced("plain text"))

    def test_multiline_implemente_instanciable_parses_synonyms(self):
        # The exact shape from the bug report, with a piped synonym so
        # we can observe the parse result in `content.synos`.
        source = (
            "// XD energie_cinetique_turbulente convection_diffusion_turbulence_multiphase "
            "energie_cinetique_turbulente NO_BRACE A keyword.\n"
            "Implemente_instanciable(Energie_cinetique_turbulente,\n"
            '                        "energie_cinetique_turbulente|kinetic_energy",\n'
            "                        Convection_diffusion_turbulence_multiphase);\n"
        )
        content = self._scan_to_content(source)
        self.assertEqual(content.synos, {"energie_cinetique_turbulente": ["kinetic_energy"]})

    def test_multiline_add_synonym_parses(self):
        source = 'Implemente_instanciable(Foo, "foo", Bar);\nAdd_synonym(Foo,\n            "barbar");\n'
        content = self._scan_to_content(source)
        self.assertEqual(content.synos, {"foo": ["barbar"]})

    def test_singleline_macro_still_works(self):
        # Regression guard: the common one-line form must be unaffected.
        source = 'Implemente_instanciable(Foo, "foo|bar", Base);\n'
        content = self._scan_to_content(source)
        self.assertEqual(content.synos, {"foo": ["bar"]})

    def test_unterminated_macro_gives_clear_error(self):
        # Never-closed paren (corrupt source): no bare ValueError —
        # a pretty_error pointing at the file must surface instead.
        source = 'Implemente_instanciable(Foo,\n        "foo",\n'
        with self.assertRaises(Exception) as ctx:
            self._scan_to_content(source)
        msg = str(ctx.exception)
        self.assertIn("fixture.cpp", msg)
        self.assertIn("Implemente_instanciable", msg)


class TestMultiLineXDAddP(unittest.TestCase):
    """An XD_ADD_P / XD_ADD_DICO tag pulls data from the `.ajouter(...)` /
    `.dictionnaire(...)` C++ call on the same line. trustify does NOT
    support that call being split across several lines (the tag is a
    `//` comment embedded mid-statement; the call may even open on an
    already-dropped prior line). Such cases must fail with one clear,
    actionable message — not the misleading "method is not called"
    wording, and not a cryptic "unbalanced parens".
    """

    HEADER = "// XD foo objet_u foo BRACE A block.\n"

    def _scan_err(self, body: str) -> str:
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            Path(d, "fixture.cpp").write_text(self.HEADER + body)
            with self.assertRaises(Exception) as ctx:
                TRAD2Content().scanSourceFiles([d])
        return str(ctx.exception)

    def test_call_opens_on_tag_line_continues_below(self):
        # Comment on the FIRST line, the .ajouter() call continues below.
        msg = self._scan_err(
            'param.ajouter("P0", &alphaE_, // XD_ADD_P rien Pressure\n                Param::REQUIRED);\n'
        )
        self.assertIn("fixture.cpp", msg)
        self.assertIn("multi-line", msg)
        self.assertIn("single line", msg)

    def test_call_opens_above_tag_on_closing_line(self):
        # The .ajouter() opens on a PRIOR (dropped) line; the tag sits on
        # the closing line. Used to say "method is not called" — wrong.
        msg = self._scan_err('param.ajouter("P0",\n                &alphaE_); // XD_ADD_P rien Pressure\n')
        self.assertIn("multi-line", msg)
        self.assertIn("single line", msg)
        self.assertNotIn("is not called", msg)

    def test_dico_call_split_across_lines(self):
        # Same hazard for XD_ADD_DICO / .dictionnaire().
        body = (
            'param.ajouter("model", &model_); // XD_ADD_P chaine(into=[]) model\n'
            'model_.dictionnaire("keps", // XD_ADD_DICO keps\n'
            "                    0);\n"
        )
        msg = self._scan_err(body)
        self.assertIn("multi-line", msg)
        self.assertIn("single line", msg)

    def test_genuine_misplacement_keeps_original_message(self):
        # A tag on a line with NO call at all (balanced, no .ajouter) is a
        # real misplacement, not a multi-line call — keep the original,
        # specific diagnostic rather than blaming multi-line.
        msg = self._scan_err("int x = 0; // XD_ADD_P rien Pressure\n")
        self.assertIn("is not called", msg)
        self.assertNotIn("multi-line", msg)


class TestXDNameValidation(unittest.TestCase):
    """A TRUST keyword/attribute name must match [A-Za-z0-9_]+. An illegal
    character (e.g. the '*' in 'y*_switch') used to surface only much
    later — as a bare ValueError in the pydantic generator (attribute
    names) or an invalid Python identifier in the generated module
    (keyword names) — with no pointer back at the offending XD tag.
    Validating at extraction time (in BuildFromTab) makes the diagnostic
    carry the real source file and line.
    """

    def _scan_err(self, source: str) -> str:
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            Path(d, "fixture.cpp").write_text(source)
            with self.assertRaises(Exception) as ctx:
                TRAD2Content().scanSourceFiles([d])
        return str(ctx.exception)

    def test_illegal_attr_name_carries_location(self):
        src = "// XD foo objet_u foo BRACE A block.\n// XD attr y*_switch entier y*_switch OPT The switch.\n"
        msg = self._scan_err(src)
        self.assertIn("fixture.cpp", msg)
        self.assertIn("y*_switch", msg)
        self.assertIn("*", msg)

    def test_illegal_keyword_name_carries_location(self):
        src = "// XD y*_switch objet_u y*_switch NO_BRACE A block.\n"
        msg = self._scan_err(src)
        self.assertIn("fixture.cpp", msg)
        self.assertIn("y*_switch", msg)

    def test_illegal_xd_add_p_param_name_carries_location(self):
        # The name comes from the C++ .ajouter() arg, not a // XD attr line.
        src = '// XD foo objet_u foo BRACE A block.\nparam.ajouter("y*_switch", &x_); // XD_ADD_P entier The switch.\n'
        msg = self._scan_err(src)
        self.assertIn("fixture.cpp", msg)
        self.assertIn("y*_switch", msg)

    def test_valid_names_still_pass(self):
        # Regression guard: legal names (digits, underscores) are untouched.
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        src = "// XD foo objet_u foo BRACE A block.\n// XD attr y_switch_2 entier y_switch_2 OPT The switch.\n"
        with tempfile.TemporaryDirectory() as d:
            Path(d, "fixture.cpp").write_text(src)
            content = TRAD2Content()
            content.scanSourceFiles([d])  # must NOT raise
        blk = next(b for b in content.data if b.name == "foo")
        self.assertEqual([a.name for a in blk.attrs], ["y_switch_2"])


class TestScanErrorAccumulation(unittest.TestCase):
    """XD extraction accumulates every error and reports them together at
    the end of the scan, instead of aborting on the first one. Warnings are
    unaffected (still emitted live via the logger)."""

    def _scan_err(self, files: dict) -> str:
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            for nam, txt in files.items():
                Path(d, nam).write_text(txt)
            with self.assertRaises(Exception) as ctx:
                TRAD2Content().scanSourceFiles([d])
        return str(ctx.exception)

    def test_multiple_errors_same_file_all_reported(self):
        # Two independent bad attrs under one valid block: the scan must not
        # stop at the first — both must appear in the single combined error.
        src = (
            "// XD foo objet_u foo BRACE A block.\n"
            "// XD attr bad*nm entier bad*nm OPT first: illegal name.\n"
            "// XD attr okname ba**d okname OPT second: illegal type.\n"
        )
        msg = self._scan_err({"fixture.cpp": src})
        self.assertIn("bad*nm", msg)  # first error
        self.assertIn("ba**d", msg)  # second error survived the first
        self.assertIn("2 error", msg)

    def test_errors_across_files_all_reported(self):
        # One error in each of two files — both surface in one report.
        msg = self._scan_err(
            {
                "a.cpp": "// XD foo objet_u foo BRACE A.\n// XD attr n1 b*d1 n1 OPT x.\n",
                "b.cpp": "// XD bar objet_u bar BRACE B.\n// XD attr n2 b*d2 n2 OPT y.\n",
            }
        )
        self.assertIn("b*d1", msg)
        self.assertIn("b*d2", msg)

    def test_clean_source_does_not_raise(self):
        import tempfile
        from pathlib import Path

        from trustify.core.trad2_utilities import TRAD2Content

        with tempfile.TemporaryDirectory() as d:
            Path(d, "ok.cpp").write_text("// XD foo objet_u foo BRACE A block.\n// XD attr p entier p OPT good.\n")
            TRAD2Content().scanSourceFiles([d])  # must NOT raise


class TestUtf8EncodingUnderCLocale(unittest.TestCase):
    """Audit 4.5: every open() in the schema-extraction path
    (TRAD2 read/write, source_locations.json read/write, C++/.xd
    source read) must pass `encoding="utf-8"`. Otherwise the
    platform-default codec applies — under LANG=C / LC_ALL=POSIX
    that falls back to ASCII, and any non-ASCII byte (a common
    accented character in a C++ description string) raises
    `UnicodeDecodeError`. Schema generation then breaks on
    locale-misconfigured machines, and the same source text hashes
    differently across locales so a cache written on machine A
    cannot be reused on machine B.

    We force the failure scenario from a subprocess: LANG=C, LC_ALL=C,
    and PYTHONUTF8=0 disable Python's UTF-8 mode fallback (PEP 540),
    leaving the interpreter at the genuinely-broken ASCII default
    that locale-misconfigured production machines would exhibit.
    """

    def _run_under_c_locale(self, script: str, extra_env: dict | None = None) -> "subprocess.CompletedProcess[str]":
        import os
        import subprocess
        import sys

        env = {**os.environ, "LANG": "C", "LC_ALL": "C", "PYTHONUTF8": "0"}
        if extra_env:
            env.update(extra_env)
        return subprocess.run(
            [sys.executable, "-c", script],
            env=env,
            capture_output=True,
            text=True,
        )

    def test_trad2_read_handles_non_ascii_under_c_locale(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            trad2 = Path(tmp) / "TRAD2"
            # Non-ASCII (é) in a description string — exactly the kind
            # of byte that breaks ASCII decoding.
            trad2.write_text(
                "ma_classe objet_u ma_classe 0 Base class décrite\n",
                encoding="utf-8",
            )
            result = self._run_under_c_locale(
                "from trustify.core.trad2_utilities import TRAD2Content; "
                f"TRAD2Content.BuildContentFromTRAD2({str(trad2)!r})",
            )
            self.assertEqual(
                result.returncode,
                0,
                f"BuildContentFromTRAD2 failed under LANG=C; the open() reading the TRAD2 "
                f"file must pass encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )

    def test_source_locations_json_read_handles_non_ascii_under_c_locale(self):
        import json
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            trad2 = Path(tmp) / "TRAD2"
            trad2.write_text(
                "ma_classe objet_u ma_classe 0 base\n",
                encoding="utf-8",
            )
            # The JSON includes a non-ASCII filesystem path (real
            # CEA HPC trees do exist with accented usernames).
            payload = {
                "version": 1,
                "entries": {"ma_classe": {"project": "trust", "path": "src/Décrit.cpp", "line": 1}},
            }
            json_path = Path(tmp) / "source_locations.json"
            json_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")

            result = self._run_under_c_locale(
                "from trustify.core.trad2_utilities import TRAD2Content; "
                f"TRAD2Content.BuildContentFromTRAD2({str(trad2)!r}, source_locations={str(json_path)!r})",
            )
            self.assertEqual(
                result.returncode,
                0,
                f"BuildContentFromTRAD2 failed under LANG=C reading source_locations.json; "
                f"the open() must pass encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )

    def test_trad2_write_handles_non_ascii_under_c_locale(self):
        # The writer must also force UTF-8: under LANG=C the default
        # ASCII codec rejects any non-ASCII char in a description.
        # The TRAD2Content is built manually so the block carries a
        # SourceLocation directly (toTRAD2 writes both the TRAD2 text
        # AND source_locations.json — the latter requires SourceLocation
        # objects, not the [path, lineno] lists produced by the loader).
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "TRAD2_out"
            script = (
                "from trustify.core.trad2_utilities import TRAD2Content, TRAD2Block; "
                "from trustify.core.source_location import SourceLocation; "
                "blk = TRAD2Block(); "
                "blk.name = 'ma_classe'; "
                "blk.name_base = 'objet_u'; "
                "blk.synos = ['ma_classe']; "
                "blk.mode = 'NO_BRACE'; "
                "blk.desc = 'D\\u00e9crit ma classe'; "
                "blk.info = SourceLocation(project='trust', path='src/X.cpp', line=1); "
                "c = TRAD2Content(); "
                "c.data = [blk]; "
                f"c.toTRAD2({str(out)!r})"
            )
            result = self._run_under_c_locale(script)
            self.assertEqual(
                result.returncode,
                0,
                f"toTRAD2 failed under LANG=C; the open() writing the TRAD2 file must pass "
                f"encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )

    def test_cpp_xd_scanner_handles_non_ascii_under_c_locale(self):
        # The C++/.xd scanner reads source files; many real TRUST .cpp
        # files have accented characters in comments or descriptions.
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            cpp = Path(tmp) / "fixture.cpp"
            cpp.write_text(
                "// XD ma_classe objet_u ma_classe 0 Décrit ma classe\n",
                encoding="utf-8",
            )
            result = self._run_under_c_locale(
                f"from trustify.core.trad2_utilities import TRAD2Content; TRAD2Content().scanOneCppFile({str(cpp)!r})",
            )
            self.assertEqual(
                result.returncode,
                0,
                f"scanOneCppFile failed under LANG=C; the open() reading C++/.xd sources "
                f"must pass encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )


if __name__ == "__main__":
    unittest.main()
