"""Focused unit tests for trustify.projects — Phase 2 of the reshaping."""

import os
import tempfile
import unittest
from pathlib import Path

from tests._baltik_fixtures import mock_workspace


class TestProjectResolution(unittest.TestCase):
    def test_resolve_projects_prepends_trust_root(self):
        from trustify.projects import resolve_projects

        result = resolve_projects(
            projects=["/abs/baltik_a", "/abs/baltik_b"],
            trust_root="/abs/trust",
        )
        self.assertEqual([p.path for p in result], [Path("/abs/trust"), Path("/abs/baltik_a"), Path("/abs/baltik_b")])

    def test_resolve_projects_no_trust_root(self):
        from trustify.projects import resolve_projects

        result = resolve_projects(projects=["/abs/standalone"], trust_root=None)
        self.assertEqual([str(p.path) for p in result], ["/abs/standalone"])

    def test_resolve_projects_empty_raises(self):
        from trustify.projects import resolve_projects

        with self.assertRaises(ValueError):
            resolve_projects(projects=None, trust_root=None)

    def test_resolve_projects_expands_user_and_envvar(self):
        from trustify.projects import resolve_projects

        os.environ["TST_BALTIK_ROOT"] = "/abs/expanded_baltik"
        try:
            result = resolve_projects(
                projects=["${TST_BALTIK_ROOT}/sub"],
                trust_root="~/trust_home",
            )
        finally:
            del os.environ["TST_BALTIK_ROOT"]
        paths = [str(p.path) for p in result]
        self.assertTrue(paths[0].endswith("/trust_home"), paths[0])
        self.assertEqual(paths[1], "/abs/expanded_baltik/sub")


class TestEffectiveProjectArgs(unittest.TestCase):
    """`effective_trust_root` / `effective_projects` apply the env-var
    fallbacks used everywhere in `api.*` and the CLI. They are the
    canonical site for `$TRUST_ROOT` and `$project_directory` lookups.

    Every test save/restores the env keys it touches — other test
    modules in this suite (notably `test_rw_full_datasets`) need a real
    `$TRUST_ROOT` in their environment.
    """

    # Lowercase env var name dictated by TRUST's build environment.
    _PD_KEY = "project_directory"

    @staticmethod
    def _patch_env(key, value):
        """Save current env[key] (or sentinel) and set/unset to value. Returns
        a callable that restores the original state.
        """
        sentinel = object()
        original = os.environ.get(key, sentinel)

        def _restore():
            if original is sentinel:
                os.environ.pop(key, None)
            else:
                os.environ[key] = original

        if value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = value
        return _restore

    def test_effective_trust_root_returns_arg_when_given(self):
        from trustify.projects import effective_trust_root

        # Explicit arg always wins over env var.
        restore = self._patch_env("TRUST_ROOT", "/from/env")
        try:
            self.assertEqual(effective_trust_root("/from/arg"), "/from/arg")
        finally:
            restore()

    def test_effective_trust_root_falls_back_to_env_var(self):
        from trustify.projects import effective_trust_root

        restore = self._patch_env("TRUST_ROOT", "/from/env")
        try:
            self.assertEqual(effective_trust_root(None), "/from/env")
        finally:
            restore()

    def test_effective_trust_root_returns_none_when_neither(self):
        from trustify.projects import effective_trust_root

        restore = self._patch_env("TRUST_ROOT", None)
        try:
            self.assertIsNone(effective_trust_root(None))
        finally:
            restore()

    def test_effective_projects_returns_arg_when_given(self):
        from trustify.projects import effective_projects

        restore = self._patch_env(self._PD_KEY, "/from/env_baltik")
        try:
            self.assertEqual(effective_projects(["/from/arg"]), ["/from/arg"])
        finally:
            restore()

    def test_effective_projects_falls_back_to_project_directory_env(self):
        """In TRUST's baltik build environment, $project_directory points
        at the active baltik root — treat it as a default `--projects`
        entry, symmetric with $TRUST_ROOT for `--trust-root`. Only fires
        when the caller passed None (i.e. "didn't specify anything").
        """
        from trustify.projects import effective_projects

        restore = self._patch_env(self._PD_KEY, "/from/env_baltik")
        try:
            self.assertEqual(effective_projects(None), ["/from/env_baltik"])
        finally:
            restore()

    def test_effective_projects_empty_list_suppresses_env_fallback(self):
        """`projects=[]` is the CLI's explicit "no baltik overlay" answer
        (e.g. when the anchor-inside-trust_root rule fired). It must NOT
        fall back to `$project_directory`, or that suppression silently
        unravels at the API layer.
        """
        from trustify.projects import effective_projects

        restore = self._patch_env(self._PD_KEY, "/from/env_baltik")
        try:
            self.assertEqual(effective_projects([]), [])
        finally:
            restore()

    def test_effective_projects_returns_none_when_neither(self):
        from trustify.projects import effective_projects

        restore = self._patch_env(self._PD_KEY, None)
        try:
            self.assertIsNone(effective_projects(None))
        finally:
            restore()


class TestConfigDiscovery(unittest.TestCase):
    def test_discover_config_walks_up_to_find_trustify_json(self):
        from trustify.projects import TrustifyConfig, discover_config

        with tempfile.TemporaryDirectory() as d:
            top = Path(d).resolve()
            cfg = top / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/some/baltik"]}}')
            sub = top / "a" / "b" / "c"
            sub.mkdir(parents=True)
            found = discover_config(sub)
            self.assertIsInstance(found, TrustifyConfig)
            self.assertEqual(found.source_path, cfg)
            self.assertEqual(found.workspace.projects, ["/some/baltik"])

    def test_discover_config_walks_up_from_file_path(self):
        """Anchor may be a file rather than a directory — used by the CLI
        when checking a dataset that lives outside the cwd."""
        from trustify.projects import discover_config

        with tempfile.TemporaryDirectory() as d:
            top = Path(d).resolve()
            cfg = top / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/some/baltik"]}}')
            sub = top / "a"
            sub.mkdir()
            data_file = sub / "foo.data"
            data_file.write_text("")
            found = discover_config(data_file)
            self.assertEqual(found.source_path, cfg)

    def test_discover_config_returns_none_when_absent(self):
        from trustify.projects import discover_config

        with tempfile.TemporaryDirectory() as d:
            # Use the tempdir itself as the anchor; it won't reach our
            # real .trustify.json because the walk stops at /.
            sub = Path(d) / "x"
            sub.mkdir()
            # The walk WILL still reach / and might find a real config;
            # if so, just assert that one was found and skip. Otherwise
            # None.
            found = discover_config(sub)
            if found is not None:
                self.skipTest(f"environment has a config at {found.source_path}")
            self.assertIsNone(found)

    def test_load_config_returns_dict(self):
        # Back-compat shim — new code uses TrustifyConfig.from_file.
        from trustify.projects import load_config

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/x", "/y"], "trust_root": "/z"}}')
            data = load_config(cfg)
            self.assertEqual(data["workspace"]["projects"], ["/x", "/y"])
            self.assertEqual(data["workspace"]["trust_root"], "/z")

    def test_load_config_rejects_invalid_keys_with_clear_error(self):
        from trustify.projects import ConfigError, load_config

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": "not-a-list"}}')
            with self.assertRaises(ConfigError) as ctx:
                load_config(cfg)
            self.assertIn("projects", str(ctx.exception))


class TestTrustifyConfigDataclass(unittest.TestCase):
    """`TrustifyConfig.from_file` parses the three recognized top-level
    sections (`workspace`, `cache`, `lsp`) into a typed structure and
    keeps the original payload in `raw` for forward compatibility.
    """

    def test_parses_all_recognized_keys(self):
        from trustify.projects import TrustifyConfig

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text(
                '{"workspace": {"projects": ["/a", "/b"], "trust_root": "/t", '
                '"auto_resolve_dependencies": false}, '
                '"cache": {"auto_trim": false}, '
                '"lsp": {"schema_dir": "/s", "enum_dedup_threshold": 8}}'
            )
            c = TrustifyConfig.from_file(cfg)
            self.assertEqual(c.workspace.projects, ["/a", "/b"])
            self.assertEqual(c.workspace.trust_root, "/t")
            self.assertFalse(c.workspace.auto_resolve_dependencies)
            self.assertFalse(c.cache.auto_trim)
            self.assertEqual(c.lsp.schema_dir, "/s")
            self.assertEqual(c.lsp.enum_dedup_threshold, 8)
            self.assertEqual(c.source_path, cfg.resolve())

    def test_defaults_when_keys_absent(self):
        from trustify.projects import TrustifyConfig

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text("{}")
            c = TrustifyConfig.from_file(cfg)
            self.assertIsNone(c.workspace.projects)
            self.assertIsNone(c.workspace.trust_root)
            self.assertTrue(c.workspace.auto_resolve_dependencies)  # default True
            self.assertTrue(c.cache.auto_trim)  # default True
            self.assertIsNone(c.lsp.schema_dir)
            self.assertIsNone(c.lsp.enum_dedup_threshold)

    def test_to_api_kwargs_returns_projects_and_trust_root(self):
        from trustify.projects import TrustifyConfig

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/a"], "trust_root": "/t"}}')
            c = TrustifyConfig.from_file(cfg)
            self.assertEqual(c.to_api_kwargs(), {"projects": ["/a"], "trust_root": "/t"})

    def test_raw_payload_preserved_for_forward_compat(self):
        """Unrecognized keys (e.g. future LSP-side ones) survive in `raw`
        so consumers can pick them up without going through this module.
        """
        from trustify.projects import TrustifyConfig

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/a"]}, "future_key": {"sub": 42}}')
            c = TrustifyConfig.from_file(cfg)
            self.assertEqual(c.raw.get("future_key"), {"sub": 42})

    def test_rejects_lsp_wrong_types(self):
        from trustify.projects import ConfigError, TrustifyConfig

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text('{"lsp": {"enum_dedup_threshold": "not-an-int"}}')
            with self.assertRaises(ConfigError) as ctx:
                TrustifyConfig.from_file(cfg)
            self.assertIn("enum_dedup_threshold", str(ctx.exception))

    def test_rejects_workspace_wrong_types(self):
        """Each `workspace.*` field carries a type — wrong types raise
        ConfigError with the offending key name in the message."""
        from trustify.projects import ConfigError, TrustifyConfig

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text('{"workspace": {"auto_resolve_dependencies": "yes please"}}')
            with self.assertRaises(ConfigError) as ctx:
                TrustifyConfig.from_file(cfg)
            self.assertIn("workspace.auto_resolve_dependencies", str(ctx.exception))


class TestValidateProjects(unittest.TestCase):
    """`validate_projects` requires every non-trust_root project to carry
    a top-level `project.cfg` baltik marker. The trust_root entry is
    always exempt — TRUST itself has no project.cfg at its root.
    """

    def test_passes_when_every_project_has_project_cfg(self):
        from trustify.projects import resolve_projects, validate_projects

        with mock_workspace() as ws:
            ba = ws.baltik("baltik_a")
            bb = ws.baltik("baltik_b")
            resolved = resolve_projects(projects=[str(ba), str(bb)], trust_root=None)
            validate_projects(resolved)  # no raise

    def test_trust_root_is_exempt(self):
        from trustify.projects import resolve_projects, validate_projects

        with mock_workspace() as ws:
            # trust_root carries no project.cfg — must be accepted.
            trust = ws.root / "trust"
            trust.mkdir()
            resolved = resolve_projects(projects=None, trust_root=str(trust))
            validate_projects(resolved)

    def test_raises_when_project_lacks_project_cfg(self):
        from trustify.projects import ConfigError, resolve_projects, validate_projects

        with mock_workspace() as ws:
            bare = ws.root / "not_a_baltik"
            bare.mkdir()  # no project.cfg
            resolved = resolve_projects(projects=[str(bare)], trust_root=None)
            with self.assertRaises(ConfigError) as ctx:
                validate_projects(resolved)
            msg = str(ctx.exception)
            self.assertIn("not_a_baltik", msg)
            self.assertIn("project.cfg", msg)
            self.assertIn("--trust-root", msg)
            self.assertIn("baltik_configure", msg)

    def test_raises_for_first_offending_project_only(self):
        """Stop at the first bad project — the user will fix one at a
        time anyway."""
        from trustify.projects import ConfigError, resolve_projects, validate_projects

        with mock_workspace() as ws:
            good = ws.baltik("good")
            bad = ws.root / "bad"
            bad.mkdir()
            resolved = resolve_projects(projects=[str(good), str(bad)], trust_root=None)
            with self.assertRaises(ConfigError) as ctx:
                validate_projects(resolved)
            # Only "bad" should be mentioned in the error.
            self.assertIn("bad", str(ctx.exception))

    def test_raises_when_baltik_lacks_src_dir(self):
        """A baltik with a `project.cfg` but no `src/` subdirectory used
        to surface as a bare `Exception("Directory '...' not found!")`
        from the source scanner, mid-pipeline. validate_projects now
        rejects this up-front with a clear ConfigError pointing at the
        offending baltik.
        """
        from trustify.projects import ConfigError, resolve_projects, validate_projects

        with mock_workspace() as ws:
            baltik = ws.root / "no_src_baltik"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname : no_src_baltik\n")
            # NB: no src/ directory created.
            resolved = resolve_projects(projects=[str(baltik)], trust_root=None)
            with self.assertRaises(ConfigError) as ctx:
                validate_projects(resolved)
            msg = str(ctx.exception)
            self.assertIn("no_src_baltik", msg)
            self.assertIn("src/", msg)
            self.assertIn("baltik_configure", msg)

    def test_passes_with_empty_src_dir(self):
        """An empty `src/` is valid — the baltik just contributes no XD
        tags. Make sure validate_projects doesn't reject it.
        """
        from trustify.projects import resolve_projects, validate_projects

        with mock_workspace() as ws:
            baltik = ws.root / "empty_src_baltik"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname : empty_src_baltik\n")
            (baltik / "src").mkdir()
            resolved = resolve_projects(projects=[str(baltik)], trust_root=None)
            validate_projects(resolved)  # no raise


class TestConfigBesideProjectCfg(unittest.TestCase):
    """`.trustify.json` that co-exists with a `project.cfg` must not
    redeclare `projects` or `trust_root` — the project is already
    auto-detected from the baltik marker. Surfaced as a ConfigError so
    the misconfiguration doesn't shadow the auto-detected baltik
    silently (the previous behaviour).
    """

    def test_rejects_projects_key_inside_baltik(self):
        from trustify.projects import ConfigError, TrustifyConfig

        with mock_workspace() as ws:
            baltik = ws.baltik("the_baltik")
            cfg = baltik / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/opt/other-baltik"]}}')
            with self.assertRaises(ConfigError) as ctx:
                TrustifyConfig.from_file(cfg)
            msg = str(ctx.exception)
            self.assertIn("project.cfg", msg)
            self.assertIn("projects", msg)
            # Must point at the actual offending file in the error.
            self.assertIn(str(cfg.resolve()), msg)

    def test_rejects_trust_root_key_inside_baltik(self):
        from trustify.projects import ConfigError, TrustifyConfig

        with mock_workspace() as ws:
            baltik = ws.baltik("the_baltik")
            cfg = baltik / ".trustify.json"
            cfg.write_text('{"workspace": {"trust_root": "/opt/some/trust"}}')
            with self.assertRaises(ConfigError) as ctx:
                TrustifyConfig.from_file(cfg)
            self.assertIn("trust_root", str(ctx.exception))

    def test_rejects_both_keys_inside_baltik(self):
        from trustify.projects import ConfigError, TrustifyConfig

        with mock_workspace() as ws:
            baltik = ws.baltik("the_baltik")
            cfg = baltik / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/x"], "trust_root": "/y"}}')
            with self.assertRaises(ConfigError) as ctx:
                TrustifyConfig.from_file(cfg)
            msg = str(ctx.exception)
            self.assertIn("projects", msg)
            self.assertIn("trust_root", msg)

    def test_lsp_only_config_inside_baltik_is_fine(self):
        """The canonical use case for a co-located .trustify.json: LSP /
        cache settings without redeclaring the project list.
        """
        from trustify.projects import TrustifyConfig

        with mock_workspace() as ws:
            baltik = ws.baltik("the_baltik")
            cfg = baltik / ".trustify.json"
            cfg.write_text('{"lsp": {"enum_dedup_threshold": 6}, "cache": {"auto_trim": false}}')
            c = TrustifyConfig.from_file(cfg)
            self.assertEqual(c.lsp.enum_dedup_threshold, 6)
            self.assertFalse(c.cache.auto_trim)

    def test_projects_key_outside_baltik_is_fine(self):
        """Same .trustify.json content, but no project.cfg next to it →
        fine. The check is strictly about coexistence."""
        from trustify.projects import TrustifyConfig

        with tempfile.TemporaryDirectory() as d:
            cfg = Path(d) / ".trustify.json"
            cfg.write_text('{"workspace": {"projects": ["/opt/baltik"], "trust_root": "/opt/trust"}}')
            c = TrustifyConfig.from_file(cfg)
            self.assertEqual(c.workspace.projects, ["/opt/baltik"])


class TestExpandDependencies(unittest.TestCase):
    """`expand_dependencies` recursively expands each project's
    `[dependencies]` from project.cfg, mirroring the behaviour of
    `baltik_configure`. Post-order ensures every dep appears before the
    baltik that declares it (overlay-correct in `scanSourceFiles`).
    """

    def test_returns_input_unchanged_when_no_dependencies(self):
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            ba = ws.baltik("ba")
            result = expand_dependencies([str(ba)])
            self.assertEqual(result, [str(ba)])

    def test_post_order_dep_before_declaring_baltik(self):
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            dep = ws.baltik("dep")
            main = ws.baltik("main", deps={"dep": dep})
            result = expand_dependencies([str(main)])
            self.assertEqual(result, [str(dep), str(main)])

    def test_transitive_expansion(self):
        """A → B → C should produce [C, B, A] (overlay: A overrides B
        overrides C in scanSourceFiles)."""
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            c = ws.baltik("c")
            b = ws.baltik("b", deps={"c": c})
            a = ws.baltik("a", deps={"b": b})
            self.assertEqual(expand_dependencies([str(a)]), [str(c), str(b), str(a)])

    def test_deduplicates_diamond_dependency(self):
        """A → B → D, A → C → D should yield [D, B, C, A] (D appears
        only once even though both B and C list it)."""
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            shared = ws.baltik("shared")
            left = ws.baltik("left", deps={"shared": shared})
            right = ws.baltik("right", deps={"shared": shared})
            top = ws.baltik("top", deps={"left": left, "right": right})
            result = expand_dependencies([str(top)])
            self.assertEqual(result.count(str(shared)), 1)
            # Ordering: shared is reached first via `left`, so it appears before
            # both left and right.
            self.assertLess(result.index(str(shared)), result.index(str(left)))
            self.assertLess(result.index(str(left)), result.index(str(top)))
            self.assertLess(result.index(str(right)), result.index(str(top)))

    def test_relative_path_resolved_against_baltik_root(self):
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            dep = ws.baltik("dep")
            main = ws.baltik("main", deps={"dep": "../dep"})
            result = expand_dependencies([str(main)])
            self.assertIn(str(dep), result)

    def test_pwd_backtick_substituted_with_baltik_root(self):
        """baltik_configure evaluates ``\\`pwd\\``` as the declaring
        baltik's directory. Mirror that here so .cfg files copy-pasted
        from a working baltik work without modification."""
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            dep = ws.baltik("dep")
            # Declare with the literal `pwd` form; trustify must
            # substitute it with `main`'s directory.
            main = ws.baltik("main", deps={"dep": "`pwd`/../dep"})
            result = expand_dependencies([str(main)])
            self.assertIn(str(dep), result)

    def test_raises_when_dep_dir_missing(self):
        from trustify.projects import ConfigError, expand_dependencies

        with mock_workspace() as ws:
            main = ws.baltik("main", deps={"ghost": str(ws.root / "nope")})
            with self.assertRaises(ConfigError) as ctx:
                expand_dependencies([str(main)])
            self.assertIn("ghost", str(ctx.exception))
            self.assertIn("does not exist", str(ctx.exception))

    def test_raises_when_dep_name_mismatch(self):
        """If main's project.cfg declares `dep_X: ./somewhere` but
        ./somewhere/project.cfg sets `[description].name = other`, that
        is rejected (mirrors baltik_configure's invalid_dependency_name_error).
        """
        from trustify.projects import ConfigError, expand_dependencies

        with mock_workspace() as ws:
            dep = ws.baltik("actual_name")
            main = ws.baltik("main", deps={"wrong_name": dep})
            with self.assertRaises(ConfigError) as ctx:
                expand_dependencies([str(main)])
            msg = str(ctx.exception)
            self.assertIn("wrong_name", msg)
            self.assertIn("actual_name", msg)

    def test_mixed_case_dep_name_preserved(self):
        """Dependency keys in `[dependencies]` must keep their original
        case (baltik_configure parses with sed/awk and matches case-
        sensitively). Regression for the configparser default that
        lowercased option keys and made `TrioCFD : ...` mismatch the
        dep's own `[description].name = TrioCFD`.
        """
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            dep = ws.baltik("TrioCFD")
            main = ws.baltik("main", deps={"TrioCFD": dep})
            result = expand_dependencies([str(main)])
            self.assertEqual(result, [str(dep), str(main)])

    def test_strips_inline_hash_comment_from_dep_path(self):
        """Mirror `baltik_configure`'s `sed 's/#.*//'`: a trailing
        `# comment` after a dependency path must not bleed into the
        path. Without `inline_comment_prefixes`, the value would be
        `<dep_dir> # main dep` and the path resolution would fail.
        """
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            dep = ws.baltik("dep")
            main = ws.baltik("main")
            (main / "project.cfg").write_text(f"[description]\nname : main\n\n[dependencies]\ndep : {dep} # main dep\n")
            result = expand_dependencies([str(main)])
            self.assertEqual(result, [str(dep), str(main)])

    def test_strips_outer_quotes_from_dep_path(self):
        """`baltik_configure`'s `eval "echo $entry"` strips surrounding
        `"..."` / `'...'` quotes. Without this, a literal-quoted path
        like `dep = "/abs/path"` would be passed through verbatim and
        fail to resolve.
        """
        from trustify.projects import expand_dependencies

        with mock_workspace() as ws:
            dep = ws.baltik("dep")
            main = ws.baltik("main")
            # Quoted dep path; trustify must strip the quotes the same
            # way baltik's `eval` does.
            (main / "project.cfg").write_text(f'[description]\nname : main\n\n[dependencies]\ndep : "{dep}"\n')
            result = expand_dependencies([str(main)])
            self.assertEqual(result, [str(dep), str(main)])

    def test_rejects_two_baltiks_with_same_description_name(self):
        """`baltik_configure`'s `check_dependencies` (in
        `baltik_dependencies_management`) errors out when the same dep
        name resolves to two distinct paths — `too_many_paths_error`.
        Trustify must mirror that: two physically distinct baltiks both
        claiming `[description].name = X` is a misconfiguration.
        """
        from trustify.projects import ConfigError, expand_dependencies

        with mock_workspace() as ws:
            # Both baltiks declare the same `[description].name`, despite
            # living at different paths. Dependency-graph-wise, "shared"
            # is now ambiguous.
            a = ws.baltik("shared")
            b_dir = ws.root / "other_shared"
            b_dir.mkdir()
            (b_dir / "project.cfg").write_text("[description]\nname : shared\n")
            with self.assertRaises(ConfigError) as ctx:
                expand_dependencies([str(a), str(b_dir)])
            msg = str(ctx.exception)
            self.assertIn("shared", msg)
            self.assertIn("too_many_paths_error", msg)

    def test_none_and_empty_pass_through(self):
        from trustify.projects import expand_dependencies

        self.assertIsNone(expand_dependencies(None))
        self.assertEqual(expand_dependencies([]), [])


class TestDetectBaltikRoot(unittest.TestCase):
    """`detect_baltik_root` walks up looking for `project.cfg` and
    returns the directory that contains it (or None). This is the
    third-priority `--projects` source in the CLI's four-signal
    precedence ladder.
    """

    def test_returns_baltik_dir_when_project_cfg_present(self):
        from trustify.projects import detect_baltik_root

        with mock_workspace() as ws:
            baltik = ws.baltik("my_baltik")
            sub = baltik / "src" / "Kernel"
            sub.mkdir(parents=True)
            self.assertEqual(detect_baltik_root(sub), baltik)

    def test_returns_none_when_no_project_cfg(self):
        from trustify.projects import detect_baltik_root

        with mock_workspace() as ws:
            sub = ws.root / "no_baltik"
            sub.mkdir()
            # The walk may still find a project.cfg further up in some
            # environments — skip in that case to stay environment-independent.
            found = detect_baltik_root(sub)
            if found is not None and not str(found).startswith(str(ws.root)):
                self.skipTest(f"environment has a project.cfg at {found}")
            self.assertIsNone(found)

    def test_accepts_file_anchor(self):
        from trustify.projects import detect_baltik_root

        with mock_workspace() as ws:
            baltik = ws.baltik("b")
            data_file = baltik / "datasets" / "foo.data"
            data_file.parent.mkdir()
            data_file.write_text("")
            self.assertEqual(detect_baltik_root(data_file), baltik)


class TestProjectNameFor(unittest.TestCase):
    def test_reads_name_from_project_cfg_description_section(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "lower-case-dir"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname       : TrioCFD\nauthor     : CEA\n")
            self.assertEqual(project_name_for(baltik), "TrioCFD")

    def test_strips_whitespace_around_name_value(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "b"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname:    MyName   \n")
            self.assertEqual(project_name_for(baltik), "MyName")

    def test_falls_back_to_basename_when_no_project_cfg(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "fallback_name"
            baltik.mkdir()
            self.assertEqual(project_name_for(baltik), "fallback_name")

    def test_falls_back_to_basename_when_description_section_missing(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "fallback_name"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[dependencies]\nname : ignored\n")
            self.assertEqual(project_name_for(baltik), "fallback_name")

    def test_falls_back_to_basename_when_name_field_missing(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "fallback_name"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nauthor : CEA\n")
            self.assertEqual(project_name_for(baltik), "fallback_name")

    def test_falls_back_to_basename_when_name_value_empty(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "fallback_name"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname :   \n")
            self.assertEqual(project_name_for(baltik), "fallback_name")

    def test_falls_back_to_basename_when_project_cfg_malformed(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "fallback_name"
            baltik.mkdir()
            # Garbage that configparser can't read.
            (baltik / "project.cfg").write_text("][not a section\n=garbage\n")
            self.assertEqual(project_name_for(baltik), "fallback_name")

    def test_accepts_string_and_path_input(self):
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "b"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname : Foo\n")
            self.assertEqual(project_name_for(str(baltik)), "Foo")
            self.assertEqual(project_name_for(baltik), "Foo")

    def test_strips_inline_hash_comment_from_name(self):
        """`baltik_configure` runs `sed 's/#.*//'` on every line of
        project.cfg, so a trailing `# comment` after a value is dropped.
        ConfigParser's default keeps it; we opt into
        `inline_comment_prefixes=("#",)` so trustify matches baltik.
        """
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "b"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nname = toto #this is the name\n")
            self.assertEqual(project_name_for(baltik), "toto")

    def test_default_section_does_not_leak_into_description(self):
        """ConfigParser's default `[DEFAULT]` section is special — its
        keys get injected into every other section. `baltik_configure`'s
        sed pipeline does NOT do this; `[DEFAULT]` is just an (unused)
        regular section. Without redirecting `default_section=`, a
        `[DEFAULT]\\nname = sneaky` would make trustify return `sneaky`
        even when `[description]` has no `name` — and baltik would error
        out.
        """
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "fallback_name"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[DEFAULT]\nname = sneaky\n\n[description]\nauthor = CEA\n")
            # No `name` in [description] → fallback to basename, NOT to
            # the [DEFAULT].name leak.
            self.assertEqual(project_name_for(baltik), "fallback_name")

    def test_uppercase_name_key_does_not_match(self):
        """`baltik_configure`'s `sed -n "/name/p"` is case-sensitive, so
        `Name = Foo` is not recognised. With `optionxform = str`, trustify
        also stores keys case-sensitively, so the lowercase lookup misses.
        Mirrors baltik's behaviour — both fall back when the key is
        miscased.
        """
        from trustify.projects import project_name_for

        with tempfile.TemporaryDirectory() as d:
            baltik = Path(d) / "fallback_name"
            baltik.mkdir()
            (baltik / "project.cfg").write_text("[description]\nName = Foo\n")
            self.assertEqual(project_name_for(baltik), "fallback_name")


class TestResolveDepPathQuoteStripping(unittest.TestCase):
    """Audit 6.9: `_resolve_dep_path` used to strip outer quotes
    whenever first and last char matched, no matter what the inner
    content was — so `"foo"bar"` (first/last both `"`, unmatched `"`
    in the middle) silently became `foo"bar`. Tighten the strip to
    require the inner content NOT to contain the same quote character.
    A well-formed `"some path"` still strips.
    """

    def test_well_formed_quoted_path_still_strips(self):
        from trustify.projects import _resolve_dep_path

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = _resolve_dep_path('"foo bar"', root)
            self.assertEqual(out.name, "foo bar")

    def test_well_formed_single_quoted_path_still_strips(self):
        from trustify.projects import _resolve_dep_path

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = _resolve_dep_path("'foo bar'", root)
            self.assertEqual(out.name, "foo bar")

    def test_unbalanced_inner_quote_is_not_stripped(self):
        from trustify.projects import _resolve_dep_path

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # `"foo"bar"` — three `"` characters; the previous "first
            # and last match → strip" produced `foo"bar` which is even
            # more malformed. Leave the value intact so the caller's
            # path-doesn't-exist failure is faithful to the input.
            out = _resolve_dep_path('"foo"bar"', root)
            self.assertEqual(out.name, '"foo"bar"')

    def test_mixed_quote_in_inner_content_does_not_block_strip(self):
        # Outer `"`, inner `'` is fine — different quote chars don't
        # interfere with each other's balance.
        from trustify.projects import _resolve_dep_path

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = _resolve_dep_path("\"o'br'\"", root)
            self.assertEqual(out.name, "o'br'")


if __name__ == "__main__":
    unittest.main()
