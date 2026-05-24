"""End-to-end CLI tests for trustify.cli — Phase 2 of the reshaping."""

import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests._baltik_fixtures import mock_workspace

# argcomplete is an opt-in extra (pip install trustify[optional]); the
# install-completion subcommand exists either way but only emits a
# usable script when the library is importable. Tests that exercise
# the happy path get conditionally skipped when it's missing.
_HAS_ARGCOMPLETE = importlib.util.find_spec("argcomplete") is not None


def _run(*args, **kwargs):
    """Run `python -m trustify.cli ...` and return CompletedProcess."""
    cmd = [sys.executable, "-m", "trustify.cli", *args]
    env = kwargs.pop("env", os.environ.copy())
    return subprocess.run(cmd, env=env, capture_output=True, text=True, **kwargs)


class TestCliEntryPoint(unittest.TestCase):
    def test_no_args_prints_help_and_exits_nonzero(self):
        result = _run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("trustify", (result.stderr + result.stdout).lower())

    def test_version_flag(self):
        result = _run("--version")
        self.assertEqual(result.returncode, 0)
        self.assertRegex(result.stdout.strip(), r"\S+")


class TestCliGenerate(unittest.TestCase):
    def test_generate_writes_files_to_out(self):
        with tempfile.TemporaryDirectory() as d:
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "generate_schema", "--out", d)
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertTrue((Path(d) / "trustify_gen.py").exists())


class TestCliCheck(unittest.TestCase):
    def test_check_returns_zero_for_valid_dataset(self):
        ds = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data")
        result = _run("--trust-root", os.environ["TRUST_ROOT"], "check", ds)
        self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")

    def test_check_missing_file_errors_fast(self):
        with tempfile.TemporaryDirectory() as d:
            missing = Path(d) / "nope.data"
            result = _run("check", str(missing))
            self.assertEqual(result.returncode, 2)
            self.assertIn("check: no such file", result.stderr)
            self.assertIn(str(missing), result.stderr)

    def test_check_skips_trustify_not_dataset_by_default(self):
        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT #\nnot valid TRUST\n")
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "check", str(data))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("SKIPPED", result.stdout)
            self.assertIn("TRUSTIFY NOT", result.stdout)

    def test_check_no_skip_forces_validation_of_marked_dataset(self):
        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT #\nnot valid TRUST\n")
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "check", "--no-skip", str(data))
            self.assertEqual(result.returncode, 1)
            self.assertIn("FAILED", result.stdout)

    def test_check_no_skip_warns_on_obsolete_marker(self):
        # A previously-marked dataset that now passes the check should
        # emit a stderr warning suggesting the marker can be removed.
        valid = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked = Path(d) / "marked.data"
            marked.write_text("# TRUSTIFY NOT #\n" + Path(valid).read_text())
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "check", "--no-skip", str(marked))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("PASSED", result.stdout)
            self.assertIn("obsolete", result.stderr.lower())
            self.assertIn("TRUSTIFY NOT", result.stderr)

    def test_batch_check_no_skip_warns_on_obsolete_marker(self):
        # batch-check should list every obsolete marker at the end,
        # to stderr, with one line per offending file.
        valid = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked_ok = Path(d) / "marked_ok.data"
            marked_ok.write_text("# TRUSTIFY NOT #\n" + Path(valid).read_text())
            marked_bad = Path(d) / "marked_bad.data"
            marked_bad.write_text("# TRUSTIFY NOT #\nnot valid TRUST\n")
            result = _run(
                "--trust-root", os.environ["TRUST_ROOT"], "batch-check", "--no-skip", str(marked_ok), str(marked_bad)
            )
            # One legitimate failure, one obsolete-marker pass.
            self.assertEqual(result.returncode, 1)
            self.assertIn("1 obsolete-marker", result.stdout)
            self.assertIn("obsolete", result.stderr.lower())
            self.assertIn(str(marked_ok), result.stderr)

    def test_batch_check_summary_orders_failed_last(self):
        # The summary block should sort entries so FAILED cases land
        # at the bottom of the output — they're the actionable items
        # and putting them last keeps them in the terminal viewport.
        # SKIPPED entries come first, FAILED entries last.
        with tempfile.TemporaryDirectory() as d:
            # Mix file names so submission order interleaves SKIPPED
            # and FAILED — without an explicit sort we'd see them
            # alternating in the summary.
            (Path(d) / "a_skip.data").write_text("# TRUSTIFY NOT #\nirrelevant\n")
            (Path(d) / "b_fail.data").write_text("not valid TRUST\n")
            (Path(d) / "c_skip.data").write_text("# TRUSTIFY NOT #\nirrelevant\n")
            (Path(d) / "d_fail.data").write_text("not valid TRUST\n")
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "batch-check", "-q", d)
            self.assertEqual(result.returncode, 1)
            # Find the positions of the per-file lines in the summary.
            out = result.stdout
            skip_positions = [out.find("a_skip.data"), out.find("c_skip.data")]
            fail_positions = [out.find("b_fail.data"), out.find("d_fail.data")]
            for p in skip_positions + fail_positions:
                self.assertNotEqual(p, -1, f"missing entry in:\n{out}")
            # Every SKIPPED line must come before every FAILED line.
            self.assertLess(max(skip_positions), min(fail_positions), f"FAILED not after SKIPPED in:\n{out}")

    def test_check_renders_justification_on_skipped_line(self):
        # `trustify check` SKIPPED line appends ` — <justification>`
        # when the tag carries one.
        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT: sed template, not real TRUST #\nnot valid TRUST\n")
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "check", str(data))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("SKIPPED", result.stdout)
            self.assertIn("sed template, not real TRUST", result.stdout)

    def test_check_renders_justification_on_obsolete_warning(self):
        # `trustify check --no-skip` WARNING includes the justification.
        valid = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked = Path(d) / "marked.data"
            marked.write_text("# TRUSTIFY NOT: kept for legacy fixture #\n" + Path(valid).read_text())
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "check", "--no-skip", str(marked))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("obsolete", result.stderr.lower())
            self.assertIn("kept for legacy fixture", result.stderr)

    def test_batch_check_renders_justification_in_summary(self):
        # `trustify batch-check` SKIPPED entries in the summary include
        # the justification.
        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT: sed template, not real TRUST #\nnot valid TRUST\n")
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "batch-check", str(data))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("[SKIPPED]", result.stdout)
            self.assertIn("sed template, not real TRUST", result.stdout)

    def test_batch_check_renders_justification_on_obsolete_and_summary_out(self):
        # `trustify batch-check --no-skip --summary-out FILE`:
        #  - the obsolete WARNING block on stderr includes the justification
        #  - the --summary-out file mirrors the same lines.
        valid = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked_ok = Path(d) / "marked_ok.data"
            marked_ok.write_text("# TRUSTIFY NOT: kept for legacy fixture #\n" + Path(valid).read_text())
            summary_out = Path(d) / "summary.txt"
            result = _run(
                "--trust-root",
                os.environ["TRUST_ROOT"],
                "batch-check",
                "--no-skip",
                str(marked_ok),
                "--summary-out",
                str(summary_out),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("[obsolete]", result.stderr)
            self.assertIn("kept for legacy fixture", result.stderr)
            self.assertTrue(summary_out.is_file())
            file_text = summary_out.read_text()
            self.assertIn("[obsolete]", file_text)
            self.assertIn("kept for legacy fixture", file_text)

    def test_batch_check_summary_out_file_has_no_ansi_codes(self):
        """`--summary-out` is a log file — it must NEVER contain ANSI
        escape sequences, even when the terminal would render them
        (colors are applied at print time, not during buffer
        construction). Subprocess tests run with a non-TTY stdout
        anyway, but force `TERM=xterm-256color` and unset `NO_COLOR`
        so a buggy implementation that colors the buffer would surface
        immediately.
        """
        with tempfile.TemporaryDirectory() as d:
            # Mix of SKIPPED + PASSED so the summary block actually has
            # content to render (FAILED entries would be ideal too, but
            # we don't have a guaranteed-broken dataset on hand).
            skipped = Path(d) / "opted_out.data"
            skipped.write_text("# TRUSTIFY NOT: sed template #\nnot valid TRUST\n")
            summary_out = Path(d) / "summary.txt"
            env = os.environ.copy()
            env.pop("NO_COLOR", None)
            env["TERM"] = "xterm-256color"
            result = _run(
                "--trust-root",
                os.environ["TRUST_ROOT"],
                "batch-check",
                str(skipped),
                "--summary-out",
                str(summary_out),
                env=env,
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertTrue(summary_out.is_file())
            file_text = summary_out.read_text()
            # No ANSI escape sequences anywhere in the file.
            self.assertNotIn("\x1b", file_text, f"summary-out file contains ANSI escape codes:\n{file_text!r}")
            # Sanity: file does carry the SKIPPED line in plain text.
            self.assertIn("[SKIPPED]", file_text)

    def test_batch_check_missing_path_errors_fast(self):
        with tempfile.TemporaryDirectory() as d:
            good = Path(d) / "good.data"
            good.write_text("dimension 2\n")
            missing = Path(d) / "nope.data"
            result = _run("batch-check", str(good), str(missing))
            self.assertEqual(result.returncode, 2)
            self.assertIn("batch-check: no such file or directory", result.stderr)
            self.assertIn(str(missing), result.stderr)
            self.assertNotIn(str(good), result.stderr)

    def test_batch_check_empty_directory_short_circuits(self):
        with tempfile.TemporaryDirectory() as d:
            result = _run("batch-check", d)
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("no .data files found", result.stdout)

    def test_batch_check_exclude_all_short_circuits(self):
        # All discovered .data files are excluded → no schema generation,
        # no error, short-circuit. Doesn't need TRUST_ROOT.
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "a.data").write_text("dimension 2\n")
            (Path(d) / "b.data").write_text("dimension 2\n")
            result = _run("batch-check", d, "--exclude", "*.data")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("no .data files found", result.stdout)


class TestResolveCheckProjectPaths(unittest.TestCase):
    """Unit tests for `batch-check-projects` path resolution: `--only` role
    filtering, `--subdir` appending, and skip-missing — without a full
    schema build (no $TRUST_ROOT needed)."""

    def _resolved(self, trust_root, primaries, projects):
        from trustify.cli import _CliResolution

        return _CliResolution(projects=projects, trust_root=trust_root, primaries=primaries)

    def _call(self, resolved, *, only=None, subdir="tests/Reference"):
        import argparse
        from types import SimpleNamespace

        from trustify.cli import _resolve_check_project_paths

        ns = SimpleNamespace(only=only, subdir=subdir)
        return {p.resolve() for p in _resolve_check_project_paths(ns, argparse.ArgumentParser(), resolved)}

    def _trip(self, root):
        """Make t/p/dep dirs each with tests/Reference; return (T, P, D, resolved)."""
        t, pp, dd = root / "t", root / "p", root / "dep"
        for x in (t, pp, dd):
            (x / "tests" / "Reference").mkdir(parents=True)
        resolved = self._resolved(str(t), [str(pp)], [str(pp.resolve()), str(dd.resolve())])
        return t, pp, dd, resolved

    def test_no_only_returns_all_with_subdir(self):
        with tempfile.TemporaryDirectory() as d:
            t, pp, dd, resolved = self._trip(Path(d))
            self.assertEqual(
                self._call(resolved),
                {
                    (t / "tests/Reference").resolve(),
                    (pp / "tests/Reference").resolve(),
                    (dd / "tests/Reference").resolve(),
                },
            )

    def test_only_dependency_filters(self):
        with tempfile.TemporaryDirectory() as d:
            _t, _p, dd, resolved = self._trip(Path(d))
            self.assertEqual(self._call(resolved, only="dependency"), {(dd / "tests/Reference").resolve()})

    def test_only_multiple_roles(self):
        with tempfile.TemporaryDirectory() as d:
            t, pp, _d, resolved = self._trip(Path(d))
            self.assertEqual(
                self._call(resolved, only="trust_root,primary"),
                {(t / "tests/Reference").resolve(), (pp / "tests/Reference").resolve()},
            )

    def test_subdir_override(self):
        with tempfile.TemporaryDirectory() as d:
            dd = Path(d) / "dep"
            (dd / "tests" / "Validation").mkdir(parents=True)
            resolved = self._resolved(None, [], [str(dd.resolve())])
            self.assertEqual(
                self._call(resolved, only="dependency", subdir="tests/Validation"),
                {(dd / "tests/Validation").resolve()},
            )

    def test_missing_subdir_is_skipped_with_note(self):
        import contextlib
        import io

        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            pp, dd = root / "p", root / "dep"
            (pp / "tests" / "Reference").mkdir(parents=True)
            dd.mkdir()  # dependency has no tests/Reference
            resolved = self._resolved(None, [str(pp)], [str(pp.resolve()), str(dd.resolve())])
            buf = io.StringIO()
            with contextlib.redirect_stderr(buf):
                got = self._call(resolved)
            self.assertEqual(got, {(pp / "tests/Reference").resolve()})
            self.assertIn("skipping", buf.getvalue())
            self.assertIn(str(dd.resolve()), buf.getvalue())

    def test_unresolvable_workspace_errors(self):
        with self.assertRaises(SystemExit):
            self._call(self._resolved(None, None, None))


class TestBatchCheckProjects(unittest.TestCase):
    @unittest.skipUnless("TRUST_ROOT" in os.environ, "needs TRUST_ROOT for schema generation")
    def test_only_dependency_checks_dep_datasets_not_primary(self):
        with mock_workspace() as ws:
            dep = ws.baltik("mydep")
            main = ws.baltik("mymain", deps={"mydep": dep})
            (dep / "tests" / "Reference").mkdir(parents=True)
            # Deliberately invalid so it surfaces as FAILED in the summary
            # (present in output regardless of verbosity), proving the dir
            # was selected and checked.
            (dep / "tests" / "Reference" / "dep_ds.data").write_text("pas_un_vrai_mot_clef\n")
            (main / "tests" / "Reference").mkdir(parents=True)
            (main / "tests" / "Reference" / "main_ds.data").write_text("pas_un_vrai_mot_clef\n")
            result = _run(
                "--trust-root",
                os.environ["TRUST_ROOT"],
                "--projects",
                str(main),
                "batch-check-projects",
                "--only",
                "dependency",
            )
            combined = result.stdout + result.stderr
            self.assertIn("dep_ds.data", combined, f"stdout={result.stdout}\nstderr={result.stderr}")
            self.assertNotIn("main_ds.data", combined)
            # The selected directories are listed up front, on stdout.
            self.assertIn("checking 1 directory:", result.stdout)
            self.assertIn(str((dep / "tests" / "Reference").resolve()), result.stdout)


def _env_without_trust_root():
    """Subprocess env with `$TRUST_ROOT` and `$project_directory`
    popped. Modernize tests run against a mock baltik only; without
    this guard the resolver's env fallback pulls in the entire TRUST
    source tree (slow and irrelevant to the mock fixture)."""
    env = os.environ.copy()
    env.pop("TRUST_ROOT", None)
    env.pop("project_directory", None)
    return env


class TestCliModernize(unittest.TestCase):
    """End-to-end tests for `trustify modernize`."""

    def test_modernize_dry_run_clean_workspace_exits_zero(self):
        with mock_workspace() as ws:
            ws.baltik(
                "clean_baltik",
                xd_blocks=["// XD clean_block objet_u clean_block BRACE description\n"],
            )
            result = _run(
                "--projects",
                str(ws._baltiks["clean_baltik"]),
                "modernize",
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertEqual(result.stdout.strip(), "")
            self.assertIn("0 files would change", result.stderr)

    def test_modernize_dry_run_pending_changes_exits_one(self):
        with mock_workspace() as ws:
            ws.baltik(
                "legacy_baltik",
                xd_blocks=["// XD legacy_block objet_u legacy_block 0 description\n"],
            )
            result = _run(
                "--projects",
                str(ws._baltiks["legacy_baltik"]),
                "modernize",
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 1, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            # Unified diff signature on stdout.
            self.assertIn("--- ", result.stdout)
            self.assertIn("+++ ", result.stdout)
            self.assertIn("NO_BRACE", result.stdout)
            self.assertIn("brace", result.stderr.lower())

    def test_modernize_apply_writes_files_exits_zero(self):
        with mock_workspace() as ws:
            ws.baltik(
                "apply_baltik",
                xd_blocks=["// XD apply_block objet_u apply_block 1 description\n"],
            )
            result = _run(
                "--projects",
                str(ws._baltiks["apply_baltik"]),
                "modernize",
                "--apply",
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("modernized", result.stderr.lower())
            xd_files = list(Path(ws._baltiks["apply_baltik"]).rglob("*.xd"))
            self.assertTrue(xd_files, "fixture should have produced at least one .xd file")
            content = xd_files[0].read_text()
            self.assertIn("BRACE", content)
            self.assertNotIn(" 1 description", content)

    def test_modernize_apply_on_clean_workspace_is_noop(self):
        with mock_workspace() as ws:
            ws.baltik(
                "noop_baltik",
                xd_blocks=["// XD noop_block objet_u noop_block BRACE description\n"],
            )
            result = _run(
                "--projects",
                str(ws._baltiks["noop_baltik"]),
                "modernize",
                "--apply",
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0)
            self.assertIn("0 files modernized", result.stderr)


class TestCliModernizeScope(unittest.TestCase):
    """Scope tests for `trustify modernize`.

    Modernize is destructive (in-place source rewrites), so it does NOT
    use the four-signal resolution. The rules live in
    `cli/__init__.py::_scope_modernize`:

      - explicit --projects: those + warning, no transitive deps, no trust_root
      - else auto-detected baltik (project.cfg walk-up / $project_directory)
      - else only --trust-root / $TRUST_ROOT: modernize TRUST
      - .trustify.json contributing -> error
      - nothing resolves -> error
    """

    _LEGACY_BLOCK = "// XD legacy_block objet_u legacy_block 1 description\n"

    def test_default_scope_uses_cwd_baltik(self):
        """No --projects + CWD inside a baltik → modernize only that baltik."""
        with mock_workspace() as ws:
            target = ws.baltik("target", xd_blocks=[self._LEGACY_BLOCK])
            other = ws.baltik("other", xd_blocks=[self._LEGACY_BLOCK])
            result = _run(
                "modernize",
                "--apply",
                cwd=str(target),
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (target / "src" / "markers_0.xd").read_text())
            # The other baltik must not have been modernized.
            self.assertIn(" 1 description", (other / "src" / "markers_0.xd").read_text())
            self.assertNotIn("warning", result.stderr.lower())

    def test_default_scope_uses_project_directory_env(self):
        """No --projects + CWD outside any baltik + $project_directory set →
        modernize only that baltik."""
        with mock_workspace() as ws:
            target = ws.baltik("env_target", xd_blocks=[self._LEGACY_BLOCK])
            other = ws.baltik("env_other", xd_blocks=[self._LEGACY_BLOCK])
            env = _env_without_trust_root()
            env["project_directory"] = str(target)
            result = _run(
                "modernize",
                "--apply",
                cwd=str(ws.root),  # outside both baltiks (no project.cfg walk-up)
                env=env,
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (target / "src" / "markers_0.xd").read_text())
            self.assertIn(" 1 description", (other / "src" / "markers_0.xd").read_text())

    def test_explicit_projects_emits_warning(self):
        """--projects A B → both modernized + warning to stderr (apply path).

        Uses distinct .xd basenames in each baltik so the scanner's
        last-project-wins basename dedup (intended for the baltik-overlay-of-
        trust use case) doesn't collapse both fixtures into one.
        """
        with mock_workspace() as ws:
            a = ws.baltik("a")
            (a / "src" / "a_markers.xd").write_text(self._LEGACY_BLOCK)
            b = ws.baltik("b")
            (b / "src" / "b_markers.xd").write_text(self._LEGACY_BLOCK)
            result = _run(
                "--projects",
                str(a),
                "--projects",
                str(b),
                "modernize",
                "--apply",
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (a / "src" / "a_markers.xd").read_text())
            self.assertIn("BRACE", (b / "src" / "b_markers.xd").read_text())
            self.assertIn("warning:", result.stderr.lower())
            self.assertIn("--projects", result.stderr)
            self.assertIn(str(a), result.stderr)
            self.assertIn(str(b), result.stderr)

    def test_explicit_projects_warning_shown_in_dry_run(self):
        """Warning is independent of --apply: dry-run also emits it."""
        with mock_workspace() as ws:
            a = ws.baltik("a", xd_blocks=[self._LEGACY_BLOCK])
            result = _run(
                "--projects",
                str(a),
                "modernize",  # no --apply
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 1)  # pending changes
            self.assertIn("warning:", result.stderr.lower())
            self.assertIn(str(a), result.stderr)

    def test_explicit_projects_skips_transitive_deps(self):
        """--projects A modernizes A even if A's project.cfg declares deps —
        modernize never expands [dependencies]."""
        with mock_workspace() as ws:
            dep = ws.baltik("dep", xd_blocks=[self._LEGACY_BLOCK])
            main = ws.baltik("main", deps={"dep": dep}, xd_blocks=[self._LEGACY_BLOCK])
            result = _run(
                "--projects",
                str(main),
                "modernize",
                "--apply",
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (main / "src" / "markers_0.xd").read_text())
            # The transitive dep was NOT walked / rewritten.
            self.assertIn(" 1 description", (dep / "src" / "markers_0.xd").read_text())

    def test_explicit_projects_ignores_trust_root(self):
        """--projects A --trust-root T → A modernized, T untouched."""
        with mock_workspace() as ws:
            a = ws.baltik("a", xd_blocks=[self._LEGACY_BLOCK])
            trust = ws.root / "trust"
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(self._LEGACY_BLOCK)
            result = _run(
                "--projects",
                str(a),
                "--trust-root",
                str(trust),
                "modernize",
                "--apply",
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (a / "src" / "markers_0.xd").read_text())
            self.assertIn(" 1 description", (trust / "src" / "markers.xd").read_text())

    def test_trust_root_only_modernizes_trust(self):
        """--trust-root T alone (no baltik, no $project_directory) → T modernized."""
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(self._LEGACY_BLOCK)
            result = _run(
                "--trust-root",
                str(trust),
                "modernize",
                "--apply",
                cwd=str(ws.root),  # not inside trust, no baltik
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (trust / "src" / "markers.xd").read_text())
            self.assertNotIn("warning", result.stderr.lower())

    def test_trust_root_env_only_modernizes_trust(self):
        """$TRUST_ROOT alone (no --trust-root, no baltik) also routes to trust."""
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(self._LEGACY_BLOCK)
            env = _env_without_trust_root()
            env["TRUST_ROOT"] = str(trust)
            result = _run(
                "modernize",
                "--apply",
                cwd=str(ws.root),
                env=env,
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (trust / "src" / "markers.xd").read_text())

    def test_baltik_wins_over_trust_root(self):
        """CWD inside baltik + --trust-root T → modernize the baltik, not T."""
        with mock_workspace() as ws:
            baltik = ws.baltik("b", xd_blocks=[self._LEGACY_BLOCK])
            trust = ws.root / "trust"
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(self._LEGACY_BLOCK)
            result = _run(
                "--trust-root",
                str(trust),
                "modernize",
                "--apply",
                cwd=str(baltik),
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("BRACE", (baltik / "src" / "markers_0.xd").read_text())
            self.assertIn(" 1 description", (trust / "src" / "markers.xd").read_text())

    def test_inside_trust_root_warns_but_modernizes_project_directory_baltik(self):
        """CWD inside $TRUST_ROOT + a `$project_directory` overlay no
        longer routes to trust: the env-named baltik is modernized (with
        a WARNING), leaving the TRUST tree untouched. The choice is left
        to the user, matching the schema resolver."""
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            sub = trust / "tests" / "foo"
            sub.mkdir(parents=True)
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(self._LEGACY_BLOCK)
            stale_baltik = ws.baltik("stale", xd_blocks=[self._LEGACY_BLOCK])
            env = _env_without_trust_root()
            env["TRUST_ROOT"] = str(trust)
            env["project_directory"] = str(stale_baltik)
            result = _run(
                "modernize",
                "--apply",
                cwd=str(sub),
                env=env,
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("WARNING", result.stderr)
            # The env-named baltik is modernized...
            self.assertIn("BRACE", (stale_baltik / "src" / "markers_0.xd").read_text())
            # ...and the surrounding TRUST tree is left untouched.
            self.assertIn(" 1 description", (trust / "src" / "markers.xd").read_text())

    def test_baltik_inside_trust_root_is_modernized_with_warning(self):
        """A baltik whose sources live inside $TRUST_ROOT (e.g. ICoCo)
        can be modernized: an auto-detected project.cfg under the TRUST
        tree targets that baltik (with a WARNING) instead of being
        suppressed in favour of the whole TRUST tree."""
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(self._LEGACY_BLOCK)
            inner = ws.baltik("icoco_like", parent=trust, xd_blocks=[self._LEGACY_BLOCK])
            env = _env_without_trust_root()
            env["TRUST_ROOT"] = str(trust)
            result = _run(
                "modernize",
                "--apply",
                cwd=str(inner),
                env=env,
            )
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("WARNING", result.stderr)
            # The inside-TRUST baltik is modernized...
            self.assertIn("BRACE", (inner / "src" / "markers_0.xd").read_text())
            # ...and the surrounding TRUST tree is left untouched.
            self.assertIn(" 1 description", (trust / "src" / "markers.xd").read_text())

    def test_trustify_json_workspace_projects_errors(self):
        """.trustify.json with workspace.projects → modernize refuses."""
        with mock_workspace() as ws:
            a = ws.baltik("a", xd_blocks=[self._LEGACY_BLOCK])
            # Drop .trustify.json next to (but NOT inside) the baltik so the
            # _check_no_project_overlap_with_baltik validation passes.
            cfg_dir = ws.root / "cfg_root"
            cfg_dir.mkdir()
            (cfg_dir / ".trustify.json").write_text(json.dumps({"workspace": {"projects": [str(a)]}}))
            result = _run(
                "modernize",
                cwd=str(cfg_dir),
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 2, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn(".trustify.json", result.stderr)
            self.assertIn("modernize", result.stderr)
            # Source still untouched.
            self.assertIn(" 1 description", (a / "src" / "markers_0.xd").read_text())

    def test_trustify_json_workspace_trust_root_errors(self):
        """.trustify.json with only workspace.trust_root also errors."""
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(self._LEGACY_BLOCK)
            cfg_dir = ws.root / "cfg_root"
            cfg_dir.mkdir()
            (cfg_dir / ".trustify.json").write_text(json.dumps({"workspace": {"trust_root": str(trust)}}))
            result = _run(
                "modernize",
                cwd=str(cfg_dir),
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 2, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn(".trustify.json", result.stderr)

    def test_nothing_resolved_errors(self):
        """No --projects, no baltik, no $project_directory, no trust_root → error."""
        with mock_workspace() as ws:
            empty = ws.root / "empty"
            empty.mkdir()
            result = _run(
                "modernize",
                cwd=str(empty),
                env=_env_without_trust_root(),
            )
            self.assertEqual(result.returncode, 2, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("modernize", result.stderr)
            self.assertIn("nothing resolved", result.stderr.lower())


class TestCliInitConfig(unittest.TestCase):
    def test_init_config_freezes_resolved_workspace(self):
        """init-config freezes the currently-resolved workspace into a
        .trustify.json, going through the same four-signal resolution
        + dependency expansion as the schema-consuming commands. The
        baltik must therefore be a real one with `project.cfg` on disk
        — the point of init-config is to record a *working* config,
        not an aspirational one.
        """
        with mock_workspace() as ws:
            baltik = ws.baltik("baltik_a")
            trust = ws.root / "trust"
            trust.mkdir()  # trust_root is exempt from project.cfg validation
            out_dir = ws.root / "out"
            out_dir.mkdir()
            result = _run("--projects", str(baltik), "--trust-root", str(trust), "init-config", str(out_dir))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            data = json.loads((out_dir / ".trustify.json").read_text())
            self.assertEqual(data["workspace"]["projects"], [str(baltik)])
            self.assertEqual(data["workspace"]["trust_root"], str(trust))
            self.assertTrue(data["workspace"]["auto_resolve_dependencies"])
            self.assertEqual(data["lsp"], {"enum_dedup_threshold": 4})

    def test_init_config_no_args_writes_full_scaffold(self):
        # Isolate from any .trustify.json discoverable upward from the
        # test cwd: run the subprocess inside the tempdir and clear
        # TRUST_ROOT so the resolver gets nothing.
        with tempfile.TemporaryDirectory() as d:
            env = os.environ.copy()
            env.pop("TRUST_ROOT", None)
            result = _run("init-config", d, cwd=d, env=env)
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            data = json.loads((Path(d) / ".trustify.json").read_text())
            self.assertEqual(data["workspace"]["trust_root"], "")
            self.assertEqual(data["workspace"]["projects"], [])
            self.assertTrue(data["workspace"]["auto_resolve_dependencies"])
            self.assertEqual(data["lsp"], {"enum_dedup_threshold": 4})


class TestCliFormat(unittest.TestCase):
    """Tests for the trustify format subcommand. Pure text manipulation
    — no schema needed."""

    SAMPLE = "dimension 2\n{\nfoo 1\n}\n"
    EXPECTED = "dimension 2\n{\n    foo 1\n}\n"

    def test_format_to_stdout(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "in.data"
            f.write_text(self.SAMPLE)
            result = _run("format", str(f))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertEqual(result.stdout, self.EXPECTED)

    def test_format_with_out_writes_file(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "in.data"
            out = Path(d) / "out.data"
            f.write_text(self.SAMPLE)
            result = _run("format", str(f), "--out", str(out))
            self.assertEqual(result.returncode, 0)
            self.assertEqual(out.read_text(), self.EXPECTED)

    def test_format_in_place(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "in.data"
            f.write_text(self.SAMPLE)
            result = _run("format", str(f), "--in-place")
            self.assertEqual(result.returncode, 0)
            self.assertEqual(f.read_text(), self.EXPECTED)

    def test_format_out_and_in_place_are_mutex(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            out = Path(d) / "out.data"
            result = _run("format", str(f), "--out", str(out), "--in-place")
            self.assertEqual(result.returncode, 2)
            self.assertIn("not allowed with argument", result.stderr)

    def test_format_missing_file_errors_fast(self):
        with tempfile.TemporaryDirectory() as d:
            missing = Path(d) / "nope.data"
            result = _run("format", str(missing))
            self.assertEqual(result.returncode, 2)
            self.assertIn("format: no such file", result.stderr)
            self.assertIn(str(missing), result.stderr)

    def test_format_rejects_multiple_files(self):
        with tempfile.TemporaryDirectory() as d:
            f1 = Path(d) / "a.data"
            f1.write_text(self.SAMPLE)
            f2 = Path(d) / "b.data"
            f2.write_text(self.SAMPLE)
            result = _run("format", str(f1), str(f2))
            self.assertEqual(result.returncode, 2)
            self.assertIn("unrecognized arguments", result.stderr)


class TestCliBatchFormat(unittest.TestCase):
    """Tests for the trustify batch-format subcommand.

    Audit 1.10: dry-run is the default; `--apply` (or
    `TRUSTIFY_FORCE_FORMAT` set to any non-empty value, matching the
    NO_COLOR convention used elsewhere) flips to in-place rewrite.
    Mirrors `modernize`'s pattern — destructive is opt-in."""

    SAMPLE = "dimension 2\n{\nfoo 1\n}\n"
    EXPECTED = "dimension 2\n{\n    foo 1\n}\n"

    def test_batch_format_default_is_dry_run(self):
        # Without --apply, the file must be UNCHANGED on disk; the
        # CLI must report a non-zero exit so a CI integration catches
        # unformatted .data files without modifying anything.
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            result = _run("batch-format", str(f))
            self.assertEqual(result.returncode, 1, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertEqual(f.read_text(), self.SAMPLE, "default mode must not modify the file")
            # And the diff is on stdout.
            self.assertIn("foo 1", result.stdout)

    def test_batch_format_default_returns_zero_when_no_changes_needed(self):
        # If every file is already canonical, dry-run still returns 0
        # (no work to do, no diffs to surface).
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.EXPECTED)
            result = _run("batch-format", str(f))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertEqual(f.read_text(), self.EXPECTED)

    def test_batch_format_apply_writes_in_place(self):
        with tempfile.TemporaryDirectory() as d:
            files = []
            for name in ("a.data", "b.data"):
                p = Path(d) / name
                p.write_text(self.SAMPLE)
                files.append(p)
            result = _run("batch-format", "--apply", *[str(p) for p in files])
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            for p in files:
                self.assertEqual(p.read_text(), self.EXPECTED)

    def test_batch_format_apply_single_file(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            result = _run("batch-format", "--apply", str(f))
            self.assertEqual(result.returncode, 0)
            self.assertEqual(f.read_text(), self.EXPECTED)

    def test_batch_format_TRUSTIFY_FORCE_FORMAT_env_enables_apply(self):
        # The baltik Makefile route: setting the env var should have
        # the same effect as passing --apply.
        import os

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            env = os.environ.copy()
            env["TRUSTIFY_FORCE_FORMAT"] = "1"
            result = _run("batch-format", str(f), env=env)
            self.assertEqual(result.returncode, 0)
            self.assertEqual(f.read_text(), self.EXPECTED)

    def test_batch_format_TRUSTIFY_FORCE_FORMAT_empty_does_not_apply(self):
        # NO_COLOR convention: empty value is treated as unset.
        import os

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            env = os.environ.copy()
            env["TRUSTIFY_FORCE_FORMAT"] = ""
            result = _run("batch-format", str(f), env=env)
            self.assertEqual(result.returncode, 1)
            self.assertEqual(f.read_text(), self.SAMPLE)

    def test_batch_format_missing_path_errors_fast(self):
        with tempfile.TemporaryDirectory() as d:
            good = Path(d) / "good.data"
            good.write_text(self.SAMPLE)
            missing = Path(d) / "nope.data"
            result = _run("batch-format", "--apply", str(good), str(missing))
            self.assertEqual(result.returncode, 2)
            self.assertIn("batch-format: no such file or directory", result.stderr)
            self.assertIn(str(missing), result.stderr)
            self.assertEqual(good.read_text(), self.SAMPLE, "good file must be untouched when peer is missing")

    def test_batch_format_rejects_out_flag(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            out = Path(d) / "out.data"
            result = _run("batch-format", "--apply", str(f), "--out", str(out))
            self.assertEqual(result.returncode, 2)
            self.assertIn("unrecognized arguments", result.stderr)

    def test_batch_format_rejects_in_place_flag(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            result = _run("batch-format", "--apply", str(f), "--in-place")
            self.assertEqual(result.returncode, 2)
            self.assertIn("unrecognized arguments", result.stderr)

    # ---- Directory recursion ----

    def _make_tree(self, root: Path) -> None:
        """Build:
        root/a.data
        root/b.data
        root/sub/c.data
        root/sub/d.data
        root/build/skip.data
        root/build/nested/deep.data
        root/sub/legacy/old.data
        root/notes.txt        (non-.data, must be ignored)
        """
        (root / "a.data").write_text(self.SAMPLE)
        (root / "b.data").write_text(self.SAMPLE)
        (root / "notes.txt").write_text("ignored")
        (root / "sub").mkdir()
        (root / "sub" / "c.data").write_text(self.SAMPLE)
        (root / "sub" / "d.data").write_text(self.SAMPLE)
        (root / "build").mkdir()
        (root / "build" / "skip.data").write_text(self.SAMPLE)
        (root / "build" / "nested").mkdir()
        (root / "build" / "nested" / "deep.data").write_text(self.SAMPLE)
        (root / "sub" / "legacy").mkdir()
        (root / "sub" / "legacy" / "old.data").write_text(self.SAMPLE)

    def _all_data_files(self, root: Path):
        return sorted(root.rglob("*.data"))

    def test_batch_format_directory_recurses(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self._make_tree(root)
            result = _run("batch-format", "--apply", str(root))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            for f in self._all_data_files(root):
                self.assertEqual(f.read_text(), self.EXPECTED, f"{f} was not reformatted")
            self.assertEqual((root / "notes.txt").read_text(), "ignored", "non-.data file must be untouched")

    def test_batch_format_mixed_file_and_dir(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self._make_tree(root)
            extra = Path(d) / "extra.data"
            extra.write_text(self.SAMPLE)
            result = _run("batch-format", "--apply", str(root / "sub"), str(extra))
            self.assertEqual(result.returncode, 0)
            self.assertEqual(extra.read_text(), self.EXPECTED)
            self.assertEqual((root / "sub" / "c.data").read_text(), self.EXPECTED)
            # File outside the listed sub-dir is left alone.
            self.assertEqual((root / "a.data").read_text(), self.SAMPLE)

    def test_batch_format_empty_directory_short_circuits(self):
        with tempfile.TemporaryDirectory() as d:
            result = _run("batch-format", d)
            self.assertEqual(result.returncode, 0)
            self.assertIn("no .data files found", result.stdout)

    def test_batch_format_skips_non_data_extensions(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "keep.data").write_text(self.SAMPLE)
            (root / "skip.txt").write_text("ignored")
            (root / "skip.dataset").write_text("ignored")
            result = _run("batch-format", "--apply", str(root))
            self.assertEqual(result.returncode, 0)
            self.assertEqual((root / "keep.data").read_text(), self.EXPECTED)
            self.assertEqual((root / "skip.txt").read_text(), "ignored")
            self.assertEqual((root / "skip.dataset").read_text(), "ignored")

    # ---- Exclude semantics ----

    def test_batch_format_exclude_by_file_basename(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self._make_tree(root)
            result = _run("batch-format", "--apply", str(root), "--exclude", "old.data")
            self.assertEqual(result.returncode, 0)
            # old.data anywhere in the tree is left alone ...
            self.assertEqual((root / "sub" / "legacy" / "old.data").read_text(), self.SAMPLE)
            # ... everything else is reformatted.
            self.assertEqual((root / "a.data").read_text(), self.EXPECTED)
            self.assertEqual((root / "sub" / "c.data").read_text(), self.EXPECTED)

    def test_batch_format_exclude_by_dir_basename_prunes(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self._make_tree(root)
            result = _run("batch-format", "--apply", str(root), "--exclude", "build")
            self.assertEqual(result.returncode, 0)
            # Every file under build/ — at any depth — is untouched.
            self.assertEqual((root / "build" / "skip.data").read_text(), self.SAMPLE)
            self.assertEqual((root / "build" / "nested" / "deep.data").read_text(), self.SAMPLE)
            # Siblings still get formatted.
            self.assertEqual((root / "a.data").read_text(), self.EXPECTED)

    def test_batch_format_exclude_by_relative_path(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self._make_tree(root)
            # Anchored path pattern — only the listed root's sub/legacy
            # matches, a hypothetical other/legacy wouldn't.
            result = _run("batch-format", "--apply", str(root), "--exclude", "sub/legacy")
            self.assertEqual(result.returncode, 0)
            self.assertEqual((root / "sub" / "legacy" / "old.data").read_text(), self.SAMPLE)
            self.assertEqual((root / "sub" / "c.data").read_text(), self.EXPECTED)

    def test_batch_format_exclude_trailing_slash_is_dir_only(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root / "old.data").write_text(self.SAMPLE)
            (root / "old").mkdir()
            (root / "old" / "inside.data").write_text(self.SAMPLE)
            result = _run("batch-format", "--apply", str(root), "--exclude", "old/")
            self.assertEqual(result.returncode, 0)
            # The directory `old/` is pruned ...
            self.assertEqual((root / "old" / "inside.data").read_text(), self.SAMPLE)
            # ... but a file literally named `old.data` is NOT a directory,
            # so the dir-only pattern doesn't match it.
            self.assertEqual((root / "old.data").read_text(), self.EXPECTED)

    def test_batch_format_explicit_file_bypasses_exclude(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            f = root / "legacy.data"
            f.write_text(self.SAMPLE)
            result = _run("batch-format", "--apply", str(f), "--exclude", "legacy.data")
            self.assertEqual(result.returncode, 0)
            self.assertEqual(f.read_text(), self.EXPECTED, "explicit file must format even if it matches --exclude")

    def test_batch_format_multiple_excludes(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            self._make_tree(root)
            result = _run("batch-format", "--apply", str(root), "--exclude", "build", "--exclude", "old.data")
            self.assertEqual(result.returncode, 0)
            self.assertEqual((root / "build" / "skip.data").read_text(), self.SAMPLE)
            self.assertEqual((root / "sub" / "legacy" / "old.data").read_text(), self.SAMPLE)
            self.assertEqual((root / "a.data").read_text(), self.EXPECTED)
            self.assertEqual((root / "sub" / "c.data").read_text(), self.EXPECTED)

    def test_batch_format_skips_symlinked_files_by_default(self):
        """Audit 2.6: symlinked files at the walk's leaves are skipped
        by default. Without this guard, batch-check / batch-format
        against a directory containing symlinks-to-elsewhere would
        silently include files outside the intended scope."""
        import os as _os

        with tempfile.TemporaryDirectory() as outer:
            outside = Path(outer) / "outside.data"
            outside.write_text(self.SAMPLE)
            with tempfile.TemporaryDirectory() as d:
                root = Path(d)
                (root / "real.data").write_text(self.SAMPLE)
                link = root / "linked.data"
                _os.symlink(outside, link)

                result = _run("batch-format", str(root))
                # Symlink target must not have been formatted.
                self.assertEqual(outside.read_text(), self.SAMPLE)
                # The walk surfaced one diff for the non-symlink file
                # only; exit 1 because that file would change.
                self.assertEqual(result.returncode, 1, f"stderr: {result.stderr}\nstdout: {result.stdout}")
                self.assertIn("linked.data", result.stderr, "expected a skip warning on stderr")
                self.assertIn("symlink", result.stderr.lower())
                # And the real file's diff IS on stdout.
                self.assertIn("real.data", result.stdout)

    def test_batch_format_follow_symlinks_includes_symlinked_files(self):
        """--follow-symlinks opts in: symlinked files at the walk's
        leaves are processed (no skip warning). The atomic write
        resolves the link, so the symlink survives and the TARGET
        gets the formatted content (matches vim / sed
        --follow-symlinks semantics)."""
        import os as _os

        with tempfile.TemporaryDirectory() as outer:
            outside = Path(outer) / "outside.data"
            outside.write_text(self.SAMPLE)
            with tempfile.TemporaryDirectory() as d:
                root = Path(d)
                link = root / "linked.data"
                _os.symlink(outside, link)

                result = _run("batch-format", "--apply", "--follow-symlinks", str(root))
                self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
                # Target was rewritten via the link; link still a link.
                self.assertEqual(outside.read_text(), self.EXPECTED)
                self.assertTrue(link.is_symlink(), "atomic write must not replace the symlink")
                # No skip warning.
                self.assertNotIn("symlink", result.stderr.lower())

    def test_batch_format_explicit_symlink_arg_writes_through_link(self):
        """An explicit symlink arg bypasses the walk-time skip — same
        precedent as the existing `--exclude` bypass for explicit CLI
        files. The atomic write follows the link to the target so the
        symlink stays intact and the underlying file is updated."""
        import os as _os

        with tempfile.TemporaryDirectory() as outer:
            outside = Path(outer) / "outside.data"
            outside.write_text(self.SAMPLE)
            with tempfile.TemporaryDirectory() as d:
                link = Path(d) / "linked.data"
                _os.symlink(outside, link)

                result = _run("batch-format", "--apply", str(link))
                self.assertEqual(result.returncode, 0)
                self.assertTrue(link.is_symlink(), "atomic write must not replace the symlink")
                self.assertEqual(outside.read_text(), self.EXPECTED)

    def test_batch_format_does_not_recurse_symlinked_directory_by_default(self):
        """Symlinked directories are not recursed by default — preserves
        the existing `os.walk(followlinks=False)` behaviour."""
        import os as _os

        with tempfile.TemporaryDirectory() as outer:
            real_subdir = Path(outer) / "real_sub"
            real_subdir.mkdir()
            (real_subdir / "deep.data").write_text(self.SAMPLE)
            with tempfile.TemporaryDirectory() as d:
                root = Path(d)
                _os.symlink(real_subdir, root / "sub_link")
                (root / "top.data").write_text(self.SAMPLE)

                result = _run("batch-format", str(root))
                # The walk only saw top.data; deep.data inside the
                # symlinked subdir was not visited.
                self.assertNotIn("deep.data", result.stdout)
                self.assertEqual((real_subdir / "deep.data").read_text(), self.SAMPLE)

    def test_batch_format_follow_symlinks_recurses_symlinked_dir(self):
        import os as _os

        with tempfile.TemporaryDirectory() as outer:
            real_subdir = Path(outer) / "real_sub"
            real_subdir.mkdir()
            (real_subdir / "deep.data").write_text(self.SAMPLE)
            with tempfile.TemporaryDirectory() as d:
                root = Path(d)
                _os.symlink(real_subdir, root / "sub_link")

                result = _run("batch-format", "--apply", "--follow-symlinks", str(root))
                self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
                # The symlinked dir's deep.data WAS visited and formatted.
                self.assertEqual((real_subdir / "deep.data").read_text(), self.EXPECTED)

    def test_batch_format_exclude_glob_metachars(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            for name in ("foo.data", "foo_backup.data", "bar.data"):
                (root / name).write_text(self.SAMPLE)
            result = _run("batch-format", "--apply", str(root), "--exclude", "*_backup.data")
            self.assertEqual(result.returncode, 0)
            self.assertEqual((root / "foo_backup.data").read_text(), self.SAMPLE)
            self.assertEqual((root / "foo.data").read_text(), self.EXPECTED)
            self.assertEqual((root / "bar.data").read_text(), self.EXPECTED)


class TestCliWorkspaceResolution(unittest.TestCase):
    """Per-path workspace discovery: the CLI walks up from each path
    argument looking for `.trustify.json` / `project.cfg`. Multi-path
    commands require all anchors to agree on both signals; the CLI
    errors out otherwise (chosen over silently picking one).

    These tests drive `batch-check`, which goes through
    `_resolve_project_args` like every other schema-consuming
    subcommand. The conflict / dep-expansion errors fire *during*
    resolution (parser.error → exit 2), well before any schema
    generation, so no real TRUST_ROOT is required.
    """

    SAMPLE = "dimension 2\n{\nfoo 1\n}\n"

    def test_two_paths_under_different_configs_errors_out(self):
        with mock_workspace() as ws:
            ws_a = ws.root / "ws_a"
            ws_b = ws.root / "ws_b"
            ws_a.mkdir()
            ws_b.mkdir()
            (ws_a / ".trustify.json").write_text('{"workspace": {"projects": ["/opt/baltik-A"]}}')
            (ws_b / ".trustify.json").write_text('{"workspace": {"projects": ["/opt/baltik-B"]}}')
            data_a = ws_a / "a.data"
            data_b = ws_b / "b.data"
            data_a.write_text(self.SAMPLE)
            data_b.write_text(self.SAMPLE)
            result = _run("batch-check", str(data_a), str(data_b))
            self.assertEqual(result.returncode, 2, f"stdout={result.stdout}\nstderr={result.stderr}")
            self.assertIn("different .trustify.json files", result.stderr)

    def test_two_paths_under_same_config_do_not_trigger_conflict(self):
        """Two data paths under the same workspace must NOT fire the
        multi-path conflict detector. Downstream schema generation
        will still fail in this test (the fake baltik has no XD
        sources, no real TRUST_ROOT), but the conflict check itself
        — the thing under test — must stay silent.
        """
        with mock_workspace() as ws:
            baltik = ws.baltik("baltik_a")
            (ws.root / ".trustify.json").write_text(f'{{"workspace": {{"projects": ["{baltik}"]}}}}')
            data_a = ws.root / "a.data"
            data_b = ws.root / "b.data"
            data_a.write_text(self.SAMPLE)
            data_b.write_text(self.SAMPLE)
            result = _run("batch-check", str(data_a), str(data_b))
            self.assertNotIn("different .trustify.json", result.stderr)
            self.assertNotIn("different baltik roots", result.stderr)

    def test_dependency_resolution_error_propagates_to_parser_error(self):
        """If a project's [dependencies] section points at a missing
        directory, the CLI should report it cleanly (parser.error,
        exit 2), not crash with a stack trace.
        """
        with mock_workspace() as ws:
            main = ws.baltik("main", deps={"ghost": "./ghost"})
            data = main / "a.data"
            data.write_text(self.SAMPLE)
            result = _run("batch-check", str(data))
            self.assertEqual(result.returncode, 2, f"stdout={result.stdout}\nstderr={result.stderr}")
            self.assertIn("ghost", result.stderr)
            self.assertIn("dependency", result.stderr.lower())

    def test_no_auto_deps_flag_skips_dependency_expansion(self):
        """With --no-auto-deps, a project whose [dependencies] points at
        a missing dir is NOT expanded — and so doesn't surface the dep
        error. Schema generation downstream may still fail for other
        reasons (empty schema vs real dataset), but the specific dep
        error must NOT appear.
        """
        with mock_workspace() as ws:
            main = ws.baltik("main", deps={"ghost": "./ghost"})
            data = main / "a.data"
            data.write_text(self.SAMPLE)
            result = _run("--no-auto-deps", "batch-check", str(data))
            self.assertNotIn("ghost", result.stderr)
            self.assertNotIn("dependency resolution failed", result.stderr.lower())

    def test_two_paths_under_different_baltik_roots_errors_out(self):
        with mock_workspace() as ws:
            ba = ws.baltik("baltik_a")
            bb = ws.baltik("baltik_b")
            data_a = ba / "a.data"
            data_b = bb / "b.data"
            data_a.write_text(self.SAMPLE)
            data_b.write_text(self.SAMPLE)
            result = _run("batch-check", str(data_a), str(data_b))
            self.assertEqual(result.returncode, 2, f"stdout={result.stdout}\nstderr={result.stderr}")
            self.assertIn("different baltik roots", result.stderr)

    def test_inside_trust_root_warns_but_keeps_stale_project_directory(self):
        """A `$project_directory` overlay while the anchor sits inside
        $TRUST_ROOT is no longer suppressed — the overlay is applied (so
        a user can deliberately check their baltik against TRUST
        datasets) but a WARNING is emitted so an accidental stale env is
        visible rather than silently dropped.
        """
        with mock_workspace() as ws, tempfile.TemporaryDirectory() as d:
            stale = ws.baltik("stale_baltik")
            trust_root = Path(d) / "trust"
            (trust_root / "src").mkdir(parents=True)
            env = os.environ.copy()
            env["project_directory"] = str(stale)
            result = _run(
                "--trust-root",
                str(trust_root),
                "projects",
                cwd=str(trust_root),
                env=env,
            )
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            self.assertIn(str(stale.resolve()), result.stdout)
            self.assertIn("WARNING", result.stderr)
            self.assertIn("inside $TRUST_ROOT", result.stderr)

    def test_resolution_summary_logs_implicit_signals(self):
        """When any signal is implicit (auto-detect, env var, transitive
        dep), trustify echoes a one-shot 'resolved workspace' block on
        stderr so users see what tree they actually got.
        """
        # Use a trust_root in a sibling tempdir so the baltik is NOT
        # under it — otherwise the anchor-inside-trust_root suppression
        # rule kicks in and the baltik auto-detection is skipped.
        with mock_workspace() as ws, tempfile.TemporaryDirectory() as trust_dir:
            main = ws.baltik("main")
            (Path(trust_dir) / "src").mkdir()
            data = main / "a.data"
            data.write_text("dimension 2\n")
            env = os.environ.copy()
            env["TRUST_ROOT"] = str(trust_dir)
            env.pop("project_directory", None)
            result = _run("check", str(data), env=env)
            self.assertIn("trustify: resolved workspace:", result.stderr)
            self.assertIn("$TRUST_ROOT env", result.stderr)
            self.assertIn("auto-detected", result.stderr)
            self.assertIn(str(main), result.stderr)

    def test_resolution_summary_silent_when_all_explicit(self):
        """When every signal came from a CLI flag (no env var, no
        auto-detect, no dep expansion), the summary stays quiet.
        """
        with mock_workspace() as ws, tempfile.TemporaryDirectory() as d:
            main = ws.baltik("main")  # no deps
            trust_root = Path(d) / "trust"
            (trust_root / "src").mkdir(parents=True)
            data = Path(d) / "a.data"
            data.write_text("dimension 2\n")
            env = os.environ.copy()
            env.pop("TRUST_ROOT", None)
            env.pop("project_directory", None)
            result = _run(
                "--trust-root",
                str(trust_root),
                "-p",
                str(main),
                "--no-auto-deps",  # otherwise the dep walk would log
                "check",
                str(data),
                env=env,
            )
            self.assertNotIn("resolved workspace", result.stderr)

    def test_warns_when_implicit_overlay_active_inside_trust_root(self):
        """When the anchor sits inside $TRUST_ROOT and an implicit baltik
        signal ($project_directory or an auto-detected project.cfg)
        supplies an overlay, trustify warns (naming the overlaid baltik)
        but proceeds — it no longer silently drops the overlay.
        """
        with mock_workspace() as ws, tempfile.TemporaryDirectory() as d:
            trust_root = Path(d) / "trust"
            (trust_root / "src").mkdir(parents=True)
            stale = ws.baltik("stale_baltik")
            env = os.environ.copy()
            env["project_directory"] = str(stale)
            # Anchor (CWD) inside trust_root → warning fires, overlay kept.
            result = _run(
                "--trust-root",
                str(trust_root),
                "generate_schema",
                "--out",
                str(Path(d) / "schema_out"),
                cwd=str(trust_root),
                env=env,
            )
            self.assertIn("WARNING", result.stderr)
            self.assertIn(str(stale.resolve()), result.stderr)
            self.assertIn("inside $TRUST_ROOT", result.stderr)

    def test_baltik_inside_trust_root_is_overlaid_with_warning(self):
        """A baltik whose sources live inside $TRUST_ROOT (e.g. ICoCo) is
        a legitimate overlay target: an auto-detected project.cfg inside
        the TRUST tree must still be used (so a dev can check whether
        their baltik's XD tags break TRUST datasets), with a WARNING
        rather than a silent drop.
        """
        with mock_workspace() as ws:
            trust_root = ws.root / "trust"
            (trust_root / "src").mkdir(parents=True)
            inner = ws.baltik("icoco_like", parent=trust_root)
            data = inner / "a.data"
            data.write_text("dimension 2\n")
            env = os.environ.copy()
            env["TRUST_ROOT"] = str(trust_root)
            env.pop("project_directory", None)
            result = _run("check", str(data), env=env)
            self.assertIn("WARNING", result.stderr)
            self.assertIn("inside $TRUST_ROOT", result.stderr)
            self.assertIn(str(inner.resolve()), result.stderr)

    def test_resolution_summary_skipped_for_projects_subcommand(self):
        """`projects` itself is a stdout dump of the resolved list;
        re-printing the same data to stderr would be noise.
        """
        with mock_workspace() as ws:
            main = ws.baltik("main")
            env = os.environ.copy()
            env["TRUST_ROOT"] = str(ws.root)
            result = _run("-p", str(main), "projects", env=env)
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            self.assertNotIn("resolved workspace", result.stderr)

    def test_inside_trust_root_keeps_explicit_projects_flag_without_warning(self):
        """`--projects` is an explicit user request: it applies even when
        the anchor sits inside $TRUST_ROOT, and the inside-$TRUST_ROOT
        warning does NOT fire (the user clearly asked for the overlay).
        """
        with mock_workspace() as ws, tempfile.TemporaryDirectory() as d:
            explicit = ws.baltik("explicit")
            trust_root = Path(d) / "trust"
            (trust_root / "src").mkdir(parents=True)
            env = os.environ.copy()
            env.pop("project_directory", None)
            result = _run(
                "--trust-root",
                str(trust_root),
                "-p",
                str(explicit),
                "projects",
                cwd=str(trust_root),
                env=env,
            )
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            self.assertIn(str(explicit.resolve()), result.stdout)
            self.assertNotIn("WARNING", result.stderr)


class TestCliProjects(unittest.TestCase):
    """`trustify projects` dumps the resolved overlay list. Default output
    is one path per line (shell-script-friendly); --json adds structured
    role labels; --only / --exclude filter by role.

    Roles: trust_root (the TRUST source tree), primary (explicitly
    listed — CLI flag / .trustify.json / auto-detect), dependency
    (pulled in transitively from project.cfg).
    """

    def test_lists_resolved_baltiks_in_overlay_order(self):
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            trust.mkdir()
            dep = ws.baltik("dep")
            main = ws.baltik("main", deps={"dep": dep})
            result = _run("--trust-root", str(trust), "-p", str(main), "projects")
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            # Order: trust_root, then deps in post-order, then primaries last.
            lines = [line for line in result.stdout.splitlines() if line]
            self.assertEqual(lines, [str(trust), str(dep), str(main)])

    def test_json_output_carries_roles(self):
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            trust.mkdir()
            dep = ws.baltik("dep")
            main = ws.baltik("main", deps={"dep": dep})
            result = _run("--trust-root", str(trust), "-p", str(main), "projects", "--json")
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            payload = json.loads(result.stdout)
            self.assertEqual(
                payload,
                [
                    {"path": str(trust), "role": "trust_root"},
                    {"path": str(dep), "role": "dependency"},
                    {"path": str(main), "role": "primary"},
                ],
            )

    def test_only_dependency_skips_trust_and_primary(self):
        """The canonical use case from the baltik scaffold Makefile:
        feed the dependency baltik paths to Doxygen without TRUST or
        the current baltik (both handled separately in the Doxyfile).
        """
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            trust.mkdir()
            dep = ws.baltik("dep")
            main = ws.baltik("main", deps={"dep": dep})
            result = _run("--trust-root", str(trust), "-p", str(main), "projects", "--only=dependency")
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            self.assertEqual(result.stdout.strip(), str(dep))

    def test_exclude_filters_out_listed_roles(self):
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            trust.mkdir()
            dep = ws.baltik("dep")
            main = ws.baltik("main", deps={"dep": dep})
            result = _run("--trust-root", str(trust), "-p", str(main), "projects", "--exclude=trust_root,primary")
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            self.assertEqual(result.stdout.strip(), str(dep))

    def test_only_and_exclude_are_mutually_exclusive(self):
        with mock_workspace() as ws:
            main = ws.baltik("main")
            result = _run("-p", str(main), "projects", "--only=primary", "--exclude=trust_root")
            self.assertEqual(result.returncode, 2)
            self.assertIn("mutually exclusive", result.stderr)

    def test_invalid_role_is_rejected(self):
        with mock_workspace() as ws:
            main = ws.baltik("main")
            result = _run("-p", str(main), "projects", "--only=nope")
            self.assertEqual(result.returncode, 2)
            self.assertIn("invalid role", result.stderr.lower())
            self.assertIn("nope", result.stderr)

    def test_errors_when_nothing_resolves(self):
        """No --projects, no --trust-root, no env, no walkup → resolution
        produces an empty list, so the command errors out (matches every
        other schema-needing command). Run inside an isolated tempdir
        with $TRUST_ROOT / $project_directory cleared.
        """
        env = os.environ.copy()
        env.pop("TRUST_ROOT", None)
        env.pop("project_directory", None)
        with tempfile.TemporaryDirectory() as d:
            result = _run("projects", cwd=d, env=env)
            # Sometimes the test machine has a .trustify.json upstream we
            # can't isolate from — skip in that case.
            if result.returncode == 0:
                self.skipTest("environment has an ambient config / detect that resolves projects")
            self.assertEqual(result.returncode, 2, f"stderr={result.stderr}")
            self.assertIn("no projects to resolve", result.stderr.lower())

    def test_empty_filter_result_is_not_an_error(self):
        """Filtering down to a role with no matches outputs nothing and
        exits 0 — only "no projects resolved at all" is an error.
        """
        with mock_workspace() as ws:
            main = ws.baltik("main")  # no deps
            # Only deps exist? There are none. Should print empty, exit 0.
            result = _run("-p", str(main), "projects", "--only=dependency")
            self.assertEqual(result.returncode, 0, f"stderr={result.stderr}")
            self.assertEqual(result.stdout.strip(), "")


class TestCliGlobalFlags(unittest.TestCase):
    """Sanity checks for global-flag placement rules."""

    def test_trust_root_after_subcommand_rejected(self):
        result = _run("check", "foo.data", "--trust-root", "/x")
        self.assertEqual(result.returncode, 2)
        self.assertIn("unrecognized arguments", result.stderr)

    def test_schema_excludes_trust_root(self):
        result = _run("--schema", "/s", "--trust-root", "/t", "generate_schema")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--schema is mutually exclusive", result.stderr)

    def test_schema_excludes_projects(self):
        result = _run("--schema", "/s", "-p", "/p", "generate_schema")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--schema is mutually exclusive", result.stderr)

    def test_schema_rejected_with_non_consumers(self):
        # `--schema` only makes sense for commands that read a generated
        # schema. Reject it on every other subcommand.
        for cmd_args in (
            ("generate_schema",),
            ("format", "x.data"),
            ("batch-format", "x.data"),
            ("init-config",),
        ):
            with self.subTest(cmd=cmd_args[0]):
                result = _run("--schema", "/s", *cmd_args)
                self.assertEqual(result.returncode, 2, f"stderr: {result.stderr}")
                self.assertIn(f"--schema cannot be used with {cmd_args[0]}", result.stderr)


_MINI_TRAD2_CLI = (
    "pb_hydraulique objet_lecture pb_hydraulique BRACE A problem.\n"
    "  attr tinit floattant t_init OPT initial time\n"
    "  attr tmax floattant tmax REQ max time\n"
    "pb_base objet_lecture pb_base BRACE abstract base.\n"
)


class TestCliKeywords(unittest.TestCase):
    """End-to-end test driving the CLI via `--schema` so we don't have to
    run the full schema generation here.
    """

    def _write_schema(self, d: Path) -> Path:
        schema_dir = d / "schema"
        schema_dir.mkdir()
        (schema_dir / "TRAD2_trustify").write_text(_MINI_TRAD2_CLI)
        return schema_dir

    def test_keywords_writes_both_files(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = self._write_schema(tmp)
            out_dir = tmp / "out"
            result = _run("--schema", str(schema_dir), "generate_keywords", "--out", str(out_dir))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertTrue((out_dir / "Keywords.txt").exists())
            self.assertTrue((out_dir / "Keywords.Vim").exists())

    def test_keywords_content_shape(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = self._write_schema(tmp)
            result = _run("--schema", str(schema_dir), "generate_keywords", "--out", str(tmp))
            self.assertEqual(result.returncode, 0)
            txt = (tmp / "Keywords.txt").read_text()
            vim = (tmp / "Keywords.Vim").read_text()
            self.assertIn("|pb_hydraulique\n", txt)
            self.assertIn("|tinit\n", txt)
            self.assertNotIn("|pb_base\n", txt)
            self.assertTrue(vim.startswith("syntax keyword TRUSTLanguageKeywords  "))
            self.assertIn("  pb_hydraulique  ", vim)


class TestCliGenerateMarkdown(unittest.TestCase):
    """CLI surface tests for `trustify generate_markdown` — covers argparse
    validation and `--out` precondition violations.

    These tests do not run the full generation (no `--trust-root` /
    schema is provided to those that don't need it)."""

    def test_generate_markdown_missing_out_arg_errors(self):
        result = _run("generate_markdown")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--out", result.stderr)

    def test_generate_markdown_out_is_file_errors(self):
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "not_a_dir"
            f.write_text("")
            result = _run("--schema", "/nonexistent", "generate_markdown", "--out", str(f))
            self.assertEqual(result.returncode, 2)
            self.assertIn("not a directory", result.stderr)

    def test_generate_markdown_out_is_nonempty_dir_errors(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "stale.md").write_text("")
            result = _run("--schema", "/nonexistent", "generate_markdown", "--out", d)
            self.assertEqual(result.returncode, 2)
            self.assertIn("not empty", result.stderr)

    def test_generate_markdown_writes_keyword_reference(self):
        if "TRUST_ROOT" not in os.environ:
            self.skipTest("TRUST_ROOT must be set")
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            result = _run("--trust-root", os.environ["TRUST_ROOT"], "generate_markdown", "--out", str(out))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertTrue((out / "keyword_reference.md").exists())


def _make_cache_entry(root, entry_id, *, provenance=None):
    """Build a fake `<root>/.cache/trustify/<entry_id>` and return its path.

    Writes all three is_cached completion sentinels: `trustify_gen.py`,
    `trustify_gen_pyd_<entry_id>.py` (digest-suffix convention enforced
    by `api.generate_schema`), and `provenance.json` (empty `{}` when
    `provenance` is None — keeps the "unknown timestamp / version"
    cases in this file working, since `_load_entry` returns None
    fields for an empty provenance dict).
    """
    import json

    cache_dir = Path(root) / ".cache" / "trustify" / entry_id
    cache_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "trustify_gen.py").write_text("# stub gen\n")
    (cache_dir / f"trustify_gen_pyd_{entry_id}.py").write_text("# stub pyd\n")
    (cache_dir / "provenance.json").write_text(
        json.dumps(provenance if provenance is not None else {}, indent=2, sort_keys=True) + "\n"
    )
    return cache_dir


def _run_with_home(home, *args, input_text=None):
    """Invoke the CLI with `HOME` pointing at `home` so the cache root
    is redirected to `<home>/.cache/trustify/`.
    """
    env = os.environ.copy()
    env["HOME"] = str(home)
    cmd = [sys.executable, "-m", "trustify.cli", *args]
    return subprocess.run(cmd, env=env, capture_output=True, text=True, input=input_text)


class TestCliCache(unittest.TestCase):
    """Tests for the `trustify cache` sub-subcommand surface."""

    def test_list_no_entries_prints_message(self):
        with tempfile.TemporaryDirectory() as home:
            result = _run_with_home(home, "cache", "list")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("No cache entries", result.stdout)

    def test_list_renders_table(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(
                home,
                "abc1234567890def",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1.2.3",
                    "projects": {"trust": "/abs/trust", "baltik_a": "/abs/b"},
                    "generation_time_seconds": 12.3,
                    "generation_steps_seconds": {
                        "scan_sources": 11.0,
                        "write_trad2": 0.1,
                        "generate_pyd_and_pars": 1.2,
                    },
                },
            )
            result = _run_with_home(home, "cache", "list")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("abc12345", result.stdout)  # 8-char id
            self.assertIn("2026-05-20 18:00:00", result.stdout)
            self.assertIn("1.2.3", result.stdout)
            self.assertIn("12.3s", result.stdout)  # total time
            # Table shows project names only, comma-separated. Absolute
            # paths and per-step timings are reserved for `cache list --json`.
            self.assertIn("trust", result.stdout)
            self.assertIn("baltik_a", result.stdout)
            self.assertNotIn("/abs/trust", result.stdout)
            self.assertNotIn("scan_sources", result.stdout)

    def test_list_renders_unknown_time_for_old_entries(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc1234567890def", provenance=None)
            result = _run_with_home(home, "cache", "list")
            self.assertEqual(result.returncode, 0)
            # TIME column should show (unknown) when no timings are stored.
            self.assertIn("(unknown)", result.stdout)

    def test_list_renders_minutes_for_long_generations(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(
                home,
                "abc1234567890def",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {},
                    "generation_time_seconds": 75.0,
                },
            )
            result = _run_with_home(home, "cache", "list")
            self.assertEqual(result.returncode, 0)
            self.assertIn("1m15s", result.stdout)

    def test_list_marks_missing_provenance_as_unknown(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc1234567890def", provenance=None)
            result = _run_with_home(home, "cache", "list")
            self.assertEqual(result.returncode, 0)
            self.assertIn("(unknown)", result.stdout)

    def test_list_json_emits_machine_readable(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(
                home,
                "abc1234567890def",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1.2.3",
                    "projects": {"trust": "/abs/trust"},
                    "generation_time_seconds": 5.5,
                    "generation_steps_seconds": {"scan_sources": 4.0, "write_trad2": 0.1, "generate_pyd_and_pars": 1.4},
                },
            )
            result = _run_with_home(home, "cache", "list", "--json")
            self.assertEqual(result.returncode, 0)
            data = json.loads(result.stdout)
            self.assertEqual(len(data), 1)
            self.assertEqual(data[0]["id"], "abc1234567890def")
            self.assertEqual(data[0]["trustify_version"], "1.2.3")
            self.assertEqual(data[0]["generation_time_seconds"], 5.5)
            self.assertEqual(
                data[0]["generation_steps_seconds"],
                {
                    "scan_sources": 4.0,
                    "write_trad2": 0.1,
                    "generate_pyd_and_pars": 1.4,
                },
            )

    def test_path_prints_absolute_dir(self):
        with tempfile.TemporaryDirectory() as home:
            entry = _make_cache_entry(home, "abc1234567890def")
            result = _run_with_home(home, "cache", "path", "abc")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertEqual(result.stdout.strip(), str(entry.resolve()))

    def test_path_ambiguous_prefix_errors(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc11111")
            _make_cache_entry(home, "abc22222")
            result = _run_with_home(home, "cache", "path", "abc")
            self.assertEqual(result.returncode, 2)
            self.assertIn("ambiguous", result.stderr)

    def test_path_unknown_prefix_errors(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc11111")
            result = _run_with_home(home, "cache", "path", "ZZ")
            self.assertEqual(result.returncode, 2)
            self.assertIn("no cache entry matches", result.stderr)

    def test_copy_creates_dest_with_filtered_contents(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc1234567890def")
            with tempfile.TemporaryDirectory() as dst_parent:
                dest = Path(dst_parent) / "out"
                result = _run_with_home(
                    home,
                    "cache",
                    "copy",
                    "abc",
                    str(dest),
                )
                self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
                self.assertTrue((dest / "trustify_gen.py").is_file())
                self.assertEqual(result.stdout.strip(), str(dest))

    def test_copy_refuses_nonempty_dest(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc1234567890def")
            with tempfile.TemporaryDirectory() as dst_parent:
                dest = Path(dst_parent) / "out"
                dest.mkdir()
                (dest / "blocker").write_text("x")
                result = _run_with_home(
                    home,
                    "cache",
                    "copy",
                    "abc",
                    str(dest),
                )
                self.assertEqual(result.returncode, 2)
                self.assertIn("not empty", result.stderr)

    def test_clean_explicit_ids_skips_prompt(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc11111")
            _make_cache_entry(home, "def22222")
            result = _run_with_home(home, "cache", "clean", "abc")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("abc11111", result.stdout)
            self.assertFalse((Path(home) / ".cache" / "trustify" / "abc11111").exists())
            self.assertTrue((Path(home) / ".cache" / "trustify" / "def22222").exists())

    def test_clean_all_with_force_skips_prompt(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc11111")
            _make_cache_entry(home, "def22222")
            result = _run_with_home(home, "cache", "clean", "--force")
            self.assertEqual(result.returncode, 0)
            self.assertFalse((Path(home) / ".cache" / "trustify" / "abc11111").exists())
            self.assertFalse((Path(home) / ".cache" / "trustify" / "def22222").exists())

    def test_clean_all_prompts_and_aborts_on_no(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc11111")
            result = _run_with_home(home, "cache", "clean", input_text="n\n")
            # Returncode = 1 on user-declined; entry untouched.
            self.assertEqual(result.returncode, 1)
            self.assertIn("Aborted", result.stdout)
            self.assertTrue((Path(home) / ".cache" / "trustify" / "abc11111").exists())

    def test_clean_all_prompts_and_proceeds_on_yes(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc11111")
            result = _run_with_home(home, "cache", "clean", input_text="y\n")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertFalse((Path(home) / ".cache" / "trustify" / "abc11111").exists())

    def test_clean_no_entries_short_circuits(self):
        with tempfile.TemporaryDirectory() as home:
            # No cache populated.
            result = _run_with_home(home, "cache", "clean", "--force")
            self.assertEqual(result.returncode, 0)
            self.assertIn("No cache entries", result.stdout)

    def test_clean_aborts_cleanly_on_eof(self):
        """Audit 5.2: `cache clean` (without --force, with an entry to
        delete) prompts via `input(...)`. Closing stdin used to raise
        `EOFError` straight through the top-level except, producing the
        confusing one-liner `trustify: EOFError: EOF when reading a
        line`. It should be treated like a user-declined prompt:
        print 'Aborted.' and exit 1, leaving the entry untouched.
        Common trigger: scripted invocations that pipe nothing into the
        CLI (`trustify cache clean </dev/null`).
        """
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(home, "abc11111")
            # `input=""` closes stdin immediately, so input() in the
            # child sees EOF and raises EOFError — exactly the
            # </dev/null scenario.
            result = _run_with_home(home, "cache", "clean", input_text="")
            self.assertEqual(result.returncode, 1, f"stderr: {result.stderr}\nstdout: {result.stdout}")
            self.assertIn("Aborted", result.stdout)
            self.assertNotIn("EOFError", result.stderr)
            self.assertTrue((Path(home) / ".cache" / "trustify" / "abc11111").exists())

    def test_cache_rejects_global_schema_flag(self):
        result = _run("--schema", "/abs/dir", "cache", "list")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--schema cannot be used with cache", result.stderr)

    def test_trim_dedupes_same_project_set(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(
                home,
                "trust_a1",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            _make_cache_entry(
                home,
                "trust_a2",
                provenance={
                    "generated_at": "2026-05-19T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            result = _run_with_home(home, "cache", "trim")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("Kept trust_a1", result.stdout)
            self.assertIn("deleted trust_a2", result.stdout)
            self.assertTrue((Path(home) / ".cache" / "trustify" / "trust_a1").exists())
            self.assertFalse((Path(home) / ".cache" / "trustify" / "trust_a2").exists())

    def test_trim_no_duplicates_short_circuits(self):
        with tempfile.TemporaryDirectory() as home:
            _make_cache_entry(
                home,
                "trust_a1",
                provenance={
                    "generated_at": "2026-05-20T18:00:00Z",
                    "trustify_version": "1",
                    "projects": {"trust": "/abs/trust"},
                },
            )
            result = _run_with_home(home, "cache", "trim")
            self.assertEqual(result.returncode, 0)
            self.assertIn("No duplicate", result.stdout)


@unittest.skipUnless(_HAS_ARGCOMPLETE, "argcomplete not installed (optional extra trustify[optional])")
class TestCliInstallCompletion(unittest.TestCase):
    """Tests for `trustify install-completion`. All tests redirect $HOME and
    $XDG_DATA_HOME/$XDG_CONFIG_HOME at a tempdir so we never touch the
    real user dotfiles.
    """

    @staticmethod
    def _run_sandboxed(home, *args):
        env = os.environ.copy()
        env["HOME"] = str(home)
        # Force XDG resolution to the tempdir so default destinations are
        # predictable and confined.
        env["XDG_DATA_HOME"] = str(Path(home) / ".local" / "share")
        env["XDG_CONFIG_HOME"] = str(Path(home) / ".config")
        # Avoid $SHELL bleeding through and changing the default shell.
        env.pop("SHELL", None)
        cmd = [sys.executable, "-m", "trustify.cli", *args]
        return subprocess.run(cmd, env=env, capture_output=True, text=True)

    def test_print_bash_emits_bash_completion_script(self):
        with tempfile.TemporaryDirectory() as home:
            result = self._run_sandboxed(home, "install-completion", "--print", "--shell", "bash")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("_python_argcomplete", result.stdout)
            self.assertIn("trustify", result.stdout)

    def test_print_fish_emits_fish_completion_script(self):
        with tempfile.TemporaryDirectory() as home:
            result = self._run_sandboxed(home, "install-completion", "--print", "--shell", "fish")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("__fish_trustify_complete", result.stdout)

    def test_print_zsh_emits_compdef_script(self):
        with tempfile.TemporaryDirectory() as home:
            result = self._run_sandboxed(home, "install-completion", "--print", "--shell", "zsh")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("compdef", result.stdout)

    def test_bash_default_dest_uses_xdg_data_home(self):
        with tempfile.TemporaryDirectory() as home:
            result = self._run_sandboxed(home, "install-completion", "--shell", "bash")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            dest = Path(home) / ".local/share/bash-completion/completions/trustify"
            self.assertTrue(dest.is_file(), f"expected file at {dest}")
            self.assertIn("_python_argcomplete", dest.read_text())
            self.assertIn(str(dest), result.stdout)

    def test_fish_default_dest_uses_xdg_config_home(self):
        with tempfile.TemporaryDirectory() as home:
            result = self._run_sandboxed(home, "install-completion", "--shell", "fish")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            dest = Path(home) / ".config/fish/completions/trustify.fish"
            self.assertTrue(dest.is_file(), f"expected file at {dest}")
            self.assertIn("__fish_trustify_complete", dest.read_text())

    def test_zsh_default_dest_writes_underscore_autoload_function(self):
        with tempfile.TemporaryDirectory() as home:
            result = self._run_sandboxed(home, "install-completion", "--shell", "zsh")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            dest = Path(home) / ".local/share/zsh/site-functions/_trustify"
            self.assertTrue(dest.is_file(), f"expected file at {dest}")
            # Zsh setup is non-trivial — the hint should tell the user what to do.
            self.assertIn("fpath", result.stdout)
            self.assertIn(str(dest.parent), result.stdout)

    def test_existing_file_without_force_errors(self):
        with tempfile.TemporaryDirectory() as home:
            dest = Path(home) / "trustify-completion"
            dest.write_text("old\n")
            result = self._run_sandboxed(home, "install-completion", "--shell", "bash", "--dest", str(dest))
            self.assertEqual(result.returncode, 2)
            self.assertIn("already exists", result.stderr)
            self.assertEqual(dest.read_text(), "old\n", "file must be untouched")

    def test_force_overwrites_existing_file(self):
        with tempfile.TemporaryDirectory() as home:
            dest = Path(home) / "trustify-completion"
            dest.write_text("old\n")
            result = self._run_sandboxed(home, "install-completion", "--shell", "bash", "--dest", str(dest), "--force")
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertIn("_python_argcomplete", dest.read_text())

    def test_custom_dest_creates_parent_dirs(self):
        with tempfile.TemporaryDirectory() as home:
            dest = Path(home) / "deep" / "nested" / "trustify"
            result = self._run_sandboxed(home, "install-completion", "--shell", "bash", "--dest", str(dest))
            self.assertEqual(result.returncode, 0, f"stderr: {result.stderr}")
            self.assertTrue(dest.is_file())


class TestCliTrustifyDebugTraceback(unittest.TestCase):
    """Audit 5.1: the CLI's top-level except collapses every uncaught
    exception to a one-line stderr message, swallowing the traceback.
    `TRUSTIFY_DEBUG=1` already enables verbose parser/generator logging
    — its documented purpose is "debug trustify itself" — but it didn't
    affect the top-level handler, so contributors hitting an internal
    crash had to drop down to `python -c` to recover the stack. The fix
    re-raises under TRUSTIFY_DEBUG so Python prints the full traceback;
    the end-user one-liner is preserved when the var is unset.

    Drives the test via a Python -c script that mocks
    `trustify.cache.iter_entries` to raise a marker RuntimeError, then
    invokes the `cache list` subcommand (no project-resolution detour).
    """

    _SCRIPT = (
        "import sys\n"
        "from unittest.mock import patch\n"
        "import trustify.cache as _cache\n"
        "import trustify.cli as _cli\n"
        "with patch.object(_cache, 'iter_entries', side_effect=RuntimeError('boom_audit_5_1')):\n"
        "    rc = _cli.main(['cache', 'list'])\n"
        "sys.exit(rc)\n"
    )

    def test_no_trustify_debug_collapses_to_one_line(self):
        # End-user behaviour preserved: one-line stderr, no traceback.
        env = {k: v for k, v in os.environ.items() if k != "TRUSTIFY_DEBUG"}
        result = subprocess.run([sys.executable, "-c", self._SCRIPT], env=env, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("RuntimeError", result.stderr)
        self.assertIn("boom_audit_5_1", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_trustify_debug_re_raises_traceback(self):
        env = {**os.environ, "TRUSTIFY_DEBUG": "1"}
        result = subprocess.run([sys.executable, "-c", self._SCRIPT], env=env, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Traceback", result.stderr)
        self.assertIn("RuntimeError", result.stderr)
        self.assertIn("boom_audit_5_1", result.stderr)


class TestCliInternalErrorAlwaysShowsTraceback(unittest.TestCase):
    """Audit 5.4: `TrustifyInternalError` signals a programmer bug or
    a schema configuration error inside trustify itself (the 5 sites
    in `base.py` that used to `raise Exception('Internal error')`).
    The CLI must surface the full traceback for these — independent
    of TRUSTIFY_DEBUG — because the stack is the only way to locate
    the bug. End-user-facing parse errors (TrustifyParseError) keep
    the one-line behaviour governed by TRUSTIFY_DEBUG.
    """

    _SCRIPT = (
        "import sys\n"
        "from unittest.mock import patch\n"
        "from trustify import TrustifyInternalError\n"
        "import trustify.cache as _cache\n"
        "import trustify.cli as _cli\n"
        "with patch.object(_cache, 'iter_entries', side_effect=TrustifyInternalError('boom_audit_5_4')):\n"
        "    rc = _cli.main(['cache', 'list'])\n"
        "sys.exit(rc)\n"
    )

    def test_internal_error_surfaces_traceback_without_debug(self):
        env = {k: v for k, v in os.environ.items() if k != "TRUSTIFY_DEBUG"}
        result = subprocess.run([sys.executable, "-c", self._SCRIPT], env=env, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Traceback", result.stderr)
        self.assertIn("TrustifyInternalError", result.stderr)
        self.assertIn("boom_audit_5_4", result.stderr)


class TestCliGeneratePdf(unittest.TestCase):
    def test_help_lists_generate_pdf(self):
        result = _run("generate_pdf", "--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("--from-markdown", result.stdout)
        self.assertIn("--keep-build", result.stdout)
        self.assertIn("--doxygen", result.stdout)

    def test_out_is_required(self):
        # --from-markdown given but no --out: argparse rejects (exit 2).
        with tempfile.TemporaryDirectory() as d:
            result = _run("generate_pdf", "--from-markdown", d)
            self.assertEqual(result.returncode, 2)
            self.assertIn("--out", result.stderr)

    def test_bad_from_markdown_dir_exits_two(self):
        # An empty dir is not a generate_markdown output -> precondition
        # (exit 2). generate_pdf probes doxygen + the LaTeX toolchain
        # BEFORE validating --from-markdown, so we stub all three on PATH
        # (fake doxygen via --doxygen; fake pdflatex/make via a PATH dir)
        # to reach the precondition deterministically, regardless of what
        # the test machine actually has installed.
        with tempfile.TemporaryDirectory() as d:
            empty = Path(d) / "empty"
            empty.mkdir()
            out = Path(d) / "manual.pdf"

            stub_dir = Path(d) / "stubbin"
            stub_dir.mkdir()
            doxygen = stub_dir / "doxygen"
            doxygen.write_text("#!/bin/sh\necho 1.16.1\n")
            doxygen.chmod(0o755)
            for tool in ("pdflatex", "make"):
                exe = stub_dir / tool
                exe.write_text("#!/bin/sh\nexit 0\n")
                exe.chmod(0o755)

            env = os.environ.copy()
            env["PATH"] = f"{stub_dir}{os.pathsep}{env.get('PATH', '')}"
            result = _run(
                "generate_pdf",
                "--out",
                str(out),
                "--from-markdown",
                str(empty),
                "--doxygen",
                str(doxygen),
                env=env,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("generate_pdf", result.stderr)


if __name__ == "__main__":
    unittest.main()
