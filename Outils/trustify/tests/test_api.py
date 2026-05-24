"""Focused unit tests for trustify.api — Phase 2 of the reshaping.

Most tests in this file invoke the real schema generation against TRUST sources.
"""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class TestApiGenerate(unittest.TestCase):
    def test_generate_against_trust_root_creates_expected_files(self):
        from trustify.api import generate_schema

        trust_root = os.environ["TRUST_ROOT"]
        with tempfile.TemporaryDirectory() as out:
            d = generate_schema(projects=None, trust_root=trust_root, out=Path(out))
            self.assertTrue((d / "trustify_gen.py").exists())
            # Pyd module name is digest-suffixed
            # (`trustify_gen_pyd_<digest>.py`) so distinct schemas can
            # coexist in one process.
            self.assertTrue(list(d.glob("trustify_gen_pyd_*.py")))
            self.assertTrue((d / "TRAD2_trustify").exists())

    def test_generate_is_idempotent_for_same_inputs(self):
        from trustify.api import generate_schema
        from trustify.cache import is_cached

        trust_root = os.environ["TRUST_ROOT"]
        with tempfile.TemporaryDirectory() as out:
            out_p = Path(out)
            generate_schema(projects=None, trust_root=trust_root, out=out_p)
            self.assertTrue(is_cached(out_p))
            mtime_first = (out_p / "trustify_gen.py").stat().st_mtime
            generate_schema(projects=None, trust_root=trust_root, out=out_p)
            mtime_second = (out_p / "trustify_gen.py").stat().st_mtime
            self.assertEqual(
                mtime_first, mtime_second, "second generate_schema(...) must reuse cached output, not regenerate"
            )


class TestApiGenerateProjectScanPath(unittest.TestCase):
    """Each project (trust_root included) must contribute `<path>/src` only."""

    def test_generate_errors_when_project_lacks_project_cfg(self):
        """A directory passed as --projects must be a proper baltik:
        validate_projects fires before any scanning happens, so a bare
        directory without project.cfg is rejected with a clear message.
        """
        from trustify.api import generate_schema
        from trustify.projects import ConfigError

        with tempfile.TemporaryDirectory() as proj:
            with self.assertRaises(ConfigError) as cm:
                generate_schema(projects=[proj], trust_root=None)
            msg = str(cm.exception)
            self.assertIn("project.cfg", msg)
            self.assertIn("baltik_configure", msg)

    def test_generate_passes_only_src_subdirs_to_scanner(self):
        from unittest import mock

        from tests._baltik_fixtures import make_baltik, mock_workspace
        from trustify.api import generate_schema
        from trustify.core import trad2_utilities

        with mock_workspace() as ws:
            tr = ws.root / "trust"
            (tr / "src").mkdir(parents=True)
            p1 = make_baltik(ws.root, "p1")
            # `make_baltik` already creates `src/` (canonical baltik
            # layout, required by `validate_projects`). Decoy dirs the
            # previous logic would have picked up — must be ignored now.
            (p1 / "build" / "src").mkdir(parents=True)
            (p1 / "extra").mkdir()
            with (
                mock.patch.object(
                    trad2_utilities.TRAD2Content,
                    "BuildFromOrgAndSources",
                    side_effect=RuntimeError("stop"),
                ) as m,
                self.assertRaises(RuntimeError),
            ):
                generate_schema(projects=[str(p1)], trust_root=str(tr))
            m.assert_called_once()
            args, _kwargs = m.call_args
            # signature: (trad2_org, src_dirs, trust_root=...)
            self.assertEqual(
                args[1],
                [str(tr / "src"), str(p1 / "src")],
                "generate_schema() must scan only `<path>/src` for each project, in [trust_root, *projects] order",
            )


class TestApiLoadDataset(unittest.TestCase):
    def test_load_dataset_with_projects_returns_pydantic_model(self):
        from trustify.api import load_dataset

        ds = load_dataset(
            os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data"),
            trust_root=os.environ["TRUST_ROOT"],
        )
        self.assertTrue(hasattr(ds, "model_dump") or hasattr(ds, "dict"))


class TestApiCheck(unittest.TestCase):
    def test_check_returns_result_with_status_attribute(self):
        from trustify.api import CheckResult, check

        result = check(
            os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data"),
            trust_root=os.environ["TRUST_ROOT"],
        )
        self.assertIsInstance(result, CheckResult)
        self.assertIn(result.status, ("PASSED", "FAILED", "SKIPPED"))
        self.assertEqual(result.status, "PASSED")

    def test_check_captures_diff_on_round_trip_mismatch(self):
        """When the write-after-parse comparison fails, the unified diff
        produced by `check_str_equality` must land in `CheckResult.diff`
        so the batch-check summary (and `--summary-out` file) can
        surface it without users re-running per-file. Regression: the
        diff used to be printed via `logger.error` and then dropped.
        """
        from unittest import mock

        from trustify.api import check

        ds_path = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data")
        # Force a round-trip mismatch by making `prune_after_end` return
        # something the serialized output cannot match. Anything that
        # differs from `toDatasetTokens()`'s output triggers the
        # check_str_equality branch we care about.
        with mock.patch(
            "trustify.core.misc_utilities.prune_after_end",
            return_value="!!! deliberately wrong expected output !!!",
        ):
            result = check(ds_path, trust_root=os.environ["TRUST_ROOT"])
        self.assertEqual(result.status, "FAILED")
        self.assertEqual(result.message, "dataset is not bit-identical after round-trip")
        self.assertTrue(result.diff, "expected CheckResult.diff to carry the captured unified diff")
        # Spot-check that the diff actually reflects the injected wrong
        # expected text (the unified-diff '+'/'-' markers must show up).
        self.assertIn("deliberately wrong expected output", result.diff)

    def test_batch_check_aggregates_individual_results(self):
        from trustify.api import BatchCheckResult, batch_check

        ds_dir = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets")
        result = batch_check(
            [os.path.join(ds_dir, "diffusion_implicite_jdd6.data")],
            trust_root=os.environ["TRUST_ROOT"],
        )
        self.assertIsInstance(result, BatchCheckResult)
        self.assertEqual(result.total, 1)
        self.assertEqual(result.passed, 1)
        self.assertEqual(result.failed, 0)

    def test_batch_check_parallel_matches_sequential(self):
        # jobs>=2 must produce the same aggregate counts AND the same
        # per-file outcomes in the same order as jobs=1.
        from trustify.api import batch_check

        ds_dir = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets")
        files = [
            os.path.join(ds_dir, "diffusion_implicite_jdd6.data"),
            os.path.join(ds_dir, "Canal_perio_VEF_2D.data"),
            os.path.join(ds_dir, "distance_paroi_jdd1.data"),
        ]
        seq = batch_check(files, trust_root=os.environ["TRUST_ROOT"])
        par = batch_check(files, trust_root=os.environ["TRUST_ROOT"], jobs=2)
        self.assertEqual(
            (seq.total, seq.passed, seq.failed, seq.skipped),
            (par.total, par.passed, par.failed, par.skipped),
        )
        self.assertEqual(
            [(r.path, r.status) for r in seq.results],
            [(r.path, r.status) for r in par.results],
        )

    def test_batch_check_rejects_negative_jobs(self):
        from trustify.api import batch_check

        with self.assertRaisesRegex(ValueError, "jobs must be >= 0"):
            batch_check([], trust_root=os.environ["TRUST_ROOT"], jobs=-1)

    def test_batch_check_jobs_zero_uses_all_cpus(self):
        # `jobs=0` is the standard CLI convention for "use all available
        # CPUs". batch_check must resolve it to os.cpu_count() before
        # handing it to the multiprocessing.Pool — verified by stubbing
        # both os.cpu_count and multiprocessing.get_context so we can
        # capture the actual N that reaches Pool().
        from unittest import mock

        from trustify.api import batch_check

        ds = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets/Canal_perio_VEF_2D.data")
        captured: dict = {}

        class _FakePool:
            def __init__(self, n):
                captured["n"] = n

            def imap(self, fn, work, chunksize=None):  # noqa: ARG002 — match Pool.imap signature
                return iter([fn(w) for w in work])

            def close(self):
                pass

            def join(self):
                pass

            def terminate(self):
                pass

        class _FakeCtx:
            def Pool(self, n):
                return _FakePool(n)

        with (
            mock.patch("os.cpu_count", return_value=7),
            mock.patch("multiprocessing.get_context", return_value=_FakeCtx()),
        ):
            batch_check([ds], trust_root=os.environ["TRUST_ROOT"], jobs=0)
        self.assertEqual(captured["n"], 7)

    def test_check_skips_dataset_with_trustify_not_marker(self):
        # Default behaviour: a dataset carrying `# TRUSTIFY NOT #`
        # short-circuits to SKIPPED without consulting the schema.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT #\nthis is not even valid TRUST\n")
            result = check(str(data), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(result.status, "SKIPPED")
            self.assertIn("TRUSTIFY NOT", result.message)

    def test_check_no_skip_forces_validation_of_marked_dataset(self):
        # no_skip=True bypasses the marker — the dataset goes through
        # the parser and fails on its (intentionally invalid) content.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT #\nthis is not even valid TRUST\n")
            result = check(str(data), trust_root=os.environ["TRUST_ROOT"], no_skip=True)
            self.assertEqual(result.status, "FAILED")

    def test_batch_check_no_skip_threads_through(self):
        # Smoke test: batch_check forwards no_skip to each per-file check.
        from trustify.api import batch_check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT #\nthis is not even valid TRUST\n")
            default = batch_check([str(data)], trust_root=os.environ["TRUST_ROOT"])
            forced = batch_check([str(data)], trust_root=os.environ["TRUST_ROOT"], no_skip=True)
            self.assertEqual(default.skipped, 1)
            self.assertEqual(default.failed, 0)
            self.assertEqual(forced.skipped, 0)
            self.assertEqual(forced.failed, 1)

    def test_check_flags_obsolete_marker_when_passing(self):
        # A passing dataset that carries `# TRUSTIFY NOT #` and is
        # checked with no_skip=True must report obsolete_marker=True
        # so callers can suggest dropping the marker.
        from trustify.api import check

        ds_dir = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets")
        valid = os.path.join(ds_dir, "diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked = Path(d) / "marked.data"
            marked.write_text("# TRUSTIFY NOT #\n" + Path(valid).read_text())
            r = check(str(marked), trust_root=os.environ["TRUST_ROOT"], no_skip=True)
            self.assertEqual(r.status, "PASSED")
            self.assertTrue(r.obsolete_marker)

    def test_check_does_not_flag_obsolete_marker_in_default_mode(self):
        # When no_skip is False, the marked dataset is SKIPPED — no
        # round-trip check runs, so obsolete_marker stays False.
        from trustify.api import check

        ds_dir = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets")
        valid = os.path.join(ds_dir, "diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked = Path(d) / "marked.data"
            marked.write_text("# TRUSTIFY NOT #\n" + Path(valid).read_text())
            r = check(str(marked), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(r.status, "SKIPPED")
            self.assertFalse(r.obsolete_marker)

    def test_check_extracts_marker_justification(self):
        # A `# TRUSTIFY NOT: <text> #` tag carries the justification on
        # the result, in addition to the SKIPPED status.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT: sed template, not real TRUST #\nnot valid TRUST\n")
            r = check(str(data), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(r.status, "SKIPPED")
            self.assertEqual(r.marker_justification, "sed template, not real TRUST")

    def test_check_bare_marker_has_empty_justification(self):
        # Backward compat: the bare sentinel still triggers SKIPPED with
        # an empty justification.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT #\nnot valid TRUST\n")
            r = check(str(data), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(r.status, "SKIPPED")
            self.assertEqual(r.marker_justification, "")

    def test_check_justification_strips_and_collapses_whitespace(self):
        # Multi-space and newline runs inside the justification collapse
        # to single spaces; leading/trailing whitespace is stripped.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "opted_out.data"
            data.write_text("# TRUSTIFY NOT:   foo    bar\n   baz #\nnot valid TRUST\n")
            r = check(str(data), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(r.marker_justification, "foo bar baz")

    def test_check_empty_justification_is_treated_as_bare(self):
        # `# TRUSTIFY NOT: #` and `# TRUSTIFY NOT:    #` are equivalent
        # to the bare form — empty marker_justification.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            for content in (
                "# TRUSTIFY NOT: #\nnot valid TRUST\n",
                "# TRUSTIFY NOT:    #\nnot valid TRUST\n",
            ):
                data = Path(d) / "opted_out.data"
                data.write_text(content)
                r = check(str(data), trust_root=os.environ["TRUST_ROOT"])
                self.assertEqual(r.status, "SKIPPED", f"failed for: {content!r}")
                self.assertEqual(r.marker_justification, "", f"failed for: {content!r}")

    def test_check_malformed_marker_missing_colon_fails(self):
        # `# TRUSTIFY NOT toto #` (missing colon between NOT and the
        # justification) is malformed. The legacy code silently ignored
        # it — the dataset would parse as if no marker was present,
        # which is a surprising correctness hole. The check must FAIL
        # with a clear "ill-formed tag" message.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "malformed.data"
            data.write_text("# TRUSTIFY NOT toto #\nnot valid TRUST\n")
            r = check(str(data), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(r.status, "FAILED")
            self.assertIn("ill-formed", r.message.lower())
            self.assertIn("TRUSTIFY NOT", r.message)

    def test_check_malformed_marker_reports_line_number(self):
        # The error message should pinpoint the line of the malformed
        # tag so the user can find it quickly.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "malformed.data"
            data.write_text("dimension 2\n\n\n# TRUSTIFY NOT bogus #\n")
            r = check(str(data), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(r.status, "FAILED")
            self.assertIn("line 4", r.message)

    def test_check_malformed_marker_alongside_well_formed_still_fails(self):
        # A file with BOTH a well-formed tag AND a malformed one must
        # fail (rather than silently honouring the well-formed one).
        # The malformed form is a typo that the user wants flagged.
        from trustify.api import check

        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "mixed.data"
            data.write_text("# TRUSTIFY NOT: legitimate skip #\n# TRUSTIFY NOT bogus typo #\nnot valid\n")
            r = check(str(data), trust_root=os.environ["TRUST_ROOT"])
            self.assertEqual(r.status, "FAILED")
            self.assertIn("ill-formed", r.message.lower())

    def test_check_no_skip_passing_carries_justification(self):
        # When `--no-skip` is set AND the parse PASSES, the result
        # carries `obsolete_marker=True` AND `marker_justification` (so
        # the CLI obsolete-WARNING line can render it).
        from trustify.api import check

        ds_dir = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets")
        valid = os.path.join(ds_dir, "diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked = Path(d) / "marked.data"
            marked.write_text("# TRUSTIFY NOT: kept for legacy fixture #\n" + Path(valid).read_text())
            r = check(str(marked), trust_root=os.environ["TRUST_ROOT"], no_skip=True)
            self.assertEqual(r.status, "PASSED")
            self.assertTrue(r.obsolete_marker)
            self.assertEqual(r.marker_justification, "kept for legacy fixture")

    def test_check_does_not_flag_obsolete_marker_for_unmarked_pass(self):
        # An unmarked passing dataset shouldn't get the obsolete flag
        # regardless of no_skip.
        from trustify.api import check

        ds_dir = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets")
        valid = os.path.join(ds_dir, "diffusion_implicite_jdd6.data")
        for no_skip in (False, True):
            r = check(valid, trust_root=os.environ["TRUST_ROOT"], no_skip=no_skip)
            self.assertEqual(r.status, "PASSED")
            self.assertFalse(r.obsolete_marker, f"unexpected obsolete flag with no_skip={no_skip}")

    def test_batch_check_obsolete_marker_count(self):
        # A mix: one marked-passing dataset (obsolete) + one marked-failing
        # dataset (legitimate marker, kept). With no_skip=True both run;
        # only the passing one is flagged obsolete.
        from trustify.api import batch_check

        ds_dir = os.path.join(os.environ["TRUST_ROOT"], "Outils/trustify/tests/datasets")
        valid = os.path.join(ds_dir, "diffusion_implicite_jdd6.data")
        with tempfile.TemporaryDirectory() as d:
            marked_ok = Path(d) / "marked_ok.data"
            marked_ok.write_text("# TRUSTIFY NOT #\n" + Path(valid).read_text())
            marked_bad = Path(d) / "marked_bad.data"
            marked_bad.write_text("# TRUSTIFY NOT #\nnot valid TRUST\n")
            r = batch_check([str(marked_ok), str(marked_bad)], trust_root=os.environ["TRUST_ROOT"], no_skip=True)
            self.assertEqual(r.passed, 1)
            self.assertEqual(r.failed, 1)
            obsolete = [c for c in r.results if c.obsolete_marker]
            self.assertEqual(len(obsolete), 1)
            self.assertEqual(obsolete[0].path, str(marked_ok))


class TestExtractMarker(unittest.TestCase):
    """Audit 1.7: `_extract_marker` used to apply its regex to the
    raw file text, so a `# TRUSTIFY NOT #` literal embedded in a
    `/* ... */` block or a `"..."` quoted string would silently opt
    the dataset out of `batch-check`. The fix gates the regex on
    standalone `#...#` comments only — `/* */` blocks and `"..."`
    strings are skipped during the walk."""

    def test_marker_inside_slash_star_block_is_ignored(self):
        from trustify.api import _extract_marker

        text = "/* user wrote # TRUSTIFY NOT # here by accident */\nreal content\n"
        has_marker, justification, err = _extract_marker(text)
        self.assertEqual(err, "")
        self.assertFalse(has_marker, "marker inside /* ... */ must not skip the dataset")
        self.assertEqual(justification, "")

    def test_marker_inside_quoted_string_is_ignored(self):
        from trustify.api import _extract_marker

        # A `"..."` string containing the literal marker text must not
        # opt the dataset out.
        text = 'system "echo # TRUSTIFY NOT # to user"\nreal content\n'
        has_marker, _, err = _extract_marker(text)
        self.assertEqual(err, "")
        self.assertFalse(has_marker, "marker inside a quoted string must not skip the dataset")

    def test_malformed_marker_inside_slash_star_is_not_an_error(self):
        # A `# TRUSTIFY NOT toto #` inside a /* */ block is just user
        # prose — must NOT be flagged as a syntax error.
        from trustify.api import _extract_marker

        text = "/* prose: # TRUSTIFY NOT toto # would be malformed if real */\n"
        _has, _just, err = _extract_marker(text)
        self.assertEqual(err, "", "malformed-prefix check must skip non-comment contexts")

    def test_real_marker_still_detected_with_embedded_text_around_it(self):
        # A genuine standalone `# TRUSTIFY NOT #` must still be
        # detected even when text elsewhere in the file contains
        # quoted or block-commented look-alikes.
        from trustify.api import _extract_marker

        text = (
            '/* decoy: # TRUSTIFY NOT # */\nsystem "another decoy: # TRUSTIFY NOT #"\n# TRUSTIFY NOT #\nreal content\n'
        )
        has_marker, _just, err = _extract_marker(text)
        self.assertEqual(err, "")
        self.assertTrue(has_marker, "real standalone marker must be detected past decoys")

    def test_bare_marker_still_detected(self):
        # Sanity: the simple bare marker is still picked up after the
        # gating refactor.
        from trustify.api import _extract_marker

        text = "# TRUSTIFY NOT #\nirrelevant body\n"
        has_marker, just, err = _extract_marker(text)
        self.assertEqual(err, "")
        self.assertTrue(has_marker)
        self.assertEqual(just, "")

    def test_inline_marker_still_detected(self):
        # Sanity: `# TRUSTIFY NOT: reason #` justification still works.
        from trustify.api import _extract_marker

        text = "# TRUSTIFY NOT: a real reason #\nbody\n"
        has_marker, just, err = _extract_marker(text)
        self.assertEqual(err, "")
        self.assertTrue(has_marker)
        self.assertEqual(just, "a real reason")

    def test_real_malformed_marker_still_flagged(self):
        # Sanity: a malformed marker OUTSIDE any other context must
        # still be reported as an error (the existing safety net).
        from trustify.api import _extract_marker

        text = "# TRUSTIFY NOT toto #\nbody\n"
        _has, _just, err = _extract_marker(text)
        self.assertNotEqual(err, "", "malformed marker outside other contexts must still error")


class TestApiInitConfig(unittest.TestCase):
    def test_init_config_writes_file_with_resolved_projects(self):
        from trustify.api import init_config

        with tempfile.TemporaryDirectory() as d:
            out = init_config(directory=d, projects=["/abs/baltik_a"], trust_root="/abs/trust")
            self.assertEqual(out, Path(d) / ".trustify.json")
            self.assertTrue(out.exists())
            import json

            data = json.loads(out.read_text())
            self.assertEqual(data["workspace"]["projects"], ["/abs/baltik_a"])
            self.assertEqual(data["workspace"]["trust_root"], "/abs/trust")

    def test_init_config_refuses_overwrite_by_default(self):
        from trustify.api import init_config

        with tempfile.TemporaryDirectory() as d:
            existing = Path(d) / ".trustify.json"
            existing.write_text("{}")
            with self.assertRaises(FileExistsError):
                init_config(directory=d, projects=["/x"])

    def test_init_config_overwrites_when_requested(self):
        from trustify.api import init_config

        with tempfile.TemporaryDirectory() as d:
            existing = Path(d) / ".trustify.json"
            existing.write_text("{}")
            init_config(directory=d, projects=["/x"], overwrite=True)
            self.assertIn("/x", existing.read_text())


_MINI_TRAD2 = (
    # Block kept: pb_hydraulique (concrete) — class name + a few attrs +
    # one chaine(into=[...]) enum attr.
    "pb_hydraulique objet_lecture pb_hydraulique|pb_hydro BRACE A problem.\n"
    "  attr tinit floattant t_init OPT initial time\n"
    "  attr tmax floattant tmax REQ max time\n"
    '  attr methode chaine(into=["explicit","implicit"]) methode REQ method\n'
    # Block kept: pb_foo — exercises the enum-value extraction with an
    # unsafe `{` value that must be filtered out.
    "pb_foo objet_lecture pb_foo BRACE another problem\n"
    '  attr aco chaine(into=["{"]) aco REQ opening brace\n'
    # Block dropped: name ends with _base.
    "pb_base objet_lecture pb_base BRACE abstract base.\n"
    "  attr secret chaine secret OPT should still be picked up via attrs\n"
)


def _write_mini_trad2(tmp_dir: Path) -> Path:
    """Write a small synthetic TRAD2_trustify file. Returns the parent dir
    so it can be passed as `schema=`.
    """
    schema_dir = tmp_dir / "schema"
    schema_dir.mkdir()
    (schema_dir / "TRAD2_trustify").write_text(_MINI_TRAD2)
    return schema_dir


class TestApiKeywords(unittest.TestCase):
    """Unit tests for `api.generate_keywords` — drive it via the `schema=` path so we
    don't need a full TRUST_ROOT generation in this test class.
    """

    def test_keywords_writes_both_files(self):
        from trustify.api import generate_keywords

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = _write_mini_trad2(tmp)
            out_dir = tmp / "out"
            generate_keywords(out=out_dir, schema=schema_dir)
            self.assertTrue((out_dir / "Keywords.txt").exists())
            self.assertTrue((out_dir / "Keywords.Vim").exists())

    def test_keywords_txt_has_pipe_prefix_and_sorted_lines(self):
        from trustify.api import generate_keywords

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = _write_mini_trad2(tmp)
            generate_keywords(out=tmp, schema=schema_dir)
            lines = (tmp / "Keywords.txt").read_text().splitlines()
            self.assertTrue(all(line.startswith("|") for line in lines))
            stripped = [line[1:] for line in lines]
            self.assertEqual(stripped, sorted(stripped))
            self.assertEqual(len(set(stripped)), len(stripped), "Keywords.txt must be deduplicated")

    def test_keywords_vim_is_single_syntax_line(self):
        from trustify.api import generate_keywords

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = _write_mini_trad2(tmp)
            generate_keywords(out=tmp, schema=schema_dir)
            vim_text = (tmp / "Keywords.Vim").read_text()
            non_blank = [ln for ln in vim_text.splitlines() if ln.strip()]
            self.assertEqual(len(non_blank), 1)
            self.assertTrue(non_blank[0].startswith("syntax keyword TRUSTLanguageKeywords  "))

    def test_keywords_include_class_attr_and_enum_values(self):
        from trustify.api import generate_keywords

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = _write_mini_trad2(tmp)
            generate_keywords(out=tmp, schema=schema_dir)
            kws = {line[1:] for line in (tmp / "Keywords.txt").read_text().splitlines()}
            # Class names
            self.assertIn("pb_hydraulique", kws)
            self.assertIn("pb_foo", kws)
            # Class synonym
            self.assertIn("pb_hydro", kws)
            # Attribute names
            self.assertIn("tinit", kws)
            self.assertIn("tmax", kws)
            self.assertIn("methode", kws)
            # Attribute synonym
            self.assertIn("t_init", kws)
            # Enum values
            self.assertIn("explicit", kws)
            self.assertIn("implicit", kws)

    def test_keywords_drop_base_and_deriv_classes(self):
        from trustify.api import generate_keywords

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = _write_mini_trad2(tmp)
            generate_keywords(out=tmp, schema=schema_dir)
            kws = {line[1:] for line in (tmp / "Keywords.txt").read_text().splitlines()}
            self.assertNotIn("pb_base", kws)
            self.assertNotIn("objet_lecture", kws, "objet_lecture is the parent class, not a block")

    def test_keywords_filter_unsafe_tokens(self):
        from trustify.api import generate_keywords

        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            schema_dir = _write_mini_trad2(tmp)
            generate_keywords(out=tmp, schema=schema_dir)
            kws = {line[1:] for line in (tmp / "Keywords.txt").read_text().splitlines()}
            self.assertNotIn("{", kws, '{ from chaine(into=["{"]) must be filtered out')

    def test_keywords_missing_trad2_raises(self):
        from trustify.api import generate_keywords

        with tempfile.TemporaryDirectory() as d:
            empty_schema = Path(d) / "empty"
            empty_schema.mkdir()
            with self.assertRaises(FileNotFoundError):
                generate_keywords(out=Path(d) / "out", schema=empty_schema)


class TestApiFormatDataset(unittest.TestCase):
    def test_format_dataset_text_returns_string(self):
        from trustify import format_dataset

        result = format_dataset("dimension 2\n")
        self.assertIsInstance(result, str)

    def test_format_dataset_file_reads_from_disk(self):
        from trustify.api import format_dataset_file

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "x.data"
            f.write_text("dimension 2\n{\nfoo 1\n}\n")
            out = format_dataset_file(f)
            self.assertIn("    foo 1", out)


class TestPickMpStartMethod(unittest.TestCase):
    """Audit 2.7: the batch-check Pool used a hardcoded
    `multiprocessing.get_context("fork")`. Fork inherits the parent's
    full state, including C-extension internal locks held by other
    threads — pydantic + the schema import sit on top of several such
    C extensions, so fork-then-imap can deadlock unpredictably.
    Python 3.14 even moves the macOS default away from fork.

    The picker prefers `forkserver` (clean child process, available
    on Linux + macOS) and falls back to `spawn` (Windows). An env
    override lets users work around any forkserver-specific issue
    without code changes."""

    @staticmethod
    def _set_env(name, value):
        previous = os.environ.get(name)

        def restore():
            if previous is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = previous

        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value
        return restore

    def test_env_override_wins(self):
        from trustify.api import _pick_mp_start_method

        restore = self._set_env("TRUSTIFY_MP_START_METHOD", "spawn")
        try:
            self.assertEqual(_pick_mp_start_method(), "spawn")
        finally:
            restore()

    def test_env_empty_string_falls_back_to_default(self):
        from trustify.api import _pick_mp_start_method

        restore = self._set_env("TRUSTIFY_MP_START_METHOD", "")
        try:
            # Empty string is "not set" — same NO_COLOR convention.
            method = _pick_mp_start_method()
            self.assertIn(method, ("forkserver", "spawn"))
            self.assertNotEqual(method, "")
        finally:
            restore()

    def test_prefers_forkserver_when_available(self):
        from unittest import mock

        from trustify.api import _pick_mp_start_method

        restore = self._set_env("TRUSTIFY_MP_START_METHOD", None)
        try:
            with mock.patch(
                "multiprocessing.get_all_start_methods",
                return_value=["fork", "forkserver", "spawn"],
            ):
                self.assertEqual(_pick_mp_start_method(), "forkserver")
        finally:
            restore()

    def test_falls_back_to_spawn_when_no_forkserver(self):
        from unittest import mock

        from trustify.api import _pick_mp_start_method

        restore = self._set_env("TRUSTIFY_MP_START_METHOD", None)
        try:
            with mock.patch(
                "multiprocessing.get_all_start_methods",
                return_value=["spawn"],
            ):
                self.assertEqual(_pick_mp_start_method(), "spawn")
        finally:
            restore()

    def test_never_picks_fork_silently(self):
        # Even on a platform that only advertises `fork` + `spawn`
        # (no forkserver), prefer `spawn` over `fork`. Fork-after-
        # C-extension-import is the entire reason for this picker.
        from unittest import mock

        from trustify.api import _pick_mp_start_method

        restore = self._set_env("TRUSTIFY_MP_START_METHOD", None)
        try:
            with mock.patch(
                "multiprocessing.get_all_start_methods",
                return_value=["fork", "spawn"],
            ):
                self.assertEqual(_pick_mp_start_method(), "spawn")
        finally:
            restore()


class TestApiBatchFormat(unittest.TestCase):
    """Audit 1.10: `batch_format` is opt-in destructive.

    Default behaviour is dry-run — returns a report carrying the
    unified diff for every file that would change, but writes
    nothing. `apply=True` switches to in-place rewrite (atomic via
    temp-file + os.replace, matching modernize). The pre-fix
    behaviour silently rewrote every file in place — combined with
    the baltik `trustify_format_tests` target that had no dry-run
    guard, a misclick on an unsaved working tree could destroy
    reference inputs.
    """

    SAMPLE = "dimension 2\n{\nfoo 1\n}\n"
    EXPECTED = "dimension 2\n{\n    foo 1\n}\n"

    def test_default_is_dry_run_no_files_modified(self):
        from trustify.api import batch_format

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            batch_format([f])
            # Source must NOT be modified on the dry-run code path.
            self.assertEqual(f.read_text(), self.SAMPLE)

    def test_default_returns_diffs_in_result(self):
        from trustify.api import batch_format

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            result = batch_format([f])
            changed = result.changed_files
            self.assertEqual(len(changed), 1)
            self.assertEqual(str(changed[0].path), str(f))
            # The diff captures the indentation we'd inject.
            self.assertIn("foo 1", changed[0].diff)
            self.assertIn("+", changed[0].diff)
            self.assertIn("-", changed[0].diff)

    def test_apply_writes_files(self):
        from trustify.api import batch_format

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            batch_format([f], apply=True)
            self.assertEqual(f.read_text(), self.EXPECTED)

    def test_already_formatted_file_produces_no_diff(self):
        from trustify.api import batch_format

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.EXPECTED)
            result = batch_format([f])
            self.assertEqual(result.changed_files, [])

    def test_apply_preserves_original_when_atomic_rename_fails(self):
        """Apply must use atomic temp+rename so a write failure
        cannot truncate the source file (audit 2.3 pattern reused)."""
        from unittest.mock import patch

        from trustify.api import batch_format

        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "a.data"
            f.write_text(self.SAMPLE)
            with (
                patch("trustify.modernize.os.replace", side_effect=OSError("simulated rename failure")),
                self.assertRaises(OSError),
            ):
                batch_format([f], apply=True)
            self.assertEqual(f.read_text(), self.SAMPLE, "source corrupted by failed apply")


class TestTopLevelReExports(unittest.TestCase):
    """Phase-4 (LSP migration): trustify's top-level namespace now exposes
    the parser + base + exception primitives the LSP uses, so callers no
    longer have to reach into `trustify.core.*`. These names previously
    only worked via the Phase-2 shim files; this test guards the public
    contract once the shims are gone (last task of this plan)."""

    def test_dataset_parser_importable_from_top_level(self):
        from trustify import Dataset_Parser

        self.assertTrue(isinstance(Dataset_Parser, type))

    def test_trust_parser_and_stream_importable_from_top_level(self):
        from trustify import TRUSTParser, TRUSTStream

        self.assertTrue(isinstance(TRUSTParser, type))
        self.assertTrue(isinstance(TRUSTStream, type))

    def test_trustify_parse_error_importable_from_top_level(self):
        from trustify import TrustifyParseError

        self.assertTrue(isinstance(TrustifyParseError, type))
        # Must be a real Exception subclass — the LSP catches it.
        self.assertTrue(issubclass(TrustifyParseError, Exception))


class TestDataFileUtf8EncodingUnderCLocale(unittest.TestCase):
    """Audit 4.6: every public api.py entry point that opens a `.data`
    file must pass `encoding="utf-8"`. The platform default codec falls
    back to ASCII under LANG=C / LC_ALL=POSIX, so any non-ASCII byte
    in a `.data` comment raises `UnicodeDecodeError` and the entire
    check / format pipeline breaks on locale-misconfigured machines.
    `modernize._atomic_write_text` already uses `encoding="utf-8"` —
    the inconsistency on the read side is the tell.

    Same subprocess pattern as the trad2_utilities tests (audit 4.5):
    LANG=C, LC_ALL=C, PYTHONUTF8=0 disables PEP 540 UTF-8 mode so the
    interpreter really exhibits the broken-ASCII default.
    """

    def _run_under_c_locale(self, script: str) -> subprocess.CompletedProcess:
        env = {**os.environ, "LANG": "C", "LC_ALL": "C", "PYTHONUTF8": "0"}
        return subprocess.run(
            [sys.executable, "-c", script],
            env=env,
            capture_output=True,
            text=True,
        )

    def test_format_dataset_file_handles_non_ascii_under_c_locale(self):
        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "x.data"
            # Plain comment — formatter doesn't need a schema, so this
            # exercises the open() in isolation.
            data.write_text("# Décrit ma classe #\n", encoding="utf-8")
            result = self._run_under_c_locale(
                f"from trustify.api import format_dataset_file; format_dataset_file({str(data)!r})"
            )
            self.assertEqual(
                result.returncode,
                0,
                f"format_dataset_file failed under LANG=C; the open() must pass "
                f"encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )

    def test_batch_format_handles_non_ascii_under_c_locale(self):
        with tempfile.TemporaryDirectory() as d:
            data = Path(d) / "x.data"
            data.write_text("# Décrit ma classe #\n", encoding="utf-8")
            result = self._run_under_c_locale(f"from trustify.api import batch_format; batch_format([{str(data)!r}])")
            self.assertEqual(
                result.returncode,
                0,
                f"batch_format failed under LANG=C; the open() must pass "
                f"encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )

    def test_load_dataset_handles_non_ascii_under_c_locale(self):
        if "TRUST_ROOT" not in os.environ:
            self.skipTest("requires TRUST_ROOT (schema needed to load .data)")
        # Take a real dataset, prepend a non-ASCII comment, save copy.
        ds_src = Path(os.environ["TRUST_ROOT"]) / "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data"
        with tempfile.TemporaryDirectory() as d:
            target = Path(d) / "diffusion_with_accent.data"
            target.write_text(
                "# Décrit ce dataset #\n" + ds_src.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            trust_root = os.environ["TRUST_ROOT"]
            result = self._run_under_c_locale(
                f"from trustify.api import load_dataset; load_dataset({str(target)!r}, trust_root={trust_root!r})"
            )
            self.assertEqual(
                result.returncode,
                0,
                f"load_dataset failed under LANG=C; the open() reading the .data file "
                f"must pass encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )

    def test_check_handles_non_ascii_under_c_locale(self):
        if "TRUST_ROOT" not in os.environ:
            self.skipTest("requires TRUST_ROOT (schema needed to check .data)")
        ds_src = Path(os.environ["TRUST_ROOT"]) / "Outils/trustify/tests/datasets/diffusion_implicite_jdd6.data"
        with tempfile.TemporaryDirectory() as d:
            target = Path(d) / "diffusion_with_accent.data"
            target.write_text(
                "# Décrit ce dataset #\n" + ds_src.read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            trust_root = os.environ["TRUST_ROOT"]
            # check() catches OSError and stuffs the message into the
            # CheckResult, but UnicodeDecodeError is not an OSError —
            # it propagates, so a non-zero returncode = bug present.
            result = self._run_under_c_locale(
                f"from trustify.api import check; check({str(target)!r}, trust_root={trust_root!r})"
            )
            self.assertEqual(
                result.returncode,
                0,
                f"check failed under LANG=C; the open() reading the .data file must pass "
                f"encoding='utf-8'.\n--- stderr ---\n{result.stderr}",
            )


class TestGenerateSchemaRefusesSharedCacheOnUnknownVersion(unittest.TestCase):
    """Audit 6.11: when `importlib.metadata.version("trustify")` fails,
    `__version__` falls back to the sentinel `"0.0.0+unknown"`. That
    sentinel went straight into the cache hash, so every broken
    install shared one cache bucket — two unrelated trustify trees
    would mistakenly serve each other's generated schema.

    With the fix, generate_schema(out=None) refuses to use the shared
    cache when the version is the unknown sentinel, raising a clear
    error pointing the user at `--out` or at fixing their install.
    generate_schema(out=<dir>) (the user-controlled bypass path)
    still works — it never touches the shared cache, so the bucket-
    sharing concern doesn't apply.
    """

    def _patch_unknown_version(self):
        """Force __version__ to the unknown sentinel for the duration
        of the test. _version() inside api.py does a late `from
        trustify import __version__` per call, so monkey-patching the
        module attribute is enough.
        """
        import trustify

        saved = trustify.__version__
        trustify.__version__ = "0.0.0+unknown"
        return saved, trustify

    def test_shared_cache_path_refuses_with_unknown_version(self):
        from trustify.api import generate_schema

        saved, trustify_mod = self._patch_unknown_version()
        try:
            with self.assertRaises(Exception) as cm:
                generate_schema(projects=None, trust_root=os.environ["TRUST_ROOT"], out=None)
            msg = str(cm.exception)
            # The message names the situation and offers concrete
            # remediation paths.
            self.assertIn("unknown", msg.lower())
            self.assertIn("--out", msg)
        finally:
            trustify_mod.__version__ = saved

    def test_explicit_out_path_still_works_with_unknown_version(self):
        # `--out=<dir>` is the documented bypass — the user owns the
        # destination, no shared-cache bucket coordination needed, so
        # the unknown-version refusal must NOT fire here.
        from trustify.api import generate_schema

        saved, trustify_mod = self._patch_unknown_version()
        try:
            with tempfile.TemporaryDirectory() as out:
                d = generate_schema(projects=None, trust_root=os.environ["TRUST_ROOT"], out=Path(out))
                self.assertTrue((d / "trustify_gen.py").exists())
        finally:
            trustify_mod.__version__ = saved


if __name__ == "__main__":
    unittest.main()
