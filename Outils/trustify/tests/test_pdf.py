"""Unit tests for trustify.pdf — the generate_pdf orchestration.

doxygen / pdflatex / make are external tools absent from CI, so every
subprocess call and PATH lookup is mocked. These tests exercise the
orchestration logic (resolution, validation, env wiring, error
surfacing), never the real toolchain.
"""

import os
import subprocess
import unittest
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock

from trustify import api, pdf


def _make_markdown_dir(root: Path) -> Path:
    """Build a minimal `generate_markdown`-shaped dir at root (created)."""
    root.mkdir(parents=True, exist_ok=True)
    (root / "keyword_reference.md").write_text("@page foo Foo\n")
    (root / "kw_foo.md").write_text("# foo\n")
    (root / "figures").mkdir(exist_ok=True)
    return root


class TestValidateMarkdownDir(unittest.TestCase):
    def test_accepts_valid_markdown_dir(self):
        with TemporaryDirectory() as d:
            md = _make_markdown_dir(Path(d))
            self.assertEqual(pdf.validate_markdown_dir(str(md)), md.resolve())

    def test_rejects_non_directory(self):
        with TemporaryDirectory() as d:
            missing = Path(d) / "nope"
            with self.assertRaises(ValueError):
                pdf.validate_markdown_dir(str(missing))

    def test_rejects_dir_without_keyword_reference(self):
        with TemporaryDirectory() as d:
            (Path(d) / "kw_foo.md").write_text("# foo\n")
            with self.assertRaises(ValueError):
                pdf.validate_markdown_dir(d)

    def test_rejects_dir_without_kw_files(self):
        with TemporaryDirectory() as d:
            (Path(d) / "keyword_reference.md").write_text("@page foo Foo\n")
            with self.assertRaises(ValueError):
                pdf.validate_markdown_dir(d)

    def test_missing_figures_is_warning_not_error(self):
        with TemporaryDirectory() as d:
            (Path(d) / "keyword_reference.md").write_text("@page foo Foo\n")
            (Path(d) / "kw_foo.md").write_text("# foo\n")
            # No figures/ dir: must still validate, only warn on stderr.
            self.assertEqual(pdf.validate_markdown_dir(d), Path(d).resolve())


class TestEmbeddedDoxyfile(unittest.TestCase):
    def test_doxyfile_is_path_free_and_uses_env_vars(self):
        self.assertIn("$(KW_INPUT_DIR)", pdf._DOXYFILE)
        self.assertIn("$(KW_OUTPUT_DIR)", pdf._DOXYFILE)
        self.assertIn("GENERATE_LATEX         = YES", pdf._DOXYFILE)


class TestResolveDoxygen(unittest.TestCase):
    @staticmethod
    def _cp(version_line):
        """A CompletedProcess whose stdout carries a doxygen version line."""
        return subprocess.CompletedProcess([], 0, stdout=version_line, stderr="")

    def test_explicit_flag_wins_and_version_accepted(self):
        with mock.patch.object(pdf.subprocess, "run", return_value=self._cp("1.16.1\n")):
            self.assertEqual(pdf.resolve_doxygen("/opt/doxygen"), Path("/opt/doxygen"))

    def test_env_var_used_when_no_explicit(self):
        with (
            mock.patch.dict(os.environ, {"TRUST_DOXYGEN_BINARY": "/env/doxygen"}, clear=False),
            mock.patch.object(pdf.subprocess, "run", return_value=self._cp("1.9.0\n")),
        ):
            self.assertEqual(pdf.resolve_doxygen(None), Path("/env/doxygen"))

    def test_path_lookup_when_no_explicit_or_env(self):
        with (
            mock.patch.dict(os.environ, {}, clear=True),
            mock.patch.object(pdf.shutil, "which", return_value="/usr/bin/doxygen"),
            mock.patch.object(pdf.subprocess, "run", return_value=self._cp("1.10\n")),
        ):
            self.assertEqual(pdf.resolve_doxygen(None), Path("/usr/bin/doxygen"))

    def test_not_found_raises(self):
        with (
            mock.patch.dict(os.environ, {}, clear=True),
            mock.patch.object(pdf.shutil, "which", return_value=None),
            self.assertRaises(RuntimeError),
        ):
            pdf.resolve_doxygen(None)

    def test_non_doxygen_binary_rejected(self):
        with (
            mock.patch.object(pdf.subprocess, "run", return_value=self._cp("not a version\n")),
            self.assertRaises(RuntimeError),
        ):
            pdf.resolve_doxygen("/bin/false")


class TestProbeLatexToolchain(unittest.TestCase):
    def test_passes_when_both_present(self):
        with mock.patch.object(pdf.shutil, "which", side_effect=lambda n: f"/usr/bin/{n}"):
            pdf.probe_latex_toolchain()  # must not raise

    def test_raises_when_pdflatex_missing(self):
        def which(name):
            return None if name == "pdflatex" else f"/usr/bin/{name}"

        with (
            mock.patch.object(pdf.shutil, "which", side_effect=which),
            self.assertRaises(RuntimeError) as cm,
        ):
            pdf.probe_latex_toolchain()
        self.assertIn("pdflatex", str(cm.exception))

    def test_raises_when_make_missing(self):
        def which(name):
            return None if name == "make" else f"/usr/bin/{name}"

        with (
            mock.patch.object(pdf.shutil, "which", side_effect=which),
            self.assertRaises(RuntimeError) as cm,
        ):
            pdf.probe_latex_toolchain()
        self.assertIn("make", str(cm.exception))


class TestRunDoxygen(unittest.TestCase):
    def test_writes_doxyfile_and_sets_env(self):
        with TemporaryDirectory() as d:
            build = Path(d)
            with mock.patch.object(pdf.subprocess, "run") as run:
                pdf.run_doxygen(Path("/usr/bin/doxygen"), Path("/in/md"), build)
            # Doxyfile written into the build dir, content == embedded twin.
            self.assertEqual((build / pdf.DOXYFILE_NAME).read_text(), pdf._DOXYFILE)
            # Env wiring + cwd captured from the recorded call.
            args, kwargs = run.call_args
            self.assertEqual(kwargs["env"]["KW_INPUT_DIR"], "/in/md")
            self.assertEqual(kwargs["env"]["KW_OUTPUT_DIR"], str(build))
            self.assertEqual(kwargs["cwd"], build)
            self.assertEqual(args[0][0], "/usr/bin/doxygen")

    def test_nonzero_exit_raises_runtimeerror(self):
        with (
            TemporaryDirectory() as d,
            mock.patch.object(pdf.subprocess, "run", side_effect=subprocess.CalledProcessError(2, ["doxygen"])),
            self.assertRaises(RuntimeError) as cm,
        ):
            pdf.run_doxygen(Path("doxygen"), Path("/in"), Path(d))
        self.assertIn("doxygen exited with status 2", str(cm.exception))


class TestRunLatex(unittest.TestCase):
    def test_success_removes_stale_pdf_and_runs_make(self):
        with TemporaryDirectory() as d:
            build = Path(d)
            latex = build / "latex"
            latex.mkdir()
            stale = latex / "refman.pdf"
            stale.write_text("old")
            with mock.patch.object(
                pdf.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, stdout="", stderr="")
            ) as run:
                pdf.run_latex(build)
            self.assertFalse(stale.exists())  # stale pdf dropped first
            args, _ = run.call_args
            self.assertEqual(args[0], ["make", "-C", str(latex)])

    def test_failure_raises_and_surfaces_log_tail(self):
        with (
            TemporaryDirectory() as d,
            mock.patch.object(
                pdf.subprocess,
                "run",
                return_value=subprocess.CompletedProcess([], 2, stdout="make says", stderr="boom"),
            ),
            self.assertRaises(RuntimeError) as cm,
        ):
            build = Path(d)
            latex = build / "latex"
            latex.mkdir()
            (latex / "refman.log").write_text("\n".join(f"line{i}" for i in range(100)))
            pdf.run_latex(build)
        self.assertIn("make", str(cm.exception))


class TestFinalizePdf(unittest.TestCase):
    def test_copies_refman_to_output(self):
        with TemporaryDirectory() as d:
            build = Path(d) / "build"
            (build / "latex").mkdir(parents=True)
            (build / "latex" / "refman.pdf").write_text("PDFDATA")
            out = Path(d) / "manual.pdf"
            pdf.finalize_pdf(build, out)
            self.assertEqual(out.read_text(), "PDFDATA")

    def test_missing_refman_raises(self):
        with TemporaryDirectory() as d:
            build = Path(d) / "build"
            (build / "latex").mkdir(parents=True)
            with self.assertRaises(RuntimeError):
                pdf.finalize_pdf(build, Path(d) / "manual.pdf")

    def test_missing_output_parent_raises(self):
        with TemporaryDirectory() as d:
            build = Path(d) / "build"
            (build / "latex").mkdir(parents=True)
            (build / "latex" / "refman.pdf").write_text("PDFDATA")
            out = Path(d) / "nope" / "manual.pdf"  # parent does not exist
            with self.assertRaises(RuntimeError):
                pdf.finalize_pdf(build, out)


class TestApiGeneratePdf(unittest.TestCase):
    """generate_pdf with all toolchain steps patched — no doxygen/latex."""

    @staticmethod
    def _patch_toolchain(stack):
        """Patch the four pdf.* steps to no-ops; return the finalize mock.

        finalize writes a stub PDF so the output path exists for the
        completion log's stat() and for content assertions.
        """
        stack.enter_context(mock.patch.object(pdf, "resolve_doxygen", return_value=Path("doxygen")))
        stack.enter_context(mock.patch.object(pdf, "probe_latex_toolchain"))
        stack.enter_context(mock.patch.object(pdf, "run_doxygen"))
        stack.enter_context(mock.patch.object(pdf, "run_latex"))
        return stack.enter_context(
            mock.patch.object(pdf, "finalize_pdf", side_effect=lambda *a: Path(a[1]).write_text("PDF"))
        )

    def test_from_markdown_skips_generation(self):
        with TemporaryDirectory() as d, ExitStack() as stack:
            md = _make_markdown_dir(Path(d) / "md")
            out = Path(d) / "manual.pdf"
            fin = self._patch_toolchain(stack)
            gm = stack.enter_context(mock.patch.object(api, "generate_markdown"))
            result = api.generate_pdf(out=str(out), from_markdown=str(md))
            gm.assert_not_called()  # reused existing markdown
            fin.assert_called_once()
            self.assertEqual(result, out.resolve())
            # finalize received the resolved output path as its 2nd positional arg.
            self.assertEqual(fin.call_args.args[1], out.resolve())

    def test_end_to_end_generates_markdown_into_tempdir(self):
        with TemporaryDirectory() as d, ExitStack() as stack:
            out = Path(d) / "manual.pdf"
            self._patch_toolchain(stack)
            gm = stack.enter_context(mock.patch.object(api, "generate_markdown"))
            result = api.generate_pdf(out=str(out), trust_root=None)
            self.assertEqual(result, out.resolve())
            gm.assert_called_once()
            md_used = Path(gm.call_args.kwargs["out"])
            # The temp markdown dir is cleaned up by default.
            self.assertFalse(md_used.exists())

    def test_emits_progress_logs_to_stderr(self):
        import io
        from contextlib import redirect_stderr

        with TemporaryDirectory() as d, ExitStack() as stack:
            md = _make_markdown_dir(Path(d) / "md")
            out = Path(d) / "manual.pdf"
            self._patch_toolchain(stack)
            stack.enter_context(mock.patch.object(api, "generate_markdown"))
            buf = io.StringIO()
            with redirect_stderr(buf):
                api.generate_pdf(out=str(out), from_markdown=str(md))
            err = buf.getvalue()
            # Phase markers + a completion line, all on stderr.
            self.assertIn("[1/3]", err)
            self.assertIn("[2/3]", err)
            self.assertIn("[3/3]", err)
            self.assertIn("doxygen", err)
            self.assertIn("LaTeX", err)
            self.assertIn("wrote", err)

    def test_missing_out_parent_raises_precondition(self):
        # The out-parent check runs first, before the toolchain probe,
        # so no step patches are needed.
        with TemporaryDirectory() as d:
            mdp = _make_markdown_dir(Path(d) / "md")
            out = Path(d) / "nope" / "manual.pdf"  # parent does not exist
            with self.assertRaises(api._PreconditionError):
                api.generate_pdf(out=str(out), from_markdown=str(mdp))

    def test_bad_from_markdown_raises_precondition(self):
        with TemporaryDirectory() as d, ExitStack() as stack:
            out = Path(d) / "manual.pdf"
            bad = Path(d) / "empty"
            bad.mkdir()  # exists but not a generate_markdown output
            stack.enter_context(mock.patch.object(pdf, "resolve_doxygen", return_value=Path("doxygen")))
            stack.enter_context(mock.patch.object(pdf, "probe_latex_toolchain"))
            with self.assertRaises(api._PreconditionError):
                api.generate_pdf(out=str(out), from_markdown=str(bad))
