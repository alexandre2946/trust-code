"""Unit tests for `trustify._color`."""

import os
import unittest
from unittest import mock


class TestColorize(unittest.TestCase):
    def test_returns_text_unchanged_when_disabled(self):
        from trustify._color import colorize

        # Force-disabled: no escape sequences anywhere in the output.
        self.assertEqual(colorize("PASSED FAILED SKIPPED", enabled=False), "PASSED FAILED SKIPPED")
        self.assertNotIn("\x1b", colorize("WARNING something", enabled=False))

    def test_wraps_known_status_keywords_when_enabled(self):
        from trustify._color import colorize

        out = colorize("[FAILED] /path/foo.data", enabled=True)
        # Token surrounded by an ANSI red prefix + reset suffix.
        self.assertIn("\x1b[31mFAILED\x1b[0m", out)
        # The brackets and the path stay intact.
        self.assertTrue(out.startswith("["))
        self.assertIn("/path/foo.data", out)

    def test_distinct_colors_per_status(self):
        from trustify._color import colorize

        passed = colorize("PASSED", enabled=True)
        failed = colorize("FAILED", enabled=True)
        skipped = colorize("SKIPPED", enabled=True)
        self.assertIn("\x1b[32m", passed)  # green
        self.assertIn("\x1b[31m", failed)  # red
        self.assertIn("\x1b[33m", skipped)  # yellow
        # All three terminate with the reset code.
        for s in (passed, failed, skipped):
            self.assertTrue(s.endswith("\x1b[0m"), s)

    def test_does_not_color_lowercase_descriptive_words(self):
        """The `Summary: X/Y test(s) passed successfully (N skipped, M
        failed)` line uses lowercase English words — they must NOT be
        colorized (case-sensitive matching).
        """
        from trustify._color import colorize

        line = "Summary: 5/7 test(s) passed successfully (1 skipped, 1 failed)"
        out = colorize(line, enabled=True)
        self.assertEqual(out, line, "lowercase descriptive words must stay plain")

    def test_colors_warning_and_obsolete(self):
        from trustify._color import colorize

        self.assertIn("\x1b[33mWARNING\x1b[0m", colorize("WARNING: stale marker", enabled=True))
        self.assertIn("\x1b[33mobsolete\x1b[0m", colorize("  [obsolete] /path", enabled=True))


class TestColorsEnabled(unittest.TestCase):
    """`colors_enabled()` policy: TTY stdout AND `NO_COLOR` env unset."""

    def test_disabled_when_stdout_not_a_tty(self):
        from trustify import _color

        with mock.patch("sys.stdout") as fake_stdout:
            fake_stdout.isatty.return_value = False
            # NO_COLOR irrelevant when isatty is False.
            with mock.patch.dict(os.environ, {}, clear=False):
                os.environ.pop("NO_COLOR", None)
                self.assertFalse(_color.colors_enabled())

    def test_disabled_when_no_color_env_set(self):
        from trustify import _color

        with mock.patch("sys.stdout") as fake_stdout:
            fake_stdout.isatty.return_value = True
            with mock.patch.dict(os.environ, {"NO_COLOR": "1"}, clear=False):
                self.assertFalse(_color.colors_enabled())

    def test_enabled_on_tty_without_no_color(self):
        from trustify import _color

        with mock.patch("sys.stdout") as fake_stdout:
            fake_stdout.isatty.return_value = True
            with mock.patch.dict(os.environ, {}, clear=False):
                os.environ.pop("NO_COLOR", None)
                self.assertTrue(_color.colors_enabled())

    def test_per_stream_decision(self):
        """`colors_enabled(stream=...)` checks that specific stream's
        TTY status — important for the logger, which writes to stderr.
        Stdout being a TTY must NOT enable coloring when stderr (the
        actual sink) is redirected to a file.
        """
        import sys

        from trustify import _color

        fake_stdout = mock.MagicMock()
        fake_stdout.isatty.return_value = True
        fake_stderr = mock.MagicMock()
        fake_stderr.isatty.return_value = False
        with (
            mock.patch.object(sys, "stdout", fake_stdout),
            mock.patch.object(sys, "stderr", fake_stderr),
            mock.patch.dict(os.environ, {}, clear=False),
        ):
            os.environ.pop("NO_COLOR", None)
            self.assertTrue(_color.colors_enabled())  # default = stdout
            self.assertFalse(_color.colors_enabled(sys.stderr))


class TestPaint(unittest.TestCase):
    def test_wraps_text_in_color_when_enabled(self):
        from trustify._color import RED, RESET, paint

        out = paint("error message", RED, enabled=True)
        self.assertEqual(out, f"{RED}error message{RESET}")

    def test_returns_text_unchanged_when_disabled(self):
        from trustify._color import RED, paint

        self.assertEqual(paint("error message", RED, enabled=False), "error message")
        self.assertNotIn("\x1b", paint("error message", RED, enabled=False))


class TestLoggerRespectsColorsEnabled(unittest.TestCase):
    """The trustify logger's `CustomFormatter` must skip ANSI codes when
    the destination (stderr) is not a TTY. Regression: the formatter
    used to embed color codes unconditionally at class-definition time,
    so log lines captured by CI / piped to a file picked up
    `\\033[1;33m`-style noise.
    """

    def test_no_ansi_when_stderr_not_a_tty(self):
        import logging

        from trustify.core.misc_utilities import CustomFormatter

        record = logging.LogRecord(
            name="trustify",
            level=logging.WARNING,
            pathname=__file__,
            lineno=42,
            msg="something happened",
            args=(),
            exc_info=None,
        )
        fake_stderr = mock.MagicMock()
        fake_stderr.isatty.return_value = False
        import sys

        with mock.patch.object(sys, "stderr", fake_stderr):
            out = CustomFormatter().format(record)
        self.assertNotIn("\x1b", out, f"expected plain text, got: {out!r}")
        self.assertNotIn("\033", out, f"expected plain text, got: {out!r}")

    def test_ansi_when_stderr_is_a_tty(self):
        import logging

        from trustify._color import YELLOW
        from trustify.core.misc_utilities import CustomFormatter

        record = logging.LogRecord(
            name="trustify",
            level=logging.WARNING,
            pathname=__file__,
            lineno=42,
            msg="something happened",
            args=(),
            exc_info=None,
        )
        fake_stderr = mock.MagicMock()
        fake_stderr.isatty.return_value = True
        import sys

        with mock.patch.object(sys, "stderr", fake_stderr), mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("NO_COLOR", None)
            out = CustomFormatter().format(record)
        self.assertIn(YELLOW, out)


if __name__ == "__main__":
    unittest.main()
