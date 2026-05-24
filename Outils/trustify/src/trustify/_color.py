"""ANSI coloring for trustify output.

Two flavours:

- `colorize(text)` — wraps the recognised status keywords (`PASSED`,
  `FAILED`, `SKIPPED`, `WARNING`, `obsolete`) anywhere they appear in
  `text`. Used by the CLI summary block and by `batch_check`'s verbose
  progress lines.
- `paint(text, color)` — wraps the entire `text` in a specific color
  code. Used by the logger formatter and a couple of error renderers
  that need precise control (e.g. only color the filename, not the
  whole message).

Both gate on `colors_enabled(stream)`: ANSI codes are emitted only when
the target stream is a TTY AND the `NO_COLOR` environment variable
(https://no-color.org) is not set. Decisions are made per call (no
cache), so a runtime change to `sys.stdout` / `sys.stderr` / the env
takes effect immediately — important for tests that monkeypatch any of
these.

Buffers built elsewhere (notably the `batch-check` summary list) are
kept ANSI-free; colorization happens at print time only. `--summary-out`
files and other log-style sinks therefore never see escape sequences.
"""

from __future__ import annotations

import os
import re
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import TextIO


# Full palette — kept compact and named after the underlying SGR codes
# so a reader can map them to `man console_codes(4)` without ceremony.
BLACK = "\033[0;30m"
RED = "\033[0;31m"
GREEN = "\033[0;32m"
BROWN = "\033[0;33m"
BLUE = "\033[0;34m"
PURPLE = "\033[0;35m"
CYAN = "\033[0;36m"
YELLOW = "\033[1;33m"
BOLD = "\033[1m"
NEGATIVE = "\033[7m"
RESET = "\033[0m"

_KEYWORD_COLORS = {
    "PASSED": "\x1b[32m",  # green
    "FAILED": "\x1b[31m",  # red
    "SKIPPED": "\x1b[33m",  # yellow
    "WARNING": "\x1b[33m",  # yellow (same hue as SKIPPED — both "attention, not failure")
    "obsolete": "\x1b[33m",  # yellow (lowercase variant in `[obsolete]` summary rows)
}
_PATTERN = re.compile(r"\b(" + "|".join(re.escape(k) for k in _KEYWORD_COLORS) + r")\b")


def colors_enabled(stream: TextIO | None = None) -> bool:
    """True when `stream` is an interactive terminal AND `NO_COLOR` is
    unset. `stream=None` defaults to `sys.stdout`; pass `sys.stderr`
    explicitly when deciding for a stderr-bound consumer (e.g. the
    logger's `StreamHandler`).

    Checked on every call — no cache — so a test that monkeypatches
    the stream or the env sees the change immediately.
    """
    if os.environ.get("NO_COLOR"):
        return False
    if stream is None:
        stream = sys.stdout
    return stream.isatty()


def paint(text: str, color: str, *, enabled: bool | None = None) -> str:
    """Wrap the entire `text` with the given color code + reset suffix.

    `enabled=None` defers to `colors_enabled()` (stdout). Pass `True` /
    `False` to force a decision regardless of the runtime stream
    state — used by unit tests, and by callers that already resolved
    a stream-specific decision (e.g. the logger formatter passing
    `colors_enabled(sys.stderr)`).
    """
    if enabled is None:
        enabled = colors_enabled()
    if not enabled:
        return text
    return f"{color}{text}{RESET}"


def colorize(text: str, *, enabled: bool | None = None) -> str:
    """Wrap recognised status tokens (`PASSED`, `FAILED`, `SKIPPED`,
    `WARNING`, `obsolete`) anywhere they appear in `text`.

    `enabled=None` defers to `colors_enabled()`. Passing `True` /
    `False` forces the decision — useful in unit tests that need to
    exercise the colored path regardless of how the test harness
    redirects stdout.
    """
    if enabled is None:
        enabled = colors_enabled()
    if not enabled:
        return text
    return _PATTERN.sub(lambda m: f"{_KEYWORD_COLORS[m.group(1)]}{m.group(1)}{RESET}", text)
