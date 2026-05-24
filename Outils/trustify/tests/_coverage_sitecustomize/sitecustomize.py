# Coverage hook for subprocess-spawning tests.
#
# Many trustify tests spawn a fresh Python via subprocess.run(
# [sys.executable, "-c", ...]) to exercise CLI behaviour, LANG=C
# fallbacks, concurrent imports, etc. Code executed in those
# subprocesses is invisible to the parent's `coverage run` unless
# each subprocess starts its own coverage tracker.
#
# `make coverage` arranges for two things:
#   - PYTHONPATH includes the directory holding THIS file, so it's
#     found and executed at every Python startup;
#   - COVERAGE_PROCESS_START points at the pyproject.toml carrying
#     the [tool.coverage.run] config.
#
# `coverage.process_startup()` is a no-op when COVERAGE_PROCESS_START
# is unset, so this file is safe to drop into any sys.path — it has
# zero effect outside coverage runs.
import coverage

coverage.process_startup()
