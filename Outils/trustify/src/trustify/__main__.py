"""Enable `python -m trustify ...` as an alias for the `trustify` console script.

The docs build invokes trustify through a specific interpreter
(`$(DOC_PYTHON) -m trustify`) so it works regardless of whether the package
was installed into a conda env, an active venv, or a `--user` site (where the
console script lands in the user bin dir, not next to the interpreter).
"""

from trustify.cli import main

raise SystemExit(main())
