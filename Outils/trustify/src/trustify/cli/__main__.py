"""Allow `python -m trustify.cli` to work alongside the console_scripts entry."""

from trustify.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
