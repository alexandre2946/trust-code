"""Runs every shared corpus pair through the Python formatter.

The corpus (at tests/corpus/) is the behavioral oracle shared with
the TypeScript test suite. Byte-identical output is required from
both implementations.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from trustify import format_dataset

CORPUS_DIR = Path(__file__).resolve().parent / "corpus"
INPUTS_DIR = CORPUS_DIR / "inputs"
EXPECTED_DIR = CORPUS_DIR / "expected"


def _corpus_cases() -> list[str]:
    if not INPUTS_DIR.is_dir():
        return []
    return sorted(p.name for p in INPUTS_DIR.glob("*.data"))


CASES = _corpus_cases()
assert CASES, f"No corpus cases found in {INPUTS_DIR}"


class TestCorpus(unittest.TestCase):
    """One test method per corpus case, generated below."""


def _make_corpus_test(case_name: str):
    def test(self):
        input_path = INPUTS_DIR / case_name
        expected_path = EXPECTED_DIR / case_name
        self.assertTrue(
            expected_path.exists(),
            f"missing expected file for corpus input '{case_name}'",
        )
        input_text = input_path.read_text(encoding="utf-8")
        expected_text = expected_path.read_text(encoding="utf-8")
        actual = format_dataset(input_text)
        self.assertEqual(actual, expected_text)

    return test


# Generate one test method per corpus case. Method names must be valid
# Python identifiers, so swap dots/dashes for underscores.
for _case_name in CASES:
    _safe = _case_name.replace(".", "_").replace("-", "_")
    setattr(TestCorpus, f"test_{_safe}", _make_corpus_test(_case_name))


if __name__ == "__main__":
    unittest.main()
