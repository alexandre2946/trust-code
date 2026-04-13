""" 
Script to convert all XD lines to use the new 'OPT', 'REQ' or 'BRACE', 'NO_BRACE' tags.
Thanks AI.
"""

import sys
import re

# Replacement maps
ATTR_MAP = {
    "0": "REQ",
    "1": "OPT",
}

BRACE_MAP = {
    "-1": "INHERITS_BRACE",
    "0": "NO_BRACE",
    "1": "BRACE",
}

COMMA_MAP = {
    "-1": "INHERITS_COMMA",
    "0": "NO_COMMA",
    "1": "COMMA",
}

def transform_line(line: str) -> str:

    # --- XD attr ---
    # // XD attr word1 word2 word3 an_int ...
    line = re.sub(
        r"(//\s*XD\s+attr\s+\S+\s+\S+\s+\S+\s+)(\d+)",
        lambda m: m.group(1) + ATTR_MAP.get(m.group(2), m.group(2)),
        line,
    )

    # --- XD listobj (handle both integers in one pass) ---
    line = re.sub(
        r"(//\s*XD\s+\S+\s+listobj\s+\S+\s+)(-?\d+)(.*?\s)(-?\d+)",
        lambda m: (
            m.group(1)
            + BRACE_MAP.get(m.group(2), m.group(2))
            + m.group(3)
            + COMMA_MAP.get(m.group(4), m.group(4))
        ),
        line,
    )

    # --- XD general (non-listobj, non-attr) ---
    line = re.sub(
        r"(//\s*XD\s+(?!attr\b)(?!\S+\s+listobj\b)\S+\s+\S+\s+\S+\s+)(-?\d+)",
        lambda m: m.group(1) + BRACE_MAP.get(m.group(2), m.group(2)),
        line,
    )

    return line


def process_file(filepath: str):
    with open(filepath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    updated_lines = [transform_line(line) for line in lines]

    with open(filepath, "w", encoding="utf-8") as f:
        f.writelines(updated_lines)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        sys.exit(1)

    process_file(sys.argv[1])
    
