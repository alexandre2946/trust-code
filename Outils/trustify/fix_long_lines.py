import sys
import os
import tempfile

MAX_LEN = 260
CONT_PREFIX = "// XD_CONT "

def split_line_wordwise(line, max_len=MAX_LEN):
    words = line.split()
    chunks = []
    current = ""

    for word in words:
        if len(current) + (1 if current else 0) + len(word) > max_len:
            if current:
                chunks.append(current)
                current = word
            else:
                # Word longer than max_len → force split
                chunks.append(word[:max_len])
                remainder = word[max_len:]
                while len(remainder) > max_len:
                    chunks.append(remainder[:max_len])
                    remainder = remainder[max_len:]
                current = remainder
        else:
            current = word if not current else current + " " + word

    if current:
        chunks.append(current)

    return chunks


def process_file_inplace(file_path):
    dir_name = os.path.dirname(file_path)

    with tempfile.NamedTemporaryFile("w", delete=False, dir=dir_name, encoding="utf-8") as tmpfile:
        temp_name = tmpfile.name

        with open(file_path, "r", encoding="utf-8") as infile:
            for line in infile:
                stripped_line = line.rstrip("\n")

                if stripped_line.lstrip().startswith("// XD") and len(stripped_line) > MAX_LEN:
                    chunks = split_line_wordwise(stripped_line)

                    # First line stays unchanged
                    tmpfile.write(chunks[0] + "\n")

                    # Continuation lines get the prefix
                    for chunk in chunks[1:]:
                        tmpfile.write(CONT_PREFIX + chunk + "\n")
                else:
                    tmpfile.write(line)

    os.replace(temp_name, file_path)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py input.txt")
        sys.exit(1)

    process_file_inplace(sys.argv[1])

