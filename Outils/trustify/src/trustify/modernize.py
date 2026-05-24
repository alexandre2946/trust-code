"""XD-tag modernization engine.

Public entry point: `modernize(projects, trust_root, apply)`.

Five rules are applied per file, in order:

1. Numeric brace flag on XD block headers -> named token
   (BRACE / NO_BRACE / INHERITS_BRACE).
2. Numeric opt flag on XD attr lines -> REQ / OPT.
3. XD lines > 120 chars (block headers, XD attr) split at word
   boundary into base + XD_CONT continuation(s).
4. XD_ADD_P (any level) canonicalised unconditionally to a
   metadata-only opener line plus a (possibly empty) XD_CONT line
   carrying the description.
5. Multi-space whitespace between structural tokens in the XD
   comment collapses to a single space (the C++ prefix is left
   verbatim). Folded into rules 3 and 4 as a side effect of
   rebuilding the line via `" ".join(...)`.

Existing XD_CONT-broken blocks are folded first so the output is
idempotent: running modernize on its own output is a no-op.
"""

from __future__ import annotations

import contextlib
import difflib
import glob
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

LINE_LENGTH_LIMIT = 120

_BRACE_MAP = {
    "-3": "BRACE",
    "-2": "NO_BRACE",
    "-1": "INHERITS_BRACE",
    "0": "NO_BRACE",
    "1": "BRACE",
}
_OPT_MAP = {"0": "REQ", "1": "OPT"}
_NAMED_BRACE = {"BRACE", "NO_BRACE", "INHERITS_BRACE"}
_NAMED_OPT = {"REQ", "OPT"}

# Match an XD opener and capture the level prefix + the tag suffix
# (so we can dispatch on the structured form: XD vs XD attr vs XD_ADD_P
# vs XD_ADD_DICO). XD_CONT lines are detected by a separate regex
# (also captured here so the per-file processor can fold continuations).
# Note: `\s*` (zero spaces accepted) so `//XD foo` is classified too —
# `apply_split_rule` then re-emits the canonical `// XD foo` form via
# the `" ".join(seg_tokens)` rebuild (audit 1.8).
_OPENER_RE = re.compile(r"//\s*([23]?)XD(?:_(ADD_P|ADD_DICO))?\s")
_CONT_RE = re.compile(r"^\s*//\s*([23]?)XD_CONT(?:\s|$)")


@dataclass
class PerFileCounts:
    brace: int = 0
    opt: int = 0
    splits: int = 0
    unsplittable: int = 0


@dataclass
class FileReport:
    path: str
    new_content: str
    diff: str  # empty when no change
    counts: PerFileCounts


@dataclass
class ModernizeReport:
    files: list[FileReport] = field(default_factory=list)
    summary: dict[str, int] = field(default_factory=lambda: {"brace": 0, "opt": 0, "splits": 0, "unsplittable": 0})

    @property
    def changed_files(self) -> list[FileReport]:
        return [f for f in self.files if f.diff]


@dataclass
class LineMeta:
    lvl: str  # "" / "2" / "3"
    kind: str  # "block" / "attr" / "add_p" / "add_dico"
    indent: str  # whitespace before "//"


def classify_line(line: str) -> LineMeta | None:
    """Classify a source line. Returns None for non-XD lines and for
    XD_CONT continuation lines (continuations are handled by the
    per-file processor, not by individual rules).
    """
    m = _OPENER_RE.search(line)
    if m is None:
        return None
    lvl = m.group(1)
    suffix = m.group(2)
    indent_match = re.match(r"^(\s*)", line)
    indent = indent_match.group(1) if indent_match else ""
    if suffix == "ADD_P":
        return LineMeta(lvl=lvl, kind="add_p", indent=indent)
    if suffix == "ADD_DICO":
        return LineMeta(lvl=lvl, kind="add_dico", indent=indent)
    # Bare XD/2XD/3XD: dispatch on the next token (attr vs block).
    after = line[m.end() :].lstrip().split(None, 1)
    if after and after[0] == "attr":
        return LineMeta(lvl=lvl, kind="attr", indent=indent)
    return LineMeta(lvl=lvl, kind="block", indent=indent)


def is_cont(line: str) -> str | None:
    """Return the level prefix ("" / "2" / "3") if `line` is an XD_CONT
    continuation, else None."""
    m = _CONT_RE.match(line)
    return m.group(1) if m else None


def _split_at_tag(line: str) -> tuple[str, list[str]]:
    """Split `line` into `(cpp_prefix, seg_tokens)`.

    - `cpp_prefix` is everything BEFORE the `//` of the XD opener. May
      be just leading whitespace for pure-comment lines, or include C++
      code for `XD_ADD_P` / `XD_ADD_DICO` lines (and for the rarer case
      of `XD attr` lines emitted from within a `param.ajouter*` call).
    - `seg_tokens` is the whitespace-split list starting at `//`.
      `seg_tokens[0]` is always `//`.

    Returns `("", [])` when the line carries no opener (defensive —
    callers already classified it).
    """
    m = _OPENER_RE.search(line)
    if m is None:
        return "", []
    tag_start = m.start()
    # Normalise the `//` prefix to a single trailing space BEFORE the
    # whitespace split, so a `//XD ...` input produces the same
    # `["//", "XD", ...]` token list as `// XD ...` does. Without this,
    # `.split()` would emit `["//XD", ...]` as a single token, shift
    # every downstream index by one, and silently skip the brace / opt
    # rewrites (audit 1.8). Also collapses any multi-space form into
    # the canonical `// ` so the rebuild emits `// XD ...`.
    body = re.sub(r"^//\s*", "// ", line[tag_start:])
    return line[:tag_start], body.split()


def apply_brace_rule(line: str, meta: LineMeta) -> tuple[str, int]:
    """Rule 1: rewrite numeric brace flag -> named token on block headers."""
    if meta.kind != "block":
        return line, 0
    cpp_prefix, seg_tokens = _split_at_tag(line)
    # seg_tokens layout:
    #   ["//", "<lvl>XD", "<name>", "<name_base>", "<syno>", "<brace_flag>", ...desc...]
    if len(seg_tokens) < 6:
        return line, 0
    flag = seg_tokens[5]
    if flag.upper() in _NAMED_BRACE:
        return line, 0
    mapped = _BRACE_MAP.get(flag)
    if mapped is None:
        return line, 0
    seg_tokens[5] = mapped
    return cpp_prefix + " ".join(seg_tokens), 1


def apply_opt_rule(line: str, meta: LineMeta) -> tuple[str, int]:
    """Rule 2: rewrite numeric opt flag -> REQ/OPT on attr lines."""
    if meta.kind != "attr":
        return line, 0
    cpp_prefix, seg_tokens = _split_at_tag(line)
    # seg_tokens layout:
    #   ["//", "<lvl>XD", "attr", "<name>", "<type>", "<syno>", "<opt>", ...desc...]
    if len(seg_tokens) < 7:
        return line, 0
    flag = seg_tokens[6]
    if flag.upper() in _NAMED_OPT:
        return line, 0
    mapped = _OPT_MAP.get(flag)
    if mapped is None:
        return line, 0
    seg_tokens[6] = mapped
    return cpp_prefix + " ".join(seg_tokens), 1


# Number of structured (non-description) tokens before the free-text
# description, per tag kind. The tokens are:
#   - block:    ["//", "<lvl>XD", "<name>", "<name_base>", "<syno>", "<flag>"]
#   - attr:     ["//", "<lvl>XD", "attr", "<name>", "<type>", "<syno>", "<opt>"]
#   - add_p:    ["//", "<lvl>XD_ADD_P", "<type>"]
#   - add_dico: ["//", "<lvl>XD_ADD_DICO", "<name>"]
# The description occupies positions [_STRUCTURED_TOKENS[kind]:].
_STRUCTURED_TOKENS = {"block": 6, "attr": 7, "add_p": 3, "add_dico": 3}


def apply_split_rule(line: str, meta: LineMeta) -> tuple[list[str], int]:
    """Rule 3 / Rule 4 / Rule 5: split, canonicalise, or collapse.

    Rule 5 (always, all kinds): consecutive whitespace inside the XD
    comment portion collapses to a single space. The C++ prefix (if
    any) is preserved verbatim — only the `// ...` payload is
    normalised. Implemented as a side effect of rebuilding the line
    via `cpp_prefix + " ".join(seg_tokens)`.

    Rule 4 (XD_ADD_P, any level): unconditionally rewrite into a
    metadata-only opener line plus a (possibly empty) XD_CONT line
    carrying the description. Never reports unsplittable — the C++
    prefix may legitimately push the opener past 120 chars.

    Rule 3 (all other XD kinds): greedy word-packer that splits lines
    > 120 chars at word boundaries. The OPENER tag and its structured
    tokens MUST stay on the first line — only the description can
    move into XD_CONT continuation lines. If the opener+structured
    prefix alone exceeds 120 chars, the line is reported as
    unsplittable.

    Returns `(output_lines, n_unsplittable)`. Callers compute the "net
    new splits" count from the output-line count delta — an idempotent
    re-run reports zero splits.
    """
    if meta.kind == "add_p":
        return _canonicalize_add_p(line, meta), 0
    cpp_prefix, seg_tokens = _split_at_tag(line)
    if not seg_tokens:
        return [line], 0
    # Rule 5: rebuild the XD comment with single-space-collapsed tokens.
    normalized = cpp_prefix + " ".join(seg_tokens)
    if len(normalized) <= LINE_LENGTH_LIMIT:
        return [normalized], 0
    prefix_count = _STRUCTURED_TOKENS[meta.kind]
    if len(seg_tokens) <= prefix_count:
        # No description text after the structured part — nothing to
        # split off. Pathological; emit the normalised line unchanged.
        return [normalized], 1
    structured = " ".join(seg_tokens[:prefix_count])
    desc_words = seg_tokens[prefix_count:]
    # The first line is `<cpp_prefix><structured>` plus a greedy
    # tail of description words. If even this base exceeds the limit
    # there is no place to put the opener — unsplittable.
    first_base = cpp_prefix + structured
    if len(first_base) > LINE_LENGTH_LIMIT:
        return [normalized], 1
    cont_prefix = f"{meta.indent}// {meta.lvl}XD_CONT"
    output: list[str] = []
    current = first_base
    on_cont = False
    for word in desc_words:
        candidate = f"{current} {word}"
        if len(candidate) <= LINE_LENGTH_LIMIT:
            current = candidate
            continue
        output.append(current)
        current = f"{cont_prefix} {word}"
        on_cont = True
        if len(current) > LINE_LENGTH_LIMIT:
            # Single description word overflows even a fresh
            # continuation. Emit the normalised line unchanged rather
            # than force-break a long identifier.
            return [normalized], 1
    output.append(current)
    if not on_cont:
        return [normalized], 0
    return output, 0


def _canonicalize_add_p(line: str, meta: LineMeta) -> list[str]:
    """Canonical form for XD_ADD_P: metadata-only opener + XD_CONT desc.

    Always emits at least one XD_CONT line — bare placeholder when the
    description is empty, otherwise one or more packed at 120 chars.
    Returns the input unchanged if structurally malformed (defensive).
    """
    cpp_prefix, seg_tokens = _split_at_tag(line)
    if not seg_tokens:
        return [line]
    prefix_count = _STRUCTURED_TOKENS["add_p"]
    if len(seg_tokens) < prefix_count:
        return [line]
    structured = " ".join(seg_tokens[:prefix_count])
    desc_words = seg_tokens[prefix_count:]
    opener = cpp_prefix + structured
    cont_prefix = f"{meta.indent}// {meta.lvl}XD_CONT"
    if not desc_words:
        return [opener, cont_prefix]
    output = [opener]
    current = cont_prefix
    for word in desc_words:
        candidate = f"{current} {word}"
        if len(candidate) <= LINE_LENGTH_LIMIT:
            current = candidate
            continue
        output.append(current)
        current = f"{cont_prefix} {word}"
    output.append(current)
    return output


def is_unsplittable(line: str) -> bool:
    """Return True iff modernize cannot fit the post-rewrite opener
    line within LINE_LENGTH_LIMIT.

    Used by the scanner to suppress the "long XD line" warning on
    lines that no amount of `XD_CONT` splitting can shorten — the
    structured prefix (opener tokens that MUST stay on line 1)
    already exceeds the limit.

    Returns False for non-XD lines (defensive — the warning gate
    already filters those, but we don't want a surprising True here).
    """
    meta = classify_line(line)
    if meta is None:
        return False
    cpp_prefix, seg_tokens = _split_at_tag(line)
    if not seg_tokens:
        return False
    prefix_count = _STRUCTURED_TOKENS[meta.kind]
    if len(seg_tokens) < prefix_count:
        return False
    structured = " ".join(seg_tokens[:prefix_count])
    return len(cpp_prefix + structured) > LINE_LENGTH_LIMIT


def _fold_block(opener: str, conts: list[str], lvl: str) -> str:
    """Fold an opener + its XD_CONT continuations into one logical line.

    Mirrors the scanner's `curr_block[1] = " ".join([curr_block[1], lin])`
    behaviour: strip the `// <lvl>XD_CONT ` prefix from each continuation
    line and concatenate with single-space separators.
    """
    if not conts:
        return opener.rstrip("\n")
    parts = [opener.rstrip("\n").rstrip()]
    for c in conts:
        c = c.rstrip("\n").lstrip()
        m = re.match(rf"^//\s*{lvl}XD_CONT\s*(.*)$", c)
        if m is None:
            # Defensive: pass through (shouldn't happen — we pre-filtered
            # continuations via is_cont upstream).
            parts.append(c)
        else:
            parts.append(m.group(1))
    return " ".join(p for p in parts if p)


def modernize_file_content(content: str) -> tuple[str, PerFileCounts]:
    """Apply all three rules to `content`. Pure function — no file IO.

    Preserves non-XD lines verbatim and keeps the trailing newline (or
    its absence) intact.
    """
    counts = PerFileCounts()
    lines = content.splitlines(keepends=False)
    trailing_nl = content.endswith("\n") if content else False
    out_lines: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        meta = classify_line(line)
        if meta is None:
            # Non-XD line OR orphan XD_CONT — pass through unchanged.
            out_lines.append(line)
            i += 1
            continue
        # Found an opener. Gather following continuations at the same level.
        cont_lines: list[str] = []
        j = i + 1
        while j < len(lines):
            cont_lvl = is_cont(lines[j])
            if cont_lvl is not None and cont_lvl == meta.lvl:
                cont_lines.append(lines[j])
                j += 1
            else:
                break
        folded = _fold_block(line, cont_lines, meta.lvl)
        folded, n_brace = apply_brace_rule(folded, meta)
        folded, n_opt = apply_opt_rule(folded, meta)
        split_lines, n_unsplittable = apply_split_rule(folded, meta)
        # Count "net new" splits: how many MORE lines the output has
        # vs the input block. An idempotent re-run (same input lines,
        # same output lines) reports zero splits.
        n_splits = max(0, len(split_lines) - (1 + len(cont_lines)))
        counts.brace += n_brace
        counts.opt += n_opt
        counts.splits += n_splits
        counts.unsplittable += n_unsplittable
        out_lines.extend(split_lines)
        i = j
    joined = "\n".join(out_lines)
    if trailing_nl:
        joined += "\n"
    return joined, counts


def _resolve_src_dirs(projects, trust_root):
    """Map projects + trust_root through the standard resolution to a
    list of `<project>/src` directories (same input as scanSourceFiles).

    Unlike the schema-building entry points, modernize does NOT apply the
    `$project_directory` / `$TRUST_ROOT` env fallbacks here. Its scope is
    decided entirely upstream by the CLI's `_scope_modernize`, which
    reads those env vars explicitly and picks exactly one target. Re-
    applying `effective_trust_root` would silently re-add the whole TRUST
    tree whenever `$TRUST_ROOT` happens to be set in the shell — widening
    an intended baltik-only rewrite (e.g. modernizing a baltik that lives
    inside TRUST) into a rewrite of all of TRUST. `resolve_projects`
    raises if nothing was passed, so a caller that supplies neither a
    project nor a trust_root gets a clear error rather than an implicit
    env-driven scope.
    """
    from trustify.projects import Project, resolve_projects, validate_projects

    resolved: list[Project] = resolve_projects(projects, trust_root)
    validate_projects(resolved)
    return [str(p.path / "src") for p in resolved]


def _walk_xd_files(src_dirs: list[str]) -> list[str]:
    """Return the sorted, deduped list of .cpp + .xd files under src_dirs.

    Symlinked files are skipped: `glob` with `recursive=True` returns
    symlinked file leaves, and modernize used to write through them —
    a baltik that symlinks shared `.xd` files from `$TRUST_ROOT` would
    mutate `$TRUST_ROOT` under its feet (audit 2.4). The skip is silent;
    real-world TRUST sources don't symlink XD-bearing files.
    """
    seen: dict[str, str] = {}
    for d in src_dirs:
        for ext in ("cpp", "xd"):
            pattern = os.path.join(d, "**", f"*.{ext}")
            for p in sorted(glob.glob(pattern, recursive=True)):
                if os.path.islink(p):
                    continue
                seen[os.path.basename(p)] = p  # last project wins on basename collision
    return sorted(seen.values())


def _atomic_write_text(path: str, content: str) -> None:
    """Write `content` to `path` atomically: stage in a sibling
    `.<name>.tmp.<pid>` file, then `os.replace` onto the final path.
    A SIGINT, disk-full, or any other failure mid-write leaves `path`
    untouched at its previous content (audit 2.3) — without this, the
    `Path.write_text` truncate-then-write sequence could corrupt the
    source on any failure between truncate and full-write.

    When `path` is a symbolic link, resolve it via `os.path.realpath`
    and write to the TARGET. The link itself stays intact and the
    target gets the new content. Without this, `os.replace(tmp,
    symlink)` would silently overwrite the link with a regular file
    and the original target would keep its old bytes — surprising,
    and a footgun for any baltik that symlinks shared `.data`
    fixtures from `$TRUST_ROOT`. Matches vim / emacs /
    `sed --follow-symlinks` semantics.

    On any exception the temp file is cleaned up.
    """
    p = Path(path)
    if p.is_symlink():
        # Resolve so the temp file lands in the TARGET's parent dir
        # (atomic rename requires same filesystem). The symlink at
        # the original path is then untouched — readers through it
        # see the new content via the now-updated target.
        p = Path(os.path.realpath(p))
    tmp = p.parent / f".{p.name}.tmp.{os.getpid()}"
    try:
        tmp.write_text(content, encoding="utf-8")
        os.replace(tmp, p)
    except BaseException:
        with contextlib.suppress(OSError):
            tmp.unlink(missing_ok=True)
        raise


def modernize(
    projects: list[str] | None = None,
    trust_root: str | None = None,
    apply: bool = False,
) -> ModernizeReport:
    """Modernize XD tags in the resolved workspace.

    Dry-run by default; pass `apply=True` to write changes in place.
    """
    src_dirs = _resolve_src_dirs(projects, trust_root)
    files = _walk_xd_files(src_dirs)
    report = ModernizeReport()
    for path in files:
        with open(path, encoding="utf-8") as f:
            content = f.read()
        new_content, counts = modernize_file_content(content)
        diff = ""
        if new_content != content:
            diff = "".join(
                difflib.unified_diff(
                    content.splitlines(keepends=True),
                    new_content.splitlines(keepends=True),
                    fromfile=path,
                    tofile=path,
                )
            )
            if apply:
                _atomic_write_text(path, new_content)
        report.files.append(FileReport(path=path, new_content=new_content, diff=diff, counts=counts))
        report.summary["brace"] += counts.brace
        report.summary["opt"] += counts.opt
        report.summary["splits"] += counts.splits
        report.summary["unsplittable"] += counts.unsplittable
    return report
