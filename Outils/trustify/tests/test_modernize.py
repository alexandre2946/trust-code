"""Unit tests for trustify.modernize — the XD-tag rewrite engine."""

import os
import unittest
from pathlib import Path


class TestBraceRule(unittest.TestCase):
    """Rule 1: numeric brace flag on XD block headers -> named token."""

    def test_each_legacy_brace_int_maps_to_named_token(self):
        from trustify.modernize import modernize_file_content

        src = (
            "// XD a objet_u a -3 d1\n"
            "// XD b objet_u b -2 d2\n"
            "// XD c objet_u c -1 d3\n"
            "// XD d objet_u d 0 d4\n"
            "// XD e objet_u e 1 d5\n"
        )
        out, counts = modernize_file_content(src)
        expected = (
            "// XD a objet_u a BRACE d1\n"
            "// XD b objet_u b NO_BRACE d2\n"
            "// XD c objet_u c INHERITS_BRACE d3\n"
            "// XD d objet_u d NO_BRACE d4\n"
            "// XD e objet_u e BRACE d5\n"
        )
        self.assertEqual(out, expected)
        self.assertEqual(counts.brace, 5)
        self.assertEqual(counts.opt, 0)
        self.assertEqual(counts.splits, 0)
        self.assertEqual(counts.unsplittable, 0)

    def test_already_named_brace_flag_is_idempotent(self):
        from trustify.modernize import modernize_file_content

        src = "// XD a objet_u a BRACE d1\n// XD b objet_u b NO_BRACE d2\n// XD c objet_u c INHERITS_BRACE d3\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(out, src)
        self.assertEqual(counts.brace, 0)

    def test_brace_rule_does_not_touch_attr_or_add_p(self):
        from trustify.modernize import modernize_file_content

        # XD attr has a 5th-position opt flag, not a 4th-position brace
        # flag — the brace rule must not pick up the `0` in the attr
        # line and mistake it for a brace flag.
        src = (
            "// XD attr nam typ syn 0 description\n"
            'param.ajouter_flag("P0", &x); // XD_ADD_P rien Pressure\n'
            'param.dictionnaire("foo"); // XD_ADD_DICO foo\n'
        )
        _out, counts = modernize_file_content(src)
        self.assertEqual(counts.brace, 0)


class TestOptRule(unittest.TestCase):
    """Rule 2: numeric opt flag on XD attr lines -> REQ/OPT."""

    def test_zero_maps_to_REQ_one_maps_to_OPT(self):
        from trustify.modernize import modernize_file_content

        src = "// XD attr a typ a 0 d1\n// XD attr b typ b 1 d2\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(
            out,
            "// XD attr a typ a REQ d1\n// XD attr b typ b OPT d2\n",
        )
        self.assertEqual(counts.opt, 2)
        self.assertEqual(counts.brace, 0)

    def test_already_named_opt_flag_is_idempotent(self):
        from trustify.modernize import modernize_file_content

        src = "// XD attr a typ a REQ d1\n// XD attr b typ b OPT d2\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(out, src)
        self.assertEqual(counts.opt, 0)

    def test_opt_rule_does_not_touch_block_header(self):
        from trustify.modernize import modernize_file_content

        # The `0` at position 4 of a block header is a brace flag, not an
        # opt flag. The opt rule must leave it alone (the brace rule will
        # handle it).
        src = "// XD a objet_u a 0 desc\n"
        _out, counts = modernize_file_content(src)
        self.assertEqual(counts.opt, 0)


class TestSplitRule(unittest.TestCase):
    """Rule 3: lines > 120 chars split at word boundary with XD_CONT."""

    def test_short_line_untouched(self):
        from trustify.modernize import modernize_file_content

        src = "// XD a objet_u a BRACE short desc\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(out, src)
        self.assertEqual(counts.splits, 0)

    def test_line_over_120_splits_at_word_boundary(self):
        from trustify.modernize import modernize_file_content

        word = "x" * 5
        words = " ".join([word] * 25)
        src = f"// XD a objet_u a BRACE {words}\n"
        self.assertGreater(len(src.rstrip()), 120)
        out, counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        for line in out_lines:
            self.assertLessEqual(len(line), 120, f"line over 120: {line!r}")
        self.assertTrue(out_lines[0].startswith("// XD a objet_u a BRACE"))
        for line in out_lines[1:]:
            self.assertTrue(line.lstrip().startswith("// XD_CONT "), line)
        self.assertGreaterEqual(counts.splits, 1)
        self.assertEqual(counts.unsplittable, 0)

    def test_level_2_uses_2XD_CONT(self):
        from trustify.modernize import modernize_file_content

        word = "y" * 5
        words = " ".join([word] * 25)
        src = f"// 2XD attr a typ a OPT {words}\n"
        out, _counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        for line in out_lines[1:]:
            self.assertTrue(line.lstrip().startswith("// 2XD_CONT "), line)

    def test_xd_add_p_tag_stays_on_first_line_with_cpp_prefix(self):
        # Reproduce a bug: a long `// XD_ADD_P` line with C++ code
        # before the comment had its XD_ADD_P tag tokenized as part of
        # the description and pushed into the continuation as
        # `// XD_CONT XD_ADD_P ...`. The opener MUST stay on line 0.
        from trustify.modernize import modernize_file_content

        src = (
            '  param.ajouter("LISGRP", &lis_face_grp_); // XD_ADD_P '
            "listchaine List of face groups to read and register in the MED file and to use later.\n"
        )
        self.assertGreater(len(src.rstrip()), 120)
        out, counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertIn("XD_ADD_P", out_lines[0])
        for line in out_lines[1:]:
            self.assertNotIn("XD_ADD_P", line, f"opener tag leaked into continuation: {line!r}")
            self.assertTrue(line.lstrip().startswith("// XD_CONT "), line)
        # Sanity: still got at least one continuation.
        self.assertGreater(len(out_lines), 1)
        self.assertGreaterEqual(counts.splits, 1)

    def test_xd_attr_opt_flag_after_cpp_prefix(self):
        # `// XD attr` lines can also have C++ code before the comment.
        # The opt-flag rule must find the flag at the right offset
        # relative to the `//` — not the start of the line.
        from trustify.modernize import modernize_file_content

        src = '  param.ajouter("nam", &x); // XD attr nam typ syn 0 a short desc\n'
        out, counts = modernize_file_content(src)
        # The `0` must have been rewritten to `REQ`.
        self.assertIn("REQ", out)
        self.assertNotIn(" 0 ", out.split("//", 1)[1])
        self.assertEqual(counts.opt, 1)

    def test_unsplittable_single_word_is_reported(self):
        from trustify.modernize import modernize_file_content

        long_word = "z" * 150
        src = f"// XD a objet_u a BRACE {long_word}\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(out, src)
        self.assertEqual(counts.splits, 0)
        self.assertEqual(counts.unsplittable, 1)


class TestAddPCanonicalForm(unittest.TestCase):
    """Rule 4: XD_ADD_P always canonicalises to a metadata-only opener
    line plus a (possibly empty) XD_CONT-anchored description, regardless
    of total line length. Idempotent."""

    def test_short_add_p_with_description_still_splits(self):
        from trustify.modernize import modernize_file_content

        src = 'param.ajouter("k", &x); // XD_ADD_P int short desc\n'
        self.assertLess(len(src.rstrip()), 120)
        out, counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertEqual(len(out_lines), 2)
        self.assertTrue(out_lines[0].endswith("// XD_ADD_P int"))
        self.assertEqual(out_lines[1].lstrip(), "// XD_CONT short desc")
        self.assertEqual(counts.splits, 1)
        self.assertEqual(counts.unsplittable, 0)

    def test_add_p_with_empty_description_emits_bare_xd_cont(self):
        from trustify.modernize import modernize_file_content

        src = 'param.ajouter_flag("k", &x); // XD_ADD_P rien\n'
        out, counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertEqual(len(out_lines), 2)
        self.assertTrue(out_lines[0].endswith("// XD_ADD_P rien"))
        self.assertEqual(out_lines[1].lstrip(), "// XD_CONT")
        self.assertEqual(counts.splits, 1)
        self.assertEqual(counts.unsplittable, 0)

    def test_long_cpp_prefix_does_not_count_as_unsplittable(self):
        from trustify.modernize import modernize_file_content

        long_arg = '"x"' * 50
        src = f"param.ajouter({long_arg}, &x); // XD_ADD_P int description here\n"
        self.assertGreater(len(src.rstrip()), 120)
        out, counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertEqual(len(out_lines), 2)
        self.assertIn("// XD_ADD_P int", out_lines[0])
        self.assertEqual(out_lines[1].lstrip(), "// XD_CONT description here")
        self.assertEqual(counts.unsplittable, 0)

    def test_level_2_and_3_add_p_canonicalize(self):
        from trustify.modernize import modernize_file_content

        src = (
            'param.ajouter("k", &x); // 2XD_ADD_P int level two desc\n'
            'param.ajouter("k", &y); // 3XD_ADD_P int level three desc\n'
        )
        out, _ = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertEqual(len(out_lines), 4)
        self.assertIn("// 2XD_ADD_P int", out_lines[0])
        self.assertEqual(out_lines[1].lstrip(), "// 2XD_CONT level two desc")
        self.assertIn("// 3XD_ADD_P int", out_lines[2])
        self.assertEqual(out_lines[3].lstrip(), "// 3XD_CONT level three desc")

    def test_existing_xd_cont_continuations_are_re_canonicalised(self):
        from trustify.modernize import modernize_file_content

        src = 'param.ajouter("k", &x); // XD_ADD_P int short desc\n// XD_CONT continued here\n'
        out, _ = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertIn("// XD_ADD_P int", out_lines[0])
        self.assertNotIn("short", out_lines[0])
        self.assertNotIn("continued", out_lines[0])
        for line in out_lines[1:]:
            self.assertTrue(line.lstrip().startswith("// XD_CONT"), line)

    def test_canonical_form_is_idempotent(self):
        from trustify.modernize import modernize_file_content

        srcs = [
            'param.ajouter("k", &x); // XD_ADD_P int short desc\n',
            'param.ajouter_flag("k", &x); // XD_ADD_P rien\n',
            'param.ajouter("k", &x); // 2XD_ADD_P int level two\n',
        ]
        for src in srcs:
            out_once, _ = modernize_file_content(src)
            out_twice, counts_twice = modernize_file_content(out_once)
            self.assertEqual(out_once, out_twice, f"non-idempotent for: {src!r}")
            self.assertEqual(counts_twice.splits, 0, f"net splits on re-run for: {src!r}")
            self.assertEqual(counts_twice.unsplittable, 0)

    def test_indentation_of_xd_cont_matches_cpp_statement(self):
        from trustify.modernize import modernize_file_content

        src = '    param.ajouter("k", &x); // XD_ADD_P int description\n'
        out, _ = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertTrue(out_lines[1].startswith("    // XD_CONT "), out_lines[1])


class TestWhitespaceCollapse(unittest.TestCase):
    """Rule 5: multi-space whitespace inside the XD comment portion
    collapses to single spaces. The C++ prefix (if any) is preserved
    verbatim — only the `// ...` payload is normalised."""

    def test_short_xd_attr_with_multi_spaces_collapses(self):
        from trustify.modernize import modernize_file_content

        src = "// XD   attr name typ name OPT short desc\n"
        out, _ = modernize_file_content(src)
        self.assertEqual(out, "// XD attr name typ name OPT short desc\n")

    def test_block_header_with_multi_spaces_collapses(self):
        from trustify.modernize import modernize_file_content

        src = "// XD  block  objet_u  block  BRACE  short desc\n"
        out, _ = modernize_file_content(src)
        self.assertEqual(out, "// XD block objet_u block BRACE short desc\n")

    def test_long_line_collapses_then_fits_under_120(self):
        from trustify.modernize import modernize_file_content

        # Real-world Solv_Optimal.cpp shape: extra whitespace alone
        # pushes the line over the limit. After collapse, fits in 120
        # without needing XD_CONT.
        src = (
            "// XD   attr pas_de_solution_initiale rien pas_de_solution_initiale OPT"
            " Resolution isn't initialized with the solution xx\n"
        )
        self.assertGreater(len(src.rstrip()), 120)
        out, counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertEqual(len(out_lines), 1, "should not need XD_CONT after collapse")
        self.assertLessEqual(len(out_lines[0]), 120)
        self.assertNotIn("XD   attr", out)
        self.assertEqual(counts.unsplittable, 0)

    def test_long_line_still_needs_split_after_collapse(self):
        from trustify.modernize import modernize_file_content

        # Even after collapse the line is > 120 -> split rule kicks in,
        # description moves to XD_CONT, all output lines are single-spaced.
        long_desc = " ".join(["word"] * 30)
        src = f"// XD  a  objet_u  a  BRACE  {long_desc}\n"
        out, _ = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertGreater(len(out_lines), 1)
        for line in out_lines:
            self.assertNotIn("  ", line.lstrip(), f"double space leaked: {line!r}")
            self.assertLessEqual(len(line), 120)

    def test_no_change_when_already_single_spaced(self):
        from trustify.modernize import modernize_file_content

        src = "// XD attr name typ name OPT short desc\n"
        out, _ = modernize_file_content(src)
        self.assertEqual(out, src)

    def test_no_space_after_slash_is_recognized_and_normalized(self):
        """Audit 1.8: `_OPENER_RE` used to require `\\s+` between `//`
        and `XD`, so a `//XD foo ...` line slipped through modernize
        entirely (never classified, never warned about). Allow zero
        spaces too — and emit the canonical single-space form on
        output so the in-tree convention stays uniform.
        """
        from trustify.modernize import modernize_file_content

        src = "//XD legacy_block objet_u legacy_block 0 a block\n"
        out, counts = modernize_file_content(src)
        # The line is recognised — its numeric brace flag must have
        # been rewritten.
        self.assertEqual(counts.brace, 1)
        # And the output is canonicalised to `// XD ...` (single space).
        self.assertEqual(out, "// XD legacy_block objet_u legacy_block NO_BRACE a block\n")

    def test_no_space_canonical_form_is_idempotent(self):
        # Running modernize on its own output is a no-op (the rebuild
        # already produced the canonical `// XD ...` form).
        from trustify.modernize import modernize_file_content

        src = "//XD legacy_block objet_u legacy_block 0 desc\n"
        once, _ = modernize_file_content(src)
        twice, counts = modernize_file_content(once)
        self.assertEqual(once, twice)
        self.assertEqual(counts.brace, 0)
        self.assertEqual(counts.splits, 0)

    def test_no_space_XD_CONT_also_recognized(self):
        # Same fix applies to XD_CONT — `//XD_CONT ...` should be
        # folded into its opener (and re-emitted with the canonical
        # single space).
        from trustify.modernize import modernize_file_content

        long_word = "x" * 100
        src = f"// XD legacy_block objet_u legacy_block 0 short\n//XD_CONT and a very long continuation {long_word}\n"
        out, _ = modernize_file_content(src)
        # The continuation must have been folded — the input two lines
        # become one logical block. The output must contain the
        # continuation's text and the canonical form.
        self.assertIn("// XD_CONT", out, "//XD_CONT was not recognised as a continuation")

    def test_no_space_works_for_level_2_block(self):
        from trustify.modernize import modernize_file_content

        src = "//2XD legacy_block objet_u legacy_block 0 a block\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(counts.brace, 1)
        self.assertEqual(out, "// 2XD legacy_block objet_u legacy_block NO_BRACE a block\n")

    def test_no_space_works_for_level_3_block(self):
        from trustify.modernize import modernize_file_content

        src = "//3XD legacy_block objet_u legacy_block 0 a block\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(counts.brace, 1)
        self.assertEqual(out, "// 3XD legacy_block objet_u legacy_block NO_BRACE a block\n")

    def test_no_space_works_for_level_2_attr(self):
        from trustify.modernize import modernize_file_content

        src = "//2XD attr legacy_opt int legacy_opt 1 an attribute\n"
        out, counts = modernize_file_content(src)
        self.assertEqual(counts.opt, 1)
        self.assertEqual(out, "// 2XD attr legacy_opt int legacy_opt OPT an attribute\n")

    def test_no_space_works_for_XD_ADD_P(self):
        # add_p is canonicalised to a metadata-only opener + XD_CONT
        # carrying the description. Confirm both lines come out with
        # the canonical `// ` prefix.
        from trustify.modernize import modernize_file_content

        src = "//XD_ADD_P int just a number\n"
        out, _ = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertEqual(out_lines[0], "// XD_ADD_P int")
        self.assertEqual(out_lines[1], "// XD_CONT just a number")

    def test_no_space_works_for_level_2_XD_ADD_P(self):
        from trustify.modernize import modernize_file_content

        src = "//2XD_ADD_P int just a number\n"
        out, _ = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertEqual(out_lines[0], "// 2XD_ADD_P int")
        self.assertEqual(out_lines[1], "// 2XD_CONT just a number")

    def test_no_space_works_for_XD_ADD_DICO(self):
        from trustify.modernize import modernize_file_content

        # XD_ADD_DICO has no rewrite rules of its own; the test
        # confirms classification + Rule-5 canonical-form rebuild.
        src = "//XD_ADD_DICO some_name\n"
        out, _ = modernize_file_content(src)
        self.assertEqual(out, "// XD_ADD_DICO some_name\n")

    def test_no_space_works_for_level_3_XD_CONT(self):
        # `_fold_block` uses a per-level regex internally — needs to
        # accept the no-space form for every level too.
        from trustify.modernize import modernize_file_content

        long_word = "x" * 100
        src = f"// 3XD legacy_block objet_u legacy_block 0 short\n//3XD_CONT a long continuation {long_word}\n"
        out, _ = modernize_file_content(src)
        self.assertIn("// 3XD_CONT", out, "//3XD_CONT was not recognised as a level-3 continuation")

    def test_cpp_prefix_whitespace_preserved(self):
        from trustify.modernize import modernize_file_content

        # Multi-spaces in the C++ prefix (e.g., between `);` and `//`)
        # are NOT touched — only the XD comment is normalised.
        src = '  param.ajouter("x", &y);   // XD attr  name  typ  name  OPT  desc\n'
        out, _ = modernize_file_content(src)
        # Triple-space between `;` and `//` is preserved.
        self.assertIn(";   //", out)
        # But `XD attr name typ name OPT desc` is single-spaced.
        self.assertIn("// XD attr name typ name OPT desc", out)

    def test_whitespace_collapse_is_idempotent(self):
        from trustify.modernize import modernize_file_content

        src = "// XD   attr  name  typ  name  OPT   desc\n"
        out_once, _ = modernize_file_content(src)
        out_twice, counts_twice = modernize_file_content(out_once)
        self.assertEqual(out_once, out_twice)
        self.assertEqual(counts_twice.splits, 0)


class TestIdempotence(unittest.TestCase):
    def test_existing_XD_CONT_block_folds_and_resplits_canonically(self):
        from trustify.modernize import modernize_file_content

        word = "w" * 5
        chunk_a = " ".join([word] * 8)
        chunk_b = " ".join([word] * 8)
        chunk_c = " ".join([word] * 9)
        src = f"// XD a objet_u a BRACE {chunk_a}\n// XD_CONT {chunk_b}\n// XD_CONT {chunk_c}\n"
        out_once, _ = modernize_file_content(src)
        out_twice, counts_twice = modernize_file_content(out_once)
        self.assertEqual(out_once, out_twice, "modernize should be idempotent")
        self.assertEqual(counts_twice.splits, 0)
        self.assertEqual(counts_twice.brace, 0)
        self.assertEqual(counts_twice.opt, 0)


class TestCombinedRules(unittest.TestCase):
    def test_brace_rename_runs_before_split_decision(self):
        from trustify.modernize import modernize_file_content

        # Build a line that is EXACTLY 120 chars with `1` and 124 chars
        # once `1` -> `BRACE` (+4 chars). Confirms renames happen first
        # and the split decision uses the post-rename length.
        prefix = "// XD a objet_u a 1 "  # 20 chars
        # Total target with `1`: 120 chars; description budget: 100 chars.
        # 16 words of "xxxxx" + 15 separator spaces = 95 chars. Pad with
        # " yyyy" (5 chars) to hit 100 exactly.
        desc = " ".join(["xxxxx"] * 16) + " yyyy"
        self.assertEqual(len(desc), 100)
        src = f"{prefix}{desc}\n"
        self.assertEqual(len(src.rstrip()), 120)
        out, counts = modernize_file_content(src)
        out_lines = out.rstrip("\n").split("\n")
        self.assertIn("BRACE", out_lines[0])
        self.assertNotIn(" 1 ", out_lines[0])
        for line in out_lines:
            self.assertLessEqual(len(line), 120)
        self.assertEqual(counts.brace, 1)
        self.assertGreaterEqual(counts.splits, 1)


class TestIsUnsplittable(unittest.TestCase):
    """The `is_unsplittable` predicate gates the scanner's long-line
    warning: True iff modernize cannot fit the post-rewrite opener
    line within 120 chars."""

    def test_xd_attr_with_long_chaine_into_type_is_unsplittable(self):
        from trustify.modernize import is_unsplittable

        values = ",".join(f'"value_{i}"' for i in range(20))
        line = f"// XD attr name chaine(into=[{values}]) syn OPT short desc"
        self.assertGreater(len(line), 120)
        self.assertTrue(is_unsplittable(line))

    def test_xd_block_header_with_very_long_name_is_unsplittable(self):
        from trustify.modernize import is_unsplittable

        long_name = "a" * 100
        line = f"// XD {long_name} objet_u {long_name} BRACE desc"
        self.assertGreater(len(line), 120)
        self.assertTrue(is_unsplittable(line))

    def test_long_description_with_short_prefix_is_splittable(self):
        from trustify.modernize import is_unsplittable

        long_desc = " ".join(["word"] * 40)
        line = f"// XD a objet_u a BRACE {long_desc}"
        self.assertGreater(len(line), 120)
        self.assertFalse(is_unsplittable(line))

    def test_xd_add_p_with_long_type_is_unsplittable(self):
        from trustify.modernize import is_unsplittable

        values = ",".join(f'"v_{i}"' for i in range(20))
        line = f"// XD_ADD_P chaine(into=[{values}]) short desc"
        self.assertGreater(len(line), 120)
        self.assertTrue(is_unsplittable(line))

    def test_xd_add_p_with_short_type_is_splittable(self):
        from trustify.modernize import is_unsplittable

        long_desc = " ".join(["word"] * 40)
        line = f"// XD_ADD_P int {long_desc}"
        self.assertGreater(len(line), 120)
        self.assertFalse(is_unsplittable(line))

    def test_non_xd_line_returns_false(self):
        from trustify.modernize import is_unsplittable

        line = "int x = 42; // ordinary comment, even if very long " + ("y" * 100)
        self.assertGreater(len(line), 120)
        self.assertFalse(is_unsplittable(line))

    def test_short_xd_line_returns_false(self):
        from trustify.modernize import is_unsplittable

        line = "// XD a objet_u a BRACE short desc"
        self.assertLess(len(line), 120)
        self.assertFalse(is_unsplittable(line))


class _NoTrustRoot:
    """Context manager that pops $TRUST_ROOT for the duration of a block.

    Kept as belt-and-suspenders: `modernize._resolve_src_dirs` no longer
    applies the `$TRUST_ROOT` env fallback (scope is decided by the CLI),
    but stripping the env keeps these unit tests robust to any future
    re-introduction and documents that the mock-baltik scope is the only
    intended input.
    """

    def __enter__(self):
        self._saved = os.environ.pop("TRUST_ROOT", None)
        return self

    def __exit__(self, *exc):
        if self._saved is not None:
            os.environ["TRUST_ROOT"] = self._saved


class TestWorkspaceWalker(unittest.TestCase):
    def test_trust_root_env_does_not_widen_baltik_scope(self):
        """Regression: with `trust_root=None` and `$TRUST_ROOT` set in the
        shell, modernize must rewrite ONLY the passed baltik — not the
        whole TRUST tree. The scope is owned by the CLI's
        `_scope_modernize`; the impl no longer re-applies the env
        fallback (which used to silently fold in all of TRUST whenever
        `$TRUST_ROOT` was set)."""
        from tests._baltik_fixtures import mock_workspace
        from trustify.modernize import modernize

        legacy = "// XD legacy_block objet_u legacy_block 0 a block with the legacy brace flag\n"
        with mock_workspace() as ws:
            trust = ws.root / "trust"
            (trust / "src").mkdir(parents=True)
            (trust / "src" / "markers.xd").write_text(legacy)
            baltik = ws.baltik("only_me", xd_blocks=[legacy])
            saved = os.environ.get("TRUST_ROOT")
            os.environ["TRUST_ROOT"] = str(trust)
            try:
                report = modernize(projects=[str(baltik)], trust_root=None, apply=True)
            finally:
                if saved is None:
                    os.environ.pop("TRUST_ROOT", None)
                else:
                    os.environ["TRUST_ROOT"] = saved
            changed = {fr.path for fr in report.changed_files}
            self.assertIn(str(baltik / "src" / "markers_0.xd"), changed)
            self.assertNotIn(str(trust / "src" / "markers.xd"), changed)
            self.assertIn("NO_BRACE", (baltik / "src" / "markers_0.xd").read_text())
            self.assertIn(" 0 a block", (trust / "src" / "markers.xd").read_text())

    def test_modernize_dry_run_reports_diff_and_counts(self):
        from tests._baltik_fixtures import mock_workspace
        from trustify.modernize import modernize

        with _NoTrustRoot(), mock_workspace() as ws:
            ws.baltik(
                "modernize_fixture",
                xd_blocks=[
                    "// XD legacy_block objet_u legacy_block 0 a block with the legacy brace flag\n",
                    "// XD attr legacy_opt int legacy_opt 1 an attribute with the legacy opt flag\n",
                    "// XD long_block objet_u long_block BRACE " + "x" * 150 + "\n",
                ],
            )
            project = str(ws._baltiks["modernize_fixture"])
            report = modernize(projects=[project], trust_root=None, apply=False)
            self.assertGreater(report.summary["brace"], 0)
            self.assertGreater(report.summary["opt"], 0)
            self.assertEqual(report.summary["unsplittable"], 1)
            for fr in report.files:
                if fr.diff:
                    on_disk = Path(fr.path).read_text()
                    self.assertNotEqual(on_disk, fr.new_content, "dry-run wrote to disk")

    def test_modernize_apply_preserves_source_when_write_fails(self):
        """Audit 2.3: modernize used to call `Path.write_text` directly,
        which truncates then writes. A SIGINT / disk-full between
        truncate and the final write would leave the source file
        truncated or partially written, with no backup. The fix is
        write-to-temp + os.replace. Simulate a rename failure and
        assert the source content survives unchanged."""
        from unittest.mock import patch

        from tests._baltik_fixtures import mock_workspace
        from trustify.modernize import modernize

        with _NoTrustRoot(), mock_workspace() as ws:
            ws.baltik(
                "atomic_fixture",
                xd_blocks=[
                    "// XD legacy_block objet_u legacy_block 0 a block with the legacy brace flag\n",
                ],
            )
            project = str(ws._baltiks["atomic_fixture"])
            src_files = list(Path(project, "src").rglob("*.xd"))
            self.assertEqual(len(src_files), 1)
            target = src_files[0]
            original = target.read_text()

            # The fix routes apply-writes through os.replace; patching
            # it forces the rename to fail AFTER the temp file is
            # written. The bug's direct `Path.write_text` would never
            # call os.replace, so the patch would be a no-op and the
            # file would still get overwritten — the assertion below
            # catches that case.
            with (
                patch("trustify.modernize.os.replace", side_effect=OSError("simulated rename failure")),
                self.assertRaises(OSError),
            ):
                modernize(projects=[project], trust_root=None, apply=True)

            self.assertEqual(
                target.read_text(),
                original,
                "modernize must not corrupt the source when the atomic write fails",
            )

    def test_modernize_apply_skips_symlinked_files(self):
        """Audit 2.4: a symlinked .xd file inside the baltik src/
        must NOT be rewritten through — that would mutate the target
        (which can live outside the workspace, e.g. shared TRUST
        sources). Skip symlinks during the walk."""
        from tests._baltik_fixtures import mock_workspace
        from trustify.modernize import modernize

        with _NoTrustRoot(), mock_workspace() as ws:
            baltik = ws.baltik("symlink_fixture")
            # File OUTSIDE the baltik src/ — must not be touched.
            outside = ws.root / "outside.xd"
            outside_content = "// XD legacy_block objet_u legacy_block 0 a block with the legacy brace flag\n"
            outside.write_text(outside_content)
            # Symlink from inside the baltik src/ to that outside file.
            link = baltik / "src" / "markers_sym.xd"
            link.symlink_to(outside)

            project = str(baltik)
            modernize(projects=[project], trust_root=None, apply=True)

            self.assertEqual(
                outside.read_text(),
                outside_content,
                "modernize followed a symlink and rewrote a file outside the workspace",
            )

    def test_modernize_apply_writes_and_is_idempotent(self):
        from tests._baltik_fixtures import mock_workspace
        from trustify.modernize import modernize

        with _NoTrustRoot(), mock_workspace() as ws:
            ws.baltik(
                "modernize_apply_fixture",
                xd_blocks=[
                    "// XD legacy_block objet_u legacy_block 0 a block with the legacy brace flag\n",
                ],
            )
            project = str(ws._baltiks["modernize_apply_fixture"])
            r1 = modernize(projects=[project], trust_root=None, apply=True)
            self.assertGreater(r1.summary["brace"], 0)
            r2 = modernize(projects=[project], trust_root=None, apply=True)
            self.assertEqual(r2.summary["brace"], 0)
            self.assertEqual(r2.summary["opt"], 0)
            self.assertEqual(r2.summary["splits"], 0)


class TestAtomicWriteTextSymlink(unittest.TestCase):
    """`_atomic_write_text` is shared by modernize and batch-format
    (via api.batch_format). When the target path is a symlink, the
    naive `os.replace(tmp, p)` overwrites the SYMLINK with a regular
    file — the link is gone and the original target file silently
    keeps its old content. Resolve the symlink so the write
    atomically updates the TARGET and the link stays intact (the
    semantic vim, emacs, sed --follow-symlinks all use)."""

    def test_atomic_write_through_symlink_updates_target_keeps_link(self):
        import os as _os
        import tempfile as _tempfile

        from trustify.modernize import _atomic_write_text

        with _tempfile.TemporaryDirectory() as outer:
            outside = Path(outer) / "outside.txt"
            outside.write_text("original\n")
            with _tempfile.TemporaryDirectory() as d:
                link = Path(d) / "linked.txt"
                _os.symlink(outside, link)

                _atomic_write_text(str(link), "new content\n")

                # The link path still resolves through a symlink.
                self.assertTrue(link.is_symlink(), "atomic write must NOT replace the symlink")
                # And both the link and the target read the new content.
                self.assertEqual(link.read_text(), "new content\n")
                self.assertEqual(outside.read_text(), "new content\n")

    def test_atomic_write_regular_file_unchanged_behaviour(self):
        # Sanity: a regular file (not a symlink) is still written
        # atomically into its own parent dir.
        import tempfile as _tempfile

        from trustify.modernize import _atomic_write_text

        with _tempfile.TemporaryDirectory() as d:
            target = Path(d) / "plain.txt"
            target.write_text("original\n")
            _atomic_write_text(str(target), "new content\n")
            self.assertEqual(target.read_text(), "new content\n")
            self.assertFalse(target.is_symlink())

    def test_atomic_write_broken_symlink_creates_target(self):
        # A symlink pointing at a non-existent target: write through it.
        # `os.path.realpath` returns the unresolved path; the write
        # creates a new file there. The symlink ends up pointing at a
        # real file.
        import os as _os
        import tempfile as _tempfile

        from trustify.modernize import _atomic_write_text

        with _tempfile.TemporaryDirectory() as outer:
            ghost = Path(outer) / "ghost.txt"  # never created
            with _tempfile.TemporaryDirectory() as d:
                link = Path(d) / "linked.txt"
                _os.symlink(ghost, link)
                self.assertFalse(ghost.exists())

                _atomic_write_text(str(link), "fresh\n")

                self.assertTrue(link.is_symlink())
                self.assertTrue(ghost.exists())
                self.assertEqual(ghost.read_text(), "fresh\n")


if __name__ == "__main__":
    unittest.main()
