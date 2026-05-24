"""End-to-end smoke test for trustify.doc.generator.

Runs `api.generate_markdown` against the real $TRUST_ROOT schema and asserts the
filesystem invariants (files exist, family-hub count ≥ family count).
No markdown content assertions — per the design spec, content shape is
explicitly out of scope for testing.

Requires $TRUST_ROOT to be set. Skip otherwise.
"""

import os
import tempfile
import unittest
from pathlib import Path


@unittest.skipUnless("TRUST_ROOT" in os.environ, "TRUST_ROOT must be set")
class TestDocSmokeFullTrust(unittest.TestCase):
    def test_generation_produces_index_and_family_hubs(self):
        from trustify.api import generate_markdown, generate_schema
        from trustify.core.misc_utilities import import_parser_module
        from trustify.doc.generator import DocGenerator

        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            generate_markdown(out=out, trust_root=os.environ["TRUST_ROOT"])

            # Index file must exist.
            self.assertTrue((out / "keyword_reference.md").is_file(), "keyword_reference.md was not written")

            # At least one Markdown family file must exist.
            all_md = list(out.glob("kw_*.md"))
            self.assertGreater(len(all_md), 0, "no kw_*.md files were emitted — generator dropped every class")

            # Family-presence invariant: for every concrete keyword family
            # (after applying the skip-list), at least one `kw_<family>*.md`
            # file must exist. Catches silent-drop regressions without
            # needing to discriminate hub vs branch by filename shape —
            # family names contain underscores in practice (e.g.
            # pb_hydraulique), so any "count underscores in stem"
            # heuristic is brittle.
            schema_dir = generate_schema(trust_root=os.environ["TRUST_ROOT"])
            schema_mod = import_parser_module(schema_dir / "trustify_gen.py")
            families = {DocGenerator._ultimate_parent_name(c) for c in schema_mod.all_constrain_base_pyd()}
            families -= {"dataset", "declaration", "bloc_comment", "objet_u", "listobj"}
            missing = [f for f in families if not list(out.glob(f"kw_{f}*.md"))]
            self.assertEqual(missing, [], f"families with no kw_*.md output: {missing}")

            # End-to-end extras chain: figures/ should exist (after the
            # migration) and contain at least one of the migrated jpegs.
            figures = out / "figures"
            self.assertTrue(figures.is_dir(), "figures/ should have been created by api.generate_markdown image copy")
            jpegs = list(figures.glob("*.jpeg"))
            self.assertGreater(
                len(jpegs), 0, "at least one .jpeg from docs/trustify/extras/ should have been copied to figures/"
            )

    def test_index_renders_full_inheritance_tree(self):
        """Top-level index must show the FULL inheritance tree, not just
        family hubs. Family roots stay @subpage (page hierarchy intact);
        every descendant becomes @ref so users can scan the complete
        tree from one place. Each family root is preceded by a level-2
        heading so doxygen's right-side TOC can navigate to it.
        """
        from trustify.api import generate_markdown

        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            generate_markdown(out=out, trust_root=os.environ["TRUST_ROOT"])
            index = (out / "keyword_reference.md").read_text()

            # Every family-hub @subpage must be followed by at least one
            # nested @ref line at a deeper indent — otherwise the tree
            # didn't expand and we've regressed to the flat list.
            self.assertIn("@subpage ", index)
            self.assertIn("  - @ref ", index, "index has no nested @ref children — tree expansion is broken")

            # Sanity check: the index must NOT emit @subpage for anything
            # below a family root. The legacy "flat list" layout was all
            # @subpage; the new layout reserves @subpage for top-level
            # entries only. Detect a regression by counting @subpage
            # entries at >0 indent.
            deep_subpages = [ln for ln in index.splitlines() if ln.startswith("  ") and "@subpage " in ln]
            self.assertEqual(
                deep_subpages,
                [],
                f"@subpage leaked below a family root (only top-level entries should be @subpage): {deep_subpages[:3]}",
            )

            # Every @subpage line must be immediately preceded (after a
            # blank line) by a level-2 heading — that's what makes the
            # family root navigable from doxygen's right-side TOC.
            lines = index.splitlines()
            heading_count = sum(1 for ln in lines if ln.startswith("## "))
            subpage_count = sum(1 for ln in lines if "@subpage " in ln and not ln.startswith("  "))
            self.assertEqual(
                heading_count,
                subpage_count,
                f"heading count ({heading_count}) does not match top-level @subpage count ({subpage_count}) "
                "— every family root must carry a ## heading for TOC navigation.",
            )

    def test_kw_ref_path_override_writes_index_outside_out_dir(self):
        from trustify.api import generate_markdown

        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "out"
            ref = Path(d) / "elsewhere" / "keyword_reference.md"
            ref.parent.mkdir(parents=True)
            generate_markdown(out=out, trust_root=os.environ["TRUST_ROOT"], kw_ref_path=ref)
            self.assertTrue(ref.is_file(), "kw_ref_path override was not honoured")
            self.assertFalse((out / "keyword_reference.md").exists(), "index was written into --out despite override")

    def test_multi_origin_does_not_misclassify_trust_keywords_as_baltik(self):
        """Regression: with TRUST + a baltik project, every TRUST keyword
        used to land under the baltik's per-origin index because
        ``_origin_of`` read ``_infoMain`` from PYD classes that don't
        carry it, falling through to ``classifier.labels[0]`` (= the
        baltik, since baltiks come first in the mapping)."""
        from tests._baltik_fixtures import mock_workspace
        from trustify.api import generate_markdown

        with mock_workspace() as ws:
            # One XD declaration so the baltik contributes a keyword the
            # smoke test can find in its per-origin index.
            baltik = ws.baltik(
                "fake_baltik",
                xd_blocks=["// XD fake_baltik_marker_kw objet_lecture fake_baltik_marker_kw NO_BRACE Marker keyword."],
            )

            out = ws.root / "out"
            generate_markdown(
                out=out,
                trust_root=os.environ["TRUST_ROOT"],
                projects=[str(baltik)],
            )

            trust_idx = out / "keyword_reference_trust.md"
            baltik_idx = out / "keyword_reference_fake_baltik.md"
            self.assertTrue(trust_idx.is_file(), "per-origin TRUST index missing in multi-origin layout")
            self.assertTrue(baltik_idx.is_file(), "per-origin baltik index missing in multi-origin layout")

            trust_text = trust_idx.read_text()
            baltik_text = baltik_idx.read_text()

            # The bug symptom: TRUST index reads "TRUST does not introduce
            # any new keywords." while every TRUST family lands under the
            # baltik. Assert the opposite — TRUST owns many families,
            # baltik owns very few.
            trust_subpages = trust_text.count("@subpage ")
            baltik_subpages = baltik_text.count("@subpage ")
            self.assertGreater(
                trust_subpages,
                10,
                f"TRUST per-origin index has too few @subpage entries ({trust_subpages}) — "
                "regression: TRUST keywords are being misclassified as the baltik origin.",
            )
            self.assertLess(
                baltik_subpages,
                trust_subpages,
                f"baltik @subpage count ({baltik_subpages}) ≥ TRUST count ({trust_subpages}) — "
                "regression: TRUST keywords are being misclassified as the baltik origin.",
            )

            # The baltik's marker keyword must appear in a baltik-prefixed
            # page (its ultimate-parent family) and must NOT appear in any
            # TRUST-prefixed page. The marker inherits from `objet_lecture`,
            # so it ends up inlined under `kw_fake_baltik_objet_lecture.md`.
            baltik_pages = [p for p in out.glob("kw_fake_baltik_*.md") if "fake_baltik_marker_kw" in p.read_text()]
            trust_pages = [p for p in out.glob("kw_trust_*.md") if "fake_baltik_marker_kw" in p.read_text()]
            self.assertTrue(
                baltik_pages,
                "baltik-defined keyword is not rendered under any baltik-prefixed page",
            )
            self.assertFalse(
                trust_pages,
                f"baltik-defined keyword leaked into TRUST-prefixed pages: {[p.name for p in trust_pages]}",
            )


if __name__ == "__main__":
    unittest.main()
