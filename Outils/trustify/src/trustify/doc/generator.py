"""Doc generator orchestrator.

`DocGenerator(schema_dir, classifier, extras).generate(out_dir, kw_ref_path=None)`
imports the generated parser module, groups classes by (origin,
family), and emits one Markdown page per family/branch/leaf plus index
pages. The single-origin layout mirrors ga/doxy's legacy file naming so
hand-written ``@ref kw_<family>`` references in TRUST's docs keep
working. ``extras`` is a ``trustify.doc.extras.ExtrasIndex`` threaded
through to the renderer for ``\\input`` / ``\\includeimage`` resolution.
"""

from pathlib import Path

from trustify.doc.origin import OriginClassifier, origin_slug

_FAMILIES_TO_SKIP = {"dataset", "declaration", "bloc_comment"}
_TYPE_MAP = {"objet_u": "Objet_U", "listobj": "list of objects"}


class DocGenerator:
    """Drive the keyword reference manual generation.

    Construction validates the generated-module path; ``generate`` does
    the actual work. Stateless between calls except for the
    ``_TRUST_BASE_CLS`` populated when the module is imported.
    """

    def __init__(self, schema_dir: Path, classifier: OriginClassifier, extras):
        """Initialise the generator.

        ``schema_dir`` is a cache entry directory that contains
        ``trustify_gen.py`` (the parser module) and its sibling
        ``trustify_gen_pyd_<digest>.py``. ``extras`` is an
        ``ExtrasIndex`` instance used by the renderer to resolve
        ``\\input`` and ``\\includeimage`` directives.
        """
        schema_dir = Path(schema_dir)
        parser_path = schema_dir / "trustify_gen.py"
        if not parser_path.is_file():
            raise FileNotFoundError(f"DocGenerator: generated parser not found at {parser_path}")
        self.schema_dir = schema_dir
        self.classifier = classifier
        self.extras = extras
        self._trust_base_cls = None
        self._schema_mod = None

    @staticmethod
    def _ultimate_parent_name(cls) -> str:
        """Return the lowercase name of `cls`'s ultimate parent.

        Exposed for the smoke test's family-count invariant. Walks
        ``__bases__`` upwards, stopping at the first class with multiple
        bases or whose only base is the TRUST root.
        """
        # Short-circuit: if `cls` IS Objet_u (the TRUST root), don't
        # walk past it into pydantic.BaseModel / object — return
        # "objet_u" directly so the caller's _TYPE_MAP skip-list
        # works.
        if cls.__name__ == "Objet_u":
            return "objet_u"
        current = cls
        while True:
            bases = current.__bases__
            if len(bases) != 1:
                return current.__name__.lower()
            parent = bases[0]
            # Heuristic: TRUST's root is Objet_u; recognise it by name to
            # avoid importing the generated module here.
            if parent.__name__ == "Objet_u":
                return current.__name__.lower()
            current = parent

    def _get_parents(self, cls):
        """Return (direct_parent, ultimate_parent) for a keyword class.

        Direct parent  = one step up the inheritance chain.
        Ultimate parent = highest ancestor below ``Objet_u``.
        Used to decide where each keyword lands in the page hierarchy.
        """
        first_parent = None

        def walk(c):
            nonlocal first_parent
            bases = c.__bases__
            if len(bases) != 1:
                if first_parent is None:
                    first_parent = c
                return c
            if bases[0] is self._trust_base_cls:
                if first_parent is None:
                    first_parent = c
                return c
            if first_parent is None:
                first_parent = bases[0]
            return walk(bases[0])

        ultimate = walk(cls)
        return first_parent, ultimate

    def _origin_of(self, cls) -> str:
        # _infoMain (a `(project, path, line)` triple post-refactor) is
        # emitted only on parser classes, not on the PYD siblings that
        # all_constrain_base_pyd() returns. Look up the parser to read
        # the triple. Falls back to the first registered label when the
        # class is a synthetic root (objet_u / listobj_impl) without a
        # parser sibling.
        try:
            parser_cls = self._schema_mod.get_parser_from_pyd(cls)
        except AttributeError:
            parser_cls = None
        info = getattr(parser_cls, "_infoMain", None) if parser_cls is not None else None
        project = info[0] if info and len(info) >= 1 else ""
        return self.classifier.classify_project(project) or self.classifier.labels[0]

    INTRO_BLURB = (
        "This manual is the complete reference for all keywords accepted "
        "by TRUST `.data` input files. Keywords are extracted "
        "automatically from `// XD` annotations embedded in the C++ "
        "source code and assembled by the trustify pipeline — the content "
        "here always reflects the current state of the source.\n\n"
        "Keywords are grouped by their base class. Each section below "
        "covers one family: the base class description is shown first, "
        "followed by all concrete keywords that inherit from it, "
        "organised recursively by intermediate parent when the hierarchy "
        "has more than one level.\n\n"
    )

    def generate(self, out_dir: Path, kw_ref_path: Path | None = None) -> None:
        """Emit the full keyword reference manual into ``out_dir``.

        Steps:
          1. Import the generated parser module (loads all classes).
          2. Load the pydantic module for Objet_u.
          3. Group concrete classes by (origin, ultimate-parent family).
          4. Emit one hub page per family + branch sub-pages for
             non-leaf intermediates, with all leaves rendered inline.
          5. Emit indices: per-origin sub-indices (multi-origin only),
             top-level ``keyword_reference.md``.
        """
        from trustify.core.misc_utilities import import_parser_module
        from trustify.doc.render import render_keyword_card

        out_dir = Path(out_dir)
        # Load the parser module — exposes class objects and is the
        # namespace for sibling lookups during rendering.
        parser_mod_path = self.schema_dir / "trustify_gen.py"
        schema_mod = import_parser_module(parser_mod_path)
        self._trust_base_cls = schema_mod.module.Objet_u
        # _origin_of() reads _infoMain via the parser-class sibling — stash
        # the schema module so it can do the PYD-to-parser lookup without
        # threading schema_mod through every call site.
        self._schema_mod = schema_mod

        # Classify every keyword class by family & direct parent.
        all_cls = schema_mod.all_constrain_base_pyd()
        top_cls, direct_parent = {}, {}
        for c in all_cls:
            first, ultimate = self._get_parents(c)
            cn = c.__name__.lower()
            upn = ultimate.__name__.lower()
            if upn in _TYPE_MAP or cn in _TYPE_MAP:
                continue
            top_cls[c] = ultimate
            direct_parent[c] = first

        # Per-class origin + per-family origin.
        kw_origin = {c: self._origin_of(c) for c in top_cls}
        family_origin = {p: self._origin_of(p) for p in set(top_cls.values())}

        # Group keywords by (origin, family).
        by_origin = {label: {} for label in self.classifier.labels}
        for c, parent in top_cls.items():
            by_origin[kw_origin[c]].setdefault(parent, []).append(c)

        is_multi = self.classifier.is_multi
        per_origin_families = self._emit_families(out_dir, by_origin, family_origin, is_multi, render_keyword_card)
        self._emit_indices(out_dir, per_origin_families, is_multi, kw_ref_path)

    def _emit_families(self, out_dir, by_origin, family_origin, is_multi, render_card):
        """Write hub/branch/leaf pages for every (origin, family) bucket.

        Returns a dict mapping origin label → list of per-family records
        (``{family_root, page_id, family_name, children_map}``) — the
        ``children_map`` is what the index renderer walks to emit a
        full inheritance tree on the top-level / per-origin index.
        """
        per_origin_families = {label: [] for label in self.classifier.labels}

        for origin in self.classifier.labels:
            families = sorted(by_origin[origin].keys(), key=lambda c: c.__name__)
            for parent_cls in families:
                par_name = parent_cls.__name__.lower()
                if par_name in _FAMILIES_TO_SKIP:
                    continue

                kw_list = sorted(by_origin[origin][parent_cls], key=lambda c: c.__name__)
                lst_set = set(kw_list)

                page_id = self._family_hub_id(parent_cls, origin, is_multi)
                page_title = f"Keywords derived from {par_name}"

                native_family = family_origin.get(parent_cls) == origin

                # Build the intermediate-node sub-tree only for native
                # families. Non-native (baltik extending TRUST) keeps the
                # layout flat: every keyword is rendered directly under the
                # family root, no abstract intermediates re-emitted.
                children_map = {}
                parent_of = {}
                if native_family:
                    for c in kw_list:
                        current = c
                        while current is not parent_cls:
                            bases = current.__bases__
                            if len(bases) != 1:
                                break
                            up = bases[0]
                            if current not in parent_of:
                                parent_of[current] = up
                            children_map.setdefault(up, set()).add(current)
                            current = up
                    children_map = {node: sorted(kids, key=lambda c: c.__name__) for node, kids in children_map.items()}
                else:
                    children_map = {parent_cls: kw_list}

                per_origin_families[origin].append(
                    {
                        "family_root": parent_cls,
                        "page_id": page_id,
                        "family_name": par_name,
                        "children_map": children_map,
                    }
                )

                self._write_subtree(
                    parent_cls,
                    parent_cls,
                    page_id,
                    page_title,
                    out_dir,
                    origin,
                    is_multi,
                    native_family,
                    children_map,
                    parent_of,
                    lst_set,
                    render_card,
                    is_hub=True,
                )

        return per_origin_families

    def _write_subtree(
        self,
        family_root,
        cls,
        page_id,
        page_title,
        out_dir,
        origin,
        is_multi,
        native_family,
        children_map,
        parent_of,
        lst_set,
        render_card,
        is_hub,
    ):
        """Emit the page for `cls` and recurse into its branch children."""
        kids = children_map.get(cls, [])
        branches = [k for k in kids if k in children_map]
        leaves = [k for k in kids if k not in children_map]

        s = f"@page {page_id} {page_title}\n\n"

        if is_hub and not native_family:
            native_id = self._family_hub_id(cls, self._origin_of(cls), is_multi)
            s += (
                f"Keywords introduced by **{origin}** that extend the "
                f"`{cls.__name__}` family. See @ref {native_id} for the "
                f"base-class documentation and the canonical keyword list.\n\n"
                f"---\n\n"
            )
        elif not is_hub or cls in lst_set:
            direct = parent_of.get(cls, cls)
            s += render_card(cls, direct, self.extras)
            if kids:
                s += "\n\n---\n\n"

        for b in branches:
            s += f"- @subpage {self._branch_id(family_root, b, origin, is_multi)}\n"

        if branches and leaves:
            s += "\n"

        for i, leaf in enumerate(leaves):
            if i > 0:
                s += "\n\n---\n\n"
            s += render_card(leaf, parent_of.get(leaf, leaf), self.extras)

        fname = (
            self._family_hub_fname(family_root, origin, is_multi)
            if is_hub
            else self._branch_fname(family_root, cls, origin, is_multi)
        )
        (out_dir / fname).write_text(s)

        for b in branches:
            self._write_subtree(
                family_root,
                b,
                self._branch_id(family_root, b, origin, is_multi),
                f"Keywords derived from {b.__name__.lower()}",
                out_dir,
                origin,
                is_multi,
                native_family,
                children_map,
                parent_of,
                lst_set,
                render_card,
                is_hub=False,
            )

    def _emit_indices(self, out_dir, per_origin_families, is_multi, kw_ref_path):
        """Write the per-origin sub-indices (multi-origin only) and the
        top-level ``keyword_reference.md`` index.

        The index that surfaces a particular set of families
        (``keyword_reference.md`` in single-origin builds, each
        ``keyword_reference_<slug>.md`` in multi-origin builds) renders
        the **full** inheritance tree under every family root, not just
        the direct list of family hubs. Family roots are emitted as
        ``@subpage`` (the canonical doxygen page hierarchy stays
        intact); every descendant becomes a ``@ref`` so users can scan
        the complete tree from one place without losing their bearings
        in nested @subpage indirection.
        """
        out_dir = Path(out_dir)
        if is_multi:
            index_s = "@page KeywordsReference Keywords Reference Manual\n\n" + self.INTRO_BLURB
            for label in self.classifier.labels:
                index_s += f"- @subpage KW_{origin_slug(label)}\n"
            for label in self.classifier.labels:
                slug = origin_slug(label)
                sub_id = f"KW_{slug}"
                sub_s = f"@page {sub_id} {label} Keyword Reference\n\n"
                families = per_origin_families[label]
                if families:
                    sub_s += f"Keyword families documented under **{label}**.\n\n"
                    for fam in families:
                        sub_s += self._render_index_tree(fam, label, is_multi)
                else:
                    sub_s += f"*{label} does not introduce any new keywords.*\n"
                (out_dir / f"keyword_reference_{slug}.md").write_text(sub_s)
        else:
            only_label = self.classifier.labels[0]
            index_s = "@page KeywordsReference Keywords Reference Manual\n\n" + self.INTRO_BLURB
            for fam in per_origin_families[only_label]:
                index_s += self._render_index_tree(fam, only_label, is_multi)

        index_path = Path(kw_ref_path) if kw_ref_path is not None else out_dir / "keyword_reference.md"
        index_path.write_text(index_s)

    def _render_index_tree(self, fam, origin, is_multi) -> str:
        """Render one family's full inheritance tree as nested bullets,
        preceded by a level-2 heading so the family root shows up in the
        right-side TOC of the doxygen HTML build.

        The family root is emitted as ``@subpage`` to its hub page so
        doxygen records the parent-child page relationship. Every
        descendant is ``@ref``: branches link to their dedicated branch
        page; leaves link to their ``kw_<name>`` card anchor on
        whichever hub/branch page renders them inline.
        """
        family_root = fam["family_root"]
        children_map = fam["children_map"]
        page_id = fam["page_id"]
        family_name = family_root.__name__.lower()

        def _walk(cls, depth):
            indent = "  " * depth
            name = cls.__name__.lower()
            if cls is family_root:
                line = f'{indent}- @subpage {page_id} "{name}"\n'
            elif cls in children_map:
                branch_id = self._branch_id(family_root, cls, origin, is_multi)
                line = f'{indent}- @ref {branch_id} "{name}"\n'
            else:
                line = f'{indent}- @ref kw_{name} "{name}"\n'
            for kid in children_map.get(cls, []):
                line += _walk(kid, depth + 1)
            return line

        return f"## {family_name}\n\n{_walk(family_root, 0)}\n"

    @staticmethod
    def _family_hub_id(parent_cls, origin, is_multi) -> str:
        par = parent_cls.__name__.lower()
        return f"KW_{origin_slug(origin)}_{par}" if is_multi else f"KW_{par}"

    @staticmethod
    def _family_hub_fname(parent_cls, origin, is_multi) -> str:
        par = parent_cls.__name__.lower()
        return f"kw_{origin_slug(origin)}_{par}.md" if is_multi else f"kw_{par}.md"

    @staticmethod
    def _branch_id(family_cls, branch_cls, origin, is_multi) -> str:
        par = family_cls.__name__.lower()
        bn = branch_cls.__name__.lower()
        return f"KW_{origin_slug(origin)}_{par}_{bn}" if is_multi else f"KW_{par}_{bn}"

    @staticmethod
    def _branch_fname(family_cls, branch_cls, origin, is_multi) -> str:
        par = family_cls.__name__.lower()
        bn = branch_cls.__name__.lower()
        return f"kw_{origin_slug(origin)}_{par}_{bn}.md" if is_multi else f"kw_{par}_{bn}.md"
