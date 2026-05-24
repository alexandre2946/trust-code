"""Runtime hacks applied to the generated parser module.

Historically the body of this file was appended verbatim to the bottom of
``trustify_gen.py`` during generation. That baked a copy of every hack into
every cache entry and prevented ``pip install -U trustify`` from refreshing
older caches. Phase 1 of the reshaping switches to a callable hook: the
generated parser module's footer calls ``apply(sys.modules[__name__])`` with
itself, and we look every name up via ``module.__dict__``. The hacks see the
generated classes; the hack source no longer rides along in the cache.

Every name is fetched defensively with ``.get`` so the hook is a no-op when
applied against a reduced fixture (e.g. ``TRAD_2_adr_simple``) that doesn't
define the full TRUST class hierarchy.
"""

import contextlib
from typing import Annotated

from trustify.core.naming import ToParserName


def apply(module):
    """Apply Phase-1-extracted hacks to ``module`` (the freshly-imported
    generated parser module). The function mutates ``module.__dict__``.
    Safe to call multiple times.
    """
    g = module.__dict__
    if g.get("Objet_u") is None:
        # No generated content at all — nothing to hack.
        return

    _apply_skip_type_reads(g, module)
    _apply_listeqn_write_hack(g)
    _apply_dynamic_liste_equations(g, module)
    _apply_problem_generic_attribute_handler(g, module)


def _apply_skip_type_reads(g, module):
    """For Milieu / Mor_eqn / Corps_postraitement / Pb_base subclasses, the
    attribute name in the dataset IS the class name, so the parser must not
    re-read the type. Counter-example for Eqn_base_Parser: in Scalaires_passifs
    we DO want the type read (see WC_multi_2_3_espece.data). Same for
    Fluide_base_Parser: when a problem class declares a polymorphic-fluid
    attribute (e.g. `fluide0` of type `fluide_base` inside Fluide_Diphasique
    / Fluide_Diphasique_IJK in TrioCFD), the dataset shape is
    `fluide0 Fluide_Incompressible { ... }` — the attribute name is NOT the
    class name, so the type discriminator must be read.
    """
    Milieu_base_Parser = g.get("Milieu_base_Parser")
    Mor_eqn_Parser = g.get("Mor_eqn_Parser")
    Corps_postraitement_Parser = g.get("Corps_postraitement_Parser")
    Pb_base_Parser = g.get("Pb_base_Parser")
    Eqn_base_Parser = g.get("Eqn_base_Parser")
    Fluide_base_Parser = g.get("Fluide_base_Parser")
    bases = tuple(
        b for b in (Milieu_base_Parser, Mor_eqn_Parser, Corps_postraitement_Parser, Pb_base_Parser) if b is not None
    )
    if not bases:
        return
    for c in module._all_constrain_base_parser:
        # Per-iteration tolerance: a bad entry in _all_constrain_base_parser
        # shouldn't abort the whole sweep; issubclass() can raise TypeError
        # when handed a non-class value. Narrowed from `except Exception`
        # (audit 5.3) so an AttributeError from a typo'd `_read_type`
        # write — or any other real bug in the hack body — surfaces
        # instead of being silently swallowed.
        try:
            if issubclass(c, bases):
                c._read_type = False
        except TypeError:  # noqa: PERF203
            pass
    if Eqn_base_Parser is not None:
        Eqn_base_Parser._read_type = True
    if Fluide_base_Parser is not None:
        Fluide_base_Parser._read_type = True


def _apply_listeqn_write_hack(g):
    """When writing back a Listeqn (a list of equations of varying derived
    type), the item type MUST be output even though it wasn't read.
    Motivating dataset: Quasi_Comp_Cond_GP_VDF_FM.data.
    """
    ListOfBase_Parser = g.get("ListOfBase_Parser")
    Listeqn_Parser = g.get("Listeqn_Parser")
    if ListOfBase_Parser is None or Listeqn_Parser is None:
        return

    def toDSToken_hack(self):
        self._outputItemType = True
        return ListOfBase_Parser.toDatasetTokens(self)

    Listeqn_Parser.toDatasetTokens = toDSToken_hack


def _apply_dynamic_liste_equations(g, module):
    """Add a hidden ``liste_equations`` field to every Problem_read_generic
    subclass so the parser can stash arbitrary trailing equations (front-
    tracking / multi-phase generic problems).
    """
    Eqn_base = g.get("Eqn_base")
    Problem_read_generic = g.get("Problem_read_generic")
    if Eqn_base is None or Problem_read_generic is None:
        return
    # _all_constrain_base_pyd depends on Objet_u being present in the module;
    # bail out cleanly if the schema isn't fully loaded yet.
    if g.get("Objet_u") is None:
        return
    to_modif = [c for c in module._all_constrain_base_pyd if issubclass(c, Problem_read_generic)]
    for c in to_modif:
        c._synonyms["liste_equations"] = []
        c_new = c.with_fields(liste_equations=(Annotated[list[Eqn_base], "Listeqn"], []))
        c_new.__doc__ = c.__doc__
        g[c.__name__] = c_new


def _apply_problem_generic_attribute_handler(g, module):
    """Override ``handleUnexpectedAttribute`` on Problem_read_generic_Parser:
    if the unknown attribute name appears in ``solved_equations``, treat it as
    an equation alias and parse the following block as the declared equation
    class — appending the resulting parser object to the hidden
    ``liste_equations`` list created by ``_apply_dynamic_liste_equations``.

    Also installs a class-level ``toDatasetTokens`` wrapper that handles
    the FTD_IJK inner-param-block write-back (audit 6.10 — previous
    version monkey-patched ``self.toDatasetTokens`` from inside the
    handler, leaking the rebind onto every instance the handler ever
    ran against and inverting the normal class-level override path).
    The wrapper dispatches on per-instance state ``_ftd_ijk_param_block``
    and falls through to the original method when the block isn't set.
    """
    Problem_read_generic_Parser = g.get("Problem_read_generic_Parser")
    Builtin_Parser = g.get("Builtin_Parser")
    TRUSTTokens = g.get("TRUSTTokens")
    Eqn_base = g.get("Eqn_base")
    if not all((Problem_read_generic_Parser, Builtin_Parser, TRUSTTokens, Eqn_base)):
        return

    def handleUnexpectedAttribute_pb_generic(self, stream, tok, nams):
        # FTD_IJK-style inner brace block: Probleme_FTD_IJK_base::readOn
        # calls `param.lire_avec_accolades(is)` after reading the equations,
        # which reads an extra `{ key val ... }` block for problem-specific
        # parameters (nom_sauvegarde, sauvegarder_xyz, nom_reprise). The
        # block is often empty, sometimes carries a few of those keys.
        # When the outer parser hits `{` here, slurp the whole inner block
        # as raw tokens; the class-level toDatasetTokens wrapper installed
        # below picks the stashed state up at write-back time.
        if tok == "{":
            opening = stream.lastReadTokens()
            captured_low = list(opening.low())
            captured_orig = list(opening.orig())
            depth = 1
            while depth > 0:
                ttk = stream.nextLow()
                tok_obj = stream.lastReadTokens()
                captured_low.extend(tok_obj.low())
                captured_orig.extend(tok_obj.orig())
                if ttk == "{":
                    depth += 1
                elif ttk == "}":
                    depth -= 1
            block_tok = TRUSTTokens(low=captured_low, orig=captured_orig)
            # Stash for write-back. The block lands in the output stream
            # right before the next attribute that the outer parser reads
            # after we return — typically `postraitement`. We don't know
            # the next attribute name yet (we are mid-loop), so remember
            # the current _attrInOrder length and let the toDatasetTokens
            # wrapper prepend the block tokens to whichever attr ends up
            # at that position.
            self._ftd_ijk_param_block = block_tok
            self._ftd_ijk_inner_pos = len(self._attrInOrder)
            return

        slv_eq = {}
        slv_eq_att = []
        with contextlib.suppress(Exception):
            slv_eq_att = self._pyd_value.solved_equations
        for two_words in slv_eq_att:
            slv_eq[two_words.mot_2] = two_words.mot_1
        if tok not in slv_eq:
            raise self.GenErr(
                stream,
                f"Unexpected attribute or equation alias '{tok}' in keyword '{nams}'",
                token=tok,
                attr=tok,
                kind="unexpected-attribute",
            ) from None
        tok_full = stream.lastReadTokens()
        self._attr_ok["liste_equations"] = True
        eq_cls = getattr(module, ToParserName(slv_eq[tok]))
        self.Dbg(f"@FUNC@ About to ReadFromTokens class '{eq_cls}' alias '{tok}'")
        obj = eq_cls.ReadFromTokens(stream)
        obj._parser._tokens["cls_nam"] = tok_full
        if self._pyd_value.liste_equations == []:
            self._attrInOrder.append("liste_equations")
            self._leafParsers["liste_equations"] = Builtin_Parser.InstanciateFromBuiltin(
                Annotated[list[Eqn_base], "Listeqn"]
            )
            self._tokens["liste_equations"] = TRUSTTokens()

            def dummyGetBrace(_br):
                return ""

            self._leafParsers["liste_equations"].getBraceTokens = dummyGetBrace
        self._pyd_value.liste_equations.append(obj)

    Problem_read_generic_Parser.handleUnexpectedAttribute = handleUnexpectedAttribute_pb_generic

    # Class-level toDatasetTokens wrapper. Dispatches on per-instance
    # `_ftd_ijk_param_block`: absent (the common case) → delegate to
    # original; present → splice the slurped block before the attr
    # that landed at _ftd_ijk_inner_pos. Audit 6.10 — eliminates the
    # instance-method monkey-patch that the handler used to perform.
    original_toDS = Problem_read_generic_Parser.toDatasetTokens

    def toDatasetTokens_with_ftd_ijk_writeback(self):
        block = getattr(self, "_ftd_ijk_param_block", None)
        if block is None:
            return original_toDS(self)
        pos = self._ftd_ijk_inner_pos
        block_low, block_orig = list(block.low()), list(block.orig())
        saved = None
        if pos < len(self._attrInOrder):
            next_attr = self._attrInOrder[pos]
            old_tok = self._tokens.get(next_attr)
            if old_tok is not None:
                saved = (next_attr, old_tok)
                self._tokens[next_attr] = TRUSTTokens(
                    low=block_low + list(old_tok.low()),
                    orig=block_orig + list(old_tok.orig()),
                )
        try:
            s = original_toDS(self)
        finally:
            if saved is not None:
                self._tokens[saved[0]] = saved[1]
        if saved is None:
            # No attr after — splice the block right before the outer
            # closing brace.
            close_tail = self._tokens.get("}", TRUSTTokens()).orig()
            s = list(s)
            n = len(close_tail)
            s = s[:-n] + block_orig + s[-n:] if n and s[-n:] == close_tail else s + block_orig
        return s

    Problem_read_generic_Parser.toDatasetTokens = toDatasetTokens_with_ftd_ijk_writeback
