/****************************************************************************
* Copyright (c) 2026, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#include <Modele_turbulence_hyd_base.h>
#include <Op_Diff_PolyMAC_CDO_Face.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Champ_Face_PolyMAC_CDO.h>
#include <Domaine_PolyMAC_CDO.h>
#include <TRUSTTab_parts.h>
#include <Probleme_base.h>
#include <Synonyme_info.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTLists.h>
#include <Dirichlet.h>
#include <EChaine.h>

Implemente_instanciable_sans_constructeur( Op_Diff_PolyMAC_CDO_Face,"Op_Diff_PolyMAC_CDO_Face|Op_Dift_PolyMAC_CDO_Face_PolyMAC_CDO", Op_Diff_PolyMAC_CDO_base );
Add_synonym(Op_Diff_PolyMAC_CDO_Face, "Op_Diff_PolyMAC_CDO_var_Face");
Add_synonym(Op_Diff_PolyMAC_CDO_Face, "Op_Dift_PolyMAC_CDO_var_Face_PolyMAC_CDO");

Sortie& Op_Diff_PolyMAC_CDO_Face::printOn(Sortie& os) const { return Op_Diff_PolyMAC_CDO_base::printOn(os); }

Entree& Op_Diff_PolyMAC_CDO_Face::readOn(Entree& is) { return Op_Diff_PolyMAC_CDO_base::readOn(is); }

Op_Diff_PolyMAC_CDO_Face::Op_Diff_PolyMAC_CDO_Face()
{
  declare_support_masse_volumique(1);
}

void Op_Diff_PolyMAC_CDO_Face::completer()
{
  Op_Diff_PolyMAC_CDO_base::completer();
  if (polymac_flica5)
    {
      bool flag = Process::nproc() == 1 && le_dom_poly_->nb_faces() < 10000;
      // Pour des petites matrices, LU PETSc plus rapide que GMRES/ILU(0) ou MUMPS
      EChaine chl(flag ? "Petsc Cholesky_lapack { quiet }" : "Petsc gmres { precond block_jacobi_ilu { level 0 } quiet rtol 1.e-14 }");
      lire_solveur(chl);
      solveur.nommer("Op_Diff_PolyMAC_CDO_Face solver");
    }

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());

  if (domaine.domaine().nb_joints() && domaine.domaine().joint(0).epaisseur() < 1)
    Cerr << "Op_Diff_PolyMAC_CDO_Face : largeur de joint insuffisante (minimum 1)!" << finl, Process::exit();

  ch.init_ra();
  domaine.init_rf();
  domaine.init_m1();
  domaine.init_m2();

  if (polymac_flica5)
    porosite_e.ref(mon_equation->milieu().porosite_elem());

  if (que_suis_je() == "Op_Diff_PolyMAC_CDO_Face") return;

  const RefObjU& modele_turbulence = equation().get_modele(TURBULENCE);
  const Modele_turbulence_hyd_base& mod_turb = ref_cast(Modele_turbulence_hyd_base, modele_turbulence.valeur());
  const Champ_Fonc_base& alpha_t = mod_turb.viscosite_turbulente();
  associer_diffusivite_turbulente(alpha_t);
}

void Op_Diff_PolyMAC_CDO_Face::dimensionner(Matrice_Morse& mat) const
{
  if (polymac_flica5)
    {
      dimensionner_bloc(mat, -1);
      return;
    }
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const IntTab& e_f = domaine.elem_faces();
  const int nf_tot = domaine.nb_faces_tot(), na_tot = dimension < 3 ? domaine.domaine().nb_som_tot() : domaine.domaine().nb_aretes_tot();

  domaine.init_m2();

  Stencil stencil(0, 2);

  //partie vitesses : m2 Rf
  for (int e = 0; e < domaine.nb_elem_tot(); e++)
    {
      int idx = 0;
      for (int i = domaine.m2d(e); i < domaine.m2d(e + 1); i++, idx++)
        {
          const int f = e_f(e, idx);
          for (int j = domaine.m2i(i); f < domaine.nb_faces() && ch.fcl()(f, 0) < 2 && j < domaine.m2i(i + 1); j++)
            {
              const int fb = e_f(e, domaine.m2j(j));
              for (int k = domaine.rfdeb(fb); k < domaine.rfdeb(fb + 1); k++)
                stencil.append_line(f, nf_tot + domaine.rfji(k));
            }
        }
    }

  //partie vorticites : Ra m2 - m1 / nu
  for (int a = 0; a < (dimension < 3 ? domaine.nb_som() : domaine.domaine().nb_aretes()); a++)
    {
      for (int i = ch.radeb(a, 0); i < ch.radeb(a + 1, 0); i++)
        stencil.append_line(nf_tot + a, ch.raji(i));

      for (int i = domaine.m1deb(a); i < domaine.m1deb(a + 1); i++)
        stencil.append_line(nf_tot + a, nf_tot + domaine.m1ji(i, 0));
    }

  tableau_trier_retirer_doublons(stencil);
  Matrix_tools::allocate_morse_matrix(nf_tot + na_tot, nf_tot + na_tot, stencil, mat);
}

void Op_Diff_PolyMAC_CDO_Face::dimensionner_bloc(Matrice_Morse& mat, const int p) const
{
  /*
   0 | 1
   ---+---
   2 | 3
   */

  if (p > 3 || p < -1)
    Process::exit("Op_Diff_PolyMAC_CDO_Face::dimensionner_bloc : invalid bloc number! p must be in [-1, 3]");

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const IntTab& e_f = domaine.elem_faces();
  int i, j, k, a, e, f, fb, nf_tot = domaine.nb_faces_tot(), na_tot = dimension < 3 ? domaine.domaine().nb_som_tot() : domaine.domaine().nb_aretes_tot(), idx;

  domaine.init_m2();

  IntTab stencil(0, 2);
  VECT(IntTab) sp(4);
  for (int q = 0; q < 4; q++)
    {
      sp[q].resize(0, 2);
    }

  //partie vitesses : m2 Rf
  for (e = 0; e < domaine.nb_elem_tot(); e++)
    for (i = domaine.m2d(e), idx = 0; i < domaine.m2d(e + 1); i++, idx++)
      for (f = e_f(e, idx), j = domaine.m2i(i); f < domaine.nb_faces() && ch.fcl()(f, 0) < 2 && j < domaine.m2i(i + 1); j++)
        for (fb = e_f(e, domaine.m2j(j)), k = domaine.rfdeb(fb); k < domaine.rfdeb(fb + 1); k++)
          {
            stencil.append_line(f, nf_tot + domaine.rfji(k));
            sp[1].append_line(f, nf_tot + domaine.rfji(k));
          }

  //partie vorticites : Ra m2 - m1 / nu
  for (a = 0; a < (dimension < 3 ? domaine.nb_som() : domaine.domaine().nb_aretes()); a++)
    {
      for (i = ch.radeb(a, 0); i < ch.radeb(a + 1, 0); i++)
        {
          stencil.append_line(nf_tot + a, ch.raji(i));
          sp[2].append_line(a, ch.raji(i));
        }
      for (i = domaine.m1deb(a); i < domaine.m1deb(a + 1); i++)
        {
          stencil.append_line(nf_tot + a, nf_tot + domaine.m1ji(i, 0));
          sp[3].append_line(a, domaine.m1ji(i, 0));
        }
    }

  if (p == -1)
    {
      tableau_trier_retirer_doublons(stencil);
      Matrix_tools::allocate_morse_matrix(nf_tot + na_tot, nf_tot + na_tot, stencil, mat);
    }
  else
    {
      tableau_trier_retirer_doublons(sp[p]);
      const int nx = (p <= 1) ? nf_tot : na_tot;
      const int ny = (p == 0 || p == 2) ? nf_tot : na_tot;
      Matrix_tools::allocate_morse_matrix(nx, ny, sp[p], mat);
    }
}

void Op_Diff_PolyMAC_CDO_Face::update_auxiliary_variables()
{
  update_auxiliary_variables(le_champ_inco->valeurs());
}

static Matrice_Morse AF, AA;
void Op_Diff_PolyMAC_CDO_Face::update_auxiliary_variables(DoubleTab& inco)
{
  if (!polymac_flica5) return;

  if (AF.nb_lignes() == 0)
    {
      dimensionner_bloc(AF, 2);
      dimensionner_bloc(AA, 3);
    }
  else
    {
      AF.clean();
      AA.clean();
    }

  DoubleTab sm(inco);
  sm = 0.;
  DoubleTab_parts inco_parts(inco);
  DoubleTab_parts sm_parts(sm);

  contribuer_a_avec(inco, AF);
  contribuer_a_avec(inco, AA);

  AF *= -1.;
  AF.ajouter_multvect(inco_parts[0], sm_parts[1]);
  if (AA.is_diagonal())
    {
      const int n = AA.nb_lignes();
      for (int i = 0; i < n; i++)
        for (auto k = AA.get_tab1()(i) - 1; k < AA.get_tab1()(i + 1) - 1; k++)
          inco_parts[1][i] = sm_parts[1][i] / AA.get_coeff()(k);
      inco_parts[1].echange_espace_virtuel();
    }
  else
    {
      set_solveur()->reinit();
      inco_parts[1] = 0.;
      set_solveur().resoudre_systeme(AA, sm_parts[1], inco_parts[1]);
    }
}

inline DoubleTab& Op_Diff_PolyMAC_CDO_Face::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  if (polymac_flica5)
    {
      const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
      const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces();
      const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
      const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
      const DoubleVect& pe = porosite_e, &ve = domaine.volumes();
      int i, j, k, e, f, fb, a, nf_tot = domaine.nb_faces_tot(), idx;

      update_nu();
      //partie vitesses : m2 Rf
      for (e = 0; e < domaine.nb_elem_tot(); e++)
        for (i = domaine.m2d(e), idx = 0; i < domaine.m2d(e + 1); i++, idx++)
          for (f = e_f(e, idx), j = domaine.m2i(i); f < domaine.nb_faces() && j < domaine.m2i(i + 1); j++)
            for (fb = e_f(e, domaine.m2j(j)), k = domaine.rfdeb(fb); k < domaine.rfdeb(fb + 1); k++)
              resu(f) -= domaine.m2c(j) * ve(e) * (e == f_e(f, 0) ? 1 : -1) * (e == f_e(fb, 0) ? 1 : -1) * pe(e) * domaine.rfci(k) * inco(nf_tot + domaine.rfji(k));

      //partie vorticites : Ra m2 - m1 / nu
      if (resu.dimension_tot(0) == nf_tot)
        return resu; //resu ne contient que la partie "faces"
      for (a = 0; a < (dimension < 3 ? domaine.nb_som() : domaine.domaine().nb_aretes()); a++)
        {
          //rotationnel : vitesses internes
          for (i = ch.radeb(a, 0); i < ch.radeb(a + 1, 0); i++)
            resu(nf_tot + a) -= ch.raci(i) * inco(ch.raji(i));
          //rotationnel : vitesses aux bords
          for (i = ch.radeb(a, 1); i < ch.radeb(a + 1, 1); i++)
            for (k = 0; k < dimension; k++)
              resu(nf_tot + a) -= ch.racf(i, k) * ref_cast(Dirichlet, cls[ch.fcl()(ch.rajf(i), 1)].valeur()).val_imp(ch.fcl()(ch.rajf(i), 2), k);
          // -m1 / nu
          for (i = domaine.m1deb(a); i < domaine.m1deb(a + 1); i++)
            resu(nf_tot + a) += domaine.m1ci(i) / (pe(domaine.m1ji(i, 1)) * nu_(domaine.m1ji(i, 1), 0)) * inco(nf_tot + domaine.m1ji(i, 0));
        }
      return resu;
    }

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();
  const int nf_tot = domaine.nb_faces_tot();

  update_nu();

  //partie vitesses : m2 Rf
  for (int e = 0; e < domaine.nb_elem_tot(); e++)
    {
      int idx = 0;
      for (int i = domaine.m2d(e); i < domaine.m2d(e + 1); i++, idx++)
        {
          const int f = e_f(e, idx);

          for (int j = domaine.m2i(i); j < domaine.m2i(i + 1); j++)
            if (f < domaine.nb_faces() && ch.fcl()(f, 0) < 2)
              {
                const int fb = e_f(e, domaine.m2j(j));

                for (int k = domaine.rfdeb(fb); k < domaine.rfdeb(fb + 1); k++)
                  resu(f) -= domaine.m2c(j) * ve(e) * (e == f_e(f, 0) ? 1 : -1) * (e == f_e(fb, 0) ? 1 : -1) * pe(e) * domaine.rfci(k) * inco(nf_tot + domaine.rfji(k));
              }
        }
    }

  //partie vorticites : Ra m2 - m1 / nu
  if (resu.dimension_tot(0) == nf_tot)
    return resu; //resu ne contient que la partie "faces"

  /* boucle aretes*/
  for (int a = 0; a < (dimension < 3 ? domaine.nb_som() : domaine.domaine().nb_aretes()); a++)
    {
      //rotationnel : vitesses internes
      for (int i = ch.radeb(a, 0); i < ch.radeb(a + 1, 0); i++)
        resu(nf_tot + a) -= ch.raci(i) * inco(ch.raji(i));

      //rotationnel : vitesses aux bords
      for (int i = ch.radeb(a, 1); i < ch.radeb(a + 1, 1); i++)
        for (int k = 0; k < dimension; k++)
          resu(nf_tot + a) -= ch.racf(i, k) * ref_cast(Dirichlet, cls[ch.fcl()(ch.rajf(i), 1)].valeur()).val_imp(ch.fcl()(ch.rajf(i), 2), k);

      // -m1 / nu
      for (int i = domaine.m1deb(a); i < domaine.m1deb(a + 1); i++)
        resu(nf_tot + a) += domaine.m1ci(i) / (pe(domaine.m1ji(i, 1)) * nu_(domaine.m1ji(i, 1), 0)) * inco(nf_tot + domaine.m1ji(i, 0));
    }
  return resu;
}

void Op_Diff_PolyMAC_CDO_Face::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  if (polymac_flica5)
    {
      const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
      const int nf_tot = domaine.nb_faces_tot(), na_tot = dimension < 3 ? domaine.domaine().nb_som_tot() : domaine.domaine().nb_aretes_tot();

      int i = -1;
      if (matrice.nb_lignes() == nf_tot && matrice.nb_colonnes() == nf_tot)
        i = 0;
      if (matrice.nb_lignes() == nf_tot && matrice.nb_colonnes() == na_tot)
        i = 1;
      if (matrice.nb_lignes() == na_tot && matrice.nb_colonnes() == nf_tot)
        i = 2;
      if (matrice.nb_lignes() == na_tot && matrice.nb_colonnes() == na_tot)
        i = 3;
      contribuer_bloc(inco, matrice, i);
      return;
    }

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();
  const int nf_tot = domaine.nb_faces_tot();

  update_nu();

  //partie vitesses : m2 Rf
  for (int e = 0; e < domaine.nb_elem_tot(); e++)
    {
      int idx = 0;

      for (int i = domaine.m2d(e); i < domaine.m2d(e + 1); i++, idx++)
        {
          const int f = e_f(e, idx);

          for (int j = domaine.m2i(i); j < domaine.m2i(i + 1); j++)
            if (f < domaine.nb_faces() && ch.fcl()(f, 0) < 2)
              {
                const int fb = e_f(e, domaine.m2j(j));

                for (int k = domaine.rfdeb(fb); k < domaine.rfdeb(fb + 1); k++)
                  matrice(f, nf_tot + domaine.rfji(k)) += domaine.m2c(j) * ve(e) * (e == f_e(f, 0) ? 1 : -1) * (e == f_e(fb, 0) ? 1 : -1) * pe(e) * domaine.rfci(k);
              }
        }
    }

  //partie vorticites : Ra m2 - m1 / nu
  for (int a = 0; a < (dimension < 3 ? domaine.nb_som() : domaine.domaine().nb_aretes()); a++)
    {
      //rotationnel : vitesses internes
      for (int i = ch.radeb(a, 0); i < ch.radeb(a + 1, 0); i++)
        {
          const int f = ch.raji(i);

          if (ch.fcl()(f, 0) < 2)
            matrice(nf_tot + a, f) += ch.raci(i);
        }

      // -m1 / nu
      for (int i = domaine.m1deb(a); i < domaine.m1deb(a + 1); i++)
        matrice(nf_tot + a, nf_tot + domaine.m1ji(i, 0)) -= domaine.m1ci(i) / (pe(domaine.m1ji(i, 1)) * nu_(domaine.m1ji(i, 1), 0));
    }
}

void Op_Diff_PolyMAC_CDO_Face::contribuer_bloc(const DoubleTab& inco, Matrice_Morse& matrice, const int ip) const
{
  if (ip > 3 || ip < -1)
    Process::exit("Op_Diff_PolyMAC_CDO_Elem::contribuer_bloc : invalid bloc number! p must be in [-1, 3]");

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const DoubleVect& pe = porosite_e, &ve = domaine.volumes();
  int i, j, k, e, f, fb, a, nf_tot = domaine.nb_faces_tot(), idx;

  update_nu();
  //partie vitesses : m2 Rf
  for (e = 0; e < domaine.nb_elem_tot(); e++)
    for (i = domaine.m2d(e), idx = 0; i < domaine.m2d(e + 1); i++, idx++)
      for (f = e_f(e, idx), j = domaine.m2i(i); f < domaine.nb_faces() && ch.fcl()(f, 0) < 2 && j < domaine.m2i(i + 1); j++)
        for (fb = e_f(e, domaine.m2j(j)), k = domaine.rfdeb(fb); k < domaine.rfdeb(fb + 1); k++)
          {
            if (ip == -1)
              matrice(f, nf_tot + domaine.rfji(k)) += domaine.m2c(j) * ve(e) * (e == f_e(f, 0) ? 1 : -1) * (e == f_e(fb, 0) ? 1 : -1) * pe(e) * domaine.rfci(k);
            else if (ip == 1)
              matrice(f, domaine.rfji(k)) += domaine.m2c(j) * ve(e) * (e == f_e(f, 0) ? 1 : -1) * (e == f_e(fb, 0) ? 1 : -1) * pe(e) * domaine.rfci(k);
          }

  //partie vorticites : Ra m2 - m1 / nu
  for (a = 0; a < (dimension < 3 ? domaine.nb_som() : domaine.domaine().nb_aretes()); a++)
    {
      //rotationnel : vitesses internes
      for (i = ch.radeb(a, 0); i < ch.radeb(a + 1, 0); i++)
        if ((ch.fcl()(f = ch.raji(i), 0) < 2) || ip > -1)
          {
            if (ip == -1)
              matrice(nf_tot + a, f) += ch.raci(i);
            else if (ip == 2)
              matrice(a, f) += ch.raci(i);
          }
      // -m1 / nu
      for (i = domaine.m1deb(a); i < domaine.m1deb(a + 1); i++)
        {
          if (ip == -1)
            matrice(nf_tot + a, nf_tot + domaine.m1ji(i, 0)) -= domaine.m1ci(i) / (pe(domaine.m1ji(i, 1)) * nu_(domaine.m1ji(i, 1), 0));
          else if (ip == 3)
            matrice(a, domaine.m1ji(i, 0)) -= domaine.m1ci(i) / (pe(domaine.m1ji(i, 1)) * nu_(domaine.m1ji(i, 1), 0));
        }
    }
}
