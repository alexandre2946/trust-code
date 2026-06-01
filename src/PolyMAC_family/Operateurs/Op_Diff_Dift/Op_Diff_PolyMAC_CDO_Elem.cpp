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

#include <Modele_turbulence_scal_base.h>
#include <Echange_contact_PolyMAC_CDO.h>
#include <Op_Diff_PolyMAC_CDO_Elem.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Champ_Elem_PolyMAC_CDO.h>
#include <Probleme_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTTab_parts.h>
#include <EChaine.h>

Implemente_instanciable_sans_constructeur( Op_Diff_PolyMAC_CDO_Elem , "Op_Diff_PolyMAC_CDO_Elem|Op_Diff_PolyMAC_CDO_var_Elem" , Op_Diff_PolyMAC_CDO_base );
Implemente_instanciable( Op_Dift_PolyMAC_CDO_Elem , "Op_Dift_PolyMAC_CDO|Op_Dift_PolyMAC_CDO_var_P0_PolyMAC_CDO" , Op_Diff_PolyMAC_CDO_Elem );
Implemente_instanciable( Op_Diff_Nonlinear_PolyMAC_CDO_Elem, "Op_Diff_nonlinear_PolyMAC_CDO_Elem|Op_Diff_nonlinear_PolyMAC_CDO_var_Elem" , Op_Diff_PolyMAC_CDO_Elem );
Implemente_instanciable( Op_Dift_Nonlinear_PolyMAC_CDO_Elem, "Op_Dift_PolyMAC_CDO_nonlinear|Op_Dift_PolyMAC_CDO_var_P0_PolyMAC_CDO_nonlinear", Op_Diff_PolyMAC_CDO_Elem );

Sortie& Op_Diff_PolyMAC_CDO_Elem::printOn(Sortie& os) const { return Op_Diff_PolyMAC_CDO_base::printOn(os); }
Sortie& Op_Dift_PolyMAC_CDO_Elem::printOn(Sortie& os) const { return Op_Diff_PolyMAC_CDO_base::printOn(os); }
Sortie& Op_Diff_Nonlinear_PolyMAC_CDO_Elem::printOn(Sortie& os) const { return Op_Diff_PolyMAC_CDO_base::printOn(os); }
Sortie& Op_Dift_Nonlinear_PolyMAC_CDO_Elem::printOn(Sortie& os) const { return Op_Diff_PolyMAC_CDO_base::printOn(os); }

Entree& Op_Diff_PolyMAC_CDO_Elem::readOn(Entree& is) { return Op_Diff_PolyMAC_CDO_base::readOn(is); }
Entree& Op_Diff_Nonlinear_PolyMAC_CDO_Elem::readOn(Entree& is) { return Op_Diff_PolyMAC_CDO_base::readOn(is); }
Entree& Op_Dift_PolyMAC_CDO_Elem::readOn(Entree& is) { return Op_Diff_PolyMAC_CDO_base::readOn(is); }
Entree& Op_Dift_Nonlinear_PolyMAC_CDO_Elem::readOn(Entree& is) { return Op_Diff_PolyMAC_CDO_base::readOn(is); }

Op_Diff_PolyMAC_CDO_Elem::Op_Diff_PolyMAC_CDO_Elem()
{
  declare_support_masse_volumique(1);
}

void Op_Diff_PolyMAC_CDO_Elem::completer()
{
  Op_Diff_PolyMAC_CDO_base::completer();
  if (polymac_flica5)
    {
      // For small systems, a direct factorization can be faster than iterative solvers.
      bool flag = Process::nproc() == 1 && le_dom_poly_->nb_elem() < 10000;
      EChaine chl(flag ? "Petsc Cholesky_lapack { quiet }" : "Petsc Cholesky { quiet }");
      lire_solveur(chl);
      solveur.nommer("Op_Diff_PolyMAC_CDO_Elem solver");
    }
  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  if (domaine.domaine().nb_joints() && domaine.domaine().joint(0).epaisseur() < 1)
    Cerr << "Op_Diff_PolyMAC_CDO_Elem : largeur de joint insuffisante (minimum 1)!" << finl, Process::exit();
  ch.fcl();
  int nb_comp = (equation().que_suis_je() == "Transport_K_Epsilon") ? 2 : ch.valeurs().line_size();
  flux_bords_.resize(domaine.premiere_face_int(), nb_comp);

  stab_ = que_suis_je().find("nonlinear") >= 0;
  if (stab_)
    {
      delta_e.resize(0, nb_comp), delta_f.resize(0, nb_comp), delta_f_int.resize(0, nb_comp, 2);
      domaine.domaine().creer_tableau_elements(delta_e), domaine.creer_tableau_faces(delta_f), domaine.creer_tableau_faces(delta_f_int);
    }
  delta_int_a_jour_ = delta_a_jour_ = (stab_ ? 0 : 1);

  if (!que_suis_je().debute_par("Op_Dift"))
    return;
  const RefObjU& modele_turbulence = equation().get_modele(TURBULENCE);
  const Modele_turbulence_scal_base& mod_turb = ref_cast(Modele_turbulence_scal_base, modele_turbulence.valeur());
  const Champ_Fonc_base& lambda_t = mod_turb.conductivite_turbulente();
  associer_diffusivite_turbulente(lambda_t);
}

void Op_Diff_PolyMAC_CDO_Elem::update_delta_int() const
{
  if (delta_int_a_jour_)
    return; //already done
  const DoubleTab& inco = equation().inconnue().valeurs();
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const IntTab& e_f = domaine.elem_faces();
  const DoubleVect& fs = domaine.face_surfaces(), &ve = domaine.volumes();
  int i, j, k, l, e, f, fb, ne_tot = domaine.nb_elem_tot(), n, N = inco.line_size();
  double fac, fac_n;

  update_nu(); //prerequisite: nu
  delta_f_int = 0;
  DoubleTrav nu_ef(e_f.dimension(1), N), de_num(N), de_den(N);
  for (e = 0; e < domaine.nb_elem_tot(); e++)
    {
      de_num = 0, de_den = 0; //numerator / denominator of delta_e
      remplir_nu_ef(e, nu_ef);
      for (i = 0, j = domaine.m2d(e); j < domaine.m2d(e + 1); i++, j++)
        {
          //'face-face' part
          for (f = e_f(e, i), k = domaine.w2i(j); k < domaine.w2i(j + 1); k++)
            for (fb = e_f(e, l = domaine.w2j(k)), fac = fs(f) * fs(fb) / ve(e) * domaine.w2c(k), n = 0; n < N; n++)
              {
                fac_n = fac * nu_ef(l, n) * (inco(ne_tot + fb, n) - inco(e, n));
                de_num(n) += fac_n, delta_f_int(f, n, 0) += fac_n;
                delta_f_int(f, n, 1) += std::fabs(inco(ne_tot + fb, n) - inco(ne_tot + f, n));
              }
          //'face-element' part
          for (n = 0; n < N; n++)
            {
              fac = std::fabs(inco(e, n) - inco(ne_tot + f, n));
              de_den(n) += fac, delta_f_int(f, n, 1) += fac;
            }
        }
      //we can compute delta_e now
      for (n = 0; n < N; n++)
        delta_e(e, n) = de_den(n) > 1e-8 ? std::fabs(de_num(n)) / de_den(n) : 0;
    }
  delta_e.echange_espace_virtuel(), delta_f_int.echange_espace_virtuel();
  delta_int_a_jour_ = 1;
}

void Op_Diff_PolyMAC_CDO_Elem::update_delta() const
{
  if (delta_a_jour_)
    return; //already done
  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  int i, f, n, N = ch.valeurs().line_size();

  //prerequisite: delta_int internally + in Echange_contact boundary conditions
  update_delta_int();
  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      ref_cast_non_const(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_coeffs();
  //final computation of delta_f
  for (f = 0; f < domaine.nb_faces_tot(); f++)
    for (n = 0; n < N; n++)
      {
        double n_d[2] = { delta_f_int(f, n, 0), delta_f_int(f, n, 1) };
        //contribution from the other problem via Echange_contact boundary conditions
        for (i = 0; ch.fcl()(f, 0) == 3 && i < 2; i++)
          n_d[i] += ref_cast(Echange_contact_PolyMAC_CDO, cls[ch.fcl()(f, 1)].valeur()).delta_int(ch.fcl()(f, 2), n, i);
        delta_f(f, n) = n_d[1] > 1e-8 ? std::fabs(n_d[0]) / n_d[1] : 0;
      }
  delta_f.echange_espace_virtuel();
  delta_a_jour_ = 1;
}

void Op_Diff_PolyMAC_CDO_Elem::dimensionner(Matrice_Morse& mat) const
{
  dimensionner_bloc(mat, -1);
}

void Op_Diff_PolyMAC_CDO_Elem::dimensionner_bloc(Matrice_Morse& mat, const int p) const
{
  /*
   0 | 1
   ---+---
   2 | 3
   */

  if (p > 3 || p < -1)
    Process::exit("Op_Diff_PolyMAC_CDO_Elem::dimensionner_bloc : invalid bloc number! p must be in [-1, 3]");

  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const IntTab& e_f = domaine.elem_faces();
  int i, j, k, l, e, f, ne_tot = domaine.nb_elem_tot(), nf_tot = domaine.nb_faces_tot(), n, N = ch.valeurs().line_size();
  domaine.init_m2();

  Stencil stencil(0, 2);
  VECT(Stencil) sp(4);
  for (int q = 0; q < 4; q++) sp[q].resize(0, 2);

  for (e = 0; e < domaine.nb_elem_tot(); e++)
    {
      //dependence on Te: diagonal -> faces around each element
      if (e < domaine.nb_elem())
        for (n = 0; n < N; n++)
          {
            stencil.append_line(N * e + n, N * e + n);
            sp[0].append_line(N * e + n, N * e + n);
          }
      for (i = 0; i < e_f.dimension(1) && (f = e_f(e, i)) >= 0; i++)
        for (n = 0; f < domaine.nb_faces() && n < N; n++)
          {
            stencil.append_line(N * (ne_tot + f) + n, N * e + n);
            sp[2].append_line(N * f + n, N * e + n);
          }

      //dependence on Tf
      for (j = 0, k = domaine.m2d(e); k < domaine.m2d(e + 1); j++, k++)
        for (f = e_f(e, j), l = domaine.w2i(k); l < domaine.w2i(k + 1); l++)
          {
            //upper blocks: divergence
            for (n = 0; e < domaine.nb_elem() && n < N; n++)
              {
                stencil.append_line(N * e + n, N * (ne_tot + e_f(e, domaine.w2j(l))) + n);
                sp[1].append_line(N * e + n, N * e_f(e, domaine.w2j(l)) + n);
              }

            //lower blocks: continuity
            for (n = 0; f < domaine.nb_faces() && n < N; n++)
              {
                stencil.append_line(N * (ne_tot + f) + n, N * (ne_tot + e_f(e, domaine.w2j(l))) + n);
                sp[3].append_line(N * f + n, N * e_f(e, domaine.w2j(l)) + n);
              }
          }
    }
  if (p == -1)
    {
      tableau_trier_retirer_doublons(stencil);
      Matrix_tools::allocate_morse_matrix(N * (ne_tot + nf_tot), N * (ne_tot + nf_tot), stencil, mat);
    }
  else
    {
      tableau_trier_retirer_doublons(sp[p]);
      const int nx = N * ((p <= 1) ? ne_tot : nf_tot);
      const int ny = N * ((p == 0 || p == 2) ? ne_tot : nf_tot);
      Matrix_tools::allocate_morse_matrix(nx, ny, sp[p], mat);
    }
}

void Op_Diff_PolyMAC_CDO_Elem::dimensionner_termes_croises(Matrice_Morse& matrice, const Probleme_base& autre_pb, int nl, int nc) const
{
  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  int i, j, k, l, f, n, N = ch.valeurs().line_size(), ne_tot = domaine.nb_elem_tot();

  Stencil stencil(0, 2);

  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      {
        const Echange_contact_PolyMAC_CDO& cl = ref_cast(Echange_contact_PolyMAC_CDO, cls[i].valeur());
        if (cl.nom_autre_pb() != autre_pb.le_nom())
          continue; //not our problem

        /* stencil */
        const Front_VF& fvf = ref_cast(Front_VF, cl.frontiere_dis());
        for (j = 0; j < cl.item.dimension(0); j++)
          for (k = 0, f = fvf.num_face(j); k < cl.item.dimension(1) && (l = cl.item(j, k)) >= 0; k++)
            for (n = 0; n < N; n++)
              stencil.append_line(N * (ne_tot + f) + n, N * l + n);
      }

  tableau_trier_retirer_doublons(stencil);
  Matrix_tools::allocate_morse_matrix(nl, nc, stencil, matrice);
}

void Op_Diff_PolyMAC_CDO_Elem::ajouter_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, DoubleTab& resu) const
{
  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  int i, j, k, l, f, n, N = ch.valeurs().line_size(), ne_tot = domaine.nb_elem_tot();

  //prerequisites: nu, delta internally + coeffs/delta in Echange_contact boundary conditions
  update_nu(), update_delta();
  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      ref_cast_non_const(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_coeffs(), ref_cast(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_delta();

  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      {
        const Echange_contact_PolyMAC_CDO& cl = ref_cast(Echange_contact_PolyMAC_CDO, cls[i].valeur());
        if (cl.nom_autre_pb() != autre_pb.le_nom())
          continue; //not our problem
        const Front_VF& fvf = ref_cast(Front_VF, cl.frontiere_dis());
        for (j = 0; j < fvf.nb_faces(); j++)
          {
            f = fvf.num_face(j);
            for (n = 0; n < N; n++)
              resu(ne_tot + f, n) -= cl.coeff(j, 0, n) * inco(ne_tot + f, n); //face's own term
            for (k = 0; k < cl.item.dimension(1) && (l = cl.item(j, k)) >= 0; k++)
              for (n = 0; n < N; n++)
                {
                  //operator
                  resu(ne_tot + f, n) -= cl.coeff(j, k + 1, n) * autre_inco(l, n);
                  //nonlinear correction
                  if (stab_)
                    resu(ne_tot + f, n) -= std::max(delta_f(f, n), cl.delta(j, k, n)) * (inco(ne_tot + f, n) - autre_inco(l, n));
                }
          }
      }
}

void Op_Diff_PolyMAC_CDO_Elem::contribuer_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, Matrice_Morse& matrice) const
{
  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  int i, j, k, l, f, n, N = ch.valeurs().line_size(), ne_tot = domaine.nb_elem_tot();

  //prerequisites: nu, delta internally + coeffs/delta in Echange_contact boundary conditions
  update_nu(), update_delta();
  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      ref_cast_non_const(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_coeffs(), ref_cast(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_delta();

  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      {
        const Echange_contact_PolyMAC_CDO& cl = ref_cast(Echange_contact_PolyMAC_CDO, cls[i].valeur());
        if (cl.nom_autre_pb() != autre_pb.le_nom())
          continue; //not our problem
        const Front_VF& fvf = ref_cast(Front_VF, cl.frontiere_dis());
        for (j = 0; j < fvf.nb_faces(); j++) //we can fill all coefficients except the face's own (filled by contribuer_a_avec)
          for (k = 0, f = fvf.num_face(j); k < cl.item.dimension(1) && (l = cl.item(j, k)) >= 0; k++)
            for (n = 0; n < N; n++)
              matrice(N * (ne_tot + f) + n, N * l + n) += cl.coeff(j, k + 1, n) - (stab_ ? std::max(delta_f(f, n), cl.delta(j, k, n)) : 0); //operator + nonlinear correction
      }
}

DoubleTab& Op_Diff_PolyMAC_CDO_Elem::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  const IntTab& e_f = domaine.elem_faces();
  const DoubleVect& fs = domaine.face_surfaces(), &ve = domaine.volumes();
  int i, j, e, f, fb, ne_tot = domaine.nb_elem_tot(), n, N = inco.line_size();
  bool elem_only = polymac_flica5 ? resu.dimension_tot(0) == ne_tot : false;
  double fac;

  //prerequisites: nu, delta internally + coeffs/delta in Echange_contact boundary conditions
  update_nu(), update_delta();
  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      ref_cast_non_const(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_coeffs(), ref_cast(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_delta();
  flux_bords_ = 0;

  DoubleTrav nu_ef(e_f.dimension(1), N), mff(N), mfe(N), mee(N);
  for (e = 0; e < ne_tot; e++)
    {
      /* operator: divergence for element rows, continuity for face rows */
      int n_f = domaine.m2d(e + 1) - domaine.m2d(e); //number of faces of element e
      for (remplir_nu_ef(e, nu_ef), mee = 0, i = 0; i < n_f; i++, mee += mfe)
        {
          for (f = e_f(e, i), j = domaine.w2i(domaine.m2d(e) + i), mfe = 0; j < domaine.w2i(domaine.m2d(e) + i + 1); j++, mfe += mff)
            {
              for (fb = e_f(e, domaine.w2j(j)), n = 0, fac = fs(f) * fs(fb) / ve(e) * domaine.w2c(j); n < N; n++)
                mff(n) = fac * nu_ef(domaine.w2j(j), n);
              if (!elem_only)
                for (n = 0; ch.fcl()(f, 0) < 6 && n < N; n++)
                  resu(ne_tot + f, n) -= mff(n) * inco(ne_tot + fb, n);
              if (!elem_only)
                for (n = 0; f < domaine.premiere_face_int() && n < N; n++)
                  flux_bords_(f, n) -= mff(n) * inco(ne_tot + fb, n);
              for (n = 0; n < N; n++)
                resu(e, n) += mff(n) * inco(ne_tot + fb, n);

              //nonlinear correction: "faces/faces" part
              if (!elem_only)
                for (n = 0; stab_ && ch.fcl()(f, 0) < 4 && n < N; n++)
                  resu(ne_tot + f, n) -= std::max(delta_f(f, n), delta_f(fb, n)) * (inco(ne_tot + f, n) - inco(ne_tot + fb, n));
            }
          if (!elem_only)
            for (n = 0; ch.fcl()(f, 0) < 6 && n < N; n++)
              resu(ne_tot + f, n) += mfe(n) * inco(e, n);
          for (n = 0; f < domaine.premiere_face_int() && n < N; n++)
            flux_bords_(f, n) += mfe(n) * inco(e, n);

          //Echange_impose_base
          if (!elem_only)
            if (ch.fcl()(f, 0) > 0 && ch.fcl()(f, 0) < 2 && f < domaine.nb_faces())
              for (n = 0; n < N; n++)
                resu(ne_tot + f, n) -= fs(f) * ref_cast(Echange_impose_base, cls[ch.fcl()(f, 1)].valeur()).h_imp(ch.fcl()(f, 2), n)
                                       * (inco(ch.fcl()(f, 0) == 1 ? ne_tot + f : e, n) - ref_cast(Echange_impose_base, cls[ch.fcl()(f, 1)].valeur()).T_ext(ch.fcl()(f, 2), n));

          //nonlinear correction: "elements/faces" and "faces/elements" parts
          for (n = 0; stab_ && ch.fcl()(f, 0) < 4 && n < N; n++) //not applied to Dirichlet or Neumann boundary conditions
            {
              double corr = std::max(delta_e(e, n), delta_f(f, n)) * (inco(e, n) - inco(ne_tot + f, n));
              resu(e, n) -= corr, resu(ne_tot + f, n) += corr;
            }
        }
      for (n = 0; n < N; n++)
        resu(e, n) -= mee(n) * inco(e, n);
    }

  return resu;
}

void Op_Diff_PolyMAC_CDO_Elem::contribuer_bloc(const DoubleTab& inco, Matrice_Morse& matrice, const int ip) const
{
  if (ip > 3 || ip < -1)
    Process::exit("Op_Diff_PolyMAC_CDO_Elem::contribuer_bloc : invalid bloc number! p must be in [-1, 3]");

  const Champ_Elem_PolyMAC_CDO& ch = ref_cast(Champ_Elem_PolyMAC_CDO, equation().inconnue());
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  const IntTab& e_f = domaine.elem_faces();
  const DoubleVect& fs = domaine.face_surfaces(), &ve = domaine.volumes();
  int i, j, k, l, e, f, fb, ne_tot = domaine.nb_elem_tot(), n, N = inco.line_size();
  double fac;

  //prerequisites: nu, delta internally + coeffs/delta in Echange_contact boundary conditions
  update_nu(), update_delta();
  for (i = 0; i < cls.size(); i++)
    if (sub_type(Echange_contact_PolyMAC_CDO, cls[i].valeur()))
      ref_cast_non_const(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_coeffs(), ref_cast(Echange_contact_PolyMAC_CDO, cls[i].valeur()).update_delta();

  /* operator: divergence for element rows, continuity for face rows */
  DoubleTrav nu_ef(e_f.dimension(1), N), mff(N), mfe(N), mee(N);
  for (e = 0; e < ne_tot; e++)
    {
      int n_f = domaine.m2d(e + 1) - domaine.m2d(e); //number of faces of element e
      for (remplir_nu_ef(e, nu_ef), mee = 0, i = 0; i < n_f; i++, mee += mfe)
        {
          for (f = e_f(e, i), j = domaine.w2i(domaine.m2d(e) + i), mfe = 0; j < domaine.w2i(domaine.m2d(e) + i + 1); j++, mfe += mff)
            {
              for (fb = e_f(e, domaine.w2j(j)), n = 0, fac = fs(f) * fs(fb) / ve(e) * domaine.w2c(j); n < N; n++)
                mff(n) = fac * nu_ef(domaine.w2j(j), n);
              for (n = 0; f < domaine.nb_faces() && ((ch.fcl()(f, 0) < 6 && ch.fcl()(fb, 0) < 6) || ip > -1) && n < N; n++)
                {
                  if (ip == -1)
                    matrice(N * (ne_tot + f) + n, N * (ne_tot + fb) + n) += mff(n);
                  else if (ip == 3)
                    matrice(N * f + n, N * fb + n) += mff(n);
                }
              for (n = 0; e < domaine.nb_elem() && ch.fcl()(fb, 0) < 6 && n < N; n++)
                {
                  if (ip == -1)
                    matrice(N * e + n, N * (ne_tot + fb) + n) -= mff(n);
                  else if (ip == 1)
                    matrice(N * e + n, N * fb + n) -= mff(n);
                }

              //nonlinear correction: "faces/faces" part
              for (n = 0; stab_ && ch.fcl()(f, 0) < 4 && f < domaine.nb_faces() && n < N; n++)
                for (k = 0, fac = std::max(delta_f(f, n), delta_f(fb, n)); k < 2; k++)
                  {
                    if (ip == -1)
                      matrice(N * (ne_tot + f) + n, N * (ne_tot + (k ? fb : f)) + n) += (k ? -1 : 1) * fac;
                    else if (ip == 3)
                      matrice(N * f + n, N * (k ? fb : f) + n) += (k ? -1 : 1) * fac;
                  }
            }
          for (n = 0; f < domaine.nb_faces() && (ch.fcl()(f, 0) < 6 || ip > -1) && n < N; n++)
            {
              if (ip == -1)
                matrice(N * (ne_tot + f) + n, N * e + n) -= mfe(n);
              else if (ip == 2)
                matrice(N * f + n, N * e + n) -= mfe(n);
            }

          //Echange_impose_base
          if (ch.fcl()(f, 0) == 1 && f < domaine.nb_faces())
            for (n = 0; n < N; n++)
              {
                if (ip == -1)
                  matrice(N * (ne_tot + f) + n, N * (ch.fcl()(f, 0) == 1 ? ne_tot + f : e) + n) += fs(f) * ref_cast(Echange_impose_base, cls[ch.fcl()(f, 1)].valeur()).h_imp(ch.fcl()(f, 2), n);
                else
                  Process::exit("Echange_impose_base et diffusion explicite : pas bon!");
              }
          else if (ch.fcl()(f, 0) == 3 && f < domaine.nb_faces()) //paroi_contact handled monolithically -> add coefficient to face from the other side
            {
              const Echange_contact_PolyMAC_CDO& cl = ref_cast(Echange_contact_PolyMAC_CDO, cls[ch.fcl()(f, 1)].valeur());
              for (j = ch.fcl()(f, 2), n = 0; n < N; n++)
                {
                  if (ip == -1)
                    matrice(N * (ne_tot + f) + n, N * (ne_tot + f) + n) += cl.coeff(j, 0, n); //coefficient of the face itself
                  else
                    Process::exit("Echange_impose_base et diffusion explicite : pas bon!");
                }
              for (k = 0; stab_ && k < cl.item.dimension(1) && cl.item(j, k) >= 0; k++)
                for (n = 0; n < N; n++) //nonlinear correction
                  {
                    if (ip == -1)
                      matrice(N * (ne_tot + f) + n, N * (ne_tot + f) + n) += std::max(delta_f(f, n), cl.delta(j, k, n));
                    else
                      Process::exit("Echange_impose_base et diffusion explicite : pas bon!");
                  }
            }

          //nonlinear correction: "elements/faces" and "faces/elements" parts
          for (n = 0; stab_ && ch.fcl()(f, 0) < 4 && n < N; n++) //not applied to Dirichlet or Neumann boundary conditions
            {
              double corr = std::max(delta_e(e, n), delta_f(f, n));
              for (k = 0; k < 2; k++)
                for (l = 0; (k ? (f < domaine.nb_faces()) : (e < domaine.nb_elem())) && l < 2; l++)
                  {
                    if (ip == -1)
                      matrice(N * (k ? ne_tot + f : e) + n, N * (l ? ne_tot + f : e) + n) += (k == l ? 1 : -1) * corr;
                    else
                      Process::exit("Echange_impose_base et diffusion explicite : pas bon!");
                  }
            }
        }
      for (n = 0; e < domaine.nb_elem() && n < N; n++)
        if (ip == -1 || (ip == 0))
          matrice(N * e + n, N * e + n) += mee(n);
    }

}

void Op_Diff_PolyMAC_CDO_Elem::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const int ne_tot = domaine.nb_elem_tot(), nf_tot = domaine.nb_faces_tot();
  int i = -1;
  if (matrice.nb_lignes() == ne_tot && matrice.nb_colonnes() == ne_tot) i = 0;
  if (matrice.nb_lignes() == ne_tot && matrice.nb_colonnes() == nf_tot) i = 1;
  if (matrice.nb_lignes() == nf_tot && matrice.nb_colonnes() == ne_tot) i = 2;
  if (matrice.nb_lignes() == nf_tot && matrice.nb_colonnes() == nf_tot) i = 3;
  contribuer_bloc(inco, matrice, i);
}

void Op_Diff_PolyMAC_CDO_Elem::update_auxiliary_variables()
{
  update_auxiliary_variables(le_champ_inco->valeurs());
}

static Matrice_Morse FE, FF;
void Op_Diff_PolyMAC_CDO_Elem::update_auxiliary_variables(DoubleTab& inco)
{
  if (!polymac_flica5) return;
// solve M_ff T_f = -M_fe T_e
  DoubleTab sm(inco);
  sm = 0.;
  DoubleTab_parts inco_parts(inco);
  DoubleTab_parts sm_parts(sm);

  // 1. sizing
  //std::clock_t start = std::clock();
  if (FF.nb_lignes() == 0)
    {
      //Matrice_Morse FE, FF;
      dimensionner_bloc(FE, 2);
      dimensionner_bloc(FF, 3);
    }
  else
    {
      FE.clean();
      FF.clean();
    }
  //Cout << "[Op_Diff_PolyMAC_CDO_Elem] Time to initialize matrix: " << (std::clock() - start) / (double) CLOCKS_PER_SEC << finl;
  // 2. filling
  contribuer_a_avec(inco, FE);
  contribuer_a_avec(inco, FF);

  // 3. solving
  FE *= -1.;
  FE.ajouter_multvect(inco_parts[0], sm_parts[1]);

  if (FF.is_diagonal())
    {
      const int n = FF.nb_lignes();
      for (int i = 0; i < n; i++)
        for (auto k = FF.get_tab1()(i) - 1; k < FF.get_tab1()(i + 1) - 1; k++)
          inco_parts[1][i] = sm_parts[1][i] / FF.get_coeff()(k);
      inco_parts[1].echange_espace_virtuel();
    }
  else
    {
      set_solveur()->reinit();
      inco_parts[1] = 0.;
      try
        {
          set_solveur().resoudre_systeme(FF, sm_parts[1], inco_parts[1]);
        }
      catch(...)
        {
          // PL: Cholesky_lapack crash on meshes with no local diffusion:
          statistics().end_count(STD_COUNTERS::system_solver,0,0);
          Nom solv("Petsc Cholesky { quiet }");
          Cerr << "Echec du solveur dans Op_Diff_PolyMAC_CDO_Elem..." << finl;
          Cerr << "On change pour un solveur plus robuste (mais plus lent): " << solv << finl;
          EChaine chl(solv);
          lire_solveur(chl);
          solveur.nommer("Op_Diff_PolyMAC_CDO_Elem solver");
          set_solveur().resoudre_systeme(FF, sm_parts[1], inco_parts[1]);
        }
    }
}
