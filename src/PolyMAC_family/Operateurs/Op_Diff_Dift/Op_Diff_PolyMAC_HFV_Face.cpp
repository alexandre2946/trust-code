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

#include <Op_Diff_PolyMAC_HFV_Face.h>
#include <Linear_algebra_tools_impl.h>
#include <Champ_Face_PolyMAC_HFV.h>
#include <Domaine_PolyMAC_HFV.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Pb_Multiphase.h>
#include <Synonyme_info.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Perf_counters.h>


Implemente_instanciable( Op_Diff_PolyMAC_HFV_Face, "Op_Diff_PolyMAC_HFV_Face|Op_Dift_PolyMAC_HFV_Face_PolyMAC_HFV", Op_Diff_PolyMAC_HFV_base );
Add_synonym(Op_Diff_PolyMAC_HFV_Face, "Op_Diff_PolyMAC_HFV_var_Face");
Add_synonym(Op_Diff_PolyMAC_HFV_Face, "Op_Dift_PolyMAC_HFV_var_Face_PolyMAC_HFV");

Sortie& Op_Diff_PolyMAC_HFV_Face::printOn(Sortie& os) const { return Op_Diff_PolyMAC_HFV_base::printOn(os); }

Entree& Op_Diff_PolyMAC_HFV_Face::readOn(Entree& is) { return Op_Diff_PolyMAC_HFV_base::readOn(is); }

void Op_Diff_PolyMAC_HFV_Face::completer()
{
  Op_Diff_PolyMAC_HFV_base::completer();
  const Domaine_PolyMAC_HFV& domaine = ref_cast(Domaine_PolyMAC_HFV, le_dom_poly_.valeur());
  Equation_base& eq = equation();
  Champ_Face_PolyMAC_HFV& ch = ref_cast(Champ_Face_PolyMAC_HFV, le_champ_inco ? le_champ_inco.valeur() : eq.inconnue());
  ch.init_auxiliary_variables(); /* add auxiliary unknowns (vorticities at edges) */
  flux_bords_.resize(domaine.premiere_face_int(), dimension * ch.valeurs().line_size());
  if (domaine.domaine().nb_joints() && domaine.domaine().joint(0).epaisseur() < 1)
    {
      Cerr << "Op_Diff_PolyMAC_HFV_Face : largeur de joint insuffisante (minimum 1)!" << finl;
      Process::exit();
    }
  porosite_e.ref(equation().milieu().porosite_elem());
  porosite_f.ref(equation().milieu().porosite_face());
  op_ext = { this };
}

double Op_Diff_PolyMAC_HFV_Face::calculer_dt_stab() const
{
  const Domaine_PolyMAC_HFV& domaine = ref_cast(Domaine_PolyMAC_HFV, le_dom_poly_.valeur());
  const IntTab& e_f = domaine.elem_faces();

  const DoubleTab& nf = domaine.face_normales(),
                   *alp = sub_type(Pb_Multiphase, equation().probleme()) ?
                          &ref_cast(Pb_Multiphase, equation().probleme()).equation_masse().inconnue().passe() : nullptr,
                          *a_r = sub_type(Pb_Multiphase, equation().probleme()) ?
                                 &ref_cast(Pb_Multiphase, equation().probleme()).equation_masse().champ_conserve().passe() :
                                 (has_champ_masse_volumique() ? &get_champ_masse_volumique().valeurs() : nullptr); /* alpha * rho product */

  const DoubleVect& pe = equation().milieu().porosite_elem(), &vf = domaine.volumes_entrelaces(), &ve = domaine.volumes();
  update_nu();

  const int N = equation().inconnue().valeurs().line_size();
  double dt = 1e10;
  DoubleTrav flux(N);

  for (int e = 0; e < domaine.nb_elem(); e++)
    {
      flux = 0.;
      const double vol = pe(e) * ve(e);

      for (int i = 0; i < e_f.dimension(1); i++)
        {
          const int f = e_f(e, i);
          if (f < 0) continue;

          for (int n = 0; n < N; n++)
            flux(n) += domaine.nu_dot(&nu_, e, n, &nf(f, 0), &nf(f, 0)) / vf(f);
        }

      for (int n = 0; n < N; n++)
        if ((!alp || (*alp)(e, n) > 0.25) && flux(n)) /* below 0.5e-6, assume evanescence handles it */
          dt = std::min(dt, vol * (a_r ? (*a_r)(e, n) : 1) / flux(n));
    }
  return Process::mp_min(dt);
}

void Op_Diff_PolyMAC_HFV_Face::dimensionner_blocs_ext(int aux_only, matrices_t matrices, const tabs_t& semi_impl) const
{
  const Champ_Face_PolyMAC_HFV& ch = ref_cast(Champ_Face_PolyMAC_HFV, le_champ_inco ? le_champ_inco.valeur() : equation().inconnue());
  const std::string& nom_inco = ch.le_nom().getString();
  if (!matrices.count(nom_inco))
    return; //no diagonal block -> nothing to do

  const Domaine_PolyMAC_HFV& domaine = ref_cast(Domaine_PolyMAC_HFV, le_dom_poly_.valeur());
  const IntTab& e_f = domaine.elem_faces(), &f_s = domaine.face_sommets(), &e_a = domaine.domaine().elem_aretes(), &fcl = ch.fcl();

  Matrice_Morse& mat = *matrices.at(nom_inco), mat2;

  update_nu();

  ConstDoubleTab_parts p_inco(ch.valeurs());

  const int N = ch.valeurs().line_size(), nf_tot = domaine.nb_faces_tot(), D = dimension,
            N_nu = nu_.line_size(), semi = (int) semi_impl.count(nom_inco);

  Stencil stencil(0, 2);

  Cerr << "Op_Diff_PolyMAC_HFV_Face::dimensionner() : ";

  /* block (faces, edges): rot [(lambda grad)^u]*/
  if (!semi && !aux_only)
    for (int f = 0; f < domaine.nb_faces(); f++)
      for (int i = 0; i < f_s.dimension(1); i++)
        {
          const int s = f_s(f, i);
          if (s < 0) continue;

          const int a = D < 3 ? s :
                        domaine.som_arete[s].at(f_s(f, i + 1 < f_s.dimension(1) && f_s(f, i + 1) >= 0 ? i + 1 : 0)); //edge index

          for (int n = 0; n < N; n++)
            stencil.append_line(N * f + n, N * (nf_tot + a) + n);
        }

  /* blocks (edges, faces) and (edges, edges): using M2 and W1 in each element */
  Matrice33 L(0, 0, 0, 0, 0, 0, 0, 0, D < 3), iL; //diffusion tensor in each element, its inverse

  DoubleTrav inu, m2, w1, v_e, v_ea; //in the format expected by domaine.nu_dot

  nu_.nb_dim() == 2 ? inu.resize(1, N) :
  nu_.nb_dim() == 3 ? inu.resize(1, N, D) :
  inu.resize(1, N, D, D);

  if (!semi)
    for (int e = 0; e < domaine.nb_elem_tot(); e++)
      {
        //diagonal or anisotropic diagonal diffusion tensor: easy inversion!
        if (nu_.nb_dim() < 4)
          {
            for (int i = 0; i < N_nu; i++)
              inu.addr()[i] = 1. / nu_.addr()[N_nu * e + i];
          }
        else
          {
            for (int n = 0; n < N; n++) //otherwise: one matrix to invert per component
              {
                for (int d = 0; d < D; d++)
                  for (int db = 0; db < D; db++)
                    L(d, db) = nu_(e, n, d, db);

                Matrice33::inverse(L, iL);

                for (int d = 0; d < D; d++)
                  for (int db = 0; db < D; db++)
                    inu(0, n, d, db) = iL(d, db);
              }
          }

        domaine.M2(&inu, e, m2);

        if (D > 2)
          domaine.W1(&nu_, e, w1, v_e, v_ea); //only in 3D: in 2D, diagonal matrix

        //block (edges, faces): by iterating over the edges of each face
        if (!aux_only)
          for (int i = 0; i < m2.dimension(0); i++)
            {
              const int f = e_f(e, i);

              for (int j = 0; j < f_s.dimension(1); j++)
                {
                  const int s = f_s(f, j);
                  if (s < 0) continue;

                  const int a = D < 3 ? s : domaine.som_arete[s].at(f_s(f, j + 1 < f_s.dimension(1) && f_s(f, j + 1) >= 0 ? j + 1 : 0)); //edge index

                  if (a < (D < 3 ? domaine.domaine().nb_som() : domaine.domaine().nb_aretes()))
                    for (int k = 0; k < m2.dimension(1); k++)
                      {
                        const int fb = e_f(e, k);
                        for (int n = 0; n < N; n++)
                          if (fcl(f, 0) == 2 || m2(i, k, n)) //if f is Symmetry, there is also a ve contribution -> full dependence
                            stencil.append_line(N * (nf_tot + a) + n, N * fb + n);
                      }
                }
            }

        //block (edges, edges): with m1 if D = 3 (otherwise handled below)
        if (D > 2)
          for (int i = 0; i < w1.dimension(0); i++)
            {
              const int a = e_a(e, i);

              if (a < domaine.domaine().nb_aretes())
                for (int j = 0; j < w1.dimension(1); j++)
                  {
                    const int ab = e_a(e, j);
                    for (int n = 0; n < N; n++)
                      if (w1(i, j, n))
                        stencil.append_line(N * (!aux_only * nf_tot + a) + n, N * (!aux_only * nf_tot + ab) + n);
                  }
            }
      }
  if (semi || D < 3)
    for (int s = 0; s < (D < 3 ? domaine.nb_som() : domaine.domaine().nb_aretes()); s++)
      for (int n = 0; n < N; n++)
        stencil.append_line(N * (!aux_only * nf_tot + s) + n, N * (!aux_only * nf_tot + s) + n);

  tableau_trier_retirer_doublons(stencil);

  Cerr << "OK" << finl;

  Matrix_tools::allocate_morse_matrix(aux_only ? p_inco[1].size_totale() : ch.valeurs().size_totale(), aux_only ? p_inco[1].size_totale() : ch.valeurs().size_totale(), stencil, mat2);

  if (mat.nb_colonnes())
    mat += mat2;
  else
    mat = mat2;
}

// adds the diffusion contribution to the right-hand side resu
// returns resu
void Op_Diff_PolyMAC_HFV_Face::ajouter_blocs_ext(int aux_only, matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_Face_PolyMAC_HFV& ch = ref_cast(Champ_Face_PolyMAC_HFV, le_champ_inco ? le_champ_inco.valeur() : equation().inconnue());
  const Conds_lim& cls = ch.domaine_Cl_dis().les_conditions_limites();
  const Domaine_PolyMAC_HFV& domaine = ref_cast(Domaine_PolyMAC_HFV, le_dom_poly_.valeur());
  const IntTab& e_f = domaine.elem_faces(), &f_s = domaine.face_sommets(), &e_a = domaine.domaine().elem_aretes(),
                &fcl = ch.fcl(), &f_e = domaine.face_voisins();

  const std::string& nom_inco = ch.le_nom().getString();
  Matrice_Morse *mat = matrices.count(nom_inco) ? matrices.at(nom_inco) : nullptr;

  update_nu();

  const int N = ch.valeurs().line_size(), nf_tot = domaine.nb_faces_tot(),  D = dimension,
            N_nu = nu_.line_size(), semi = (int) semi_impl.count(nom_inco);

  double vecz[3] = { 0, 0, 1 }, v_cl[3];
  const double t = equation().schema_temps().temps_courant();

  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv(), &xs = domaine.domaine().coord_sommets(),
                   &xa = D < 3 ? xs : domaine.xa(), &ta = domaine.ta(), &nf = domaine.face_normales(),
                    &inco = semi_impl.count(nom_inco) ? semi_impl.at(nom_inco) : ch.valeurs();

  const DoubleVect& la = domaine.longueur_aretes(), &vf = domaine.volumes_entrelaces(), &fs = domaine.face_surfaces(), &ve = domaine.volumes();

  /* what to do with auxiliary variables? */
  if (aux_only)  /* 1) assembling the auxiliary variable resolution system itself */
    use_aux_ = 0;
  else if (mat && !semi) /* 2) fully implicit: no need for mat_aux / var_aux */
    {
      t_last_aux_ = t;
      use_aux_ = 0;
    }
  else if (t_last_aux_ < t) /* 3) first step at this time in semi-implicit: compute auxiliary variables at t and store in var_aux */
    update_aux(t);

  ConstDoubleTab_parts p_inco(inco); /* two parts of the unknown */

  const DoubleTab& omega = use_aux_ ? var_aux : p_inco[1]; /* auxiliary variables can be either in inco/semi_impl (case 1) or in var_aux (case 2) */

  /* block (faces, edges): rot [(lambda grad)^u]*/
  if (!aux_only)
    for (int f = 0; f < domaine.nb_faces(); f++)
      for (int i = 0; i < f_s.dimension(1); i++)
        {
          const int s = f_s(f, i);
          if (s < 0) continue;

          const int a = D < 3 ? s : domaine.som_arete[s].at(f_s(f, i + 1 < f_s.dimension(1) && f_s(f, i + 1) >= 0 ? i + 1 : 0)); //edge index

          auto vec = domaine.cross(3, D, D < 3 ? vecz : &ta(a, 0), &xv(f, 0), nullptr, &xa(a, 0));

          const int sgn = domaine.dot(&nf(f, 0), &vec[0]) > 0 ? 1 : -1; //edge-face orientation

          for (int n = 0; n < N; n++)
            secmem(f, n) -= sgn * vf(f) / fs(f) * (D < 3 ? 1 : la(a)) * omega(a, n);

          if (mat && !semi)
            for (int n = 0; n < N; n++)
              (*mat)(N * f + n, N * (nf_tot + a) + n) += sgn * vf(f) / fs(f) * (D < 3 ? 1 : la(a));
        }

  /* blocks (edges, faces) and (edges, edges): using M2 and W1 in each element */
  Matrice33 L(0, 0, 0, 0, 0, 0, 0, 0, D < 3), iL; //diffusion tensor in each element, its inverse and its square
  DoubleTrav dL(N), inu, m2, w1, v_e, v_ea; //determinant, inverse (in the format expected by domaine.nu_dot), matrices M2(iL) / W1(L)

  if (nu_.nb_dim() == 2)
    inu.resize(1, N);
  else if (nu_.nb_dim() == 3)
    inu.resize(1, N, D);
  else
    inu.resize(1, N, D, D);

  if (!aux_only && mat && semi)
    {
      for (int a = 0; a < xa.dimension(0); a++)
        for (int n = 0; n < N; n++) /* semi-implicit: equalities w_a^+ = var_aux */
          {
            secmem(nf_tot + a, n) += omega(a, n) - ch.valeurs()(nf_tot + a, n);
            (*mat)(N * (nf_tot + a) + n, N * (nf_tot + a) + n)++;
          }
    }
  else if (mat && !semi)
    for (int e = 0; e < domaine.nb_elem_tot(); e++) /* implicit: true equations */
      {
        //diagonal or anisotropic diagonal diffusion tensor: easy inversions!
        if (nu_.nb_dim() == 2)
          {
            for (int n = 0; n < N; n++)
              {
                inu(0, n) = 1. / nu_(e, n);
                dL(n) = std::pow(nu_(e, n), D);
              }
          }
        else if (nu_.nb_dim() == 3)
          {
            for (int n = 0; n < N; n++)
              {
                dL(n) = 1.;
                for (int d = 0; d < D; d++)
                  {
                    inu(0, n, d) = 1. / nu_(e, n, d);
                    dL(n) *= nu_(e, n, d);
                  }
              }
          }

        if (nu_.nb_dim() < 4)
          {
            for (int i = 0; i < N_nu; i++)
              inu.addr()[i] = 1. / nu_.addr()[N_nu * e + i];
          }
        else
          {
            for (int n = 0; n < N; n++) //otherwise: one matrix to invert per component
              {
                for (int d = 0; d < D; d++)
                  for (int db = 0; db < D; db++)
                    L(d, db) = nu_(e, n, d, db);

                dL(n) = Matrice33::inverse(L, iL); //returns the determinant!

                for (int d = 0; d < D; d++)
                  for (int db = 0; db < D; db++)
                    inu(0, n, d, db) = iL(d, db);
              }
          }

        domaine.M2(&inu, e, m2);
        if (D > 2)
          domaine.W1(&nu_, e, w1, v_e, v_ea); //only in 3D: in 2D, diagonal matrix

        //block (edges, faces): by iterating over the edges of each face
        for (int i = 0; i < m2.dimension(0); i++)
          {
            const int f = e_f(e, i);
            for (int j = 0; j < f_s.dimension(1); j++)
              {
                const int s = f_s(f, j);
                if (s < 0) continue;

                const int a = D < 3 ? s : domaine.som_arete[s].at(f_s(f, j + 1 < f_s.dimension(1) && f_s(f, j + 1) >= 0 ? j + 1 : 0)); // edge index

                if (a < (D < 3 ? domaine.domaine().nb_som() : domaine.domaine().nb_aretes()))
                  {
                    auto vec = domaine.cross(3, D, D < 3 ? vecz : &ta(a, 0), &xv(f, 0), nullptr, &xa(a, 0));

                    const int sgn = (e == f_e(f, 0) ? 1 : -1) * domaine.dot(&nf(f, 0), &vec[0]) > 0 ? 1 : -1; //edge-face orientation (in the outward direction of e)

                    for (int k = 0; k < m2.dimension(1); k++)
                      {
                        const int fb = e_f(e, k);
                        for (int n = 0; n < N; n++) //part x_e -> x_f with m2 + part x_f -> x_a with v_e if Neumann/Symmetry boundary
                          {
                            const double coeff = m2(i, k, n) + (fcl(f, 0) == 2 ? domaine.nu_dot(&inu, 0, n, &xa(a, 0), &xv(fb, 0), &xv(f, 0), &xp(e, 0)) / ve(e) : 0);

                            if (!coeff)
                              continue;

                            secmem(!aux_only * nf_tot + a, n) += sgn * coeff * inco(fb, n) * fs(fb) * (e == f_e(fb, 0) ? 1 : -1);

                            if (!aux_only)
                              (*mat)(N * (nf_tot + a) + n, N * fb + n) -= sgn * coeff * fs(fb) * (e == f_e(fb, 0) ? 1 : -1);
                          }
                      }

                    if (fcl(f, 0) == 3 && sub_type(Dirichlet, cls[fcl(f, 1)].valeur()))
                      for (int n = 0; n < N; n++) //if Dirichlet boundary: part x_f -> x_a with the velocity given by the BC
                        {
                          for (int d = 0; d < D; d++)
                            v_cl[d] = ref_cast(Dirichlet, cls[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), N * d + n); //v impose

                          secmem(!aux_only * nf_tot + a, n) += sgn * domaine.nu_dot(&inu, 0, n, &xa(a, 0), v_cl, &xv(f, 0));
                        }
                  }
              }
          }

        //block (edges, edges): using m1 if D = 3, diagonal * surface if D = 2
        if (D == 2)
          {
            for (int i = 0; i < e_f.dimension(1); i++)
              {
                const int f = e_f(e, i);
                if (f < 0) continue;

                for (int j = 0; j < 2; j++)
                  {
                    const int s = f_s(f, j);
                    if (s < domaine.nb_som())
                      {
                        auto vec = domaine.cross(D, D, &xv(f, 0), &xs(s, 0), &xp(e, 0), &xp(e, 0));

                        const double surf = std::abs(vec[2]) / 2.; //area of triangle (e, f, s)

                        for (int n = 0; n < N; n++)
                          secmem(!aux_only * nf_tot + s, n) -= surf * inco(nf_tot + s, n) / dL(n);

                        for (int n = 0; n < N; n++)
                          (*mat)(N * (!aux_only * nf_tot + s) + n, N * (!aux_only * nf_tot + s) + n) += surf / dL(n);
                      }
                  }
              }
          }
        else
          for (int i = 0; i < w1.dimension(0); i++)
            {
              const int a = e_a(e, i);

              if (a < domaine.domaine().nb_aretes())
                for (int j = 0; j < w1.dimension(1); j++)
                  {
                    const int ab = e_a(e, j);
                    for (int n = 0; n < N; n++)
                      if (w1(i, j, n))
                        {
                          secmem(!aux_only * nf_tot + a, n) -= w1(i, j, n) * la(ab) * inco(nf_tot + ab, n) / dL(n);
                          (*mat)(N * (!aux_only * nf_tot + a) + n, N * (!aux_only * nf_tot + ab) + n) += w1(i, j, n) * la(ab) / dL(n);
                        }
                  }
            }
      }
}

