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

#include <Op_Conv_EF_Stab_PolyMAC_HFV_Face.h>
#include <Champ_Face_PolyMAC_HFV.h>
#include <Dirichlet_homogene.h>
#include <Masse_ajoutee_base.h>
#include <Schema_Temps_base.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Domaine_Poly_base.h>
#include <Pb_Multiphase.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTLists.h>
#include <Dirichlet.h>
#include <Param.h>
#include <cmath>

Implemente_instanciable( Op_Conv_EF_Stab_PolyMAC_HFV_Face, "Op_Conv_EF_Stab_PolyMAC_HFV_Face", Op_Conv_EF_Stab_PolyMAC_CDO_Face );
Implemente_instanciable( Op_Conv_Amont_PolyMAC_HFV_Face, "Op_Conv_Amont_PolyMAC_HFV_Face", Op_Conv_EF_Stab_PolyMAC_HFV_Face );
Implemente_instanciable( Op_Conv_Centre_PolyMAC_HFV_Face, "Op_Conv_Centre_PolyMAC_HFV_Face", Op_Conv_EF_Stab_PolyMAC_HFV_Face );

// XD Op_Conv_EF_Stab_PolyMAC_HFV_Face interprete Op_Conv_EF_Stab_PolyMAC_HFV_Face BRACE Class
// XD_CONT Op_Conv_EF_Stab_PolyMAC_HFV_Face

Sortie& Op_Conv_EF_Stab_PolyMAC_HFV_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_CDO_base::printOn(os); }
Sortie& Op_Conv_Amont_PolyMAC_HFV_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_CDO_base::printOn(os); }
Sortie& Op_Conv_Centre_PolyMAC_HFV_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_CDO_base::printOn(os); }

Entree& Op_Conv_EF_Stab_PolyMAC_HFV_Face::readOn(Entree& is) { return Op_Conv_EF_Stab_PolyMAC_CDO_Face::readOn(is); }

Entree& Op_Conv_Amont_PolyMAC_HFV_Face::readOn(Entree& is)
{
  alpha_ = 1.0;
  return Op_Conv_PolyMAC_CDO_base::readOn(is);
}

Entree& Op_Conv_Centre_PolyMAC_HFV_Face::readOn(Entree& is)
{
  alpha_ = 0.0;
  return Op_Conv_PolyMAC_CDO_base::readOn(is);
}

double Op_Conv_EF_Stab_PolyMAC_HFV_Face::calculer_dt_stab() const
{
  double dt = 1e10;
  const Domaine_Poly_base& domaine = le_dom_poly_.valeur();
  const DoubleVect& fs = domaine.face_surfaces(), &pf = equation().milieu().porosite_face(), &ve = domaine.volumes(), &pe = equation().milieu().porosite_elem();
  const DoubleTab& vit = vitesse_->valeurs(),
                   *alp = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()).equation_masse().inconnue().passe() : nullptr;
  const IntTab& e_f = domaine.elem_faces(), &f_e = domaine.face_voisins();
  const int N = vit.line_size();
  DoubleTrav flux(N); //sum of fluxes pf * |f| * vf, minimum volume of elements/faces affected by this flux

  for (int e = 0; e < domaine.nb_elem(); e++)
    {
      // Compute the effective volume of the element
      const double vol = pe(e) * ve(e);
      flux = 0.;

      // Loop over faces associated with the element
      for (int i = 0; i < e_f.dimension(1); i++)
        {
          int f = e_f(e, i);
          if (f < 0) continue; // non-existent face

          for (int n = 0; n < N; n++)
            {
              // Add incoming flux for component n: only incoming fluxes count
              double flux_f = pf(f) * fs(f) * std::max((e == f_e(f, 1) ? 1 : -1) * vit(f, n), 0.);
              flux(n) += flux_f;
            }
        }

      // Compute the time step for each component n
      for (int n = 0; n < N; n++)
        if ((!alp || (*alp)(e, n) > 1e-3) && std::abs(flux(n)) > 1e-12)
          dt = std::min(dt, vol / flux(n));
    }

  return Process::mp_min(dt);
}

void Op_Conv_EF_Stab_PolyMAC_HFV_Face::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_Poly_base& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_HFV& ch = ref_cast(Champ_Face_PolyMAC_HFV, equation().inconnue());
  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces(), &fcl = ch.fcl(), &equiv = domaine.equiv();
  const DoubleTab& nf = domaine.face_normales(), &inco = ch.valeurs(), &xp = domaine.xp(), &xv = domaine.xv();
  const DoubleVect& fs = domaine.face_surfaces(), &ve = domaine.volumes();

  const std::string& nom_inco = ch.le_nom().getString();

  if (!matrices.count(nom_inco) || semi_impl.count(nom_inco))
    return; //no diagonal block or semi-implicit -> nothing to do

  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;
  const Masse_ajoutee_base *corr = pbm && pbm->has_correlation("masse_ajoutee") ? &ref_cast(Masse_ajoutee_base, pbm->get_correlation("masse_ajoutee")) : nullptr;
  Matrice_Morse& mat = *matrices.at(nom_inco), mat2;

  const int N = equation().inconnue().valeurs().line_size();

  Stencil stencil(0, 2);

  /* This block acts on elements only; the diagonal of the matrix is omitted. */
  for (int f = 0; f < domaine.nb_faces_tot(); f++)
    {
      // Check if the face is internal or satisfies boundary conditions
      if (f_e(f, 0) >= 0 && (f_e(f, 1) >= 0 || fcl(f, 0) == 3))
        {
          // Loop over elements associated with this face
          for (int i = 0; i < 2 ; i++)
            {
              const int e = f_e(f, i);
              if (e < 0) continue; // virtual elem

              for (int j = 0; j < 2 ; j++)
                {
                  const int eb = f_e(f, j);
                  // Loop over faces connected to the current element
                  if (eb < 0) continue;

                  for (int k = 0; k < e_f.dimension(1); k++)
                    {
                      const int fb = e_f(e, k);
                      if (fb <0) continue;

                      if (fb < domaine.nb_faces())
                        {
                          int fc = equiv(f, i, k);
                          // Case where a face equivalence exists
                          if (fc >= 0)
                            {
                              for (int n = 0; n < N; n++)
                                for (int m = (corr ? 0 : n); m < (corr ? N : n + 1); m++)
                                  stencil.append_line(N * fb + n, N * fc + m);
                            }
                          // Case without equivalence: contributions between element faces
                          else if (f_e(f, 1) >= 0)
                            {
                              for (int l = 0; l < e_f.dimension(1); l++)
                                {
                                  fc = e_f(eb, l);
                                  if (fc < 0) continue;

                                  const double critere = fs(fc) * domaine.dot(&xv(fc, 0), &nf(fb, 0), &xp(eb, 0));
                                  if (std::abs(critere) > 1e-6 * ve(eb) * fs(fb))
                                    for (int n = 0; n < N; n++)
                                      for (int m = (corr ? 0 : n); m < (corr ? N : n + 1); m++)
                                        stencil.append_line(N * fb + n, N * fc + m);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

  // Sort and remove duplicates from the stencil
  tableau_trier_retirer_doublons(stencil);

  // Allocate a sparse matrix based on the stencil
  Matrix_tools::allocate_morse_matrix(inco.size_totale(), inco.size_totale(), stencil, mat2);

  // Add mat2 to the existing matrix or initialize 'mat'
  if (mat.nb_colonnes())
    mat += mat2;
  else
    mat = mat2;
}

// adds the convection contribution to the right-hand side resu
// returns resu
void Op_Conv_EF_Stab_PolyMAC_HFV_Face::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_Poly_base& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_HFV& ch = ref_cast(Champ_Face_PolyMAC_HFV, equation().inconnue());
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();
  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces(), &fcl = ch.fcl(), &equiv = domaine.equiv();
  const DoubleTab& vit = ch.passe(), &nf = domaine.face_normales(), &vfd = domaine.volumes_entrelaces_dir(), &xp = domaine.xp(), &xv = domaine.xv();
  const DoubleVect& fs = domaine.face_surfaces(), &pe = porosite_e, &pf = porosite_f, &ve = domaine.volumes();

  /* a_r : alpha_rho product if Pb_Multiphase -> via semi-implicit, or by retrieving the conserved field from the mass equation */
  const std::string& nom_inco = ch.le_nom().getString();
  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr;
  const Masse_ajoutee_base *corr = pbm && pbm->has_correlation("masse_ajoutee") ? &ref_cast(Masse_ajoutee_base, pbm->get_correlation("masse_ajoutee")) : nullptr;
  const DoubleTab& inco = semi_impl.count(nom_inco) ? semi_impl.at(nom_inco) : ch.valeurs(),
                   *a_r = !pbm ? nullptr : semi_impl.count("alpha_rho") ? &semi_impl.at("alpha_rho") : &pbm->equation_masse().champ_conserve().valeurs(),
                    *alp = pbm ? &pbm->equation_masse().inconnue().passe() : nullptr,
                     &rho = equation().milieu().masse_volumique().passe();

  Matrice_Morse *mat = matrices.count(nom_inco) && !semi_impl.count(nom_inco) ? matrices.at(nom_inco) : nullptr;

  const int N = inco.line_size(), D = dimension;

  DoubleTrav dfac(2, N, N), masse(N, N);

  // Loop over all faces of the domain
  for (int f = 0; f < domaine.nb_faces_tot(); f++)
    {
      if (f_e(f, 0) >= 0 && (f_e(f, 1) >= 0 || fcl(f, 0) == 1 || fcl(f, 0) == 3))
        {
          // Compute face contributions
          dfac = 0.;
          for (int i = 0; i < 2; i++)
            {
              // Diagonal mass with correction if needed
              masse = 0.;
              int e = f_e(f, (f_e(f, i) >= 0) ? i : 0);
              for (int n = 0; n < N; n++)
                masse(n, n) = a_r ? (*a_r)(e, n) : 1.;

              if (corr)
                corr->ajouter(&(*alp)(e, 0), &rho(e, 0), masse);

              // Contribution a dfac
              e = f_e(f, i);
              int eb = f_e(f, i);
              for (int n = 0; n < N; n++)
                for (int m = 0; m < N; m++)
                  {
                    const double signe = (vit(f, m) * (i ? -1 : 1) >= 0) ? 1. : (vit(f, m) ? -1. : 0.);
                    dfac((fcl(f, 0) == 1) ? 0 : i, n, m) += fs(f) * vit(f, m) * pe((eb >= 0) ? eb : f_e(f, 0)) * masse(n, m) * (1. + signe * alpha_) / 2;
                  }
            }

          // Contributions to matrices and right-hand side
          for (int i = 0; i < 2 ; i++)
            {
              const int e = f_e(f, i);
              if (e < 0) continue;

              for (int k = 0; k < e_f.dimension(1) ; k++)
                {
                  const int fb = e_f(e, k);
                  if (fb < 0) continue;

                  if (fb < domaine.nb_faces())
                    {
                      int fc = equiv(f, i, k);
                      // Equivalence case: source face -> target face
                      if (fc >= 0 || f_e(f, 1) < 0)
                        {
                          for (int j = 0; j < 2; j++)
                            {
                              int eb = f_e(f, j);
                              int fd = (j == i) ? fb : fc; // source face or element

                              //multiplier to convert from vf to ve
                              double mult = (fd < 0 || domaine.dot(&nf(fb, 0), &nf(fd, 0)) > 0) ? 1 : -1;
                              mult *= (fd >= 0) ? pf(fd) / pe(eb) : 1;

                              for (int n = 0; n < N; n++)
                                for (int m = 0; m < N; m++)
                                  {
                                    if (dfac(j, n, m))
                                      {
                                        double fac = (i ? -1 : 1) * vfd(fb, e != f_e(fb, 0)) * dfac(j, n, m) / ve(e);

                                        // Update the right-hand side
                                        if (fd >= 0)
                                          secmem(fb, n) -= fac * mult * inco(fd, m);
                                        else
                                          {
                                            // Dirichlet boundary condition
                                            for (int d = 0; d < D; d++)
                                              secmem(fb, n) -= fac * nf(fb, d) / fs(fb) * ref_cast(Dirichlet, cls[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), N * d + m);
                                          }
                                        if (!incompressible_)
                                          secmem(fb, n) += fac * inco(fb, m);

                                        // Update the matrix
                                        if (mat)
                                          {
                                            if (fd >= 0)
                                              (*mat)(N * fb + n, N * fd + m) += fac * mult;

                                            if (!incompressible_)
                                              (*mat)(N * fb + n, N * fb + m) -= fac;
                                          }
                                      }
                                  }
                            }
                        }
                      // No equivalence: n_f * element operator
                      else
                        {
                          for (int j = 0; j < 2; j++)
                            {
                              const int eb = f_e(f, j);
                              for (int l = 0; l < e_f.dimension(1) ; l++)
                                {
                                  fc = e_f(eb, l);
                                  if (fc < 0) continue;

                                  double critere = fs(fc) * domaine.dot(&xv(fc, 0), &nf(fb, 0), &xp(eb, 0));
                                  if (std::abs(critere) > 1e-6 * ve(eb) * fs(fb))
                                    {
                                      for (int n = 0; n < N; n++)
                                        for (int m = 0; m < N; m++)
                                          if (dfac(j, n, m))
                                            {
                                              double fac = (i ? -1 : 1) * vfd(fb, e != f_e(fb, 0)) * dfac(j, n, m) * ((eb == f_e(fc, 0)) ? 1 : -1) * critere / (ve(e) * ve(eb) * fs(fb));
                                              secmem(fb, n) -= fac * inco(fc, m);
                                              if (mat && fac)
                                                (*mat)(N * fb + n, N * fc + m) += fac;
                                            }
                                    }
                                }
                              // Correction part if 'comp'
                              if (!incompressible_)
                                {
                                  for (int l = 0; l < e_f.dimension(1) ; l++)
                                    {
                                      fc = e_f(e, l);
                                      if (fc < 0) continue;

                                      double critere = fs(fc) * domaine.dot(&xv(fc, 0), &nf(fb, 0), &xp(e, 0));
                                      if (std::abs(critere) > 1e-6 * ve(e) * fs(fb))
                                        {
                                          for (int n = 0; n < N; n++)
                                            for (int m = 0; m < N; m++)
                                              {
                                                if (dfac(j, n, m))
                                                  {
                                                    double fac = (i ? -1 : 1) * vfd(fb, e != f_e(fb, 0)) * dfac(j, n, m) * ((e == f_e(fc, 0)) ? 1 : -1) * critere / (ve(e) * ve(e) * fs(fb));
                                                    secmem(fb, n) += fac * inco(fc, m);
                                                    if (mat && fac)
                                                      (*mat)(N * fb + n, N * fc + m) -= fac;
                                                  }
                                              }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
