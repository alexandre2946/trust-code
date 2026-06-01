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

#include <Op_VEF_Face.h>
#include <Matrice_Morse.h>
#include <Equation_base.h>
#include <Sortie.h>
#include <Probleme_base.h>

#include <Champ_Uniforme.h>
#include <Schema_Temps_base.h>
#include <Milieu_base.h>
#include <Operateur_base.h>
#include <Operateur_Diff_base.h>
#include <Op_Conv_VEF_base.h>
#include <EcrFicPartage.h>
#include <SFichier.h>
#include <Process.h>
#include <Matrice_Morse_Diag.h>
#include <TRUSTTrav.h>
#include <Dirichlet_homogene.h>
#include <Periodique.h>
#include <Symetrie.h>

/*! @brief Dimensioning of the matrix that will receive the coefficients from convection and diffusion for the face-based case.
 *
 *  This matrix has a Morse matrix structure.
 *  We start by computing the sizes of arrays tab1 and tab2.
 *
 */

void Op_VEF_Face::dimensionner(const Domaine_VEF& le_dom, const Domaine_Cl_VEF& le_dom_cl, Matrice_Morse& la_matrice) const
{
  // Dimensioning of the matrix that will receive the coefficients from
  // convection and diffusion for the face-based case.
  // This matrix has a Morse matrix structure.
  // We start by computing the sizes of arrays tab1 and tab2.
  // To do so, we need to find the neighbouring faces of the current face.

  int nfin = le_dom.nb_faces_tot();
  int nb_faces_elem = le_dom.domaine().nb_faces_elem();
  const int nb_comp = le_dom_cl.equation().inconnue().valeurs().line_size();
  la_matrice.dimensionner(nfin * nb_comp, nfin * nb_comp, 0);

  auto& tab1 = la_matrice.get_set_tab1();
  auto& tab2 = la_matrice.get_set_tab2();
  auto& coeff = la_matrice.get_set_coeff();
  coeff = 0;

  const IntTab& elem_faces = le_dom.elem_faces();
  const IntTab& face_voisins = le_dom.face_voisins();

  // For each face, we associate an integer array and a list of reals:
  // voisines[i] = {j such that j>i and M(i,j) is non-zero }

  //  IntVect rang_voisin(nfin*nb_comp);
  IntTrav rang_voisin(nfin * nb_comp);
  rang_voisin = nb_comp;
  // Process all faces
  int j;
  ToDo_Kokkos("Port with kokkos ? It will be called once...");
  for (int num_face = 0; num_face < nfin; num_face++)
    {
      int elem1 = face_voisins(num_face, 0);
      int elem2 = face_voisins(num_face, 1);

      for (int i = 0; i < nb_faces_elem; i++)
        {
          if ((j = elem_faces(elem1, i)) != num_face)
            {
              for (int k = 0; k < nb_comp; k++)
                {
                  rang_voisin(num_face * nb_comp + k) += nb_comp;
                }
            }
          if (elem2 != -1)
            if ((j = elem_faces(elem2, i)) != num_face)
              {
                for (int k = 0; k < nb_comp; k++)
                  {
                    rang_voisin(num_face * nb_comp + k) += nb_comp;
                  }
              }
        }
    }

  // The neighbouring faces of num_face have now been counted;
  // dimension tab1 and tab2 to the number of faces

  tab1(0) = 1;
  for (int num_face = 0; num_face < nfin; num_face++)
    {
      for (int k = 0; k < nb_comp; k++)
        {
          tab1(num_face * nb_comp + 1 + k) = rang_voisin(num_face * nb_comp + k) + tab1(num_face * nb_comp + k);
        }
    }
  la_matrice.dimensionner(nfin * nb_comp, tab1(nfin * nb_comp) - 1);

  for (int num_face = 0; num_face < nfin; num_face++)
    {
      for (int k = 0; k < nb_comp; k++)
        {
          for (int kk = 0; kk < nb_comp; kk++)
            {
              int modulo = (k + kk) % nb_comp;
              tab2[tab1[num_face * nb_comp + k] - 1 + kk] = num_face * nb_comp + 1 + modulo;
            }
          rang_voisin[num_face * nb_comp + k] = (int)(tab1[num_face * nb_comp + k] + nb_comp - 1);
        }
    }

  // Process all faces
  for (int num_face = 0; num_face < nfin; num_face++)
    {
      int elem1 = face_voisins(num_face, 0);
      int elem2 = face_voisins(num_face, 1);

      for (int i = 0; i < nb_faces_elem; i++)
        {
          if ((j = elem_faces(elem1, i)) != num_face)
            {
              for (int k = 0; k < nb_comp; k++)
                {
                  for (int kk = 0; kk < nb_comp; kk++)
                    {
                      int modulo = (k + kk) % nb_comp;
                      tab2[rang_voisin[num_face * nb_comp + k] + kk] = j * nb_comp + 1 + modulo;
                    }
                  rang_voisin[num_face * nb_comp + k] += nb_comp;
                }
            }
          if (elem2 != -1)
            {
              if ((j = elem_faces(elem2, i)) != num_face)
                {
                  for (int k = 0; k < nb_comp; k++)
                    {
                      for (int kk = 0; kk < nb_comp; kk++)
                        {
                          int modulo = (k + kk) % nb_comp;
                          tab2[rang_voisin[num_face * nb_comp + k] + kk] = j * nb_comp + 1 + modulo;
                        }
                      rang_voisin[num_face * nb_comp + k] += nb_comp;
                    }
                }
            }
        }
    }
}

/*! @brief Modify the matrix coefficients and the right-hand side for Dirichlet boundary conditions.
 *
 */

void Op_VEF_Face::modifier_pour_Cl(const Domaine_VEF& le_dom, const Domaine_Cl_VEF& le_dom_cl, Matrice_Morse& la_matrice, DoubleTab& tab_secmem) const
{
  // Dimensioning of the matrix that will receive the coefficients from
  // convection and diffusion for the face-based case.
  // This matrix has a Morse matrix structure.
  // We start by computing the sizes of arrays tab1 and tab2.
  const Conds_lim& les_cl = le_dom_cl.les_conditions_limites();
  const DoubleTab& champ_inconnue = le_dom_cl.equation().inconnue().valeurs();
  const int nb_comp = champ_inconnue.line_size();
  ArrOfDouble normale(nb_comp);
  for (const auto &itr : les_cl)
    {
      const Cond_lim_base& la_cl = itr.valeur();
      const Front_VF& la_front_dis = ref_cast(Front_VF, la_cl.frontiere_dis());
      int nfaces = la_front_dis.nb_faces_tot();
      if (sub_type(Dirichlet, la_cl) || sub_type(Dirichlet_homogene, la_cl))
        {
          bool has_val_imp = sub_type(Dirichlet, la_cl);
          CDoubleTabView val_imp;
          if (has_val_imp) val_imp = ref_cast(Dirichlet, la_cl).tab_val_imp().view_ro();
          auto tab1 = la_matrice.get_tab1().view_ro();
          CIntArrView num_face = la_front_dis.num_face().view_ro();
          DoubleArrView coeff = la_matrice.get_set_coeff().view_wo();
          DoubleTabView secmem = tab_secmem.view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nfaces, KOKKOS_LAMBDA(const int ind_face)
          {
            int face = num_face(ind_face);
            for (int comp = 0; comp < nb_comp; comp++)
              {
                auto idiag = tab1[face * nb_comp + comp] - 1;
                coeff[idiag] = 1;
                // for the neighbours
                auto nbvois = tab1[face * nb_comp + 1 + comp] - tab1[face * nb_comp + comp];
                for (auto k = 1; k < nbvois; k++)
                  coeff[idiag + k] = 0;
                // for the right-hand side
                int j = nb_comp == 1 ? 0 : comp;
                secmem(face, j) = has_val_imp ? val_imp(ind_face, j) : 0.;
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else if (sub_type(Symetrie, la_cl) && le_dom_cl.equation().inconnue().nature_du_champ() == vectoriel)
        {
          const auto& tab1 = la_matrice.get_tab1();
          const auto& tab2 = la_matrice.get_tab2();
          const DoubleTab& face_normales = le_dom.face_normales();
          ArrOfDouble somme(la_matrice.nb_colonnes()); // Dimensioned to the maximum size
          ToDo_Kokkos("critical");
          for (int ind_face = 0; ind_face < nfaces; ind_face++)
            {
              int face = la_front_dis.num_face(ind_face);
              double max_coef = 0;
              int ind_max = -1;
              double n2 = 0;
              for (int comp = 0; comp < nb_comp; comp++)
                {
                  normale[comp] = face_normales(face, comp);
                  if (std::fabs(normale[comp]) > std::fabs(max_coef))
                    {
                      max_coef = normale[comp];
                      ind_max = comp;
                    }
                  n2 += normale[comp] * normale[comp];
                }
              normale /= sqrt(n2);
              max_coef = normale[ind_max];

              // Start by recomputing secmem=secmem-A*present to allow modifying A (and project at the same time)
              auto nb_coeff_ligne = tab1[face * nb_comp + 1] - tab1[face * nb_comp];
              for (auto k = 0; k < nb_coeff_ligne; k++)
                {
                  for (int comp = 0; comp < nb_comp; comp++)
                    {
                      int j = tab2[tab1[face * nb_comp + comp] - 1 + k] - 1;
                      //assert(j!=(face*nb_comp+comp));
                      //if ((j>=(face*nb_comp))&&(j<(face*nb_comp+nb_comp)))
                      const double coef_ij = la_matrice(face * nb_comp + comp, j);
                      int face2 = j / nb_comp;
                      int comp2 = j - face2 * nb_comp;
                      tab_secmem(face, comp) -= coef_ij * champ_inconnue(face2, comp2);
                    }
                }
              double somme_b = 0;

              for (int comp = 0; comp < nb_comp; comp++)
                somme_b += tab_secmem(face, comp) * normale[comp];

              // subtract secmem.n * n
              for (int comp = 0; comp < nb_comp; comp++)
                tab_secmem(face, comp) -= somme_b * normale[comp];

              // the same diagonal must be set everywhere, take the average
              double ref = 0;
              for (int comp = 0; comp < nb_comp; comp++)
                {
                  int j0 = face * nb_comp + comp;
                  ref += la_matrice(j0, j0);
                }
              ref /= nb_comp;

              for (int comp = 0; comp < nb_comp; comp++)
                {
                  int j0 = face * nb_comp + comp;
                  double rap = ref / la_matrice(j0, j0);
                  for (auto k = 0; k < nb_coeff_ligne; k++)
                    {
                      int j = tab2[tab1[j0] - 1 + k] - 1;
                      la_matrice(j0, j) *= rap;
                    }
                  assert(est_egal(la_matrice(j0, j0), ref));
                }
              // zero all off-diagonal coefficients of the block
              //
              for (auto k = 1; k < nb_coeff_ligne; k++)
                {
                  for (int comp = 0; comp < nb_comp; comp++)
                    {
                      int j = tab2[tab1[face * nb_comp + comp] - 1 + k] - 1;
                      assert(j != (face * nb_comp + comp));
                      if ((j >= (face * nb_comp)) && (j < (face * nb_comp + nb_comp)))
                        la_matrice(face * nb_comp + comp, j) = 0;
                    }
                }

              // for the off-diagonal blocks ensure that Aij.ni=0
              //ArrOfDouble somme(nb_coeff_ligne);
              for (auto k = 0; k < nb_coeff_ligne; k++)
                {
                  somme[(int)k] = 0;
                  int j = tab2[tab1[face * nb_comp] - 1 + k] - 1;

                  // the coefficient j must exist on all nb_comp rows
                  double dsomme = 0;
                  for (int comp = 0; comp < nb_comp; comp++)
                    dsomme += la_matrice(face * nb_comp + comp, j) * normale[comp];

                  // subtract sum*ni

                  for (int comp = 0; comp < nb_comp; comp++)
                    // only modify coefficients not involving u(face,comp)
                    if ((j < (face * nb_comp)) || (j >= (face * nb_comp + nb_comp)))
                      la_matrice(face * nb_comp + comp, j) -= (dsomme) * normale[comp];
                }
              // Finally recompute secmem=secmem+A*champ_inconnue (A has been heavily modified)
              for (auto k = 0; k < nb_coeff_ligne; k++)
                {
                  for (int comp = 0; comp < nb_comp; comp++)
                    {
                      int j = tab2[tab1[face * nb_comp + comp] - 1 + k] - 1;
                      int face2 = j / nb_comp;
                      int comp2 = j - face2 * nb_comp;
                      const double coef_ij = la_matrice(face * nb_comp + comp, j);
                      tab_secmem(face, comp) += coef_ij * champ_inconnue(face2, comp2);
                    }
                }
              {
                // verification
                double somme_c = 0;
                for (int comp = 0; comp < nb_comp; comp++)
                  somme_c += tab_secmem(face, comp) * normale[comp];
                // subtract secmem.n * n
                for (int comp = 0; comp < nb_comp; comp++)
                  tab_secmem(face, comp) -= somme_c * normale[comp];
              }
            }
        }
    }
}

void Op_VEF_Face::modifier_flux(const Operateur_base& op) const
{
  controle_modifier_flux_ = 1;
  DoubleTab& flux_bords_ = op.flux_bords();
  if (flux_bords_.nb_dim() != 2)
    return;
  const Probleme_base& pb = op.equation().probleme();

  const Domaine_VEF& le_dom_vef = ref_cast(Domaine_VEF, op.equation().domaine_dis());
  int nb_compo = flux_bords_.dimension(1);
  // Multiply the boundary flux by rho*Cp unless it is a diffusion operator with conductivity as the field
  if (op.equation().inconnue().le_nom() == "temperature" && !( sub_type(Operateur_Diff_base,op) && ref_cast(Operateur_Diff_base,op).diffusivite().le_nom() == "conductivite"))
    {
      const Champ_base& rho = op.equation().milieu().masse_volumique();
      const Champ_Don_base& Cp = op.equation().milieu().capacite_calorifique();
      bool rho_uniforme = sub_type(Champ_Uniforme,rho);
      bool cp_uniforme = sub_type(Champ_Uniforme,Cp);
      bool is_rho_u = (sub_type(Op_Conv_VEF_base, op) && ref_cast(Op_Conv_VEF_base,op).vitesse().le_nom() == "rho_u") ? true : false;
      const int nb_faces_bords = le_dom_vef.nb_faces_bord();
      CDoubleArrView rho_face = static_cast<const DoubleVect&>(rho.valeurs()).view_ro();
      CDoubleArrView Cp_face = static_cast<const DoubleVect&>(Cp.valeurs()).view_ro();
      DoubleArrView flux_bords = static_cast<DoubleVect&>(flux_bords_).view_rw();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_faces_bords, KOKKOS_LAMBDA(
                             const int face)
      {
        // if in QC temperature mode and div(rhou * T) has been computed,
        // do not multiply by rho again
        flux_bords(face) *= (is_rho_u ? 1 : rho_face(rho_uniforme?0:face)) * Cp_face(cp_uniforme?0:face);
      });
      end_gpu_timer(__KERNEL_NAME__);
    }

  // Multiply by rho if incompressible Navier-Stokes
  Nom nom_eqn = op.equation().que_suis_je();
  if (nom_eqn.debute_par("Navier_Stokes") && pb.milieu().que_suis_je() == "Fluide_Incompressible")
    {
      const Champ_base& rho = op.equation().milieu().masse_volumique();
      if (sub_type(Champ_Uniforme, rho))
        {
          double coef = rho.valeurs()(0, 0);
          int nb_faces_bord = le_dom_vef.nb_faces_bord();
          DoubleTabView flux_bords = flux_bords_.view_rw();

          Kokkos::parallel_for(
            start_gpu_timer(__KERNEL_NAME__),
            range_2D({0,0}, {nb_faces_bord,nb_compo}),
            KOKKOS_LAMBDA (int face, int k)
          {
            flux_bords(face, k) *= coef;
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
}

/*! @brief Print the face fluxes of a VEF operator (e.g. diffusion, convection).
 *
 */
int Op_VEF_Face::impr(Sortie& os, const Operateur_base& op) const
{
  const Domaine_VEF& le_dom_vef = ref_cast(Domaine_VEF, op.equation().domaine_dis());
  const DoubleTab& flux_bords_ = op.flux_bords();
  if (flux_bords_.nb_dim() != 2)
    {
      Cout << "Printing of fluxes is not implemented for the operator " << op.que_suis_je() << finl;
      return 1;
    }
  if (controle_modifier_flux_ == 0)
    if (max_abs_array(flux_bords_) != 0)
      {
        Cerr << op.que_suis_je() << " appelle  Op_VEF_Face::impr sans avoir appeler  Op_VEF_Face::modifier_flux, on arrete tout " << finl;
        Process::exit();
      }
  int nb_compo = flux_bords_.dimension(1);
  const Probleme_base& pb = op.equation().probleme();
  const Schema_Temps_base& sch = pb.schema_temps();
  // Print moments only if requested and if treating the velocity diffusion operator
  int impr_mom = 0;
  if (le_dom_vef.domaine().moments_a_imprimer() && sub_type(Operateur_Diff_base, op) && op.equation().inconnue().le_nom() == "vitesse")
    impr_mom = 1;

  const int impr_sum = (le_dom_vef.domaine().bords_a_imprimer_sum().est_vide() ? 0 : 1);
  const int impr_bord = (le_dom_vef.domaine().bords_a_imprimer().est_vide() ? 0 : 1);

  // Compute moments
  DoubleTab xgr;
  if (impr_mom)
    xgr = le_dom_vef.calculer_xgr();

  // Loop over the boundaries to sum the fluxes per boundary into the flux_bord array
  DoubleVect bilan(nb_compo);
  bilan = 0;
  int nb_cl = le_dom_vef.nb_front_Cl();
  DoubleTrav flux_bords(4, nb_cl, nb_compo);
  flux_bords = 0.;
  /*
   flux_bord(k)             ->     flux_bords(0,num_cl,k)
   flux_bord_perio1(k)      ->     flux_bords(1,num_cl,k)
   flux_bord_perio2(k)      ->     flux_bords(2,num_cl,k)
   moment(k)                ->     flux_bords(3,num_cl,k)
   */
  for (int num_cl = 0; num_cl < nb_cl; num_cl++)
    {
      const Cond_lim& la_cl = op.equation().domaine_Cl_dis().les_conditions_limites(num_cl);
      const Front_VF& frontiere_dis = ref_cast(Front_VF, la_cl->frontiere_dis());
      int ndeb = frontiere_dis.num_premiere_face();
      int nfin = ndeb + frontiere_dis.nb_faces();
      int perio = (sub_type(Periodique,la_cl.valeur()) ? 1 : 0);
      for (int face = ndeb; face < nfin; face++)
        {
          for (int k = 0; k < nb_compo; k++)
            {
              flux_bords(0, num_cl, k) += flux_bords_(face, k);
              if (perio)
                {
                  if (face < (ndeb + frontiere_dis.nb_faces() / 2))
                    flux_bords(1, num_cl, k) += flux_bords_(face, k);
                  else
                    flux_bords(2, num_cl, k) += flux_bords_(face, k);
                }
            }
          if (impr_mom)
            {
              // Compute the moment exerted by the fluid on the boundary (OM x F)
              if (Objet_U::dimension == 2)
                flux_bords(3, num_cl, 0) += flux_bords_(face, 1) * xgr(face, 0) - flux_bords_(face, 0) * xgr(face, 1);
              else
                {
                  flux_bords(3, num_cl, 0) += flux_bords_(face, 2) * xgr(face, 1) - flux_bords_(face, 1) * xgr(face, 2);
                  flux_bords(3, num_cl, 1) += flux_bords_(face, 0) * xgr(face, 2) - flux_bords_(face, 2) * xgr(face, 0);
                  flux_bords(3, num_cl, 2) += flux_bords_(face, 1) * xgr(face, 0) - flux_bords_(face, 0) * xgr(face, 1);
                }
            }
        } // end for face
    } // end for num_cl

  // Sum the contributions from each processor
  Process::mp_sum_for_each_item(flux_bords);

  // Write to files
  if (Process::je_suis_maitre())
    {
      op.ouvrir_fichier(Flux, "", 1);
      op.ouvrir_fichier(Flux_moment, "moment", impr_mom);
      op.ouvrir_fichier(Flux_sum, "sum", impr_sum);

      // Write time
      Flux.add_col(sch.temps_courant());
      if (impr_mom)
        Flux_moment.add_col(sch.temps_courant());
      if (impr_sum)
        Flux_sum.add_col(sch.temps_courant());

      // Write flux on boundaries
      for (int num_cl = 0; num_cl < nb_cl; num_cl++)
        {
          const Frontiere_dis_base& la_fr = op.equation().domaine_Cl_dis().les_conditions_limites(num_cl)->frontiere_dis();
          const Cond_lim& la_cl = op.equation().domaine_Cl_dis().les_conditions_limites(num_cl);
          int perio = (sub_type(Periodique,la_cl.valeur()) ? 1 : 0);
          for (int k = 0; k < nb_compo; k++)
            {
              if (perio)
                {
                  Flux.add_col(flux_bords(1, num_cl, k));
                  Flux.add_col(flux_bords(2, num_cl, k));
                }
              else
                Flux.add_col(flux_bords(0, num_cl, k));
              if (impr_mom)
                Flux_moment.add_col(flux_bords(3, num_cl, k));
              if (le_dom_vef.domaine().bords_a_imprimer_sum().contient(la_fr.le_nom()))
                Flux_sum.add_col(flux_bords(0, num_cl, k));

              // Sum the fluxes from all boundaries into the bilan array
              bilan(k) += flux_bords(0, num_cl, k);
            }
        }

      // Print the balance values and go to a new line
      for (int k = 0; k < nb_compo; k++)
        Flux.add_col(bilan(k));
      Flux << finl;
      if (impr_mom)
        Flux_moment << finl;
      if (impr_sum)
        Flux_sum << finl;
    }

  const LIST(Nom) &Liste_bords_a_imprimer = le_dom_vef.domaine().bords_a_imprimer();
  if (!Liste_bords_a_imprimer.est_vide())
    {
      EcrFicPartage Flux_face;
      op.ouvrir_fichier_partage(Flux_face, "", impr_bord);
      // Print per face if requested
      for (int num_cl = 0; num_cl < nb_cl; num_cl++)
        {
          const Frontiere_dis_base& la_fr = op.equation().domaine_Cl_dis().les_conditions_limites(num_cl)->frontiere_dis();
          const Cond_lim& la_cl = op.equation().domaine_Cl_dis().les_conditions_limites(num_cl);
          const Front_VF& frontiere_dis = ref_cast(Front_VF, la_cl->frontiere_dis());
          int ndeb = frontiere_dis.num_premiere_face();
          int nfin = ndeb + frontiere_dis.nb_faces();
          // Print per face
          if (Liste_bords_a_imprimer.contient(la_fr.le_nom()))
            {
              Flux_face << "# Flux par face sur " << la_fr.le_nom() << " au temps ";
              sch.imprimer_temps_courant(Flux_face);
              Flux_face << " : " << finl;
              const DoubleTab& xv = le_dom_vef.xv();
              for (int face = ndeb; face < nfin; face++)
                {
                  if (Objet_U::dimension == 2)
                    Flux_face << "# Face a x= " << xv(face, 0) << " y= " << xv(face, 1);
                  else if (Objet_U::dimension == 3)
                    Flux_face << "# Face a x= " << xv(face, 0) << " y= " << xv(face, 1) << " z= " << xv(face, 2);
                  for (int k = 0; k < nb_compo; k++)
                    Flux_face << " surface_face(m2)= " << le_dom_vef.face_surfaces(face) << " flux_par_surface(W/m2)= " << flux_bords_(face, k) / le_dom_vef.face_surfaces(face) << " flux(W)= "
                              << flux_bords_(face, k);
                  Flux_face << finl;
                }
              Flux_face.syncfile();
            }
        }
    }
  return 1;
}

/////////////////////////////////////////
// Method for the implicit scheme
/////////////////////////////////////////
void modif_matrice_pour_periodique_avant_contribuer(Matrice_Morse& matrice_morse, const Equation_base& eqn)
{
  const int nb_comp = eqn.inconnue().valeurs().line_size();
  const Domaine_Cl_dis_base& domaine_Cl_VEF = eqn.domaine_Cl_dis();
  const Domaine_VF& domaine_VEF = ref_cast(Domaine_VF, eqn.domaine_dis());
  int nb_bords = domaine_VEF.nb_front_Cl();
  int nb_faces_elem = domaine_VEF.elem_faces().dimension(1);
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique, la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          // only iterate over half the periodic faces;
          // the result will be copied to the associated face at the end
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces() / 2;
          CIntTabView elem_faces = domaine_VEF.elem_faces().view_ro();
          CIntTabView face_voisins = domaine_VEF.face_voisins().view_ro();
          CIntArrView face_associee = la_cl_perio.face_associee().view_ro();
          Matrice_Morse_View matrice;
          matrice.set(matrice_morse);
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(num1, num2),
                               KOKKOS_LAMBDA(const int num_face)
          {
            for (int dir = 0; dir < 2; dir++)
              {
                int elem1 = face_voisins(num_face, dir);
                int fac_asso = face_associee(num_face - num1) + num1;
                for (int i = 0; i < nb_faces_elem; i++)
                  {
                    int j = elem_faces(elem1, i);
                    for (int nc = 0; nc < nb_comp; nc++)
                      {
                        int n0 = num_face * nb_comp + nc;
                        int n0perio = fac_asso * nb_comp + nc;
                        if (((j == num_face) || (j == fac_asso)))
                          {
                            if (dir == 0)
                              {
                                assert(matrice(n0, n0perio) == 0);
                                assert(matrice(n0perio, n0) == 0);
                                assert(matrice(n0, n0) == matrice(n0perio, n0perio));
                                double coeff = (matrice(n0, n0)) / 2.;
                                matrice.store(n0, n0, coeff);
                                matrice.store(n0perio, n0perio, coeff);
                              }
                          }
                        else
                          {
                            for (int nc2 = 0; nc2 < nb_comp; nc2++)
                              {
                                int j20 = j * nb_comp + nc2;
                                assert(matrice(n0, j20) == matrice(n0perio, j20));
                                double coeff = (matrice(n0, j20) / 2.);
                                matrice.store(n0, j20, coeff);
                                matrice.store(n0perio, j20, coeff);
                              }
                          }
                      }
                  }
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      //  matrice.imprimer(Cerr);
    }
}

void modif_matrice_pour_periodique_apres_contribuer(Matrice_Morse& matrice_morse, const Equation_base& eqn)
{
  const int nb_comp = eqn.inconnue().valeurs().line_size();
  const Domaine_Cl_dis_base& domaine_Cl_VEF = eqn.domaine_Cl_dis();
  const Domaine_VF& domaine_VEF = ref_cast(Domaine_VF, eqn.domaine_dis());
  int nb_bords = domaine_VEF.nb_front_Cl();
  int nb_faces_elem = domaine_VEF.elem_faces().dimension(1);
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique, la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          // only iterate over half the periodic faces;
          // the result will be copied to the associated face at the end
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces() / 2;
          CIntTabView elem_faces = domaine_VEF.elem_faces().view_ro();
          CIntTabView face_voisins = domaine_VEF.face_voisins().view_ro();
          CIntArrView face_associee = la_cl_perio.face_associee().view_ro();
          Matrice_Morse_View matrice;
          matrice.set(matrice_morse);
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(num1, num2),
                               KOKKOS_LAMBDA(const int num_face)
          {
            for (int dir = 0; dir < 2; dir++)
              {
                int elem1 = face_voisins(num_face, dir);
                int fac_asso = face_associee(num_face - num1) + num1;
                for (int i = 0; i < nb_faces_elem; i++)
                  {
                    int j = elem_faces(elem1, i);
                    //if (j!=num_face)
                    for (int nc = 0; nc < nb_comp; nc++)
                      {
                        int n0 = num_face * nb_comp + nc;
                        //    int j0=j*nb_comp+nc;
                        int n0perio = fac_asso * nb_comp + nc;
                        if (((j == num_face) || (j == fac_asso)))
                          {
                            if (dir == 0)
                              {
                                for (int nc2 = 0; nc2 < nb_comp; nc2++)
                                  {
                                    int j0 = num_face * nb_comp + nc2;
                                    int j0perio = fac_asso * nb_comp + nc2;
                                    matrice.atomic_add(n0, j0, matrice(n0, j0perio));
                                    matrice.store(n0, j0perio, 0);
                                    matrice.atomic_add(n0perio, j0perio, matrice(n0perio, j0));
                                    matrice.store(n0perio, j0, 0);
                                    double coeff = (matrice(n0, j0) +
                                                    matrice(n0perio, j0perio));
                                    matrice.store(n0, j0, coeff);
                                    if (nc != nc2)
                                      {
                                        matrice.store(n0perio, j0, coeff);
                                        matrice.store(n0perio, j0perio, 0);
                                      }
                                    else
                                      matrice.store(n0perio, j0perio, coeff);
                                  }
                              }
                          }
                        else
                          {
                            for (int nc2 = 0; nc2 < nb_comp; nc2++)
                              {
                                int j20 = j * nb_comp + nc2;
                                double coeff = (matrice(n0, j20) + matrice(n0perio, j20));
                                matrice.store(n0, j20, coeff);
                                matrice.store(n0perio, j20, coeff);
                              }
                          }
                      }
                  }
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      //  matrice.imprimer(Cerr);
    }
}

/*! @brief Divides the coefficients on the periodic face rows by 2 in preparation for applying modifier_matrice_pour_periodique_apres_contribuer, which will sum the 2 rows of the associated periodic faces.
 *
 */

void Op_VEF_Face::modifier_matrice_pour_periodique_avant_contribuer(Matrice_Morse& matrice, const Equation_base& eqn) const
{
  // if matrice_morse_diag, no contribution n0 n0perio
  if (sub_type(Matrice_Morse_Diag, matrice))
    return;
  modif_matrice_pour_periodique_avant_contribuer(matrice, eqn);
}
/*! @brief Sums the 2 rows of the associated periodic faces, allowing computations in the code without needing to find the associated face.
 *
 *  Only half the periodic faces are iterated over in contribuer_a_avec (in general).
 *
 */

void Op_VEF_Face::modifier_matrice_pour_periodique_apres_contribuer(Matrice_Morse& matrice, const Equation_base& eqn) const
{
  // if matrice_morse_diag, no contribution n0 n0perio
  if (sub_type(Matrice_Morse_Diag, matrice))
    return;

  modif_matrice_pour_periodique_apres_contribuer(matrice, eqn);

  // verify that the matrix is properly periodic
#ifndef  NDEBUG

  const int nb_comp = eqn.inconnue().valeurs().line_size();

  const Domaine_Cl_dis_base& domaine_Cl_VEF = eqn.domaine_Cl_dis();
  const Domaine_VF& domaine_VEF = ref_cast(Domaine_VF, eqn.domaine_dis());
  int nb_bords = domaine_VEF.nb_front_Cl();

  const auto& tab1 = matrice.get_tab1();
  const auto& tab2 = matrice.get_tab2();

  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      int num1 = le_bord.num_premiere_face();
      //      int num2 = num1 + le_bord.nb_faces();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          int fac_asso;
          // only iterate over half the periodic faces;
          // the result will be copied to the associated face at the end
          int num2 = num1 + le_bord.nb_faces() / 2;
          for (int num_face = num1; num_face < num2; num_face++)
            for (int nc = 0; nc < nb_comp; nc++)
              {
                fac_asso = la_cl_perio.face_associee(num_face - num1) + num1;
                int n0 = num_face * nb_comp + nc;
                int n0perio = fac_asso * nb_comp + nc;
                // verify that the 2 rows are identical (except for the diagonal entry which is not in the same position)
                for (auto j = tab1[n0] - 1; j < tab1[n0 + 1] - 1; j++)
                  {
                    int c = tab2[j] - 1;
                    if ((c != n0) && (c != n0perio))
                      {
                        double test = matrice(n0, c) - matrice(n0perio, c);
                        if (test != 0)
                          {
                            Cerr << "Pb non-periodic matrix face" << num_face << " component " << nc << " column " << c << finl;
                            Cerr << " diff " << test << " coef1 " << matrice(n0, c) << " coef2 " << matrice(n0perio, c) << finl;
                            Process::exit();
                          }
                      }
                  }
                if ((matrice(n0, n0perio) != 0) || (matrice(n0perio, n0) != 0))
                  {
                    Cerr << "Pb non-periodic matrix face" << num_face << " component " << nc << finl;
                    Cerr << " non-zero coef" << matrice(n0, n0perio) << " " << matrice(n0perio, n0) << finl;
                    Process::exit();
                  }
                if (matrice(n0, n0) != matrice(n0perio, n0perio))
                  {
                    double test = matrice(n0, n0perio) - matrice(n0perio, n0);
                    Cerr << "Pb non-periodic matrix face" << num_face << " component " << nc << finl;
                    Cerr << " diff " << test << " differing coef " << matrice(n0, n0) << " " << matrice(n0perio, n0perio) << finl;
                    Process::exit();
                  }
              }
        }
    }

  /*
   Matrice_Morse es(matrice);
   modif_matrice_pour_periodique_avant_contribuer(es,eqn);
   modif_matrice_pour_periodique_apres_contribuer(es,eqn);
   es.coeff_-=(matrice.coeff_);
   Cerr<<" error after modifier_matrice_pour_periodique_apres_contribuer"<< mp_max_abs_vect(es.coeff_)<<finl;
   assert(mp_max_abs_vect(es.coeff_)<1e-9);
   */

#endif
}
