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

#include <Op_Conv_Coloc_Elem_base.h>
#include <Neumann_paroi_flux_nul.h>
#include <Milieu_composite_Euler.h>
#include <Sortie_supersonique.h>
#include <Entree_supersonique.h>
#include <Conservation_Euler.h>
#include <Champ_Inc_P0_base.h>
#include <Fluide_reel_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Momentum_Euler.h>
#include <Domaine_Coloc.h>
#include <Pb_Euler.h>
#include <array>

Implemente_base(Op_Conv_Coloc_Elem_base,"Op_Conv_Coloc_Elem_base",Op_Conv_Coloc_base);

Sortie& Op_Conv_Coloc_Elem_base::printOn(Sortie& os) const { return Op_Conv_Coloc_base::printOn(os); }
Entree& Op_Conv_Coloc_Elem_base::readOn(Entree& is) { Op_Conv_Coloc_base::readOn(is);  return is; }

void Op_Conv_Coloc_Elem_base::Riemann_solver(DoubleTab& num_flux) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const Pb_Euler& pb = ref_cast(Pb_Euler, equation().probleme());
  const DoubleTab& w = le_champ_inco->valeurs();
  const Conds_lim& cls = le_dcl_coloc_->les_conditions_limites();
  const IntTab& f_e = domaine.face_voisins();
  const IntTab& fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();
  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  const Conservation_Euler& eq = ref_cast(Conservation_Euler, equation());
  const DoubleTab& vit_n = pb.equation_qdm().vitesse_normale();
  const DoubleTab& p = pb.equation_qdm().pression().valeurs();
  const int nb_phases = pb.nb_phases();

  // compute left/right fluxes
  DoubleTrav flux_l(domaine.nb_faces(), num_flux.line_size()), flux_r(domaine.nb_faces(), num_flux.line_size());
  eq.compute_fluxes_on_all_faces(flux_l, flux_r);

  for (int f = 0; f < domaine.nb_faces(); f++)
    {
      if (fcl(f, 0) == 0) // face interne
        scheme(num_flux, f, flux_l, flux_r);
      else
        {
          const int e = f_e(f, 0) >= 0 ? f_e(f, 0) : f_e(f, 1); //pas besoin
          // assert(f_e(f,0) >= 0 && vit_n (f,0)!=123.123); //pas besoin

          //tableaux utilitaires sur les CLs : fcl(f, .) = (type de la CL, no de la CL, indice dans la CL)
          //types de CL : 0 -> pas de CL
          //              1 -> Echange_externe_impose
          //              2 -> Echange_global_impose
          //              3 -> Echange_contact_Coloc
          //              4 -> Neumann_paroi
          //              5 -> Neumann_val_ext ou Neumann_homogene ou Symetrie
          //              6 -> Dirichlet
          //              7 -> Dirichlet_homogene

          std::array<double, 3> normal { 0., 0., 0. };
          for (int d = 0; d < Objet_U::dimension; d++)
            normal[d] = domaine.face_normales(f, d) / domaine.face_surfaces(f);

          if (sub_type(Sortie_supersonique, cls[fcl(f, 1)].valeur())) //Neumann_val_ext : 5
            {
              for (int n = 0; n < nb_phases; n++)
                num_flux(f, n) = eq.flux_bord(w(e, n), vit_n(f, n), alpha(e, n) * p(e, n));
            }
          else if (sub_type(Neumann_paroi_flux_nul, cls[fcl(f, 1)].valeur())) //Neumann_homogene : 5
            {
              for (int n = 0; n < nb_phases; n++)
                num_flux(f, n) = 0;
            }
          else if (sub_type(Entree_supersonique, cls[fcl(f, 1)].valeur())) // Dirichlet  : 6
            {
              const Conds_lim& cls_alpha = pb.equation_fraction().domaine_Cl_dis().les_conditions_limites();
              const Conds_lim& cls_rho = pb.equation_masse().domaine_Cl_dis().les_conditions_limites();
              const Conds_lim& cls_qdm = pb.equation_qdm().domaine_Cl_dis().les_conditions_limites();
              const Conds_lim& cls_p = pb.equation_energie().domaine_Cl_dis().les_conditions_limites();

              for (int n = 0; n < nb_phases; n++)
                {
                  std::array<double, 3> vitesse_bord { 0., 0., 0. };
                  double vitesse_normale_bord = 0;
                  double norme_vitesse = 0;
                  const double alpha_bord = ref_cast(Dirichlet, cls_alpha[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n);
                  const double p_bord = ref_cast(Dirichlet, cls_p[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n);
                  const double rho_bord = ref_cast(Dirichlet, cls_rho[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n);

                  for (int d = 0; d < Objet_U::dimension; d++)
                    {
                      vitesse_bord[d] = ref_cast(Dirichlet, cls_qdm[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n + nb_phases * d);
                      vitesse_normale_bord += vitesse_bord[d] * normal[d];
                      norme_vitesse += vitesse_bord[d] * vitesse_bord[d];
                    }
                  /* TODO a refaire */
                  const Fluide_reel_base& phase = ref_cast(Fluide_reel_base, ref_cast(Milieu_composite_Euler,eq.milieu()).get_fluid(n));
                  const double inco_bord = (!(sub_type(Energy_Euler, equation()))) ? alpha_bord * rho_bord : alpha_bord * phase.init_energie_tot(rho_bord, norme_vitesse, p_bord);
                  num_flux(f, n) = eq.flux_bord(inco_bord, vitesse_normale_bord, alpha_bord * p_bord);
                }
            }
          else
            {
              Cerr << "The BC of type " << fcl(f, 0) << "for the equation " << eq.que_suis_je() << " is not available \n";
              Process::exit();
            }
        }
    }
}
