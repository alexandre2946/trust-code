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

#include <Op_Conv_Coloc_Vect_base.h>
#include <Milieu_composite_Euler.h>
#include <Sortie_supersonique.h>
#include <Entree_supersonique.h>
#include <Conservation_Euler.h>
#include <Champ_Inc_P0_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Fluide_reel_base.h>
#include <Momentum_Euler.h>
#include <Domaine_Coloc.h>
#include <Pb_Euler.h>
#include <array>

Implemente_base(Op_Conv_Coloc_Vect_base,"Op_Conv_Coloc_Vect_base",Op_Conv_Coloc_base);

Sortie& Op_Conv_Coloc_Vect_base::printOn(Sortie& os) const { return Operateur_Conv_base::printOn(os); }
Entree& Op_Conv_Coloc_Vect_base::readOn(Entree& is) { Operateur_Conv_base::readOn(is);  return is; }

void Op_Conv_Coloc_Vect_base::Riemann_solver(DoubleTab& num_flux) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const IntTab& f_e = domaine.face_voisins();
  const IntTab& fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();
  const Conds_lim& cls = le_dcl_coloc_->les_conditions_limites();
  const Momentum_Euler& eq = ref_cast(Momentum_Euler, equation());
  const DoubleTab& vit_n = eq.vitesse_normale();
  const DoubleTab& vit = eq.vitesse().valeurs();
  const DoubleTab& p = eq.pression().valeurs();
  const Pb_Euler& pb = ref_cast(Pb_Euler, equation().probleme());
  const int nb_phase = pb.nb_phases();
  const DoubleTab& rho = pb.equation_masse().densite().valeurs();
  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();

  // compute left/right fluxes on internal faces
  DoubleTrav flux_l(domaine.nb_faces(), num_flux.line_size()), flux_r(domaine.nb_faces(), num_flux.line_size());
  eq.compute_fluxes_on_all_faces(flux_l, flux_r);
  scheme(num_flux, flux_l, flux_r);

  // Boundary faces treatement
  for (int f = 0; f < domaine.nb_faces(); f++)
    {
      if (fcl(f, 0) != 0)
        {
          //tableaux de correspondance lies aux CLs : fcl(f, .) = { type de CL, num de la CL, indice de la face dans la CL }
          //types de CL : 0 -> pas de CL
          //              1 -> Neumann
          //              2 -> Navier ou symetrie
          //              3 -> Dirichlet ou Neumann_homogene
          //              4 -> Dirichlet_homogene
          //              5 -> Periodique

          const int e = f_e(f, 0) >= 0 ? f_e(f, 0) : f_e(f, 1); //pas besoin
          assert(f_e(f, 0) >= 0 && vit_n(f, 0) != 123.123); //pas besoin

          std::array<double, 3> normal { 0., 0., 0. };
          for (int d = 0; d < Objet_U::dimension; d++)
            normal[d] = domaine.face_normales(f, d) / domaine.face_surfaces(f);

          if (sub_type(Sortie_supersonique, cls[fcl(f, 1)].valeur())) // 1 -> Neumann
            {
              for (int n = 0; n < nb_phase; n++)
                {
                  const double p_bord = p(e, n);
                  const double rho_bord = rho(e, n);
                  const double alpha_bord = alpha(e, n);

                  for (int d = 0; d < Objet_U::dimension; d++)
                    {
                      num_flux(f, n + nb_phase * d) = normal[d] * p_bord + vit(e, n + nb_phase * d) * rho_bord * vit_n(f, n);
                      num_flux(f, n + nb_phase * d) *= alpha_bord;
                    }
                }
            }
          else if (sub_type(Entree_fluide_vitesse_imposee, cls[fcl(f, 1)].valeur())) // 3 -> Dirichlet
            {
              const Conds_lim& cls_alpha = pb.equation_fraction().domaine_Cl_dis().les_conditions_limites();
              const Conds_lim& cls_rho = pb.equation_masse().domaine_Cl_dis().les_conditions_limites();
              const Conds_lim& cls_p = pb.equation_energie().domaine_Cl_dis().les_conditions_limites();

              for (int n = 0; n < nb_phase; n++)
                {
                  std::array<double, 3> vitesse_bord { 0., 0., 0. };
                  double vitesse_normale_bord = 0;

                  const double rho_bord = ref_cast(Dirichlet, cls_rho[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n);
                  const double p_bord = ref_cast(Dirichlet, cls_p[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n);
                  const double alpha_bord = ref_cast(Dirichlet, cls_alpha[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n);

                  for (int d = 0; d < Objet_U::dimension; d++)
                    {
                      vitesse_bord[d] = ref_cast(Dirichlet, cls[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), n + nb_phase * d);
                      vitesse_normale_bord += vitesse_bord[d] * normal[d];
                    }

                  for (int d = 0; d < Objet_U::dimension; d++)
                    {
                      num_flux(f, n + nb_phase * d) = normal[d] * p_bord + vitesse_bord[d] * rho_bord * vitesse_normale_bord;
                      num_flux(f, n + nb_phase * d) *= alpha_bord;
                    }
                }
            }
          else if ( sub_type(Symetrie,cls[fcl(f, 1)].valeur()) && !sub_type(Sortie_supersonique, cls[fcl(f, 1)].valeur()))
            {
              //Slip wall : u_n=0
              for (int n = 0; n < nb_phase; n++)
                {
                  std::array<double, 3> vitesse_bord { 0., 0., 0. };

                  const double p_bord = p(e, n);
                  const double rho_bord = rho(e, n);
                  const double alpha_bord = alpha(e, n);
                  double vitesse_normale_bord = 0;

                  for (int d = 0; d < Objet_U::dimension; d++)
                    {
                      num_flux(f, n + nb_phase * d) = normal[d] * p_bord + vitesse_bord[d] * rho_bord * vitesse_normale_bord;
                      num_flux(f, n + nb_phase * d) *= alpha_bord;
                    }
                }
            }
          else
            {
              Cerr << " La CL de type " << fcl(f, 0) << " pour l'equation " << eq.que_suis_je() << " n est pas diponible .....\n";
              Process::exit();
            }

        }
    }
}
