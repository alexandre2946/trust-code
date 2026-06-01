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

#include <Op_NConserv_HLL_Coloc_Vect.h>
#include <Interface_Baer_Nunziato.h>
#include <Conservation_Euler_base.h>
#include <Milieu_composite_Euler.h>
#include <Coloc_Operator_tools.h>
#include <Sortie_supersonique.h>
#include <Champ_Inc_P0_base.h>
#include <Fluide_reel_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Momentum_Euler.h>
#include <Domaine_Coloc.h>
#include <Dirichlet.h>
#include <Pb_Euler.h>
#include <array>

Implemente_instanciable(Op_NConserv_HLL_Coloc_Vect, "Op_NConserv_HLL_Coloc_Vect", Op_NConserv_Coloc_base);

Sortie& Op_NConserv_HLL_Coloc_Vect::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_HLL_Coloc_Vect::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is); return is;}

void Op_NConserv_HLL_Coloc_Vect::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();
  assert(secmem.line_size() / Objet_U::dimension == 2);

  const int N = domaine.nb_faces();
  DoubleTrav num_flux_left(N, Objet_U::dimension);
  DoubleTrav num_flux_right(N, Objet_U::dimension);

  Abgral_scheme(num_flux_left, num_flux_right);

  for (int f = 0; f < N; f++)
    for (int i = 0; i < 2; i++)
      {
        const int e = f_e(f, i);
        if (e >= 0 && e < domaine.nb_elem())
          {
            double val = (i ? num_flux_right(f, 0) : num_flux_left(f, 0)) * fs(f);
            secmem(e, 0) -= val;
            secmem(e, 1) += val;
            val = (i ? num_flux_right(f, 1) : num_flux_left(f, 1)) * fs(f);
            secmem(e, 2) -= val;
            secmem(e, 3) += val;
          }
      }
}

void Op_NConserv_HLL_Coloc_Vect::Abgral_scheme(DoubleTab& num_flux_left, DoubleTab& num_flux_right) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const IntTab& f_e = domaine.face_voisins();
  const IntTab& fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();

  const Momentum_Euler& eq = ref_cast(Momentum_Euler, equation());
  const DoubleTab& vit_n = eq.vitesse_normale();
  const DoubleTab& p = eq.pression().valeurs();

  const Pb_Euler& pb = ref_cast(Pb_Euler, equation().probleme());
  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  const DoubleTab& c = pb.equation_qdm().vitesse_son();

  const Interface_Baer_Nunziato& interface = ref_cast(Interface_Baer_Nunziato, ref_cast(Milieu_composite_Euler,pb.milieu()).interface_phase());
  const int n = interface.id_phase_vitesse_inter();
  const int m = interface.id_phase_pression_inter();
  const int nb_phase = pb.nb_phases();

  flux_bords_.resize(domaine.nb_faces_bord(), num_flux_left.line_size());
  flux_bords_ = 0.;

  for (int f = 0; f < domaine.nb_faces(); f++)
    {
      const int el = f_e(f, 0), er = f_e(f, 1);
      if (fcl(f, 0) == 0) // faces internes
        {
          double Sm = 0., Sp = 0., un_l = 0.;
          compute_non_conservative_hll_left_bounds(vit_n, c, f, el, er, m, n, nb_phase, Sm, Sp, un_l);
          for (int d = 0; d < Objet_U::dimension; d++)
            {
              double n_d = domaine.face_normales(f, d) / domaine.face_surfaces(f);
              num_flux_left(f, d) = (Sp * alpha(el, 0) - Sm * alpha(er, 0)) * p(el, m) * n_d;
              num_flux_left(f, d) /= -(Sp - Sm);
            }

          compute_non_conservative_hll_right_bounds(vit_n, c, f, el, er, m, n, nb_phase, Sm, Sp, un_l);

          for (int d = 0; d < Objet_U::dimension; d++)
            {
              double n_d = -domaine.face_normales(f, d) / domaine.face_surfaces(f);
              num_flux_right(f, d) = (Sp * alpha(er, 0) - Sm * alpha(el, 0)) * p(er, m) * n_d;
              num_flux_right(f, d) /= -(Sp - Sm);
            }
        }
      else // faces bords
        {
          //lookup arrays for boundary conditions: fcl(f, .) = { BC type, BC index, face index within BC }
          //BC types: 0 -> no BC
          //          1 -> Neumann
          //          2 -> Navier or symmetry
          //          3 -> Dirichlet or Neumann_homogene
          //          4 -> Dirichlet_homogene
          //          5 -> Periodique

          assert(er < 0 && el >= 0 && vit_n(f, 0) != -123.123);
          const int e = el;

          const Conds_lim& cls_qdm = pb.equation_qdm().domaine_Cl_dis().les_conditions_limites();
          const Conds_lim& cls_alpha = pb.equation_fraction().domaine_Cl_dis().les_conditions_limites();

          std::array<double, 3> normal { 0., 0., 0. };
          for (int d = 0; d < Objet_U::dimension; d++)
            normal[d] = domaine.face_normales(f, d) / domaine.face_surfaces(f);

          if (sub_type(Sortie_supersonique, cls_qdm[fcl(f, 1)].valeur()))
            {
              const double alpha_bord = alpha(e, 0);
              for (int d = 0; d < Objet_U::dimension; d++)
                {
                  num_flux_left(f, d) = -alpha_bord * p(e, m) * normal[d];
                  flux_bords_(f, d) = num_flux_left(f, d) * domaine.face_surfaces(f);
                }
            }
          else if (sub_type(Dirichlet, cls_qdm[fcl(f, 1)].valeur()))
            {
              const double alpha_bord = ref_cast(Dirichlet, cls_alpha[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), 0);
              for (int d = 0; d < Objet_U::dimension; d++)
                {
                  num_flux_left(f, d) = -alpha_bord * p(e, m) * normal[d];
                  flux_bords_(f, d) = num_flux_left(f, d) * domaine.face_surfaces(f);
                }
            }
          else if ( sub_type(Symetrie,cls_qdm[fcl(f, 1)].valeur()) && !sub_type(Sortie_supersonique, cls_qdm[fcl(f, 1)].valeur()))
            {
              //Slip wall : u_n=-u_n
              const double alpha_bord = alpha(e, 0);
              for (int d = 0; d < Objet_U::dimension; d++)
                {
                  num_flux_left(f, d) = -alpha_bord * p(e, m) * normal[d];
                  flux_bords_(f, d) = num_flux_left(f, d) * domaine.face_surfaces(f);
                }
            }
          else
            {
              Cerr << " The BC of type " << fcl(f, 0) << " for the equation " << eq.que_suis_je() << " is not available .....\n";
              Process::exit();
            }
        }
    }
}
