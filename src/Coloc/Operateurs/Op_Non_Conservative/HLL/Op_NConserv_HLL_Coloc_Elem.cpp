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

#include <Op_NConserv_HLL_Coloc_Elem.h>
#include <Interface_Baer_Nunziato.h>
#include <Milieu_composite_Euler.h>
#include <Neumann_paroi_flux_nul.h>
#include <Coloc_Operator_tools.h>
#include <Entree_supersonique.h>
#include <Sortie_supersonique.h>
#include <Conservation_Euler.h>
#include <Champ_Inc_P0_base.h>
#include <Fluide_reel_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Momentum_Euler.h>
#include <Domaine_Coloc.h>
#include <Pb_Euler.h>

Implemente_instanciable(Op_NConserv_HLL_Coloc_Elem, "Op_NConserv_HLL_Coloc_Elem", Op_NConserv_Coloc_Elem_base);

Sortie& Op_NConserv_HLL_Coloc_Elem::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_HLL_Coloc_Elem::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is); return is;}

void Op_NConserv_HLL_Coloc_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();
  assert(secmem.line_size() == 2);
  const int nb_faces = domaine.nb_faces();
  DoubleTrav num_flux_left(nb_faces), num_flux_right(nb_faces);

  Abgral_scheme(num_flux_left, num_flux_right);

  for (int f = 0; f < nb_faces; f++)
    for (int i = 0; i < 2; i++)
      {
        int e = f_e(f, i);
        if (e >= 0 && e < domaine.nb_elem())
          {
            secmem(e, 0) -= (i ? num_flux_right(f) : num_flux_left(f)) * fs(f);
            secmem(e, 1) += (i ? num_flux_right(f) : num_flux_left(f)) * fs(f);
          }
      }
}

void Op_NConserv_HLL_Coloc_Elem::Abgral_scheme(DoubleTab& num_flux_left, DoubleTab& num_flux_right) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const IntTab& f_e = domaine.face_voisins();
  const IntTab& fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();

  const Pb_Euler& pb = ref_cast(Pb_Euler, equation().probleme());
  const Conservation_Euler& eq = ref_cast(Conservation_Euler, equation());

  const Interface_Baer_Nunziato& interface = ref_cast(Interface_Baer_Nunziato, ref_cast(Milieu_composite_Euler,pb.milieu()).interface_phase());
  const int n = interface.id_phase_vitesse_inter();
  const int m = interface.id_phase_pression_inter();
  const int nb_phases = pb.nb_phases();

  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  const DoubleTab& vit_n = pb.equation_qdm().vitesse_normale();
  const DoubleTab& p = pb.equation_qdm().pression().valeurs();
  const DoubleTab& c = pb.equation_qdm().vitesse_son();

  const Conds_lim& cls = equation().domaine_Cl_dis().les_conditions_limites();
  const Conds_lim& cls_alpha = pb.equation_fraction().domaine_Cl_dis().les_conditions_limites();

  // faces internes
  if (sub_type(Fraction_Euler, equation()))
    {
      for (int f = 0; f < domaine.nb_faces(); f++)
        if (fcl(f, 0) == 0)
          {
            const int el = f_e(f, 0), er = f_e(f, 1);
            double Sm = 0., Sp = 0., un_l = 0.;
            compute_non_conservative_hll_left_bounds(vit_n, c, f, el, er, m, n, nb_phases, Sm, Sp, un_l);

            num_flux_left(f) = (Sp * alpha(el, 0) - Sm * alpha(er, 0)) * un_l + Sp * Sm * (alpha(er, 0) - alpha(el, 0));
            num_flux_left(f) /= (Sp - Sm);

            compute_non_conservative_hll_right_bounds(vit_n, c, f, el, er, m, n, nb_phases, Sm, Sp, un_l);

            num_flux_right(f) = (Sp * alpha(er, 0) - Sm * alpha(el, 0)) * un_l + Sp * Sm * (alpha(el, 0) - alpha(er, 0));
            num_flux_right(f) /= (Sp - Sm);
          }
    }
  else if (sub_type(Energy_Euler, equation()))
    {
      for (int f = 0; f < domaine.nb_faces(); f++)
        if (fcl(f, 0) == 0)
          {
            const int el = f_e(f, 0), er = f_e(f, 1);
            double Sm = 0., Sp = 0., un_l = 0.;
            compute_non_conservative_hll_left_bounds(vit_n, c, f, el, er, m, n, nb_phases, Sm, Sp, un_l);

            num_flux_left(f) = (Sp * alpha(el, 0) - Sm * alpha(er, 0)) * un_l * p(el, m);
            num_flux_left(f) /= -(Sp - Sm);

            compute_non_conservative_hll_right_bounds(vit_n, c, f, el, er, m, n, nb_phases, Sm, Sp, un_l);

            num_flux_right(f) = (Sp * alpha(er, 0) - Sm * alpha(el, 0)) * un_l * p(er, m);
            num_flux_right(f) /= -(Sp - Sm);
          }
    }
  else
    {
      Cerr << "Op_NConserv_HLL_Coloc_Elem should not be used for equation " << equation().que_suis_je() << finl;
      Process::exit();
    }

  // faces bords
  for (int f = 0; f < domaine.nb_faces(); f++)
    {
      if (fcl(f, 0) != 0)
        {
          const int e = f_e(f, 0) >= 0 ? f_e(f, 0) : f_e(f, 1); //pas besoin
          assert(f_e(f, 0) >= 0 && vit_n(f, 0) != 123.123); //pas besoin
          num_flux_right(f) = 123.123; //pas besoin

          if (sub_type(Sortie_supersonique, cls[fcl(f, 1)].valeur()))
            {
              num_flux_left(f) = eq.termes_NonConservatif(alpha(e, 0), vit_n(f, n), p(e, m));
            }
          else if (sub_type(Entree_supersonique, cls[fcl(f, 1)].valeur())) // Dirichlet  : 6
            {
              const double alpha_bord = ref_cast(Dirichlet, cls_alpha[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), 0);
              num_flux_left(f) = eq.termes_NonConservatif(alpha_bord, vit_n(f, n), p(e, m));
            }
          else if (sub_type(Neumann_paroi_flux_nul, cls[fcl(f, 1)].valeur())) //Neumann_homogene : 5
            {
              num_flux_left(f) = eq.termes_NonConservatif(alpha(e, 0), vit_n(f, n), p(e, m));
            }
          else
            {
              Cerr << "The BC of type " << fcl(f, 0) << "for the equation " << eq.que_suis_je() << " is not available \n";
              Process::exit();
            }
        }
    }
}
