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

#include <Fluide_Incompressible.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Milieu_composite_Euler.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Pb_Multiphase.h>
#include <Interprete.h>
#include <Domaine_VF.h>
#include<Fluide_reel_base.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>

Implemente_instanciable(Milieu_composite_Euler, "Milieu_composite_Euler", Milieu_composite);
// XD Milieu_composite_Euler Milieu_composite Milieu_composite_Euler -1 Composite medium made of several sub mediums.

Sortie& Milieu_composite_Euler::printOn(Sortie& os) const { return os; }
Entree& Milieu_composite_Euler::readOn(Entree& is) { return Milieu_composite::readOn(is); }

void Milieu_composite_Euler::discretiser(const Probleme_base& pb, const  Discretisation_base& dis)
{
  Milieu_composite::discretiser(pb, dis);
  res_en_T_ = true;
  inter_lu_->assoscier_pb(pb);
}

void Milieu_composite_Euler::init_energie_tot(DoubleTab& alpha_energie_tot_jdd) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler, equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& alpha = equation("alpha").inconnue().valeurs();
  const DoubleTab& p = qdm.pression().valeurs();
  const int Nb_phase = (int) fluides_.size();
  const int Nb_elem = qdm.domaine_dis().nb_elem_tot();

  for (int n = 0; n < Nb_phase; n++)
    {
      const DoubleTab& U = qdm.vitesse_phase(n).valeurs();
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base, get_fluid(n));
      for (int i = 0; i < Nb_elem; i++)
        {
          double nom_u2 = 0;
          for (int d = 0; d < dimension; d++)
            nom_u2 += U(i, d) * U(i, d);
          alpha_energie_tot_jdd(i, n) = alpha(i, n) * phase.init_energie_tot(rho(i, n), nom_u2, p(i, n));
        }
    }
}

void Milieu_composite_Euler::calculer_pression(DoubleTab& p) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler, equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& alpha = equation("alpha").inconnue().valeurs();
  const DoubleTab& alpha_rhoE = equation("alpha_energie_tot").inconnue().valeurs();
  DoubleTab rhoE = alpha_rhoE;
  tab_divide_any_shape(rhoE, alpha); // @suppress("Function cannot be resolved")
  const int Nb_phase = (int) fluides_.size();
  const int Nb_elem = qdm.domaine_dis().nb_elem_tot();

  for (int n = 0; n < Nb_phase; n++)
    {
      const DoubleTab& U = qdm.vitesse_phase(n).valeurs();
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base, get_fluid(n));
      for (int i = 0; i < Nb_elem; i++)
        {
          double nom_u2 = 0;
          for (int d = 0; d < dimension; d++)
            nom_u2 += U(i, d) * U(i, d);
          p(i, n) = phase.calculer_pression(rho(i, n), nom_u2, rhoE(i, n));
          assert(p(i, n) > 0);
        }
    }
}

void Milieu_composite_Euler::calculer_vitesse_son(DoubleTab& c) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler, equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& p = qdm.pression().valeurs();

  const int Nb_phase = (int) fluides_.size();
  const int Nb_elem = qdm.domaine_dis().nb_elem_tot();
  for (int n = 0; n < Nb_phase; n++)
    {
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base, get_fluid(n));
      for (int i = 0; i < Nb_elem; i++)
        c(i, n) = phase.calculer_vitesse_son(rho(i, n), p(i, n));
    }
}
