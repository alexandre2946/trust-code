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

#include <Interface_Baer_Nunziato.h>
#include <Milieu_composite_Euler.h>
#include <Discretisation_base.h>
#include <Fluide_reel_base.h>
#include <Champ_Uniforme.h>
#include <Momentum_Euler.h>
#include <Pb_Euler.h>

Implemente_instanciable(Milieu_composite_Euler, "Milieu_composite_Euler", Milieu_composite);
// XD Milieu_composite_Euler Milieu_composite Milieu_composite_Euler INHERITS_BRACE Composite medium made of several sub mediums.

Sortie& Milieu_composite_Euler::printOn(Sortie& os) const { return os; }
Entree& Milieu_composite_Euler::readOn(Entree& is)
{
  Milieu_composite::readOn(is);

  // XXX pour le moment on force ca .. a retirer apres
  if (has_saturation_)
    Process::exit("We dont accept at present a saturation object in Milieu_composite_Euler ... But we will soon !\n");

  if (static_cast<int>(fluides_.size()) == 1 && inter_lu_)
    {
      Cerr << "Error while reading Milieu_composite_Euler !!" << finl;
      Cerr << "You are simulating a Single-Phase Euler problem. No need to define an interface !" << finl;
      Cerr << "Please remove the object " << inter_lu_->que_suis_je() << " from your data file !" << finl;
      Process::exit();
    }

  if (inter_lu_ && !sub_type(Interface_Baer_Nunziato, inter_lu_.valeur()))
    Process::exit("We dont accept at present an interface object with a type different than Interface_Baer_Nunziato !\n");

  return is;
}

void Milieu_composite_Euler::discretiser(const Probleme_base& pb, const  Discretisation_base& dis)
{
  Cerr << "Composite Euler medium discretization" << finl;
  // on discretise seulement la porosite
  Milieu_base::discretiser_porosite(pb,dis);

  res_en_T_ = true;

  if (inter_lu_)
    {
      inter_lu_->assoscier_pb(pb);
      inter_lu_->completer();
    }
}

void Milieu_composite_Euler::init_energie_tot(DoubleTab& alpha_energie_tot_jdd) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler, equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& alpha = equation("alpha").inconnue().valeurs();
  const DoubleTab& p = qdm.pression().valeurs();
  const int Nb_phase = static_cast<int>(fluides_.size());
  const int Nb_elem = qdm.domaine_dis().nb_elem();

  for (int n = 0; n < Nb_phase; n++)
    {
      const DoubleTab& vit_phase = qdm.vitesse_phase(n).valeurs();
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base, get_fluid(n));

      for (int i = 0; i < Nb_elem; i++)
        {
          double nom_u2 = 0;
          for (int d = 0; d < Objet_U::dimension; d++)
            nom_u2 += vit_phase(i, d) * vit_phase(i, d);

          alpha_energie_tot_jdd(i, n) = alpha(i, n) * phase.init_energie_tot(rho(i, n), nom_u2, p(i, n));
        }
    }

  alpha_energie_tot_jdd.echange_espace_virtuel();
}

void Milieu_composite_Euler::calculer_pression(DoubleTab& p) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler, equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& alpha = equation("alpha").inconnue().valeurs();
  const DoubleTab& alpha_rhoE = equation("alpha_energie_tot").inconnue().valeurs();

  DoubleTrav rhoE(alpha_rhoE);
  rhoE = alpha_rhoE; // XXX
  tab_divide_any_shape(rhoE, alpha); // @suppress("Function cannot be resolved")

  const int Nb_phase = static_cast<int>(fluides_.size());
  const int Nb_elem = qdm.domaine_dis().nb_elem();

  for (int n = 0; n < Nb_phase; n++)
    {
      const DoubleTab& vit_phase = qdm.vitesse_phase(n).valeurs();
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base, get_fluid(n));

      for (int e = 0; e < Nb_elem; e++)
        {
          double nom_u2 = 0;
          for (int d = 0; d < Objet_U::dimension; d++)
            nom_u2 += vit_phase(e, d) * vit_phase(e, d);

          p(e, n) = phase.calculer_pression(rho(e, n), nom_u2, rhoE(e, n));
          assert(p(e, n) > 0);
        }
    }

  p.echange_espace_virtuel();
}

void Milieu_composite_Euler::calculer_vitesse_son(DoubleTab& c) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler, equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& p = qdm.pression().valeurs();

  const int Nb_phase = static_cast<int>(fluides_.size());
  const int Nb_elem = qdm.domaine_dis().nb_elem();
  assert(c.dimension(0) == Nb_elem && c.line_size() == Nb_phase);
  assert(c.dimension_tot(0) == qdm.domaine_dis().nb_elem_tot());

  for (int n = 0; n < Nb_phase; n++)
    {
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base, get_fluid(n));
      for (int e = 0; e < Nb_elem; e++)
        c(e, n) = phase.calculer_vitesse_son(rho(e, n), p(e, n));
    }

  c.echange_espace_virtuel();
}
