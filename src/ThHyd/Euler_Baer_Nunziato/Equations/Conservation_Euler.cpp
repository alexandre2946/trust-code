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

#include <EcritureLectureSpecial.h>
#include <Pb_Euler.h>
#include <Conservation_Euler.h>
#include <Champ_Uniforme.h>
#include <Matrice_Morse.h>
#include <Discret_Thyd.h>
#include <Fluide_base.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Domaine.h>
#include <EChaine.h>
#include <Param.h>

Implemente_instanciable(Conservation_Euler, "Conservation_Euler", Convection_Diffusion_std);

Sortie& Conservation_Euler::printOn(Sortie& is) const
{
  return Equation_base::printOn(is);
}

Entree& Conservation_Euler::readOn(Entree& is)
{
  assert(l_inco_ch_.non_nul());
  assert(le_fluide_.non_nul());
  champs_compris_.ajoute_champ(l_inco_ch_);
  Equation_base::readOn(is);
  return is;
}

int Conservation_Euler::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      is >> terme_diffusif;
      terme_diffusif.associer_diffusivite_pour_pas_de_temps(diffusivite_pour_pas_de_temps());
      return 1;
    }
  else if (mot == "convection")
    {
      Cerr << "Reading and typing of the convection operator : " << finl;
      const Champ_base& ch_vitesse_transportante = vitesse_pour_transport();
      associer_vitesse(ch_vitesse_transportante);
      terme_convectif.associer_vitesse(ch_vitesse_transportante);
      is >> terme_convectif;
      terme_convectif.associer_eqn(*this);
      return 1;
    }
  else if (mot == "termes_non_conservatifs|non_conservative_terms")
    {
      Cerr << "Reading and typing of the non_conservative_terms operator : " << finl;
      is >> terme_nconserv_;
      terme_nconserv_.associer_eqn(*this);
      return 1;
    }
  else
    return Equation_base::lire_motcle_non_standard(mot, is);
}

void Conservation_Euler::associer_milieu_base(const Milieu_base& un_milieu) //ok
{
  const Fluide_base& un_fluide = ref_cast(Fluide_base, un_milieu);
  associer_fluide(un_fluide);
}

const Milieu_base& Conservation_Euler::milieu() const
{
  return fluide();
}

Milieu_base& Conservation_Euler::milieu()
{
  return fluide();
}

const Fluide_base& Conservation_Euler::fluide() const
{
  if (le_fluide_.est_nul())
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

Fluide_base& Conservation_Euler::fluide()
{
  assert(le_fluide_.non_nul());
  return le_fluide_.valeur();
}

void Conservation_Euler::associer_fluide(const Fluide_base& un_fluide)
{
  assert(sub_type(Fluide_base,un_fluide));
  le_fluide_ = ref_cast(Fluide_base, un_fluide);
}

