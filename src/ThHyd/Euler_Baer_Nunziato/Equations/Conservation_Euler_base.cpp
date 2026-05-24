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

#include <Conservation_Euler_base.h>
#include <Pb_Euler.h>
#include <EChaine.h>

Implemente_base(Conservation_Euler_base, "Conservation_Euler_base", Convection_Diffusion_std);
// XD cons_euler eqn_base cons_euler INHERITS_BRACE Base class equation for a multi-phase Euler conservation equations
// XD attr termes_non_conservatifs bloc_op_non_conservativtifs non_conservative_terms OPT Keyword to alter the
// XD_CONT non-conservative scheme.

// XD bloc_op_non_conservativtifs objet_lecture nul NO_BRACE not_set
// XD attr aco chaine(into=["{"]) aco REQ Opening curly bracket.
// XD attr operateur op_non_conservativtifs_deriv operateur REQ not_set
// XD attr acof chaine(into=["}"]) acof REQ Closing curly bracket.

Sortie& Conservation_Euler_base::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Conservation_Euler_base::readOn(Entree& is)
{
  assert(l_inco_ch_ && le_fluide_);
  Equation_base::readOn(is);

  if (!sub_type(Pb_Euler, probleme()))
    {
      Cerr << "Equation " << que_suis_je() << " can only used with a pb of type Pb_Euler not " << probleme().que_suis_je() << " !!" << finl;
      Process::exit();
    }
  return is;
}

int Conservation_Euler_base::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "convection")
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

void Conservation_Euler_base::add_missing_nconserv_op()
{
  assert(sub_type(Pb_Euler, probleme()));
  // si monophasique et termes_non_conservatifs non-lu ... on type sans error !
  const bool is_single_phase = (ref_cast(Pb_Euler, probleme()).nb_phases() == 1);
  if (is_single_phase && !terme_nconserv_)
    {
      EChaine str(" { negligeable } ");
      str >> terme_nconserv_;
      terme_nconserv_.associer_eqn(*this);
    }
  else if (!is_single_phase && !terme_nconserv_)
    {
      Cerr << "Error while reading " << que_suis_je() << " !!! non_conservative_terms operator is not read although it is required !!! Fix your data file ..." << finl;
      Process::exit();
    }
}

void Conservation_Euler_base::associer_milieu_base(const Milieu_base& un_milieu)
{
  associer_fluide(ref_cast(Fluide_base, un_milieu));
}

void Conservation_Euler_base::associer_fluide(const Fluide_base& un_fluide)
{
  assert(sub_type(Fluide_base,un_fluide));
  le_fluide_ = ref_cast(Fluide_base, un_fluide);
}
