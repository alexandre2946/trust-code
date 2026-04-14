/****************************************************************************
* Copyright (c) 2022, CEA
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
//#include <Pb_Multiphase_HEM.h>
#include <Pb_Euler.h>
#include <Fraction_Euler.h>
#include <Champ_Uniforme.h>
#include <Matrice_Morse.h>
#include <Discret_Thyd.h>
#include <Fluide_base.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Domaine.h>
#include <EChaine.h>
#include <Param.h>
#include <Momentum_Euler.h>
//#include <Op_Conv_Rusanov_Coloc_Elem.h>

Implemente_instanciable(Fraction_Euler, "Fraction_Euler", Conservation_Euler);

Sortie& Fraction_Euler::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Fraction_Euler::readOn(Entree& is)
{
  Conservation_Euler::readOn(is);
  terme_nconserv_.associer_eqn(*this);
  return is;
}

const Operateur& Fraction_Euler::operateur(int i) const
{
  if (i)
    {
      Cerr << que_suis_je() <<" : wrong operator number " << i << finl;
      Process::exit();
    }
  return terme_nconserv_;
}

Operateur& Fraction_Euler::operateur(int i)
{
  if (i)
    {
      Cerr << que_suis_je() <<" : wrong operator number " << i << finl;
      Process::exit();
    }
  return terme_nconserv_;
}

void Fraction_Euler::discretiser()
{
  int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  double temps = schema_temps().temps_courant();
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  Cerr << "Volume fraction discretization" << finl;

  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());
  dis.discretiser_champ("temperature",domaine_dis(),"alpha","sans_dimension", pb.nb_phases(),nb_valeurs_temp,temps,l_inco_ch_);
  l_inco_ch_->fixer_nature_du_champ(pb.nb_phases() == 1 ? scalaire : pb.nb_phases() == dimension ? vectoriel : multi_scalaire); //pfft
  for (int i = 0; i < pb.nb_phases(); i++)
    l_inco_ch_->fixer_nom_compo(i, Nom("alpha_") + pb.nom_phase(i));
  champs_compris_.ajoute_champ(l_inco_ch_);
  Equation_base::discretiser();
  Cerr << "Fraction_Euler::discretiser() ok" << finl;
}


void Fraction_Euler::set_param(Param& param)
{
  Equation_base::set_param(param);
  //param.ajouter_non_std("diffusion",(this));
  //param.ajouter_non_std("convection",(this));
  param.ajouter_non_std("termes_non_conservatifs",(this));
  param.ajouter_non_std("non_conservative_terms",(this));

}
