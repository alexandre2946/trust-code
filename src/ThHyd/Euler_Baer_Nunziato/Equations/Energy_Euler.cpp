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

#include <Momentum_Euler.h>
#include <Energy_Euler.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <Pb_Euler.h>
#include <Domaine.h>
#include <Param.h>

Implemente_instanciable(Energy_Euler, "Energy_Euler|Energie_Euler", Conservation_Euler_base);
// XD energy_euler cons_euler energie_euler INHERITS_BRACE Internal energy conservation equation for a multi-phase Euler problem where the unknown is the temperature

Sortie& Energy_Euler::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Energy_Euler::readOn(Entree& is)
{
  Conservation_Euler_base::readOn(is);
  add_missing_nconserv_op();

  terme_convectif.set_fichier("Convection_energie");
  terme_convectif.set_description("Convective heat transfer rate=Integral(-h*u*ndS) [W] if SI units used");

  terme_nconserv_.set_fichier("Non_conservative_energie");
  terme_nconserv_.set_description("Conribution of non_conservative operator in energy equation");

  return is;
}

void Energy_Euler::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("termes_non_conservatifs|non_conservative_terms", (this));
  param.ajouter_non_std("convection", (this), Param::REQUIRED);
}

void Energy_Euler::discretiser()
{
  Cerr << "Energy_Euler discretization" << finl;
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());

  const double temps = schema_temps().temps_courant();
  const int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  const int N = pb.nb_phases();

  dis.discretiser_champ("temperature", domaine_dis(), "alpha_energie_tot", "J/m3", N, nb_valeurs_temp, temps, l_inco_ch_);
  l_inco_ch_->fixer_nature_du_champ(N == 1 ? scalaire : multi_scalaire);

  for (int i = 0; i < N; i++)
    l_inco_ch_->fixer_nom_compo(i, Nom("alpha_energie_tot_") + pb.nom_phase(i));

  champs_compris_.ajoute_champ(l_inco_ch_);

  Equation_base::discretiser();

  Cerr << "Energy_Euler discretization ==> ok" << finl;
}

// on surcharge pour pas effecter energie_tot et pour lire mais rien faire !!
Entree& Energy_Euler::lire_cond_init(Entree& is)
{
  Cerr << "Reading of initial conditions\n";
  Nom nom;
  Motcle motlu;
  is >> nom;
  motlu = nom;
  if (motlu != Motcle("{"))
    {
      Cerr << "We expected a { while reading " << que_suis_je() << finl;
      Cerr << "and not : " << nom << finl;
      exit();
    }
  is >> nom;
  motlu = nom;
  if (motlu != Motcle(inconnue().le_nom()) && motlu != "Energie_tot")
    {
      Cerr << nom << " is not the name of the unknown " << inconnue().le_nom() << finl;
      exit();
    }
  OWN_PTR(Champ_Don_base) ch_init;
  is >> ch_init;

  const int nb_comp = ch_init->nb_comp();
  verifie_ch_init_nb_comp(inconnue(), nb_comp);

  is >> nom;
  motlu = nom;
  if (motlu != Motcle("}"))
    {
      Cerr << "We expected a } while reading " << que_suis_je() << finl;
      Cerr << "and not : " << nom << finl;
      exit();
    }
  return is;
}

const Operateur& Energy_Euler::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_convectif;
    case 1:
      return terme_nconserv_;
    default:
      Cerr << que_suis_je() << " : wrong operator number " << i << finl;
      Process::exit();
    }
  // Pour les compilos!!
  return terme_convectif;
}

Operateur& Energy_Euler::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_convectif;
    case 1:
      return terme_nconserv_;
    default:
      Cerr << que_suis_je() << " : wrong operator number " << i << finl;
      Process::exit();
    }
  // Pour les compilos!!
  return terme_convectif;
}
