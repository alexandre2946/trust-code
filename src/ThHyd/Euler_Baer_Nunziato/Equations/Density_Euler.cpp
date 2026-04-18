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
#include <Momentum_Euler.h>
#include <Density_Euler.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <Pb_Euler.h>
#include <Domaine.h>
#include <Param.h>

Implemente_instanciable(Density_Euler, "Masse_Euler|Density_Euler", Conservation_Euler);
// XD masse_euler cons_euler density_euler -1 Mass consevation equation for a multi-phase Euler problem where the unknown is the alpha (void fraction)

Sortie& Density_Euler::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Density_Euler::readOn(Entree& is)
{
  Conservation_Euler::readOn(is);
  assert(densite_.non_nul());
  terme_convectif.associer_eqn(*this);
  terme_convectif.set_fichier("Debit");
  terme_convectif.set_description((Nom) "Mass flow rate=Integral(-rho*u*ndS) [kg/s] if SI units used");
  return is;
}

void Density_Euler::discretiser()
{
  int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  double temps = schema_temps().temps_courant();
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());

  Cerr << "Density_Euler discretization" << finl;

  dis.discretiser_champ("temperature", domaine_dis(), "alpha_rho", "kg/m3", pb.nb_phases(), nb_valeurs_temp, temps, l_inco_ch_);
  l_inco_ch_->fixer_nature_du_champ(pb.nb_phases() == 1 ? scalaire : pb.nb_phases() == dimension ? vectoriel : multi_scalaire); //pfft
  for (int i = 0; i < pb.nb_phases(); i++)
    l_inco_ch_->fixer_nom_compo(i, Nom("alpha_rho_") + pb.nom_phase(i));
  champs_compris_.ajoute_champ(l_inco_ch_);

  dis.discretiser_champ("temperature", domaine_dis(), "densite", "kg/m3", pb.nb_phases(), nb_valeurs_temp, temps, densite_);
  l_inco_ch_->fixer_nature_du_champ(pb.nb_phases() == 1 ? scalaire : pb.nb_phases() == dimension ? vectoriel : multi_scalaire); //pfft
  for (int i = 0; i < pb.nb_phases(); i++)
    densite_->fixer_nom_compo(i, Nom("densite_") + pb.nom_phase(i));
  champs_compris_.ajoute_champ(densite_);

  Equation_base::discretiser();
  Cerr << "Density_Euler::discretiser() ok" << finl;
}

Entree& Density_Euler::lire_cond_init(Entree& is)
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
  if (motlu != Motcle(densite().le_nom()))
    {
      Cerr << que_suis_je() << " : expected " << densite().le_nom() << " instead of " << nom << finl;
      exit();
    }
  OWN_PTR(Champ_Don_base) ch_init;
  is >> ch_init;
  const int nb_comp = ch_init->nb_comp();
  verifie_ch_init_nb_comp(densite(), nb_comp);

  densite().affecter(ch_init.valeur());
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

void Density_Euler::init_alpha_rho()
{
  const DoubleTab& alpha = ref_cast(Pb_Euler,probleme()).equation_fraction().inconnue().valeurs();
  const DoubleTab& rho = densite().valeurs();
  DoubleTab& alpha_rho = inconnue().valeurs();
  alpha_rho = rho;
  tab_multiply_any_shape(alpha_rho, alpha);
}

void Density_Euler::mettre_a_jour_champs_conserves(double temps, int reset)
{
  Equation_base::mettre_a_jour_champs_conserves(temps);
  const DoubleTab& alpha = ref_cast(Pb_Euler,probleme()).equation_fraction().inconnue().valeurs();
  const DoubleTab& alpha_rho = inconnue().valeurs();
  DoubleTab& rho = densite().valeurs();
  rho = alpha_rho;
  tab_divide_any_shape(rho, alpha);
}

DoubleTab Density_Euler::flux(const int f, const int left_or_right) const
{
  //left_or_right = 0 : left et 1 right;
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const DoubleTab& vit_normale = ref_cast(Momentum_Euler, probleme().equation(0)).vitesse_normale();
  const DoubleTab& alpha_rho = inconnue().valeurs();
  const int e = dom.face_voisins(f, left_or_right);
  const int nb_phases = ref_cast(Pb_Euler,probleme()).nb_phases();
  DoubleTrav flux_(nb_phases);
  for (int n = 0; n < nb_phases; n++)
    flux_(n) = alpha_rho(e, n) * vit_normale(f, n + left_or_right * nb_phases);
  return flux_;
}

const Operateur& Density_Euler::operateur(int i) const
{
  if (i)
    {
      Cerr << que_suis_je() << " : wrong operator number " << i << finl;
      Process::exit();
    }
  return terme_convectif;
}

Operateur& Density_Euler::operateur(int i)
{
  if (i)
    {
      Cerr << que_suis_je() << " : wrong operator number " << i << finl;
      Process::exit();
    }
  return terme_convectif;
}

void Density_Euler::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("convection", (this));
}
