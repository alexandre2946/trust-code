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

#include <Fraction_Euler.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <Pb_Euler.h>
#include <Domaine.h>
#include <Param.h>

Implemente_instanciable(Fraction_Euler, "Fraction_Euler", Conservation_Euler_base);
// XD fraction_euler cons_euler fraction_euler INHERITS_BRACE Void fraction conservation equation for a multi-phase
// XD_CONT Euler problem where the unknown is the temperature

Sortie& Fraction_Euler::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Fraction_Euler::readOn(Entree& is)
{
  Conservation_Euler_base::readOn(is);
  add_missing_nconserv_op();
  verifier_somme_alpha();

  terme_nconserv_.set_fichier("Non_conservative_fraction");
  terme_nconserv_.set_description("Conribution of non_conservative operator in fraction equation");

  return is;
}

void Fraction_Euler::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("termes_non_conservatifs|non_conservative_terms", (this));
}

void Fraction_Euler::discretiser()
{
  Cerr << "Fraction_Euler discretization" << finl;
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());

  const double temps = schema_temps().temps_courant();
  const int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  const int N = pb.nb_phases();

  dis.discretiser_champ("temperature", domaine_dis(), "alpha", "sans_dimension", N, nb_valeurs_temp, temps, l_inco_ch_);
  l_inco_ch_->fixer_nature_du_champ(N == 1 ? scalaire : multi_scalaire);

  for (int i = 0; i < N; i++)
    l_inco_ch_->fixer_nom_compo(i, Nom("alpha_") + pb.nom_phase(i));

  champs_compris_.ajoute_champ(l_inco_ch_);

  Equation_base::discretiser();
  Cerr << "Fraction_Euler discretization ==> ok" << finl;
}

const Operateur& Fraction_Euler::operateur(int i) const
{
  if (i)
    {
      Cerr << que_suis_je() << " : wrong operator number " << i << finl;
      Process::exit();
    }
  return terme_nconserv_;
}

Operateur& Fraction_Euler::operateur(int i)
{
  if (i)
    {
      Cerr << que_suis_je() << " : wrong operator number " << i << finl;
      Process::exit();
    }
  return terme_nconserv_;
}

void Fraction_Euler::mettre_a_jour(double temps)
{
  Conservation_Euler_base::mettre_a_jour(temps);
#ifndef NDEBUG
  verifier_somme_alpha();
#endif
}

void Fraction_Euler::verifier_somme_alpha()
{
  Cerr << "Fraction_Euler::verifier_somme_alpha() ..." ;
  const DoubleTab& vals = l_inco_ch_->valeurs();
  const int ne = vals.dimension(0), nl = vals.line_size();
  DoubleVect vals_somme(ne);
  vals_somme = 0.;

  for (int i = 0; i < ne; i++)
    for (int j = 0; j < nl; j++)
      vals_somme(i) += vals(i, j);

  const double min_a = mp_min_vect(vals_somme), max_a = mp_max_vect(vals_somme);

  if (min_a < 1. - 1.e-12 || max_a > 1. + 1.e-12)
    {
      Cerr << " KO !!! " << finl;
      Cerr << "WHAT ?? The sum of the void fraction (per cell) is not 1 !!!! You should do something !" << finl;
      Process::exit();
    }
  else
    Cerr << " OK ! " << finl;
}
