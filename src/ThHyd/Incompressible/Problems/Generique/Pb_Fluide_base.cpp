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

#include <Cond_lim_rayo_milieu_transp.h>
#include <Pb_Fluide_base.h>
#include <Fluide_base.h>

Implemente_base(Pb_Fluide_base, "Pb_Fluide_base", Probleme_base);

Sortie& Pb_Fluide_base::printOn(Sortie& os) const { return Probleme_base::printOn(os); }
Entree& Pb_Fluide_base::readOn(Entree& is) { return Probleme_base::readOn(is); }

Entree& Pb_Fluide_base::lire_radiation_models(Entree& is, Motcle& mot)
{
  assert (mot == "Modele_rayonnement_milieu_transparent" || mot == "Transparent_medium_radiation_model");

  // TODO FIXME do better here ...
  if (!(Motcle(que_suis_je()).debute_par("Pb_HYDRAULIQUE") || Motcle(que_suis_je()).debute_par("Pb_THERMOHYDRAULIQUE") ))
    {
      Cerr << "The transparent medium radiation model is not yet tested with a problem of type " << que_suis_je() << finl;
      Cerr << "Please contact the TRUST team." << finl;
      Process::exit();
    }

  // set flag is_rad_transp_med_ in Fluide_base
  bool flag_set = false;
  for (auto& itr : le_milieu_)
    if (sub_type(Fluide_base, itr.valeur()))
      {
        ref_cast(Fluide_base, itr.valeur()).set_rayo_transp_flag();
        flag_set = true;
        break;
      }

  if (!flag_set)
    {
      Cerr << "Using a transparent medium radiation model with a problem of type " << que_suis_je() << " that dont have a fluid medium is forbidden !!!" << finl;
      Process::exit();
    }

  // si bon on type et on lit !
  mod_rayo_transp_.typer(mot.getChar());
  is >> mod_rayo_transp_.valeur();
  mod_rayo_transp_->associer_pb_fluide_rayo(*this);

  return is;
}

void Pb_Fluide_base::completer()
{
  Probleme_base::completer();

  if (mod_rayo_transp_.non_nul())
    mod_rayo_transp_->completer();
}

int Pb_Fluide_base::postraiter(int force)
{
  int ok = Probleme_base::postraiter(force);

  if (!ok)
    return 0;

  if (mod_rayo_transp_.non_nul())
    mod_rayo_transp_->postraiter();

  return ok;
}

void Pb_Fluide_base::validateTimeStep()
{
  Probleme_base::validateTimeStep();

  if (mod_rayo_transp_.non_nul())
    mod_rayo_transp_->mettre_a_jour(presentTime());
}

int Pb_Fluide_base::expression_predefini(const Motcle& motlu, Nom& expression)
{
  if (motlu=="ENERGIE_CINETIQUE_TOTALE")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " energie_cinetique_totale } ";
      return 1;
    }
  else if (motlu=="ENERGIE_CINETIQUE_ELEM")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " energie_cinetique_elem } ";
      return 1;
    }
  else if (motlu=="VISCOUS_FORCE_X")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " viscous_force_X } ";
      return 1;
    }
  else if (motlu=="VISCOUS_FORCE_Y")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " viscous_force_y } ";
      return 1;
    }
  else if (motlu=="VISCOUS_FORCE_Z")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " viscous_force_z } ";
      return 1;
    }
  else if (motlu=="VISCOUS_FORCE")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " viscous_force } ";
      return 1;
    }
  else if (motlu=="PRESSURE_FORCE_X")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " pressure_force_x } ";
      return 1;
    }
  else if (motlu=="PRESSURE_FORCE_Y")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " pressure_force_y } ";
      return 1;
    }
  else if (motlu=="PRESSURE_FORCE_Z")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " pressure_force_z } ";
      return 1;
    }
  else if (motlu=="PRESSURE_FORCE")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " pressure_force } ";
      return 1;
    }
  else if (motlu=="TOTAL_FORCE_X")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " total_force_x } ";
      return 1;
    }
  else if (motlu=="TOTAL_FORCE_Y")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " total_force_y } ";
      return 1;
    }
  else if (motlu=="TOTAL_FORCE_Z")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " total_force_z } ";
      return 1;
    }
  else if (motlu=="TOTAL_FORCE")
    {
      expression = "predefini { pb_champ ";
      expression += le_nom();
      expression += " total_force } ";
      return 1;
    }
  else
    return Probleme_base::expression_predefini(motlu,expression);
}
