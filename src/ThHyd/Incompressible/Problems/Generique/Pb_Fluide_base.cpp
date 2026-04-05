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
#include <Pb_Couple_rayo_transp.h>
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

  // test si c'est un pb couple de type Pb_Couple_rayo_transp ou pas !
  if (is_coupled() && !sub_type(Pb_Couple_rayo_transp, pbc_.valeur()))
    {
      Cerr << "You asked for using a transparent medium radiation model with a coupled problem of type " << pbc_->que_suis_je() << finl;
      Cerr << "Please update your data file by using instead a coupled problem of type Pb_Couple_rayo_transp ..." << finl;
      Process::exit();
    }

  // set flag is_rad_transp_med_ in Fluide_base
  bool flag_set = false;
  for (auto& itr : le_milieu_)
    {
      if (sub_type(Fluide_base, itr.valeur()))
        {
          ref_cast(Fluide_base, itr.valeur()).set_rayo_transp_flag();
          flag_set = true;
          break;
        }
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

void Pb_Fluide_base::preparer_calcul()
{
  if (mod_rayo_transp_.non_nul() && !is_coupled()) // sinon c'est fait dans Pb_Couple_rayo_transp ...
    assoscier_rayo_model_CL();

  Probleme_base::preparer_calcul();
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

void Pb_Fluide_base::assoscier_rayo_model_CL()
{
  if (mod_rayo_transp_.est_nul()) return; /* rien a faire */

  // TODO FIXME for now we suppose that we have only one radiation model in the coupled pb ...
  // see test in Pb_Couple_rayo_transp::initialize
  const int nb_pbs = is_coupled() ? pbc_->nb_problemes() : 1;
  for (int l = 0; l < nb_pbs; l++)
    {
      Probleme_base& le_pb = is_coupled() ? ref_cast(Probleme_base, pbc_->probleme(l)) : *this;

      for (int j = 0; j < le_pb.nombre_d_equations(); j++)
        {
          Domaine_Cl_dis_base& la_zcl = le_pb.equation(j).domaine_Cl_dis();
          for (int num_cl = 0; num_cl < la_zcl.nb_cond_lim(); num_cl++)
            {
              Cond_lim_base& la_cl = la_zcl.les_conditions_limites(num_cl).valeur();

              Cond_lim_rayo_milieu_transp *la_cl_rayo;
              if (la_cl.is_bc_rayo_milieu_transp(la_cl_rayo))
                la_cl_rayo->associer_modele_rayo(mod_rayo_transp_.valeur());
            }
        }
    }

  mod_rayo_transp_->preparer_calcul();
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
