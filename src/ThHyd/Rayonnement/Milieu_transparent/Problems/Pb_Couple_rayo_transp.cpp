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

#include <Modele_rayo_transp.h>
#include <Pb_Couple_rayo_transp.h>
#include <Paroi_rayo_transp.h>
#include <Probleme_base.h>
#include <Fluide_base.h>

Implemente_instanciable(Pb_Couple_rayo_transp, "Pb_Couple_rayo_transp", Probleme_Couple);

Entree& Pb_Couple_rayo_transp::readOn(Entree& is) { return is; }

Sortie& Pb_Couple_rayo_transp::printOn(Sortie& os) const { return Probleme_Couple::printOn(os); }

void Pb_Couple_rayo_transp::initialize()
{
  assoscier_rayo_model_CL();
  le_modele_de_rayo_->preparer_calcul();
  Probleme_Couple::initialize();
}

void Pb_Couple_rayo_transp::associer_modele_rayo_transp(const Modele_rayo_transp& mod)
{
  if (le_modele_de_rayo_.non_nul())
    {
      Cerr << "Error in Pb_Couple_rayo_transp::associer_modele_rayo_transp. It seems that you have another model associated to the problem " << le_nom() << finl;
      Process::exit();
    }
  le_modele_de_rayo_ = mod;
}

int Pb_Couple_rayo_transp::postraiter(int force)
{
  int ok = Probleme_Couple::postraiter(force);
  if (!ok)
    return 0;

  return le_modele_de_rayo_->postraiter();
}

void Pb_Couple_rayo_transp::validateTimeStep()
{
  Probleme_Couple::validateTimeStep();
  le_modele_de_rayo_->mettre_a_jour(presentTime());
}

void Pb_Couple_rayo_transp::assoscier_rayo_model_CL()
{
  for (int l = 0; l < nb_problemes(); l++)
    {
      Probleme_base& le_pb = ref_cast(Probleme_base, probleme(l));

      for (int j = 0; j < le_pb.nombre_d_equations(); j++)
        {
          Domaine_Cl_dis_base& la_zcl = le_pb.equation(j).domaine_Cl_dis();
          for (int num_cl = 0; num_cl < la_zcl.nb_cond_lim(); num_cl++)
            {
              Cond_lim_base& la_cl = la_zcl.les_conditions_limites(num_cl).valeur();

              Cond_lim_rayo_milieu_transp *la_cl_rayo;
              if (la_cl.is_bc_rayo_milieu_transp(la_cl_rayo))
                la_cl_rayo->associer_modele_rayo(le_modele_de_rayo_.valeur());
            }
        }
    }
}
