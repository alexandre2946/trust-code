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

#include <Temperature_imposee_paroi_rayo_transp.h>
#include <Champ_front_contact_rayo_transp_VEF.h>
#include <Pb_Fluide_base.h>
#include <Equation_base.h>
#include <Domaine_VEF.h>

Implemente_instanciable(Temperature_imposee_paroi_rayo_transp, "Paroi_temperature_imposee_rayo_transp", Temperature_imposee_paroi);

Sortie& Temperature_imposee_paroi_rayo_transp::printOn(Sortie& is) const { return is; }

Entree& Temperature_imposee_paroi_rayo_transp::readOn(Entree& s) { return Temperature_imposee_paroi::readOn(s); }

int Temperature_imposee_paroi_rayo_transp::initialiser(double temps)
{
  assert(le_modele_rayo_.est_nul());

  // on recupere le modele rayo seulement si pb fluide et rayo ... !
  const Probleme_base& this_pb = domaine_Cl_dis().equation().probleme();
  if (sub_type(Pb_Fluide_base, this_pb)) // sinon Pb_conduction par exemple ;)
    {
      if (this_pb.milieu().is_rayo_transp())
        {
          le_modele_rayo_ = ref_cast(Pb_Fluide_base, this_pb).get_mod_rayo_transp();

          if (le_modele_rayo_->nom_pb_rayonnant() != this_pb.le_nom())
            error_pb_name(que_suis_je(), this_pb.le_nom(), le_modele_rayo_->nom_pb_rayonnant());
        }
      else
        error_non_rad_bc(que_suis_je(), this_pb.le_nom(), frontiere_dis().frontiere().le_nom(), "temperature_imposee_paroi");
    }

  return Temperature_imposee_paroi::initialiser(temps);
}

void Temperature_imposee_paroi_rayo_transp::completer()
{
  Temperature_imposee_paroi::completer();
  preparer_surface(frontiere_dis(), domaine_Cl_dis());
}

void Temperature_imposee_paroi_rayo_transp::calculer_Teta_i(double temps)
{
  if (sub_type(Champ_front_contact_rayo_transp_VEF, le_champ_front.valeur()))
    {
      Champ_front_contact_rayo_transp_VEF& Ch_contact = ref_cast(Champ_front_contact_rayo_transp_VEF, le_champ_front.valeur());
      Ch_contact.calculer_temperature_bord(temps);
    }
  else
    {
      // La temperature de paroi etant directement donnee par le champ_front associe a la condition a la limite, il n'y a rien a calculer ici
    }

  const Front_VF& front_vf = ref_cast(Front_VF, frontiere_dis());
  for (int numfa = 0; numfa < front_vf.nb_faces(); numfa++)
    teta_i_[numfa] = val_imp(numfa);
}
