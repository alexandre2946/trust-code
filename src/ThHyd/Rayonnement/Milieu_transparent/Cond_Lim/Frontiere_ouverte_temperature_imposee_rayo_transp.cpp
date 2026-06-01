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

#include <Frontiere_ouverte_temperature_imposee_rayo_transp.h>
#include <Pb_Fluide_base.h>
#include <Front_VF.h>

Implemente_instanciable(Frontiere_ouverte_temperature_imposee_rayo_transp, "Frontiere_ouverte_temperature_imposee_rayo_transp", Entree_fluide_temperature_imposee);

// XD frontiere_ouverte_temperature_imposee_rayo_transp frontiere_ouverte_temperature_imposee frontiere_ouverte_temperature_imposee_rayo_transp INHERITS_BRACE Radiation imposed temperature condition at the open boundary called bord (edge) (in the case of fluid inlet). This condition must be associated with an imposed inlet velocity condition. The imposed temperature value is expressed in C or K.

Sortie& Frontiere_ouverte_temperature_imposee_rayo_transp::printOn(Sortie& is) const { return is; }

Entree& Frontiere_ouverte_temperature_imposee_rayo_transp::readOn(Entree& s) { return Entree_fluide_temperature_imposee::readOn(s); }

int Frontiere_ouverte_temperature_imposee_rayo_transp::initialiser(double temps)
{
  assert(!le_modele_rayo_);

  // retrieve the radiation model ... !
  const Probleme_base& this_pb = domaine_Cl_dis().equation().probleme();
  if (sub_type(Pb_Fluide_base, this_pb))
    {
      if (this_pb.milieu().is_rayo_transp())
        {
          le_modele_rayo_ = ref_cast(Pb_Fluide_base, this_pb).get_mod_rayo_transp();

          if (le_modele_rayo_->nom_pb_rayonnant() != this_pb.le_nom())
            error_pb_name(que_suis_je(), this_pb.le_nom(), le_modele_rayo_->nom_pb_rayonnant());
        }
      else
        error_non_rad_bc(que_suis_je(), this_pb.le_nom(), frontiere_dis().frontiere().le_nom(), "frontiere_ouverte_temperature_imposee");
    }
  else
    {
      Cerr << "The BC " << que_suis_je() << " should be associated to a fluid problem and not a a one of type " << this_pb.que_suis_je() << " !!!" << finl;
      Process::exit();
    }

  return Entree_fluide_temperature_imposee::initialiser(temps);
}

void Frontiere_ouverte_temperature_imposee_rayo_transp::completer()
{
  Entree_fluide_temperature_imposee::completer();
  preparer_surface(frontiere_dis(), domaine_Cl_dis());
}

void Frontiere_ouverte_temperature_imposee_rayo_transp::mettre_a_jour(double temps)
{
  Entree_fluide_temperature_imposee::mettre_a_jour(temps);
  calculer_Teta_i();
}

void Frontiere_ouverte_temperature_imposee_rayo_transp::calculer_Teta_i()
{
  const Front_VF& front_vf = ref_cast(Front_VF, frontiere_dis());
  for (int numfa = 0; numfa < front_vf.nb_faces(); numfa++)
    teta_i_[numfa] = val_imp(numfa);
}
