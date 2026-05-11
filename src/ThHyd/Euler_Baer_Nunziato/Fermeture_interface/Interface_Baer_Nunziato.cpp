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

#include <Interface_Baer_Nunziato.h>
#include <Pb_Euler.h>
#include <Param.h>

Implemente_instanciable(Interface_Baer_Nunziato, "Interface_Baer_Nunziato", Interface_base);
// XD Interface_Baer_Nunziato Interface_base Interface_Baer_Nunziato INHERITS_BRACE Interface Baer Nunziato class

Sortie& Interface_Baer_Nunziato::printOn(Sortie& os) const { return os; }
Entree& Interface_Baer_Nunziato::readOn(Entree& is)
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Interface_Baer_Nunziato::set_param(Param& param) const
{
  param.ajouter("velocity_from_phase|phase_vitesse", &nom_phase_vitesse_ , Param::REQUIRED); // XD_ADD_P chaine Name of phase used to take the velocity at the interface. The pressure is automatically taken from the other phase
}

void Interface_Baer_Nunziato::completer()
{
  assert (pb_);
  if (!sub_type(Pb_Euler, pb_.valeur()))
    {
      Cerr << "Interface_Baer_Nunziato should only be used with a problem of type Pb_Euler, not " << pb_->que_suis_je() << " !!!!" << finl;
      Process::exit();
    }

  if (!sub_type(Milieu_composite_Euler, pb_->milieu()))
    {
      Cerr << "Milieu_composite_Euler should only be used with a problem of type Pb_Euler, not " << pb_->milieu().que_suis_je() << " !!!!" << finl;
      Process::exit();
    }

  const Milieu_composite_Euler& mil = ref_cast(Milieu_composite_Euler, pb_->milieu());
  const auto& noms_phases = mil.noms_phases();

  if (noms_phases.size() != 2)
    {
      Cerr << "Error in Interface_Baer_Nunziato : the class is actually coded only for 2 phase problem. We detect the following phases :\n " << noms_phases << finl;
      Process::exit();
    }

  for (int i = 0; i < noms_phases.size(); i++)
    {
      if (nom_phase_vitesse_ == noms_phases[i])
        {
          id_vitesse_interface_ = i;
          id_pression_interface_ = 1 - i; // the other is for pressure
          break;
        }
    }

  if (id_vitesse_interface_ == -123 && id_pression_interface_ == -123)
    {
      Cerr << "Error in Interface_Baer_Nunziato : the phase name " << nom_phase_vitesse_ << " is not found in the provided phases :\n " << noms_phases << finl;
      Process::exit();
    }
}
