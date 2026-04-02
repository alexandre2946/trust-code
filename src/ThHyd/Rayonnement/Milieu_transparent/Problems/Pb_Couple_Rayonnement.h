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

#ifndef Pb_Couple_Rayonnement_included
#define Pb_Couple_Rayonnement_included

#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Probleme_Couple.h>
#include <TRUST_Ref.h>

class Cond_lim_base;
class Schema_Temps_base;
class Discretisation_base;

class Pb_Couple_Rayonnement: public Probleme_Couple
{
  Declare_instanciable(Pb_Couple_Rayonnement);
public:
  int associer_(Objet_U&) override;
  int postraiter(int force = 1) override;

  void completer();
  void validateTimeStep() override;
  void initialize() override;

  inline Modele_Rayonnement_Milieu_Transparent& le_modele_rayo() { return le_modele_de_rayo_.valeur(); }

protected:
  OBS_PTR(Modele_Rayonnement_Milieu_Transparent) le_modele_de_rayo_;
};

#endif /* Pb_Couple_Rayonnement_included */
