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

#ifndef Pb_Couple_rayo_semi_transp_included
#define Pb_Couple_rayo_semi_transp_included

#include <Schema_Temps_base.h>
#include <Probleme_Couple.h>
#include <Domaine.h>
#include <TRUST_Ref.h>

class Cond_lim_rayo_semi_transp;
class Modele_rayo_semi_transp;
class Cond_lim_base;

/*! @brief classe derivee de Probleme_Couple Cette classe couple, outre des Probleme_base, un modele de
 *
 *      rayonnement.
 *      Seule l'etape d'initialisation impose actuellement que le
 *      modele de rayonnement soit differencie, et donc l'existence de
 *      la classe Pb_Couple_rayo_semi_transp.
 *
 *
 * @sa Probleme_Couple Modele_rayo_semi_transp
 */
class Pb_Couple_rayo_semi_transp: public Probleme_Couple
{
  Declare_instanciable(Pb_Couple_rayo_semi_transp);
public:

  void initialize() override;

  void le_modele_rayo_associe(const Modele_rayo_semi_transp&);
  int associer_(Objet_U&) override;

  inline Modele_rayo_semi_transp& modele() { return le_modele_.valeur(); }
  inline const Modele_rayo_semi_transp& modele() const { return le_modele_.valeur(); }

protected:
  DerObjU der_domaine_clone;
  OBS_PTR(Modele_rayo_semi_transp) le_modele_;
};

#endif /* Pb_Couple_rayo_semi_transp_included */

