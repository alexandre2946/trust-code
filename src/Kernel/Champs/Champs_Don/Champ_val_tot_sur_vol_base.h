/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef Champ_val_tot_sur_vol_base_included
#define Champ_val_tot_sur_vol_base_included

#include <Champ_Uniforme_Morceaux.h>

class Domaine_dis_base;
class Domaine_Cl_dis_base;

/*! @brief Champ_val_tot_sur_vol_base Base class derived from Champ_Uniforme_Morceaux representing fields
 *
 *      whose evaluation is expressed as val_lue_loc/(Somme_vol_poro_loc)
 *
 *      val_lue_loc denotes the value read for a spatial localization (subdomain or default domain)
 *      VDF case: Somme_vol_poro_loc is the sum of vol_element*poro_volumique
 *                for elements contained in localization loc
 *      VEF case: Somme_vol_poro_loc is the sum of vol_entrelaces*poro_surface
 *                for faces contained in localization loc
 *      Somme_vol_poro_loc is evaluated by the eval_contrib_loc() method of derived classes
 *
 *  The user syntax to follow is:
 *                Valeur_totale_sur_volume nom_domaine nb_comp { defaut val_lue_dom ...domainei val_lue_szi ... }
 *                nom_domaine : name of the computational domain
 *                  nb_comp     : number of field components
 *                  val_lue_dom : value read for the default domain
 *                  val_lue_szi : value read for subdomain domainei
 *
 *
 *
 *
 */
class Champ_val_tot_sur_vol_base : public Champ_Uniforme_Morceaux
{

  Declare_base(Champ_val_tot_sur_vol_base);

public:

  void evaluer(const Domaine_dis_base& zdis,const Domaine_Cl_dis_base& zcldis);
  virtual DoubleVect& eval_contrib_loc(const Domaine_dis_base& zdis,const Domaine_Cl_dis_base& zcldis,
                                       DoubleVect& vol)=0;
};

#endif
