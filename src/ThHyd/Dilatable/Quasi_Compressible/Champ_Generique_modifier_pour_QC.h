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


#ifndef Champ_Generique_modifier_pour_QC_included
#define Champ_Generique_modifier_pour_QC_included

#include <Champ_Gen_de_Champs_Gen.h>
#include <TRUST_Ref.h>

class Milieu_base;

/*! @brief class Champ_Generique_modifier_pour_QC OWN_PTR(Champ_base) intended to post-process a field of a quasi-compressible problem
 *
 *  that we want to multiply or divide by the density (rho).
 *  The class carries a REF to the medium which must be of type Fluide_Quasi_Compressible.
 *
 */

//// Data file syntax to follow
//
// "field_name" modifier_pour_QC {
//                 [division]
//                source "generic_field_type" { ... source ref_Champ { Pb_champ "pb_name" "discrete_field_name" } }
//               }
// "field_name" set by user will be the name of the generic field
// "division" activates the division of the field by rho otherwise the field is multiplied by rho
// "generic_field_type" type of a generic field

class Champ_Generique_modifier_pour_QC : public Champ_Gen_de_Champs_Gen
{
  Declare_instanciable_sans_constructeur(Champ_Generique_modifier_pour_QC);

public:

  Champ_Generique_modifier_pour_QC();
  void set_param(Param& param) const override;
  void completer(const Postraitement_base& post) override;
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  void nommer_source() override;

protected:

  bool diviser_ = false;
  OBS_PTR(Milieu_base) mon_milieu_;
};

#endif

