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


#ifndef Champ_Generique_Reduction_0D_included
#define Champ_Generique_Reduction_0D_included

#include <Champ_Gen_de_Champs_Gen.h>

/*! @brief class Champ_Generique_Reduction_0D
 *
 */

// Field intended to post-process a field reduced to 0D dimension
// We construct a field taking the reduced value at every point in space
// The class carries the type of method to perform the reduction
//// Data file syntax to follow
//
// "field_name" Reduction_0D { method "method_type"
//                source "generic_field_type" { ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name" } }
//                [ source_name "source_name" ]
//                [ reference_source "reference_source" ]
//                [ reference_sources "comma_separated_list" ]
//                [ sources "list_generic_field" ]
//               }
// "field_name" set by the user will be the name of the generic field
// "generic_field_type" type of a generic field
// "method_type" indicates the type of reduction requested
// The possible options are available in $TRUST_ROOT/src/Kernel/Champs/Champ_Generique_Reduction_0D.cpp

class Champ_Generique_Reduction_0D : public Champ_Gen_de_Champs_Gen
{

  Declare_instanciable(Champ_Generique_Reduction_0D);

public:

  void set_param(Param& param) const override;
  void completer(const Postraitement_base& post) override;
  const Noms get_property(const Motcle& query) const override;
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  void nommer_source() override;
  void extraire(double& val_extraites, const DoubleVect& val_source, const bool basis_function, const int composante_VDF=-1) const;
  const Motcle get_directive_pour_discr() const override;

protected:

  Motcle methode_; //Type of reduction: min, max, moyenne ou somme
  int numero_proc_ = -10; // number of the proc containing the leftmost cell
  int numero_elem_ = -10; // local number of the leftmost cell on numero_proc_
  mutable DoubleVect volume_controle_; //Work array

private:
  mutable OWN_PTR(Champ_Fonc_base) espace_stockage_;
  mutable OWN_PTR(Champ_base) source_espace_stockage_;
  mutable DoubleVect un_;
};

#endif

