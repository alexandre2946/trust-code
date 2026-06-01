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


#ifndef Champ_Generique_Predefini_included
#define Champ_Generique_Predefini_included

#include <Champ_Gen_de_Champs_Gen.h>


/*! @brief class Champ_Generique_Predefini
 *
 */

// Field intended to encapsulate a generic field whose expression is predefined
// so that users have shortcuts in terms of syntax
//
//// Data file syntax to follow
//
// "field_name" Predefini { Pb_champ "pb_name" "field_name_to_create" }
//
// "field_name" set by the user will be the name of the constructed generic field (field_)
// "field_name_to_create" type of generic field to create (ex: kinetic_energy)


class Champ_Generique_Predefini : public Champ_Gen_de_Champs_Gen
{

  Declare_instanciable(Champ_Generique_Predefini);

public:

  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void completer(const Postraitement_base& post) override;
  const Champ_Generique_base&  get_source(int i) const override;
  const Noms get_property(const Motcle& query) const override;
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  void nommer(const Nom&) override;
  const Nom& get_nom_post() const override;
  void nommer_source() override;
  Nom construit_expression();

protected:

  Nom type_champ_;         //Type of predefined field to read (ex: energie_cinetique)
  Nom nom_pb_;
  OWN_PTR(Champ_Generique_base) champ_;  //The generic field predefined by type_champ_

};

#endif

