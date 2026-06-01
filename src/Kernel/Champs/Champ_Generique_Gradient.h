/****************************************************************************
* Copyright (c) 2024, CEA
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


#ifndef Champ_Generique_Gradient_included
#define Champ_Generique_Gradient_included

#include <Champ_Generique_Operateur_base.h>
#include <Operateur_Grad.h>

/*! @brief class Champ_Generique_Gradient OWN_PTR(Champ_base) intended to post-process the gradient of a generic field
 *
 *  The class carries a statistical operator "gradient"
 *
 */

//// Data file syntax to follow
//
// "field_name" Gradient {
//                source "generic_field_type" { ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name" } }
//               }
// "field_name" set by the user will be the name of the generic field
// "generic_field_type" type of a generic field
// This type of field implies that the source field has boundary conditions
// Its application is restricted to certain discrete fields (pressure VDF and VEF or temperature in VEF)

class Champ_Generique_Gradient : public Champ_Generique_Operateur_base
{

  Declare_instanciable(Champ_Generique_Gradient);

public:

  const Noms get_property(const Motcle& query) const override;
  Entity  get_localisation(const int index = -1) const override;
  const   Motcle             get_directive_pour_discr() const override;
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;

  inline const Operateur_base& Operateur() const override;
  inline Operateur_base& Operateur() override;
  void completer(const Postraitement_base& post) override;
  void nommer_source() override;

protected:

  Operateur_Grad Op_Grad_;

};

inline const Operateur_base& Champ_Generique_Gradient::Operateur() const
{
  return Op_Grad_.valeur();
}

inline Operateur_base& Champ_Generique_Gradient::Operateur()
{
  return Op_Grad_.valeur();
}

#endif
