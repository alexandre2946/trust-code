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

#ifndef Champ_Generique_Extraction_included
#define Champ_Generique_Extraction_included

#include <Champ_Gen_de_Champs_Gen.h>


/*! @brief A generic field that performs the extraction of a field on a boundary
 *
 */

//// Data file syntax to follow
//
// "field_name" Extraction { domain "domain_name" boundary_name "boundary_name" [ method ] "method_type"
//                source "generic_field_type" { ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name" } }
//               }
// "field_name" set by the user will be the name of the generic field
// "domain_name" name of the domain to which the boundary belongs
// "boundary_name" name of the boundary on which we want to perform the extraction
// "method_type" type of method to perform the extraction
//                  ("trace" default method or "boundary_field" to extract the_boundary_field)
// "generic_field_type" type of a generic field

class Champ_Generique_Extraction : public Champ_Gen_de_Champs_Gen
{
  Declare_instanciable_sans_constructeur(Champ_Generique_Extraction);

public:
  Champ_Generique_Extraction();
  void set_param(Param& param) const override;
  Entity  get_localisation(const int index = -1) const override;
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Noms        get_property(const Motcle& query) const override;
  void nommer_source() override;
  void completer(const Postraitement_base& post) override;
  const Domaine& get_ref_domain() const override;
  void get_copy_domain(Domaine&) const override;
  const Domaine_dis_base& get_ref_domaine_dis_base() const override;
  void discretiser_domaine();
  const  Motcle  get_directive_pour_discr() const override;

protected :
  Nom dom_extrac_;              // Name of the extraction domain
  Nom nom_fr_;                  // Name of the boundary on which the extraction is performed
  Nom methode_;                 // Type of method for extraction ("trace" or "champ_frontiere")
  OBS_PTR(Domaine) domaine_;        // Reference to the extraction domain
  OBS_PTR(Domaine_dis_base) le_dom_dis;  // The discretized domain - real owner is Domaine_dis_cache
};

#endif
