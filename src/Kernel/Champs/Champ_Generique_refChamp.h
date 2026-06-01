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

#ifndef Champ_Generique_refChamp_included
#define Champ_Generique_refChamp_included

#include <Postraitement.h>

/*! @brief Special field class that encapsulates a reference to a volume field of TRUST of type Champ_base.
 *
 *   The class does not manage the memory of Champ_base.
 *
 */

// Syntax to follow in the data file
//
// "field_name" refChamp { Pb_champ "pb_name" "discrete_field_name" }
// "field_name" set by the user will be the name of the generic field
// "pb_name" name of the problem to which the discrete field belongs
// "discrete_field_name" name of the target field or one of its components

class Champ_Generique_refChamp : public Champ_Generique_base
{
  Declare_instanciable(Champ_Generique_refChamp);
public:

  void initialize(const Champ_base& champ);
  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int            get_dimension() const override;
  void              get_property_names(Motcles& list) const override;
  const Noms        get_property(const Motcle& query) const override;
  Entity            get_localisation(const int index = -1) const override;

  const DoubleTab&  get_ref_values() const override;
  void              get_copy_values(DoubleTab&) const override;
  void              get_xyz_values(const DoubleTab& coords, DoubleTab& values, ArrOfBit& validity_flag) const override;

  const Domaine_Cl_dis_base&  get_ref_zcl_dis_base() const override;

  const DoubleTab&  get_ref_coordinates() const override;
  void              get_copy_coordinates(DoubleTab&) const override;
  const IntTab&     get_ref_connectivity(Entity index1, Entity index2) const override;
  void              get_copy_connectivity(Entity index1, Entity index2, IntTab&) const override;

  const Probleme_base& get_ref_pb_base() const override;
  void              reset() override;
  void completer(const Postraitement_base& post) override;
  void              mettre_a_jour(double temps) override;

  //Particular get_champ() that does not use a storage space
  //Updates the discrete field if it is a computed field of the problem
  //and returns the reference
  const Champ_base& get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base& get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  virtual const Champ_base& get_ref_champ_base() const; //returns the reference to the encapsulated field

  virtual void set_ref_champ(const Champ_base&);

  double get_time() const override;
  const  Motcle            get_directive_pour_discr() const override;
  void set_nom_champ(const Motcle&);
  void nommer_source(const Postraitement_base& post);
  int get_info_type_post() const override;

  //The compo_ attribute of Champ_Generique_refChamp is only filled for CGNS
  //instances created by macro, in order to reproduce component names for face fields
  const Noms& fixer_noms_compo(const Noms& noms) override;
  const Noms& fixer_noms_synonyms(const Noms& noms) override;

protected:
  Noms compo_,syno_;

  OBS_PTR(Champ_base) ref_champ_;
  mutable OWN_PTR(Champ_base) ptr_champ_; /* XXX Elie Saikali : not sure what to do here */

  //temporary - check usefulness
  Motcle localisation_;

  Nom nom_champ_;

};
#endif
