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

#ifndef Champ_Generique_Interpolation_included
#define Champ_Generique_Interpolation_included

#include <Champ_Gen_de_Champs_Gen.h>


/*! @brief A generic field constructed as an interpolation of another generic field (interpolation at vertices or elements).
 *
 *   The interpolation will be performed on a domain that can be
 *   the computational domain or a different domain specified by the user
 *
 */

//// Data file syntax to follow
//
// "field_name" Interpolation { [ domain ] "interp_domain"
//                localisation "loc"
//                source "generic_field_type" { ...source ref_Champ { Pb_champ "pb_name" "discrete_field_name" } }
//               }
// "field_name" set by the user will be the name of the generic field
// "interp_domain" name of the interpolation domain in case it differs from the computational domain (optional)
// "loc" designates the interpolation localisation "elem", "elem_dg", "faces" or "som"
// "generic_field_type" type of a generic field

class Champ_Generique_Interpolation : public Champ_Gen_de_Champs_Gen
{
  Declare_instanciable(Champ_Generique_Interpolation);
public:
  void reset() override;
  void set_param(Param& param) const override;
  //int lire_motcle_non_standard(const Motcle&, Entree&);
  virtual int     set_localisation(const Motcle& localisation, int exit_on_error = 1);
  virtual int     set_methode(const Motcle& methode, int exit_on_error = 1);
  virtual int     set_domaine(const Nom& nom_domaine, int exit_on_error = 1);
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  virtual const Champ_base&  get_champ_with_calculer_champ_post() const;

  const DoubleTab&  get_ref_values() const override;
  void              get_copy_values(DoubleTab&) const override;
  //virtual void  get_xyz_values(const DoubleTab & coords, DoubleTab & values, ArrOfBit & validity_flag) const;

  Entity  get_localisation(const int index = -1) const override;
  const Domaine&      get_ref_domain() const override;
  void get_copy_domain(Domaine&) const override;

  const DoubleTab&  get_ref_coordinates() const override;
  void              get_copy_coordinates(DoubleTab&) const override;
  //virtual const IntTab&   get_ref_connectivity(Entity index1, Entity index2) const;
  //virtual void            get_copy_connectivity(Entity index1, Entity index2, IntTab &) const;

  const Noms        get_property(const Motcle& query) const override;
  void nommer_source() override;
  void completer(const Postraitement_base& post) override;
  const Domaine_dis_base& get_ref_domaine_dis_base() const override;
  const   Motcle            get_directive_pour_discr() const override;
  void discretiser_domaine();

  const Noms& fixer_noms_compo(const Noms& noms) override;
  const Noms& fixer_noms_synonyms(const Noms& noms) override;
  //The compo_ attribute of Champ_Generique_Interpolation is only filled for Champ_Generique_Interpolation
  //instances created by macro, in order to reproduce component names in lml files
  Noms compo_,syno_;

private:
  Motcle            localisation_;                 // interpolation localisation: elem, som
  Motcle            methode_;                      // calculer_champ_post, etc...
  Nom               nom_domaine_lu_;               // Name of the domain read
  OBS_PTR(Domaine)      domaine_;                      // domain on which we want to interpolate the field (native domain if reference is null)
  OBS_PTR(Domaine_dis_base)  le_dom_dis;                    // filled if interpolation domain differs from native domain. A REF since Domaine_dis_cache owns the memory
  // ex: Sonde uses valeur_aux...() which requires a discretized domain
  int optimisation_sous_maillage_,optimisation_demande_;
  ArrOfInt renumerotation_maillage_;
  mutable OWN_PTR(Champ_Fonc_base) espace_stockage_;
  mutable OWN_PTR(Champ_base) espace_stockage_source_;
};

#endif
