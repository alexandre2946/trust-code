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

#ifndef Discretisation_base_included
#define Discretisation_base_included

#include <Domaine_forward.h>
#include <Champ_base.h> // For Nature_du_champ

#include <TRUST_Ref.h>
#include <Reorder_Mesh.h>

class Champ_Fonc_Tabule;
class Schema_Temps_base;
class Domaine_dis_base;
class Champ_Fonc_base;
class Champ_Inc_base;
class Champ_Don_base;
class Probleme_base;
class Equation_base;
class Champ_base;
class Motcle;

/*! @brief class Discretisation_base This class represents a spatial discretization scheme, which
 *
 *      will be associated with a problem. Discretisation_base is the
 *      abstract class that is the basis of the hierarchy of
 *      spatial discretizations.
 *
 * @sa Probleme_base, Abstract class from which all spatial discretizations must derive., Abstract method:, void domaine_Cl_dis(Domaine_dis_base& , Domaine_Cl_dis_base& ) const
 */
class Discretisation_base : public Objet_U
{
  Declare_base(Discretisation_base);

public :

  // MODIF ELI LAUCOIN (7/06/2007):
  // Adding a pure virtual method that specifies the mode of assembling contributions to the residual of operators and sources.
  //
  // By default, in TRUST, we go through contribuer_a_avec, then through ajouter. In the porous modules MTMS (and EF?), we rather go through contribuer_a_avec and
  // contribuer_au_second_membre, because we do not use a priori a convection scheme (like Quick) for which the Jacobian is not exact.
  //
  // The default value returned is VIA_AJOUTER, so as not to disturb the normal behavior of Trio-U classes outside the kernel.
  //
  enum type_calcul_du_residu { VIA_CONTRIBUER_AU_SECOND_MEMBRE = 0, VIA_AJOUTER = 1 };
  inline virtual type_calcul_du_residu codage_du_calcul_du_residu() const { return VIA_AJOUTER; }

  void associer_domaine(const Domaine& dom);

  virtual Nom domaine_cl_dis_type() const = 0;

  virtual void discretiser_variables() const;
  virtual Domaine_dis_base& discretiser() const;

  // Creation of scalar or vector fields (essentially call to the general method, do not overload these methods, they are only here for convenience)
  void discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, const Nom& nom, const Nom& unite, int nb_comp, int nb_pas_dt, double temps, OWN_PTR(Champ_Inc_base)& champ,
                         const Nom& sous_type=NOM_VIDE) const;
  void discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, const Nom& nom, const Nom& unite, int nb_comp, double temps, OWN_PTR(Champ_Fonc_base)& champ) const;
  void discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, const Nom& nom, const Nom& unite, int nb_comp, double temps, OWN_PTR(Champ_Don_base)& champ) const;

  // Creation of general fields (possibly multiscalar):
  // * These methods must be overloaded.
  // * Each method includes the demande_description keyword, which causes the display of all understood directives (and calls the ancestor with the same keyword).
  // * If scalar or vector field, the first name and first unit are used
  virtual void discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& nom, const Noms& unite, int nb_comp, int nb_pas_dt, double temps,
                                 OWN_PTR(Champ_Inc_base)& champ, const Nom& sous_type=NOM_VIDE) const;
  virtual void discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& nom, const Noms& unite, int nb_comp, double temps, OWN_PTR(Champ_Fonc_base)& champ) const;
  virtual void discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& nom, const Noms& unite, int nb_comp, double temps, OWN_PTR(Champ_Don_base)& champ) const;

  void nommer_completer_champ_physique(const Domaine_dis_base& domaine_vdf, const Nom& nom_champ, const Nom& unite, Champ_base& champ, const Probleme_base& pbi) const;
  int verifie_sous_type(Nom& type, const Nom& sous_type, const Motcle& directive) const;

  void volume_maille(const Schema_Temps_base& sch, const Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const;
  void mesh_numbering(const Schema_Temps_base& sch, const Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const;
  virtual void residu(const Domaine_dis_base&, const Champ_Inc_base&, OWN_PTR(Champ_Fonc_base)&) const;

  static void creer_champ(OWN_PTR(Champ_Inc_base)& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, int nb_pas_dt, double temps,
                          const Nom& directive=NOM_VIDE, const Nom& nom_discretisation=NOM_VIDE);
  static void creer_champ(OWN_PTR(Champ_Fonc_base)& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, double temps, const Nom& directive = NOM_VIDE,
                          const Nom& nom_discretisation=NOM_VIDE);
  static void creer_champ(OWN_PTR(Champ_Don_base)& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, double temps, const Nom& directive = NOM_VIDE,
                          const Nom& nom_discretisation=NOM_VIDE);

  virtual Nom get_name_of_type_for(const Nom& class_operateur, const Nom& type_operteur,const Equation_base& eqn, const OBS_PTR(Champ_base)& champ_supp =OBS_PTR(Champ_base)()) const;

  // usefull methods to detect discretization
  virtual bool is_ef() const { return false; }
  virtual bool is_dg() const { return false; }
  virtual bool is_vdf() const { return false; }
  virtual bool is_vef() const { return false; }
  virtual bool is_PolyMAC_CDO() const { return false; }
  virtual bool is_PolyMAC_MPFA() const { return false; }
  virtual bool is_PolyMAC_HFV() const { return false; }
  virtual bool is_poly_family() const { return false; }
  virtual bool is_coloc() const { return false; }

  const Reorder_Mesh& get_reorder() const { return reorder_; }

protected:
  static const Motcle DEMANDE_DESCRIPTION;
  static const Nom NOM_VIDE;
  OBS_PTR(Domaine) le_domaine_;
  Reorder_Mesh reorder_;   ///< Helper object to renumber entities (nodes, elems, faces) if requested

private:
  void test_demande_description(const Motcle& , const Nom&) const;
  static void champ_fixer_membres_communs(Champ_base& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, double temps);

  virtual void modifier_champ_tabule(const Domaine_dis_base& domaine_dis,Champ_Fonc_Tabule& ch_tab,const VECT(OBS_PTR(Champ_base))& ch_inc) const ;
};

using Discretisation = TRUST_Deriv<Discretisation_base>;

#endif /* Discretisation_base_included */
