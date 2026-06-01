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

#ifndef Domaine_EF_included
#define Domaine_EF_included

#include <Elem_EF_base.h>
#include <TRUST_Deriv.h>
#include <Domaine_VF.h>

class Geometrie;

/*! @brief class Domaine_EF
 *
 * @brief Instantiable class derived from Domaine_VF.
 *  	This class holds the geometric information required by the
 *  	Finite Element (FE) method (Crouzeix-Raviart element).
 *  	The class stores a number of pieces of information about faces.
 *  	The set of faces also includes boundary faces and
 *       joint faces. Faces are split into 2 categories:
 *            - non-standard faces: on a joint, on a boundary, or internal
 *              but belonging to a boundary element
 *            - standard faces: internal faces not belonging to a boundary element
 *       This distinction corresponds to the treatment of boundary conditions:
 *       standard faces do not "see" the boundary conditions.
 *       The full set of faces is numbered as follows:
 *            - faces on a Domaine_joint appear first
 *     	       (in the order of the les_joints vector)
 *    	     - faces on a Domaine_bord appear next
 * 	       (in the order of the les_bords vector)
 *   	     - non-standard internal faces appear next
 *            - standard internal faces appear last
 *       All non-standard faces requiring special treatment are thus grouped first.
 *       Two element types are distinguished:
 *            - non-standard elements: they have at least one boundary face
 *            - standard elements: they have no boundary face
 *       Standard (resp. non-standard) elements are not stored consecutively
 *       in the Domaine object. The rang_elem_non_std array is used to
 *       selectively access one or the other type of element.
 *
 */

class Domaine_EF: public Domaine_VF
{
  Declare_instanciable(Domaine_EF);
public:
  void typer_elem(Domaine& domaine_geom) override;
  void discretiser() override;
  void swap(int, int, int);
  void modifier_pour_Cl(const Conds_lim&) override;

  inline const Elem_EF_base& type_elem() const { return type_elem_.valeur(); }
  inline int nb_elem_Cl() const { return nb_elem() - nb_elem_std_; }
  inline int nb_faces_joint() const { return 0; /*    return nb_faces_joint_;    A FAIRE */ }
  inline int nb_faces_std() const { return nb_faces_std_; }
  inline int nb_elem_std() const { return nb_elem_std_; }
  inline double carre_pas_du_maillage() const { return h_carre; }
  inline double carre_pas_maille(int i) const { return h_carre_(i); }
  inline IntVect& rang_elem_non_std() { return rang_elem_non_std_; }

  inline const IntVect& rang_elem_non_std() const { return rang_elem_non_std_; }
  void calculer_volumes_entrelaces();
  void calculer_volumes_sommets(const Domaine_Cl_dis_base& zcl);
  virtual void calculer_IPhi(const Domaine_Cl_dis_base& zcl);
  virtual void calculer_Bij(DoubleTab& bij_);
  virtual void calculer_Bij_gen(DoubleTab& bij);
  void calculer_Bij() { calculer_Bij(Bij_); }

  //  inline const DoubleVect& volumes_sommets() const { return volumes_sommets_; }
  inline const DoubleVect& volumes_thilde() const { return volumes_thilde_; }
  inline const DoubleVect& volumes_sommets_thilde() const { return volumes_sommets_thilde_; }
  inline const DoubleVect& porosite_sommet() const { return porosite_sommets_; }
  inline DoubleVect& set_porosite_sommet() { return porosite_sommets_; }

  void calculer_h_carre();
  void calculer_porosites_sommets();
  inline const DoubleTab& Bij() const { return Bij_; }
  inline const DoubleTab& Bij_thilde() const { return Bij_thilde_; }
  inline const DoubleTab& IPhi() const { return IPhi_; }
  inline const DoubleTab& IPhi_thilde() const { return IPhi_thilde_; }

  virtual void verifie_compatibilite_domaine();

protected:
  DoubleTab IPhi_, IPhi_thilde_;

private:
  DoubleVect porosite_sommets_, volumes_sommets_thilde_, volumes_thilde_;
  //  OWN_PTR(Champ_Don_base) champ_porosite_sommets_,champ_porosite_lu_;
  DoubleTab Bij_, Bij_thilde_;                         // storage of the Bije matrices

  double h_carre = 1.e30;			 // square of the mesh step
  DoubleVect h_carre_;			// square of the mesh step per cell
  OWN_PTR(Elem_EF_base) type_elem_;                  // type of the discretisation element

  Sortie& ecrit(Sortie& os) const;
  IntVect orientation_;
};

#endif /* Domaine_EF_included */
