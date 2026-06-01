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

#ifndef Faces_builder_included
#define Faces_builder_included

#include <Bords.h>
#include <Bords_Internes.h>
#include <Poly_geom_base.h>
#include <Joint.h>
#include <NettoieNoeuds.h>
#include <Polyedriser.h>
#include <Raccord_base.h>
#include <Raccords.h>
#include <TRUST_Ref.h>
#include <TRUSTTab.h>
#include <Domaine_forward.h>
#include <Static_Int_Lists.h>
#include <NettoieNoeuds.h>
/*! @brief Helper class for building the faces of a domain.
 * (used only to create the arrays of real faces)
 *
 */

template <typename _SIZE_>
class Faces_builder_32_64
{
public:
  using int_t = _SIZE_;
  using Domaine_t = Domaine_32_64<_SIZE_>;
  using Static_Int_Lists_t = Static_Int_Lists_32_64<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;
  using Faces_t = Faces_32_64<_SIZE_>;
  using Groupes_Faces_t = Groupes_Faces_32_64<_SIZE_>;
  using Groupe_Faces_t = Groupe_Faces_32_64<_SIZE_>;
  using Poly_geom_base_t = Poly_geom_base_32_64<_SIZE_>;
  using Bords_t = Bords_32_64<_SIZE_>;
  using Frontiere_t = Frontiere_32_64<_SIZE_>;
  using Raccords_t = Raccords_32_64<_SIZE_>;
  using Bords_Internes_t = Bords_Internes_32_64<_SIZE_>;
  using Joints_t = Joints_32_64<_SIZE_>;
  using Joint_t = Joint_32_64<_SIZE_>;
  using NettoieNoeuds_t =NettoieNoeuds_32_64<_SIZE_>;
  using Faces_builder_t = Faces_builder_32_64<_SIZE_>;

  Faces_builder_32_64();
  void reset();
  void creer_faces_reeles(Domaine_t& domaine, const Static_Int_Lists_t& connect_som_elem, Faces_t& les_faces, IntTab_t& elem_faces);

  // This static method is also used by 64b objects like Raffiner_Simplexes:
  static int chercher_face_element(const IntTab_t& elem_som, const IntTab& faces_element_ref,
                                   const  SmallArrOfTID_t& une_face, const int_t elem);

private:
  static int_t ajouter_une_face(const SmallArrOfTID_t& une_face, const int_t elem0, const int_t elem1, IntTab_t& faces_sommets, IntTab_t& faces_voisins);

  int chercher_face_element(const SmallArrOfTID_t& une_face, const int_t elem) const;

  void check_erreur_faces(const char *message, const ArrOfInt_t& liste_faces) const;
  void creer_faces_frontiere(const int_t nb_voisins_attendus, Frontiere_t& frontiere, IntTab_t& faces_sommets, IntTab_t& faces_voisins, IntTab_t& elem_faces) const;
  void creer_faces_internes(IntTab_t& faces_sommets, IntTab_t& elem_faces, IntTab_t& faces_voisins) const;
  void identification_groupe_faces(Groupe_Faces_t& groupe_int, const IntTab_t& elem_faces) const;

  const IntTab_t& les_elements() const { return *les_elements_ptr_; }
  const Static_Int_Lists_t& connectivite_som_elem() const { return *connectivite_som_elem_ptr_; }
  const IntTab& faces_element_reference(int_t elem_t) const; // returns the faces for an element; for polyhedra this depends on the element index.


  // All the following members are initialized in creer_faces_reeles:

  // Shortcut to the domain elements
  const IntTab_t * les_elements_ptr_;
  // The element-vertex connectivity (for each vertex, list of
  // adjacent elements, including virtual vertices and elements)
  const Static_Int_Lists_t * connectivite_som_elem_ptr_;
  /*! @brief Faces of the reference element (see elem_geom_base::get_tab_faces_sommets_locaux)
   *
   */
  IntTab faces_element_reference_old_;
  int_t is_polyedre_;
  // for check_erreur_faces:
  OBS_PTR(Domaine_t)   ref_domaine_;
  OBS_PTR(IntTab_t) faces_sommets_;
  OBS_PTR(IntTab_t) face_elem_;
};


using Faces_builder = Faces_builder_32_64<int>;
using Faces_builder_64 = Faces_builder_32_64<trustIdType>;

#endif
