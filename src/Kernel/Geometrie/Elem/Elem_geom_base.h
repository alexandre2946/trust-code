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

#ifndef Elem_geom_base_included
#define Elem_geom_base_included

#include <TRUSTTabs_forward.h>
#include <Faces.h>

template<class _SIZE_> class Domaine_32_64;

/*! @brief Class Elem_geom_base This class is the base class for the definition of elements
 *
 *     geometric constituting a mesh (i.e. a geometric Domain)
 *     A geometric element has an associated Domain to which it provides
 *     the basic routines for manipulating its elements.
 *     (A domain has only one type of geometric element)
 *
 * @sa Hexaedre Prisme Rectangle Segment Tetraedre Triangle, Domaine, Abstract class, Abstract methods:, int face_sommet(int i, int j) const, int nb_som() const, int nb_faces(int=0) const, int nb_som_face(int=0) const, int est_structure() const, const Nom& nom_lml() const, int contient(const ArrOfDouble&, int ) const, int contient(const ArrOfInt&, int ) const, Type_Face type_face(int=0) const, void calculer_volumes(DoubleVect& ) const, void calculer_normales(const IntTab& , DoubleTab& ) const
 */
template <typename _SIZE_>
class Elem_geom_base_32_64 : public Objet_U
{

  Declare_base_32_64(Elem_geom_base_32_64);

public:

  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using IntVect_t = IntVect_T<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;
  using ArrOfDouble_t= ArrOfDouble_T<_SIZE_>;
  using DoubleVect_t = DoubleVect_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;
  using Domaine_t = Domaine_32_64<_SIZE_>;
  using Faces_t = Faces_32_64<_SIZE_>;


  virtual void creer_faces_elem(Faces_t& ,int_t ,Type_Face ) const;
  inline void creer_faces_elem(Faces_t& ,int_t ) const;
  /// Returns the number of the j-th vertex of the i-th face of the element
  virtual int face_sommet(int i, int j) const=0;
  inline void associer_domaine(const Domaine_32_64<int_t>& dom) { mon_dom = dom; }
  /// Nb of vertices for the element
  virtual int nb_som() const=0;
  /// Nb of faces for the element
  virtual int nb_faces(int=0) const=0;
  /// Nb of vertices for one face of the element
  virtual int nb_som_face(int=0) const=0;
  virtual bool est_structure() const=0;
  virtual const Nom& nom_lml() const=0;
  /// DOes the element 'elem' contains the point 'pos'
  virtual int contient(const ArrOfDouble& pos, int_t elem) const=0;
  /// Returns 1 if the vertices specified by the parameter "pos" are the vertices of the element "element" of the associated domain
  virtual int contient(const SmallArrOfTID_t& soms, int_t elem) const=0;
  /// Nb of face types of the elemnt (for example 2 for a prism)
  virtual int nb_type_face() const;
  virtual int num_face(int face, Type_Face& type) const;
  /// Type of the face of the element - face_typ < nb_type_face()
  virtual Type_Face type_face(int face_typ=0) const=0;
  /// Compute all centers of mass of all elements in the domain
  virtual void calculer_centres_gravite(DoubleTab_t& ) const ;
  virtual void reordonner() { }
  /// Compute vols of all elements in the domain
  virtual void calculer_volumes(DoubleVect_t& vols) const =0;
  virtual void calculer_normales(const IntTab_t& faces_sommets , DoubleTab_t& face_normales) const;
  virtual int get_tab_faces_sommets_locaux(IntTab& faces_som_local) const;
  virtual void get_tab_aretes_sommets_locaux(IntTab& aretes_som_local) const;

protected:
  OBS_PTR(Domaine_t) mon_dom;
};

/*! @brief Create the faces of the element of the specified domain.
 *
 * @param (Faces_t& faces) the faces of the elements to create
 * @param (int elem) the number of the element of the domain whose faces we want to create
 */
template <typename _SIZE_>
inline void Elem_geom_base_32_64<_SIZE_>::creer_faces_elem(Faces_t& faces, _SIZE_ elem) const
{
  assert(nb_type_face() == 1);
  creer_faces_elem(faces, elem, type_face());
}

using Elem_geom_base = Elem_geom_base_32_64<int>;
using Elem_geom_base_64 = Elem_geom_base_32_64<trustIdType>;

#endif
