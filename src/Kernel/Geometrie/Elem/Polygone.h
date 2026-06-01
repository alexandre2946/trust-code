/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Polygone_included
#define Polygone_included

#include <Poly_geom_base.h>

/*! @brief Class Polygone represents the geometric element Polygone.
 *
 *  A polyhedron is an element defined by its faces of type Type_Face::polygone_3D
 *
 * @sa Poly_geom_base Elem_geom
 */
template <typename _SIZE_>
class Polygone_32_64  : public Poly_geom_base_32_64<_SIZE_>
{

  Declare_instanciable_32_64(Polygone_32_64);

public :
  // Template classes with different template parameter do not see each other - build_reduced() needs this.
  template <class T> friend class Polygone_32_64;

  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;
  using DoubleVect_t = DoubleVect_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;
  using Domaine_t = Domaine_32_64<_SIZE_>;


  inline int face_sommet(int i, int j) const override;
  inline int nb_som() const override;
  inline int nb_faces(int=0) const override;
  inline int nb_som_face(int=0) const override;
  const Nom& nom_lml() const override;

  int contient(const ArrOfDouble& pos, int_t elem) const override;
  int contient(const SmallArrOfTID_t& soms, int_t elem) const override;
  inline Type_Face type_face(int=0) const override;
  void calculer_volumes(DoubleVect_t& vols) const override;
  void calculer_centres_gravite(DoubleTab_t& xp) const override;

  void calculer_un_centre_gravite(const int_t elem, DoubleVect& xp) const override ;

  int get_tab_faces_sommets_locaux(IntTab& faces_som_local) const override;
  int get_tab_faces_sommets_locaux(IntTab& faces_som_local, int_t elem) const override;
  void affecte_connectivite_numero_global(const ArrOfInt_t& FacesIndex,const ArrOfInt_t& PolygonIndex,IntTab_t& les_elems);

  int get_nb_som_elem_max() const;
  inline int get_nb_face_elem_max() const { return nb_face_elem_max_ ;  }

  _SIZE_ get_somme_nb_faces_elem() const override;

  inline const ArrOfInt_t& getPolygonIndex() const        { return PolygonIndex_; }
  inline const ArrOfInt_t& getElemIndex() const override  { return PolygonIndex_; }
  void rebuild_index();
  void build_reduced(OWN_PTR(Elem_geom_base_32_64<int>)& type_elem, const ArrOfInt_t& elems_sous_part) const override;
  void compute_virtual_index() override;

protected:
  // Members herited from top classes:
  using Objet_U::dimension;
  using Elem_geom_base_32_64<_SIZE_>::mon_dom;
  // FacesIndex_[f] first vertex index of face f in Nodes
  using Poly_geom_base_32_64<_SIZE_>::FacesIndex_;
  using Poly_geom_base_32_64<_SIZE_>::nb_som_elem_max_;
  using Poly_geom_base_32_64<_SIZE_>::nb_face_elem_max_;

  /*! PolygonIndex_[ele] gives the index of the first face of element 'ele'
   *  See FaceIndex_ in base class.
   */
  ArrOfInt_t PolygonIndex_;
};


/*! @brief Returns the index of the j-th vertex of the i-th face of the element.
 *
 * @param (int i) a face index
 * @param (int j) a vertex index
 * @return (int) the index of the j-th vertex of the i-th face
 */
template <typename _SIZE_>
inline int Polygone_32_64<_SIZE_>::face_sommet(int face, int sommet) const
{
  BLOQUE;
  return -1;
}


/*! @brief Returns the number of vertices of a Polygone element.
 *
 * @return (int) the maximum number of vertices of a polyhedron
 */
template <typename _SIZE_>
inline int Polygone_32_64<_SIZE_>::nb_som() const
{
  return get_nb_som_elem_max();
}


/*! @brief Returns the number of faces of the specified type that the geometric element has.
 *
 *     A Polygone has 1 type of face: polygon_3D
 *
 * @param (int i) the face type
 */
template <typename _SIZE_>
inline int Polygone_32_64<_SIZE_>::nb_faces(int i) const
{
  assert(i==0);
  switch(i)
    {
    case 0:
      return get_nb_face_elem_max();
    default :
      Cerr << "Error, a polyhedron has 1 type of faces and not " << i << finl;
      Process::exit();
      return -1;
    }
}


/*! @brief Returns the maximum number of vertices of the faces of the specified type.
 *
 * @param (int i) the face type
 * @return (int) the number of vertices of faces of type i
 */
template <typename _SIZE_>
inline int Polygone_32_64<_SIZE_>::nb_som_face(int i) const
{
  assert(i==0);
  return 2;
}

/*! @brief Returns the i-th face type.
 *
 * A polyhedron has 2 face types: quadrangle and triangle
 *
 * @param (int i) the rank of the face type to return
 * @return (Type_Face) the type of face i
 */
template <typename _SIZE_>
inline Type_Face Polygone_32_64<_SIZE_>::type_face(int i) const
{
  assert(i<=0);
  switch(i)
    {
    case 0:
      return Type_Face::segment_2D;
    default :
      Cerr << "Error, a polyhedron has 1 type of faces and not " << i << finl;
      Process::exit();
      return Type_Face::quadrangle_3D;
    }
}


using Polygone = Polygone_32_64<int>;
using Polygone_64 = Polygone_32_64<trustIdType>;

#endif
