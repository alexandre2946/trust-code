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

#ifndef Point_included
#define Point_included

#include <Elem_geom_base.h>

/*! @brief Class Point: represents the geometric element point.
 *
 *     A point has 1 face and 1 vertex.
 *
 * @sa Elem_geom_base Elem_geom, .CONTRAINTES, .INVARIANTS, .HTML, .EPS
 */
template <typename _SIZE_>
class Point_32_64 : public Elem_geom_base_32_64<_SIZE_>
{

  Declare_instanciable_32_64(Point_32_64);

public :

  using int_t = _SIZE_;
  using IntTab_t = IntTab_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;
  using DoubleVect_t = DoubleVect_T<_SIZE_>;
  using Domaine_t = Domaine_32_64<_SIZE_>;

  inline int face_sommet(int i, int j) const override;
  inline int face_sommet0(int i) const;

  inline int nb_som() const override { return 1; }
  inline int nb_faces(int=0) const override;
  inline int nb_som_face(int=0) const override;
  inline bool est_structure() const override { return true; }
  const Nom& nom_lml() const override;
  int contient(const ArrOfDouble& pos, int_t elem) const override;
  int contient(const SmallArrOfTID_t& soms, int_t elem) const override;
  inline Type_Face type_face(int=0) const override;
  void reordonner() override { }
  void calculer_volumes(DoubleVect_t& vols) const override { vols = 1.; }
  int get_tab_faces_sommets_locaux(IntTab& faces_som_local) const override;
};


/*! @brief Returns the index of the j-th vertex of the i-th face of the element.
 *
 * @param i a face index
 * @param j a vertex index
 * @return the index of the j-th vertex of the i-th face
 */
template <typename _SIZE_>
inline int Point_32_64<_SIZE_>::face_sommet(int i, int j) const
{
  assert(i<1);
  switch(i)
    {
    case 0:
      return face_sommet0(j);
    default :
      return -1;
    }
}


/*! @brief Returns the index of the i-th vertex of face 0. NOTE: THE CODE DOES NOT DO WHAT ONE WOULD EXPECT.
 *
 * @param i the index of the vertex to return
 * @return always returns 0 (if i=0)
 */
template <typename _SIZE_>
inline int Point_32_64<_SIZE_>::face_sommet0(int i) const
{
  assert(i==0);
  return 0;
}


/*! @brief Returns the number of faces of the specified type that the geometric element has.
 *
 * @param i the face type
 * @return the number of faces of type i
 */
template <typename _SIZE_>
inline int Point_32_64<_SIZE_>::nb_faces(int i) const
{
  assert(i==0);
  return 1;
}


/*! @brief Returns the number of vertices of faces of the specified type.
 *
 *     NOTE: THE CODE DOES NOT DO WHAT ONE WOULD EXPECT.
 *
 * @param i the face type
 * @return the number of vertices of faces of type i
 */
template <typename _SIZE_>
inline int Point_32_64<_SIZE_>::nb_som_face(int i) const
{
  assert(i==0);
  return 1;
}


/*! @brief Returns the i-th face type.
 *
 * A point has only one face type.
 *
 * @param i the rank of the face type to return
 * @return a face type
 */
template <typename _SIZE_>
inline Type_Face Point_32_64<_SIZE_>::type_face(int i) const
{
  assert(i==0);
  return Type_Face::point_1D;
}

using Point = Point_32_64<int>;
using Point_64 = Point_32_64<trustIdType>;


#endif
