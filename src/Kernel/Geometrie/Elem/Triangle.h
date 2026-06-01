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

#ifndef Triangle_included
#define Triangle_included

#include <Elem_geom_base.h>
#include <TRUSTTab.h>

/*! @brief Class Triangle: represents the geometric element Triangle.
 *
 *     A triangle has 3 faces, 3 vertices, and a single face type with
 *     2 vertices per face.
 *
 * @sa Elem_geom_base Elem_geom
 */
template <typename _SIZE_>
class Triangle_32_64 : public Elem_geom_base_32_64<_SIZE_>
{

  Declare_instanciable_32_64(Triangle_32_64);

public :

  using int_t = _SIZE_;
  using IntTab_t = IntTab_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;
  using DoubleVect_t = DoubleVect_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;
  using Domaine_t = Domaine_32_64<_SIZE_>;


  inline int face_sommet(int i, int j) const override;
  inline int face_sommet0(int i) const;
  inline int face_sommet1(int i) const;
  inline int face_sommet2(int i) const;
  inline int nb_som() const override { return 3; }
  inline int nb_faces(int = 0) const override;
  inline int nb_som_face(int = 0) const override;
  inline bool est_structure() const override { return false; }
  const Nom& nom_lml() const override;
  int contient(const ArrOfDouble& pos, int_t elem) const override;
  int contient(const SmallArrOfTID_t& soms, int_t elem) const override;
  inline Type_Face type_face(int=0) const override;
  void calculer_volumes(DoubleVect_t& vols) const override;
  void calculer_normales(const IntTab_t& faces_sommets , DoubleTab_t& face_normales) const override;
  int get_tab_faces_sommets_locaux(IntTab& faces_som_local) const override;

protected:
  // Members herited from top classes:
  using Objet_U::dimension;
  using Elem_geom_base_32_64<_SIZE_>::mon_dom;

};

/*! Computes the area of a triangle in 2D or 3D.
 * Parameters passed:
 * pos: array containing the coordinates of the 3 vertices of the triangle
 */
inline double aire_triangle(const DoubleTab& pos)
{
  double x0 = pos(0,0);
  double y0 = pos(0,1);
  double x1 = pos(1,0);
  double y1 = pos(1,1);
  double x2 = pos(2,0);
  double y2 = pos(2,1);
  double dz = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);
  if (Objet_U::dimension==2)
    return 0.5 * std::fabs(dz);
  else
    {
      double z0 = pos(0,2);
      double z1 = pos(1,2);
      double z2 = pos(2,2);
      double dx = (y1-y0)*(z2-z0) - (y2-y0)*(z1-z0);
      double dy = (z1-z0)*(x2-x0) - (x1-x0)*(z2-z0);
      return 0.5 * sqrt(dx*dx+dy*dy+dz*dz);
    }
}

/*! @brief Returns the index of the j-th vertex of the i-th face of the element.
 *
 * @param i a face index
 * @param j a vertex index
 * @return the index of the j-th vertex of the i-th face
 */
template <typename _SIZE_>
inline int Triangle_32_64<_SIZE_>::face_sommet(int i, int j) const
{
  assert(i<3);
  switch(i)
    {
    case 0:
      return face_sommet0(j);
    case 1:
      return face_sommet1(j);
    case 2:
      return face_sommet2(j);
    default :
      return -1;
    }
}

/*! @brief Returns the number of faces of the specified type that the geometric element has.
 *
 * @param i the face type
 * @return the number of faces of type i
 */
template <typename _SIZE_>
inline int Triangle_32_64<_SIZE_>::nb_faces(int i) const
{
  assert(i==0);
  return 3;
}


/*! @brief Returns the number of vertices of faces of the specified type.
 *
 * @param i the face type
 * @return the number of vertices of faces of type i
 */
template <typename _SIZE_>
inline int Triangle_32_64<_SIZE_>::nb_som_face(int i) const
{
  assert(i==0);
  return 2;
}


/*! @brief Returns the index of the i-th vertex of face 0.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 0
 */
template <typename _SIZE_>
inline int Triangle_32_64<_SIZE_>::face_sommet0(int i) const
{
  assert(i>=0);
  assert(i<2);
  return (i+1);
}


/*! @brief Returns the index of the i-th vertex of face 1.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 1
 */
template <typename _SIZE_>
inline int Triangle_32_64<_SIZE_>::face_sommet1(int i) const
{
  assert(i>=0);
  assert(i<2);
  return ((i+2)%3);
}


/*! @brief Returns the index of the i-th vertex of face 2.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 2
 */
template <typename _SIZE_>
inline int Triangle_32_64<_SIZE_>::face_sommet2(int i) const
{
  assert(i>=0);
  assert(i<2);
  return i;
}


/*! @brief Returns the i-th face type.
 *
 * A triangle has only one face type.
 *
 * @param i the rank of the face type to return
 * @return a face type
 */
template <typename _SIZE_>
inline Type_Face Triangle_32_64<_SIZE_>::type_face(int i) const
{
  assert(i==0);
  return Type_Face::segment_2D;
}


using Triangle = Triangle_32_64<int>;
using Triangle_64 = Triangle_32_64<trustIdType>;

#endif
