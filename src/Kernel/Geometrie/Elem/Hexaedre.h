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

#ifndef Hexaedre_included
#define Hexaedre_included

#include <Elem_geom_base.h>

/*! @brief Class Hexaedre: represents a geometric element with 6 faces, 8 vertices, and
 *
 *     a single face type (Type_Face::quadrangle_3D) with 4 vertices.
 *
 * @sa Hexaedre_axi Elem_geom_base Elem_geom, A hexahedron is a 3D element
 */
template <typename _SIZE_>
class Hexaedre_32_64  : public Elem_geom_base_32_64<_SIZE_>
{

  Declare_instanciable_32_64(Hexaedre_32_64);

public :
  using int_t = _SIZE_;
  using IntTab_t = IntTab_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;
  using DoubleVect_t = DoubleVect_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;
  using Domaine_t = Domaine_32_64<_SIZE_>;


  int face_sommet(int i, int j) const override;
  inline int face_sommet0(int i) const;
  inline int face_sommet1(int i) const;
  inline int face_sommet2(int i) const;
  inline int face_sommet3(int i) const;
  inline int face_sommet4(int i) const;
  inline int face_sommet5(int i) const;
  inline int nb_som() const override { return 8; }
  inline int nb_faces(int=0) const override;
  inline int nb_som_face(int=0) const override;
  inline bool est_structure() const override { return true; };
  const Nom& nom_lml() const override;
  int contient(const ArrOfDouble& pos, int_t elem) const override;
  int contient(const SmallArrOfTID_t& soms, int_t elem) const override;
  inline Type_Face type_face(int=0) const override;
  void reordonner() override ;
  int reordonner_elem();
  void calculer_volumes(DoubleVect_t& vols) const override;
  void calculer_normales(const IntTab_t& faces_sommets , DoubleTab_t& face_normales) const override;
  int get_tab_faces_sommets_locaux(IntTab& faces_som_local) const override;
};

/*! @brief Returns the number of faces of the specified type that the geometric element has.
 *
 * @param i the face type
 * @return the number of faces of type i
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::nb_faces(int i) const
{
  assert(i==0);
  return 6;
}


/*! @brief Returns the number of vertices of the faces of the specified type.
 *
 * @param i the face type
 * @return the number of vertices of faces of type i
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::nb_som_face(int i) const
{
  assert(i==0);
  return 4;
}


/*! @brief Returns the index of the i-th vertex of face 0.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 0
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::face_sommet0(int i) const
{
  //face_sommet0(0)=0;
  //face_sommet0(1)=2;
  //face_sommet0(2)=4;
  //face_sommet0(3)=6;
  assert(i>=0);
  assert(i<4);
  return 2*i;
}


/*! @brief Returns the index of the i-th vertex of face 1.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 1
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::face_sommet1(int i) const
{
  //face_sommet1(0)=0;
  //face_sommet1(1)=1;
  //face_sommet1(2)=4;
  //face_sommet1(3)=5;
  assert(i>=0);
  assert(i<4);
  return 2*(i/2)+i;
}


/*! @brief Returns the index of the i-th vertex of face 2.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 2
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::face_sommet2(int i) const
{
  //face_sommet2(0)=0;
  //face_sommet2(1)=1;
  //face_sommet2(2)=2;
  //face_sommet2(3)=3;
  assert(i>=0);
  assert(i<4);
  return i;
}


/*! @brief Returns the index of the i-th vertex of face 3.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 3
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::face_sommet3(int i) const
{
  //face_sommet3(0)=1;
  //face_sommet3(1)=3;
  //face_sommet3(2)=5;
  //face_sommet3(3)=7;
  assert(i>=0);
  assert(i<4);
  return 2*i+1;
}


/*! @brief Returns the index of the i-th vertex of face 4.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 4
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::face_sommet4(int i) const
{
  //face_sommet4(0)=2;
  //face_sommet4(1)=3;
  //face_sommet4(2)=6;
  //face_sommet4(3)=7;
  assert(i>=0);
  assert(i<4);
  return 2*(i/2)+i+2;
}


/*! @brief Returns the index of the i-th vertex of face 5.
 *
 * @param i the index of the vertex to return
 * @return the index of the i-th vertex of face 5
 */
template<typename _SIZE_>
inline int Hexaedre_32_64<_SIZE_>::face_sommet5(int i) const
{
  //face_sommet4(0)=4;
  //face_sommet4(1)=5;
  //face_sommet4(2)=6;
  //face_sommet4(3)=7;
  assert(i>=0);
  assert(i<4);
  return i+4;
}


/*! @brief Returns the i-th face type.
 *
 * A hexahedron has only one face type.
 *
 * @param i the rank of the face type to return
 * @return a face type
 */
template<typename _SIZE_>
inline Type_Face Hexaedre_32_64<_SIZE_>::type_face(int i) const
{
  assert(i==0);
  return Type_Face::quadrangle_3D;
}

using Hexaedre = Hexaedre_32_64<int>;
using Hexaedre_64 = Hexaedre_32_64<trustIdType>;

#endif
