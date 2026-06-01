/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Triangle.h>
#include <Domaine.h>

Implemente_instanciable_32_64(Triangle_32_64,"Triangle",Elem_geom_base_32_64<_T_>);

/*! @brief Faces of the reference triangle: 3 faces of two vertices each.
 *
 * Face i is the face opposite to vertex i
 * (see get_tab_faces_sommets_locaux).
 *
 */
static int faces_sommets_triangle[3][2] =
{
  { 1, 2 },
  { 2, 0 },
  { 0, 1 }
};

template <typename _SIZE_>
Sortie& Triangle_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  return s;
}

template <typename _SIZE_>
Entree& Triangle_32_64<_SIZE_>::readOn(Entree& s )
{
  return s;
}


/*! @brief Returns the LML name of a triangle = "PRISM6".
 *
 * @return Always equal to "PRISM6" (or "TRIANGLE_3D" in 3D).
 */
template <typename _SIZE_>
const Nom& Triangle_32_64<_SIZE_>::nom_lml() const
{
  static Nom nom="PRISM6";
  if (dimension==3) nom="TRIANGLE_3D";
  return nom;
}


/*! @brief Returns 1 if element "ielem" of the domain associated with this geometric element contains the point
 *
 * with coordinates specified by "pos". Returns 0 otherwise.
 *
 * @param pos Coordinates of the point to locate.
 * @param ielem Index of the domain element in which to search for the point.
 * @return 1 if the point belongs to element "ielem", 0 otherwise.
 */
template <typename _SIZE_>
int Triangle_32_64<_SIZE_>::contient(const ArrOfDouble& pos, int_t ielem) const
{
  assert(pos.size_array()==2);
  const Domaine_t& dom=this->mon_dom.valeur();
  assert(ielem<dom.nb_elem_tot());
  int_t som0 = dom.sommet_elem(ielem,0);
  int_t som1 = dom.sommet_elem(ielem,1);
  int_t som2 = dom.sommet_elem(ielem,2);
  assert((som0>=0) && (som0<dom.nb_som_tot()));
  assert((som1>=0) && (som1<dom.nb_som_tot()));
  assert((som2>=0) && (som2<dom.nb_som_tot()));
  double prod,p0,p1,p2;

  // First check if the point is one of the triangle vertices.
  // GF: this test is removed to be consistent with Tetraedre::contient and to avoid issues in Champ_implementation_P1::form_function, which does not have this test.
  /*
    if( (est_egal(dom.coord(som0,0),pos(0)) && est_egal(dom.coord(som0,1),pos(1)))
    || (est_egal(dom.coord(som1,0),pos(0)) && est_egal(dom.coord(som1,1),pos(1)))
    || (est_egal(dom.coord(som2,0),pos(0)) && est_egal(dom.coord(som2,1),pos(1))) )
    return 1;

  */
  // Note: vertices are stored in arbitrary order.
  // Determine the orientation (counter-clockwise or clockwise) for the vertex numbering:
  // Compute prod = 01 cross 02 along z
  // prod > 0 : counter-clockwise
  // prod < 0 : clockwise
  prod = (dom.coord(som1,0)-dom.coord(som0,0))*(dom.coord(som2,1)-dom.coord(som0,1))
         - (dom.coord(som1,1)-dom.coord(som0,1))*(dom.coord(som2,0)-dom.coord(som0,0));
  double signe;
  if (prod >= 0)
    signe = 1;
  else
    signe = -1;
  // Compute p0 = 0M cross 1M along z
  p0 = (pos[0]-dom.coord(som0,0))*(pos[1]-dom.coord(som1,1))
       - (pos[1]-dom.coord(som0,1))*(pos[0]-dom.coord(som1,0));
  p0 *= signe;
  // Compute p1 = 1M cross 2M along z
  p1 = (pos[0]-dom.coord(som1,0))*(pos[1]-dom.coord(som2,1))
       - (pos[1]-dom.coord(som1,1))*(pos[0]-dom.coord(som2,0));
  p1 *= signe;
  // Compute p2 = 2M cross 0M along z
  p2 = (pos[0]-dom.coord(som2,0))*(pos[1]-dom.coord(som0,1))
       - (pos[1]-dom.coord(som2,1))*(pos[0]-dom.coord(som0,0));
  p2 *= signe;
  double epsilon=std::fabs(prod)*Objet_U::precision_geom;
  if ((p0>-epsilon) && (p1>-epsilon) && (p2>-epsilon))
    return 1;
  else
    return 0;
}


/*! @brief Returns 1 if the vertices specified by "som" are the vertices of element "element"
 *
 * in the domain associated with this geometric element. Returns 0 otherwise.
 *
 * @param som Vertex indices to compare with those of element "element".
 * @param element Index of the domain element whose vertices are to be compared.
 * @return 1 if the specified vertices are those of the given element, 0 otherwise.
 */
template <typename _SIZE_>
int Triangle_32_64<_SIZE_>::contient(const SmallArrOfTID_t& som, int_t element ) const
{
  const Domaine_t& domaine=this->mon_dom.valeur();
  if((domaine.sommet_elem(element,0)==som[0])&&
      (domaine.sommet_elem(element,1)==som[1])&&
      (domaine.sommet_elem(element,2)==som[2]))
    return 1;
  else
    return 0;
}

/*! @brief Computes the volumes (areas) of the elements of the associated domain.
 *
 * @param volumes Vector to fill with the volumes of domain elements.
 */
template <typename _SIZE_>
void Triangle_32_64<_SIZE_>::calculer_volumes(DoubleVect_t& volumes) const
{
  const Domaine_t& domaine=this->mon_dom.valeur();
  const DoubleTab_t& coord = domaine.coord_sommets();
  DoubleTab pos(3,dimension);
  int_t size_tot = domaine.nb_elem_tot();
  assert(volumes.size_totale()==size_tot);
  for (int_t num_poly=0; num_poly<size_tot; num_poly++)
    {
      for (int i=0; i<3; i++)
        {
          int_t Si = domaine.sommet_elem(num_poly,i);
          for (int j=0; j<dimension; j++)
            pos(i,j) = coord(Si,j);
        }
      volumes[num_poly] = aire_triangle(pos);
    }
}

/*! @brief Computes the face normals of the elements of the associated domain.
 *
 * @param Face_sommets Vertex indices of the faces in the domain vertex list.
 * @param face_normales Output array to fill with face normals.
 */
template <typename _SIZE_>
void Triangle_32_64<_SIZE_>::calculer_normales(const IntTab_t& Face_sommets, DoubleTab_t& face_normales) const
{
  const Domaine_t& domaine_geom = this->mon_dom.valeur();
  const DoubleTab_t& les_coords = domaine_geom.coord_sommets();
  int_t nbfaces = Face_sommets.dimension(0);
  double x1,y1;
  int_t n0,n1;
  for (int_t numface=0; numface<nbfaces; numface++)
    {
      n0 = Face_sommets(numface,0);
      n1 = Face_sommets(numface,1);
      x1 = les_coords(n0,0)-les_coords(n1,0);
      y1 = les_coords(n0,1)-les_coords(n1,1);
      face_normales(numface,0) = -y1;
      face_normales(numface,1) = x1;
    }
}

/*! @brief See ElemGeomBase::get_tab_faces_sommets_locaux.
 *
 */
template <typename _SIZE_>
int Triangle_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local) const
{
  faces_som_local.resize(3,2);
  for (int i=0; i<3; i++)
    for (int j=0; j<2; j++)
      faces_som_local(i,j) = faces_sommets_triangle[i][j];
  return 1;
}


template class Triangle_32_64<int>;
#if INT_is_64_ == 2
template class Triangle_32_64<trustIdType>;
#endif

