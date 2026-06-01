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

#include <Tetraedre.h>
#include <Domaine.h>
#include <Linear_algebra_tools_impl.h>
#include <algorithm>
using std::swap;

Implemente_instanciable_32_64(Tetraedre_32_64,"Tetraedre",Elem_geom_base_32_64<_T_>);

static int faces_sommets_tetra[4][3] =
{
  { 1, 2, 3 },
  { 2, 3, 0 },
  { 3, 0, 1 },
  { 0, 1, 2 }
};

template <typename _SIZE_>
Sortie& Tetraedre_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  return s;
}

template <typename _SIZE_>
Entree& Tetraedre_32_64<_SIZE_>::readOn(Entree& s )
{
  return s;
}


/*! @brief Returns the LML name of a tetrahedron = "TETRA4".
 *
 * @return Always equal to "TETRA4".
 */
template <typename _SIZE_>
const Nom& Tetraedre_32_64<_SIZE_>::nom_lml() const
{
  static Nom nom="TETRA4";
  return nom;
}


namespace
{
/*! @brief tests if 2 points are on the same side of a plane defined by three points
*
* Takes the coordinates of all points involved as arguments (5 points, so 15 arguments)
* The order is X, Y, Z coord of a point, then next point
*
* The first nine arguments are for the points defining the plane (X/Y/Z 0 to 2)
*
* The next 6 describe the point for which we want to test they are on the same side (X3/Y3/Z3 and Mx/My/Mz)
*
* The use case is testing if a point M is inside a tetrahedra.
* To do that, call this function 4 times in a row while cycling the first 4 points,
* which must correspond to the 4 vertexes of the tetrahedra,
* as done in function Tetraedre_32_64<_SIZE_>::contient
* (hence the names of the arguments, which may be confusing for a different use case)
*
*
*
* @return 1 if the point belongs to the tetrahedron, 0 otherwise.
*/
inline bool is_on_same_side_of_plane(const double& X0, const double& Y0, const double& Z0,
                                     const double& X1, const double& Y1, const double& Z1,
                                     const double& X2, const double& Y2, const double& Z2,
                                     const double& X3, const double& Y3, const double& Z3,
                                     const double& Mx, const double& My, const double& Mz
                                    )
{

  // computes the normal vector of the plane
  double xn = (Y1 - Y0) * (Z2 - Z0) - (Z1 - Z0) * (Y2 - Y0);
  double yn = (Z1 - Z0) * (X2 - X0) - (X1 - X0) * (Z2 - Z0);
  double zn = (X1 - X0) * (Y2 - Y0) - (Y1 - Y0) * (X2 - X0);

  // computes the scalar product between normal vector and a vector from the plane to each of the points we want to test
  double prod1 = xn * (X3 - X0) + yn * (Y3 - Y0) + zn * (Z3 - Z0);
  double prod2 = xn * (Mx - X0) + yn * (My - Y0) + zn * (Mz - Z0);

  // if scalar products have the same sign, the points are on the same sides
  // we allow a slight tolerance for points very close to the plane
  if (prod1 * prod2 < 0 && std::fabs(prod2)>std::fabs(prod1)*Objet_U::precision_geom)
    {
      return false;
    }
  else
    {
      return true;
    }
}
}

/*! @brief Returns 1 if element "ielem" of the domain associated with this geometric element contains the point with coordinates "pos". Returns 0 otherwise.
 *
 * @param pos Coordinates of the point to locate.
 * @param ielem Index of the domain element in which to search for the point.
 * @return 1 if the point belongs to element "ielem", 0 otherwise.
 */
template <typename _SIZE_>
int Tetraedre_32_64<_SIZE_>::contient(const ArrOfDouble& pos, int_t ielem) const
{
  // 29/01/2010 CPU optimisation of this method (50% faster) by PL
  assert(pos.size_array()==3);
  const Domaine_t& domaine=mon_dom.valeur();
  const DoubleTab_t& coord=domaine.coord_sommets();

  int_t som0 = domaine.sommet_elem(ielem,0);
  int_t som1 = domaine.sommet_elem(ielem,1);
  int_t som2 = domaine.sommet_elem(ielem,2);
  int_t som3 = domaine.sommet_elem(ielem,3);
  double X0 = coord(som0,0);
  double Y0 = coord(som0,1);
  double Z0 = coord(som0,2);
  double X1 = coord(som1,0);
  double Y1 = coord(som1,1);
  double Z1 = coord(som1,2);
  double X2 = coord(som2,0);
  double Y2 = coord(som2,1);
  double Z2 = coord(som2,2);
  double X3 = coord(som3,0);
  double Y3 = coord(som3,1);
  double Z3 = coord(som3,2);

  // Here we used to test if the point was one of the vertexes of the tetra using est_egal
  // probably not worth it, happened in 0.03% of test cases according to gcov
  // must mean we rarely lookup for a point of the mesh using this function


  // However, it might be worth to check a simpler distance to center of tetra first
  // Using some bound at which we are certain the point is outside (one that is easier to compute than circumradius preferably, don't know if that exists)
  // depending on usage of this function, may avoid testing on each face in a lot of cases

  // Now we do the real work
  // For a point to be inside a tetra, for each face made of three of the 4 vertexes
  // the point must be on the same side as the fourth vertex
  // We test that with the function is_on_same_side_of_plane defined in this file


  // test som3 and pos are on same side
  if (not is_on_same_side_of_plane(X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, pos[0], pos[1], pos[2]))
    {
      return false;
    }

  // test som2 and pos are on same side
  if (not is_on_same_side_of_plane(X3, Y3, Z3, X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, pos[0], pos[1], pos[2]))
    {
      return false;
    }

  // test som1 and pos are on same side
  if (not is_on_same_side_of_plane(X2, Y2, Z2, X3, Y3, Z3, X0, Y0, Z0, X1, Y1, Z1, pos[0], pos[1], pos[2]))
    {
      return false;
    }

  // test som0 and pos are on same side
  if (not is_on_same_side_of_plane(X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X0, Y0, Z0, pos[0], pos[1], pos[2]))
    {
      return false;
    }

  return true;

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
int Tetraedre_32_64<_SIZE_>::contient(const SmallArrOfTID_t& som, int_t element ) const
{
  const Domaine_t& domaine=mon_dom.valeur();
  if((domaine.sommet_elem(element,0)==som[0])&&
      (domaine.sommet_elem(element,1)==som[1])&&
      (domaine.sommet_elem(element,1)==som[2])&&
      (domaine.sommet_elem(element,1)==som[3]))
    return 1;
  else
    return 0;
}

/*! @brief Computes the volumes of the elements of the associated domain.
 *
 * @param tab_volumes Vector to fill with the volumes of domain elements.
 */
template <typename _SIZE_>
void Tetraedre_32_64<_SIZE_>::calculer_volumes(DoubleVect_t& tab_volumes) const
{
  const Domaine_t& domaine=mon_dom.valeur();

  int_t size_tot = domaine.nb_elem_tot();
  assert(tab_volumes.size_totale()==size_tot);
  ConstView<_SIZE_,2> les_Polys = domaine.les_elems().view_ro();
  CDoubleTabView coord = domaine.coord_sommets().view_ro();
  auto volumes = tab_volumes.view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), size_tot, KOKKOS_LAMBDA(const int_t num_poly)
  {
    int_t s0 = les_Polys(num_poly, 0);
    int_t s1 = les_Polys(num_poly, 1);
    int_t s2 = les_Polys(num_poly, 2);
    int_t s3 = les_Polys(num_poly, 3);
    double x0 = coord(s0, 0), y0 = coord(s0, 1), z0 = coord(s0, 2);
    double x1 = coord(s1, 0), y1 = coord(s1, 1), z1 = coord(s1, 2);
    double x2 = coord(s2, 0), y2 = coord(s2, 1), z2 = coord(s2, 2);
    double x3 = coord(s3, 0), y3 = coord(s3, 1), z3 = coord(s3, 2);
    volumes(num_poly) = Kokkos::fabs((x1-x0)*((y2-y0)*(z3-z0)-(y3-y0)*(z2-z0))-
                                     (x2-x0)*((y1-y0)*(z3-z0)-(y3-y0)*(z1-z0))+
                                     (x3-x0)*((y1-y0)*(z2-z0)-(y2-y0)*(z1-z0)))/6;
  });
  end_gpu_timer(__KERNEL_NAME__);
}


/*! @brief Computes the face normals of the elements of the associated domain.
 *
 * @param Face_sommets Vertex indices of the faces in the domain vertex list.
 * @param face_normales Output array to fill with face normals.
 */
template <typename _SIZE_>
void Tetraedre_32_64<_SIZE_>::calculer_normales(const IntTab_t& Face_sommets, DoubleTab_t& face_normales) const
{
  const Domaine_t& domaine_geom = mon_dom.valeur();
  const DoubleTab_t& les_coords = domaine_geom.coord_sommets();
  int_t nbfaces = Face_sommets.dimension(0);
  for (int_t numface=0; numface<nbfaces; numface++)
    {

      int_t n0 = Face_sommets(numface,0);
      int_t n1 = Face_sommets(numface,1);
      int_t n2 = Face_sommets(numface,2);

      double x1 = les_coords(n0,0) - les_coords(n1,0);
      double y1 = les_coords(n0,1) - les_coords(n1,1);
      double z1 = les_coords(n0,2) - les_coords(n1,2);

      double x2 = les_coords(n2,0) - les_coords(n1,0);
      double y2 = les_coords(n2,1) - les_coords(n1,1);
      double z2 = les_coords(n2,2) - les_coords(n1,2);

      face_normales(numface,0) = (y1*z2 - y2*z1)/2;
      face_normales(numface,1) = (-x1*z2 + x2*z1)/2;
      face_normales(numface,2) = (x1*y2 - x2*y1)/2;
    }
}

/*! @brief See ElemGeomBase::get_tab_faces_sommets_locaux.
 */
template <typename _SIZE_>
int Tetraedre_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local) const
{
  // a tetrahedron has four faces of three vertices each
  faces_som_local.resize(4,3);
  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      faces_som_local(i,j) = faces_sommets_tetra[i][j];
  return 1;
}

template <typename _SIZE_>
void Tetraedre_32_64<_SIZE_>::get_tab_aretes_sommets_locaux(IntTab& tab) const
{
  // a tetrahedron has six edges of two vertices each
  tab.resize(6, 2);
  int count = 0;
  // one edge between each pair of tetra vertices: n * (n-1) / 2 edges with n=4
  for (int i = 0; i < 3; i++)
    {
      for (int j = i + 1; j < 4; j++)
        {
          tab(count, 0) = i;
          tab(count, 1) = j;
          count++;
        }
    }
  assert(count == 6);
}


///*! Computes the barycentric coordinate in a tetrahedron corresponding to a
// * Cartesian coordinate "point". Note: if "point" is outside the tetra, one or more
// * barycentric coordinates will be negative.
// * polys is the tetrahedron connectivity table (vertex indices),
// * coords is the vertex coordinate table,
// * le_poly is the tetrahedron index whose barycentric coordinates are to be computed.
// *
// * The result is stored in coord_bary (weights of the first three vertices,
// * the fourth being implicitly 1 minus the sum of the other three).
// * If epsilon is non-zero, the return value is the uncertainty on the barycentric
// * coordinates for an uncertainty epsilon on the Cartesian coordinates.
// * (computed in Linfini norm, i.e. the max error over each component)
// */
//template <typename _SIZE_>
//double Tetraedre_32_64<_SIZE_>::coord_bary(const IntTab& polys, const DoubleTab& coords,
//                             const Vecteur3& point, int le_poly, Vecteur3& coord_bary, double epsilon)
//{
//  Matrice33 m;
//  Vecteur3 origine;
//  matrice_base_tetraedre(polys, coords, le_poly, m, origine);
//  Matrice33 inverse_m;
//  Matrice33::inverse(m, inverse_m);
//  Vecteur3 v(point-origine);
//  Matrice33::produit(inverse_m, v, coord_bary);
//
//  double resu;
//  if (epsilon > 0.)
//    {
//      // An error epsilon on the "point" coordinate results in an error on coord_bary:
//      double norm = inverse_m.norme_Linfini();
//      resu = norm * epsilon;
//    }
//  else
//    {
//      resu = 0.;
//    }
//  return resu;
//}


template class Tetraedre_32_64<int>;
#if INT_is_64_ == 2
template class Tetraedre_32_64<trustIdType>;
#endif

