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

#include <Prisme.h>
#include <Domaine.h>

//      5
//     /|
//    3----4
//    | |  |
//    | 2  |
//    |/   |
//    0----1

static int faces_sommets_prisme[5][4] =
{
  { 0, 1, 3, 4 },
  { 0, 2, 3, 5 },
  { 1, 2, 4, 5 },
  { 0, 1, 2, -1 },
  { 3, 4, 5, -1 }
};


Implemente_instanciable_32_64(Prisme_32_64,"Prisme",Elem_geom_base_32_64<_T_>);

template <typename _SIZE_>
Sortie& Prisme_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  return s;
}

template <typename _SIZE_>
Entree& Prisme_32_64<_SIZE_>::readOn(Entree& s )
{
  return s;
}


/*! @brief Reorders the vertices of the Prism.
 *
 * DOES NOTHING: TO BE CODED
 *
 */
template <typename _SIZE_>
void Prisme_32_64<_SIZE_>::reordonner()
{
  Cerr << "Prisme_32_64<_SIZE_>::reordonner must be coded " << finl;
  Process::exit();
  // to be coded: restore conforming node numbering.
}


/*! @brief Returns the LML name of a prism = "PRISM6".
 *
 * @return always equal to "PRISM6"
 */
template <typename _SIZE_>
const Nom& Prisme_32_64<_SIZE_>::nom_lml() const
{
  static Nom nom="PRISM6";
  return nom;
}


/*! @brief DOES NOTHING: TO BE CODED; always returns 0.
 *
 * @brief Returns 1 if element "element" of the domain associated with this geometric element contains the point
 *               with coordinates specified by parameter "pos".
 *     Returns 0 otherwise.
 *
 * @param pos coordinates of the point to locate
 * @param ielem the element index in the domain in which the point is searched.
 * @return 1 if the point belongs to element "element", 0 otherwise
 */
template <typename _SIZE_>
int Prisme_32_64<_SIZE_>::contient(const ArrOfDouble& pos, int_t ielem ) const
{
  // to be coded:
  // does the prism with element number contain the point with
  // coordinates pos?
  // adaptation of the tetrahedron contains method
  assert(pos.size_array()==3);
  const Domaine_t& domaine=this->mon_dom.valeur();
  const Domaine_t& dom=domaine;
  double prod1,prod2,xn,yn,zn;
  int_t som0, som1, som2, som3,som4,som5;

  // First check if the sought point is not one of the
  // vertices of the triangle
  som0 = domaine.sommet_elem(ielem,0);
  som1 = domaine.sommet_elem(ielem,1);
  som2 = domaine.sommet_elem(ielem,2);
  som3 = domaine.sommet_elem(ielem,3);
  som4 = domaine.sommet_elem(ielem,4);
  som5 = domaine.sommet_elem(ielem,5);
  if( ( est_egal(dom.coord(som0,0),pos[0]) && est_egal(dom.coord(som0,1),pos[1]) && est_egal(dom.coord(som0,2),pos[2]) )
      || (est_egal(dom.coord(som1,0),pos[0]) && est_egal(dom.coord(som1,1),pos[1]) && est_egal(dom.coord(som1,2),pos[2]))
      || (est_egal(dom.coord(som2,0),pos[0]) && est_egal(dom.coord(som2,1),pos[1]) && est_egal(dom.coord(som2,2),pos[2]))
      || (est_egal(dom.coord(som3,0),pos[0]) && est_egal(dom.coord(som3,1),pos[1]) && est_egal(dom.coord(som3,2),pos[2]))
      || (est_egal(dom.coord(som4,0),pos[0]) && est_egal(dom.coord(som4,1),pos[1]) && est_egal(dom.coord(som4,2),pos[2]))
      || (est_egal(dom.coord(som5,0),pos[0]) && est_egal(dom.coord(som5,1),pos[1]) && est_egal(dom.coord(som5,2),pos[2]))
    )
    return 1;

  for (int j=0; j<5; j++)
    {
      switch(j)
        {
        case 0 :
          som0 = domaine.sommet_elem(ielem,0);
          som1 = domaine.sommet_elem(ielem,1);
          som2 = domaine.sommet_elem(ielem,3);
          som3 = domaine.sommet_elem(ielem,2);
          break;
        case 1 :
          som0 = domaine.sommet_elem(ielem,0);
          som1 = domaine.sommet_elem(ielem,2);
          som2 = domaine.sommet_elem(ielem,3);
          som3 = domaine.sommet_elem(ielem,4);
          break;
        case 2 :
          som0 = domaine.sommet_elem(ielem,1);
          som1 = domaine.sommet_elem(ielem,2);
          som2 = domaine.sommet_elem(ielem,4);
          som3 = domaine.sommet_elem(ielem,0);
          break;
        case 3 :
          som0 = domaine.sommet_elem(ielem,0);
          som1 = domaine.sommet_elem(ielem,1);
          som2 = domaine.sommet_elem(ielem,2);
          som3 = domaine.sommet_elem(ielem,3);
          break;
        case 4 :
          som0 = domaine.sommet_elem(ielem,3);
          som1 = domaine.sommet_elem(ielem,4);
          som2 = domaine.sommet_elem(ielem,5);
          som3 = domaine.sommet_elem(ielem,0);
          break;
        }

      // Algorithm: vertex 3 and point M must for j=0 to 3 be on the same side
      // as the plane formed by points som0,som1,som2.
      // compute the normal to the plane som0,som1,som2:
      xn = (dom.coord(som1,1)-dom.coord(som0,1))*(dom.coord(som2,2)-dom.coord(som0,2))
           - (dom.coord(som1,2)-dom.coord(som0,2))*(dom.coord(som2,1)-dom.coord(som0,1));
      yn = (dom.coord(som1,2)-dom.coord(som0,2))*(dom.coord(som2,0)-dom.coord(som0,0))
           - (dom.coord(som1,0)-dom.coord(som0,0))*(dom.coord(som2,2)-dom.coord(som0,2));
      zn = (dom.coord(som1,0)-dom.coord(som0,0))*(dom.coord(som2,1)-dom.coord(som0,1))
           - (dom.coord(som1,1)-dom.coord(som0,1))*(dom.coord(som2,0)-dom.coord(som0,0));
      prod1 = xn * ( dom.coord(som3,0) - dom.coord(som0,0) )
              + yn * ( dom.coord(som3,1) - dom.coord(som0,1) )
              + zn * ( dom.coord(som3,2) - dom.coord(som0,2) );
      prod2 = xn * ( pos[0] - dom.coord(som0,0) )
              + yn * ( pos[1] - dom.coord(som0,1) )
              + zn * ( pos[2] - dom.coord(som0,2) );
      // If the point is on the plane (prod2 nearly zero): cannot conclude...
      if (prod1*prod2 < 0 && std::fabs(prod2)>std::fabs(prod1)*Objet_U::precision_geom) return 0;
    }
  return 1;
}


/*! @brief DOES NOTHING: TO BE CODED, always returns 0. Returns 1 if the vertices specified by parameter "pos"
 *
 *     are the vertices of element "element" of the domain associated with
 *     the geometric element.
 *
 * @param (IntVect& pos) the vertex indices to compare with those of element "element"
 * @param (int element) the index of the domain element whose vertices are to be compared
 * @return (int) 1 if the vertices passed as parameter are those of the specified element, 0 otherwise
 */
template <typename _SIZE_>
int Prisme_32_64<_SIZE_>::contient(const SmallArrOfTID_t& pos, int_t element ) const
{
  // to be coded:
  Process::exit();
  return 0;
}


/*! @brief DOES NOTHING: TO BE CODED. Computes the volumes of the elements of the associated domain.
 *
 * @param (DoubleVect& volumes) the vector containing the volume values of the domain elements
 */
template <typename _SIZE_>
void Prisme_32_64<_SIZE_>::calculer_volumes(DoubleVect_t& volumes) const
{
  const Domaine_t& domaine=this->mon_dom.valeur();
  const IntTab_t& elem=domaine.les_elems();
  const DoubleTab_t& coord=domaine.coord_sommets();
  int_t size_tot = domaine.nb_elem_tot();
  assert(volumes.size_totale()==size_tot);
  for (int_t num_poly=0; num_poly<size_tot; num_poly++)
    {
      int_t i1=elem(num_poly,0),
            i2=elem(num_poly,1),
            i3=elem(num_poly,2),
            i4=elem(num_poly,3),
            i5=elem(num_poly,4),
            i6=elem(num_poly,5);

      // retrieved from MEDMEM_Formulae.hxx
      double a1 = (coord(i2,0)-coord(i3,0))/2.0, a2 = (coord(i2,1)-coord(i3,1))/2.0, a3 = (coord(i2,2)-coord(i3,2))/2.0;
      double b1 = (coord(i5,0)-coord(i6,0))/2.0, b2 = (coord(i5,1)-coord(i6,1))/2.0, b3 = (coord(i5,2)-coord(i6,2))/2.0;
      double c1 = (coord(i4,0)-coord(i1,0))/2.0, c2 = (coord(i4,1)-coord(i1,1))/2.0, c3 = (coord(i4,2)-coord(i1,2))/2.0;
      double d1 = (coord(i5,0)-coord(i2,0))/2.0, d2 = (coord(i5,1)-coord(i2,1))/2.0, d3 = (coord(i5,2)-coord(i2,2))/2.0;
      double e1 = (coord(i6,0)-coord(i3,0))/2.0, e2 = (coord(i6,1)-coord(i3,1))/2.0, e3 = (coord(i6,2)-coord(i3,2))/2.0;
      double f1 = (coord(i1,0)-coord(i3,0))/2.0, f2 = (coord(i1,1)-coord(i3,1))/2.0, f3 = (coord(i1,2)-coord(i3,2))/2.0;
      double h1 = (coord(i4,0)-coord(i6,0))/2.0, h2 = (coord(i4,1)-coord(i6,1))/2.0, h3 = (coord(i4,2)-coord(i6,2))/2.0;

      double A = a1*c2*f3 - a1*c3*f2 - a2*c1*f3 + a2*c3*f1 +
                 a3*c1*f2 - a3*c2*f1;
      double B = b1*c2*h3 - b1*c3*h2 - b2*c1*h3 + b2*c3*h1 +
                 b3*c1*h2 - b3*c2*h1;
      double C = (a1*c2*h3 + b1*c2*f3) - (a1*c3*h2 + b1*c3*f2) -
                 (a2*c1*h3 + b2*c1*f3) + (a2*c3*h1 + b2*c3*f1) +
                 (a3*c1*h2 + b3*c1*f2) - (a3*c2*h1 + b3*c2*f1);
      double D = a1*d2*f3 - a1*d3*f2 - a2*d1*f3 + a2*d3*f1 +
                 a3*d1*f2 - a3*d2*f1;
      double E = b1*d2*h3 - b1*d3*h2 - b2*d1*h3 + b2*d3*h1 +
                 b3*d1*h2 - b3*d2*h1;
      double F = (a1*d2*h3 + b1*d2*f3) - (a1*d3*h2 + b1*d3*f2) -
                 (a2*d1*h3 + b2*d1*f3) + (a2*d3*h1 + b2*d3*f1) +
                 (a3*d1*h2 + b3*d1*f2) - (a3*d2*h1 + b3*d2*f1);
      double G = a1*e2*f3 - a1*e3*f2 - a2*e1*f3 + a2*e3*f1 +
                 a3*e1*f2 - a3*e2*f1;
      double H = b1*e2*h3 - b1*e3*h2 - b2*e1*h3 + b2*e3*h1 +
                 b3*e1*h2 - b3*e2*h1;
      double P = (a1*e2*h3 + b1*e2*f3) - (a1*e3*h2 + b1*e3*f2) -
                 (a2*e1*h3 + b2*e1*f3) + (a2*e3*h1 + b2*e3*f1) +
                 (a3*e1*h2 + b3*e1*f2) - (a3*e2*h1 + b3*e2*f1);
      volumes[num_poly]=std::fabs(-2.0*(2.0*(A + B + D + E + G + H) + C + F + P)/9.0);

    }

}

/*! @brief Fills the array faces_som_local(i,j) giving, for 0 <= i < nb_faces() and 0 <= j < nb_som_face(i), the local vertex index
 *
 *   on the element.
 *   We have 0 <= faces_sommets_locaux(i,j) < nb_som()
 *  If all faces of the element do not have the same number of vertices, the number
 *  of columns of the array is the largest number of vertices, and the unused entries
 *  of the array are set to -1.
 *  Returns 1 if all faces have the same number of elements, 0 otherwise.
 *
 */
template <typename _SIZE_>
int Prisme_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local) const
{
  faces_som_local.resize(5,4);
  for (int i=0; i<5; i++)
    for (int j=0; j<4; j++)
      faces_som_local(i,j) = faces_sommets_prisme[i][j];
  return 1;
}


template class Prisme_32_64<int>;
#if INT_is_64_ == 2
template class Prisme_32_64<trustIdType>;
#endif

