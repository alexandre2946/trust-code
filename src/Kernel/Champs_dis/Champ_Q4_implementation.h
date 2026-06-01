/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef Champ_Q4_implementation_included
#define Champ_Q4_implementation_included

#include <Champ_implementation_divers.h>
#include <TRUSTTab.h>

class Champ_Q4_implementation: public Champ_implementation_divers
{
public:

  DoubleVect& valeur_a_elem(const DoubleVect& position, DoubleVect& val, int le_poly) const override;
  double valeur_a_elem_compo(const DoubleVect& position, int le_poly, int ncomp) const override;
  DoubleTab& valeur_aux_elems(const DoubleTab& positions, const IntVect& les_polys, DoubleTab& valeurs) const override;
  DoubleVect& valeur_aux_elems_compo(const DoubleTab& positions, const IntVect& les_polys, DoubleVect& valeurs, int ncomp) const override;
  DoubleTab& valeur_aux_sommets(const Domaine&, DoubleTab&) const override;
  DoubleVect& valeur_aux_sommets_compo(const Domaine&, DoubleVect&, int) const override;
  DoubleTab& remplir_coord_noeuds(DoubleTab& positions) const override;
  int remplir_coord_noeuds_et_polys(DoubleTab& positions, IntVect& polys) const override;
  int imprime_Q4(Sortie&, int) const;
};

/*! @brief Computes the barycentric coordinate of a point (x,y) with respect to the specified vertex of a triangle or rectangle (an element).
 *
 *     This computation concerns a 2D point.
 *
 * @param (IntTab& polys) array containing the indices of the elements with respect to which the barycentric coordinate is to be computed. polys(i,0) is the index of vertex 0 of element i in the coordinate array (coord).
 * @param (DoubleTab& coord) the coordinates of the vertices for which the barycentric coordinates are to be computed.
 * @param (double x) the first Cartesian coordinate of the point whose barycentric coordinates are to be computed
 * @param (double y) the second Cartesian coordinate of the point whose barycentric coordinates are to be computed
 * @param (int le_poly) the index of the element (in the polys array) with respect to which the barycentric coordinate will be computed.
 * @param (int i) the index of the vertex with respect to which the barycentric coordinate is wanted.
 * @return (double) the barycentric coordinate of point (x,y) with respect to the specified vertex (i) in the specified element (le_poly)
 * @throws arithmetic error, null denominator
 * @throws computation error, invalid barycentric coordinate
 */
inline double coord_barycentrique(const IntTab& polys, const DoubleTab& coord, double x, double y, int le_poly, int i)
{
  int som0, som1, som2;
  int nb_som_elem = polys.dimension(1);
  //Selection of barycentric coordinate computation depending on element type
  //Triangle case
  if (nb_som_elem == 3)
    {
      switch(i)
        {
        case 0:
          {
            som0 = polys(le_poly, 0);
            som1 = polys(le_poly, 1);
            som2 = polys(le_poly, 2);
            break;
          }
        case 1:
          {
            som0 = polys(le_poly, 1);
            som1 = polys(le_poly, 2);
            som2 = polys(le_poly, 0);
            break;
          }
        case 2:
          {
            som0 = polys(le_poly, 2);
            som1 = polys(le_poly, 0);
            som2 = polys(le_poly, 1);
            break;
          }
        default:
          {
            som0 = -1;
            som1 = -1;
            som2 = -1;
            Cerr << "Error in Champ_P1::coord_barycentrique : " << finl;
            Cerr << "A triangle does not have " << i << "nodes " << finl;
            Process::exit();
          }

        }
      double den = (coord(som2, 0) - coord(som1, 0)) * (coord(som0, 1) - coord(som1, 1)) - (coord(som2, 1) - coord(som1, 1)) * (coord(som0, 0) - coord(som1, 0));

      double num = (coord(som2, 0) - coord(som1, 0)) * (y - coord(som1, 1)) - (coord(som2, 1) - coord(som1, 1)) * (x - coord(som1, 0));

      assert(den != 0.);
      double coord_bary = num / den;
      if ((coord_bary < -Objet_U::precision_geom) || (coord_bary > 1. + Objet_U::precision_geom))
        {
          Cerr << "WARNING: The barycentric coordinate of point :" << finl;
          Cerr << "x= " << x << " y=" << y << finl;
          Cerr << "is not between 0 and 1 : " << coord_bary << finl;
          Cerr << "On the element " << le_poly << " of the processor " << Process::me() << finl;
        }
      return coord_bary;
    }
  //Rectangle case
  else if (nb_som_elem == 4)
    {
      double alpha_x = 0.;
      double alpha_y = 0.;
      double delta_x, delta_y;
//      int som0,som1,som2;
      double x0, y0;

      switch(i)
        {
        case 0:
          {
            alpha_x = -1.;
            alpha_y = -1;

            break;
          }
        case 1:
          {
            alpha_x = 1.;
            alpha_y = -1.;

            break;
          }
        case 2:
          {
            alpha_x = -1.;
            alpha_y = 1.;

            break;
          }
        case 3:
          {
            alpha_x = 1.;
            alpha_y = 1.;

            break;
          }
        default:
          {
            Cerr << "Error in Champ_P1::coord_barycentrique : " << finl;
            Process::exit();
          }

        }

      som0 = polys(le_poly, 0);
      som1 = polys(le_poly, 1);
      som2 = polys(le_poly, 2);
      delta_x = coord(som1, 0) - coord(som0, 0);
      delta_y = coord(som2, 1) - coord(som0, 1);
      x0 = coord(som0, 0) + delta_x / 2.;
      y0 = coord(som0, 1) + delta_y / 2.;
      double coord_bary = 0.25 * (1. + 2. * alpha_x * (x - x0) / delta_x) * (1. + 2. * alpha_y * (y - y0) / delta_y);

      if ((coord_bary < -Objet_U::precision_geom) || (coord_bary > 1. + Objet_U::precision_geom))
        {
          Cerr << "WARNING: The barycentric coordinate of point :" << finl;
          Cerr << "x= " << x << " y=" << y << finl;
          Cerr << "is not between 0 and 1 : " << coord_bary << finl;
          Cerr << "On the element " << le_poly << " of the processor " << Process::me() << finl;
        }
      return coord_bary;
    }

  Cerr << "The number of nodes by element " << nb_som_elem << " does not correspond to a treated situation." << finl;
  Process::exit();
  return 0.;
}

/*! @brief Computes the barycentric coordinate of a point (x,y,z) with respect to the specified vertex of a tetrahedron or hexahedron (an element).
 *
 *     This computation concerns a 3D point.
 *
 * @param (IntTab& polys) array containing the indices of the elements with respect to which the barycentric coordinate is to be computed. polys(i,0) is the index of vertex 0 of element i in the coordinate array (coord).
 * @param (DoubleTab& coord) the coordinates of the vertices for which the barycentric coordinates are to be computed.
 * @param (double x) the first Cartesian coordinate of the point whose barycentric coordinates are to be computed
 * @param (double y) the second Cartesian coordinate of the point whose barycentric coordinates are to be computed
 * @param (double z) the third Cartesian coordinate of the point whose barycentric coordinates are to be computed
 * @param (int le_poly) the index of the element (in the polys array) with respect to which the barycentric coordinate will be computed.
 * @param (int i) the index of the vertex with respect to which the barycentric coordinate is wanted.
 * @return (double) the barycentric coordinate of point (x,y,z) with respect to the specified vertex (i) in the specified element (le_poly)
 * @throws a tetrahedron does not have more than 4 vertices
 * @throws a hexahedron does not have more than 8 vertices
 * @throws arithmetic error, null denominator
 * @throws computation error, invalid barycentric coordinate
 */
inline double coord_barycentrique(const IntTab& polys, const DoubleTab& coord, double x, double y, double z, int le_poly, int i)
{
  int som0, som1, som2, som3;
  int nb_som_elem = polys.dimension(1);
  //Selection of barycentric coordinate computation depending on element type
  //Tetrahedron case
  if (nb_som_elem == 4)
    {
      switch(i)
        {
        case 0:
          {
            som0 = polys(le_poly, 0);
            som1 = polys(le_poly, 1);
            som2 = polys(le_poly, 2);
            som3 = polys(le_poly, 3);
            break;
          }
        case 1:
          {
            som0 = polys(le_poly, 1);
            som1 = polys(le_poly, 2);
            som2 = polys(le_poly, 3);
            som3 = polys(le_poly, 0);
            break;
          }
        case 2:
          {
            som0 = polys(le_poly, 2);
            som1 = polys(le_poly, 3);
            som2 = polys(le_poly, 0);
            som3 = polys(le_poly, 1);
            break;
          }
        case 3:
          {
            som0 = polys(le_poly, 3);
            som1 = polys(le_poly, 0);
            som2 = polys(le_poly, 1);
            som3 = polys(le_poly, 2);
            break;
          }
        default:
          {
            som0 = -1;
            som1 = -1;
            som2 = -1;
            som3 = -1;
            Cerr << "Error in Champ_P1::coord_barycentrique : " << finl;
            Cerr << "A tetrahedron does not have " << i << "nodes " << finl;
            Process::exit();
          }
        }

      double xp = (coord(som2, 1) - coord(som1, 1)) * (coord(som0, 2) - coord(som1, 2)) - (coord(som2, 2) - coord(som1, 2)) * (coord(som0, 1) - coord(som1, 1));
      double yp = (coord(som2, 2) - coord(som1, 2)) * (coord(som0, 0) - coord(som1, 0)) - (coord(som2, 0) - coord(som1, 0)) * (coord(som0, 2) - coord(som1, 2));
      double zp = (coord(som2, 0) - coord(som1, 0)) * (coord(som0, 1) - coord(som1, 1)) - (coord(som2, 1) - coord(som1, 1)) * (coord(som0, 0) - coord(som1, 0));
      double den = xp * (coord(som3, 0) - coord(som1, 0)) + yp * (coord(som3, 1) - coord(som1, 1)) + zp * (coord(som3, 2) - coord(som1, 2));

      xp = (coord(som2, 1) - coord(som1, 1)) * (z - coord(som1, 2)) - (coord(som2, 2) - coord(som1, 2)) * (y - coord(som1, 1));
      yp = (coord(som2, 2) - coord(som1, 2)) * (x - coord(som1, 0)) - (coord(som2, 0) - coord(som1, 0)) * (z - coord(som1, 2));
      zp = (coord(som2, 0) - coord(som1, 0)) * (y - coord(som1, 1)) - (coord(som2, 1) - coord(som1, 1)) * (x - coord(som1, 0));
      double num = xp * (coord(som3, 0) - coord(som1, 0)) + yp * (coord(som3, 1) - coord(som1, 1)) + zp * (coord(som3, 2) - coord(som1, 2));

      assert(den != 0.);
      double coord_bary = num / den;
      if ((coord_bary < -Objet_U::precision_geom) || (coord_bary > 1. + Objet_U::precision_geom))
        {
          Cerr << "WARNING: The barycentric coordinate of point :" << finl;
          Cerr << "x= " << x << " y=" << y << " z=" << z << finl;
          Cerr << "is not between 0 and 1 : " << coord_bary << finl;
          Cerr << "On the element " << le_poly << " of the processor " << Process::me() << finl;
        }
      return coord_bary;
    }
  //Hexahedron case
  else if (nb_som_elem == 8)
    {
      double alpha_x = 0.;
      double alpha_y = 0.;
      double alpha_z = 0.;
      double delta_x, delta_y, delta_z;
      //int som0,som1,som2;
      double x0, y0, z0;

      switch(i)
        {
        case 0:
          {
            alpha_x = -1.;
            alpha_y = -1.;
            alpha_z = -1.;
            break;
          }
        case 1:
          {
            alpha_x = 1.;
            alpha_y = -1.;
            alpha_z = -1.;
            break;
          }
        case 2:
          {
            alpha_x = -1.;
            alpha_y = 1.;
            alpha_z = -1.;
            break;
          }
        case 3:
          {
            alpha_x = 1.;
            alpha_y = 1.;
            alpha_z = -1.;
            break;
          }
        case 4:
          {
            alpha_x = -1.;
            alpha_y = -1;
            alpha_z = 1.;
            break;
          }
        case 5:
          {
            alpha_x = 1.;
            alpha_y = -1.;
            alpha_z = 1.;
            break;
          }
        case 6:
          {
            alpha_x = -1.;
            alpha_y = 1.;
            alpha_z = 1.;
            break;
          }
        case 7:
          {
            alpha_x = 1.;
            alpha_y = 1.;
            alpha_z = 1.;
            break;
          }
        default:
          {
            Cerr << "Error in Champ_P1::coord_barycentrique : " << finl;
            Process::exit();
          }
        }

      som0 = polys(le_poly, 0);
      som1 = polys(le_poly, 1);
      som2 = polys(le_poly, 2);
      som3 = polys(le_poly, 4);

      delta_x = coord(som1, 0) - coord(som0, 0);
      delta_y = coord(som2, 1) - coord(som0, 1);
      delta_z = coord(som3, 2) - coord(som0, 2);

      x0 = coord(som0, 0) + delta_x / 2.;
      y0 = coord(som0, 1) + delta_y / 2.;
      z0 = coord(som0, 2) + delta_z / 2.;

      double coord_bary = (1. / 8.) * (1. + 2. * alpha_x * (x - x0) / delta_x) * (1. + 2. * alpha_y * (y - y0) / delta_y) * (1. + 2. * alpha_z * (z - z0) / delta_z);

      if ((coord_bary < -Objet_U::precision_geom) || (coord_bary > 1. + Objet_U::precision_geom))
        {
          Cerr << "WARNING: The barycentric coordinate of point :" << finl;
          Cerr << "x= " << x << " y=" << y << " z=" << z << finl;
          Cerr << "is not between 0 and 1 : " << coord_bary << finl;
          Cerr << "On the element " << le_poly << " of the processor " << Process::me() << finl;
        }
      return coord_bary;
    }

  Cerr << "The number of nodes by element " << nb_som_elem << " does not correspond to a treated situation." << finl;
  Process::exit();
  return 0.;
}

#endif
