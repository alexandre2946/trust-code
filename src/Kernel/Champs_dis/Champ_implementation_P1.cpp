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

#include <Champ_implementation_P1.h>
#include <Domaine_Poly_base.h>
#include <Octree_Double.h>
#include <Domaine.h>

double Champ_implementation_P1::form_function(const ArrOfDouble& position, const IntTab& les_elems, const DoubleTab& nodes, ArrOfInt& index, int cell, int ddl) const
{
  int nb_nodes_per_cell = les_elems.dimension(1);
  assert(ddl < nb_nodes_per_cell);

  for (int i = 0; i < nb_nodes_per_cell; i++)
    index[i] = les_elems(cell, (i + ddl) % nb_nodes_per_cell);

  if (nb_nodes_per_cell == 2)
    {
      double num = 0.;
      double den = 0.;
      for (int d = 0; d < Objet_U::dimension; d++)
        {
          num += (position[d] - nodes(index[0], d)) * (position[d] - nodes(index[0], d));
          den += (nodes(index[1], d) - nodes(index[0], d)) * (nodes(index[1], d) - nodes(index[0], d));
        }
      double res = 1. - sqrt(num / den);
      // Cerr<<"ii "<<res<<" "<<index[0]<<" "<<ddl<<" "<<(ddl)%nb_nodes_per_cell<<finl;
      return res;
    }
  double num = 0.;
  double den = 0.;

  if (Objet_U::dimension == 2)
    {
      den = (nodes(index[2], 0) - nodes(index[1], 0)) * (nodes(index[0], 1) - nodes(index[1], 1)) - (nodes(index[2], 1) - nodes(index[1], 1)) * (nodes(index[0], 0) - nodes(index[1], 0));

      num = (nodes(index[2], 0) - nodes(index[1], 0)) * (position[1] - nodes(index[1], 1)) - (nodes(index[2], 1) - nodes(index[1], 1)) * (position[0] - nodes(index[1], 0));
    }
  else if (Objet_U::dimension == 3)
    {
      double xp = (nodes(index[2], 1) - nodes(index[1], 1)) * (nodes(index[0], 2) - nodes(index[1], 2)) - (nodes(index[2], 2) - nodes(index[1], 2)) * (nodes(index[0], 1) - nodes(index[1], 1));

      double yp = (nodes(index[2], 2) - nodes(index[1], 2)) * (nodes(index[0], 0) - nodes(index[1], 0)) - (nodes(index[2], 0) - nodes(index[1], 0)) * (nodes(index[0], 2) - nodes(index[1], 2));

      double zp = (nodes(index[2], 0) - nodes(index[1], 0)) * (nodes(index[0], 1) - nodes(index[1], 1)) - (nodes(index[2], 1) - nodes(index[1], 1)) * (nodes(index[0], 0) - nodes(index[1], 0));

      den = xp * (nodes(index[3], 0) - nodes(index[1], 0)) + yp * (nodes(index[3], 1) - nodes(index[1], 1)) + zp * (nodes(index[3], 2) - nodes(index[1], 2));

      xp = (nodes(index[2], 1) - nodes(index[1], 1)) * (position[2] - nodes(index[1], 2)) - (nodes(index[2], 2) - nodes(index[1], 2)) * (position[1] - nodes(index[1], 1));

      yp = (nodes(index[2], 2) - nodes(index[1], 2)) * (position[0] - nodes(index[1], 0)) - (nodes(index[2], 0) - nodes(index[1], 0)) * (position[2] - nodes(index[1], 2));

      zp = (nodes(index[2], 0) - nodes(index[1], 0)) * (position[1] - nodes(index[1], 1)) - (nodes(index[2], 1) - nodes(index[1], 1)) * (position[0] - nodes(index[1], 0));

      num = xp * (nodes(index[3], 0) - nodes(index[1], 0)) + yp * (nodes(index[3], 1) - nodes(index[1], 1)) + zp * (nodes(index[3], 2) - nodes(index[1], 2));
    }
  else
    {
      Cerr << "Error in Champ_implementation_P1::form_function : Invalid dimension" << finl;
      Process::exit();
    }

  assert(den != 0.);
  double result = num / den;

  if ((result < -Objet_U::precision_geom) || (result > 1. + Objet_U::precision_geom))
    {
      Cerr << "WARNING: The barycentric coordinate of point :" << finl;
      Cerr << "x= " << position[0] << " y=" << position[1];
      if (Objet_U::dimension == 3)
        {
          Cerr << " z=" << position[3];
        }
      Cerr << finl;
      Cerr << "is not between 0 and 1 : " << result << finl;
      Cerr << "On the element " << cell << " of the processor " << Process::me() << finl;

#ifndef NDEBUG
      Process::exit();
#endif

    }

  return result;
}
void Champ_implementation_P1::value_interpolation(const DoubleTab& positions, const ArrOfInt& cells, const DoubleTab& values, DoubleTab& resu, int ncomp) const
{
  const Domaine& domaine = get_domaine_geom();
  const Domaine_Poly_base *zpoly = sub_type(Domaine_Poly_base, get_domaine_dis()) ? &ref_cast(Domaine_Poly_base, get_domaine_dis()) : nullptr;
  const IntTab& les_elems = domaine.les_elems();
  const DoubleTab& nodes = domaine.les_sommets();
  const int nb_nodes_per_cell = domaine.nb_som_elem(), N = resu.line_size();
  ArrOfInt index(nb_nodes_per_cell);
  ArrOfDouble position(Objet_U::dimension);
  resu = 0;
  if (zpoly) //polyhedron -> volume-weighted interpolation
    {
      const DoubleTab& v_es = zpoly->vol_elem_som();
      const DoubleVect& ve = zpoly->volumes();
      const IntTab& es_d = zpoly->elem_som_d();
      for (int ic = 0; ic < cells.size_array(); ic++)
        {
          int cell = cells[ic];
          if (cell < 0)
            continue;
          assert(cell >= 0);
          assert(cell < les_elems.dimension_tot(0));
          if (ncomp != -1)
            for (int j = 0, k = es_d(cell); k < es_d(cell + 1); j++, k++)
              resu(ic) += values(les_elems(cell, j), ncomp) * v_es(k) / ve(cell);
          else
            {
              assert(values.line_size() == N);
              for (int j = 0, k = es_d(cell); k < es_d(cell + 1); k++)
                for (int n = 0, s = les_elems(cell, j); n < N; n++)
                  resu(ic, n) += values(s, n) * v_es(k) / ve(cell);
            }

        }
    }
  else
    for (int ic = 0; ic < cells.size_array(); ic++)
      {
        int cell = cells[ic];
        if (cell < 0)
          continue;
        for (int k = 0; k < Objet_U::dimension; k++)
          position[k] = positions(ic, k);

        assert(cell >= 0);
        assert(cell < les_elems.dimension_tot(0));
        if (ncomp != -1)
          {
            for (int j = 0; j < nb_nodes_per_cell; j++)
              {
                int node = les_elems(cell, j);
                resu(ic) += values(node, ncomp) * form_function(position, les_elems, nodes, index, cell, j);
              }
          }
        else
          {

            assert(values.line_size() == N);
            for (int j = 0; j < nb_nodes_per_cell; j++)
              {
                double weight = form_function(position, les_elems, nodes, index, cell, j);
                int node = les_elems(cell, j);
                for (int k = 0; k < N; k++)
                  resu(ic, k) += values(node, k) * weight;
              }
          }
      }
}

/*! @brief Initializes the array of values at the vertices of domain dom from values read from "input".
 *
 * (the array is sized, associated with the parallel structure, and filled with real and virtual values)
 *   The file must have the following format (n is the number of nodal values
 *   stored, x, y, z are the coordinates of the vertices, compo1... are the values
 *   of the components.
 *   n  (int)
 *   x y [z] compo1 [compo2 [compo3 ... ]]   (type double)
 *
 */
void Champ_implementation_P1::init_from_file(DoubleTab& val, const Domaine& dom, int nb_comp, double tolerance, Entree& input)
{
  val.resize(0, nb_comp);
  dom.creer_tableau_sommets(val, RESIZE_OPTIONS::NOCOPY_NOINIT);

  // Build an octree with the domain vertices:
  const DoubleTab& coord = dom.coord_sommets();
  Octree_Double octree;
  octree.build_nodes(coord, 0 /* do not include virtual nodes */);

  // Read values from input
  int nb_val_lues;
  input >> nb_val_lues;
  const int dim = coord.dimension(1);
  const int tmp_size = dim + nb_comp;
  ArrOfDouble tmp(tmp_size);
  ArrOfDouble node_coord; // points to the dim first elements of "tmp" (coordinates of the node)
  node_coord.ref_array(tmp, 0 /* start index */, dim /* size */);
  ArrOfInt node_list;


  ArrOfInt count(dom.nb_som()); // number of times this coordinate has been found

  for (int i_val = 0; i_val < nb_val_lues; i_val++)
    {
      input.get(tmp.addr(), tmp_size);
      // first call returns more points (some might be at a larger distance)
      octree.search_elements_box(node_coord, tolerance, node_list);
      octree.search_nodes_close_to(node_coord, coord, node_list, tolerance);
      const int n = node_list.size_array();
      if (n > 1)
        {
          Cerr << "Error in Champ_som_lu::readOn: point " << node_coord << " corresponds to " << node_list.size_array() << " nodes in the geometry within the specified tolerance" << finl;
          Process::exit();
        }
      if (n == 1)
        {
          const int node_index = node_list[0];
          count[node_index]++;
          for (int i = 0; i < nb_comp; i++)
            val(node_index, i) = tmp[dim + i];
        }
      else
        {
          // This vertex is not on this processor...
        }
    }

  if (dom.nb_som() > 0 && min_array(count) == 0)
    {
      Cerr << "Error in Champ_som_lu::readOn: some coordinates were not found in the file" << finl;
      Process::exit();
    }
  val.echange_espace_virtuel();
}

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
double coord_barycentrique_P1(const IntTab& polys, const DoubleTab& coord, double x, double y, int le_poly, int i)
{
  int nb_som_elem = polys.dimension(1);
  //Selection of barycentric coordinate computation depending on element type
  //Triangle case
  if (nb_som_elem == 3)
    return coord_barycentrique_P1_triangle(polys, coord, x, y, le_poly, i);
  //Rectangle case
  else if (nb_som_elem == 4)
    return coord_barycentrique_P1_rectangle(polys, coord, x, y, le_poly, i);
  Cerr << "The number of nodes by element " << nb_som_elem << " does not correspond to a treated situation in the coord_barycentrique_P1 function." << finl;
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
double coord_barycentrique_P1(const IntTab& polys, const DoubleTab& coord, double x, double y, double z, int le_poly, int i)
{
  int nb_som_elem = polys.dimension(1);
  //Selection of barycentric coordinate computation depending on element type
  //Tetrahedron case
  if (nb_som_elem == 4)
    return coord_barycentrique_P1_tetraedre(polys, coord, x, y, z, le_poly, i);
  else if (nb_som_elem == 8)
    return coord_barycentrique_P1_hexaedre(polys, coord, x, y, z, le_poly, i);
  Cerr << "The number of nodes by element " << nb_som_elem << " does not correspond to a treated situation in the coord_barycentrique_P1 function." << finl;
  Process::exit();
  return 0.;
}
