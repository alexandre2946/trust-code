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

#include <Reordonner_faces_periodiques.h>
#include <Connectivite_som_elem.h>
#include <Partitionneur_base.h>
#include <communications.h>
#include <Array_tools.h>
#include <TRUSTLists.h>
#include <ArrOfBit.h>
#include <Domaine.h>
#include <Param.h>

Implemente_base_32_64(Partitionneur_base_32_64,"Partitionneur_base",Objet_U_With_Params);
// XD partitionneur_deriv objet_u partitionneur_deriv INHERITS_BRACE not_set
// XD attr nb_parts entier nb_parts OPT The number of non empty parts that must be generated (generally equal to the
// XD_CONT number of processors in the parallel run).

template <typename _SIZE_>
Sortie& Partitionneur_base_32_64<_SIZE_>::printOn(Sortie& os) const
{
  exit();
  return os;
}

template <typename _SIZE_>
void Partitionneur_base_32_64<_SIZE_>::set_param(Param& param) const
{
}

/*! @brief Corrects the partition so that element 0 of the initial domain is on the first sub-domain of the partition.
 *
 *   The first sub-domain and the one containing element 0 are swapped.
 *
 */
template <typename _SIZE_>
void Partitionneur_base_32_64<_SIZE_>::corriger_elem0_sur_proc0(BigIntVect_& elem_part)
{
  Cerr << "Correction of the splitting to put the element 0 on processor 0." << finl;
  int pe_to_xchange =  elem_part[0];
  envoyer_broadcast(pe_to_xchange, 0);
  if (pe_to_xchange == 0)
    {
      Cerr << " No correction to be made" << finl;
      return;
    }

  Cerr << " Exchange of parts 0 and " << pe_to_xchange << finl;
  const int_t n = elem_part.size_reelle();
  for (int_t i = 0; i < n; i++)
    {
      const int pe = elem_part[i];
      if (pe == 0)
        elem_part[i] = pe_to_xchange;
      else if (pe == pe_to_xchange)
        elem_part[i] = 0;
    }

}

namespace
{
/*! @brief Builds (size and content) the elements array with, for each face of the given boundary, the index of the adjacent domain element.
 *
 * @param som_elem Vertex-element connectivity for the domain, computed using construire_connectivite_som_elem.
 * @param faces The faces of the boundary to process (for each face, indices of the vertices).
 * @param nom_faces A boundary name to print in case of error.
 * @param elements The array to fill.
 */
template <typename _SIZE_>
void chercher_elems_voisins_faces(const Static_Int_Lists_32_64<_SIZE_>& som_elem,
                                  const IntTab_T<_SIZE_>& faces,
                                  const Nom& nom_faces,
                                  ArrOfInt_T<_SIZE_>& elements)
{
  using int_t = _SIZE_;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;

  const int_t nb_faces = faces.dimension(0);
  elements.resize_array(nb_faces);
  if (nb_faces == 0)
    return;
  const int nb_som_faces = faces.dimension_int(1);
  SmallArrOfTID_t une_face(nb_som_faces);
  SmallArrOfTID_t voisins;
  for (int_t i = 0; i < nb_faces; i++)
    {
      for (int j = 0; j < nb_som_faces; j++)
        une_face[j] = faces(i, j);
      find_adjacent_elements(som_elem, une_face, voisins);
      const int nb_voisins = voisins.size_array();
      if (nb_voisins != 1)
        {
          Cerr << "Error in chercher_elems_voisins_faces : the face " << i
               << "\n of boundary " << nom_faces << " has " << nb_voisins
               << " neighboring elements with the indices : " << voisins << finl;
          Process::exit();
        }
      elements[i] = voisins[0];
    }
}
}

/*! @brief Computes a connectivity graph between elements connected by periodic faces.
 *
 * If element i is a neighbour of element j through a periodic face, then there exists
 *   k such that graph(i,k)==j and there exists k2 such that graph(j,k2)==i.
 *
 * @param domaine The domain to process.
 * @param som_elem Vertex-element connectivity for the given domain. WARNING: periodic boundary faces are assumed to be ordered according to the periodic boundary convention. See check_faces_periodiques().
 * @param my_offset Global element offset for this process.
 * @param graph Where the result is stored. Return value: number of elements in the graph (equal to the number of periodic faces).
 */
template <typename _SIZE_>
typename Partitionneur_base_32_64<_SIZE_>::int_t
Partitionneur_base_32_64<_SIZE_>::calculer_graphe_connexions_periodiques(const Domaine_t& domaine,
                                                                         const Static_Int_Lists_t& som_elem, const int_t my_offset,
                                                                         Static_Int_Lists_t& graph)
{
  const Noms& liste_bords_periodiques = domaine.bords_perio();
  const int_t nb_elem = domaine.nb_elem();

  // For each element, how many periodic faces does it have?
  ArrOfInt_t nb_faces_perio(nb_elem);
  // List of correspondences element0 <=> element1
  // between the element adjacent to a face and the element adjacent to the opposite periodic face
  IntTab_t correspondances(0,2);

  // First step: fill nb_faces_perio and correspondances
  // Loop over periodic boundaries
  const int nb_bords = domaine.nb_bords();
  for (int i_bord = 0; i_bord < nb_bords; i_bord++)
    {
      const Bord_t& bord = domaine.bord(i_bord);
      if (!liste_bords_periodiques.contient_(bord.le_nom()))
        continue;
      Cerr << " Checking of the boundary " << bord.le_nom();
      {
        ArrOfDouble delta;
        ArrOfDouble erreur;
        const int ok = Reordonner_faces_periodiques_32_64<_SIZE_>::check_faces_periodiques(domaine.bord(i_bord), delta, erreur);
        const int d = delta.size_array();
        Cerr << " Delta = ";
        for (int i = 0; i < d; i++) Cerr << delta[i] << " ";
        Cerr << " Error = ";
        for (int i = 0; i < d; i++) Cerr << erreur[i] << " ";
        Cerr << finl;
        if (!ok)
          {
            Cerr << "You need to use the Declarer_bord_perio keyword on the periodic boundaries." << finl;
            Cerr << "See the reference manual to use this keyword on your data file." << finl;
            exit();
          }
      }
      ArrOfInt_t elems_voisins;
      chercher_elems_voisins_faces<_SIZE_>(som_elem, bord.faces().les_sommets(), bord.le_nom(), elems_voisins);

      // Loop over the periodic boundary faces, two by two.
      // It is assumed that faces appear in the order:
      //  first all faces from one end of the domain,
      //  then, in the same order, the faces from the other end.
      assert(bord.nb_faces() % 2 == 0); // Even count, necessarily
      const int_t nb_faces = bord.nb_faces() / 2;
      for (int_t i = 0; i < nb_faces; i++)
        {
          // Indices of the two elements "neighbouring" through the periodic face:
          int_t elem0 = elems_voisins[i];
          int_t elem1 = elems_voisins[i+nb_faces]; // Index of the corresponding periodic face

          ++nb_faces_perio[elem0 - my_offset*(elem0 >= my_offset)];
          ++nb_faces_perio[elem1 - my_offset*(elem1 >= my_offset)];
          if (elem0 == elem1)
            {
              Cerr << "Error in calculer_correspondance_faces_perio: the faces " << i
                   << " and " << i + nb_faces
                   << "\n of the boundary " << bord.le_nom()
                   << " are neighbors of the same element " << elem0 << finl;
            }
          const int_t n = correspondances.dimension(0);
          correspondances.resize(n+1, 2);
          correspondances(n, 0) = elem0;
          correspondances(n, 1) = elem1;
        }
    }

  // Second step:
  // Build "graph" from the correspondances array.
  graph.set_list_sizes(nb_faces_perio);
  // Reuse the nb_faces_perio array to store the number
  // of elements already filled in each list:
  nb_faces_perio = 0;
  const int_t n = correspondances.dimension(0);
  for (int_t i = 0; i < n; i++)
    {
      const int_t elem0 = correspondances(i, 0),
                  elem1 = correspondances(i, 1);

      const int_t j0 = nb_faces_perio[elem0 - my_offset*(elem0 >= my_offset)]++;
      graph.set_value(elem0 - my_offset*(elem0 >= my_offset), j0, elem1);
      const int_t j1 = nb_faces_perio[elem1 - my_offset*(elem1 >= my_offset)]++;
      graph.set_value(elem1 - my_offset*(elem1 >= my_offset), j1, elem0);
    }
  Cerr << " There is " << n*2 << " periodic connections." << finl;
  return n * 2;
}

/*! @brief Modifies elem_part to ensure the following properties: 1) Elements that have a boundary vertex are associated
 *
 *      with a processor that owns an adjacent boundary face.
 *   2) If a processor owns a real periodic vertex, it necessarily
 *      also owns the associated renum_som_perio (hence an element
 *      that has this vertex and a periodic face).
 *   This property is essential for periodicity (existence
 *   of renum_som_perio for all vertices).
 *   For other boundaries, this correction may be unnecessary, but
 *   this is not certain. Without this correction, isolated boundary
 *   vertices can exist (a processor owns a boundary vertex but no face).
 *   If boundary vertices are searched by scanning boundary faces,
 *   the result is wrong. With this correction, that algorithm is correct.
 *
 */
template <typename _SIZE_>
typename Partitionneur_base_32_64<_SIZE_>::int_t
Partitionneur_base_32_64<_SIZE_>::corriger_sommets_bord(const Domaine_t& domaine,
                                                        const ArrOfInt_t& renum_som_perio,
                                                        const Static_Int_Lists_t& som_elem,
                                                        BigIntVect_& elem_part)
{
  using ArrOfBit_t = ArrOfBit_32_64<_SIZE_>;

  const int_t nb_som_tot = domaine.nb_som_tot();
  const int_t nb_elem = domaine.nb_elem();
  const int_t nb_elem_tot = domaine.nb_elem_tot();
  const Noms& liste_bords_perio = domaine.bords_perio();

  // First step:
  // Mark boundary vertices:
  ArrOfBit_t sommet_bord(nb_som_tot);
  ArrOfBit_t sommet_bord_perio(nb_som_tot);
  sommet_bord = 0;
  sommet_bord_perio = 0;
  // element_bord indicates whether the element is adjacent to a boundary face
  ArrOfBit_t element_bord(nb_elem_tot);
  // element_bord_perio indicates whether the element is adjacent to a periodic boundary face
  ArrOfBit_t element_bord_perio(nb_elem_tot);
  element_bord = 0;
  element_bord_perio = 0;

  const int nb_bords = domaine.nb_bords();
  for (int i_bord = 0; i_bord < nb_bords; i_bord++)
    {
      const Bord_t& bord = domaine.bord(i_bord);
      const IntTab_t& faces_sommets = bord.faces().les_sommets();
      const int_t nb_faces_bord = faces_sommets.dimension(0);
      const int nb_som_face = faces_sommets.dimension_int(1);
      ArrOfInt_t elems_voisins;
      const bool is_perio = (liste_bords_perio.contient_(bord.le_nom()));

      chercher_elems_voisins_faces(som_elem, faces_sommets, bord.le_nom(), elems_voisins);

      // For each element neighbouring the boundary faces, mark that element.
      // Associate to each boundary vertex the index of the smallest part
      // that contains an adjacent boundary face.
      for (int_t i_face = 0; i_face < nb_faces_bord; i_face++)
        {
          const int_t elem_voisin = elems_voisins[i_face];
          element_bord.setbit(elem_voisin);
          if (is_perio)
            element_bord_perio.setbit(elem_voisin);
          for (int i_som = 0; i_som < nb_som_face; i_som++)
            {
              const int_t som = faces_sommets(i_face, i_som);
              sommet_bord.setbit(som);
              if (is_perio)
                sommet_bord_perio.setbit(som);
            }
        }
    }

  // Second step:
  // Loop over elements that have no boundary face but do have
  //  boundary vertices:
  //  For each boundary vertex of the element, build the list of
  //  "allowed parts":
  //    If the renum_som_perio vertex is adjacent to a periodic boundary face,
  //    the allowed parts are those that own a periodic face adjacent to the vertex.
  //    Otherwise, they are the parts containing an adjacent boundary face.
  //  Compute the intersection of these lists and if the element does not belong
  //  to an allowed part, choose one and assign it.
  int_t count = 0;
  {
    BigArrOfInt_ parties_autorisees, tmp;

    const IntTab_t& elements = domaine.les_elems();
    const int nb_som_elem = elements.dimension_int(1);
    for (int_t elem = 0; elem < nb_elem; elem++)
      {
        // Does the element have a periodic vertex?
        bool has_som_perio = false;
        int isom;
        for (isom = 0; isom < nb_som_elem; isom++)
          if (sommet_bord_perio[elements(elem, isom)])
            has_som_perio = true;
        // Loop over element vertices:
        parties_autorisees.resize_array(0);
        int nb_sommets_bord = 0; // Number of boundary vertices of the element
        for (isom = 0; isom < nb_som_elem; isom++)
          {
            const int_t som = elements(elem, isom);
            // Process only periodic boundary vertices if the element has a periodic boundary,
            // otherwise process only boundary vertices.
            if (has_som_perio && !sommet_bord_perio[som])
              continue;
            if (!sommet_bord[som])
              continue;

            nb_sommets_bord++;
            // Put in tmp the list of parts associated with adjacent boundary elements
            {
              tmp.resize_array(0);
              const int_t renum_som = renum_som_perio[som];
              const int_t n = som_elem.get_list_size(renum_som);
              for (int_t i = 0; i < n; i++)
                {
                  const int_t elem2 = som_elem(renum_som, i);
                  int_t test;
                  if (has_som_perio)
                    test = element_bord_perio[elem2];
                  else
                    test = element_bord[elem2];
                  if (test)
                    {
                      const int p = elem_part[elem2];
                      tmp.append_array(p);
                    }
                }
              array_trier_retirer_doublons(tmp);
            }
            // Compute the intersection between tmp and parties_autorisees
            if (parties_autorisees.size_array() > 0)
              array_calculer_intersection(parties_autorisees, tmp);
            else
              parties_autorisees = tmp;
          }
        if (nb_sommets_bord > 0)
          {
            // Does the element belong to an allowed part?
            const int_t n = parties_autorisees.size_array();
            if (n == 0)
              {
                Cerr << "Partitionneur_base_32_64<_SIZE_>::corriger_sommets_bord : Error. No part authorized because of the periodicity for the mesh element " << elem << finl;
                elem_part[elem] = -1;
              }
            else
              {
                int_t i;
                const int p = elem_part[elem];
                for (i = 0; i < n; i++)
                  if (parties_autorisees[i] == p)
                    break;
                if (i >= n)
                  {
                    // Another part must be assigned:
                    elem_part[elem] = parties_autorisees[0];
                    count++;
                  }
              }
          }
      }
    Cerr << "Partitionneur_base_32_64<_SIZE_>::corriger_sommets_bord : " << count << " modified elements" << finl;
  }
  return count;
}

/*! @brief Applies corrections to elem_part so that multi-periodicity is correct:
 *
 *   If a vertex belongs to several periodic boundaries,
 *   all adjacent elements are assigned to the same processor.
 *
 */
template <typename _SIZE_>
typename Partitionneur_base_32_64<_SIZE_>::int_t
Partitionneur_base_32_64<_SIZE_>::corriger_multiperiodique(const Domaine_t& domaine,
                                                           const ArrOfInt_t& renum_som_perio,
                                                           const Static_Int_Lists_t& som_elem,
                                                           BigIntVect_& elem_part)
{
  const int_t nb_som = domaine.nb_som();
  const int_t nb_elem = domaine.nb_elem();
  const Noms& liste_bords_perio = domaine.bords_perio();

  // For each periodic vertex, select a part to which it belongs
  // (the smallest among the parts of the adjacent periodic elements)
  // Initialised to -1
  BigArrOfInt_ partie_associee(nb_som);
  partie_associee= -1;
  // For each vertex, to which periodic boundary(ies) does it belong?
  // (sum of 2^n where n is the index of the boundary name in the periodic boundary list)
  // Initialised to 0
  BigArrOfInt_ marqueur_bord(nb_som);

  int deux_puissance_i_bord = 1;
  for (auto& itr : liste_bords_perio)
    {
      const Bord_t& bord = domaine.bord(itr);
      const IntTab_t& faces_sommets = bord.faces().les_sommets();
      const int_t nb_faces_bord = faces_sommets.dimension(0);
      const int nb_som_face = faces_sommets.dimension_int(1);
      ArrOfInt_t elems_voisins;
      chercher_elems_voisins_faces<_SIZE_>(som_elem, faces_sommets, bord.le_nom(), elems_voisins);
      for (int_t i_face = 0; i_face < nb_faces_bord; i_face++)
        {
          const int_t elem_voisin = elems_voisins[i_face];
          const int part = elem_part[elem_voisin];
          for (int i_som = 0; i_som < nb_som_face; i_som++)
            {
              const int_t som = faces_sommets(i_face, i_som);
              const int old_part = partie_associee[som];
              if (old_part < 0 || old_part > part)
                partie_associee[som] = part;
              marqueur_bord[som] |= deux_puissance_i_bord; // Bitwise OR operator
            }
        }
      deux_puissance_i_bord *= 2;
      if (deux_puissance_i_bord > 65536)
        {
          // In any case, more than 3 periodic boundaries is highly suspicious...
          Cerr << "Error in Partitionneur_base_32_64<_SIZE_>::corriger_multiperiodique : there is too many periodic boundaries." << finl;
          exit();
        }
    }

  // Transform marqueur_bord:
  // 1 if the vertex belongs to several periodic boundaries,
  // 0 otherwise
  for (int_t sommet = 0; sommet < nb_som; sommet++)
    {
      // Count the number of bits set to 1 in the marker
      const int marq = marqueur_bord[sommet];
      int n = 0;
      for (int x = 1; x < marq; x = x * 2)
        {
          if (marq & x) // bitwise AND
            n++;
        }
      // Marker is 1 if the vertex belongs to several boundaries
      marqueur_bord[sommet] = (n > 1);
    }
  // Second step: assign elements adjacent to a multi-periodic vertex
  // to the part associated with the renum_som_perio vertex of that vertex.
  const IntTab_t& les_elems = domaine.les_elems();
  const int nb_som_elem = les_elems.dimension_int(1);
  int_t count = 0;
  for (int_t elem = 0; elem < nb_elem; elem++)
    {
      // This int will be -1 if no vertex of the element is periodic,
      // otherwise it is the smallest of the parts associated with the periodic vertices
      int new_part = -1;
      for (int isom = 0; isom < nb_som_elem; isom++)
        {
          const int_t sommet = les_elems(elem, isom);
          if (marqueur_bord[sommet])
            {
              // This vertex belongs to several periodic boundaries
              const int_t renum = renum_som_perio[sommet];
              const int part = partie_associee[renum];
              assert(marqueur_bord[renum]);
              assert(part > -1);
              if (new_part < 0 || new_part > part)
                new_part = part;
            }
        }
      if (new_part >= 0 && new_part != elem_part[elem])
        {
          count++;
          elem_part[elem] = new_part;
        }
    }
  Cerr << "Partitionneur_base_32_64<_SIZE_>::corriger_multiperiodique : " << count << " modified elements" << finl;

  // Add another criteria to fix some problem when building renum_som_perio in parallel:
  // Check every periodic nodes are surrounded by elements on the same parts
  SmallArrOfTID_T<_SIZE_> node(2);
  int_t another_count=0;
  bool err=false;  // never changed later??
  for (int_t som=0; som<nb_som; som++)
    {
      // Periodic nodes:
      node[0]=som;
      node[1]=renum_som_perio[som];
      if (node[0]!=node[1])
        {
          // Loop on each element surrounding the nodes
          IntLists part(2);
          for (int i=0; i<2; i++)
            {
              for (int_t i_elem=0; i_elem<som_elem.get_list_size(node[i]); i_elem++)
                {
                  int_t elem = som_elem(node[i],i_elem);
                  part[i].add_if_not(elem_part[elem]);
                }
            }
          // Check if same number:
          if (part[0].size()!=part[1].size())
            {
              Cerr << "Warning: Not the same number of parts around the periodic nodes " << node[0] << " and " << node[1] << " !" << finl;
              Cerr << "We try to fix:" << finl;
              int smaller = ( part[0].size() < part[1].size() ? 0 : 1);
              int bigger  = 1 - smaller;
              int first_part = part[smaller][0]; // Take arbitrary the first part of the smallest list
              for (int i=0; i<part[bigger].size(); i++)
                {
                  int i_part = part[bigger][i];
                  if (!part[smaller].contient(i_part))
                    {
                      // i_part -> first_part on all elements surrounding the node with more parts:
                      for (int_t i_elem=0; i_elem<som_elem.get_list_size(node[bigger]); i_elem++)
                        {
                          int_t elem = som_elem(node[bigger],i_elem);
                          if (elem_part[elem] == i_part)
                            {
                              elem_part[elem] = first_part;
                              another_count++;
                              Cerr << "Element " << elem << " moved from part " << i_part << " to " << first_part << finl;
                            }
                        }
                    }
                }
            }
          else
            {
              // Check if same parts around the 2 nodes:
              for (int i=0; i<part[0].size(); i++)
                {
                  int i_part = part[0][i];
                  if (!part[1].contient(i_part))
                    {
                      // Implement an algorithm as just above ?
                      Cerr << "Warning: different parts around the periodic nodes " << node[0] << " and " << node[1] << " !" << finl;
                      //for (int j=0;j<2;j++)
                      //   for (int k=0;k<part[j].size();k++)
                      //      Cerr << "node " << j << " part " << part[j][k] << finl;
                      Cerr << "We try to fix:" << finl;
                      // Look for the j_part not in the part[1] list:
                      for (int j=0; j<part[1].size(); j++)
                        {
                          int j_part = part[1][j];
                          if (!part[0].contient(j_part))
                            {
                              // Choice between i_part and j_part:
                              // We take i_part arbitrary so j_part -> i_part
                              for (int_t i_elem=0; i_elem<som_elem.get_list_size(node[1]); i_elem++)
                                {
                                  int_t elem = som_elem(node[1],i_elem);
                                  if (elem_part[elem] == j_part)
                                    {
                                      elem_part[elem] = i_part;
                                      Cerr << "Element " << elem << " moved from part " << j_part << " to " << i_part << finl;
                                      another_count++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  if (err)
    {
      Cerr << "Error in Partitionneur_base_32_64<_SIZE_>::corriger_multiperiodique" << finl;
      Cerr << "It will create possible problems for creation of renum_som_perio array." << finl;
      Cerr << "Contact TRUST support or change the partition options." << finl;
      exit();
    }
  count+=another_count;
  Cerr << "Partitionneur_base_32_64<_SIZE_>::corriger_multiperiodique : plus " << another_count << " another modified elements" << finl;
  return count;
}

/*! @brief Corrects elem_part so that element i is on the same partition elem_part[i] as all elements connected to it in the graph
 *
 *   (elements with indices graph_elements_perio(i, j) for all j).
 *
 * @param graph_elements_perio Graph computed by calculer_graphe_connexions_periodiques.
 * @param som_elem Vertex-element connectivity.
 * @param domaine The domain.
 * @param elem_part For each element, which part it belongs to. Return value: number of elements whose partition was corrected.
 */
template <typename _SIZE_>
typename Partitionneur_base_32_64<_SIZE_>::int_t
Partitionneur_base_32_64<_SIZE_>::corriger_bords_avec_graphe(const Static_Int_Lists_t& graph_elements_perio,
                                                             const Static_Int_Lists_t& som_elem,
                                                             const Domaine_t& domaine,
                                                             BigIntVect_& elem_part)
{
  // Algorithm: loop over all elements in order.
  //  For each element, assign to all linked elements the part to which the current
  //  element belongs. Since the graph is symmetric, if an element has already been
  //  processed, nothing is changed on subsequent passes. One pass is sufficient.
  const Noms& liste_bords_periodiques = domaine.bords_perio();
  const int_t n = graph_elements_perio.get_nb_lists(); //elem_part.size_array();
  //assert(n == graph_elements_perio.get_nb_lists());
  int_t count = 0;
  for (int_t i = 0; i < n; i++)
    {
      const int_t m = graph_elements_perio.get_list_size(i);
      const int part = elem_part[i];
      for (int_t j = 0; j < m; j++)
        {
          const int_t elem2 = graph_elements_perio(i,j);
          if (elem_part[elem2] != part)
            {
              elem_part[elem2] = part;
              count++;
            }
        }
    }
  Cerr << "Partitionneur_base_32_64<_SIZE_>::corriger_bords_avec_graphe : " << count
       << " modified periodic elements" << finl;

  const int_t nb_sommets_reels = domaine.nb_som();
  ArrOfInt_t renum_som_perio(nb_sommets_reels);
  // Initialise the renum_som_perio array
  for (int_t i = 0; i < nb_sommets_reels; i++)
    renum_som_perio[i] = i;
  bool parallel_algo = Process::is_parallel();
  Reordonner_faces_periodiques_32_64<_SIZE_>::renum_som_perio(domaine, renum_som_perio,
                                                              parallel_algo /* no virtual space in sequential */);

  if (liste_bords_periodiques.size() > 1)
    count += corriger_multiperiodique(domaine, renum_som_perio, som_elem, elem_part);
  count += corriger_sommets_bord(domaine, renum_som_perio, som_elem, elem_part);
  return count;
}

/*! @brief Computes the periodic element connectivity graphs and calls corriger_periodique_avec_graphe.
 *
 * (Method to use when the connectivity graph is not yet available;
 *    if the graph is already at hand, call corriger_periodique_avec_graphe directly.)
 *
 */
template <typename _SIZE_>
void Partitionneur_base_32_64<_SIZE_>::corriger_bords_avec_liste(const Domaine_t& dom,
                                                                 const int_t my_offset,
                                                                 BigIntVect_& elem_part)
{
  Cerr << "Correction of the splitting for the periodicity" << finl;
  Static_Int_Lists_t som_elem;
  Cerr << " Construction of the connectivity som_elem" << finl;
  construire_connectivite_som_elem(dom.nb_som_tot(),
                                   dom.les_elems(),
                                   som_elem,
                                   1 /* include virtual elements */);
  Cerr << " Construction of graph connectivity for periodic elements" << finl;
  Static_Int_Lists_t graph_elements_perio;
  calculer_graphe_connexions_periodiques(dom,
                                         som_elem,
                                         my_offset,
                                         graph_elements_perio);
  const int_t count = corriger_bords_avec_graphe(graph_elements_perio,
                                                 som_elem,
                                                 dom,
                                                 elem_part);
  Cerr << "corriger_bords_avec_liste : we have corrected " << count << " elements all in all." << finl;
}



template class Partitionneur_base_32_64<int>;
#if INT_is_64_ == 2
template class Partitionneur_base_32_64<trustIdType>;
#endif


