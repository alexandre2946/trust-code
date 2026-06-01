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

#include <Connex_components.h>
#include <communications.h>
#include <TRUSTTab.h>
#include <ArrOfBit.h>

/*! @brief Computes the connected sets by faces of non-"marked" elements (elements are connected to each other via a symmetric graph
 *
 *   passing through the faces).
 *   A connected domain portion has a unique number 0 <= i < N
 *   and is delimited either by a boundary, or by a neighbor element "marked"
 *   by num_compo[elem] = -1.
 *   This method is sequential (may be called on a single processor)
 *
 * @param (elem_faces)
 * @param (faces_elem)
 * @param (num_compo)
 */
int search_connex_components_local(const IntTab& elem_faces, const IntTab& faces_elem, IntVect& num_compo)
{
  const int nbelem = num_compo.size_totale();
  const int nb_voisins = elem_faces.dimension(1);
  assert(elem_faces.dimension_tot(0) == nbelem);
  {
    int i;
    for (i = 0; i < nbelem; i++)
      if (num_compo[i] != -1)
        num_compo[i] = -2;
  }
  int start_element = 0;
  int num_compo_courant = 0;
  ArrOfInt liste_elems;

  ArrOfInt tmp_liste;

  do
    {
      // Find the next element not yet assigned to a connected component
      while (start_element < nbelem && num_compo[start_element] >= -1)
        start_element++;
      if (start_element == nbelem)
        break;
      // Search for the elements of the connected component starting from this element
      liste_elems.resize_array(1);
      liste_elems[0] = start_element;
      num_compo[start_element] = num_compo_courant;
      while (liste_elems.size_array() > 0)
        {
          tmp_liste.resize_array(0);
          const int liste_elems_size = liste_elems.size_array();
          for (int i_elem = 0; i_elem < liste_elems_size; i_elem++)
            {
              const int elem = liste_elems[i_elem];
              // Add the unassigned neighbors of this element to the list to
              // be processed in the next step
              for (int j = 0; j < nb_voisins; j++)
                {
                  const int face = elem_faces(elem, j);
                  const int voisin = faces_elem(face, 0) + faces_elem(face, 1) - elem;
                  if (voisin >= 0)
                    {
                      const int num = num_compo[voisin];
                      if (num == -2)
                        {
                          num_compo[voisin] = num_compo_courant;
                          tmp_liste.append_array(voisin);
                        }
                    }
                }
            }
          liste_elems = tmp_liste;
        }
      num_compo_courant++;
    }
  while (1);
  // Returns the number of local connected components found
  return num_compo_courant;
}

/*! @brief Searches for the connected components of a local (non-distributed across processors) non-symmetric graph.
 *
 * @param (graph)
 * @param (connex_components)
 */
int compute_graph_connex_components(const IntTab& graph, ArrOfInt& connex_components)
{
  // connex_components must already have the correct size on entry!
  const int nb_sommets = connex_components.size_array();

  // renum_data defines linked lists of "vertex" numbers belonging to
  //  the same connected component.
  // renum_data(i,0) = number of the first "vertex" in the list to which i belongs
  // renum_data(i,1) = number of the next "vertex" in the list
  IntTab renum_data(nb_sommets, 2);
  // At the start, each vertex is alone in a list:
  int i_sommet;
  for (i_sommet = 0; i_sommet < nb_sommets; i_sommet++)
    {
      renum_data(i_sommet, 0) = i_sommet;
      renum_data(i_sommet, 1) = -1; // end of list
    }
  const int nbcouples = graph.dimension(0);
  for (int i_couple = 0; i_couple < nbcouples; i_couple++)
    {
      const int compo1 = graph(i_couple, 0); // the smaller one
      const int compo2 = graph(i_couple, 1); // the larger one
      assert(compo1 != compo2);
      // If the two components are already in the same list,
      // do nothing.
      if (renum_data(compo1, 0) == renum_data(compo2, 0))
        continue;
      // Merge list1 containing compo1 and list2 containing compo2:
      // 1) find the end of the first list
      int fin_liste1 = compo1;
      for (;;)
        {
          const int next = renum_data(fin_liste1, 1);
          if (next < 0)
            break;
          fin_liste1 = next;
        }
      // 2) append list2 at the end of list1:
      const int debut_liste2 = renum_data(compo2, 0);
      renum_data(fin_liste1, 1) = debut_liste2;
      // 2) update the beginning of list for list2:
      i_sommet = debut_liste2;
      const int debut_liste1 = renum_data(compo1, 0);
      do
        {
          renum_data(i_sommet, 0) = debut_liste1;
          i_sommet = renum_data(i_sommet, 1);
        }
      while (i_sommet >= 0);
    }

  // Create a contiguous numbering for the components:
  // Next number to assign
  int count = 0;
  connex_components = -1;
  for (i_sommet = 0; i_sommet < nb_sommets; i_sommet++)
    {
      if (connex_components[i_sommet] < 0)
        {
          // vertex not yet processed
          // Assign a new number to all vertices of the connected component
          // to which i_sommet belongs:
          for (int i = renum_data(i_sommet, 0); i >= 0; i = renum_data(i, 1))
            connex_components[i] = count;
          // New number for the next component
          count++;
        }
    }
  // Return the number of connected components found
  return count;
}

/*! @brief Searches for the connected components of a set of elements distributed across all processors.
 *
 * This method is parallel and must be called at the
 *   same time on all processors.
 *
 * @param (num_compo)
 * @param (nb_local_components)
 */
int compute_global_connex_components(IntVect& num_compo, int nb_local_components)
{
  const int nbelem = num_compo.size();
  const int nbelem_tot = num_compo.size_totale();
  //int i;

  // Transform local connected component indices into a global index
  // (a shift is added to the global indices using mppartial_sum())
  const int decalage = static_cast<int>(Process::mppartial_sum(nb_local_components)); // compo number are never huge
  const int nb_total_components = static_cast<int>(Process::mp_sum(nb_local_components));
  for (int i = 0; i < nbelem_tot; i++)
    if (num_compo[i] >= 0)
      num_compo[i] += decalage;

  // To find correspondences between a local component number and a
  // number of the same component on the neighboring processor, we create a copy of
  // the num_compo array on which we perform an echange_espace_virtuel(). Thus,
  // in the virtual entries of the array, num_compo holds the number of the
  // local component and copie_compo holds the number of that same component on
  // the processor that owns the element. These two numbers therefore designate
  // the same connected component.
  IntVect copie_compo(num_compo);
  copie_compo.echange_espace_virtuel();

  // Search for equivalences between local component numbers and
  // neighboring component numbers. We build a graph whose
  // edges connect equivalent components.
  // Marker array for equivalences already found.
  // Dimensions = nb local components * nb total components
  //  (to avoid counting the same component more than once).
  ArrOfBit markers(nb_local_components * nb_total_components);
  markers = 0;
  // Correspondence table between local and remote connected components
  IntTab graph;

  int graph_size = 0;
  // Iterate over virtual elements only
  for (int i = nbelem; i < nbelem_tot; i++)
    {
      int compo = num_compo[i];
      if (compo < 0)
        continue;
      int compo2 = copie_compo[i];
      // Index of the pair compo2/compo in the markers array
      // The num_compo array must contain only local components:
      assert(compo >= decalage && compo - decalage < nb_local_components);
      // compo2 is necessarily a remote component.
      assert(compo2 < decalage || compo2 - decalage >= nb_local_components);
      const int index = (compo - decalage) * nb_total_components + compo2;
      if (!markers.testsetbit(index))
        {
          graph.resize(graph_size+1, 2);
          // Put the smaller component number in column 0:
          if (compo2 < compo)
            {
              int tmp = compo;
              compo = compo2;
              compo2 = tmp;
            }
          graph(graph_size, 0) = compo;
          graph(graph_size, 1) = compo2;
          graph_size++;
        }
    }

  ArrOfInt renum;
  if (Process::je_suis_maitre())
    {
      // Receive graph portions from other processors
      IntTab tmp;
      const int nproc = Process::nproc();
      int pe;
      for (pe = 1; pe < nproc; pe++)
        {
          recevoir(tmp, pe, 54 /* tag */);
          const int n2 = tmp.dimension(0);
          graph.resize(graph_size + n2, 2);
          for (int i = 0; i < n2; i++)
            {
              graph(graph_size, 0) = tmp(i, 0);
              graph(graph_size, 1) = tmp(i, 1);
              graph_size++;
            }
        }
      // Compute the connected components of the graph
      renum.resize_array(nb_total_components);
      const int n = compute_graph_connex_components(graph, renum);
      Process::Journal() << "compute_global_connex_components: nb_components=" << n << finl;
    }
  else
    {
      // Send the local graph to processor 0
      envoyer(graph, 0, 54 /* tag */);
    }

  // Receive the connected components
  envoyer_broadcast(renum, 0 /* source processor */);

  // Renumber the components in num_compo
  for (int i = 0; i < nbelem_tot; i++)
    {
      const int x = num_compo[i];
      if (x >= 0)
        {
          const int new_x = renum[x];
          num_compo[i] = new_x;
        }
    }
  // Verification: if we do a virtual space exchange,
  //  this should not change the connected component numbers!

  int nb_components = 0;
  // All processors hold the same renum array, so all compute
  //  the same maximum!
  if (renum.size_array() > 0)
    nb_components = max_array(renum) + 1;
  return nb_components;
}

