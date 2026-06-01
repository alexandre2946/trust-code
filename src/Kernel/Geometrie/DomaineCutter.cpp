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

#include <DomaineCutter.h>
#include <ArrOfBit.h>
#include <Domaine.h>
#include <Connectivite_som_elem.h>
#include <SFichierBin.h>
#include <Array_tools.h>
#include <Scatter.h>
#include <TRUSTArrays.h>
#include <Sous_Domaine.h>
#include <Sparskit.h>
#include <Poly_geom_base.h>
#include <Sortie_Brute.h>
#include <TRUSTVect.h>
#include <FichierHDFPar.h>
#include <communications.h>

Implemente_instanciable_32_64(DomaineCutter_32_64,"DomaineCutter",Objet_U);

template<typename _SIZE_>
Sortie& DomaineCutter_32_64<_SIZE_>::printOn(Sortie& os) const
{
  Process::exit("Error : DomaineCutter_32_64<_SIZE_>::printOn should not be used.");
  return os;
}

template<typename _SIZE_>
Entree& DomaineCutter_32_64<_SIZE_>::readOn(Entree& s)
{
  Process::exit("Error : DomaineCutter_32_64<_SIZE_>::readOn should not be used.");
  return s;
}


namespace  // anonymous namespace
{

/*! @brief Creation of the vertex list of sub-domain "partie".
 *
 * This is the set of vertices of elements belonging to this sub-domain.
 *  Only real elements are processed.
 *
 * @param (nb_sommets) number of vertices of the global domain
 * @param (les_elems)
 * @param (elem_part) partitioning array (for each element i of the global domain, elem_part[i] is the number of the sub-domain to which it is assigned)
 * @param (partie) the number of the sub-domain to build
 * @param (liste_sommets) on output: list of vertices of the sub-domain: liste_sommets[i] is the index in domaine_globale of the i-th vertex of the sub-domain. Indices are sorted in ascending order.
 * @param (liste_inverse_sommets) on output: sized to nb_sommets and initialized. liste_inverse_sommet[i] is the index of the vertex in the sub-domain, or -1 if vertex i is not in the sub-domain)
 */
template<typename _SIZE_>
void construire_liste_sommets_sousdomaine(const _SIZE_ nb_sommets,
                                          const IntTab_T<_SIZE_>& les_elems,
                                          const ArrOfInt_T<_SIZE_>& liste_elements,
                                          const int i_part,
                                          const Static_Int_Lists_32_64<_SIZE_> *som_raccord,
                                          SmallArrOfTID_T<_SIZE_>& liste_sommets,
                                          BigArrOfInt_T<_SIZE_>& liste_inverse_sommets)
{
  using int_t = _SIZE_;

  const int nb_elem_part = Process::check_int_overflow(liste_elements.size_array());
  const int_t nb_sommets_par_element = les_elems.dimension(1);

  // Algorithm: iterate over elements; for elements of the part,
  // mark the vertices of the element with a flag.
  // Then iterate over vertices, and those whose flag is set
  // are added to the list.

  // First, count the vertices of the part and
  // fill drapeau_sommet
  ArrOfBit_32_64<_SIZE_> drapeau_sommet(nb_sommets);

  drapeau_sommet = 0;
  // Number of vertices of part "part"
  int nb_sommets_part = 0;
  for (int i_elem = 0; i_elem < nb_elem_part; i_elem++)
    {
      const int_t elem = liste_elements[i_elem];
      for (int j = 0; j < nb_sommets_par_element; j++)
        {
          int_t sommet = les_elems(elem, j);
          if (sommet>-1)
            {
              int bit = drapeau_sommet.testsetbit(sommet);
              // If the flag was not set, this is one more vertex
              if (! bit)
                nb_sommets_part++;
            }
        }
    }
  //vertices to add due to som_raccord
  if (som_raccord)
    for (int s = 0; s < som_raccord->get_nb_lists(); s++)
      for (int_t i = 0; i < som_raccord->get_list_size(s); i++)
        if ((*som_raccord)(s, i) == i_part && !drapeau_sommet.testsetbit(s)) //the vertex is requested by this proc
          nb_sommets_part++; //if we don't already have it, add it

  // Fill liste_sommets and liste_inverse_sommets
  liste_sommets.resize_array(0); // Forget previous values
  liste_sommets.resize_array(nb_sommets_part);

  liste_inverse_sommets.resize_array(0); // Forget previous values
  liste_inverse_sommets.resize_array(nb_sommets,RESIZE_OPTIONS::NOCOPY_NOINIT);
  liste_inverse_sommets = -1;

  int n = 0;
  for (int_t i = 0; i < nb_sommets; i++)
    {
      if (drapeau_sommet[i])
        {
          liste_sommets[n] = i;
          liste_inverse_sommets[i] = n;
          n++;
        }
    }
}

/*! Fill the vertex coordinate array of a part
 * from the vertex coordinates of the complete domain (sommets_glob)
 * and the list of vertices of the part (liste_sommets).
 * sommets_loc is created as follows:
 *   sommets_loc.dimension(0) = liste_sommets.size_array()
 *   sommets_loc.dimension(1) = sommets_glob.dimension(1)
 *
 * Parameter:     sommets_glob
 * Meaning: the coordinates of the vertices of the global domain
 * Parameter:     liste_sommets
 * Meaning: the indices of the vertices of sommets_glob to copy
 *          into sommets_loc.
 */
template<typename _SIZE_>
void remplir_coordsommets_sous_domaine(const DoubleTab_T<_SIZE_>& sommets_glob,
                                       const SmallArrOfTID_T<_SIZE_>& liste_sommets,
                                       DoubleTab& sommets_loc)
{
  using int_t = _SIZE_;
  const int nb_som = liste_sommets.size_array();
  const int dim = sommets_glob.dimension_int(1);
  sommets_loc.resize(nb_som, dim);

  for (int i = 0; i < nb_som; i++)
    {
      int_t indice = liste_sommets[i];
      for (int j = 0; j < dim; j++)
        sommets_loc(i, j) = sommets_glob(indice, j);
    }
}

/*! Build the elems array of elements of part "part".
 * elems(i, j) will be the new number of vertex j of element i
 * in part. The liste_inverse is used to obtain the
 * new vertex numbers.
 * The elements of the local domain are created in ascending order of
 * their index in the global domain.
 */
template<typename _SIZE_>
void construire_elems_sous_domaine(const IntTab_T<_SIZE_>&    elems_domaine_globale,
                                   const ArrOfInt_T<_SIZE_>& liste_elements,
                                   const BigArrOfInt_T<_SIZE_>& liste_inverse_sommets,
                                   IntTab&    elems_domaine_locale,
                                   BigArrOfInt_T<_SIZE_>& liste_inverse_elements)
{
  using int_t = _SIZE_;

  const int_t nb_elem_tot                = elems_domaine_globale.dimension_tot(0);
  const int nb_sommets_par_element = elems_domaine_globale.dimension_int(1);

  liste_inverse_elements.resize(nb_elem_tot,RESIZE_OPTIONS::NOCOPY_NOINIT);
  liste_inverse_elements = -1;

  // First pass: count the number of elements in the part
  const int nb_elem_part = Process::check_int_overflow(liste_elements.size_array());
  elems_domaine_locale.resize(nb_elem_part, nb_sommets_par_element);

  // Second pass: fill the array
  for (int i_elem = 0; i_elem < nb_elem_part; i_elem++)
    {
      int_t elem = liste_elements[i_elem];

      // New element number in the part
      liste_inverse_elements[elem] = i_elem;
      // Copy the element with vertex numbers translated to
      // local numbers in the sub-domain
      for (int i = 0; i < nb_sommets_par_element; i++)
        {
          int_t sommet = elems_domaine_globale(elem, i);
          if (sommet<0)
            {
              assert(sommet == -1); // Polytopes
              elems_domaine_locale(i_elem, i) = -1;
            }
          else
            {
              int new_num = liste_inverse_sommets[sommet];
              assert(new_num >= 0);
              elems_domaine_locale(i_elem, i) = new_num;
            }
        }
    }
}

/*! @brief For a list of "faces" of the global domain, count the number of faces included in part "part" and copy them into the structure
 *
 *  faces_partie by replacing vertex numbers with local numbers
 *  in the sub-domain.
 *  The order of faces is preserved (if a face appears before another in
 *  the list of the complete domain and both are in the sub-part,
 *  their order is preserved). This is important for periodicity
 *  (assumption that there is a correspondence between face i and face i+n/2 of the
 *  periodic boundary).
 *  Note: the condition for a face to be included is "the face belongs
 *   to an element of the neighboring part". The condition "the vertices of the faces are
 *   joint vertices" is NOT sufficient.
 */
template<typename _SIZE_>
void construire_liste_faces_sous_domaine(const ArrOfInt_T<_SIZE_>& elements_voisins,
                                         const BigIntVect_T<_SIZE_>& elem_part,
                                         const int     partie,
                                         const IntTab_T<_SIZE_>& faces_sommets,
                                         const BigArrOfInt_T<_SIZE_>& liste_inverse_sommets,
                                         IntTab& faces_sommets_partie)
{
  using int_t = _SIZE_;

  const int_t nb_faces = faces_sommets.dimension(0);
  const int nb_sommets_par_face = faces_sommets.dimension_int(1);

  assert(elements_voisins.size_array() == nb_faces);
  // List of faces from the faces_sommets array to include in the sub-domain
  ArrOfInt_T<_SIZE_> liste_faces;

  // First pass: look for faces to include
  {
    for (int_t i = 0; i < nb_faces; i++)
      {
        int_t elem = elements_voisins[i];
        const int part = elem_part[elem];
        if (part == partie)
          liste_faces.append_array(i);
      }
  }

  const int nb_faces_part = Process::check_int_overflow(liste_faces.size_array());
  faces_sommets_partie.resize(nb_faces_part, nb_sommets_par_face);

  // Second pass: store the faces and convert vertex numbers
  //  to local numbers in the sub-domain.
  for (int i = 0; i < nb_faces_part; i++)
    {
      const int_t i_face = liste_faces[i];
      for (int j = 0; j < nb_sommets_par_face; j++)
        {
          int_t sommet = faces_sommets(i_face, j);
          if (sommet>-1)
            {
              int new_num = liste_inverse_sommets[sommet];
              assert(new_num >= 0);
              faces_sommets_partie(i, j) = new_num;
            }
          else
            faces_sommets_partie(i, j) = -1;
        }
    }
}

} // end anonymous NS

/*! @brief For each PE listed in the "voisins" array, if a joint with that pe does not yet exist in the domain,
 * add a joint and initialize it with "epaisseur".
 */
template <typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::ajouter_joints(Domaine32& domaine, const ArrOfInt& voisins) const
{
  Joints& joints = domaine.faces_joint();

  const int n = voisins.size_array();
  for (int i = 0; i < n; i++)
    {
      const int pe = voisins[i];
      int j = 0;
      const int nb_joints = joints.size();
      for (; j < nb_joints; j++)
        if (joints[j].PEvoisin() == pe)
          break;
      if (j == nb_joints)
        {
          Cerr << " Adding of a new joint : " << pe << finl;
          Joint& joint = joints.add(Joint());
          joint.nommer(Nom("Joint_")+Nom(pe));
          joint.associer_domaine(domaine);
          joint.affecte_epaisseur(epaisseur_joint_);
          joint.affecte_PEvoisin(pe);
          // These joints will have no common vertices. Set the initialization flag to 1.
          joint.set_joint_item(JOINT_ITEM::SOMMET).set_items_communs();
        }
    }
}


/*! @brief Starting from a list of vertices (liste_sommets_depart), iterate over neighboring elements of these vertices (if epaisseur <= 1),
 *
 *   then the neighbors of these elements (sharing a vertex with the element) if epaisseur <= 2,
 *   then the neighbors of the neighbors if epaisseur <= 3, etc.
 *   Elements belonging to "partie_a_ignorer" are not visited.
 *   The indices of the visited elements are stored in liste_elements_trouves.
 *   This method was written for construire_elements_distants_ssdom()
 */
template <typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::parcourir_epaisseurs_elements(SmallArrOfTID_t liste_sommets_depart, /* by value since we will modify it */
                                                                const int partie_a_ignorer, /* do not traverse elements of this partition */
                                                                SmallArrOfTID_t& liste_elements_trouves) const
{
  const Domaine_t& domaine          = ref_domaine_.valeur();
  const IntTab_t& elements    = domaine.les_elems();
  const BigIntVect_t& elem_part = ref_elem_part_.valeur();
  const int_t nb_sommets = som_elem_.get_nb_lists();
  const int_t nb_elements = elements.dimension_tot(0);
  const int nb_som_elem = elements.dimension_int(1);

  liste_elements_trouves.resize_array(0);

  ArrOfBit_32_64<_SIZE_> sommets_parcourus(nb_sommets);
  ArrOfBit_32_64<_SIZE_> elements_parcourus(nb_elements);
  SmallArrOfTID_t new_liste;

  sommets_parcourus = 0;
  elements_parcourus = 0;

  {
    // Mark the starting vertices and build a list of unique vertices
    const int sz_liste = liste_sommets_depart.size_array();
    for (int i = 0; i < sz_liste; i++)
      {
        const int_t sommet = liste_sommets_depart[i];
        if (sommet>-1)
          if (!sommets_parcourus.testsetbit(sommet))
            new_liste.append_array(sommet);
      }
  }

  // Loop over successive layers of elements around the joint vertices.
  for (int ep = 1; ep <= epaisseur_joint_; ep++)
    {
      liste_sommets_depart = new_liste;
      new_liste.resize_array(0);
      const int sz_liste = liste_sommets_depart.size_array();
      for (int i_sommet = 0; i_sommet < sz_liste; i_sommet++)
        {
          const int_t sommet = liste_sommets_depart[i_sommet];
          // Iterate over elements neighboring this vertex:
          const int nb_elems_voisins = Process::check_int_overflow(som_elem_.get_list_size(sommet));
          for (int i_elem = 0; i_elem < nb_elems_voisins; i_elem++)
            {
              const int_t elem = som_elem_(sommet, i_elem);
              if (elements_parcourus[elem])
                continue; // This element has already been visited
              // Domain number of the neighboring element:
              const int part = elem_part[elem];
              if (part == partie_a_ignorer)
                continue; // Only process elements from other parts

              liste_elements_trouves.append_array(elem);
              elements_parcourus.setbit(elem);

              // Add to the next list of vertices to process the vertices of this element:
              for (int i = 0; i < nb_som_elem; i++)
                {
                  const int_t sommet2 = elements(elem, i);
                  if (sommet2>-1)
                    {
                      if (sommets_parcourus.testsetbit(sommet2) == 0)
                        new_liste.append_array(sommet2);
                    }
                }
            }
        }
    }
}


/*! For the following three functions:
 * Fill the boundaries of the part associated with a processor.
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_faces_bords_ssdom(const BigArrOfInt_T<_SIZE_>& liste_inverse_sommets,
                                                               const int partie,
                                                               Domaine32& domaine_partie) const
{
  const Domaine_t& domaine = ref_domaine_.valeur();
  int i_fr = 0;
  ArrOfInt_t elements_voisins;
  for(const auto& itr: domaine.faces_bord())
    {
      const Frontiere_t& frontiere = itr;
      Frontiere& front_partie = domaine_partie.faces_bord().add(Bord());
      front_partie.nommer(frontiere.le_nom());
      front_partie.associer_domaine(domaine_partie);
      front_partie.faces().typer(frontiere.faces().type_face());
      voisins_bords_.copy_list_to_array(i_fr, elements_voisins);
      construire_liste_faces_sous_domaine(elements_voisins,
                                          ref_elem_part_.valeur(),
                                          partie,
                                          frontiere.faces().les_sommets(),
                                          liste_inverse_sommets,
                                          front_partie.faces().les_sommets());
      i_fr++;
    }
}

/*! @brief Builds the connections of the sub-partition from the connections of the global domain.
 *
 *   Raccord_local_homogene objects are transformed into Raccord_distant_homogene objects.
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_faces_raccords_ssdom(const BigArrOfInt_t& liste_inverse_sommets,
                                                                  const int partie,
                                                                  Domaine32& domaine_partie) const
{
  const Domaine_t& domaine = ref_domaine_.valeur();
  int i_fr = domaine.nb_bords();
  ArrOfInt_t elements_voisins;
  for(const auto& itr: domaine.faces_raccord())
    {
      const Frontiere_t& frontiere = itr;
      Raccord& raccord_partie = domaine_partie.faces_raccord().add(Raccord());
      Nom type_raccord = frontiere.que_suis_je();
      // Strip away a potential "_64" suffix since we are now building a reduced 32b dom.
      type_raccord.prefix("_64");
      if (type_raccord == "Raccord_local_homogene")
        type_raccord = "Raccord_distant_homogene";
      raccord_partie.typer(type_raccord);
      Frontiere& front_partie = raccord_partie.valeur();
      front_partie.nommer(frontiere.le_nom());
      front_partie.associer_domaine(domaine_partie);
      front_partie.faces().typer(frontiere.faces().type_face());
      voisins_bords_.copy_list_to_array(i_fr, elements_voisins);
      construire_liste_faces_sous_domaine(elements_voisins,
                                          ref_elem_part_.valeur(),
                                          partie,
                                          frontiere.faces().les_sommets(),
                                          liste_inverse_sommets,
                                          front_partie.faces().les_sommets());
      i_fr++;
    }
}

template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_frontieres_internes_ssdom(const BigArrOfInt_t& liste_inverse_sommets,
                                                                       const int partie,
                                                                       Domaine32& domaine_partie) const
{
  // Reminder: internal boundaries are "frontiers" inside the domain
  // (for example a zero-thickness plate in the flow)
  const Domaine_t& domaine = ref_domaine_.valeur();
  int i_fr = domaine.nb_bords()+domaine.nb_raccords();
  ArrOfInt_t elements_voisins;
  for(const auto& itr: domaine.bords_int())
    {
      const Frontiere_t& frontiere = itr;
      Frontiere& front_partie = domaine_partie.bords_int().add(Bord_Interne());
      front_partie.nommer(frontiere.le_nom());
      front_partie.associer_domaine(domaine_partie);
      front_partie.faces().typer(frontiere.faces().type_face());
      voisins_bords_.copy_list_to_array(i_fr, elements_voisins);
      construire_liste_faces_sous_domaine(elements_voisins,
                                          ref_elem_part_.valeur(),
                                          partie,
                                          frontiere.faces().les_sommets(),
                                          liste_inverse_sommets,
                                          front_partie.faces().les_sommets());
      i_fr++;
    }
}

template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_groupe_faces_ssdom(const BigArrOfInt_t& liste_inverse_sommets,
                                                                const int partie,
                                                                Domaine32& domaine_partie) const
{
  using Groupe_Faces_t = Groupe_Faces_32_64<_SIZE_>;

  const Domaine_t& domaine = ref_domaine_.valeur();
  Static_Int_Lists_t voisins;
  const int nb_groupe_faces = domaine.nb_groupes_faces();
  ArrOfInt_T<_SIZE_> nb_faces(nb_groupe_faces);
  for (int i = 0; i < nb_groupe_faces; i++)
    nb_faces[i] = domaine.groupe_faces(i).nb_faces();
  voisins.set_list_sizes(nb_faces);
  const int nb_som_face = domaine.frontiere(0).faces().nb_som_faces();
  SmallArrOfTID_t une_face(nb_som_face), elements_voisins;

  const BigIntVect_t& elem_part = ref_elem_part_.valeur();

  for (int grp = 0; grp < nb_groupe_faces; grp++)
    {
      const Groupe_Faces_t& groupe_faces = domaine.groupe_faces(grp);
      Groupe_Faces& groupe_partie = domaine_partie.groupes_faces().add(Groupe_Faces());
      groupe_partie.nommer(groupe_faces.le_nom());
      groupe_partie.associer_domaine(domaine_partie);
      groupe_partie.faces().typer(groupe_faces.faces().type_face());
      const IntTab_t& faces_sommets = groupe_faces.faces().les_sommets();
      IntTab& faces_sommets_partie = groupe_partie.faces().les_sommets();

      ArrOfInt_t liste_faces;
      // First pass: look for faces to include
      const IntTab_t& faces = groupe_faces.les_sommets_des_faces();
      const int_t n = nb_faces[grp];
      for (int_t j = 0; j < n; j++)
        {
          for (int k = 0; k < nb_som_face; k++)
            une_face[k] = faces(j, k);
          find_adjacent_elements(som_elem_, une_face, elements_voisins);

          if (elements_voisins.size_array() == 1 && elem_part[elements_voisins[0]] == partie) liste_faces.append_array(j);

          if (elements_voisins.size_array() == 2)
            if (elem_part[elements_voisins[0]] == partie || elem_part[elements_voisins[1]] == partie)
              liste_faces.append_array(j);
        }

      const int nb_faces_part = Process::check_int_overflow(liste_faces.size_array());
      faces_sommets_partie.resize(nb_faces_part, nb_som_face);

      // Second pass: store the faces and convert vertex numbers
      //  to local numbers in the sub-domain.
      for (int i = 0; i < nb_faces_part; i++)
        {
          const int_t i_face = liste_faces[i];
          for (int j = 0; j < nb_som_face; j++)
            {
              int_t sommet = faces_sommets(i_face, j);
              if (sommet>-1)
                {
                  int new_num = liste_inverse_sommets[sommet];
                  assert(new_num >= 0);
                  faces_sommets_partie(i, j) = new_num;
                }
              else
                faces_sommets_partie(i, j) = -1;
            }
        }
    }
}


/*! @brief Compute and fill domaine_partie.
 *
 * joint(i).set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants()
 *    (list of distant elements).
 *    New joints may be created.
 *  History: method coded in January 2006 by B.Mathieu. The original belief was that
 *   determining distant elements during scatter would be too complex. Eventually a way
 *   was found. As a result this method is no longer used in normal operation. In case of
 *   problems, it can be re-enabled (option 'print_more_info' in Decouper):
 *   in theory it gives exactly the same result as the algorithm in Scatter::calculer_espace_distant_elements.
 *   Note: it does nothing special for periodic boundaries...
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_elements_distants_ssdom(const int     partie,
                                                                     const SmallArrOfTID_t& liste_sommets,
                                                                     const BigArrOfInt_t& liste_inverse_elements,
                                                                     Domaine32& domaine_partie) const
{
  const Domaine_t& domaine          = ref_domaine_.valeur();
  const IntTab_t& elements    = domaine.les_elems();
  const BigIntVect_t& elem_part = ref_elem_part_.valeur();
  const int_t nb_elem = elements.dimension(0);
  const int_t nb_som_elem = elements.dimension(1);

  // First step: look for virtual elements of the sub-part
  // (elements of other parts located within the joint thickness starting
  //  from the joint vertices).
  SmallArrOfTID_t elements_virtuels;
  {
    // Build an array containing all the vertices of all the joints
    SmallArrOfTID_t sommets_joint;

    const int nb_joints = domaine_partie.nb_joints();
    for (int i_joint = 0; i_joint < nb_joints; i_joint++)
      {
        const ArrOfInt& sommets = domaine_partie.joint(i_joint).joint_item(JOINT_ITEM::SOMMET).items_communs();
        const int n = sommets.size_array();
        for (int i = 0; i < n; i++)
          {
            const int sommet_local = sommets[i];
            const int_t sommet_global = liste_sommets[sommet_local];
            sommets_joint.append_array(sommet_global);
          }
      }

    parcourir_epaisseurs_elements(sommets_joint, // Starting vertices
                                  partie, /* ignore elements of the local partition */
                                  elements_virtuels);
  }
  const int nb_elements_virtuels = elements_virtuels.size_array();
  // Build the list of neighboring domain numbers and add missing joints:
  ArrOfInt parties_voisines;
  {
    for (int i = 0; i < nb_elements_virtuels; i++)
      {
        const int_t elem = elements_virtuels[i];
        const int part = elem_part[elem];
        parties_voisines.append_array(part);
      }
    array_trier_retirer_doublons(parties_voisines);
    ajouter_joints(domaine_partie, parties_voisines);
  }

  // For each neighboring part, look for distant elements:
  const int nb_parties_voisines = parties_voisines.size_array();
  SmallArrOfTID_t sommets_depart;

  for (int i_part = 0; i_part < nb_parties_voisines; i_part++)
    {
      const int partie_voisine = parties_voisines[i_part];
      sommets_depart.resize_array(0);
      // List of vertices of the virtual elements of the neighboring part:
      for (int i_elem = 0; i_elem < nb_elements_virtuels; i_elem++)
        {
          const int_t elem = elements_virtuels[i_elem];
          if (elem_part(elem) == partie_voisine)
            {
              for (int i = 0; i < nb_som_elem; i++)
                {
                  const int_t som = elements(elem, i);
                  sommets_depart.append_array(som);
                }
            }
        }

      // Search for elements neighboring this list:
      SmallArrOfTID_t elements_distants;
      parcourir_epaisseurs_elements(sommets_depart,
                                    partie_voisine, /* ignore elements of the neighboring part */
                                    elements_distants);

      // Remove from the list the elements that are not in the local part:
      {
        int count = 0;
        const int n = elements_distants.size_array();
        for (int i = 0; i < n; i++)
          {
            const int_t elem = elements_distants[i];
            if (elem_part[elem] == partie && elem < nb_elem)
              elements_distants[count++] = elem;
          }
        elements_distants.resize_array(count);
      }
      // Translate element indices to local indices on the sub-domain.
      const int n = elements_distants.size_array();
      ArrOfInt elem_dist_loc(n);
      {
        for (int i = 0; i < n; i++)
          {
            const int_t elem = elements_distants[i];
            const int renum = liste_inverse_elements[elem];
            elem_dist_loc[i] = renum;
          }
      }
      // Sort elements in ascending order of indices
      elements_distants.ordonne_array();
      // Store the list of distant elements in the joint
      Joint&     joint = domaine_partie.joint_of_pe(partie_voisine);
      joint.set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants() = elem_dist_loc;
    }
}

/*! @brief Create joints and build vertex lists and arrays for all joints of
 *
 *   part. The vertices of the joint with "PEvoisin" are
 *   the vertices of part "part" that are also a vertex of an
 *   element belonging to PEvoisin.
 *   Joints appear in the list in ascending order of PEs.
 *   Note:  The faces to be created next necessarily use
 *          vertices found here.
 *   Property: The vertices of the joint are sorted in ascending order
 *               of their local number, hence in ascending order of
 *               their global number (see Remplir_Numeros_Sommets)
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_sommets_joints_ssdom(const SmallArrOfTID_t& liste_sommets,
                                                                  const BigArrOfInt_t& liste_inverse_sommets,
                                                                  const int partie,
                                                                  const Static_Int_Lists_t* som_raccord,
                                                                  Domaine32& domaine_partie) const
{
  Joints& les_joints = domaine_partie.faces_joint();

  const int parts = nb_parties_;
  const BigIntVect_t& elem_part = ref_elem_part_.valeur();

  // List of joint vertices (across all parts; a vertex may appear
  // multiple times in the array, at most once per neighboring pe).
  // This is the global vertex number.
  ArrsOfInt_T<_SIZE_> joints_sommets(parts);

  // Algorithm: loop over the vertices of the part.
  // For each vertex of sub-domain "p", look for the PEs to which
  // the neighboring elements of the vertex belong, and, if there is a neighboring PE
  // different from "p", add this vertex to the joints with those PEs.

  const int nb_sommets = liste_sommets.size_array();

  // i = local index of the vertex in the sub-domain
  for (int i_sommet = 0; i_sommet < nb_sommets; i_sommet++)
    {
      // Global vertex number
      int_t sommet = liste_sommets[i_sommet];
      // Loop over the elements neighboring vertex i
      int nb_elements_voisins = Process::check_int_overflow(som_elem_.get_list_size(sommet));
      for (int i = 0; i < nb_elements_voisins; i++)
        {
          int_t elem = som_elem_(sommet, i);
          int PEvoisin = elem_part[elem];
          // Should this vertex be added to the joint vertices with PEvoisin?
          if (PEvoisin != partie)
            joints_sommets[PEvoisin].append_array(sommet);
        }
      //loop over procs connected to the vertex by a connection
      if (som_raccord)
        for (int_t i = 0; i < som_raccord->get_list_size(sommet); i++)
          {
            int part_num = Process::check_int_overflow((*som_raccord)(sommet, i));
            if (part_num != partie)
              joints_sommets[part_num].append_array(sommet);
          }
    }

  // Create joints: in ascending order of the neighboring PE
  for (int PEvoisin = 0; PEvoisin < parts; PEvoisin++)
    {
      ArrOfInt_t& sommets = joints_sommets[PEvoisin];
      if (sommets.size_array() > 0)
        {
          // Sort joint vertices in ascending order of global index
          // (so all lists on neighboring processors will be sorted
          //  in the same way) and remove duplicates
          array_trier_retirer_doublons(sommets);
          const int nb_sommets2 = Process::check_int_overflow(sommets.size_array());

          // Create a new joint
          Joint& joint = les_joints.add(Joint());
          Nom nom_joint("Joint_");
          Nom nom_numero(PEvoisin);
          nom_joint+=nom_numero;
          joint.nommer(nom_joint);
          joint.affecte_epaisseur(epaisseur_joint_);
          joint.associer_domaine(domaine_partie);
          joint.affecte_PEvoisin(PEvoisin);
          ArrOfInt& sommets_locaux = joint.set_joint_item(JOINT_ITEM::SOMMET).set_items_communs();
          sommets_locaux.resize_array(nb_sommets2);
          // Fill the array (convert to local vertex index)
          for (int i = 0; i < nb_sommets2; i++)
            {
              const int_t indice_global = sommets[i];
              const int indice_local = liste_inverse_sommets[indice_global];
              assert(indice_local >= 0);
              sommets_locaux[i] = indice_local;
            }
        }
    }
}

/*! Search for joint faces (faces adjacent to two elements belonging
 * to two different processors).
 * It is assumed that the joint vertices have already been built.
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_faces_joints_ssdom(const int partie,
                                                                const DomaineCutter_Correspondance_t& correspondance,
                                                                Domaine32& domaine_partie) const
{
  const int nb_sommets_ssdom = domaine_partie.nb_som();
  const int nb_elem_ssdom = domaine_partie.nb_elem();
  const SmallArrOfTID_t& liste_sommets = correspondance.liste_sommets_;
  const BigArrOfInt_t& liste_inverse_elements = correspondance.liste_inverse_elements_;

  // First step: identify "joint elements" (elements of "partie"
  //  having a joint vertex), and mark the joint vertices
  SmallArrOfTID_t liste_elements_joint;
  ArrOfBit drapeaux_sommets_joints(nb_sommets_ssdom);
  {
    // For now, use joints as const...
    const Joints& joints_partie = domaine_partie.faces_joint();
    drapeaux_sommets_joints = 0;

    // Markers: for each element of the sub-domain, is it already in the list?
    ArrOfBit drapeaux_elements(nb_elem_ssdom);
    drapeaux_elements = 0;
    // Loop over all joint vertices
    const int nb_joints = joints_partie.size();
    for (int i_joint = 0; i_joint < nb_joints; i_joint++)
      {
        const Joint& joint = joints_partie[i_joint];
        const ArrOfInt& sommets_du_joint = joint.joint_item(JOINT_ITEM::SOMMET).items_communs();
        const int n = sommets_du_joint.size_array();
        for (int i = 0; i < n; i++)
          {
            const int i_sommet_local = sommets_du_joint[i];
            const int_t i_sommet_global = liste_sommets[i_sommet_local];
            const int nb_elem_voisins = Process::check_int_overflow(som_elem_.get_list_size(i_sommet_global));
            // Mark the joint vertex in the sub-domain
            drapeaux_sommets_joints.setbit(i_sommet_local);
            for (int j = 0; j < nb_elem_voisins; j++)
              {
                const int_t i_elem_global = som_elem_(i_sommet_global, j);
                const int i_elem_local = liste_inverse_elements[i_elem_global];
                // Add the element to the joint element list if it is in my part
                // and if it is not yet in the list.
                if (i_elem_local >= 0)
                  if (! drapeaux_elements.testsetbit(i_elem_local))
                    liste_elements_joint.append_array(i_elem_global);
              }
          }
      }
  }

  // Second step: for each joint element, build its faces and check
  // if the face has a neighboring element in another part.
  // If so, it is a joint face.

  // Array of vertex indices of all joint faces
  //  (across all joints, vertex indices in the sub-domain)
  IntTab faces_joints;
  // For each joint face, the number of the neighboring sub-domain through this face
  ArrOfInt faces_pe_voisins;
  // For each processor, number of joint faces with that proc.
  // (initialized to zero)
  ArrOfInt nb_faces_joints(nb_parties_);
  {
    // Indices of face vertices on the reference element
    IntTab faces_element_reference;
    const Elem_geom_base& type_elem = domaine_partie.type_elem().valeur();
    int is_regular = type_elem.get_tab_faces_sommets_locaux(faces_element_reference);
    // Overload for polytope meshes
    if (!is_regular)
      ref_cast(Poly_geom_base,type_elem).get_tab_faces_sommets_locaux(faces_element_reference,0);

    int nb_faces_elem = faces_element_reference.dimension(0);
    const int nb_sommets_par_face = faces_element_reference.dimension(1);
    faces_joints.resize(0, nb_sommets_par_face); // Voir *suite*
    const BigArrOfInt_t& liste_inverse_sommets = correspondance.liste_inverse_sommets_;
    const BigIntVect_t& elem_part = ref_elem_part_.valeur();
    // Vertices of the elements of the global mesh
    const Domaine_t& domaine = ref_domaine_.valeur();
    const IntTab_t& elem_som = domaine.les_elems();
    SmallArrOfTID_t une_face(nb_sommets_par_face);
    SmallArrOfTID_t elements_voisins;

    // Loop over joint elements
    const int nb_elem_joints = liste_elements_joint.size_array();
    for (int i_elem_joint = 0; i_elem_joint < nb_elem_joints; i_elem_joint++)
      {
        // Index of the element in the global domain
        const int_t i_elem_global = liste_elements_joint[i_elem_joint];
        // Loop over the faces of the element
        if (!is_regular)
          {
            int loc_idx = liste_inverse_elements[i_elem_global];
            ref_cast(Poly_geom_base,type_elem).get_tab_faces_sommets_locaux(faces_element_reference, loc_idx);
            nb_faces_elem       = faces_element_reference.dimension(0);
            while ( faces_element_reference(nb_faces_elem-1,0) == -1)
              nb_faces_elem--;
          }
        for (int i_face = 0; i_face < nb_faces_elem; i_face++)
          {
            // Build the face
            int face_ok = 1;
            for (int i = 0; i < nb_sommets_par_face; i++)
              {
                const int i_ref = faces_element_reference(i_face, i);
                if (i_ref<0)
                  une_face[i] = i_ref;
                else
                  {
                    const int_t i_som = elem_som(i_elem_global, i_ref);
                    une_face[i] = i_som;
                    // Index of the vertex in the sub-domain

                    if (i_som>-1)
                      {
                        const int i_som_local = liste_inverse_sommets[i_som];
                        if(i_som_local>0)
                          {
                            // Is the vertex on a joint?
                            if (drapeaux_sommets_joints[i_som_local] == 0)
                              face_ok = 0; // No => this face is not on a joint
                          }
                      }
                  }
              }
            // First test to quickly reject the face
            // if not all its vertices are joint vertices.
            if (face_ok)
              {
                // Search for neighboring elements of this face
                find_adjacent_elements(som_elem_, une_face, elements_voisins);
                // Search for a neighboring element that is not in my part
                const int nb_elem_voisins = elements_voisins.size_array();
                int PEvoisin = -1, i=0;
                for (; i < nb_elem_voisins; i++)
                  {
                    const int_t i_elem_glob = elements_voisins[i];
                    PEvoisin = elem_part[i_elem_glob];
                    if (PEvoisin != partie)
                      break;
                  }
                if (i < nb_elem_voisins)
                  {
                    // I have a joint face, add it to the array
                    faces_pe_voisins.append_array(PEvoisin);
                    const int n = faces_joints.dimension(0);
                    faces_joints.resize(n+1, nb_sommets_par_face);
                    for (int j = 0; j < nb_sommets_par_face; j++)
                      {
                        const int_t i_som_global = une_face[j];
                        if (i_som_global>-1)
                          {
                            const int i_som_local = liste_inverse_sommets[i_som_global];
                            faces_joints(n, j) = i_som_local;
                          }
                        else
                          faces_joints(n, j) = -1;
                      }
                    nb_faces_joints[PEvoisin]++;
                  }
              }
          }
      }
  }

  // Third step: add the faces to the joints
  {
    // The type of joint faces
    const Domaine_t& domaine_globale = ref_domaine_.valeur();
    const Type_Face& type_face_joint = domaine_globale.type_elem()->type_face();
    // Modify the joints:
    Joints& joints_partie = domaine_partie.faces_joint();
    const int nb_joints = joints_partie.size();
    const int nb_faces_joints_tot = faces_joints.dimension(0);
    // *continued*: 2D array even if there are no faces
    const int nb_sommets_par_face = faces_joints.dimension(1);

    for (int i_joint = 0; i_joint < nb_joints; i_joint++)
      {
        // Face type and memory allocation
        Joint& joint = joints_partie[i_joint];
        joint.typer_faces(type_face_joint);
        joint.faces().les_sommets().resize(0, nb_sommets_par_face);
        const int PE_voisin = joint.PEvoisin();
        const int nb_faces = nb_faces_joints[PE_voisin];
        joint.dimensionner(nb_faces);
        // Fill the face array
        IntTab& faces = joint.faces().les_sommets();
        int k = 0;
        for (int i = 0; i < nb_faces_joints_tot; i++)
          {
            const int face_pe = faces_pe_voisins[i];
            if (face_pe == PE_voisin)
              {
                for (int j = 0; j < nb_sommets_par_face; j++)
                  faces(k, j) = faces_joints(i,j);
                k++;
              }
          }
      }
  }
}

/*! @brief cancels all references and empties the arrays
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::reset()
{
  ref_domaine_.reset();
  ref_elem_part_.reset();
  nb_parties_ = -1;
  epaisseur_joint_ = -1;
  som_elem_.reset();
}

namespace // anonymous NS
{

template<typename _SIZE_>
void calculer_listes_elements_sous_domaines(const BigIntVect_T<_SIZE_>& elem_part,
                                            const int nb_parts,
                                            const _SIZE_ nbelem,
                                            Static_Int_Lists_32_64<_SIZE_>& liste_elems_sous_domaines)
{
  using int_t = _SIZE_;

  ArrOfInt_T<_SIZE_> sizes(nb_parts);
  // First pass: counting
  for (int_t i = 0; i < nbelem; i++)
    {
      const int part = elem_part[i];
      sizes[part]++;
    }
  liste_elems_sous_domaines.set_list_sizes(sizes);
  // Second pass: filling
  sizes = 0;
  for (int_t i = 0; i < nbelem; i++)
    {
      const int part = elem_part[i];
      const int_t index = sizes[part]++;
      liste_elems_sous_domaines.set_value(part, index, i);
    }
}

template<typename _SIZE_>
void calculer_elements_voisins_bords(const Domaine_32_64<_SIZE_>& domaine,
                                     const Static_Int_Lists_32_64<_SIZE_>& som_elem,
                                     Static_Int_Lists_32_64<_SIZE_>& voisins,const BigIntVect_T<_SIZE_>& elem_part,
                                     const bool permissif, Noms& bords_a_pb_)
{
  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;

  const int nb_front = domaine.nb_front_Cl();
  ArrOfInt_t nb_faces(nb_front);
  for (int i = 0; i < nb_front; i++)
    nb_faces[i] = domaine.frontiere(i).nb_faces();
  voisins.set_list_sizes(nb_faces);
  const int nb_som_face = domaine.frontiere(0).faces().nb_som_faces();
  SmallArrOfTID_t une_face(nb_som_face);
  SmallArrOfTID_t elems_voisins;

  for (int i = 0; i < nb_front; i++)
    {
      int drap=0;
      const IntTab_T<_SIZE_>& faces = domaine.frontiere(i).les_sommets_des_faces();
      const int_t n = nb_faces[i];
      for (int_t j = 0; j < n; j++)
        {
          for (int k = 0; k < nb_som_face; k++)
            une_face[k] = faces(j, k);
          find_adjacent_elements(som_elem, une_face, elems_voisins);
          const int n_voisins = elems_voisins.size_array();
          if (n_voisins != 1)
            {
              if (drap==0)
                {
                  if (permissif)
                    Cerr << "Message from DomaineCutter.cpp : calculer_elements_voisins_bords\n"
                         << " boundary face "<<  domaine.frontiere(i).le_nom()<<" with " << n_voisins << " neighbors." << finl;
                  else
                    Cerr << "Error in DomaineCutter.cpp : calculer_elements_voisins_bords\n"
                         << " boundary face "<<  domaine.frontiere(i).le_nom()<<" with " << n_voisins << " neighbors." << finl;
                  bords_a_pb_.add( domaine.frontiere(i).le_nom());
                  if (elems_voisins.size_array()==0) Process::exit();
                }
              drap=1;
              if (elem_part[elems_voisins[1]]==1)
                elems_voisins[0]=elems_voisins[1];
              if (!permissif)
                Process::exit();
            }
          voisins.set_value(i, j, elems_voisins[0]);
        }
    }
}

/*! @brief Build the name of the ".
 *
 * Domaines" file for a given proc and a domain. If partie == -1 a single filename is returned. For example
 *           DOM.Zones
 *    instead of
 *           DOM_0001.Zones
 */
void construire_nom_fichier_sous_domaine(const Nom& basename,
                                         const int partie, const int nb_parties_,
                                         const int original_proc,
                                         Nom& fichier)
{
  fichier = basename;

  if (partie > 100000)
    {
      Cerr << "Error while generating filename: nb_parties_ too large" << finl;
      Process::exit();
    }
  char s[30];
  // single file name for all procs (HDF5)
  // the number of parts is still included in the file name
  // (make_PAR.data needs to know how many domaines there are)
  if (partie < 0)
    {
      if (nb_parties_ > 10000)
        snprintf(s, 30, "_p%05d.Zones", (int)nb_parties_);
      else
        snprintf(s, 30, "_p%04d.Zones",(int) nb_parties_);
    }
  else
    {
      if (nb_parties_ > 10000)
        {
          if(original_proc < 0)
            snprintf(s, 30, "_%05d.Zones", (int)partie);
          else
            snprintf(s, 30, "_%05d_%d.Zones", (int)partie, (int)original_proc);
        }
      else
        {
          if(original_proc < 0)
            snprintf(s, 30, "_%04d.Zones", (int)partie);
          else
            snprintf(s, 30, "_%04d_%d.Zones", (int)partie, (int)original_proc);
        }
    }
  fichier += Nom(s);
}

}// end anonymous NS


/*! @brief Prepares the data structures for building sub-domains based on a partitioning provided in elem_part.
 *
 * @param (domaine_global) the domain to partition (must have a domain). A reference to this domain is stored (it must remain valid).
 * @param (elem_part) for each element, the number of the sub-domain to which it belongs, with 0 <= elem_part[i] < nb_parts. Note: a reference to this array is stored, not a local copy. The array must continue to exist until DomaineCutter is no longer used.
 * @param (nb_parts) total number of sub-domains
 * @param (les_faces)
 * @param (epaisseur_joint)
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::initialiser(const Domaine_t&   domaine_global,
                                              const BigIntVect_t& elem_part,
                                              const int     nb_parts,
                                              const int     epaisseur_joint,
                                              const bool permissif)
{
  assert(nb_parts >= 0);
  assert(elem_part.size_array() == domaine_global.nb_elem_tot());
  assert(max_array(elem_part) < nb_parts);
  if (min_array(elem_part)<0)
    {
      Cerr << "Error in DomaineCutter_32_64<_SIZE_>::initialiser" << finl;
      Cerr << "The built partition has a default." << finl;
      Cerr << "Try to change the partition tool or try changing some partitioning options." << finl;
      Cerr << "Contact TRUST support if you are still unsuccessful with your tries." << finl;
      Process::exit();
    }
  reset();

  ref_domaine_ = domaine_global;
  ref_elem_part_ = elem_part;
  nb_parties_ = nb_parts;
  epaisseur_joint_ = epaisseur_joint;

  const IntTab_t& elems = domaine_global.les_elems();
  const int_t nb_som = domaine_global.nb_som_tot();
  construire_connectivite_som_elem(nb_som, elems, som_elem_,
                                   1 /* include virtual elements */);

  Cerr << "Search of neighboring elements of boundary faces" << finl;
  calculer_elements_voisins_bords(domaine_global, som_elem_, voisins_bords_,elem_part,permissif,bords_a_pb_);
  Cerr << "Construction of lists of elements per subdomain" << finl;
  calculer_listes_elements_sous_domaines(elem_part, nb_parts, elems.dimension(0), liste_elems_sous_domaines_);
  //calculer_listes_elements_sous_domaines(elem_part, nb_parts, liste_elems_sous_domaines_);
}

/*! @brief Fills the "correspondance" structure and the "sous_domaine" for part "part".
 *
 * Look for the vertices belonging to the sub-part,
 *   compute the renumbering arrays for vertices and elements between global
 *   and local numbers (correspondance structure), fill the sub_domaine structures
 *   (vertices, elements, boundaries, joints, etc).
 *   Note: sous_domaine must be an empty object (not containing a domain)
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::construire_sous_domaine(const int part, DomaineCutter_Correspondance_t& correspondance,
                                                          Domaine32& sous_domain, const Static_Int_Lists_t *som_raccord) const
{
  using Poly_geom_base_t = Poly_geom_base_32_64<_SIZE_>;

  // The object must be initialized:
  assert(nb_parties_ >= 0);
  // sous_domaine empty
  assert(sous_domain.nb_elem() == 0);
  // Valid part number
  assert(part >= 0 && part < nb_parties_);

  correspondance.partie_ = part;

  const Domaine_t& domain = ref_domaine_.valeur();

  ArrOfInt_t elements_sous_partie;
  liste_elems_sous_domaines_.copy_list_to_array(part, elements_sous_partie);

  // Preparation of sous_domaine
  // sous_domaine.reset(); /* reset does not exist yet... */
  sous_domain.nommer(domain.le_nom());
  if (sub_type(Poly_geom_base_t, domain.type_elem().valeur()))
    {
      const Poly_geom_base_t& dom_as_poly = ref_cast(Poly_geom_base_t,domain.type_elem().valeur());
      dom_as_poly.build_reduced(sous_domain.type_elem(), elements_sous_partie);
    }
  else
    {
      // Funny! We need to build the 32b elem type from its 64b version ...!
      Nom el_nam = domain.type_elem()->que_suis_je();
      el_nam.prefix("_64");
      sous_domain.typer(el_nam);
    }
  sous_domain.type_elem()->associer_domaine(sous_domain);

  // Copy names of periodic boundaries:
  sous_domain.bords_perio() = domain.bords_perio();

  construire_liste_sommets_sousdomaine(domain.nb_som_tot(), domain.les_elems(), elements_sous_partie, part,
                                       som_raccord, correspondance.liste_sommets_ /* write */,
                                       correspondance.liste_inverse_sommets_ /* write */);

  remplir_coordsommets_sous_domaine(domain.coord_sommets(), correspondance.liste_sommets_, sous_domain.les_sommets() /* write */);
  construire_elems_sous_domaine(domain.les_elems(), elements_sous_partie, correspondance.liste_inverse_sommets_, sous_domain.les_elems() /* write */,
                                correspondance.liste_inverse_elements_ /* write */);

  const SmallArrOfTID_t& l_som = correspondance.liste_sommets_;
  const BigArrOfInt_t& l_inv_som = correspondance.liste_inverse_sommets_;
  construire_faces_bords_ssdom(l_inv_som, part, sous_domain);
  construire_faces_raccords_ssdom(l_inv_som, part, sous_domain);
  construire_frontieres_internes_ssdom(l_inv_som, part, sous_domain);
  construire_groupe_faces_ssdom(l_inv_som, part, sous_domain);
  construire_sommets_joints_ssdom(l_som, l_inv_som, part, som_raccord, sous_domain);
  construire_faces_joints_ssdom(part, correspondance, sous_domain);

  //if som_raccord is used (for DecouperMulti), then construire_elements_distants_ssdom()
  //can lead to the creation of empty joints -> it must not be called
  int compute_items_distants = Decouper_t::print_more_infos_ && !som_raccord; // To print NbElemDist informations
  if (compute_items_distants)
    {
      // This sequential algorithm is not used in normal operation, unless
      // testing the parallel algorithm in Scatter.cpp.
      // see CHECK_ALGO_ESPACE_VIRTUEL in Scatter.cpp (Benoit Mathieu) and/or option 'print_more_info' in Decouper
      construire_elements_distants_ssdom(part, correspondance.liste_sommets_, correspondance.liste_inverse_elements_, sous_domain);
    }
  else
    {
      // Initialize empty joints (otherwise assert in debug mode since joints are not initialized)
      for (int ij = 0; ij < sous_domain.nb_joints(); ij++)
        {
          Joint& joint = sous_domain.joint(ij);
          joint.set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants();
        }
    }
  Scatter::trier_les_joints(sous_domain.faces_joint());
}


template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::writeData(const Domaine& sous_domaine, Sortie& os) const
{
  os << sous_domaine;
}

/*! @brief Generation of all sub-domains of the computation and writing to disk of the files basename_000n.
 *
 * Domains for 0 <= n < nb_parties_.
 *   If "sous-domaines" are defined (in the domain.ss_domaines() field),
 *   a file is also generated per sub-domain.
 *
 *   WARNING: this method might change elem_part!!! (if 'reorder' option was specified)
 */
template<typename _SIZE_>
void DomaineCutter_32_64<_SIZE_>::ecrire_domaines(const Nom& basename, const DomainesFileOutputType format,
                                                  const int reorder, const Static_Int_Lists_t* som_raccord)
{
  assert(nb_parties_ >= 0);
  BigIntVect_t& elem_part = ref_elem_part_.valeur();
  const Domaine_t& domaine = ref_domaine_.valeur();
  const int_t nbelem = domaine.nb_elem();
  DomaineCutter_Correspondance_t dc_correspondance;

  // Needed for HDF5 Domaines output:
#ifdef MPI_
  FichierHDFPar fic_hdf;
#else
  FichierHDF fic_hdf;
#endif
  Nom nom_fichier_hdf5(basename);
  nom_fichier_hdf5+=Nom(".Zones");
  nom_fichier_hdf5 = nom_fichier_hdf5.nom_me(nb_parties_, "p", 1);

  // Build temp arrays to eventually reorder the partition numbering
  ArrOfInt ia(nb_parties_+1), ja;

  int nnz=0;
  ia[0]=1;

  //To detect my parts (when running in parallel)
  ArrOfInt myDomaines(nb_parties_);
  myDomaines = 0;
  ArrsOfInt otherProcDomaines(Process::nproc());

  //if some domains are splitted between multiple procs,
  //we assign consecutive indices to each of its fragment
  //(reading the .Zones files during Scatter will be more efficient)
  // Possible values for domaines_index[part]:
  // -2    : means that part is detained by multiple procs but not by me
  // -1    : means that part is detained by a single proc
  // i >=0 : means that my proc detains the i-th fragment of part
  ArrOfInt domaines_index(nb_parties_);
  domaines_index = -2;

  Cerr << "Generation of " << nb_parties_ << " parts:" << finl;
  IntVect EdgeCut(nb_parties_);
  ArrsOfInt Neighbours(nb_parties_);
  // 2 loops if reorder=1
  for (int loop=0; loop<1+reorder; loop++)
    {
      if (reorder)
        {
          Cerr << "====================================" << finl;
          if (loop==0)
            Cerr << "FIRST PASS, analyzing the partition:" << finl;
          else
            Cerr << "SECOND PASS, after reordering the partition:" << finl;
          Cerr << "====================================" << finl;
        }

      if(loop == reorder)
        {
          // check to see which part is shared between multiple processors:
          // 1- everyone sends the number of the parts that belong to them to the master process
          // 2- master process will assign a unique positive number to each fragment of a shared domain
          // 3- if a part is owned by a single process, it is indicated with the index -1
          // 4- the master process scatters the indices to all the proc
          ArrsOfInt domaines_indices(Process::nproc());

          for(int_t i=0; i < nbelem; i++)
            {
              const int part = elem_part[i];
              myDomaines[part] = 1;
            }

          for(int p=0; p<Process::nproc(); p++)
            {
              if(p==0)
                otherProcDomaines[p] = myDomaines;
              else
                {
                  otherProcDomaines[p].resize_array(nb_parties_);
                  otherProcDomaines[p] = 0;
                }
            }
          if(Process::je_suis_maitre())
            {
              for(int p=0; p<Process::nproc(); p++)
                {
                  domaines_indices[p].resize_array(nb_parties_);
                  domaines_indices[p] = -2;
                  if(p!=0)
                    recevoir(otherProcDomaines[p], p, 0, p+2001);
                }

              for(int part=0; part<nb_parties_; part++)
                {
                  int s = 0;
                  for(int proc=0; proc < Process::nproc(); proc++)
                    {
                      if(otherProcDomaines[proc][part])
                        domaines_indices[proc][part] = s++;
                    }

                  if(s<=1)
                    {
                      if(s==0)   //empty part: master process will write it
                        myDomaines[part] = 1;

                      //part is detained by a single proc
                      for(int proc=0; proc < Process::nproc(); proc++)
                        domaines_indices[proc][part] = -1;
                    }
                }

              for(int p=0; p<Process::nproc(); p++)
                {
                  if(p==0)
                    domaines_index = domaines_indices[p];
                  else
                    envoyer(domaines_indices[p], 0, p, p+2002);
                }
            }
          else
            {
              envoyer(myDomaines, Process::me(), 0, Process::me()+2001);
              recevoir(domaines_index, 0, Process::me(), Process::me()+2002);
            }

          if (format == DomainesFileOutputType::HDF5_SINGLE)  // create HDF5 file only once!
            {
              fic_hdf.create(nom_fichier_hdf5);
              if(Process::is_parallel())
                {
                  // creating datasets
                  Noms dataset_names;
                  if(Process::je_suis_maitre())
                    {
                      for(int part=0; part<nb_parties_; part++)
                        {
                          if(domaines_index[part] == -1)
                            {
                              std::string dname = "/zone_"  + std::to_string(part);
                              Nom dataset_name(dname);
                              dataset_names.add(dataset_name);
                            }
                          else
                            {
                              for(int proc=0; proc < Process::nproc(); proc++)
                                {
                                  if(domaines_indices[proc][part] >=0)
                                    {
                                      std::string dname = "/zone_"  + std::to_string(part) + "_" + std::to_string(domaines_indices[proc][part]);
                                      Nom dataset_name(dname);
                                      dataset_names.add(dataset_name);
                                    }
                                }
                            }
                        }
                    }

                  // estimation of an upper bound of the datasets' size
                  int ipart = 0;
                  while(!myDomaines[ipart]) ipart++;
                  Domaine32 dom_tmp;
                  construire_sous_domaine(ipart, dc_correspondance, dom_tmp);
                  Sortie_Brute os_tmp;
                  writeData(dom_tmp, os_tmp);
                  double sz_ = (double)os_tmp.get_size();
                  sz_ *= 1.5;
                  sz_ = Process::mp_max(sz_);
                  envoyer_broadcast(dataset_names,0);
                  long sz_l = lround(sz_);
                  fic_hdf.create_datasets(dataset_names, sz_l);
                }
            }
        }
      for (int i_part = 0; i_part < nb_parties_; i_part++)
        {
          if(!myDomaines[i_part])
            continue;

          assert(domaines_index[i_part] > -2);
          Cerr << " Construction of part number " << i_part << finl;
          if(domaines_index[i_part] >= 0)
            Cerr << "This part is shared between multiple processors" << finl;

          Domaine32 sous_domaine;
          construire_sous_domaine(i_part, dc_correspondance, sous_domaine, som_raccord);
          // Print some information...
          {
            const Joints& joints = sous_domaine.faces_joint();
            const int nb_joints = joints.size();
            Cerr << "  Number of nodes    : " << sous_domaine.nb_som() << finl;
            Cerr << "  Number of elements : " << sous_domaine.nb_elem() << finl;
            Cerr << "  Number of joints   : " << nb_joints << finl;
            //            char s[200];
            int nbfaces_total=0;
            int nbelemdist_total=0;
            int nbsom_total=0;
            Neighbours[i_part].resize_array(0);
            Neighbours[i_part].resize_array(nb_joints);
            for (int i = 0; i < nb_joints; i++)
              {
                const Joint& joint = joints[i];
                const int pe = joint.PEvoisin();
                Neighbours[i_part][i] = pe;
                const int nbsom = joint.joint_item(JOINT_ITEM::SOMMET).items_communs().size_array();
                nbsom_total+=nbsom;
                const int nbfaces = joint.nb_faces();
                nbfaces_total+=nbfaces;
                const int nbelemdist = joint.joint_item(JOINT_ITEM::ELEMENT).items_distants().size_array();
                nbelemdist_total+=nbelemdist;
                Cerr<< "  Joint "<<i<<" PeVoisin "<<pe<<" NbSommets "<<nbsom<<" NbFaces "<<nbfaces;
                if (Decouper_t::print_more_infos_) Cerr <<" NbElemDist "<<nbelemdist;
                Cerr<<finl;
              }
            Cerr<<"               Total:    NbSommets "<<nbsom_total<<" NbFaces "<<nbfaces_total;
            if (Decouper_t::print_more_infos_) Cerr<<" NbElemDist "<<nbelemdist_total;
            EdgeCut(i_part)=nbfaces_total;
            Cerr<<finl;

          }
          if (reorder && loop==0)
            {
              const Joints& joints = sous_domaine.faces_joint();
              const int nb_joints = joints.size();
              // Resize ja:
              ja.resize_array(nnz+nb_joints+1);
              // Fortran numbering:
              ja[nnz]=i_part+1;
              nnz++;
              ia[i_part+1]=ia[i_part]+1;
              for (int i = 0; i < nb_joints; i++)
                {
                  const int pe = joints[i].PEvoisin();
                  ja[nnz] = pe+1;
                  nnz++;
                  ia[i_part+1]++;
                }
            }
          else
            {
              // Write .Zones file(s):
              if (format == DomainesFileOutputType::BINARY_MULTIPLE)
                {
                  Nom nom_fichier(basename);
                  nom_fichier+=Nom(".Zones");
                  //nom_fichier = nom_fichier.nom_me(i_part);
                  construire_nom_fichier_sous_domaine(basename,i_part, nb_parties_, domaines_index[i_part], nom_fichier);
                  Cerr << "Writing part " << i_part << " into the "
                       << (format == DomainesFileOutputType::BINARY_MULTIPLE ? "binary" : "ascii")
                       << " file " << nom_fichier << finl;
                  SFichier os;
                  os.set_bin(1);

                  //
                  // Even for big computations, Zones files can always be written in 32b.
                  //
                  os.set_64b(false);

                  const int ok = os.ouvrir(nom_fichier);
                  if (!ok)
                    Process::exit("DomaineCutter_32_64<_SIZE_>::ecrire_domaines : Error while opening file!");
                  writeData(sous_domaine, os);
                }
              else if (format == DomainesFileOutputType::HDF5_SINGLE)
                {
                  Sortie_Brute os_hdf;
                  writeData(sous_domaine, os_hdf);

                  std::string dname = "/zone_" + std::to_string(i_part);
                  if(domaines_index[i_part] >=0)
                    dname += std::string("_") + std::to_string(domaines_index[i_part]);
                  Nom datasetname(dname);
                  if(Process::is_parallel())
                    fic_hdf.fill_dataset(datasetname, os_hdf);
                  else
                    fic_hdf.create_and_fill_dataset_SW(datasetname, os_hdf);
                }
              else
                {
                  Cerr << "DomaineCutter_32_64<_SIZE_>::ecrire_domaines : Unsupported output file type!" << finl;
                  Process::exit(1);
                }
            }

          // Write sub-domain .ssz files
          const LIST(OBS_PTR(Sous_Domaine_t)) & liste_sous_domaines = domaine.ss_domaines();
          const int nb_sous_domaines = liste_sous_domaines.size();
          if (nb_sous_domaines>0)
            {
              Cerr << " Writing of files .ssz ..." << finl;
              // One file per sub-domain, each file contains all parts...
              // If writing the first part, erase the file
              IOS_OPEN_MODE mode;
              if (i_part == 0)
                mode = ios::out; // First part, overwrite the file
              else
                mode = (ios::out | ios::app); // Append mode

              // Loop over sub-domains
              const BigArrOfInt_t& liste_inverse_elements = dc_correspondance.liste_inverse_elements_;

              for (int i_sous_domaine = 0; i_sous_domaine < nb_sous_domaines; i_sous_domaine++)
                {
                  // Indices of the sub-domain elements that are in the sub-domain:
                  const Sous_Domaine_t& sous_dom = liste_sous_domaines[i_sous_domaine];
                  ArrOfInt elements;
                  {
                    for (int_t i = 0; i < sous_dom.nb_elem_tot(); i++)
                      {
                        const int_t i_elem_global = sous_dom[i];
                        const int i_elem_local = liste_inverse_elements[i_elem_global];
                        if (i_elem_local >= 0)
                          elements.append_array(i_elem_local);
                      }
                  }
                  Nom nom_fichier(sous_dom.le_nom() + Nom(".ssz"));
                  Cerr << " Subarea " << i_sous_domaine
                       << "  file_name " << nom_fichier
                       << "  nb_elements in this part " << elements.size_array()
                       << finl;
                  // The sub-domain file is always in ascii
                  SFichier os;
                  const int ok = os.ouvrir(nom_fichier, mode);
                  if (!ok)
                    {
                      Cerr << "DomaineCutter_32_64<_SIZE_>::ecrire_domaines :\n Error while opening file"
                           << nom_fichier << finl;
                      exit();
                    }
                  os << elements;
                }
            }
        }
      if (reorder && loop==0)
        {
          // Reduce the bandwith of a the matrix connectivity between parts:
          // ToDo force not to change partition 0!
          ArrOfInt riord(nb_parties_+1);
          ArrOfInt levels(nb_parties_+1);
          ArrOfInt mask(nb_parties_+1);
          int maskval = 1;
          for (int i_part=0 ; i_part<nb_parties_+1; i_part++)
            mask[i_part] = maskval;
          int init=1;
          int nlev;
          F77NAME(PERPHN)(&nb_parties_, ja.addr(), ia.addr(), &init, mask.addr(), &maskval, &nlev, riord.addr(), levels.addr());

          // Renumber the elem_part array:
          ArrOfInt renum(nb_parties_);
          for (int i_part=0; i_part < nb_parties_; i_part++)
            renum[riord[i_part]-1]=i_part;
          int_t size=elem_part.size_reelle();
          for (int_t i=0; i<size; i++)
            elem_part[i]=renum[elem_part[i]];
          elem_part.echange_espace_virtuel();
          // Rebuild the liste_elems_sous_domaines_
          liste_elems_sous_domaines_.reset();
          calculer_listes_elements_sous_domaines(elem_part, nb_parties_, nbelem, liste_elems_sous_domaines_);
        }

      /*
      // Print statistics exactly as Fluent:
      // http://combust.hit.edu.cn:8080/fluent/Fluent60_help/html/ug/node984.htm
      Cerr << ">> Partitions:" << finl;
      Cerr << "
               P   Cells I-Cells Cell Ratio  Faces I-Faces Face Ratio Neighbors" << finl;

               0     134      10      0.075    217      10      0.046         1
               1     137      19      0.139    222      19      0.086         2
               2     134      19      0.142    218      19      0.087         2
               3     137      10      0.073    223      10      0.045         1
      Cerr << "------" << finl;
      Cerr << "Partition count             = " << nb_parties_ << finl;
      Cerr << "Cell variation              = (134 - 138)
      Cerr << "Mean cell variation         = (  -1.1% -    1.1%)
      Cerr << "Intercell variation         = (10 - 19)
      Cerr << "Intercell ratio variation   = (   7.3% -   14.2%);
      Cerr << "Global intercell ratio      =   10.7%
      Cerr << "Face variation              = (217 - 223)
      Cerr << "Interface variation         = (10 - 19)
      Cerr << "Interface ratio variation   = (   4.5% -    8.7%)
      Cerr << "Global interface ratio      =    3.4%
      Cerr << "Neighbor variation          = (1 - 2) */
    }

  // if my part is shared with other procs, then my whole part may contain joints that I don't have
  // calling mp_sum wouldn't be correct though, because the shared joints would be counted several times
  // so we need to gather neighbours of each part
  if (Decouper_t::print_more_infos_)   // involves communication: do it only if requested
    {
      if(Process::je_suis_maitre())
        {
          for(int proc=1; proc<Process::nproc(); proc++)
            {
              for(int i_part=0; i_part<nb_parties_; i_part++)
                {
                  if(otherProcDomaines[proc][i_part])
                    {
                      ArrOfInt tmp_neighbours;
                      recevoir(tmp_neighbours, proc, 0, proc+2003);
                      for(int i=0; i<tmp_neighbours.size_array(); i++)
                        Neighbours[i_part].append_array(tmp_neighbours[i]);
                    }
                }
              ArrOfInt tmp_edge_cut(nb_parties_);
              tmp_edge_cut = 0;
              recevoir(tmp_edge_cut, proc, 0, proc+2008);

              for(int i_part=0; i_part<nb_parties_; i_part++)
                EdgeCut(i_part) += tmp_edge_cut[i_part];
            }

          Cout << "\nQuality of partitioning --------------------------------------------" << finl;
          int total_edge_cut = 0;
          for (int i_part=0; i_part<nb_parties_; i_part++)
            total_edge_cut+=EdgeCut(i_part);
          Cout << "Total number of edge-cut (faces shared by processes) : " << total_edge_cut << finl;
          if (total_edge_cut>0)
            {
              double mean_edgecut_domaine = total_edge_cut / nb_parties_;
              int max_edgecut_domaine = local_min_vect(EdgeCut);
              int min_edgecut_domaine = local_max_vect(EdgeCut);

              double load_imbalance = double(max_edgecut_domaine / mean_edgecut_domaine);
              Cout << "Number of edge-cut per Domaine (min/mean/max) : " << min_edgecut_domaine  << " / "
                   << (int) (mean_edgecut_domaine) << " / " << max_edgecut_domaine  << " Load imbalance: " << load_imbalance
                   << "\n" << finl;
            }

          int mean_neighbours = 0;
          int min_neighbours = nb_parties_;
          int max_neighbours = 0;
          for (int i_part = 0; i_part < nb_parties_; i_part++)
            {
              array_trier_retirer_doublons(Neighbours[i_part]);
              int nb_neighbours = Neighbours[i_part].size_array();
              mean_neighbours += nb_neighbours;
              if(nb_neighbours < min_neighbours)   min_neighbours = nb_neighbours;
              if(nb_neighbours > max_neighbours)   max_neighbours = nb_neighbours;
            }

          mean_neighbours/=nb_parties_;
          if (mean_neighbours>0)
            {
              double load_imbalance = double(max_neighbours / mean_neighbours);
              Cout << "Number of neighbours per Domaine (min/mean/max) : " << min_neighbours << " / "
                   << (int) (mean_neighbours) << " / " <<max_neighbours << " Load imbalance: " << load_imbalance
                   << "\n" << finl;
            }
        }
      else
        {
          for(int i_part=0; i_part < nb_parties_; i_part++)
            if(myDomaines[i_part])
              envoyer(Neighbours[i_part], Process::me(), 0, Process::me()+2003);
          envoyer(EdgeCut, Process::me(), 0, Process::me()+2008);
        }
    }

  if (format == DomainesFileOutputType::HDF5_SINGLE)
    fic_hdf.close();
}

template class DomaineCutter_32_64<int>;
#if INT_is_64_ == 2
template class DomaineCutter_32_64<trustIdType>;
#endif
