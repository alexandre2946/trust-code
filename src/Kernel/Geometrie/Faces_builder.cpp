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

#include <Connectivite_som_elem.h>
#include <EcrFicCollecteBin.h>
#include <Elem_geom_base.h>
#include <Poly_geom_base.h>
#include <communications.h>
#include <NettoieNoeuds.h>
#include <Faces_builder.h>
#include <Domaine.h>
#include <Scatter.h>
#include <stdio.h>
#include <vector>
#include <array>
#include <map>

template <typename _SIZE_>
Faces_builder_32_64<_SIZE_>::Faces_builder_32_64() :
  les_elements_ptr_(0),
  connectivite_som_elem_ptr_(0),
  is_polyedre_(-1)
{
}

template <typename _SIZE_>
void Faces_builder_32_64<_SIZE_>::reset()
{
  les_elements_ptr_ = 0;
  connectivite_som_elem_ptr_ = 0;
  faces_element_reference_old_.reset();
  ref_domaine_.reset();
  faces_sommets_.reset();
  face_elem_.reset();
}

/*! @brief From the description of the domain elements and boundaries (borders, connections, face groups, and joints):
 *
 *   Fills the following structures:
 *   - for each domain boundary: fixer_num_premiere_face
 *   - les_faces.faces_sommets (real faces)
 *   - les_faces.faces_voisins (real faces)
 *   - elem_faces              (for the real faces of real elements)
 *        (elem_faces is initialised with size nb_elem_reels x nb_faces_par_elem)
 *   - joints.items_communs(FACE)
 *
 * @param domaine The domain whose faces are being built.
 * @param connect_som_elem Vertex-to-element connectivity.
 * @param les_faces The faces object to fill.
 * @param elem_faces The element-to-face connectivity array to fill.
 */
template <typename _SIZE_>
void Faces_builder_32_64<_SIZE_>::creer_faces_reeles(Domaine_t& domaine,
                                                     const Static_Int_Lists_t& connect_som_elem,
                                                     Faces_t&   les_faces,
                                                     IntTab_t& elem_faces)
{
  les_elements_ptr_ = & domaine.les_elems();

  connectivite_som_elem_ptr_ = & connect_som_elem;
  // The connectivity must include virtual vertices
  assert(connect_som_elem.get_nb_lists() == domaine.nb_som_tot());

  // Fill the reference-element face table

  is_polyedre_=0;
  if (sub_type(Poly_geom_base,domaine.type_elem().valeur()))
    {
      is_polyedre_=1;
    }
  else
    domaine.type_elem()->get_tab_faces_sommets_locaux(faces_element_reference_old_);
  // Array of size (nb_faces, nb_vertices_per_face),
  // giving for each face the indices of its vertices in the domain.
  // Vertex ordering follows the reference element, for the neighboring
  // element of the face with the smallest index.
  IntTab_t& faces_sommets = les_faces.les_sommets();

  // Array of size (nb_faces, 2) containing for each face
  // the indices of the two neighboring elements. If "i_face" has only one neighbor,
  // faces_voisins_(i_face, 1) = -1;
  IntTab_t& faces_voisins = les_faces.voisins();

  // Initialise references used in check_erreur_faces
  faces_sommets_ = faces_sommets;
  face_elem_ = faces_voisins;
  ref_domaine_ = domaine;

  // Element-to-face array:
  //  dimension(0) = number of elements,
  //  dimension(1) = number of faces per element
  //  elem_faces(i,j) = index of face j of element i in the
  //                    faces_sommets and faces_voisins arrays
  //   (element faces are in the order given by faces_element_reference)
  //  appropriate remote and virtual spaces for elements
  const int_t nb_elements          = les_elements().dimension(0);
  const int nb_faces_par_element = faces_element_reference(0).dimension(0);
  elem_faces.resize(nb_elements, nb_faces_par_element);
  elem_faces = -1;

  const int nb_sommets_par_face = faces_element_reference(0).dimension(1);
  // Each face is added with resize(n+1,...), so smart_resize is used:
  // Compute the theoretical number of faces:
  const int_t nb_faces_front = domaine.nb_faces_frontiere() + domaine.nb_faces_joint();
  int_t nb_faces_prevision = (nb_elements * nb_faces_par_element + nb_faces_front) / 2;
  if (is_polyedre_)
    {
      // all faces are already known....
      const Poly_geom_base_t& poly=ref_cast(Poly_geom_base_t,ref_domaine_->type_elem().valeur());
      nb_faces_prevision=(poly.get_somme_nb_faces_elem()+ nb_faces_front) / 2;;
    }
  // Pre-allocate memory for the expected number of faces to avoid repeated
  // reallocations (see set_smart_resize)

  faces_sommets.resize(nb_faces_prevision, nb_sommets_par_face);
  faces_sommets.resize(0, nb_sommets_par_face);

  faces_voisins.resize(nb_faces_prevision, 2);
  faces_voisins.resize(0, 2);

  // ******** Boundary processing **********
  //  note: "num_premiere_face" is initialised for boundaries here!

  // Create boundary faces
  {
    Bords_t& bords = domaine.faces_bord();
    const int n = bords.size();
    for (int i = 0; i < n; i++)
      {
        Frontiere_t& frontiere = bords[i];

        creer_faces_frontiere(1, /* one neighboring element per face */
                              frontiere,
                              faces_sommets,
                              faces_voisins,
                              elem_faces);
      }
  }
// Connections (Raccords)
  {
    Raccords_t& raccords = domaine.faces_raccord();
    const int n = raccords.size();
    for (int i = 0; i < n; i++)
      {
        Frontiere_t& frontiere = raccords[i].valeur();
        creer_faces_frontiere(1, /* one neighboring element per face */
                              frontiere,
                              faces_sommets,
                              faces_voisins,
                              elem_faces);
      }
  }

// Internal boundary faces
  {
    Bords_Internes_t& faces_int = domaine.bords_int();
    const int n = faces_int.size();
    for (int i = 0; i < n; i++)
      {
        Frontiere_t& frontiere = faces_int[i];
        creer_faces_frontiere(2, /* two neighboring elements per face */
                              frontiere,
                              faces_sommets,
                              faces_voisins,
                              elem_faces);
      }

    // Duplicate internal faces: for each face that has two neighbors,
    // create a second identical face with the second neighbor,
    // clear the second neighbor of the original face, and update
    // the neighbor face of the second neighbor:
    if (n > 0)
      {
        Cerr << "Faces_builder_32_64<_SIZE_>::creer_faces_reeles not coded for the internal faces of boundary" << finl;
        Process::exit();
        // To be done based on the old version of domaine2... and needs testing!
      }
  }

// Joint faces
  {
    Joints_t& joints = domaine.faces_joint();
    const int n = joints.size();
    for (int i = 0; i < n; i++)
      {
        Frontiere_t& frontiere = joints[i];
        creer_faces_frontiere(2, /* neighboring elements per face */
                              frontiere,
                              faces_sommets,
                              faces_voisins,
                              elem_faces);
        // Fill items_communs(FACE)
        // Joint faces are in the same order locally and on the neighboring domain.
        Joint_t& joint = joints[i];
        ArrOfInt_t& indices_faces =
          joint.set_joint_item(JOINT_ITEM::FACE).set_items_communs();
        const int_t nb_faces  = joint.nb_faces();
        indices_faces.resize_array(nb_faces);
        const int_t num_premiere_face = joint.num_premiere_face();
        for (int_t i2 = 0; i2 < nb_faces; i2++)
          indices_faces[i2] = num_premiere_face + i2;
      }
  }

// *********************************************
// Internal faces

  creer_faces_internes(faces_sommets,
                       elem_faces,
                       faces_voisins);


// Face group identification
  {
    Groupes_Faces_t& groupes_faces = domaine.groupes_faces();
    const int n = groupes_faces.size();
    for (int i = 0; i < n; i++)
      {
        Groupe_Faces_t& groupe_faces = groupes_faces[i];
        identification_groupe_faces(groupe_faces,
                                    elem_faces);
      }
  }
// *********************************************
// Done: verify that the actual number of faces matches the predicted number
  if (faces_sommets.dimension(0) != nb_faces_prevision)
    {
      Cerr << "Error in Faces_builder_32_64<_SIZE_>::creer_faces_reeles:\n"
           << " number of faces does not match predicted number of faces.\n"
           << " (problem with faces_bords_internes ?)" << finl;
      Process::exit();
    }

// Reset the smart_resize attribute of the faces_sommets and faces_voisins arrays.


// Reset class attributes
  reset();
}

/*! @brief Helper method for creer_faces_frontiere and creer_faces_internes.
 *
 * If the list is non-empty on at least one processor, prints an error message and calls exit().
 *
 * @param message Error message to display.
 * @param liste_faces List of faces that caused the error.
 */
template <typename _SIZE_>
void Faces_builder_32_64<_SIZE_>::check_erreur_faces(const char * message,
                                                     const ArrOfInt_t& liste_faces) const
{
  const int nmax = 100;
  int_t n = liste_faces.size_array();
  if (n > 0)
    {
      Cerr << "==========================" << finl;
      Cerr << "Error!" << finl << message
           << "\nSee log file of this PE for detailed info."
           << finl;
      Sortie& J = Process::Journal();
      J <<  "Error in Faces_builder_32_64<_SIZE_>::creer_faces_*\n"
        << message << finl;
      if (n > nmax)
        {
          J << "Too many faces to display (" << n << ") display only " << nmax << " first faces" << finl;
          n = nmax;
        }
      int_t i;
      J << "Display format:\n"
        << " facenumber = face index in faces_sommet array\n"
        << " som1..som4 = node index\n"
        << " elem1 elem2 = neighbouring element number\n"
        << "facenumber som1 (x1 y1 z1) som2 (x2 y2 z2) [som3 (x3 y3 z3)...] elem1 elem2" << finl;
      char s[1000];
      const DoubleTab_t& coord = ref_domaine_->coord_sommets();
      const IntTab_t&     faces = faces_sommets_.valeur();
      const IntTab_t&     face_elem = face_elem_.valeur();
      const int dim = Objet_U::dimension;
      const int_t nb_som_faces = faces.dimension(1);
      for (i = 0; i < n; i++)
        {
          char *sptr = s;
          const int_t iface = liste_faces[i];
          sptr += snprintf(sptr, 100, "%4ld ",(long) iface);
          for (int j = 0; j < nb_som_faces; j++)
            {
              const int_t isom = faces(iface,j);
              sptr += snprintf(sptr, 100, "%5ld(", (long)isom);
              for (int k = 0; k < dim; k++)
                if (isom!=-1)
                  sptr += snprintf(sptr, 100, "%10.6f", coord(isom, k));
              sptr += snprintf(sptr, 100, ")");
            }
          sptr += snprintf(sptr, 100, "%4ld %4ld", (long)face_elem(iface,0),(long) face_elem(iface,1));
          J << s << finl;
        }
      NettoieNoeuds_t::verifie_noeuds(ref_domaine_.valeur());
      Process::exit();
    }
}

/*! @brief Adds a real face to faces_sommets and faces_voisins.
 *
 */
template <typename _SIZE_>
_SIZE_ Faces_builder_32_64<_SIZE_>::ajouter_une_face(const SmallArrOfTID_t& une_face,
                                                     const int_t elem0,
                                                     const int_t elem1,
                                                     IntTab_t& faces_sommets,
                                                     IntTab_t& faces_voisins)
{
  int i;
  const int_t num_new_face        = faces_sommets.dimension(0);
  const int nb_sommets_par_face = (int)faces_sommets.dimension(1);
  const int_t new_size = num_new_face + 1;

  assert(une_face.size_array() == nb_sommets_par_face);
  faces_sommets.resize(new_size, nb_sommets_par_face);
  for (i = 0; i < nb_sommets_par_face; i++)
    faces_sommets(num_new_face, i) = une_face[i];

  faces_voisins.resize(new_size, 2);
  faces_voisins(num_new_face, 0) = elem0;
  faces_voisins(num_new_face, 1) = elem1;

  return num_new_face;
}

template <typename _SIZE_>
int Faces_builder_32_64<_SIZE_>::chercher_face_element(const IntTab_t&    elem_som,
                                                       const IntTab& faces_element_ref,
                                                       const SmallArrOfTID_t& une_face,
                                                       const int_t    elem)
{
  const int nb_faces_element = (int)faces_element_ref.dimension(0);
  const int nb_sommets_par_face = (int)faces_element_ref.dimension(1);

  int i_face, i_som2, i_som;
  for (i_face = 0; i_face < nb_faces_element; i_face++)
    {
      for (i_som = 0; i_som < nb_sommets_par_face; i_som++)
        {
          const int sommet_elem_ref = faces_element_ref(i_face, i_som);
          int_t sommet_domaine ;
          if (sommet_elem_ref==-1)
            sommet_domaine=-1;
          else
            sommet_domaine = elem_som(elem, sommet_elem_ref);
          for (i_som2 = 0; i_som2 < nb_sommets_par_face; i_som2++)
            if (une_face[i_som2] == sommet_domaine) // if vertex found, stop
              break;
          if (i_som2 == nb_sommets_par_face) // if vertex not found, stop
            break;
        }
      if (i_som == nb_sommets_par_face) // if all vertices have been found, stop
        break;
    }
  if (i_face == nb_faces_element) // if face not found
    return -1;
  else
    return i_face;
}

template <typename _SIZE_>
const IntTab& Faces_builder_32_64<_SIZE_>::faces_element_reference(int_t elem) const
{
  if (is_polyedre_==1)
    {
      const Poly_geom_base_t& poly =ref_cast(Poly_geom_base_t,ref_domaine_->type_elem().valeur());
      IntTab& elem_ref_mod=ref_cast_non_const(IntTab,faces_element_reference_old_);
      poly.get_tab_faces_sommets_locaux(elem_ref_mod,elem);

      //abort();
      //return faces_element_reference(0);
    }
  return faces_element_reference_old_;
}


/*! @brief Helper method: assumes "une_face" contains the vertex indices of a face of the element with index "elem" in the domain.
 *
 *   Searches for the number of this face on the reference element.
 *   If the vertices do not correspond to any face of the element, returns -1.
 *
 */
template <typename _SIZE_>
int Faces_builder_32_64<_SIZE_>::chercher_face_element(const SmallArrOfTID_t& une_face,
                                                       const int_t     elem) const
{
  const IntTab_t& elem_som                = les_elements();
  const IntTab& faces_element_ref       = faces_element_reference(elem);
  int i_face = chercher_face_element(elem_som, faces_element_ref, une_face, elem);
  return i_face;
}

/*! @brief Inserts the faces of the given boundary into the three arrays, after the faces already present in faces_sommets.
 *
 *   Fills:
 *    frontiere.num_premiere_face
 *   Completes:
 *    faces_sommets
 *    elem_faces
 *    faces_voisins
 *
 */
template <typename _SIZE_>
void Faces_builder_32_64<_SIZE_>::creer_faces_frontiere(const int_t nb_voisins_attendus,
                                                        Frontiere_t&   frontiere,
                                                        IntTab_t& faces_sommets,
                                                        IntTab_t& faces_voisins,
                                                        IntTab_t& elem_faces) const
{
  assert(nb_voisins_attendus == 1 || nb_voisins_attendus == 2);

  const Static_Int_Lists_t& som_elem   = connectivite_som_elem();
  const int  nb_sommets_par_face  = faces_element_reference(0).dimension(0) ? faces_element_reference(0).dimension(1) : 3;
  const int_t   num_premiere_face    = faces_sommets.dimension(0);
  const int_t   nb_elem_reels        = elem_faces.dimension(0);
  frontiere.fixer_num_premiere_face(num_premiere_face);

  const Faces_t&   faces_frontiere  = frontiere.faces();
  const IntTab_t& sommets_faces_fr = faces_frontiere.les_sommets();
  const int_t   nb_faces         = faces_frontiere.nb_faces();
  SmallArrOfTID_t une_face(nb_sommets_par_face);
  SmallArrOfTID_t  voisins;

  ArrOfInt_t liste_faces_erreur0;

  ArrOfInt_t liste_faces_erreur1;

  ArrOfInt_t liste_faces_erreur2;

  ArrOfInt_t liste_faces_erreur3;

  constexpr bool STOP_FIRST_ERR = false; // set this to true in Debug to stop gdb at the right place.

  int i_face;
  for (i_face = 0; i_face < nb_faces; i_face++)
    {
      {
        int nb_sommets_par_face_fr= (int)sommets_faces_fr.dimension(1);
        for (int i = 0; i < std::min(nb_sommets_par_face, nb_sommets_par_face_fr); i++)
          une_face[i] = sommets_faces_fr(i_face, i);
        for (int i = std::min(nb_sommets_par_face, nb_sommets_par_face_fr); i < nb_sommets_par_face; i++)
          une_face[i] = -1;
      }
      // What are the neighboring elements of this face?
      find_adjacent_elements(som_elem, une_face, voisins);
      const int_t nb_voisins = voisins.size_array();
      const int_t elem0 = (nb_voisins > 0) ? voisins[0] : -1;
      const int_t elem1 = (nb_voisins > 1) ? voisins[1] : -1;
      const int_t indice_face =
        ajouter_une_face(une_face, elem0, elem1, faces_sommets, faces_voisins);

      switch(nb_voisins)
        {
        case 0:
          {
            // Error: the face has no neighbor
            liste_faces_erreur0.append_array(indice_face);
            if(STOP_FIRST_ERR) Process::exit("A least one face has no neighbor!");
            break;
          }
        case 1:
        case 2:
          {
            if (nb_voisins_attendus == nb_voisins)
              {
                int i_voisin;
                for (i_voisin = 0; i_voisin < nb_voisins; i_voisin++)
                  {
                    const int_t elem = voisins[i_voisin];
                    // What is the face of the element?
                    const int i_face_elem = chercher_face_element(une_face, elem);
                    if (i_face_elem >= 0)
                      {
                        // If it is a real element, associate the face
                        if (elem < nb_elem_reels)
                          {
                            if (elem_faces(elem, i_face_elem) < 0)
                              elem_faces(elem, i_face_elem) = indice_face;
                            else
                              {
                                // Error: this face already exists (in this or another boundary)
                                liste_faces_erreur3.append_array(indice_face);
                                if(STOP_FIRST_ERR) Process::exit("A face already exists! Was found twice!");
                              }
                          }
                      }
                    else
                      {
                        // Error: the face does not belong to the element.
                        liste_faces_erreur0.append_array(indice_face);
                        if(STOP_FIRST_ERR) Process::exit("A face does not belong to any element!");
                      }
                  }
              }
            else
              {
                // Error: unexpected number of neighbors.
                liste_faces_erreur1.append_array(indice_face);
                if(STOP_FIRST_ERR) Process::exit("A face has an unexpected number of neighbors!");
              }
            break;
          }
        default:
          // Error: more than two neighbors, which should not happen.
          liste_faces_erreur2.append_array(indice_face);
          if(STOP_FIRST_ERR) Process::exit("A face has more than 2 neighbors!");
        }
    }
  Nom msg;
  msg = "Boundary \"";
  msg += frontiere.le_nom();
  msg += "\" contains faces which do not belong to any element.";
  check_erreur_faces(msg, liste_faces_erreur0);

  msg = "Boundary \"";
  msg += frontiere.le_nom();
  msg += "\" contains faces that belong to ";
  msg += Nom(3-nb_voisins_attendus);
  msg += " elements.\n";
  switch(nb_voisins_attendus)
    {
    case 1:
      msg += "These faces should have only 1 neighbouring element.";
      break;
    case 2:
      msg += "These faces should have 2 neighbouring elements.";
      break;
    default:
      msg = "Internal error.";
    }
  if (sub_type(Joint, frontiere))
    {
      // Two possible error sources: the joint faces are incorrect,
      // or the domain does not contain virtual elements (at a minimum the domain
      // must contain the virtual elements neighboring the joint faces).
      msg += "(Error in a Joint object: internal error in the mesh splitter or scatter ? )\n";
    }
  check_erreur_faces(msg, liste_faces_erreur1);

  msg = "Boundary \"";
  msg += frontiere.le_nom();
  msg += "\" contains faces that belong to more than 2 elements.\n";
  check_erreur_faces(msg, liste_faces_erreur2);

  msg = "Boundary \"";
  msg += frontiere.le_nom();
  msg += "\" contains faces that already exist in another boundary or in this one.\n";
  check_erreur_faces(msg, liste_faces_erreur3);
}

/*! @brief Construction of the internal faces of the domain (faces with two neighbors that are not "faces_bord_internes").
 *
 *   Joint faces have already been created.
 *
 */
template <typename _SIZE_>
void Faces_builder_32_64<_SIZE_>::creer_faces_internes(IntTab_t& faces_sommets,
                                                       IntTab_t& elem_faces,
                                                       IntTab_t& faces_voisins) const
{
  const IntTab_t& elem_som             = les_elements();
  const Static_Int_Lists_t& som_elem   = connectivite_som_elem();
  //  const IntTab_t & faces_elem_ref       = faces_element_reference();
  const int_t   nb_elem              = elem_som.dimension(0);
  const int   nb_faces_par_element = faces_element_reference(0).dimension(0);
  const int   nb_sommets_par_face  = nb_faces_par_element ? faces_element_reference(0).dimension(1) : 3;

  // Temporary array storing the vertex indices of the face being processed
  SmallArrOfTID_t une_face(nb_sommets_par_face);
  // Temporary array (list of neighboring elements of a face)
  SmallArrOfTID_t voisins;

  // List of faces with only one neighbor not listed in boundary faces (errors):
  ArrOfInt_t liste_faces_frontiere_non_declarees;

  ArrOfInt_t liste_faces_joint_non_declarees;

  // List of faces with a connectivity error (more than
  // two neighboring elements, or connection to vertices that are
  // not on any face of the element:
  ArrOfInt_t liste_faces_erreurs_connectivite;

  constexpr bool STOP_FIRST_ERR = false; // set this to true in Debug to stop gdb at the right place.

  // Loop over elements
  int_t i_elem;
  for (i_elem = 0; i_elem < nb_elem; i_elem++)
    {
      int i_face;
      // Loop over the faces of the element
      for (i_face = 0; i_face < nb_faces_par_element; i_face++)
        {

          // Index of this face in the faces_sommets array.
          // It is -1 if the face has not yet been created.
          int_t indice_face = elem_faces(i_elem, i_face);

          // Compute the vertex indices of the face in the domain:
          int i;
          // Note: this call must stay here...
          const IntTab& faces_elem_ref       = faces_element_reference(i_elem);

          for (i = 0; i < nb_sommets_par_face; i++)
            {
              // index of the vertex on the reference element
              const int i_som_ref = faces_elem_ref(i_face, i);
              // index of the vertex in the domain
              if (i_som_ref==-1)
                une_face[i] = -1;
              else
                {
                  const int_t i_som = elem_som(i_elem, i_som_ref);
                  une_face[i] = i_som;
                }
            }
          if (une_face[0]==-1)
            {
              // dummy face, do nothing
              elem_faces(i_elem, i_face) = -1;
            }
          else
            {
              // Search for neighboring elements of this face.
              // The "voisins" array is sorted in ascending order.
              find_adjacent_elements(som_elem, une_face, voisins);

              const int_t nb_voisins = voisins.size_array();
              assert (nb_voisins > 0); // There should be at least i_elem !!! (or else we have a face made of -1);

              if (nb_voisins == 1)   // ***** The face has 1 neighbor ********
                {

                  assert(voisins[0] == i_elem); // The neighboring element must be i_elem
                  // A face with only one neighboring element is a boundary face.
                  if (indice_face >= 0)
                    {
                      // Ok, this is normal; boundary faces have already been processed
                    }
                  else
                    {
                      // Error: the face does not yet exist. It should have been
                      // created from the boundaries (creer_faces_frontiere)
                      indice_face = ajouter_une_face(une_face, i_elem, -1,
                                                     faces_sommets, faces_voisins);
                      liste_faces_frontiere_non_declarees.append_array(indice_face);
                      if(STOP_FIRST_ERR) Process::exit("Non declared face!");
                    }

                }
              else if (nb_voisins == 2)     // ***** The face has 2 neighbors ********
                {

                  const int_t elem0 = voisins[0];
                  const int_t elem1 = voisins[1];
                  assert(elem0 < elem1);
                  if (indice_face >= 0)
                    {
                      // The face has already been created.
                    }
                  else
                    {
                      // The face does not yet exist.
                      if (elem0 == i_elem)
                        {
                          // Neighbors are sorted: elem0 < elem1
                          // so this is the first time this face is encountered in the
                          // element loop.
                          indice_face = ajouter_une_face(une_face, elem0, elem1,
                                                         faces_sommets, faces_voisins);

                          // Where is this face on the neighboring element?
                          const int i_face_elem1 = chercher_face_element(une_face, elem1);
                          if (i_face_elem1 >= 0)
                            {
                              if (elem1 < nb_elem) // Is the neighboring element real?
                                elem_faces(elem1, i_face_elem1) = indice_face;
                            }
                          else
                            {
                              // Error: the face vertices belong to elem1
                              // but are not on any face of that element. Mesh connectivity error.
                              liste_faces_erreurs_connectivite.append_array(indice_face);
                              if(STOP_FIRST_ERR) Process::exit("Connectivity issue with face!");
                            }
                          if (elem1 >= nb_elem)
                            {
                              // Error: the neighbor is a virtual element; this face
                              // should be in the joint faces and thus already created.
                              liste_faces_joint_non_declarees.append_array(indice_face);
                              if(STOP_FIRST_ERR) Process::exit("Pb with face: its neighbor is virtual! Should not happen here.");
                            }
                        }
                      else
                        {
                          assert(elem1 == i_elem);
                          indice_face = ajouter_une_face(une_face, elem0, elem1,
                                                         faces_sommets, faces_voisins);
                          // We should have already created this face since it is a neighbor of elem0,
                          // which has already been processed (smaller index). If we reach here,
                          // the vertices of "une_face" belong to elem0 but are not on any face
                          // of that element. This is a connectivity error.
                          liste_faces_erreurs_connectivite.append_array(indice_face);
                          if(STOP_FIRST_ERR) Process::exit("Pb with face: connectivity error.");
                        }
                    }

                }
              else                        // ***** The face has > 2 neighbours ********
                {
                  if (indice_face < 0)
                    {
                      const int_t elem0 = voisins[0];
                      const int_t elem1 = voisins[1];
                      indice_face = ajouter_une_face(une_face, elem0, elem1,
                                                     faces_sommets, faces_voisins);
                    }
                  liste_faces_erreurs_connectivite.append_array(indice_face);
                  if(STOP_FIRST_ERR) Process::exit("Pb with face: connectivity error 2.");
                }

              // If the face did not exist, it has been created and its index stored in indice_face.
              // Otherwise, the index of the existing face has been found.
              assert(indice_face >= 0);
              elem_faces(i_elem, i_face) = indice_face; /* WRITE elem_faces */
            }
        }
    }

  // Error handling:
  {
    const char * const msg1 = "We found faces which belong to one element/cell only and are not declared in any boundary ! You forgot to define at least one boundary in your mesh. Fix your mesh.\n";
    const char * const msg2 = "Joint faces are incomplete: internal error in the mesh splitter\n";
    const char * const msg3 = "Connectivity error in the mesh elements. Possible errors:\n- one face of one element belongs to more than 2 elements\n- two element have at least 3 common nodes but these nodes are not faces of these elements\n";
    check_erreur_faces(msg1, liste_faces_frontiere_non_declarees);
    check_erreur_faces(msg2, liste_faces_joint_non_declarees);
    check_erreur_faces(msg3, liste_faces_erreurs_connectivite);
  }
}

/*! @brief Identification of the face groups specified in the domain.
 *
 *   Fills the indices_faces array of a specific face group.
 *
 */
template <typename _SIZE_>
void Faces_builder_32_64<_SIZE_>::identification_groupe_faces(Groupe_Faces_t& groupe_faces,
                                                              const IntTab_t& elem_faces) const
{
  const Static_Int_Lists_t& som_elem   = connectivite_som_elem();
  const int   nb_sommets_par_face  = faces_element_reference(0).dimension(0) ? faces_element_reference(0).dimension(1) : 3;

  const Faces_t&   faces_specifiees  = groupe_faces.faces();
  const IntTab_t& sommets_faces_fr = faces_specifiees.les_sommets();
  const int_t   nb_faces         = faces_specifiees.nb_faces();
  ArrOfInt_t& indices_faces = groupe_faces.get_indices_faces();
  indices_faces.resize_array(nb_faces);

  SmallArrOfTID_t une_face(nb_sommets_par_face);
  SmallArrOfTID_t voisins;

  ArrOfInt_t liste_faces_erreur0;

  ArrOfInt_t liste_faces_erreur1;


  for (int i_face = 0; i_face < nb_faces; i_face++)
    {
      {
        int nb_sommets_par_face_fr= (int)sommets_faces_fr.dimension(1);
        for (int i = 0; i < std::min(nb_sommets_par_face, nb_sommets_par_face_fr); i++)
          une_face[i] = sommets_faces_fr(i_face, i);
        for (int i = std::min(nb_sommets_par_face, nb_sommets_par_face_fr); i < nb_sommets_par_face; i++)
          une_face[i] = -1;
      }
      // What are the neighboring elements of this face?
      find_adjacent_elements(som_elem, une_face, voisins);
      const int_t nb_voisins = voisins.size_array();

      switch(nb_voisins)
        {
        case 0:
          {
            // Error: the face has no neighbor
            liste_faces_erreur0.append_array(i_face);
            break;
          }
        case 1:
        case 2:
          {
            const int_t elem = voisins[0];
            // Which face of the element is it?
            const int i_face_elem = chercher_face_element(une_face, elem);

            if (i_face_elem >= 0)
              // What is the index of the face
              indices_faces[i_face] = elem_faces(elem,i_face_elem);
            break;
          }
        default:
          // Error: more than two neighbors, which should not happen.
          liste_faces_erreur1.append_array(i_face);
        }
    }

  Nom msg;
  msg = "Group of Faces \"";
  msg += groupe_faces.le_nom();
  msg += "\" contains faces which do not belong to any element or not virtual element.";
  check_erreur_faces(msg, liste_faces_erreur0);

  msg = "Group of Faces \"";
  msg += groupe_faces.le_nom();
  msg += "\" contains faces that belong to more than 2 elements.\n";
  check_erreur_faces(msg, liste_faces_erreur1);
}

template class Faces_builder_32_64<int>;
//#if INT_is_64_ == 2
template class Faces_builder_32_64<trustIdType>;
//#endif
