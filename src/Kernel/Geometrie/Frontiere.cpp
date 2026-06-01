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

#include <MD_Vector_tools.h>
#include <Frontiere.h>
#include <TRUSTTab.h>

Implemente_base_32_64(Frontiere_32_64,"Frontiere",Objet_U);
// XD bord_base objet_lecture bord_base INHERITS_BRACE Basic class for block sides. Block sides that are neither edges
// XD_CONT nor connectors are not specified. The duplicate nodes of two blocks in contact are automatically recognized
// XD_CONT and deleted.

/*! @brief Reads the specifications of a boundary from an input stream.
 *
 *     Reads:
 *        the name
 *        the faces
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 */
template <typename _SIZE_>
Entree& Frontiere_32_64<_SIZE_>::readOn(Entree& is)
{
  is >> nom;
  return is >> les_faces;
}


/*! @brief Writes the boundary to an output stream.
 *
 * Writes:
 *       the name of the boundary
 *       the faces
 *
 * @param (Sortie& os) an output stream
 * @return (Sortie&) the modified output stream
 */
template <typename _SIZE_>
Sortie& Frontiere_32_64<_SIZE_>::printOn(Sortie& os) const
{
  os << nom << finl;
  return os << les_faces;
}

/*! @brief Associates the boundary to the domain it belongs to.
 *
 * @param (Domaine& un_domaine) the domain to associate with the boundary
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::associer_domaine(const Domaine_t& un_domaine)
{
  le_dom=un_domaine;
  les_faces.associer_domaine(un_domaine);
}

/*! @brief Gives a name to the boundary.
 *
 * @param (Nom& name) the name to give to the boundary
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::nommer(const Nom& name)
{
  nom=name;
}

/*! @brief Adds one or more faces to the boundary; the face(s) are specified by an array
 *
 *     containing the vertex indices.
 *
 * @param (IntTab& sommets) array containing the vertex indices of the faces to add
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::ajouter_faces(const IntTab_t& sommets)
{
  les_faces.ajouter(sommets);
}

/*! @brief Sets the type of the boundary faces.
 *
 * @param (Motcle& typ) the geometric type of the faces
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::typer_faces(const Motcle& typ)
{
  les_faces.typer(typ);
}

/*! @brief Sets the type of the boundary faces.
 *
 * @param (Type_Face& typ) the geometric type of the faces
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::typer_faces(const Type_Face& typ)
{
  les_faces.typer(typ);
}

/*! @brief Returns the vertices of the boundary faces.
 *
 * @return (IntTab&) the array containing the vertex indices of the boundary faces
 */
template <typename _SIZE_>
typename Frontiere_32_64<_SIZE_>::IntTab_t& Frontiere_32_64<_SIZE_>::les_sommets_des_faces()
{
  return les_faces.les_sommets();
}
/*! @brief Returns the vertices of the boundary faces (const version).
 *
 * @return (IntTab&) the array containing the vertex indices of the boundary faces
 */
template <typename _SIZE_>
const typename Frontiere_32_64<_SIZE_>::IntTab_t& Frontiere_32_64<_SIZE_>::les_sommets_des_faces() const
{
  return les_faces.les_sommets();
}

/*! @brief Renumbers the nodes (vertices) of the faces.
 *
 * The node with index k becomes the node with index Les_Nums[k].
 *
 * @param (IntVect& Les_Nums) the renumbering vector: new_vertex[i] = Les_Nums[old_vertex[i]]
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::renum(const IntVect_t& Les_Nums)
{
  IntTab_t& les_sommets=faces().les_sommets();
  int_t nb_som = les_sommets.dimension(0);
  int nb_som_dim = les_sommets.dimension_int(1);
  for(int_t i=0; i<nb_som; i++)
    for(int j=0; j<nb_som_dim; j++)
      les_sommets(i,j)=Les_Nums[les_sommets(i,j)];
}

/*! @brief Returns the domain associated with the boundary (const version).
 *
 * @return (Domaine&) the domain associated with the boundary
 */
template <typename _SIZE_>
const typename Frontiere_32_64<_SIZE_>::Domaine_t& Frontiere_32_64<_SIZE_>::domaine() const
{
  return le_dom.valeur();
}

/*! @brief Returns the domain associated with the boundary.
 *
 * @return (Domaine&) the domain associated with the boundary
 */
template <typename _SIZE_>
typename Frontiere_32_64<_SIZE_>::Domaine_t& Frontiere_32_64<_SIZE_>::domaine()
{
  return le_dom.valeur();
}


/*! @brief Adds the vertices (and faces) of the given boundary to this object (Frontiere_32_64).
 *
 * @param (Frontiere_32_64& front) the boundary to "add" to this object
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::add(const Frontiere_32_64& front)
{
  const Faces_t& a_ajouter=front.faces();
  int_t nbf1=les_faces.nb_faces();
  int_t nbf2=a_ajouter.nb_faces();
  //max to treat the case where my front is empty
  int nbs=std::max(les_faces.nb_som_faces(), a_ajouter.nb_som_faces());
  int_t face;
  les_faces.les_sommets().resize(nbf1+nbf2, nbs);
  for(face=0; face<nbf2; face++)
    for(int som=0; som<nbs; som++)
      les_faces.sommet(nbf1+face, som)=a_ajouter.sommet(face, som);

  if(a_ajouter.voisins().nb_dim() == 1)
    return;

  int nb_voisins = a_ajouter.voisins().dimension_int(1);
  les_faces.voisins().resize(nbf1+nbf2, nb_voisins);
  for(face=0; face<nbf2; face++)
    for(int voisin=0; voisin<nb_voisins; voisin++)
      les_faces.voisin(nbf1+face, voisin)=a_ajouter.voisin(face, voisin);
}

/*! @brief Creates an array with one "row" per face of this boundary.
 *
 * See MD_Vector_tools::creer_tableau_distribue()
 */
template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::creer_tableau_faces(Array_base& v, RESIZE_OPTIONS opt) const
{
  const MD_Vector& md = les_faces.les_sommets().get_md_vector();
  MD_Vector_tools::creer_tableau_distribue(md, v, opt);
}

/*! @brief Returns the trace on the boundary of the element-based array y.
 *
 */
template <>
void Frontiere_32_64<int>::trace_elem_local(const DoubleTab& y, DoubleTab& x) const
{
  const int size = nb_faces();
  int nb_compo_ = y.line_size();

  // Resize x if not already done
  if (x.size_array() == 0 && size != 0)
    x.resize(size, nb_compo_);
  else if (x.dimension(0) != size || nb_compo_ != x.line_size())
    {
      Cerr << "Call to Frontiere_32_64<int>::trace_elem with a DoubleTab x not located on boundary faces or don't have the same number of components" << finl;
      Process::exit();
    }
  for (int i = 0; i < size; i++)
    {
      int elem = faces().voisin(i, 0);
      if (elem == -1)
        elem = faces().voisin(i, 1);

      for(int j = 0; j < nb_compo_; j++)
        x(i, j) = y(elem, j);
    }
}

/*! @brief Returns the trace on the boundary of the node-based array y.
 *
 */
template <>
void Frontiere_32_64<int>::trace_som_local(const DoubleTab& y, DoubleTab& x) const
{
  const IntTab& som_face = les_sommets_des_faces();
  const int size = nb_faces();
  int nb_compo_ = y.line_size();
  const int nsomfa = som_face.dimension_int(1);

  // Resize x if not already done
  if (x.size_array() == 0 && size != 0)
    x.resize(size, nb_compo_);
  else if (x.dimension(0) != size || nb_compo_ != x.line_size())
    {
      Cerr << "Call to Frontiere_32_64<int>::trace_elem with a DoubleTab x not located on boundary faces or don't have the same number of components" << finl;
      Process::exit();
    }

  for (int i = 0; i < size; i++)
    for(int j=0; j<nb_compo_; j++)
      {
        double s = 0.;
        for (int isom = 0; isom < nsomfa; isom++)
          s += y(som_face(i, isom), j);

        s /= nsomfa;
        x(i, j) = s;
      }
}

/*! @brief Returns the trace on the boundary of the face-based array y.
 *
 */
template <>
void Frontiere_32_64<int>::trace_face_local(const DoubleVect& y, DoubleVect& x) const
{
  Cerr << "Frontiere_32_64<int>::trace_face(const DoubleVect_t& y, DoubleVect_t& x) const" << finl;
  Cerr << "not coded yet." << finl;
  Process::exit();
}

/*! @brief Returns the trace on the boundary of the face-based array y (DoubleTab version).
 */
template <>
void Frontiere_32_64<int>::trace_face_local(const DoubleTab& y, DoubleTab& x) const
{
  int size = nb_faces();
  int M = y.line_size(), N = x.line_size();
  assert(x.dimension(0)==size);
  for (int i = 0; i < size; i++)
    {
      int face = num_premiere_face() + i;
      for (int n = 0; n < N; n++)
        x.addr()[N * i + n] = y.addr()[M * face + n];
    }
}

template <>
void Frontiere_32_64<int>::trace_som_distant(const DoubleTab&, DoubleTab&) const
{
  Cerr<<que_suis_je()<<"::trace_som_distant not implemented "<<finl;
  Process::exit();
}

template <>
void Frontiere_32_64<int>::trace_elem_distant(const DoubleTab&, DoubleTab&) const
{
  Cerr<<que_suis_je()<<"::trace_elem_distant not implemented "<<finl;
  Process::exit();
}

template <>
void Frontiere_32_64<int>::trace_face_distant(const DoubleTab&, DoubleTab&) const
{
  Cerr<<que_suis_je()<<"::trace_face_distant not implemented "<<finl;
  Process::exit();
}

template <>
void Frontiere_32_64<int>::trace_face_distant(const DoubleVect&, DoubleVect&) const
{
  Cerr<<que_suis_je()<<"::trace_face_distant not implemented "<<finl;
  Process::exit();
}



// 64 bit versions should never be called:

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_elem_local(const DoubleTab& y, DoubleTab& x) const
{
  assert(false);
  throw;
}

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_som_local(const DoubleTab& y, DoubleTab& x) const
{
  assert(false);
  throw;
}

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_face_local(const DoubleVect& y, DoubleVect& x) const
{
  assert(false);
  throw;
}

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_face_local(const DoubleTab& y, DoubleTab& x) const
{
  assert(false);
  throw;
}

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_som_distant(const DoubleTab&, DoubleTab&) const
{
  assert(false);
  throw;
}

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_elem_distant(const DoubleTab&, DoubleTab&) const
{
  assert(false);
  throw;
}

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_face_distant(const DoubleTab&, DoubleTab&) const
{
  assert(false);
  throw;
}

template <typename _SIZE_>
void Frontiere_32_64<_SIZE_>::trace_face_distant(const DoubleVect&, DoubleVect&) const
{
  assert(false);
  throw;
}



template class Frontiere_32_64<int>;
#if INT_is_64_ == 2
template class Frontiere_32_64<trustIdType>;
#endif


