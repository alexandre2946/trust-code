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

#include <Elem_geom_base.h>
#include <Domaine.h>

Implemente_base_32_64(Elem_geom_base_32_64,"Elem_geom_base",Objet_U);

/*! @brief Does nothing.
 *
 * @param s An output stream.
 * @return The output stream.
 */
template<class _SIZE_>
Sortie& Elem_geom_base_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  return s;
}

/*! @brief Does nothing.
 *
 * @param s An input stream.
 * @return The input stream.
 */
template<class _SIZE_>
Entree& Elem_geom_base_32_64<_SIZE_>::readOn(Entree& s )
{
  return s;
}


/*! @brief Returns the face parameter if the specified face type matches that of the geometric element.
 *
 * @param face Face index.
 * @param type A face type.
 * @throws Bad face type specified.
 */
template <typename _SIZE_>
int Elem_geom_base_32_64<_SIZE_>::num_face(int face, Type_Face& type) const
{
  assert(type==type_face());
  return face;
}

/*! @brief Creates the faces of the specified geometric element of the domain using the given face type.
 *
 * @param les_faces Faces object to fill.
 * @param num_elem Index of the element whose faces are to be created.
 * @param type Face type to create.
 */
template <typename _SIZE_>
void Elem_geom_base_32_64<_SIZE_>::creer_faces_elem(Faces_t& les_faces ,
                                                    _SIZE_ num_elem,
                                                    Type_Face type) const
{
  int type_id=0;
  for(; ((type_id<nb_type_face())
         &&(type!=type_face(type_id))); type_id++)
    ;
  const IntTab_t& les_Polys = mon_dom->les_elems();
  assert(les_Polys.dimension_tot(0) > num_elem);
  int face, i;
  les_faces.dimensionner(nb_faces(type_id));
  les_faces.associer_domaine(mon_dom.valeur());

  for(face=0; face<nb_faces(type_id); face++)
    {
      int face_id=num_face(face, type);
      for (i=0; i<nb_som_face(type_id); i++)
        les_faces.sommet(face_id,i) = les_Polys(num_elem,face_sommet(face_id, i));
      les_faces.completer(face_id, num_elem);
    }
}


/*! @brief Computes the centers of mass of all elements in the domain associated with this geometric element.
 *
 * @param tab_xp Array to fill with the coordinates of the centers of mass.
 */
template <typename _SIZE_>
void Elem_geom_base_32_64<_SIZE_>::calculer_centres_gravite(DoubleTab_t& tab_xp) const
{
  int_t nb_elem;
  if(tab_xp.dimension(0)==0)
    {
      nb_elem = mon_dom->nb_elem_tot();
      tab_xp.resize(nb_elem,dimension);
    }
  else
    nb_elem=tab_xp.dimension(0);

  int nb_som_elem = nb_som();
  int dim = Objet_U::dimension;
  // ToDo create a type in View_Types.h ?
  ConstView<_SIZE_,2> les_Polys = mon_dom->les_elems().view_ro();
  //CIntTabView les_Polys = mon_dom->les_elems().view_ro();
  CDoubleTabView coord = mon_dom->coord_sommets().view_ro();
  DoubleTabView xp = tab_xp.view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_elem, KOKKOS_LAMBDA(const int_t num_elem)
  {
    for (int i = 0; i < dim; i++)
      xp(num_elem, i) = 0;
    int nb_som_reel = nb_som_elem;
    while (les_Polys(num_elem, nb_som_reel - 1) == -1) nb_som_reel--;
    for (int s = 0; s < nb_som_reel; s++)
      {
        int_t num_som = les_Polys(num_elem, s);
        for (int i = 0; i < dim; i++)
          xp(num_elem, i) += coord(num_som, i) / nb_som_reel;
      }
  });
  end_gpu_timer(__KERNEL_NAME__);
}

/*! @brief Exits with an error. This method is not pure virtual for convenience.
 *
 * @param faces_sommets Vertex indices of the faces.
 * @param face_normales Array to fill with face normals.
 */
template <typename _SIZE_>
void Elem_geom_base_32_64<_SIZE_>::calculer_normales(const IntTab_t& faces_sommets , DoubleTab_t& face_normales) const
{
  Cerr << "calculer_normales method is not coded for an element " << finl;
  Cerr << "of type que_suis_je() " << finl;
  exit();
}

/*! @brief Returns the number of face types of the geometric element.
 *
 * For example, a prism (class Prisme) has 2 face types: triangle and quadrangle.
 *
 * @return Always returns 1 for this base implementation.
 */
template <typename _SIZE_>
int Elem_geom_base_32_64<_SIZE_>::nb_type_face() const
{
  return 1;
}

/*! @brief Fills faces_som_local(i,j) giving for 0 <= i < nb_faces() and 0 <= j < nb_som_face(i) the local vertex index on the element.
 *
 * We have 0 <= faces_sommets_locaux(i,j) < nb_som().
 * If faces do not all have the same number of vertices, the number of columns equals the maximum,
 * and unused entries are set to -1.
 * Returns 1 if all faces have the same number of vertices, 0 otherwise.
 *
 * @param faces_som_local Table to fill with local face-vertex indices.
 * @return 1 if all faces have the same vertex count, 0 otherwise.
 */
template <typename _SIZE_>
int Elem_geom_base_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local) const
{
  Cerr << "Elem_geom_base::faces_sommets_locaux : error.\n"
       << " Method not implemented for the object "
       << que_suis_je() << finl;
  exit();
  return 0;
}

/*! @brief Same as Elem_geom_base::get_tab_faces_sommets_locaux but for edges: aretes_som_local.
 *
 * aretes_som_local.dimension(0) = number of edges on the reference element.
 * aretes_som_local.dimension(1) = 2 (number of vertices per edge).
 * aretes_som_local(i,j) = index of a vertex of the element (0 <= n < nb_vertices_per_element).
 *
 * @param aretes_som_local Table to fill with local edge-vertex indices.
 */
template <typename _SIZE_>
void Elem_geom_base_32_64<_SIZE_>::get_tab_aretes_sommets_locaux(IntTab& aretes_som_local) const
{
  Cerr << "Elem_geom_base::aretes_sommets_locaux : error.\n"
       << " Method not implemented for the object "
       << que_suis_je() << finl;
  exit();
}


template class Elem_geom_base_32_64<int>;
#if INT_is_64_ == 2
template class Elem_geom_base_32_64<trustIdType>;
#endif

