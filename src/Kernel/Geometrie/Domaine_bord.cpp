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
#include <Domaine_bord.h>
#include <Hexaedre.h>

Implemente_instanciable_32_64(Domaine_bord_32_64,"Domaine_bord",Domaine_32_64<_T_>);

/*! @brief for now exit()
 *
 */
template <typename _SIZE_>
Entree& Domaine_bord_32_64<_SIZE_>::readOn(Entree& is)
{
  this->exit();
  return is;
}

/*! @brief for now exit()
 *
 */
template <typename _SIZE_>
Sortie& Domaine_bord_32_64<_SIZE_>::printOn(Sortie& os) const
{
  this->exit();
  return os;
}

/*! @brief builds the domain by calling extraire_domaine_bord()
 *
 */
template <typename _SIZE_>
void Domaine_bord_32_64<_SIZE_>::construire_domaine_bord(const Domaine_bord_32_64<_SIZE_>::Domaine_t& source, const Nom& nom_bord)
{
  domaine_source_ = source;
  bord_source_ = nom_bord;
  extraire_domaine_bord(source, nom_bord, *this, renum_som_);
}

/*! @brief returns a reference to the source domain
 *
 */
template <typename _SIZE_>
const Domaine_32_64<_SIZE_>& Domaine_bord_32_64<_SIZE_>::get_domaine_source() const
{
  return domaine_source_;
}

/*! @brief returns the name of the source boundary
 *
 */
template <typename _SIZE_>
const Nom& Domaine_bord_32_64<_SIZE_>::get_nom_bord_source() const
{
  return bord_source_;
}

/*! @brief returns renum_som (for each vertex of domaine_bord, index of the same vertex in the domain)
 *
 */
template <typename _SIZE_>
const ArrOfInt_T<_SIZE_>& Domaine_bord_32_64<_SIZE_>::get_renum_som() const
{
  return renum_som_;
}

/*! @brief method to convert a face type to an element type (to be moved to the Faces class ?)
 *
 */
template <typename _SIZE_>
void type_face_to_type_elem(const Elem_geom_base_32_64<_SIZE_>& type_elem, const Type_Face& type_face, Motcle& type_elem_face)
{
  switch(type_face)
    {
    case Type_Face::vide_0D:
      type_elem_face = "??";
      break;
    case Type_Face::point_1D:
      type_elem_face = "??";
      break;
    case Type_Face::segment_2D:
      type_elem_face = "segment";
      break;
    case Type_Face::segment_2D_axi:
      type_elem_face = "segment";
      break;
    case Type_Face::triangle_3D:
      type_elem_face = "triangle";
      break;
    case Type_Face::quadrilatere_2D_axi:
      type_elem_face = "quadrangle_VEF";
      break;
    case Type_Face::quadrangle_3D:
      type_elem_face = (sub_type(Hexaedre,type_elem)?"rectangle":"quadrangle_VEF");
      break;
    case Type_Face::quadrangle_3D_axi:
      type_elem_face = "quadrangle_VEF";
      break;
    default:
      type_elem_face = "??";
    }
}

/*! @brief fills the domain "dest" with the vertices and faces of the boundary "nom_bord" of the domain "src".
 *
 * The vertices of the dest domain are only the vertices that are on a boundary face.
 *   The renum_som array is sized to dest.nb_som() and filled as follows:
 *   renum_som[i] is the index in the domain "src" of vertex i of the domain "dest".
 *
 */
template <typename _SIZE_>
void Domaine_bord_32_64<_SIZE_>::extraire_domaine_bord(const Domaine_t& src,
                                                       const Nom& nom_bord,
                                                       Domaine_t& dest,
                                                       ArrOfInt_t& renum_som)
{
  if (Process::is_parallel())
    {
      Cerr << "extraire_domaine_bord in parallel: the domain created will not have a distributed structure\n"
           << " (this will be done one day... ask to B.Mathieu)" << finl;
    }

  // The destination domain must be empty:
  assert(dest.nb_elem() == 0);
  // Domain initialization:
  // Choose a name for the domain
  dest.nommer(src.le_nom() + Nom("_") + nom_bord);
  // Element type of the dest domain:
  Motcle type_elem;
  type_face_to_type_elem(src.type_elem().valeur(), src.type_elem()->type_face(), type_elem);
  const std::string suff = !std::is_same<_SIZE_, int>::value ? "_64" : "";
  Nom type_elem_64 = type_elem + suff;
  dest.type_elem().typer(type_elem_64);
  dest.type_elem()->associer_domaine(dest);

  const Frontiere_t& front = src.frontiere(nom_bord);
  const int_t nb_faces = front.faces().nb_faces();
  const int nb_som_face = front.faces().nb_som_faces();
  const IntTab_t& faces_src = front.faces().les_sommets();
  IntTab_t& elem_dest = dest.les_elems();
  elem_dest.resize(nb_faces, nb_som_face);
  renum_som.reset();

  // renum_inverse: for each vertex of the source domain, its index in the destination domain:
  ArrOfInt_t renum_inverse(src.nb_som());
  renum_inverse= -1;
  int_t nb_som_dest = 0;
  for (int_t i = 0; i < nb_faces; i++)
    {
      for (int j = 0; j < nb_som_face; j++)
        {
          const int_t som = faces_src(i, j);
          // If the vertex has not yet been encountered, assign it an index in the dest domain:
          if (renum_inverse[som] < 0)
            {
              renum_som.append_array(som);
              renum_inverse[som] = nb_som_dest++;
            }
          elem_dest(i, j) = renum_inverse[som];
        }
    }
  // Copy of vertices used in the destination domain
  DoubleTab_t& som_dest = dest.les_sommets();
  const DoubleTab_t& som_src = src.les_sommets();
  const int dim = static_cast<int>(som_src.dimension(1));
  som_dest.resize(nb_som_dest, dim);
  for (int_t i = 0; i < nb_som_dest; i++)
    {
      const int_t som = renum_som[i];
      for (int j = 0; j < dim; j++)
        som_dest(i, j) = som_src(som, j);
    }

  // TODO if needed: initialize the vertex joint to have the common items,
  //  and others if necessary.
}

template class Domaine_bord_32_64<int>;
#if INT_is_64_ == 2
template class Domaine_bord_32_64<trustIdType>;
#endif

