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

#include <Bords.h>

Implemente_instanciable_32_64(Bords_32_64, "Bords", LIST(Bord_32_64<_T_>));
// XD list_bord listobj list_bord BRACE bord_base NO_COMMA The block sides.

template <typename _SIZE_>
Sortie& Bords_32_64<_SIZE_>::printOn(Sortie& os) const { return LIST(Bord_32_64<_SIZE_>)::printOn(os); }

template <typename _SIZE_>
Entree& Bords_32_64<_SIZE_>::readOn(Entree& is) { return LIST(Bord_32_64<_SIZE_>)::readOn(is); }

/*! @brief Associates a domain to all boundaries in the list.
 *
 * @param (Domaine& un_domaine) the domain to associate to the boundaries in the list
 */
template <typename _SIZE_>
void Bords_32_64<_SIZE_>::associer_domaine(const Domaine_t& un_domaine)
{
  for (auto &itr : *this) itr.associer_domaine(un_domaine);
}

/*! @brief Returns the total number of faces of all boundaries in the list.
 *
 * @return (int) the total number of faces of all boundaries in the list
 */
template <typename _SIZE_>
typename Bords_32_64<_SIZE_>::int_t Bords_32_64<_SIZE_>::nb_faces() const
{
  int_t nombre = 0;

  for (const auto &itr : *this) nombre += itr.nb_faces();

  return nombre;
}

/*! @brief Returns the total number of faces of the specified type, for all boundaries in the list.
 *
 * @param (Type_Face type) the type of faces to count
 * @return (int) the total number of faces of the specified type, for all boundaries in the list
 */
template <typename _SIZE_>
typename Bords_32_64<_SIZE_>::int_t Bords_32_64<_SIZE_>::nb_faces(Type_Face type) const
{
  int_t nombre = 0;

  for (const auto &itr : *this)
    if (type == itr.faces().type_face()) nombre += itr.nb_faces();

  return nombre;
}


template class Bords_32_64<int>;
#if INT_is_64_ == 2
template class Bords_32_64<trustIdType>;
#endif


