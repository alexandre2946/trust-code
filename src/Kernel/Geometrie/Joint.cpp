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

#include <Joint.h>
#include <Domaine.h>

Implemente_instanciable_32_64(Joint_32_64,"Joint_32_64",Frontiere);

// ***************************************************************
//  Implementation of the Joint_32_64 class
// ***************************************************************

/*! @brief Writes the joint to an output stream.
 *
 * We write:
 *       - the boundary
 *       - the neighboring PE
 *       - the thickness
 *       - the vertices
 *
 */
template <typename _SIZE_>
Sortie& Joint_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  Frontiere_t::printOn(s) ;
  s << "format_joint_250507" << finl;
  s << PEvoisin_ << finl;
  s << epaisseur_ << finl;
  s << joint_item(JOINT_ITEM::SOMMET).items_communs() << finl;
  s << joint_item(JOINT_ITEM::ELEMENT).items_distants() << finl;
  return s ;
}

/*! @brief Reads a joint from an input stream.
 *
 */
template <typename _SIZE_>
Entree& Joint_32_64<_SIZE_>::readOn(Entree& s)
{
  Frontiere_t::readOn(s) ;
  Nom format;
  s >> format;
  if (format!="format_joint_250507")
    {
      if (Process::je_suis_maitre())
        {
          Cerr << "The format of .Zones is not recognized because probably too old." << finl;
          Cerr << "Split your mesh with an executable which is more recent" << finl;
          Cerr << "than the version 1.5.2 build 240507." << finl;
        }
      Process::exit();
    }
  s >> PEvoisin_;
  s >> epaisseur_;
  s >> set_joint_item(JOINT_ITEM::SOMMET).set_items_communs();
  s >> set_joint_item(JOINT_ITEM::ELEMENT).set_items_distants();
  return s ;
}

/*! @brief
 *
 */
template <typename _SIZE_>
void Joint_32_64<_SIZE_>::dimensionner(int i)
{
  this->faces().dimensionner(i);
  this->faces().voisins() = -1;
  this->faces().les_sommets() = -1;
}


/*! @brief Adds faces to the boundary (to the joint). See Frontiere::ajouter_faces(const IntTab&).
 *
 * @param sommets Array containing the vertex indices of the faces to add.
 */
template <typename _SIZE_>
void Joint_32_64<_SIZE_>::ajouter_faces(const IntTab_t& sommets)
{
  this->faces().ajouter(sommets);
  this->faces().voisins() = -1;
}

/*! @brief Returns the joint information for a given geometric item type, for filling the structures.
 *
 * These structures are generally filled by the Scatter class. Notable exceptions:
 * Domaine_VDF and Domaine_VF for face renumbering.
 *
 * @param item The geometric item type.
 * @return Reference to the corresponding Joint_Items structure.
 */
template <typename _SIZE_>
typename Joint_32_64<_SIZE_>::Joint_Items_t& Joint_32_64<_SIZE_>::set_joint_item(JOINT_ITEM item)
{
  switch(item)
    {
    case JOINT_ITEM::SOMMET:
      return joint_sommets_;
    case JOINT_ITEM::ELEMENT:
      return joint_elements_;
    case JOINT_ITEM::FACE:
      return joint_faces_;
    case JOINT_ITEM::ARETE:
      return joint_aretes_;
    case JOINT_ITEM::FACE_FRONT:
      return joint_faces_front_;
    default:
      Cerr << "Error in Joint_32_64<_SIZE_>::set_joint_item, invalid item number : " << (int)item << finl;
      Process::exit();
    }
  return joint_sommets_; // never arrive here
}

/*! @brief Returns the joint information for the requested type (read-only).
 *
 * @param item The geometric item type.
 * @return Const reference to the corresponding Joint_Items structure.
 */
template <typename _SIZE_>
const typename Joint_32_64<_SIZE_>::Joint_Items_t& Joint_32_64<_SIZE_>::joint_item(JOINT_ITEM item) const
{
  switch(item)
    {
    case JOINT_ITEM::SOMMET:
      return joint_sommets_;
    case JOINT_ITEM::ELEMENT:
      return joint_elements_;
    case JOINT_ITEM::FACE:
      return joint_faces_;
    case JOINT_ITEM::ARETE:
      return joint_aretes_;
    case JOINT_ITEM::FACE_FRONT:
      return joint_faces_front_;
    default:
      Cerr << "Error in Joint_32_64<_SIZE_>::set_joint_item, invalid item number : " << (int)item << finl;
      Process::exit();
    }
  return joint_sommets_; // never arrive here
}


template class Joint_32_64<int>;
#if INT_is_64_ == 2
template class Joint_32_64<trustIdType>;
#endif


