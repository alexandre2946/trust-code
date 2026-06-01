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

#ifndef Joint_included
#define Joint_included

#include <Frontiere.h>
#include <Joint_Items.h>

enum class JOINT_ITEM { SOMMET, ELEMENT, FACE, ARETE, FACE_FRONT };

/*! @brief The Joint class is a Frontiere that contains the joint faces and vertices with the
 *  neighboring domain PEvoisin() (for meshes distributed in parallel).
 *
 * It additionally carries in the Joint_Items the information needed to build
 * distributed arrays indexed by geometric items (vertices, faces, elements, edges, boundary faces).
 *
 * @sa Scatter Joint_Item
 */
template <typename _SIZE_>
class Joint_32_64 : public Frontiere_32_64<_SIZE_>
{

  Declare_instanciable_32_64(Joint_32_64);

public:

  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;

  using Joint_Items_t = Joint_Items_32_64<_SIZE_>;
  using Frontiere_t = Frontiere_32_64<_SIZE_>;

  void affecte_PEvoisin(int num) { PEvoisin_ = num; }
  void affecte_epaisseur(int ep) { epaisseur_ = ep; }
  int PEvoisin() const { return PEvoisin_; }
  int epaisseur() const { return epaisseur_; }

  void dimensionner(int);
  void ajouter_faces(const IntTab_t&);

  // Accessors for compatibility with the previous version
  // (to be removed soon)
  const IntTab_t&    renum_virt_loc() const { return joint_item(JOINT_ITEM::SOMMET).renum_items_communs(); }
  const ArrOfInt_t& esp_dist_elems() const   { return joint_item(JOINT_ITEM::ELEMENT).items_distants(); }
  const ArrOfInt_t& esp_dist_sommets() const { return joint_item(JOINT_ITEM::SOMMET).items_distants(); }
  const ArrOfInt_t& esp_dist_faces() const   { return joint_item(JOINT_ITEM::FACE).items_distants(); }

  // New interface to access joint data:
  Joint_Items_t&        set_joint_item(JOINT_ITEM type);
  const Joint_Items_t& joint_item(JOINT_ITEM type) const;

private:
  // Index of the neighboring PE corresponding to this joint
  int PEvoisin_=-1;

  // Joint thickness, used by Scatter::calculer_espace_distant_elements
  int epaisseur_=-1;

  Joint_Items_t joint_sommets_;
  Joint_Items_t joint_elements_;
  Joint_Items_t joint_faces_;
  Joint_Items_t joint_faces_front_;
  Joint_Items_t joint_aretes_;
};

using Joint = Joint_32_64<int>;
using Joint_64 = Joint_32_64<trustIdType>;

#endif
