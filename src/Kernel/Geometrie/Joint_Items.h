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

#ifndef Joint_Items_included
#define Joint_Items_included

#include <TRUSTTab.h>

/*! @brief Joint_Items holds the parallel distribution information for a particular geometric item type
 *  with a particular neighboring domain (item = vertex, element, face, etc.)
 *
 * These structures are initialised in Scatter and are then used, for example, to create a
 * distributed array indexed by geometric item indices.
 *
 * @sa class Joint
 */
template <typename _SIZE_>
class Joint_Items_32_64
{
public:

  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;

  Joint_Items_32_64();
  void reset();

  // To use these accessors, the structures must have been previously
  // initialised with set_xxx.
  int nb_items_reels() const;
  const ArrOfInt_t& items_communs() const {  return items_communs_; }
  const ArrOfInt_t& items_distants() const;
  int nb_items_virtuels() const;
  const IntTab_t& renum_items_communs() const;

  // Initialisation methods for the structures
  void set_nb_items_reels(int n);
  ArrOfInt_t& set_items_communs();
  ArrOfInt_t& set_items_distants();
  void set_nb_items_virtuels(int n);
  IntTab_t& set_renum_items_communs();

private:
  // Number of real items (allows building a distributed array)
  // using only the joint information.
  int nb_items_reels_;

  // List of items shared with the neighboring domain (the list is
  // ordered in the same way on the local domain and on the neighboring
  // domain => items_communs[i] on joint_j of domain_k represents the same
  // geometric entity as items_communs[i] on joint_k of domain_j)
  ArrOfInt_t items_communs_;

  // List of remote items to send to the neighboring domain
  // (the order of items in this list determines the order of appearance
  // of these items in the virtual space of the neighbor)
  ArrOfInt_t items_distants_;

  // Number of virtual items received from the neighboring domain.
  //  we have "nb_items_virtuels_ on joint_j of domain_k"
  //     = "items_distants.size_array() on joint_k of domain_j"
  int nb_items_virtuels_;

  // Correspondence between the local index of a shared item and the index
  // of the same item on the neighboring domain:
  // column 0 = index on the neighboring domain,
  // column 1 = index on the local domain
  // dimension(0) is equal to items_communs.size_array()
  // The order of items in the array is not specified
  IntTab_t renum_items_communs_;

  // What has been initialised?
  int flags_init_;
};

using Joint_Items = Joint_Items_32_64<int>;
using Joint_Items_64 = Joint_Items_32_64<trustIdType>;

#endif
