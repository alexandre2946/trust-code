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

#ifndef DomaineCutter_Correspondance_included
#define DomaineCutter_Correspondance_included

#include <TRUSTArray.h>

/*! @brief Helper structure holding the correspondence between vertex and element indices of the global domain and a sub-domain.
 *
 *   The fields of this class are filled by DomaineCutter::construire_sous_domaine()
 *
 */
template <typename _SIZE_>
class DomaineCutter_Correspondance_32_64
{
public:
  using SmallArrOfTID_t = SmallArrOfTID_T<_SIZE_>;  // a small number of big values -> typically an array of global indices
  using BigArrOfInt_t = BigArrOfInt_T<_SIZE_>;      // a big number of small values -> typically a (huge) array of local indices

  // The index of the constructed sub-domain
  int partie_ = -1;
  // Local/global vertex correspondence:
  //   global_vertex_index = liste_sommets[local_index]
  SmallArrOfTID_t liste_sommets_;
  // Global/local vertex correspondence:
  //   local_index = liste_inverse_sommets_[global_vertex_index]
  //   equals -1 if the vertex is not in the part
  BigArrOfInt_t liste_inverse_sommets_;
  // Global/local element correspondence:
  //   local_index = liste_inverse_elements_[global_element_index]
  BigArrOfInt_t liste_inverse_elements_;
};

#endif
