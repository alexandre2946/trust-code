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

#ifndef Reorder_Mesh_included
#define Reorder_Mesh_included

#include <TRUSTTabs_forward.h>
#include <Objet_U.h>
#include <Entree.h>

/**! @brief Various methods of reordering Morton / Hilbert
 * See wiki: https://en.wikipedia.org/wiki/Z-order_curve
 */
enum class Reorder_Algo { None, Morton, Hilbert };

/**! @brief Reorder_Mesh allows the user to trigger the renumbering of the mesh entities.
 *
 * Renumbering can be done according to a Morton or Hilbert scheme (see Morton curve or Hilbert
 * curves on Wikipedia) thus improving the data co-localisation, in the following sense:
 *
 * - when iterating through the faces of the mesh
 * - and then retrieving the adjacent elements (indirection via face_voisin typically)
 * - ensures that the elements that are thus scanned are themselves not too wide spread in memory
 *
 * This helps improving memory access performance.
 * See .cpp file for more explanations on how this works.
 *
 * This object is used as a discretisation option - see class Discret_Thyd, member reorder_
 */
class Reorder_Mesh: public Objet_U
{
  Declare_instanciable(Reorder_Mesh);

public:
  template<typename _SIZE_>
  void compute_renumbering(const DoubleTab_T<_SIZE_>& points, ArrOfInt_T<_SIZE_>& renum) const;

  template<typename _SIZE_>
  void dump_to_file(const DoubleTab_T<_SIZE_>& points, const std::string& filename) const;

  // Accessors
  Reorder_Algo algo() const { return algo_; }
  bool skip_nodes() const { return no_nodes_; }
  bool skip_elems() const { return no_elems_; }
  bool skip_faces() const { return no_faces_; }
  bool is_dump() const { return dump_; }

protected:

  Reorder_Algo algo_ = Reorder_Algo::None;  ///< Reordering algorithm
  bool no_nodes_ = false;                   ///< Whether to skip vertices in reordering
  bool no_elems_ = false;                   ///< Whether to skip elements in reordering
  bool no_faces_ = false;                   ///< Whether to skip faces in reordering
  bool dump_ = false;                       ///< Whether to dump previous and new positions into text files.
};


#endif
