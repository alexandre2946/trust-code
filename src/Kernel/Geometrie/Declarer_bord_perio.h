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
#ifndef Declarer_bord_perio_included
#define Declarer_bord_perio_included

#include <Interprete_geometrique_base.h>
#include <Connectivite_som_elem.h>
#include <Domaine_forward.h>
#include <Domaine_bord.h>
#include <Octree_Double.h>
#include <ArrOfBit.h>

/*! @brief This interpreter must be used in sequential (before mesh splitting) if the opposite vertices of a periodic boundary are not perfectly aligned.
 *
 *   (case of certain tetrahedral meshes where the mesher is overly constrained by the CAD geometry).
 *   It attempts to move the vertices to align them.
 *
 *  This interpreter corrects periodic boundaries to conform to TRUST requirements:
 *    - reorder the faces of the periodic boundary so that face i+n/2 is opposite to face i,
 *      with all faces [0 .. n/2-1] on one side and [n/2 .. n-1] on the other side
 *    - move the vertices of periodic faces if needed (if the CAD geometry is incorrect)
 *  Syntax:
 *   Declarer_bord_perio {
 *      domaine NOMDOMAINE
 *      bord    NOMBORDPERIO
 *      [ direction DIMENSION dx dy [ dz ] ]
 *      [ fichier_post BASENAME ]
 *   }
 */
template <typename _SIZE_>
class Declarer_bord_perio_32_64 : public Interprete_geometrique_base_32_64<_SIZE_>
{
  Declare_instanciable_32_64(Declarer_bord_perio_32_64);
public:
  using int_t = _SIZE_;
  using IntTab_t = IntTab_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using ArrOfDouble_t = ArrOfDouble_T<_SIZE_>;

  using ArrOfBit_t = ArrOfBit_32_64<_SIZE_>;
  using Octree_Double_t = Octree_Double_32_64<_SIZE_>;
  using Bord_t = Bord_32_64<_SIZE_>;
  using Domaine_bord_t = Domaine_bord_32_64<_SIZE_>;
  using Domaine_t = Domaine_32_64<_SIZE_>;

  Entree& interpreter_(Entree& is) override;
  void adapt_som_and_faces();

  Nom& nom_bord() { return nom_bord_; }
  const Nom& nom_bord() const { return nom_bord_; }

protected:
  void corriger_coordonnees_sommets_perio();

  Nom nom_bord_;
  ArrOfDouble direction_perio_;
  Nom nom_fichier_post_;
};

using Declarer_bord_perio = Declarer_bord_perio_32_64<int>;
using Declarer_bord_perio_64 = Declarer_bord_perio_32_64<trustIdType>;

#endif
