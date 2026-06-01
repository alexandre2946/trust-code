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

#ifndef Domaine_PolyMAC_MPFA_included
#define Domaine_PolyMAC_MPFA_included

#include <Domaine_PolyMAC_HFV.h>

class Domaine_PolyMAC_MPFA : public Domaine_PolyMAC_HFV
{
  Declare_instanciable(Domaine_PolyMAC_MPFA);
public :
  void discretiser() override;

  //stencil of the face gradient fgrad: fsten_eb([fsten_d(f), fsten_d(f + 1)[]) -> elements e, boundary faces ne_tot + f
  void init_stencils() const;
  mutable IntTab fsten_d, fsten_eb;

  //for a scalar field T at elements, interpolates [n_f.nu.grad T]_f
  //while exactly preserving fields satisfying [nu grad T]_e = const.
  //Inputs : N             : number of components
  //         is_p          : 1 if treating the pressure field (swap Neumann / Dirichlet)
  //         cls           : boundary conditions
  //         fcl(f, 0/1/2) : boundary condition data (type, index, local index) (cf. Champ_{P0,Face}_PolyMAC_MPFA)
  //         nu(e, n, ..)  : diffusivity at elements (optional)
  //         som_ext       : list of vertices to skip (e.g. direct treatment of Echange_Contact in Op_Diff_PolyMAC_MPFA_Elem)
  //         virt          : 1 if also wanting flux at virtual faces
  //         full_stencil  : 1 if wanting the full stencil (for dimensionner())
  //Outputs: phif_d(f, 0/1)                       : indices in phif_{e,c} / phif_{pe,pc} of the flux at f in [phif_d(f, 0/1), phif_d(f + 1, 0/1)[
  //         phif_e(i), phif_c(i, n, c)           : local indices/coefficients (no Echange_contact) and diagonal (independent components)
  void fgrad(int N, int is_p, const Conds_lim& cls, const IntTab& fcl, const DoubleTab *nu, const IntTab *som_ext,
             int virt, int full_stencil, IntTab& phif_d, IntTab& phif_e, DoubleTab& phif_c) const;

  //MD_Vectors for Champ_Face_PolyMAC_MPFA (faces + d x elems)
  MD_Vector mdv_ch_face;

protected:
  mutable int first_fgrad_ = 1; //to print the "MPFA-O MPFA-O(h) VFSYM" message only once per computation
};

#endif /* Domaine_PolyMAC_MPFA_included */
