/****************************************************************************
* Copyright (c) 2023, CEA
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

#ifndef Domaine_PolyMAC_HFV_included
#define Domaine_PolyMAC_HFV_included

#include <Domaine_PolyMAC_CDO.h>
#include <TRUSTLists.h>
#include <Conds_lim.h>

class Domaine_PolyMAC_HFV : public Domaine_PolyMAC_CDO
{
  Declare_instanciable(Domaine_PolyMAC_HFV);
public :
  void discretiser() override;
  void calculer_volumes_entrelaces() override;

  //for each element, normal * dual area associated with each of its edges (oriented like the edge)
  const DoubleTab& surf_elem_arete() const;

  //local matrices per element (Hodge operators) for performing interpolations:
  void M1(const DoubleTab *nu, int e, DoubleTab& m1) const; //dual edge normals -> edge tangentials: (nu|a|t_a.v)   = m1 (S_ea.v)
  void W1(const DoubleTab *nu, int e, DoubleTab& w1, DoubleTab& v_e, DoubleTab& v_ea) const; //edge tangentials -> dual edge normals: (nu S_ea.v)    = w1 (|a|t_a.v)
  //options for the nu tensor:
  //null -> nu = Id below
  //isotropic -> nu(n_e, N) with n_e = 1 (constant tensor) / nb_elem_tot() (per-element tensor), and N a number of components
  //anisotropic -> nu(n_e, N, D) (diagonal anisotropic) or nu(n_e, N, D, D) (full anisotropic)
  //the output matrix is of size (n_f, n_f, N) (for M2/W2) or (n_a, n_a, N) (for M1/W1)

private:
  mutable DoubleTab surf_elem_arete_; // normal * dual area vector
};

#endif /* Domaine_PolyMAC_HFV_included */
