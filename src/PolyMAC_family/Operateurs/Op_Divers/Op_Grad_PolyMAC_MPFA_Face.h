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

#ifndef Op_Grad_PolyMAC_MPFA_Face_included
#define Op_Grad_PolyMAC_MPFA_Face_included

#include <Op_Grad_PolyMAC_HFV_Face.h>
#include <cfloat>

/*! @brief class Op_Grad_PolyMAC_MPFA_Face
 *
 * @brief Gradient operator for PolyMAC_MPFA discretization, computing the gradient of a
 *        Champ_Elem_PolyMAC_MPFA field (pressure).
 *
 * @sa Operateur_Grad_base
 */
class Op_Grad_PolyMAC_MPFA_Face: public Op_Grad_PolyMAC_HFV_Face
{
  Declare_instanciable(Op_Grad_PolyMAC_MPFA_Face);
public:
  void completer() override;

  /* interface {dimensionner,ajouter}_blocs -> cf Equation_base.h */
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = { }) const override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = { }) const override;

  void update_grad(int full_stencil = 0) const;

  /* public for use by Assembleur_P_PolyMAC_MPFA : [grad p]_f */
  mutable IntTab fgrad_d, fgrad_e;
  mutable DoubleTab fgrad_c;

protected:
  mutable double last_gradp_ = -DBL_MAX; //last time used to interpolate grad p (set to DBL_MAX if grad p has not been re-interpolated)
  mutable std::map<int, std::map<int, double>> dgp_gb_;
  mutable std::map<int, std::map<int, double>> dgb_v_;
  mutable std::vector<std::vector<std::pair<int, double>>> dgf_pe_;
  mutable std::vector<std::vector<std::pair<int, double>>> dgf_gb_;
};

#endif /* Op_Grad_PolyMAC_MPFA_Face_included */
