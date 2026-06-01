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

#ifndef Op_Diff_PolyMAC_HFV_base_included
#define Op_Diff_PolyMAC_HFV_base_included

#include <Op_Diff_PolyMAC_CDO_Gen_base.h>

/*! @brief class Op_Diff_PolyMAC_HFV_base
 *
 *  Base class for PolyMAC_HFV diffusion operators.
 *
 *
 */
class Op_Diff_PolyMAC_HFV_base: public Op_Diff_PolyMAC_CDO_Gen_base
{
  Declare_base(Op_Diff_PolyMAC_HFV_base);
public:
  void completer() override;

  void mettre_a_jour(double t) override;

  /* extended versions of dimensionner/ajouter_blocs allowing treatment of auxiliary variables only */
  virtual void dimensionner_blocs_ext(int aux_only, matrices_t matrices, const tabs_t& semi_impl = { }) const = 0;
  virtual void ajouter_blocs_ext(int aux_only, matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = { }) const = 0;

  /* implementations of dimensionner/ajouter_blocs from these methods */
  int has_interface_blocs() const override { return 1; }

  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = { }) const override
  {
    dimensionner_blocs_ext(0, matrices, semi_impl);
  }

  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = { }) const override
  {
    ajouter_blocs_ext(0, matrices, secmem, semi_impl);
  }

  void update_nu() const override; //update

protected:
  double t_last_nu_ = -1e10; //to detect when nu and auxiliary variables need to be recomputed

  /* management of auxiliary variables in semi-implicit mode */
  void update_aux(double t) const;
  mutable double t_last_aux_ = -1e10; /* last time at which they were computed */
  mutable int use_aux_ = 0; /* whether auxiliary variables are stored in var_aux */
  mutable Matrice_Bloc mat_aux; /* system to solve: mat.var_aux = secmem */
  mutable DoubleTab var_aux;
  mutable SolveurSys solv_aux; //solver
};

#endif /* Op_Diff_PolyMAC_HFV_base_included */
