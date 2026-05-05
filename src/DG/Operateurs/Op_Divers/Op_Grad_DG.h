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

#ifndef Op_Grad_DG_included
#define Op_Grad_DG_included

#include <Domaine_DG.h>
#include <Operateur_Grad.h>
#include <TRUST_Ref.h>

class Domaine_Cl_DG;

/**
 * @brief DG gradient operator acting on a DG pressure field to produce a velocity-space residual.
 *
 * This class implements the DG discretization of the pressure gradient operator used in
 * the incompressible Navier-Stokes momentum equation. It assembles the coupling between
 * pressure DOFs and velocity DOFs through the interface_blocs mechanism.
 *
 * The weak formulation integrates -grad(p_h) against velocity test functions v_h:
 *  - *Volume term*:        integral of grad(p_h) . v_h  (element-wise, by parts)
 *  - *Internal face term*: integral of {{p_h}} * [v_h . n]  (average-jump coupling)
 *
 * This formulation is the transpose of the divergence operator assembled by Op_Div_DG
 * in the continuous sense. However, Op_Grad_DG and Op_Div_DG use independent
 * implementations: Op_Grad_DG integrates grad(phi_p) . phi_v on element volumes and
 * {{phi_p}} * [phi_v . n] on faces, while Op_Div_DG integrates phi_p * div(phi_v) and
 * [phi_v . n] * {{phi_p}}. Both produce the same matrix up to a sign when the two
 * basis sets are compatible, but they are kept separate to allow independent tuning.
 *
 * @sa Op_Div_DG, Operateur_Grad_base
 */
class Op_Grad_DG: public Operateur_Grad_base
{
  Declare_instanciable(Op_Grad_DG);
public:

  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;

  inline int has_interface_blocs() const override { return 1; }
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override;

  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = { }) const override;
  DoubleTab& calculer(const DoubleTab&, DoubleTab&) const override;
  int impr(Sortie& os) const override;

protected:
  OBS_PTR(Domaine_DG) le_dom_DG;
  OBS_PTR(Domaine_Cl_DG) le_dcl_DG;
};

#endif /* Op_Grad_DG_included */
