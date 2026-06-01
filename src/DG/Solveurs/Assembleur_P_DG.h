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

#ifndef Assembleur_P_DG_included
#define Assembleur_P_DG_included

#include <Matrice_Morse_Sym.h>
#include <Assembleur_base.h>
#include <Domaine_DG.h>
#include <TRUST_Ref.h>

class Domaine_Cl_DG;

/**
 * @brief Assembles the pressure Laplacian matrix for the DG incompressible Navier-Stokes solver.
 *
 * This class builds the global pressure matrix arising from the SIP (Symmetric Interior
 * Penalty) DG discretization of the pressure Poisson problem. It is the pressure-space
 * counterpart of the diffusion operator assembly in Op_Diff_DG_Elem, but operates
 * directly on the pressure unknown without diffusivity weighting (nu = 1).
 *
 * The assembled matrix corresponds to the bilinear form:
 *   a(p_h, q_h) = sum_T  integral_T grad(p_h).grad(q_h)
 *               - sum_f  0.5 * integral_f { grad(p_h) }.n * [q_h]   (consistency)
 *               - sum_f  0.5 * integral_f [p_h] * { grad(q_h) }.n   (symmetry)
 *               + sum_f  (eta_F/h_T) * integral_f [p_h] * [q_h]     (penalty)
 *
 * where the sums run over all elements T and all faces f (internal and Dirichlet boundary).
 *
 * The matrix is stored as a Matrice_Morse (unsymmetric storage, though the assembled
 * system is symmetric by construction). The sparsity pattern couples each pressure DOF
 * to all pressure DOFs in the same element and in all face-neighbouring elements.
 *
 * Additionally, the class:
 *  - Stores a velocity-reconstruction matrix rec (Matrice_Morse) used by
 *    corriger_vitesses() to apply the pressure correction -grad(dP) to the velocity.
 *  - Handles pressure referencing via modifier_solution(), which pins the minimum
 *    pressure to zero when no Dirichlet pressure condition is imposed (has_P_ref == 0).
 *  - Provides a stub for quasi-compressible flows (assembler_QC) which is not yet
 *    fully implemented and aborts at runtime.
 *
 * @sa Op_Diff_DG_Elem, Op_Grad_DG, Assembleur_base
 */
class Assembleur_P_DG: public Assembleur_base
{
  Declare_instanciable(Assembleur_P_DG);
public:
  void associer_domaine_dis_base(const Domaine_dis_base&) override;
  void associer_domaine_cl_dis_base(const Domaine_Cl_dis_base&) override;
  const Domaine_dis_base& domaine_dis_base() const override;
  const Domaine_Cl_dis_base& domaine_Cl_dis_base() const override;

  int assembler(Matrice&) override;
  int assembler_rho_variable(Matrice&, const Champ_Don_base& rho) override;
  int assembler_QC(const DoubleTab&, Matrice&) override;
  int assembler_mat(Matrice&, const DoubleVect&, int incr_pression, int resoudre_en_u) override;

  int modifier_secmem(DoubleTab&) override;
  int modifier_solution(DoubleTab&) override;

  void completer(const Equation_base&) override;
  inline const Equation_base& equation() const { return mon_equation.valeur(); }

  /* corrige les vitesses pour une correction en pression donnee de type (-Cp, Cv) */
  void corriger_vitesses(const DoubleTab& dP, DoubleTab& dv) const override
  {
    rec.ajouter_multvect(dP, dv);
    dv.echange_espace_virtuel();
  }

protected:
  OBS_PTR(Equation_base) mon_equation;
  OBS_PTR(Domaine_DG) le_dom_dg_;
  OBS_PTR(Domaine_Cl_DG) le_dom_Cl_dg_;

  DoubleTab les_coeff_pression;

  int has_P_ref = 0;
  int stencil_done = 0;

  Matrice_Morse rec; //for reconstructing the velocities
};

#endif /* Assembleur_P_DG_included */
