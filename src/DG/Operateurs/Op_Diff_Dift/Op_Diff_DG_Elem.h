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

#ifndef Op_Diff_DG_Elem_included
#define Op_Diff_DG_Elem_included

#include <Op_Diff_DG_base.h>
class Matrice_Morse;

/**
 * @brief Concrete DG diffusion operator acting on element-based unknowns.
 *
 * This class implements the Symmetric Interior Penalty (SIP) Galerkin discretization
 * of the diffusion operator for scalar and vector fields discretized on DG elements.
 * It assembles both the global stiffness matrix and the right-hand side contribution
 * through the interface_blocs mechanism (has_interface_blocs() returns 1).
 *
 * The SIP formulation introduces three families of face integrals at every internal face
 * shared by elements T0 and T1:
 *  - **Consistency term**: -0.5 * integral of { nu * grad(phi_i) } . n * [phi_j]
 *  - **Symmetry term**:    -0.5 * integral of { nu * grad(phi_j) } . n * [phi_i]  (transpose of consistency)
 *  - **Penalty term**:      gamma * (eta_F / h_T) * integral of [phi_i] * [phi_j]
 *
 * where gamma = 2*nu_0*nu_1/(nu_0+nu_1) is the harmonic mean of the diffusivities on
 * each side, eta_F is the element-local penalty coefficient, and h_T is the minimum
 * characteristic size of the two adjacent elements.
 *
 * At boundary faces the same three terms are applied to a single element, with the
 * boundary condition value contributing to the right-hand side through
 * contribuer_au_second_membre(). Supported boundary conditions are:
 *  - Neumann / Neumann_paroi: flux is added directly to the RHS.
 *  - Dirichlet / Dirichlet_homogene: enforced weakly via SIP penalty + symmetry terms.
 *
 * The stencil (sparsity pattern) of the assembled matrix couples every element to all
 * its face-neighbours, as computed by dimensionner() and dimensionner_blocs().
 *
 * Cross-problem coupling (dimensionner_termes_croises, ajouter_termes_croises,
 * contribuer_termes_croises) is declared but not yet implemented (throws at runtime).
 *
 * @note The ordering of degrees of freedom inside an element follows:
 *       phi0.ex, phi1.ex, ..., phi_{nb_bfunc-1}.ex,
 *       phi0.ey, phi1.ey, ..., phi_{nb_bfunc-1}.ey, ...
 *       for vector fields, and phi0, phi1, ..., phi_{nb_bfunc-1} for scalar fields.
 */
class Op_Diff_DG_Elem: public Op_Diff_DG_base
{
  Declare_instanciable( Op_Diff_DG_Elem );

public:
  virtual void calculer_flux_bord(const DoubleTab& inco) const = delete; //TODO DG a calculer dans interface_blocs

  void modifier_pour_Cl(Matrice_Morse& la_matrice, DoubleTab& secmem) const override { }
  void completer() override;

  void dimensionner(Matrice_Morse& mat) const override;

  inline int has_interface_blocs() const override { return 1; }
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;

  void dimensionner_termes_croises(Matrice_Morse&, const Probleme_base& autre_pb, int nl, int nc) const override;
  void ajouter_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, DoubleTab& resu) const override;
  void contribuer_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, Matrice_Morse& matrice) const override;
  void contribuer_au_second_membre(DoubleTab& resu ) const override;

};


#endif /* Op_Diff_DG_Elem_included */
