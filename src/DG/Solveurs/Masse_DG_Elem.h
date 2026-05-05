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

#ifndef Masse_DG_Elem_included
#define Masse_DG_Elem_included

#include <Masse_DG_base.h>

/**
 * @brief Concrete DG mass operator for element-based unknowns.
 *
 * This class provides the implementation of appliquer_impl(), which applies the
 * inverse mass matrix M^{-1} to a right-hand side vector. It is used for explicit
 * time stepping or for any operation requiring M^{-1} * v element-by-element.
 *
 * Two strategies are used depending on the basis orthonormalization flag:
 *  - **Orthonormal basis**: M is diagonal with entries equal to element volumes,
 *    so M^{-1} * v reduces to a simple element-wise division by volume via
 *    tab_divide_any_shape().
 *  - **Non-orthonormal basis**: M^{-1} is computed locally per element via
 *    BasisFunction::eval_invMassMatrix() and applied by a dense matrix-vector
 *    product for each spatial component independently.
 *
 * @sa Masse_DG_base, BasisFunction::eval_invMassMatrix
 */
class Masse_DG_Elem: public Masse_DG_base
{
  Declare_instanciable(Masse_DG_Elem);
public:
  DoubleTab& appliquer_impl(DoubleTab&) const override;
};

#endif /* Masse_DG_Elem_included */
