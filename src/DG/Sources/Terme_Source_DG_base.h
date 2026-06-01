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

#ifndef Terme_Source_DG_base_included
#define Terme_Source_DG_base_included

#include <Iterateur_Source_base.h>
#include <Source_base.h>

/**
 * @brief Abstract base class for all DG source terms.
 *
 * This class provides the minimal common infrastructure for source terms in the DG
 * framework. It owns an iterator (iter_) of type Iterateur_Source_base which encapsulates
 * the element-wise loop and delegates the actual source evaluation to a concrete
 * Evaluateur_Source. The iterator pattern allows the assembly strategy (e.g., looping
 * over elements) to be reused across different source term types.
 *
 * The interface_blocs mechanism is always active (has_interface_blocs() returns 1):
 *  - dimensionner_blocs() is a no-op since source terms add only RHS contributions
 *    (no matrix coupling).
 *  - ajouter_blocs() delegates to iter_->ajouter(secmem), which loops over elements
 *    and accumulates the evaluator output into the right-hand side.
 *
 * Derived classes must supply a concrete Iterateur_Source_base at construction time
 * and typically override associer_domaines() and associer_pb() to wire up the
 * evaluator with the domain geometry and physical fields.
 *
 * @sa Terme_Puissance_Thermique_DG_base, Terme_Source_Qdm_Elem_DG, Source_base
 */
class Terme_Source_DG_base : public Source_base
{
  Declare_base(Terme_Source_DG_base);
public:
  Terme_Source_DG_base(const Iterateur_Source_base& iter_base) { iter_ = iter_base; }

  int has_interface_blocs() const override { return 1; }
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override { }; //nothing to size
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override { iter_->ajouter(secmem); }
  void completer() override;

protected:
  OWN_PTR(Iterateur_Source_base) iter_;
};

#endif /* Terme_Source_DG_base_included */
