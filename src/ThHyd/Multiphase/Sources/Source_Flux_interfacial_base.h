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

#ifndef Source_Flux_interfacial_base_included
#define Source_Flux_interfacial_base_included

#include <Sources_Multiphase_base.h>
#include <Correlation_base.h>
#include <TRUST_Ref.h>

/*! @brief Interfacial flux source term (implemented in PolyMAC_HFV) of the form:
 *
 *     F_{kl} = - F_{lk} = - C_{kl} (u_k - u_l).
 *     The coefficient C_{kl} is computed by the Coefficient_Flux_interfacial_base hierarchy.
 *
 * @sa Source_base
 */
class Source_Flux_interfacial_base : public Sources_Multiphase_base
{
  Declare_base(Source_Flux_interfacial_base);
public :
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;
  void mettre_a_jour(double temps) override;
  void completer() override;

  /* wall-to-interface flux (per cell) and its derivatives */
  /* these arrays are only accessible in the energy equation term and are used by others */
  /* they can be filled by the diffusion operator (wall flux), a fuel rod thermal module, etc. */
  DoubleTab& qpi() const;    // qpi(e, k, l) : power deposited in element e from phase k to interface (k, l); unit: W
  DoubleTab& dT_qpi() const; // dT_qpi(e, k, l, n) : its derivative w.r.t. T[n]
  DoubleTab& da_qpi() const; // da_qpi(e, k, l, n) : its derivative w.r.t. alpha[n]
  DoubleTab& dp_qpi() const; // dp_qpi(e, k, l)    : its derivative w.r.t. p

private:
  mutable DoubleTab qpi_, dT_qpi_, da_qpi_, dp_qpi_;
  OBS_PTR(Correlation_base) correlation_; // correlation providing the interfacial flux coefficient
  int is_turb_ = 0;
  double dv_min = -1.;
  double mod2grp = -1. ;
};

#endif /* Source_Flux_interfacial_base_included */
