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

#ifndef Source_Portance_interfaciale_base_included
#define Source_Portance_interfaciale_base_included

#include <Sources_Multiphase_base.h>
#include <Correlation_base.h>
#include <Champs_Fonc.h>

/*! @brief Interfacial lift force source term of the form:
 *
 *     F_{n_l} = - F_{k} = C_{n_l, k} (u_k - u_n_l) x rot(u_n_l), where phase
 *     n_l is the carrier liquid phase and k is a gas phase.
 *     The coefficient C_{n_l, k} is computed by the Portance_interfaciale_base hierarchy.
 *
 * @sa Source_base
 */
class Source_Portance_interfaciale_base : public Sources_Multiphase_base
{
  Declare_base(Source_Portance_interfaciale_base);
public :
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override { /* Do nothing : 100% explicit */ }
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override = 0;
  const Correlation_base& correlation() const { return correlation_; }
  void creer_champ(const Motcle& motlu) override;
  void completer() override;

protected:
  OWN_PTR(Correlation_base) correlation_; // correlation providing the interfacial lift coefficient
  int n_l = -1; // liquid phase index
  double beta_ = 1. ; // To adjust the force in .data
  double g_ = 9.81;
  OWN_PTR(Champ_Fonc_base)  wobble; // postreatment
  OWN_PTR(Champ_Fonc_base)  C_lift; // postreatment
};

#endif /* Source_Portance_interfaciale_base_included */
