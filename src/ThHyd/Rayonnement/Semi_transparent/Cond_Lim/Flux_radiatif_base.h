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

#ifndef Flux_radiatif_base_included
#define Flux_radiatif_base_included

#include <Neumann_paroi.h>
class Equation_base;

class Flux_radiatif_base : public Neumann_paroi
{
  Declare_base(Flux_radiatif_base);

public :
  inline int compatible_avec_eqn(const Equation_base&) const override { return 1; }
  inline Champ_front_base& emissivite() { return emissivite_.valeur(); }
  inline const Champ_front_base& emissivite() const { return emissivite_.valeur(); }
  inline double& A() { return A_; }
  inline const double& A() const { return A_; }
  inline Champ_front_base& flux_radiatif() { return flux_radiatif_.valeur(); }
  inline const Champ_front_base& flux_radiatif() const { return flux_radiatif_.valeur(); }

  virtual void calculer_flux_radiatif(const Equation_base& eq_temp)=0;
  void completer() override;

  double flux_impose(int i) const override;
  double flux_impose(int i,int j) const override;

protected :
  double A_ = -123.;
  OWN_PTR(Champ_front_base) emissivite_;
  OWN_PTR(Champ_front_base) flux_radiatif_;
};

#endif /* Flux_radiatif_base_included */
