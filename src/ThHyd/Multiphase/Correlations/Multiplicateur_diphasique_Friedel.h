/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef Multiplicateur_diphasique_Friedel_included
#define Multiplicateur_diphasique_Friedel_included
#include <Multiplicateur_diphasique_base.h>

/*! @brief Two-phase multiplier using the Friedel correlation:
 *
 *     - applied to the liquid phase for alpha < alpha_min
 *     - applied to the vapor phase for alpha > alpha_max
 *
 *
 */

class Multiplicateur_diphasique_Friedel : public Multiplicateur_diphasique_base
{
  Declare_instanciable(Multiplicateur_diphasique_Friedel);
public:
  void coefficient(const double *alpha, const double *rho, const double *v, const double *f,
                   const double *mu, const double Dh, const double gamma, const double *Fk,
                   const double Fm, DoubleTab& coeff) const override;
protected:
  double alpha_min_ = 1, alpha_max_ = 1.1;
  int n_l = -1, n_g = -1, min_lottes_flinn_ = 0, min_sensas_ = 0; //indices of the friction phases (liquid, gas), minimum taken with Lottes-Flinn
};

#endif
