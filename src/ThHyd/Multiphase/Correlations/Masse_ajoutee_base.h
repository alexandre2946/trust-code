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

#ifndef Masse_ajoutee_base_included
#define Masse_ajoutee_base_included

#include <TRUSTTabs_forward.h>
#include <Correlation_base.h>

/*! @brief Added-mass correlation of the form:
 *       alpha_k rho_k Dv_k / Dt -> alpha_k rho_k Dv_k / Dt + sum_l ma(k, l) Dv_l / Dt
 *       This class defines a compute function with:
 *     inputs:
 *         alpha[n]  -> volume fraction of phase n
 *         rho[n]    -> density of phase n
 *
 *     input / output:
 *        a_r(k, l)   -> to be placed in the momentum equation (default: alpha(k) * rho(k) for phase k)
 *
 *     NB: the alpha passed to the correlation is the old alpha: no derivative to compute
 */

class Masse_ajoutee_base : public Correlation_base
{
  Declare_base(Masse_ajoutee_base);
public:
  virtual void ajouter(const double *alpha, const double *rho, DoubleTab& a_r  ) const = 0;
  virtual void coefficient(const double *alpha, const double *rho, DoubleTab& coeff) const {Process::exit(que_suis_je() + " : you must define a coefficient function for added mass !");};
  virtual void ajouter_inj(const double *flux_alpha, const double *alpha, const double *rho, DoubleTab& f_a_r) const = 0;
  virtual void coeff(const DoubleTab& alpha, const DoubleTab& rho, DoubleTab& coeff) const  = 0 ;

protected:
  /*! @brief Finds the continuous liquid phase index in a multiphase problem. */
  int find_liquid_phase() const;

  double limiter_liquid_ = 0.5 ; // Maximum percentage of the liquid that can be entrained by the bubbles
};

#endif
