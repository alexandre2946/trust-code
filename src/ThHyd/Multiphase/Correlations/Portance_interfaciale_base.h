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

#ifndef Portance_interfaciale_base_included
#define Portance_interfaciale_base_included
#include <Correlation_base.h>
#include <TRUSTTab.h>

/*! @brief Base class for interfacial lift force operators of the form:
 *
 *       F_{0l} = - F_{l0} = C_{0l} (u_l - u_0) x rot(u_0), where phase
 *       0 is the carrier phase and l != 0 is any other phase.
 *       This class defines a coefficient function C_{kl} depending on:
 *         alpha, p, T -> unknowns (one value per phase each)
 *         rho, mu, sigma -> physical properties (same)
 *         ndv(k, l) -> ||v_k - v_l||, filled for k < l
 *     output:
 *         coeff(k, l, 0/1) -> coefficient C_{kl} and its derivative w.r.t. ndv(k, l), filled for k < l
 *
 *
 */

class Portance_interfaciale_base : public Correlation_base
{
  Declare_base(Portance_interfaciale_base);
public:
  struct input_t
  {
    double dh;            // hydraulic diameter
    DoubleTab alpha;  // alpha[n] : void fraction of phase n
    DoubleTab T;      // T[n]     : temperature of phase n
    DoubleTab p;             // pressure
    DoubleTab nv;     // nv[k, l] : norm of ||v_k - v_l||
    DoubleTab mu;     // mu[n]         : dynamic viscosity of phase n
    DoubleTab rho;    // rho[n]        : density of phase n
    DoubleTab sigma;  // sigma[ind_trav]: surface tension sigma(ind_trav), ind_trav = (n*(N-1)-(n-1)*(n)/2) + (m-n-1)
    DoubleTab k_turb; // k_turb[n]     : turbulent kinetic energy of phase n
    DoubleTab d_bulles;//d_bulles[n]   : bubble diameter of phase n
    int e;                // element index
  };
  /* output values */
  struct output_t
  {
    DoubleTab Cl;    //Cl(k, l)       : lift coefficient between phases k and l
  };

  virtual void coefficient(const input_t& input, output_t& output) const  = 0;

protected:
  /*! @brief Finds the continuous liquid phase index in a multiphase problem. */
  int find_liquid_phase() const;
};

#endif
