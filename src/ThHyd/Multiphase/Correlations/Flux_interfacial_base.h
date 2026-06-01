/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Flux_interfacial_base_included
#define Flux_interfacial_base_included

#include <Correlation_base.h>
#include <TRUSTTabs.h>

/*! @brief Base class for interfacial heat flux correlations of the form:
 *
 *       Phi_{kl} = h_{kl}(T_k - T_l)
 *       This class defines a flux function with:
 *     inputs:
 *         D_h       -> hydraulic diameter
 *         alpha[n]  -> void fraction of phase n
 *         T[n]      -> temperature of phase n
 *         p         -> pressure
 *         nv[N * n + k]     -> norm of the velocity of phase n if k == n, norm of v_k-v_n otherwise
 *         lambda[n], mu[n], rho[n], Cp[n] -> various physical properties of phase n
 *
 *     outputs:
 *        hi(k, l)    -> heat transfer coefficient between phase k and the interface with phase l (hi(l, k) != hi(k, l) !)
 *     dT_hi(k, l, n) -> derivative of hi(k, l) w.r.t. T[n]
 *     da_hi(k, l, n) -> derivative of hi(k, l) w.r.t. a[n]
 *     dp_hi(k, l)    -> derivative of hi(k, l) w.r.t. p
 *
 *
 */
class Flux_interfacial_base : public Correlation_base
{
  Declare_base(Flux_interfacial_base);
public:
  /* input parameters */
  struct input_t
  {
    double dh;            // hydraulic diameter
    const double *alpha;  // alpha[n] : void fraction of phase n
    const double *T;      // T[n]     : temperature of phase n
    const double *T_passe;// T_passe[n]: temperature of phase n at the previous iteration
    double p;             // pressure
    const double *nv;     // nv[N * k + l] : norm of ||v_k - v_l||
    const double *lambda; // lambda[n]     : thermal conductivity of phase n
    const double *mu;     // mu[n]         : dynamic viscosity of phase n
    const double *rho;    // rho[n]        : density of phase n
    const double *Cp;     // CP[n]         : heat capacity of phase n
    const double *Lvap;   // Lvap[ind_trav]  : latent heat of phase change n=>k, ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1)
    const double *dP_Lvap;//dP_Lvap[ind_trav]: pressure derivative of latent heat, ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1)
    const double *h;      // h[n]          : enthalpy of phase n
    const double *dP_h;   // dP_h[n]       : pressure derivative of enthalpy of phase n
    const double *dT_h;   // dT_h[n]       : temperature derivative of enthalpy of phase n
    const double *d_bulles;//d_bulles[n]   : bubble diameter of phase n
    const double *k_turb; // k_turb[n]     : turbulent kinetic energy of phase n
    const double *nut;    // nut[n]        : turbulent viscosity of phase n
    const double *sigma;  //sigma[ind_trav]: surface tension sigma(ind_trav), ind_trav = (n*(N-1)-(n-1)*(n)/2) + (m-n-1)
    const double *Tsat;   // Tsat[ind_trav]: saturation temperature for phase change n <=> k
    const double *dP_Tsat;//dP_Tsat[ind_trav]: pressure derivative of saturation temperature, ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1)
    DoubleTab v;          // v(n, d)       : velocity of phase n in direction d
    int e;                // element index
  };
  /* output values */
  struct output_t
  {
    DoubleTab hi;    //hi(k, l)       : heat transfer coefficient between phase k and the interface with phase l (hi(l, k) != hi(k, l) !)
    DoubleTab dT_hi; //dT_hi(k, l, n) : derivative of hi(k, l) w.r.t. T[n]
    DoubleTab da_hi; //da_hi(k, l, n) : derivative of hi(k, l) w.r.t. a[n]
    DoubleTab dp_hi; //dp_hi(k, l)    : derivative of hi(k, l) w.r.t. p
  };
  virtual void coeffs(const input_t& input, output_t& output) const = 0;
  double dv_min() const {return dv_min_;};
  double dv_min_ = 0.01;
};

#endif
