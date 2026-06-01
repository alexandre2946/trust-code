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

#ifndef Flux_parietal_base_included
#define Flux_parietal_base_included

#include <TRUSTTabs_forward.h>
#include <Correlation_base.h>

/*! @brief Base class for wall heat flux correlations of the form:
 *
 *         sensible heat flux  : q_{p}(k)     = F(alpha_f, p, T_f, T_p, v_f, D_h, D_ch)
 *         latent heat flux    : q_{pi}(k, l) = F(alpha_f, p, T_f, T_p, v_f, D_h, D_ch)
 *           (e.g. nucleate boiling: Gamma_{kl} = q_{pi}(k, l) / Lvap)
 *       This class defines two functions q_pk, q_pi.
 */

class Flux_parietal_base : public Correlation_base
{
  Declare_base(Flux_parietal_base);
public:

  /* input parameters */
  struct input_t
  {
    int N;                // number of phases
    int f;                // face number
    double y;             // distance between the face and the center of gravity of the cell
    double D_h;           // hydraulic diameter
    double D_ch;          // heated hydraulic diameter
    double p;             // pressure
    double Tp;            // wall temperature (single value)
    const double *alpha;  // alpha[n]  -> void fraction of phase n
    const double *T;      // T[n]      -> temperature of phase n
    const double *v;      // v[n]      -> norm of the velocity of phase n
    const double *lambda; // lambda[n] -> thermal conductivity of phase n
    const double *mu;     // mu[n]     -> dynamic viscosity of phase n
    const double *rho;    // rho[n]    -> density of phase n
    const double *Cp;     // Cp[n]     -> heat capacity of phase n
    const double *Lvap;   //Lvap[(k*(N-1)-(k-1)*(k)/2) + (l-k-1)]     : latent heat of phase change n <=> k
    const double *Sigma;  //Sigma[(k*(N-1)-(k-1)*(k)/2) + (l-k-1)]    : surface tension between phases n <=> k
    const double *Tsat;   //Tsat[(k*(N-1)-(k-1)*(k)/2) + (l-k-1)]     : saturation temperature for phase change n <=> k
  };
  /* output values */
  struct output_t
  {
    DoubleTab *qpk = nullptr;     // (*qpk)(n)           -> heat flux toward phase n
    DoubleTab *da_qpk = nullptr;  // (*da_qpk)(n, m)     -> derivative w.r.t. alpha_m
    DoubleTab *dp_qpk = nullptr;  // (*dp_qpk)(n)        -> derivative w.r.t. p
    DoubleTab *dv_qpk = nullptr;  // (*dv_qpk)(n, m)     -> derivative w.r.t. v[m]
    DoubleTab *dTf_qpk = nullptr; // (*dTf_qpk)(n, m)    -> derivative w.r.t. T[m]
    DoubleTab *dTp_qpk = nullptr; // (*dTp_qpk)(n)       -> derivative w.r.t. Tp
    DoubleTab *qpi = nullptr;     // (*qpi)(k, l)        -> heat flux supplied to the phase change from k to l (to fill for k < l)
    DoubleTab *da_qpi = nullptr;  // (*da_qpi)(k, l, m)  -> derivative w.r.t. alpha_m
    DoubleTab *dp_qpi = nullptr;  // (*dp_qpi)(k, l)     -> derivative w.r.t. p
    DoubleTab *dv_qpi = nullptr;  // (*dv_qpi)(k, l,m)   -> derivative w.r.t. v[m]
    DoubleTab *dTf_qpi = nullptr; // (*dTf_qpi)(k, l, m) -> derivative w.r.t. T[m]
    DoubleTab *dTp_qpi = nullptr; // (*dTp_qpi)(k, l)    -> derivative w.r.t. Tp
    DoubleTab *d_nuc = nullptr;   // (*d_nucleation)(k)  -> nucleation diameter of phase k
    int *nonlinear = nullptr;     // nonlinear           -> set to 1 if q_pk / q_pi is nonlinear in Tp / Tf; do not modify otherwise
  };
  virtual void qp(const input_t& input, output_t& output) const = 0;
  /* 1 if T[n] must be provided at the wall, 0 if it must be provided at the cell center */
  virtual int T_at_wall() const = 0;
  virtual int calculates_bubble_nucleation_diameter() const {return 0;};
  virtual int needs_saturation() const {return 0;};
};

#endif
