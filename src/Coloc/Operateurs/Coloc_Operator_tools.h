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

#ifndef Coloc_Operator_tools_included
#define Coloc_Operator_tools_included

#include <algorithm>
#include <cmath>

inline void compute_hll_bounds(const DoubleTab& vit_n, const DoubleTab& c, const int f, const int el, const int er, const int nb_phases, double& Sm, double& Sp)
{
  Sm = 0.;
  Sp = 0.;
  for (int n = 0; n < nb_phases; n++)
    {
      const double un_l = vit_n(f, n);
      const double un_r = vit_n(f, n + nb_phases);
      const double c_l = c(el, n), c_r = c(er, n);
      const double Sm_k = std::min(un_l - c_l, un_r - c_r);
      const double Sp_k = std::max(un_l + c_l, un_r + c_r);

      Sm = std::min(Sm, std::min(0.0, Sm_k));
      Sp = std::max(Sp, std::max(0.0, Sp_k));
    }
}

inline double compute_rusanov_speed(const DoubleTab& vit_n, const DoubleTab& c, const int f, const int el, const int er, const int n, const int nb_phases)
{
  const double un_l = vit_n(f, n);
  const double un_r = vit_n(f, n + nb_phases);
  const double c_l = c(el, n), c_r = c(er, n);
  double s = std::max(fabs(un_l - c_l), fabs(un_l + c_l));

  s = std::max(s, std::max(fabs(un_r - c_r), fabs(un_r + c_r)));
  return s;
}

inline void compute_non_conservative_hll_left_bounds(const DoubleTab& vit_n, const DoubleTab& c, const int f, const int el, const int er,
                                                     const int m, const int n, const int nb_phases, double& Sm, double& Sp, double& un)
{
  int k = m; // index pressure
  double un_l = vit_n(f, k);
  double un_r = vit_n(f, k + nb_phases);
  double c_l = c(el, k);
  double c_r = c(er, k);
  const double Sm1 = std::min(un_l - c_l, un_r - c_r);
  const double Sp1 = std::max(un_l + c_l, un_r + c_r);

  k = n; // index velocity
  un_l = vit_n(f, k);
  un_r = vit_n(f, k + nb_phases);
  c_l = c(el, k);
  c_r = c(er, k);
  const double Sm2 = std::min(un_l - c_l, un_r - c_r);
  const double Sp2 = std::max(un_l + c_l, un_r + c_r);

  Sm = std::min(0.0, std::min(Sm1, Sm2));
  Sp = std::max(0.0, std::max(Sp1, Sp2));
  un = un_l;
}

inline void compute_non_conservative_hll_right_bounds(const DoubleTab& vit_n, const DoubleTab& c, const int f, const int el, const int er,
                                                      const int m, const int n, const int nb_phases, double& Sm, double& Sp, double& un)
{
  int k = m; // index pressure
  double un_r = -vit_n(f, k);
  double un_l = -vit_n(f, k + nb_phases);
  double c_l = c(er, k);
  double c_r = c(el, k);
  const double Sm1 = std::min(un_l - c_l, un_r - c_r);
  const double Sp1 = std::max(un_l + c_l, un_r + c_r);

  k = n;  // index velocity
  un_r = -vit_n(f, k);
  un_l = -vit_n(f, k + nb_phases);
  c_l = c(er, k);
  c_r = c(el, k);
  const double Sm2 = std::min(un_l - c_l, un_r - c_r);
  const double Sp2 = std::max(un_l + c_l, un_r + c_r);

  Sm = std::min(0.0, std::min(Sm1, Sm2));
  Sp = std::max(0.0, std::max(Sp1, Sp2));
  un = un_l;
}

#endif /*Coloc_Operator_tools_included*/
