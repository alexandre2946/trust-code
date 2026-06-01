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

#ifndef Sch_CN_EX_iteratif_included
#define Sch_CN_EX_iteratif_included

#include <Sch_CN_iteratif.h>
/*! @brief Extended iterative Crank-Nicolson scheme with additional stabilisation tricks.
 *
 *      Extends Sch_CN_iteratif to remain stable beyond its natural stability domain (facsec < 2).
 *
 *      An iteration damping factor omega is introduced. It improves stability but degrades solution quality:
 *      at low iteration counts, time derivatives are underestimated and conservation laws may not be satisfied.
 *
 *      To allow larger time steps, equations other than Navier-Stokes are solved by advancing n explicit
 *      Euler sub-steps. These n sub-steps are recomputed at each iteration of Sch_CN_iteratif.
 *
 *      This scheme is suited to large-scale industrial LES hydraulic simulations with solid-coupled thermics.
 *
 * @sa Sch_CN_iteratif
 */
class Sch_CN_EX_iteratif : public Sch_CN_iteratif
{
  Declare_instanciable(Sch_CN_EX_iteratif);

public :

  void set_param(Param& titi) const override;
  void mettre_a_jour_dt_stab() override;

protected :

  bool iterateTimeStepOnEquation(int i,bool& converged) override;
  virtual bool iterateTimeStepOnNS(int i,bool& converged);
  virtual bool iterateTimeStepOnOther(int i,bool& converged);

  double omega=0.1;
  void ajuster_facsec(type_convergence cv) override;

};

#endif
