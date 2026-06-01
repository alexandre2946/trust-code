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

#ifndef EDO_Pression_th_VEF_included
#define EDO_Pression_th_VEF_included

#include <EDO_Pression_th_base.h>

/*! @brief class EDO_Pression_th_VEF This class represents the ODE on the pressure associated with the computation scheme for
 *
 *      weakly compressible fluids, and relative to
 *      VEF type discretization.
 *
 * @sa Fluide_Quasi_Compressible EDO_Pression_th_base
 */

class EDO_Pression_th_VEF: public EDO_Pression_th_base
{
  Declare_base(EDO_Pression_th_VEF);

public :
  void completer() override;
  void calculer_grad(const DoubleTab&,DoubleTab&);
  double masse_totale(double P,const DoubleTab& T) override;
  double masse_totale(const DoubleTab& P,const DoubleTab& T) override;
};

#endif /* EDO_Pression_th_VEF_included */
