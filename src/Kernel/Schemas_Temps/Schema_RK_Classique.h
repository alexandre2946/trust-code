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

#ifndef Schema_RK_Classique_included
#define Schema_RK_Classique_included

#include <TRUSTSchema_RK.h>

/// \cond DO_NOT_DOCUMENT
class Schema_RK_Classique
{ };
/// \endcond

/*! @brief : class RK2_Classique This class represents a classical second-order Runge-Kutta time scheme:
 *
 *      k1 = h * f(y0)
 *      k2 = h * f(y0 + 0.5 * k1)
 *      y1 = y0 + k2
 *
 */
class RK2_Classique: public TRUSTSchema_RK<Ordre_RK::DEUX_CLASSIQUE>
{
  Declare_instanciable(RK2_Classique);
};

/*! @brief : class RK3_Classique This class represents a classical third-order Runge-Kutta time scheme:
 *
 *      k1 = h * f(y0)
 *      k2 = h * f(y0 + 0.5 * k1)
 *      k3 = h * f(y0 - k1 + 2 * k2)
 *      y1 = y0 + 1/6 * (k1 + 4 * k2 + k3)
 *
 */
class RK3_Classique: public TRUSTSchema_RK<Ordre_RK::TROIS_CLASSIQUE>
{
  Declare_instanciable(RK3_Classique);
};

/*! @brief : class RK4_Classique This class represents a classical fourth-order Runge-Kutta time scheme:
 *
 *      k1 = h * f(y0)
 *      k2 = h * f(y0 + 0.5 * k1)
 *      k3 = h * f(y0 + 0.5 * k2)
 *      k4 = h * f(y0 + k3)
 *      y1 = y0 + 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
 *
 */
class RK4_Classique: public TRUSTSchema_RK<Ordre_RK::QUATRE_CLASSIQUE>
{
  Declare_instanciable(RK4_Classique);
};

/*! @brief : class RK4_Classique_3_8 This class represents a classical fourth-order Runge-Kutta time scheme with 3/8 rule coefficients:
 *
 *      k1 = h * f(y0)
 *      k2 = h * f(y0 + 1/3 * k1)
 *      k3 = h * f(y0 - 1/3 * k1 + k2)
 *      k4 = h * f(y0 + k1 - k2 + k3)
 *      y1 = y0 + 1/8 * (k1 + 3 * k2 + 3 * k3 + k4)
 *
 */
class RK4_Classique_3_8: public TRUSTSchema_RK<Ordre_RK::QUATRE_CLASSIQUE_3_8>
{
  Declare_instanciable(RK4_Classique_3_8);
};

#endif /* Schema_RK_Classique_included */
