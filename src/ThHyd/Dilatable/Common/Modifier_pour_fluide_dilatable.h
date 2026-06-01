/****************************************************************************
* Copyright (c) 2023, CEA
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

#ifndef Modifier_pour_fluide_dilatable
#define Modifier_pour_fluide_dilatable

/*! @brief Functions intended to multiply or divide an array of values by a vector.
 *
 * These functions are used here to multiply or divide an array by rho when the medium is a Fluide_Dilatable_base.
 *
 * correction_nut_et_cisaillement_paroi_si_qc performs the conversion of the turbulent kinematic viscosity nu_t to
 *    the turbulent dynamic viscosity mu_t for compressible equations.
 *      - This conversion should only be done when using a subgrid-scale model of LES type, since in that case nu_t is returned.
 *      - In RANS simulations, since the correct turbulent quantities are entered multiplied by rho, a mu_t is already obtained
 *         and no conversion is needed.
 *
 */

#include <TRUSTTabs_forward.h>

class Fluide_Dilatable_base;
class Modele_turbulence_hyd_base;
class Milieu_base;

void multiplier_diviser_rho(DoubleVect& tab, const Fluide_Dilatable_base& le_fluide, int diviser = 0);
void diviser_par_rho_si_dilatable(DoubleVect& val, const Milieu_base& mil);
void multiplier_par_rho_si_dilatable(DoubleVect& val, const Milieu_base& mil);
void correction_nut_et_cisaillement_paroi_si_qc(Modele_turbulence_hyd_base& mod);

#endif /* Modifier_pour_fluide_dilatable */
