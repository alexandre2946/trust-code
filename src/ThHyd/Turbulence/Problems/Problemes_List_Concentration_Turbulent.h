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

#ifndef Problemes_List_Concentration_Turbulent_included
#define Problemes_List_Concentration_Turbulent_included

#include <Convection_Diffusion_Concentration_Turbulent.h>
#include <TRUSTProblem_List_Concentration_Gen.h>
#include <Pb_Thermohydraulique_Turbulent.h>
#include <Pb_Hydraulique_Turbulent.h>

/// \cond DO_NOT_DOCUMENT
class Problemes_List_Concentration_Turbulent
{ /* for source check */ };
/// \endcond

/*! @brief Turbulent thermohydraulics problem with a list of concentration equations.
 *
 *      Couples:
 *      - Turbulent Navier-Stokes equations for an incompressible fluid
 *      - Turbulent energy equation
 *      - Turbulent convection-diffusion equation(s) for one or more concentrations
 *
 * @sa Pb_Fluide_base
 */
class Pb_Thermohydraulique_List_Concentration_Turbulent: public TRUSTProblem_List_Concentration_Gen<Pb_Thermohydraulique_Turbulent, Convection_Diffusion_Concentration_Turbulent, Constituant>
{
  Declare_instanciable(Pb_Thermohydraulique_List_Concentration_Turbulent);
};

/*! @brief Turbulent hydraulics problem with a list of species-transport equations (one or more constituents).
 *
 *      Couples turbulent Navier-Stokes equations for an incompressible fluid
 *      with turbulent convection-diffusion equations for concentrations.
 *
 * @sa Pb_Fluide_base
 */
class Pb_Hydraulique_List_Concentration_Turbulent: public TRUSTProblem_List_Concentration_Gen<Pb_Hydraulique_Turbulent, Convection_Diffusion_Concentration_Turbulent, Constituant>
{
  Declare_instanciable(Pb_Hydraulique_List_Concentration_Turbulent);
};

#endif /* Problemes_List_Concentration_included */
