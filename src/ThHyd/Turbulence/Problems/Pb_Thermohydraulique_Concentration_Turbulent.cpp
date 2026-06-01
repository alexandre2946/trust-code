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

#include <Pb_Thermohydraulique_Concentration_Turbulent.h>
#include <Fluide_Incompressible.h>
#include <Constituant.h>

Implemente_instanciable(Pb_Thermohydraulique_Concentration_Turbulent, "Pb_Thermohydraulique_Concentration_Turbulent", Pb_Fluide_base);
// XD pb_thermohydraulique_concentration_turbulent Pb_base pb_thermohydraulique_concentration_turbulent INHERITS_BRACE
// XD_CONT Resolution of Navier-Stokes/energy/multiple constituent transport equations, with turbulence modelling.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_turbulent navier_stokes_turbulent navier_stokes_turbulent OPT Navier-Stokes equations as well
// XD_CONT as the associated turbulence model equations.
// XD attr convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent OPT Constituent transport equations (concentration diffusion convection) as well as the associated turbulence model equations.
// XD attr convection_diffusion_temperature_turbulent convection_diffusion_temperature_turbulent convection_diffusion_temperature_turbulent OPT Energy equation (temperature diffusion convection) as well as the associated turbulence model equations.

Sortie& Pb_Thermohydraulique_Concentration_Turbulent::printOn(Sortie& os) const { return Pb_Fluide_base::printOn(os); }

Entree& Pb_Thermohydraulique_Concentration_Turbulent::readOn(Entree& is) { return Pb_Fluide_base::readOn(is); }

/*! @brief Returns the number of equations.
 *
 *     Returns 3 because there are 3 equations in a turbulent
 *     thermohydraulic problem with concentration:
 *         - the turbulent Navier-Stokes equation
 *         - the turbulent energy equation
 *         - a turbulent convection-diffusion equation
 *
 * @return Number of equations (3).
 */
int Pb_Thermohydraulique_Concentration_Turbulent::nombre_d_equations() const
{
  return 3;
}

const Equation_base& Pb_Thermohydraulique_Concentration_Turbulent::equation(int i) const
{
  if (!(i == 0 || i == 1 || i == 2))
    {
      Cerr << "\nError in Pb_Thermohydraulique_Concentration_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else if (i == 1)
    return eq_thermique;
  else
    return eq_concentration;
}

/*! @brief Returns the hydraulic equation of type Navier_Stokes_Turbulent if i=0, returns the thermal equation of type Convection_Diffusion_Temperature_Turbulent if i=1, returns the concentration equation of type Convection_Diffusion_Concentration_Turbulent if i=2.
 *
 * @param i Index of the equation to return.
 * @return The equation corresponding to the given index.
 */
Equation_base& Pb_Thermohydraulique_Concentration_Turbulent::equation(int i)
{
  if (!(i == 0 || i == 1 || i == 2))
    {
      Cerr << "\nError in Pb_Thermohydraulique_Concentration_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else if (i == 1)
    return eq_thermique;
  else
    return eq_concentration;
}

/*! @brief Associates a medium to the problem.
 *
 * Depending on the medium type:
 *       - Fluide_Incompressible: associated to the hydraulic equation and the energy equation
 *       - Constituant: associated to the convection-diffusion equation
 *     Any other medium type causes an error.
 *
 * @param mil Physical medium to associate with the problem.
 * @throws If the medium is not of the correct physical type.
 */
void Pb_Thermohydraulique_Concentration_Turbulent::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_Incompressible, mil))
    {
      eq_hydraulique.associer_milieu_base(mil);
      eq_thermique.associer_milieu_base(mil);
    }
  else if (sub_type(Constituant, mil))
    eq_concentration.associer_milieu_base(mil);
  else
    {
      Cerr << "A medium of type " << mil.que_suis_je() << " cannot be associated with " << finl;
      Cerr << "a problem of type Pb_Thermohydraulique_Concentration_Turbulent " << finl;
      exit();
    }
}
