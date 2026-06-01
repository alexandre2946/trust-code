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

#include <Pb_Hydraulique_Concentration_Turbulent.h>
#include <Fluide_Incompressible.h>
#include <Constituant.h>

Implemente_instanciable(Pb_Hydraulique_Concentration_Turbulent, "Pb_Hydraulique_Concentration_Turbulent", Pb_Fluide_base);
// XD pb_hydraulique_concentration_turbulent Pb_base pb_hydraulique_concentration_turbulent INHERITS_BRACE Resolution of
// XD_CONT Navier-Stokes/multiple constituent transport equations, with turbulence modelling.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_turbulent navier_stokes_turbulent navier_stokes_turbulent OPT Navier-Stokes equations as well
// XD_CONT as the associated turbulence model equations.
// XD attr convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent OPT Constituent transport equations (concentration diffusion convection) as well as the associated turbulence model equations.

Sortie& Pb_Hydraulique_Concentration_Turbulent::printOn(Sortie& os) const { return Pb_Fluide_base::printOn(os); }

Entree& Pb_Hydraulique_Concentration_Turbulent::readOn(Entree& is) { return Pb_Fluide_base::readOn(is); }

/*! @brief Returns the number of equations.
 *
 *     Returns 2 because there are 2 equations in a turbulent hydraulic
 *     problem with transport:
 *       - the turbulent Navier-Stokes equation
 *       - a convection-diffusion equation (possibly vectorial) with turbulence
 *
 * @return Number of equations (2).
 */
int Pb_Hydraulique_Concentration_Turbulent::nombre_d_equations() const
{
  return 2;
}

const Equation_base& Pb_Hydraulique_Concentration_Turbulent::equation(int i) const
{
  if (!(i == 0 || i == 1))
    {
      Cerr << "\nError in Pb_Hydraulique_Concentration_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_concentration;
}

/*! @brief Returns the hydraulic equation of type Navier_Stokes_Turbulent if i=0, returns the convection-diffusion equation of type
 *
 *     Convection_Diffusion_Concentration_Turbulent if i=1
 *     (the convection-diffusion equation may be vectorial).
 *
 * @param i Index of the equation to return.
 * @return The equation corresponding to the given index.
 */
Equation_base& Pb_Hydraulique_Concentration_Turbulent::equation(int i)
{
  if (!(i == 0 || i == 1))
    {
      Cerr << "\nError in Pb_Hydraulique_Concentration_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_concentration;
}

/*! @brief Associates a medium to the problem.
 *
 * Depending on the medium type:
 *       - Fluide_Incompressible: associated to the hydraulic equation
 *       - Constituant: associated to the convection-diffusion equation
 *     Any other medium type causes an error.
 *
 * @param mil Physical medium to associate with the problem.
 * @throws If the medium is not of the correct physical type.
 */
void Pb_Hydraulique_Concentration_Turbulent::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_Incompressible, mil))
    eq_hydraulique.associer_milieu_base(mil);
  else if (sub_type(Constituant, mil))
    eq_concentration.associer_milieu_base(mil);
  else
    {
      Cerr << "A medium of type " << mil.que_suis_je() << " cannot be associated with " << finl;
      Cerr << "a problem of type Pb_Hydraulique_Concentration_Turbulent " << finl;
      exit();
    }
}
