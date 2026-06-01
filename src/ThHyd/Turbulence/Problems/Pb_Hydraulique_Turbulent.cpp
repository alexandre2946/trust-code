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

#include <Pb_Hydraulique_Turbulent.h>
#include <Fluide_Incompressible.h>

Implemente_instanciable(Pb_Hydraulique_Turbulent, "Pb_Hydraulique_Turbulent", Pb_Fluide_base);
// XD pb_hydraulique_turbulent Pb_base pb_hydraulique_turbulent INHERITS_BRACE Resolution of Navier-Stokes equations
// XD_CONT with turbulence modelling.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr navier_stokes_turbulent navier_stokes_turbulent navier_stokes_turbulent REQ Navier-Stokes equations as well
// XD_CONT as the associated turbulence model equations.

Sortie& Pb_Hydraulique_Turbulent::printOn(Sortie& os) const { return Pb_Fluide_base::printOn(os); }

Entree& Pb_Hydraulique_Turbulent::readOn(Entree& is) { return Pb_Fluide_base::readOn(is); }

/*! @brief Returns the number of equations.
 *
 *     Returns 1 because there is only 1 equation in a turbulent
 *     hydraulic problem: the Navier-Stokes equation with turbulence.
 *
 * @return Number of equations (1).
 */
int Pb_Hydraulique_Turbulent::nombre_d_equations() const
{
  return 1;
}

const Equation_base& Pb_Hydraulique_Turbulent::equation(int i) const
{
  if (!(i == 0))
    {
      Cerr << "\nError in Pb_Hydraulique_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}

/*! @brief Returns the hydraulic equation of type Navier_Stokes_Turbulent if i=0, exits otherwise.
 *
 * @param i Index of the equation to return.
 * @return The hydraulic equation of type Navier_Stokes_Turbulent.
 */
Equation_base& Pb_Hydraulique_Turbulent::equation(int i)
{
  if (!(i == 0))
    {
      Cerr << "\nError in Pb_Hydraulique_Turbulent::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}

/*! @brief Associates the medium to the problem.
 *
 * The medium must be of type incompressible fluid.
 *
 * @param mil Physical medium to associate with the problem.
 * @throws If the medium is not of type Fluide_Incompressible.
 */
void Pb_Hydraulique_Turbulent::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_Incompressible, mil))
    eq_hydraulique.associer_milieu_base(mil);
  else
    {
      Cerr << "A medium of type " << mil.que_suis_je() << " cannot be associated with " << finl;
      Cerr << "a problem of type Pb_Hydraulique_Turbulent " << finl;
      exit();
    }
}
