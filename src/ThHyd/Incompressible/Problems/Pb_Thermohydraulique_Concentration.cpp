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

#include <Pb_Thermohydraulique_Concentration.h>
#include <Fluide_Incompressible.h>
#include <Constituant.h>
#include <Verif_Cl.h>

Implemente_instanciable(Pb_Thermohydraulique_Concentration, "Pb_Thermohydraulique_Concentration", Pb_Thermohydraulique);
// XD pb_thermohydraulique_concentration Pb_base pb_thermohydraulique_concentration INHERITS_BRACE Resolution of
// XD_CONT Navier-Stokes/energy/multiple constituent transport equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr convection_diffusion_concentration convection_diffusion_concentration convection_diffusion_concentration OPT
// XD_CONT Constituent transport equations (concentration diffusion convection).
// XD attr convection_diffusion_temperature convection_diffusion_temperature convection_diffusion_temperature OPT Energy
// XD_CONT equation (temperature diffusion convection).

Sortie& Pb_Thermohydraulique_Concentration::printOn(Sortie& os) const { return Pb_Thermohydraulique::printOn(os); }
Entree& Pb_Thermohydraulique_Concentration::readOn(Entree& is) { return Pb_Thermohydraulique::readOn(is); }

/*! @brief Returns the hydraulic equation of type Navier_Stokes_std if i=0, returns the thermal equation of type
 *
 *     Convection_Diffusion_Temperature if i=1,
 *     returns the concentration equation of type
 *     Convection_Diffusion_Concentration if i=2
 *     (const version)
 *
 * @param i the index of the equation to return
 * @return the equation corresponding to the index
 */
const Equation_base& Pb_Thermohydraulique_Concentration::equation(int i) const
{
  if (i == 2) return eq_concentration;
  return Pb_Thermohydraulique::equation(i);
}

/*! @brief Returns the hydraulic equation of type Navier_Stokes_std if i=0, returns the thermal equation of type
 *
 *     Convection_Diffusion_Temperature if i=1,
 *     returns the concentration equation of type
 *     Convection_Diffusion_Concentration if i=2
 *
 * @param i the index of the equation to return
 * @return the equation corresponding to the index
 */
Equation_base& Pb_Thermohydraulique_Concentration::equation(int i)
{
  if (i == 2) return eq_concentration;
  return Pb_Thermohydraulique::equation(i);
}

/*! @brief Associates a medium to the problem. If the medium is of type:
 *
 *       - Fluide_Incompressible, it will be associated with the hydraulic equation
 *         and the energy equation.
 *       - Constituant, it will be associated with the convection-diffusion equation
 *     Any other medium type causes an error.
 *
 * @param mil the physical medium to associate with the problem
 * @throws wrong type of physical medium
 */
void Pb_Thermohydraulique_Concentration::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Constituant, mil))
    eq_concentration.associer_milieu_base(mil);
  else
    Pb_Thermohydraulique::associer_milieu_base(mil);
}

/*! @brief Tests the compatibility of the convection-diffusion and hydraulic equations.
 *
 * The test is performed on the discretized boundary conditions of each equation.
 *     Calls the library functions:
 *       tester_compatibilite_hydr_thermique(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *       tester_compatibilite_hydr_concentration(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return propagated return code
 */
int Pb_Thermohydraulique_Concentration::verifier()
{
  Pb_Thermohydraulique::verifier();
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_co = eq_concentration.domaine_Cl_dis();
  return tester_compatibilite_hydr_concentration(domaine_Cl_hydr, domaine_Cl_co);
}
