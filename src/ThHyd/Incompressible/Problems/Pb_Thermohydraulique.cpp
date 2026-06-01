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

#include <Pb_Thermohydraulique.h>
#include <Fluide_Ostwald.h>
#include <Verif_Cl.h>
#include <Champ_Uniforme.h>

Implemente_instanciable(Pb_Thermohydraulique, "Pb_Thermohydraulique", Pb_Hydraulique);
// XD pb_thermohydraulique Pb_base pb_thermohydraulique INHERITS_BRACE Resolution of thermohydraulic problem.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible OPT The fluid medium associated with the
// XD_CONT problem (only one possibility).
// XD attr fluide_ostwald fluide_ostwald fluide_ostwald OPT The fluid medium associated with the problem (only one
// XD_CONT possibility).
// XD attr fluide_sodium_liquide fluide_sodium_liquide fluide_sodium_liquide OPT The fluid medium associated with the
// XD_CONT problem (only one possibility).
// XD attr fluide_sodium_gaz fluide_sodium_gaz fluide_sodium_gaz OPT The fluid medium associated with the problem (only
// XD_CONT one possibility).
// XD attr correlations bloc_lecture correlations OPT List of correlations used in specific source terms (i.e.
// XD_CONT interfacial flux, interfacial friction, ...)
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr convection_diffusion_temperature convection_diffusion_temperature convection_diffusion_temperature OPT Energy
// XD_CONT equation (temperature diffusion convection).

Sortie& Pb_Thermohydraulique::printOn(Sortie& os) const { return Pb_Hydraulique::printOn(os); }
Entree& Pb_Thermohydraulique::readOn(Entree& is) { return Pb_Hydraulique::readOn(is); }

/*! @brief Returns the hydraulic equation of type Navier_Stokes_std if i=0, returns the thermal equation of type
 *
 *     Convection_Diffusion_Temperature if i=1
 *     (const version)
 *
 * @param i the index of the equation to return
 * @return the equation corresponding to the index
 */
const Equation_base& Pb_Thermohydraulique::equation(int i) const
{
  if (i == 0) return eq_hydraulique;
  else if (i == 1) return eq_thermique;
  else if (i < 2 + eq_opt_.size() && i > 1) return eq_opt_[i - 2].valeur();
  else
    {
      Cerr << "Pb_Thermohydraulique::equation() : Wrong equation number" << i << "!" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}

/*! @brief Returns the hydraulic equation of type Navier_Stokes_std if i=0, returns the thermal equation of type
 *
 *     Convection_Diffusion_Temperature if i=1
 *
 * @param i the index of the equation to return
 * @return the equation corresponding to the index
 */
Equation_base& Pb_Thermohydraulique::equation(int i)
{
  if (i == 0) return eq_hydraulique;
  else if (i == 1) return eq_thermique;
  else if (i < 2 + eq_opt_.size() && i > 1) return eq_opt_[i - 2].valeur();
  else
    {
      Cerr << "Pb_Thermohydraulique::equation() : Wrong equation number" << i << "!" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}

/*! @brief Associates the medium to the problem. The medium must be of type incompressible fluid.
 *
 * @param mil the physical medium to associate with the problem
 * @throws wrong type of physical medium
 */
void Pb_Thermohydraulique::associer_milieu_base(const Milieu_base& mil)
{
  Pb_Hydraulique::associer_milieu_base(mil);
  if (sub_type(Fluide_base,mil) && ref_cast(Fluide_base, mil).is_incompressible())
    eq_thermique.associer_milieu_base(mil);
  else if (sub_type(Fluide_Ostwald,mil))
    eq_thermique.associer_milieu_base(mil);
}

/*! @brief Tests the compatibility of the thermal and hydraulic equations.
 *
 * The test is performed on the discretized boundary conditions of each equation.
 *     Calls the library function:
 *       tester_compatibilite_hydr_thermique(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return propagated return code
 */
int Pb_Thermohydraulique::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_th = eq_thermique.domaine_Cl_dis();
  return tester_compatibilite_hydr_thermique(domaine_Cl_hydr,domaine_Cl_th);
}



