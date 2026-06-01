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

#include <Pb_Hydraulique_Melange_Binaire_WC.h>

Implemente_instanciable(Pb_Hydraulique_Melange_Binaire_WC,"Pb_Hydraulique_Melange_Binaire_WC",Pb_WC_base);
// XD pb_hydraulique_melange_binaire_WC Pb_base pb_hydraulique_melange_binaire_WC INHERITS_BRACE Resolution of a binary
// XD_CONT mixture problem for a weakly-compressible fluid with an iso-thermal condition. NL2 Keywords for the unknowns
// XD_CONT other than pressure, velocity, fraction_massique are : NL2 masse_volumique : density NL2 pression : reduced
// XD_CONT pressure NL2 pression_tot : total pressure NL2 pression_hydro : hydro-static pressure NL2 pression_eos :
// XD_CONT pressure used in state equation.
// XD attr fluide_weakly_compressible fluide_weakly_compressible fluide_weakly_compressible REQ The fluid medium
// XD_CONT associated with the problem.
// XD attr navier_stokes_WC navier_stokes_WC navier_stokes_WC REQ Navier-Stokes equation for a weakly-compressible
// XD_CONT fluid.
// XD attr convection_diffusion_espece_binaire_WC convection_diffusion_espece_binaire_WC convection_diffusion_espece_binaire_WC REQ Species conservation equation for a binary weakly-compressible fluid.

Sortie& Pb_Hydraulique_Melange_Binaire_WC::printOn(Sortie& os) const { return Probleme_base::printOn(os); }

Entree& Pb_Hydraulique_Melange_Binaire_WC::readOn(Entree& is) { return Probleme_base::readOn(is); }

/*! @brief Returns 2 because there are 2 equations: Navier_Stokes_WC and Convection_Diffusion_Espece_Binaire_WC.
 *
 */
int Pb_Hydraulique_Melange_Binaire_WC::nombre_d_equations() const { return 2; }

/*! @brief Returns the hydraulic equation of type Navier_Stokes_WC if i=0. Returns the mass fraction convection/diffusion equation of type Convection_Diffusion_Espece_Binaire_WC if i=1.
 *
 */
const Equation_base& Pb_Hydraulique_Melange_Binaire_WC::equation(int i) const
{
  return equation_impl(i,eq_hydraulique,eq_frac_mass);
}

Equation_base& Pb_Hydraulique_Melange_Binaire_WC::equation(int i)
{
  return equation_impl(i,eq_hydraulique,eq_frac_mass);
}

/*! @brief Tests the compatibility of the mass fraction and hydraulic equations.
 *
 */
int Pb_Hydraulique_Melange_Binaire_WC::verifier()
{
  return verifier_impl(eq_hydraulique,eq_frac_mass,false /* is_thermal */);
}
