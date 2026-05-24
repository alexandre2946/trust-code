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

#include <Problemes_List_Concentration.h>

Implemente_instanciable(Pb_Thermohydraulique_List_Concentration, "Pb_Thermohydraulique_List_Concentration", TRUSTProblem_List_Concentration_Gen<Pb_Thermohydraulique>);
Sortie& Pb_Thermohydraulique_List_Concentration::printOn(Sortie& os) const { return TRUSTProblem_List_Concentration_Gen<Pb_Thermohydraulique>::printOn(os); }
Entree& Pb_Thermohydraulique_List_Concentration::readOn(Entree& is) { return TRUSTProblem_List_Concentration_Gen<Pb_Thermohydraulique>::readOn(is); }

Implemente_instanciable(Pb_Hydraulique_List_Concentration, "Pb_Hydraulique_List_Concentration", TRUSTProblem_List_Concentration_Gen<Pb_Hydraulique>);
Sortie& Pb_Hydraulique_List_Concentration::printOn(Sortie& os) const { return TRUSTProblem_List_Concentration_Gen<Pb_Hydraulique>::printOn(os); }
Entree& Pb_Hydraulique_List_Concentration::readOn(Entree& is) { return TRUSTProblem_List_Concentration_Gen<Pb_Hydraulique>::readOn(is); }

// XD Pb_Thermohydraulique_List_Concentration Pb_base Pb_Thermohydraulique_List_Concentration INHERITS_BRACE Resolution
// XD_CONT of Navier-Stokes/energy/multiple constituent transport equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr convection_diffusion_temperature convection_diffusion_temperature convection_diffusion_temperature OPT Energy
// XD_CONT equation (temperature diffusion convection).
// XD attr list_equations listeqn list_equations REQ convection_diffusion_concentration equations. The unknown of the
// XD_CONT concentration equation number N is named concentrationN. This keyword is used to define initial conditions
// XD_CONT and the post processing fields. This kind of problem is very useful to test in only one data file (and then
// XD_CONT only one calculation) different schemes or different boundary conditions for the scalar transport equation.

// XD Pb_Hydraulique_List_Concentration Pb_base Pb_Hydraulique_List_Concentration INHERITS_BRACE Resolution of
// XD_CONT Navier-Stokes/multiple constituent transport equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr list_equations listeqn list_equations REQ convection_diffusion_concentration equations. The unknown of the
// XD_CONT concentration equation number N is named concentrationN. This keyword is used to define initial conditions
// XD_CONT and the post processing fields. This kind of problem is very useful to test in only one data file (and then
// XD_CONT only one calculation) different schemes or different boundary conditions for the scalar transport equation.
