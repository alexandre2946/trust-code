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

#include <Problemes_Scalaires_Passifs.h>

Implemente_instanciable(Pb_Thermohydraulique_Concentration_Scalaires_Passifs,"Pb_Thermohydraulique_Concentration_Scalaires_Passifs",TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Concentration>);
Sortie& Pb_Thermohydraulique_Concentration_Scalaires_Passifs::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Concentration>::printOn(os); }
Entree& Pb_Thermohydraulique_Concentration_Scalaires_Passifs::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Concentration>::readOn(is); }

Implemente_instanciable(Pb_Hydraulique_Concentration_Scalaires_Passifs,"Pb_Hydraulique_Concentration_Scalaires_Passifs",TRUSTProblem_sup_eqns<Pb_Hydraulique_Concentration>);
Sortie& Pb_Hydraulique_Concentration_Scalaires_Passifs::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Hydraulique_Concentration>::printOn(os); }
Entree& Pb_Hydraulique_Concentration_Scalaires_Passifs::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Hydraulique_Concentration>::readOn(is); }

Implemente_instanciable(Pb_Thermohydraulique_Scalaires_Passifs,"Pb_Thermohydraulique_Scalaires_Passifs",TRUSTProblem_sup_eqns<Pb_Thermohydraulique>);
Sortie& Pb_Thermohydraulique_Scalaires_Passifs::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique>::printOn(os); }
Entree& Pb_Thermohydraulique_Scalaires_Passifs::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique>::readOn(is); }

Implemente_instanciable(Pb_Conduction_Scalaires_Passifs,"Pb_Conduction_Scalaires_Passifs",TRUSTProblem_sup_eqns<Pb_Conduction>);
Sortie& Pb_Conduction_Scalaires_Passifs::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Conduction>::printOn(os); }
Entree& Pb_Conduction_Scalaires_Passifs::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Conduction>::readOn(is); }

// XD listeqn listobj nul BRACE eqn_base NO_COMMA List of equations.

// XD pb_thermohydraulique_concentration_scalaires_passifs Pb_base pb_thermohydraulique_concentration_scalaires_passifs INHERITS_BRACE Resolution of Navier-Stokes/energy/multiple constituent transport equations, with the additional passive scalar equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr convection_diffusion_concentration convection_diffusion_concentration convection_diffusion_concentration OPT
// XD_CONT Constituent transport equations (concentration diffusion convection).
// XD attr convection_diffusion_temperature convection_diffusion_temperature convection_diffusion_temperature OPT Energy
// XD_CONT equations (temperature diffusion convection).
// XD attr equations_scalaires_passifs listeqn equations_scalaires_passifs REQ Passive scalar equations. The unknowns of
// XD_CONT the passive scalar equation number N are named temperatureN or concentrationN or fraction_massiqueN. This
// XD_CONT keyword is used to define initial conditions and the post processing fields. This kind of problem is very
// XD_CONT useful to test in only one data file (and then only one calculation) different schemes or different boundary
// XD_CONT conditions for the scalar transport equation.

// XD pb_thermohydraulique_scalaires_passifs Pb_base pb_thermohydraulique_scalaires_passifs INHERITS_BRACE Resolution of
// XD_CONT thermohydraulic problem, with the additional passive scalar equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr convection_diffusion_temperature convection_diffusion_temperature convection_diffusion_temperature OPT Energy
// XD_CONT equations (temperature diffusion convection).
// XD attr equations_scalaires_passifs listeqn equations_scalaires_passifs REQ Passive scalar equations. The unknowns of
// XD_CONT the passive scalar equation number N are named temperatureN or concentrationN or fraction_massiqueN. This
// XD_CONT keyword is used to define initial conditions and the post processing fields. This kind of problem is very
// XD_CONT useful to test in only one data file (and then only one calculation) different schemes or different boundary
// XD_CONT conditions for the scalar transport equation.

// XD pb_hydraulique_concentration_scalaires_passifs Pb_base pb_hydraulique_concentration_scalaires_passifs INHERITS_BRACE Resolution of Navier-Stokes/multiple constituent transport equations with the additional passive scalar equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr convection_diffusion_concentration convection_diffusion_concentration convection_diffusion_concentration OPT
// XD_CONT Constituent transport equations (concentration diffusion convection).
// XD attr equations_scalaires_passifs listeqn equations_scalaires_passifs REQ Passive scalar equations. The unknowns of
// XD_CONT the passive scalar equation number N are named temperatureN or concentrationN or fraction_massiqueN. This
// XD_CONT keyword is used to define initial conditions and the post processing fields. This kind of problem is very
// XD_CONT useful to test in only one data file (and then only one calculation) different schemes or different boundary
// XD_CONT conditions for the scalar transport equation.

