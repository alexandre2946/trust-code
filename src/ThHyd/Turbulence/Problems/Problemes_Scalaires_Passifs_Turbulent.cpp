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

#include <Problemes_Scalaires_Passifs_Turbulent.h>

Implemente_instanciable(Pb_Thermohydraulique_Turbulent_Scalaires_Passifs,"Pb_Thermohydraulique_Turbulent_Scalaires_Passifs",TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Turbulent>);
Sortie& Pb_Thermohydraulique_Turbulent_Scalaires_Passifs::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Turbulent>::printOn(os); }
Entree& Pb_Thermohydraulique_Turbulent_Scalaires_Passifs::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Turbulent>::readOn(is); }

Implemente_instanciable(Pb_Thermohydraulique_Especes_Turbulent_QC,"Pb_Thermohydraulique_Especes_Turbulent_QC",TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Turbulent_QC>);
Sortie& Pb_Thermohydraulique_Especes_Turbulent_QC::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Turbulent_QC>::printOn(os); }
Entree& Pb_Thermohydraulique_Especes_Turbulent_QC::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Turbulent_QC>::readOn(is); }

Implemente_instanciable(Pb_Hydraulique_Concentration_Turbulent_Scalaires_Passifs,"Pb_Hydraulique_Concentration_Turbulent_Scalaires_Passifs",TRUSTProblem_sup_eqns<Pb_Hydraulique_Concentration_Turbulent>);
Sortie& Pb_Hydraulique_Concentration_Turbulent_Scalaires_Passifs::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Hydraulique_Concentration_Turbulent>::printOn(os); }
Entree& Pb_Hydraulique_Concentration_Turbulent_Scalaires_Passifs::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Hydraulique_Concentration_Turbulent>::readOn(is); }

Implemente_instanciable(Pb_Thermohydraulique_Concentration_Turbulent_Scalaires_Passifs,"Pb_Thermohydraulique_Concentration_Turbulent_Scalaires_Passifs",TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Concentration_Turbulent>);
Sortie& Pb_Thermohydraulique_Concentration_Turbulent_Scalaires_Passifs::printOn(Sortie& os) const { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Concentration_Turbulent>::printOn(os); }
Entree& Pb_Thermohydraulique_Concentration_Turbulent_Scalaires_Passifs::readOn(Entree& is) { return TRUSTProblem_sup_eqns<Pb_Thermohydraulique_Concentration_Turbulent>::readOn(is); }

// XD pb_thermohydraulique_turbulent_scalaires_passifs Pb_base pb_thermohydraulique_turbulent_scalaires_passifs INHERITS_BRACE Resolution of thermohydraulic problem, with turbulence modelling and with the additional passive scalar equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_turbulent navier_stokes_turbulent navier_stokes_turbulent OPT Navier-Stokes equations as well
// XD_CONT as the associated turbulence model equations.
// XD attr convection_diffusion_temperature_turbulent convection_diffusion_temperature_turbulent convection_diffusion_temperature_turbulent OPT Energy equations (temperature diffusion convection) as well as the associated turbulence model equations.
// XD attr equations_scalaires_passifs listeqn equations_scalaires_passifs REQ Passive scalar equations. The unknowns of
// XD_CONT the passive scalar equation number N are named temperatureN or concentrationN or fraction_massiqueN. This
// XD_CONT keyword is used to define initial conditions and the post processing fields. This kind of problem is very
// XD_CONT useful to test in only one data file (and then only one calculation) different schemes or different boundary
// XD_CONT conditions for the scalar transport equation.

// XD pb_thermohydraulique_especes_turbulent_qc Pb_base pb_thermohydraulique_especes_turbulent_qc INHERITS_BRACE
// XD_CONT Resolution of turbulent thermohydraulic problem under low Mach number with passive scalar equations.
// XD attr fluide_quasi_compressible fluide_quasi_compressible fluide_quasi_compressible REQ The fluid medium associated
// XD_CONT with the problem.
// XD attr navier_stokes_turbulent_qc navier_stokes_turbulent_qc navier_stokes_turbulent_qc REQ Navier-Stokes equations
// XD_CONT under low Mach number as well as the associated turbulence model equations.
// XD attr convection_diffusion_chaleur_turbulent_qc convection_diffusion_chaleur_turbulent_qc convection_diffusion_chaleur_turbulent_qc REQ Energy equation under low Mach number as well as the associated turbulence model equations.
// XD attr equations_scalaires_passifs listeqn equations_scalaires_passifs REQ Passive scalar equations. The unknowns of
// XD_CONT the passive scalar equation number N are named temperatureN or concentrationN or fraction_massiqueN. This
// XD_CONT keyword is used to define initial conditions and the post processing fields. This kind of problem is very
// XD_CONT useful to test in only one data file (and then only one calculation) different schemes or different boundary
// XD_CONT conditions for the scalar transport equation.

// XD pb_hydraulique_concentration_turbulent_scalaires_passifs Pb_base pb_hydraulique_concentration_turbulent_scalaires_passifs INHERITS_BRACE Resolution of Navier-Stokes/multiple constituent transport equations, with turbulence modelling and with the additional passive scalar equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_turbulent navier_stokes_turbulent navier_stokes_turbulent OPT Navier-Stokes equations as well
// XD_CONT as the associated turbulence model equations.
// XD attr convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent OPT Constituent transport equations (concentration diffusion convection) as well as the associated turbulence model equations.
// XD attr equations_scalaires_passifs listeqn equations_scalaires_passifs REQ Passive scalar equations. The unknowns of
// XD_CONT the passive scalar equation number N are named temperatureN or concentrationN or fraction_massiqueN. This
// XD_CONT keyword is used to define initial conditions and the post processing fields. This kind of problem is very
// XD_CONT useful to test in only one data file (and then only one calculation) different schemes or different boundary
// XD_CONT conditions for the scalar transport equation.

// XD pb_thermohydraulique_concentration_turbulent_scalaires_passifs Pb_base pb_thermohydraulique_concentration_turbulent_scalaires_passifs INHERITS_BRACE Resolution of Navier-Stokes/energy/multiple constituent transport equations, with turbulence modelling and with the additional passive scalar equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_turbulent navier_stokes_turbulent navier_stokes_turbulent OPT Navier-Stokes equations as well
// XD_CONT as the associated turbulence model equations.
// XD attr convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent convection_diffusion_concentration_turbulent OPT Constituent transport equations (concentration diffusion convection) as well as the associated turbulence model equations.
// XD attr convection_diffusion_temperature_turbulent convection_diffusion_temperature_turbulent convection_diffusion_temperature_turbulent OPT Energy equations (temperature diffusion convection) as well as the associated turbulence model equations.
// XD attr equations_scalaires_passifs listeqn equations_scalaires_passifs REQ Passive scalar equations. The unknowns of
// XD_CONT the passive scalar equation number N are named temperatureN or concentrationN or fraction_massiqueN. This
// XD_CONT keyword is used to define initial conditions and the post processing fields. This kind of problem is very
// XD_CONT useful to test in only one data file (and then only one calculation) different schemes or different boundary
// XD_CONT conditions for the scalar transport equation.
