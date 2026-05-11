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

#include <Interpolation_IBM_thermal_wall_law.h>

Implemente_instanciable( Interpolation_IBM_thermal_wall_law, "Interpolation_IBM_thermal_wall_law|IBM_thermal_wall_law", Interpolation_IBM_elem_fluid ) ;
// XD interpolation_ibm_thermal_wall_law interpolation_ibm_elem_fluid ibm_thermal_wall_law BRACE Immersed Boundary Method (IBM): Interpolation thermal wall law.

Sortie& Interpolation_IBM_thermal_wall_law::printOn( Sortie& os ) const
{
  return Interpolation_IBM_elem_fluid::printOn( os );
}

Entree& Interpolation_IBM_thermal_wall_law::readOn( Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Interpolation_IBM_thermal_wall_law::set_param(Param& param) const
{
  Interpolation_IBM_elem_fluid::set_param( param );
  param.ajouter("formulation_Tplus",&formulation_Tp_,Param::OPTIONAL);  // XD_ADD_P entier Choix formulation calcul T+ (Kader si rien)
  param.ajouter("boundary_type",&boundary_type_,Param::REQUIRED); // XD_ADD_P entier Choix type de cond limite
  param.ajouter("T_inlet",&T_inlet_,Param::REQUIRED); // XD_ADD_P double Demande la température d'écoulement d'entrée moyenne
  param.ajouter("Prandlt_mol",&Prandlt_mol_,Param::REQUIRED); // XD_ADD_P double Demande le Prandlt du fluide
}

double Interpolation_IBM_thermal_wall_law::Kader(double yplus, double Prandlt)
{
// Calcul de Beta
  double a = 3.85*pow(Prandlt,1/3) -1.3;
  double b = 2.12 * log(Prandlt);
  double Beta = pow(a,2) + b;
// Calcul de Gamma
  a = 0.01*pow(yplus*Prandlt,4);
  b = 1 + 5*yplus + pow(Prandlt,3);
  double Gamma = a/b;
// Calcul de theta+
  a = yplus * exp(- Gamma) * Prandlt;
  b = 2.12 * log(1 + yplus) + Beta;
  double c = exp (- 1 / Gamma);
  double thetaplus = a + b*c;

  return thetaplus;
}
