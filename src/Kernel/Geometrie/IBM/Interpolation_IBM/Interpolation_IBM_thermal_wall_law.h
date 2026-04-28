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

#ifndef Interpolation_IBM_thermal_wall_law_included
#define Interpolation_IBM_thermal_wall_law_included

#include <Interpolation_IBM_elem_fluid.h>

class Interpolation_IBM_thermal_wall_law : public Interpolation_IBM_elem_fluid
{
  Declare_instanciable( Interpolation_IBM_thermal_wall_law ) ;
public :
  double Kader(double, double);

  inline int get_formulation_Tp() { return formulation_Tp_; }
  inline int get_boundary_type() { return boundary_type_; }
  inline double get_T_inlet() { return T_inlet_; }
  inline double get_Prandlt_mol() { return Prandlt_mol_; }

  void set_param(Param&) const override;

protected :
  int formulation_Tp_ = 0; // Choix loi de theta+
  int boundary_type_ = 0; // Choix du type de condition limite (0 = temp imp., 1 = flux imp.)
  double T_inlet_ = 0.0;
  double Prandlt_mol_ = 0.0;
};

#endif /* Interpolation_IBM_thermal_wall_law_included */
