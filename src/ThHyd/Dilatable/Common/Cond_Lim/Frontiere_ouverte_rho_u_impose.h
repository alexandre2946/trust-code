/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Frontiere_ouverte_rho_u_impose_included
#define Frontiere_ouverte_rho_u_impose_included

#include <Dirichlet_entree_fluide_leaves.h>
#include <TRUST_Ref.h>

class Fluide_Dilatable_base;

/*! @brief @brief Open boundary on which the mass flux rho.U is imposed instead of the velocity U.
 *
 * The velocity is computed by dividing by rho(n+1) found in the fluid at the time of the val_imp call.
 *
 * @sa Dirichlet_entree_fluide
 */
class Frontiere_ouverte_rho_u_impose  : public Entree_fluide_vitesse_imposee_libre
{
  Declare_instanciable(Frontiere_ouverte_rho_u_impose);
public :
  void completer() override;
  int compatible_avec_eqn(const Equation_base&) const override;
  double val_imp_au_temps(double temps, int i) const override;
  double val_imp_au_temps(double temps, int i, int j) const override;
  const DoubleTab& tab_val_imp(double temps=DMAXFLOAT) const override;

protected :
  OBS_PTR(Fluide_Dilatable_base) le_fluide;
};

#endif /* Frontiere_ouverte_rho_u_impose_included */
