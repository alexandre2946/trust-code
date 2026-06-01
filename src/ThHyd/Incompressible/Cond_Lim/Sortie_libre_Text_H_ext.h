/****************************************************************************
* Copyright (c) 2024, CEA
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

#ifndef Sortie_libre_Text_H_ext_included
#define Sortie_libre_Text_H_ext_included

#include <TRUST_Ref.h>
#include <Neumann.h>

/*! @brief Neumann_sortie_libre This class represents an open boundary without imposed velocity.
 *
 *     For Navier-Stokes equations, pressure must be imposed on such a boundary.
 *     To handle the hydraulics, the class Sortie_libre_pression_imposee is derived from Neumann_sortie_libre.
 *     Boundary conditions of type Neumann_sortie_libre or derived types result in zero diffusive fluxes.
 *     However, the treatment of convective fluxes requires knowing the convected field outside the boundary in case of fluid re-entry.
 *     This is why the class carries an OWN_PTR(Champ_front_base) (member le_champ_ext).
 *     In the computation operators, boundary conditions of type Neumann_sortie_libre and derived types will be treated identically.
 *
 * @sa Neumann
 */
class Sortie_libre_Text_H_ext: public Neumann
{
  Declare_instanciable(Sortie_libre_Text_H_ext);
public:
  virtual double val_ext(int i) const;
  virtual double val_ext(int i, int j) const;
  inline void bascule_cond_lim_en_enthalpie() { type_cond_lim = 1; }
  inline void bascule_cond_lim_en_temperature() { type_cond_lim = 0; }
  double flux_impose(int i) const override;
  double flux_impose(int i, int j) const override;

protected:
  OWN_PTR(Champ_front_base) le_champ_Text, le_champ_hext;
  int type_cond_lim = -1;
};

#endif
