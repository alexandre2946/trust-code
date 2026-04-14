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

#ifndef Energy_Euler_included
#define Energy_Euler_included

#include <Conservation_Euler.h>
#include <Milieu_composite_Euler.h>
class Energy_Euler : public Conservation_Euler
{
  Declare_instanciable(Energy_Euler);
public :
  void discretiser() override;
  int verif_Cl() const override { return 1; } // TODO
  void set_param(Param& param) const override;
  int nombre_d_operateurs() const override { return 2; }
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  Entree& lire_cond_init(Entree& is) override;
  double flux_bord(const double& rhoE_bord, const double& vit_n_bord, const double& p_bord) const override;
  DoubleTab flux(const int& f, const int& left_or_right) const override;
  inline double termes_NonConservatif(const double& alpha_bord, const double& vitesse_normale_interieur, const double& p_inter) const override
  {
    return -alpha_bord * vitesse_normale_interieur * p_inter;
  }

  inline void init_energie_tot()
  {
    const Milieu_composite_Euler& mil = ref_cast(Milieu_composite_Euler, milieu());
    mil.init_energie_tot(inconnue().valeurs());
  }

};
#endif /* Energy_Euler_included */
