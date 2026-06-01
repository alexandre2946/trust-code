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

#ifndef Convection_Diffusion_Chaleur_QC_included
#define Convection_Diffusion_Chaleur_QC_included

#include <Convection_Diffusion_Chaleur_Fluide_Dilatable_base.h>

/*! @brief Particular case of Convection_Diffusion_Chaleur_Fluide_Dilatable_base for a quasi-compressible fluid, where the transported scalar is temperature (ideal gases) or enthalpy (real gases).
 *
 * Generalisation of Convection_Diffusion_Temperature for real gases.
 *
 * @sa Conv_Diffusion_std Convection_Diffusion_Temperature
 */

class Convection_Diffusion_Chaleur_QC : public Convection_Diffusion_Chaleur_Fluide_Dilatable_base
{
  Declare_instanciable_sans_constructeur(Convection_Diffusion_Chaleur_QC);

public :
  Convection_Diffusion_Chaleur_QC();
  void set_param(Param& titi) const override;
  void calculer_div_u_ou_div_rhou(DoubleTab& res) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void mettre_a_jour(double) override;
  int preparer_calcul() override;
  const Champ_base& vitesse_pour_transport() const override;

  // Inline methods
  inline bool is_generic() const override { return mode_convection_ == 2 ? true : false;}

protected :
  double TMIN_ = std::numeric_limits<double>::quiet_NaN();
  double TMAX_ = std::numeric_limits<double>::quiet_NaN();
  int mode_convection_; // 0: by divergence of u, 1: by conv(u), 2: by conv(rho*u)
};

#endif /* Convection_Diffusion_Chaleur_QC_included */
