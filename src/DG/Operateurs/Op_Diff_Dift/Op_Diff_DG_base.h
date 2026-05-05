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

#ifndef Op_Diff_DG_base_included
#define Op_Diff_DG_base_included

#include <Op_Diff_Turbulent_base.h>
#include <Operateur_Diff_base.h>
#include <Domaine_DG.h>
#include <TRUST_Ref.h>
#include <SFichier.h>
#include <Champ_Uniforme.h>

class Domaine_Cl_DG;
/**
 * @brief This class provides the common infrastructure shared by all DG diffusion operators in TRUST.
 * It handles the association with the DG domain and its boundary conditions, the management of
 * an effective diffusivity field (nu_), and the computation of a stable explicit time step.
 *
 * The stable time step calculation (calculer_dt_stab) follows a diffusive CFL criterion:
 *   dt ~ h^2 / (2 * dim * alpha_max)
 * with special handling for:
 *  - Variable-density flows: rho * h^2 / (2 * dim * nu), with VDF-like and general mesh branches.
 *  - Robin / external-heat-exchange boundary conditions (Echange_externe_impose): the effective
 *    diffusivity is scaled by the Biot number when Bi > 1, making the criterion more conservative.
 *
 * Derived classes are responsible for implementing the actual flux assembly (ajouter()) according
 * to their specific DG formulation.
 */
class Op_Diff_DG_base: public Operateur_Diff_base, public Op_Diff_Turbulent_base
{
  Declare_base(Op_Diff_DG_base);
public:
  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;

  double calculer_dt_stab() const override;

  void associer_diffusivite(const Champ_base& diffu) override
  {
    diffusivite_ = diffu;
    is_var_ = sub_type(Champ_Uniforme, diffu) ? 0 : 1;
    is_aniso_ = (diffu.nb_comp() > 1);
  }

  void completer() override;
  const Champ_base& diffusivite() const override { return diffusivite_.valeur(); }
  void mettre_a_jour(double t) override
  {
    Operateur_base::mettre_a_jour(t);
    nu_a_jour_ = 0;
  }

  void update_nu() const; //met a jour nu
  inline double nu(int i, int compo) const { return nu_(is_var_ * i, compo); }

  DoubleTab& calculer(const DoubleTab&, DoubleTab&) const override;
  int impr(Sortie& os) const override;

protected:
  OBS_PTR(Domaine_DG) le_dom_dg_;
  OBS_PTR(Domaine_Cl_DG) la_zcl_dg_;
  mutable SFichier Flux, Flux_moment, Flux_sum; // Fichiers .out

  OBS_PTR(Champ_base) diffusivite_;
  mutable int nu_a_jour_ = 0; //si on doit mettre a jour nu
  mutable DoubleTab nu_;

  bool is_var_; //if the diffusivity is Uniforme or heterogeneous
  bool is_aniso_; //if the diffusivity is anisotropic

};


#endif /* Op_Diff_DG_base_included */
