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

#ifndef QDM_Multiphase_included
#define QDM_Multiphase_included

#include <Operateur_Evanescence.h>
#include <Navier_Stokes_std.h>
#include <vector>

/*! @brief Carries the terms of the momentum equation for multiphase flow without turbulence modelling.
 *
 *     The quasi-compressible fluid assumption is used.
 *     Under these assumptions, the following form of the Navier-Stokes equations is used:
 *        DU/dt = div(viscous term) - gradP/rho + sources/rho
 *        div U = W
 *     where DU/dt : material derivative of velocity
 *           rho   : density
 *     Note: the implementation allows individual terms (viscous, convective, source) to be neglected.
 *     The unknown is the velocity field.
 *
 * @sa Equation_base Pb_Thermohydraulique_QC Navier_Stokes_std
 */
class QDM_Multiphase : public Navier_Stokes_std
{
  Declare_instanciable(QDM_Multiphase);

public :

  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void completer() override;
  void discretiser_vitesse() override;
  void discretiser_grad_p() override;
  void mettre_a_jour(double temps) override;
  int impr(Sortie& os) const override
  {
    return Equation_base::impr(os); //idem
  }
  bool initTimeStep(double dt) override;
  void abortTimeStep() override;

  void dimensionner_matrice_sans_mem(Matrice_Morse& matrice) override;

  /*
    interface {dimensionner,assembler}_blocs
    specifics: evanescence is taken into account (last)
  */
  int has_interface_blocs() const override;
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override;
  const Champ_Don_base& diffusivite_pour_transport() const override;
  const Champ_base& diffusivite_pour_pas_de_temps() const override;
  const Champ_base& vitesse_pour_transport() const override;
  void assembler_blocs_avec_inertie(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) override;
  void creer_champ(const Motcle& motlu) override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;

  double alpha_res() const ;

protected:
  Entree& lire_cond_init(Entree&) override; // to read the pressure
  int preparer_calcul() override; // calls Equation_base method

  std::vector<OWN_PTR(Champ_Inc_base)> vit_phases_; // velocity fields for each phase
  Motcles noms_vit_phases_; // their names

  std::vector<OWN_PTR(Champ_Fonc_base)> grad_vit_phases_; // velocity gradients for each phase
  Motcles noms_grad_vit_phases_; // their names

  Operateur_Evanescence evanescence_;
};

#endif /* QDM_Multiphase_included */
