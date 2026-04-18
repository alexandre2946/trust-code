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

#ifndef Momentum_Euler_included
#define Momentum_Euler_included

#include <Operateur_Evanescence.h>
#include <Operateur_NConserv.h>
#include <Navier_Stokes_std.h>
#include <vector>

/*! @brief classe Momentum_Euler Cette classe porte les termes de l'equation de la dynamique
 *
 *     pour un fluide sans modelisation de la turbulence.
 *     On suppose l'hypothese de fluide quasi compressible.
 *     Sous ces hypotheses, on utilise la forme suivante des equations de
 *     Navier_Stokes:
 *        DU/dt = div(terme visqueux) - gradP/rho + sources/rho
 *        div U = W
 *     avec DU/dt : derivee particulaire de la vitesse
 *          rho   : masse volumique
 *     Rq : l'implementation de la classe permet bien sur de negliger
 *          certains termes de l'equation (le terme visqueux, le terme
 *          convectif, tel ou tel terme source).
 *     L'inconnue est le champ de vitesse.
 *
 * @sa Equation_base Pb_Thermohydraulique_QC Navier_Stokes_std
 */
class Momentum_Euler : public Navier_Stokes_std
{
  Declare_instanciable(Momentum_Euler);
public :
  void discretiser() override;
  int nombre_d_operateurs() const override {return 3; }
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void completer() override;
  void discretiser_vitesse() override;
  int sauvegarder(Sortie& os) const override;
  void discretiser_grad_p() override;
  void mettre_a_jour(double temps) override;
  void  mettre_a_jour_p_c();
  DoubleTab& corriger_derivee_expl(DoubleTab& derivee) override;
  DoubleTab& corriger_derivee_impl(DoubleTab& derivee) override;
  int impr(Sortie& os) const override
  {
    return Equation_base::impr(os); //idem
  }
  bool initTimeStep(double dt) override;
  void abortTimeStep() override;

  inline const Champ_Inc_base& vitesse() const override { return la_vitesse.valeur(); }
  inline Champ_Inc_base& vitesse() override { return la_vitesse.valeur(); }

  virtual const Champ_Inc_base& vitesse_phase(const int i) const;
  virtual Champ_Inc_base& vitesse_phase(const int i) ;

  int has_interface_blocs() const override;
  //void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override;
  const Champ_Don_base& diffusivite_pour_transport() const override;
  const Champ_base& diffusivite_pour_pas_de_temps() const override;
  const Champ_base& vitesse_pour_transport() const override;
  //void assembler_blocs_avec_inertie(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) override { };
  void creer_champ(const Motcle& motlu) override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  double calculer_pas_de_temps() const override ;
  double alpha_res() const ;
  const DoubleTab& vitesse_son() const {return vitesse_son_;}
  DoubleTab& vitesse_son() {return vitesse_son_;}
  void calculer_vitesse_normale() ;
  void calculer_vitesse();
  const DoubleTab& vitesse_normale() const {return vitesse_normale_;};
  DoubleTab& vitesse_normale() {return vitesse_normale_;};
  const Champ_Inc_base& inconnue() const override { return l_inco_ch_.valeur();}
  Champ_Inc_base& inconnue() override { return l_inco_ch_.valeur();}
  void init_alpha_rho_u();
  virtual DoubleTab flux_(const int f, const int left_or_right) const;
  void mettre_a_jour_champs_conserves(double temps, int reset) override;

protected:
  OWN_PTR(Champ_Inc_base) l_inco_ch_;
  DoubleTab vitesse_son_, vitesse_normale_;
  Entree& lire_cond_init(Entree&) override; //pour lire la pression
  int preparer_calcul() override; //appelle la methode de Equation_base

  std::vector<OWN_PTR(Champ_Inc_base)> vit_phases_; //vitesses de chaque phase
  Motcles noms_vit_phases_; //leurs noms

  Operateur_NConserv terme_nconserv_;
};

#endif /* Momentum_Euler_included */
