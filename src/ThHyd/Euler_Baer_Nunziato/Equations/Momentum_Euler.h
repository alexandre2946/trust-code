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

#include <Operateur_NConserv.h>
#include <Navier_Stokes_std.h>

class Momentum_Euler : public Navier_Stokes_std
{
  Declare_instanciable(Momentum_Euler);
public :

  // overload to do nothing
  void discretiser_grad_p() override { }
  void discretiser_vitesse() override { }

  // overload to skip mother class
  void mettre_a_jour(double temps) override { Equation_base::mettre_a_jour(temps); }
  int has_interface_blocs() const override { return Equation_base::has_interface_blocs(); }
  int impr(Sortie& os) const override { return Equation_base::impr(os); }
  int preparer_calcul() override;

  // overload to change behavior
  Entree& lire_cond_init(Entree&) override; //pour lire la pression
  void discretiser() override;
  void set_param(Param& param) const override;
  void completer() override;
  void abortTimeStep() override;
  void creer_champ(const Motcle& motlu) override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  void mettre_a_jour_champs_conserves(double temps, int reset) override;
  void verifie_ch_init_nb_comp(const Champ_Inc_base& ch_ref, const int nb_comp) const override;

  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int sauvegarder(Sortie& os) const override;
  double calculer_pas_de_temps() const override ;

  bool initTimeStep(double dt) override;
  DoubleTab& corriger_derivee_expl(DoubleTab& derivee) override { return derivee; }
  DoubleTab& corriger_derivee_impl(DoubleTab& derivee) override { return derivee; }
  const Champ_Don_base& diffusivite_pour_transport() const override;
  const Champ_base& diffusivite_pour_pas_de_temps() const override;

  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;

  // new methods
  void mettre_a_jour_p_c();
  void calculer_vitesse_normale();
  void calculer_vitesse();
  void init_alpha_rho_u();

  double alpha_res() const ;
  const Champ_Inc_base& vitesse_phase(const int i) const;
  Champ_Inc_base& vitesse_phase(const int i) ;

  // inline methods
  inline int nombre_d_operateurs() const override { return 2; }
  inline const DoubleTab& vitesse_son() const { return vitesse_son_; }
  inline DoubleTab& vitesse_son() { return vitesse_son_; }
  inline const DoubleTab& vitesse_normale() const { return vitesse_normale_; } ;
  inline DoubleTab& vitesse_normale() { return vitesse_normale_; } ;
  inline const Champ_Inc_base& inconnue() const override { return l_inco_ch_.valeur(); }
  inline Champ_Inc_base& inconnue() override { return l_inco_ch_.valeur(); }
  inline const Champ_Inc_base& vitesse() const override { return la_vitesse.valeur(); }
  inline Champ_Inc_base& vitesse() override { return la_vitesse.valeur(); }

protected:
  OWN_PTR(Champ_Inc_base) l_inco_ch_;
  DoubleTab vitesse_son_, vitesse_normale_;

  std::vector<OWN_PTR(Champ_Inc_base)> vit_phases_; //vitesses de chaque phase
  Motcles noms_vit_phases_; //leurs noms

  Operateur_NConserv terme_nconserv_;
};

#endif /* Momentum_Euler_included */
