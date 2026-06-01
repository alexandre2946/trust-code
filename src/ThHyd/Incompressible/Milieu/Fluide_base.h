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

#ifndef Fluide_base_included
#define Fluide_base_included

#include <Milieu_base.h>
#include <TRUST_Ref.h>

class Champ_base;

/*! @brief Base class for an incompressible fluid and its properties:
 *         - kinematic viscosity (mu)
 *         - dynamic viscosity   (nu)
 *         - density             (rho)
 *         - diffusivity         (alpha)
 *         - conductivity        (lambda)
 *         - specific heat       (Cp)
 *         - constituent thermal expansion coefficient (beta_co)
 *
 * @sa Milieu_base
 */
class Fluide_base : public Milieu_base
{
  Declare_instanciable(Fluide_base);
public :
  void set_param(Param& param) const override;
  void verifier_coherence_champs(int& err, Nom& message) override;
  bool initTimeStep(double dt) override;
  void mettre_a_jour(double) override;
  int initialiser(const double temps) override;
  void creer_champs_non_lus() override;
  void discretiser(const Probleme_base& pb, const Discretisation_base& dis) override;
  virtual void set_h0_T0(double h0, double T0);
  virtual int is_incompressible() const { return 0; }

  const Champ_base& energie_interne() const;
  Champ_base& energie_interne();
  const Champ_base& enthalpie() const;
  Champ_base& enthalpie();
  const Champ_base& temperature_multiphase() const;
  Champ_base& temperature_multiphase();

  inline const Champ_Don_base& viscosite_cinematique() const { return ch_nu_.valeur(); }
  inline Champ_Don_base& viscosite_cinematique() { return ch_nu_.valeur(); }
  inline const Champ_Don_base& viscosite_dynamique() const { return ch_mu_.valeur(); }
  inline Champ_Don_base& viscosite_dynamique() { return ch_mu_.valeur(); }
  bool has_viscosite_dynamique() const { return bool(ch_mu_); }

  // Returns the constituent expansion coefficient, beta_co.
  inline const Champ_Don_base& beta_c() const { return ch_beta_co_.valeur(); }
  inline Champ_Don_base& beta_c() { return ch_beta_co_.valeur(); }
  bool has_beta_c() const { return bool(ch_beta_co_); }

  // Returns the fluid absorption coefficient
  inline Champ_Don_base& kappa() { return coeff_absorption_.valeur(); }
  inline const Champ_Don_base& kappa() const { return coeff_absorption_.valeur(); }
  bool has_kappa() const { return bool(coeff_absorption_); }

  // Returns the fluid refractive index
  inline Champ_Don_base& indice() { return indice_refraction_.valeur(); }
  inline const Champ_Don_base& indice() const { return indice_refraction_.valeur(); }

  // Returns the radiation penetration depth in the fluid defined as l = 1/(3*kappa)
  inline Champ_Don_base& longueur_rayo() { return longueur_rayo_.valeur(); }
  inline const Champ_Don_base& longueur_rayo() const { return longueur_rayo_.valeur(); }
  void typer_longeur_rayo(const Nom& typ) { longueur_rayo_.typer(typ); }

  inline bool is_rayo_semi_transp() const override { return (coeff_absorption_ && indice_refraction_); }
  inline bool is_rayo_transp() const override { return is_rad_transp_med_; }
  inline void set_rayo_transp_flag() { is_rad_transp_med_ = true; }
  inline bool is_longueur_rayo_discretised() const { return bool(longueur_rayo_); }

protected :
  void creer_e_int() const; // creation sur demande de e_int / h
  void creer_temperature_multiphase() const; // only if Energie_Multiphase_Enthalpie
  void calculer_temperature_multiphase() const; // only if Energie_Multiphase_Enthalpie

  mutable int e_int_auto_ = 0; //1 if e_int was created automatically
  static void calculer_e_int(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv); // default computation function

  mutable OWN_PTR(Champ_base) ch_e_int_, ch_h_ou_T_; //pour la creation sur demande : h is Energie_Multiphase et T si Energie_Multiphase_Enthalpie
  OWN_PTR(Champ_Don_base) ch_mu_, ch_nu_, ch_beta_co_;
  double h0_ = 0, T0_ = 0;

  bool is_rad_transp_med_ = false; // transparent radiating fluid

  // Parameters of the semi-transparent radiating fluid
  OWN_PTR(Champ_Don_base) coeff_absorption_, indice_refraction_;

  // Characteristic radiation penetration depth in the semi-transparent medium defined as l = 1/(3*kappa)
  OWN_PTR(Champ_Don_base) longueur_rayo_;

  void creer_nu();
  virtual void calculer_nu();
};

#endif /* Fluide_base_included */
