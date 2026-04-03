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

#ifndef Pb_rayo_semi_transp_included
#define Pb_rayo_semi_transp_included

#include <Equation_rayonnement_base.h>
#include <Probleme_base.h>

/*! @brief Le Pb_rayo_semi_transp est un Probleme_base qui a 4 particularites : * Son equation doit etre typee en fonction de la dicretisation.
 *
 *     Cela impose de differer certaines initialisations jusqu'a
 *     connaitre la discretisation utilisee.
 *   * Il partage son domaine avec un probleme de type hydraulique
 *   * Il n'y a qu'une seule valeur temporelle (futur=present).
 *     Il faudrait en faire un probleme independant du temps.
 *   * Il conserve une ref sur le probleme hydraulique. Cette ref est utilisee de
 *     maniere intensive.
 *
 *
 * @sa Pb_Couple_rayo_semi_transp Equation_rayonnement_base
 */
class Pb_rayo_semi_transp: public Probleme_base
{
  Declare_instanciable(Pb_rayo_semi_transp);
public:
  void terminate() override { finir(); }

  double computeTimeStep(bool& stop) const override
  {
    stop=false;
    return DMAXFLOAT;
  }

  bool initTimeStep(double dt) override;
  bool iterateTimeStep(bool& converged) override;
  void validateTimeStep() override;

  void completer() override { }
  int nombre_d_equations() const override { return 1; }

  double calculer_pas_de_temps() const override  {  return DMAXFLOAT;  }

  // Cette methode ne doivent pas servir : on passe par l'interface de Problem
  void mettre_a_jour(double temps) override { Process::exit(); }

  void preparer_calcul() override;
  void discretiser(Discretisation_base&) override;
  void associer_sch_tps_base(const Schema_Temps_base&) override;

  Champ_Inc_base& put_irradience();
  const Champ_front_base& flux_radiatif(const Nom& nom_bord) const;
  void calculer_flux_radiatif();

  void discretise_longueur_rayo();
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  void creer_champ(const Motcle& motlu) override { }

  inline Probleme_base& probleme_fluide() { return pb_fluide_.valeur(); }
  inline const Probleme_base& probleme_fluide() const { return pb_fluide_.valeur(); }
  inline const double& valeur_sigma() const { return sigma_; }
  inline void associer_probleme_fluide(Probleme_base& Pb) { pb_fluide_ = Pb; }
  void typer_lire_milieu(Entree& is) override { /* Do nothing */ }

  inline const Equation_base& equation(int i) const override
  {
    assert(i==0);
    return eq_rayo_;
  }
  inline Equation_base& equation(int i) override
  {
    assert(i==0);
    return eq_rayo_;
  }
  inline const Equation_base& get_equation_by_name(const Nom& un_nom) const override
  {
    assert(Motcle(un_nom)==Motcle("Eq_rayo_semi_transp"));
    return eq_rayo_;
  }

  inline Equation_base& getset_equation_by_name(const Nom& un_nom) override
  {
    assert(Motcle(un_nom)==Motcle("Eq_rayo_semi_transp"));
    return eq_rayo_;
  }

  inline Equation_rayonnement_base& eq_rayo()
  {
    assert(eq_rayo_.non_nul());
    return eq_rayo_.valeur();
  }

  inline const Equation_rayonnement_base& eq_rayo() const
  {
    assert(eq_rayo_.non_nul());
    return eq_rayo_.valeur();
  }

protected :
  OBS_PTR(Probleme_base) pb_fluide_;
  OWN_PTR(Equation_rayonnement_base) eq_rayo_;
  static constexpr double sigma_ = 5.67e-8;
};

#endif /* Pb_rayo_semi_transp_included */
