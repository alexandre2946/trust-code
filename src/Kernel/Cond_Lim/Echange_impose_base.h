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

#ifndef Echange_impose_base_included
#define Echange_impose_base_included

#include <Cond_lim_base.h>
#include <TRUSTTab.h>

/*! @brief Echange_impose_base: This boundary condition is used only for the energy equation.
 *
 *     It corresponds to imposing a heat exchange with the exterior
 *     of the domain by imposing an external temperature T_ext and an
 *     exchange coefficient h_imp.
 *     The flux term computed from the pair (h_imp, T_ext) is written:
 *                            h_t(T_ext - T_interior)*Surf
 *                           where h_t : global exchange coefficient.
 *     It appears on the right-hand side of the energy equation.
 *     Either the user provides an exchange coefficient corresponding
 *     only to the wall, in which case the program will compute the diffusion
 *     over the half-cell near the wall, or the user provides a global exchange
 *     coefficient that directly accounts for both the above.
 *     The two derived classes Echange_externe_impose and Echange_global_impose
 *     represent these two possibilities.
 *
 * @sa Cond_lim_base Echange_externe_impose Echange_global_impose
 */
class Echange_impose_base : public Cond_lim_base
{
  Declare_base_sans_constructeur(Echange_impose_base);
public:

  inline bool has_emissivite() const { return bool(emissivite_); }
  inline bool has_h_imp() const { return bool(h_imp_); }

  const DoubleTab& tab_T_ext(double temps=DMAXFLOAT) const;
  const DoubleTab& tab_h_imp(double temps=DMAXFLOAT) const;
  const DoubleTab& tab_emissivite(double temps=DMAXFLOAT) const;
  virtual double T_ext(int num) const;
  virtual double T_ext(int num,int k) const;
  virtual double h_imp(int num) const;
  virtual double h_imp(int num,int k) const;
  double emissivite(int num) const;
  double emissivite(int num,int k) const;

  /*! @brief Returns the T_ext field of temperature imposed at the boundary.
   *
   * @return (Champ_front_base&) the T_ext field of temperature imposed at the boundary
   */
  inline virtual Champ_front_base& T_ext() { return le_champ_front.valeur(); }
  inline virtual const Champ_front_base& T_ext() const { return le_champ_front.valeur(); }

  inline virtual Champ_front_base& h_imp() { assert (has_h_imp()); return h_imp_.valeur(); }
  inline virtual const Champ_front_base& h_imp() const { assert (has_h_imp()); return h_imp_.valeur(); }

  inline Champ_front_base& emissivite() {  assert (has_emissivite()); return emissivite_.valeur(); }
  inline const Champ_front_base& emissivite() const  {  assert (has_emissivite()); return emissivite_.valeur(); }

  // Used in the CAL of flux computation for wall laws
  virtual void liste_faces_loi_paroi(IntTab&) { }

  void mettre_a_jour(double ) override;
  int initialiser(double temps) override;
  int a_mettre_a_jour_ss_pas_dt() override { return 1; }

  // added method to avoid operating directly on champ_front
  void set_temps_defaut(double temps) override;
  void fixer_nb_valeurs_temporelles(int nb_cases) override;
  //
  void changer_temps_futur(double temps,int i) override;
  int avancer(double temps) override;
  int reculer(double temps) override;
  void associer_fr_dis_base(const Frontiere_dis_base& ) override ;

  virtual bool has_h_imp_grad() const { return false; }
  virtual double h_imp_grad(int num) const { Process::exit(que_suis_je()+ " : h_imp_grad must be overloaded !" ) ; return -1.e10 ;};
  virtual double h_imp_grad(int num,int k) const  { Process::exit(que_suis_je()+ " : h_imp_grad must be overloaded !") ; return -1.e10 ;};

protected :
  OWN_PTR(Champ_front_base) h_imp_, emissivite_ /* si Echange_externe_radiatif */;
private:
  // Stores all values on the faces (useful for GPU):
  mutable DoubleTab text_;
  mutable DoubleTab himp_;
  mutable DoubleTab eps_;
};

#endif /* Echange_impose_base_included */
