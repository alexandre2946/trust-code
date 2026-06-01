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

#ifndef Cond_lim_base_included
#define Cond_lim_base_included

#include <Champ_front_base.h>
#include <TRUST_Ref.h>
#include <Motcle.h>
#include <vector>

class Cond_lim_rayo_milieu_transp;
class Cond_lim_rayo_semi_transp;
class Discretisation_base;
class Domaine_Cl_dis_base;
class Equation_base;

/*! @brief class Cond_lim_base Base class for the hierarchy of classes that represent the different boundary conditions (Dirichlet, Neumann ...).
 *
 *      A boundary condition object serves to define, for a given equation, the boundary conditions to apply on a boundary of a domain.
 *      Each Cond_lim_base object contains a reference to the Domaine_Cl_dis_base object it is part of.
 *      Each object also contains an OWN_PTR(Champ_front_base) object containing the values to impose on the boundary.
 *
 * @sa Cond_lim Domaine_Cl_dis_base Frontiere_dis_base, Abstract class from which all objects representing boundary conditions must derive.,
 *     Abstract method:, int compatible_avec_eqn(const Equation_base&) const
 */
class Cond_lim_base : public Objet_U
{
  Declare_base(Cond_lim_base);
public:
  virtual void completer();
  virtual int initialiser(double temps);
  virtual void mettre_a_jour(double temps);
  virtual void resetTime(double time);
  virtual void calculer_coeffs_echange(double temps);
  virtual void verifie_ch_init_nb_comp() const;
  virtual inline Frontiere_dis_base& frontiere_dis();
  virtual inline const Frontiere_dis_base& frontiere_dis() const;
  virtual void associer_fr_dis_base(const Frontiere_dis_base&);
  inline Domaine_Cl_dis_base& domaine_Cl_dis();
  inline const Domaine_Cl_dis_base& domaine_Cl_dis() const;
  virtual void associer_domaine_cl_dis_base(const Domaine_Cl_dis_base&);
  inline Champ_front_base& champ_front();
  inline const Champ_front_base& champ_front() const;

  virtual void set_temps_defaut(double temps);
  virtual void fixer_nb_valeurs_temporelles(int nb_cases);

  virtual void champ_front(int, DoubleVect&) const;
  virtual int compatible_avec_eqn(const Equation_base&) const ;
  virtual int compatible_avec_discr(const Discretisation_base&) const;
  virtual void injecter_dans_champ_inc(const Champ_Inc_base&) const;

  virtual int a_mettre_a_jour_ss_pas_dt();

  // method to set the modifier_val_impl flag
  inline void set_modifier_val_imp(int);
  virtual void changer_temps_futur(double temps, int i);
  virtual int avancer(double temps);
  virtual int reculer(double temps);

  // virtual methods for radiation BCs! Note: Cond_lim_rayo_milieu_transp and Cond_lim_rayo_semi_transp do not derive from Objet_U
  virtual bool is_bc_rayo_milieu_transp(Cond_lim_rayo_milieu_transp*& la_cl_rayo)
  {
    return false; /* not radiation by default! */
  }

  virtual bool is_bc_rayo_semi_transp(Cond_lim_rayo_semi_transp*& la_cl_rayo)
  {
    return false; /* not radiation by default! */
  }

protected:
  std::vector<Motcle> app_domains;
  std::vector<Nom> supp_discs;
  OWN_PTR(Champ_front_base) le_champ_front;
  OBS_PTR(Domaine_Cl_dis_base) mon_dom_cl_dis;
  void err_pas_compatible(const Equation_base&) const;
  void err_pas_compatible(const Discretisation_base&) const;

  // flag to indicate whether or not to modify the imposed value on the boundary condition
  int modifier_val_imp = 0;
};

/*! @brief Returns the discretized boundary to which the boundary conditions apply.
 *
 * @return (Frontiere_dis_base&) the discretized boundary to which the boundary conditions are associated
 */
inline Frontiere_dis_base& Cond_lim_base::frontiere_dis()
{
  return le_champ_front->frontiere_dis();
}

/*! @brief Returns the discretized boundary to which the boundary conditions apply.
 *
 *     (const version)
 *
 * @return (Frontiere_dis_base&) the discretized boundary to which the boundary conditions are associated
 */
inline const Frontiere_dis_base& Cond_lim_base::frontiere_dis() const
{
  return le_champ_front->frontiere_dis();
}

/*! @brief Returns the domain of discretized boundary conditions to which the object belongs.
 *
 * @return (Domaine_Cl_dis_base&) the domain of discretized boundary conditions to which the object belongs
 */
inline Domaine_Cl_dis_base& Cond_lim_base::domaine_Cl_dis()
{
  return mon_dom_cl_dis.valeur();
}

/*! @brief Returns the domain of discretized boundary conditions to which the object belongs.
 *
 *     (const version)
 *
 * @return (Domaine_Cl_dis_base&) the domain of discretized boundary conditions to which the object belongs
 */
inline const Domaine_Cl_dis_base& Cond_lim_base::domaine_Cl_dis() const
{
  return mon_dom_cl_dis.valeur();
}

inline Champ_front_base& Cond_lim_base::champ_front()
{
  return le_champ_front;
}

inline const Champ_front_base& Cond_lim_base::champ_front() const
{
  return le_champ_front;
}

/*! @brief Sets the modifier_val_imp flag to the given value: - if drap == 1: modifier_val_imp=1
 *
 *     - otherwise    : modifier_val_imp=0
 *     This flag allows the BC to know whether it should return the stored value
 *     as-is, or whether it should translate it for the calling object.
 *     It is up to the BC to then decide what translation to perform.
 *     See application in Temperature_imposee_paroi_H, where the BC returns the enthalpy
 *     by default (modifier_val_imp=1), or the temperature otherwise (case of the diffusion operator).
 *     The flag is set in equation.derivee_en_temps_inco
 *     according to the operator that will be called.
 *
 * @param (drap) value to assign to the flag
 */
inline void Cond_lim_base::set_modifier_val_imp(int drap)
{
  // set the flag to 0 or 1
  modifier_val_imp = (drap==1);
}

#endif /* Cond_lim_base_included */
