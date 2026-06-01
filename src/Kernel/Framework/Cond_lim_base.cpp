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

#include <Discretisation_base.h>
#include <Domaine_Cl_dis_base.h>
#include <Cond_lim_base.h>
#include <Equation_base.h>

Implemente_base(Cond_lim_base, "Cond_lim_base", Objet_U);
// XD condlim_base objet_u condlim_base NO_BRACE Basic class of boundary conditions.

Sortie& Cond_lim_base::printOn(Sortie& s) const { return s << le_champ_front; }

Entree& Cond_lim_base::readOn(Entree& s) { return s >> le_champ_front; }

/*! @brief DOES NOTHING must be overridden in derived classes
 *
 */
void Cond_lim_base::completer()
{
  champ_front().completer();
}

int Cond_lim_base::compatible_avec_eqn(const Equation_base& eqn) const
{
  if (app_domains.size() == 0)
    {
      Cerr << "You call Cond_lim_base::compatible_avec_eqn but the std::vector app_domains is not filled ! Check your readOn !!" << finl;
      Process::exit();
    }

  Motcle dom_app = eqn.domaine_application();
  for (const auto &itr : app_domains)
    if (itr == dom_app) return 1;

  err_pas_compatible(eqn);
  return 0;
}

/*! @brief Changes the i-th future time of the BC.
 *
 */
void Cond_lim_base::changer_temps_futur(double temps, int i)
{
  champ_front().changer_temps_futur(temps, i);
}

/*! @brief Rotates the wheel of the BC
 *
 */
int Cond_lim_base::avancer(double temps)
{
  return champ_front().avancer(temps);
}

/*! @brief Rotates the wheel of the BC
 *
 */
int Cond_lim_base::reculer(double temps)
{
  return champ_front().reculer(temps);
}

/*! @brief Initialization at the beginning of the calculation.
 *
 * Must be called before any calculate_exchange_coefficients or update
 *     Unlike the update methods, the
 *     initialize methods of BCs cannot depend on the outside
 *     (it may not be initialized itself)
 *
 * @return (0 in case of error, 1 otherwise.)
 */
int Cond_lim_base::initialiser(double temps)
{
  return le_champ_front->initialiser(temps, domaine_Cl_dis().inconnue());
}

/*! @brief Performs a time update of the boundary condition.
 *
 * @param (double temps) the time step for update
 */
void Cond_lim_base::mettre_a_jour(double temps)
{
  le_champ_front->mettre_a_jour(temps);
}

/* @brief Reset current time for the boundary condition.
 *
 * A BC is regarded as an input, so here this is equivalent to a 'mettre_a_jour'
 */
void Cond_lim_base::resetTime(double time)
{
  mettre_a_jour(time);
}

/*! @brief Indicates whether this boundary condition must be updated during sub-time steps of a time scheme such as RK.
 *
 *   By default it returns 0 to indicate that no update
 *   is necessary; it must be overloaded to return 1 if needed
 *   (example: Echange_impose_base)
 *
 * @param (double temps) the time step for update
 */
int Cond_lim_base::a_mettre_a_jour_ss_pas_dt()
{
  return 0;
}

/*! @brief Computation of exchange coefficients for coupling via Champ_front_contact_VEF.
 *
 *     These computations are local to the problem and depend only on
 *     the unknown. They must therefore be done every time the
 *     unknown is modified. They are available externally and
 *     stored in the BCs.
 *     WEC: Champ_front_contact_VEF should disappear and this
 *     method with it!!!
 *
 * @param (double temps) the time step for update
 */
void Cond_lim_base::calculer_coeffs_echange(double temps)
{
  le_champ_front->calculer_coeffs_echange(temps);
}

/*! @brief Calls the verification of the field read through the equation for which the boundary condition is considered.
 *
 *  The method is overloaded in cases where the user must
 *  specify the boundary field.
 *
 */
void Cond_lim_base::verifie_ch_init_nb_comp() const
{

}

/*! @brief Associates the boundary with the object.
 *
 * The Frontiere_dis_base object is in fact associated with the member
 *     OWN_PTR(Champ_front_base) of the Cond_lim_base object that represents the field of boundary
 *     conditions imposed on the boundary.
 *
 * @param (Frontiere_dis_base& fr) the boundary on which the boundary conditions are imposed
 */
void Cond_lim_base::associer_fr_dis_base(const Frontiere_dis_base& fr)
{
  assert(le_champ_front);
  le_champ_front->associer_fr_dis_base(fr);
  modifier_val_imp = 0;
}

/*! @brief Associates the Domaine_Cl_dis_base (domain of discretized boundary conditions) with the object.
 *
 * This Domaine_Cl_dis_base stores (references) all the boundary
 *     conditions relative to a geometric domain.
 *
 * @param (Domaine_Cl_dis_base& zcl) a domain of discretized boundary conditions to which the Cond_lim_base object refers
 */
void Cond_lim_base::associer_domaine_cl_dis_base(const Domaine_Cl_dis_base& zcl)
{
  mon_dom_cl_dis = zcl;
  le_champ_front->verifier(*this);
}

/*! @brief Returns 1 if the boundary condition is compatible with the discretization passed as parameter.
 *
 * @param (Discretisation_base&) the discretization with which we want to verify compatibility
 */
int Cond_lim_base::compatible_avec_discr(const Discretisation_base& discr) const
{
  if (supp_discs.size() == 0) return 1;
  else
    {
      Nom type_discr = discr.que_suis_je();
      for (const auto &itr : supp_discs)
        if (itr == type_discr) return 1;

      err_pas_compatible(discr);
      return 0;
    }
}

/*! @brief This method is called when the boundary condition is not compatible with the equation on which we try
 *
 *     to apply it.
 *
 * @param (Equation_base& eqn) the equation with which the boundary condition is incompatible
 */
void Cond_lim_base::err_pas_compatible(const Equation_base& eqn) const
{
  Cerr << "The boundary condition " << que_suis_je() << " can't apply to " << finl << "the equation of kind " << eqn.que_suis_je() << finl;
  exit();
}

/*! @brief This method is called when the boundary condition is not compatible with the discretization on which we try
 *
 *     to apply it.
 *
 * @param (Discretisation_base& discr) the discretization with which the boundary condition is incompatible
 */
void Cond_lim_base::err_pas_compatible(const Discretisation_base& discr) const
{
  Cerr << "The boundary condition " << que_suis_je() << " can't be used with " << finl << "the discretization of kind " << discr.que_suis_je() << finl;
  exit();
}

void Cond_lim_base::champ_front(int face, DoubleVect& var) const
{
  le_champ_front->valeurs_face(face, var);
}

void Cond_lim_base::injecter_dans_champ_inc(const Champ_Inc_base&) const
{
  Cerr << "Cond_lim_base::injecter_dans_champ_inc()" << finl;
  Cerr << "this method does nothing and must be overloaded " << finl;
  Cerr << "Contact TRUST support." << finl;
  exit();
}
/*! @brief Changes the i-th future time of the BC.
 *
 */
void Cond_lim_base::set_temps_defaut(double temps)
{
  champ_front().set_temps_defaut(temps);
}
/*! @brief Called by Conds_lim::completer. Calls cha_front_base::fixer_nb_valeurs_temporelles.
 *
 */
void Cond_lim_base::fixer_nb_valeurs_temporelles(int nb_cases)
{
  champ_front().fixer_nb_valeurs_temporelles(nb_cases);
}
