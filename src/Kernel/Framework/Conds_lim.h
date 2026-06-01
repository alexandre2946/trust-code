/****************************************************************************
* Copyright (c) 2024, CEA
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

#ifndef Conds_lim_included
#define Conds_lim_included

#include <TRUST_Vector.h>

#include <Cond_lim.h>

/*! @brief class Conds_lim This class represents a vector of boundary conditions.
 *
 *      An object of this type is carried by each Domaine_Cl_dis_base associated
 *      with an equation. A class representing a vector of objects is
 *      declared through the macro VECT(classe_X)
 *
 * @sa Cond_lim Domaine_Cl_dis_base
 */
class Conds_lim : public VECT(Cond_lim)
{
  Declare_instanciable(Conds_lim);
public:

  inline int initialiser(double temps);
  inline void mettre_a_jour(double temps);
  inline void resetTime(double temps);
  inline void calculer_coeffs_echange(double temps);
  void completer(const Domaine_dis_base&);
  inline int compatible_avec_eqn(const Equation_base&) const;
  inline int compatible_avec_discr(const Discretisation_base&) const;

  inline void set_modifier_val_imp(int);
};


inline int Conds_lim::initialiser(double temps)
{
  int ok = 1;
  for (auto &itr : *this)
    ok = ok && itr->initialiser(temps);
  return ok;
}

/*! @brief Time update of all boundary conditions in the vector.
 *
 * @param (double temps) the time step for update
 */
inline void Conds_lim::mettre_a_jour(double temps)
{
  for (auto& itr : *this) itr->mettre_a_jour(temps);
}

inline void Conds_lim::resetTime(double temps)
{
  for (auto& itr : *this) itr->resetTime(temps);
}

/*! @brief Calculation of exchange coefficients for all boundary conditions in the vector.
 *
 * @param (double temps) the time step for update
 */
inline void Conds_lim::calculer_coeffs_echange(double temps)
{
  for (auto& itr : *this) itr->calculer_coeffs_echange(temps);
}

/*! @brief Returns whether ALL boundary conditions in the vector are compatible with the equation passed as parameter.
 *
 * @param (Equation_base& eqn) the equation with which we will check compatibility
 * @return (int) 1 if all boundary conditions are compatible with the equation, 0 otherwise.
 */
inline int Conds_lim::compatible_avec_eqn(const Equation_base& eqn) const
{
  int ok = 1;
  for (auto& itr : *this) ok *= itr->compatible_avec_eqn(eqn);
  return ok;
}

/*! @brief Returns whether ALL boundary conditions in the vector are compatible with the discretization passed as parameter.
 *
 * @param the discretization with which we will check compatibility
 * @return (int) 1 if all boundary conditions are compatible with the discretization, 0 otherwise.
 */
inline int Conds_lim::compatible_avec_discr(const Discretisation_base& dis) const
{
  int ok = 1;
  for (auto& itr : *this) ok *= itr->compatible_avec_discr(dis);
  return ok;
}

/*! @brief Sets the modifier_val_imp flag for all boundary conditions in the vector.
 *
 * @param (int drap) the value to assign to the flag
 */
inline void Conds_lim::set_modifier_val_imp(int drap)
{
  for (auto& itr : *this) itr->set_modifier_val_imp(drap);
}

#endif /* Conds_lim_included */
