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

#ifndef Champ_front_base_included
#define Champ_front_base_included

#include <Field_base.h>
#include <Champ_Proto.h>
#include <TRUST_Ref.h>
#include <Roue.h>

class Frontiere_dis_base;
class Champ_Inc_base;
class Domaine_dis_base;
class Cond_lim_base;

/*! @brief class Champ_front_base Base class for the hierarchy of boundary fields.
 *
 * A
 *      Champ_front_base object defines a field on the boundary of a
 *      domain. A Champ_front_base type object will be associated with
 *      each boundary condition.
 *      Champ_front_base derives from Champ_Proto, in order to have an interface
 *      conforming to all fields, and to inherit common operations on
 *      Fields.
 *      The two main methods are initialiser and
 *      mettre_a_jour. These are the only two that can modify the
 *      field values.
 *      The initialiser method is called once at the beginning of
 *      the calculation. It must not depend on data external to
 *      the equation carrying the BC (indeed, nothing guarantees that this
 *      data is initialized). On the other hand, the values of
 *      the unknown on which the BC depends can be useful for
 *      initialization => the unknown is passed as a parameter in
 *      read-only mode.
 *      The mettre_a_jour method is called at the beginning of each time step
 *      or sub-time-step, it can use data
 *      external to the equation. It is the responsibility of the algorithm to ensure
 *      that this data is relevant...
 *      In the case of stationary fields, the mettre_a_jour method
 *      has nothing to do and the values are filled once and
 *      for all by the initialiser method.
 *      In the case of unsteady fields, there are several values in time
 *      and each can be updated.
 *      Champ_front_base are divided into:
 *       * Champ_front_uniforme, constant in time and space
 *       * Champ_front_instationnaire_base, uniform in space but variable in time
 *       * Champ_front_var, variable in space.
 *      Champ_front_var are then classified according to whether they are
 *      stationary or unsteady.
 *
 *      The values are stored in a wheel of DoubleTab.
 *      If the field is uniform in space, the DoubleTabs are
 *      dimensioned to 1.
 *      If it is stationary, the wheel has only one temporal value
 *      and the time assigned to it is meaningless.
 *      If it is unsteady, the temporal values are those
 *      of the unknown of the equation to which the field refers.
 *
 *
 * @sa Champ_Proto Frontiere_dis_base, Abstract class, Abstract method:, Champ_front_base& affecter_(const Champ_front_base& ch)
 */
class Champ_front_base : public Field_base, public Champ_Proto
{
  Declare_base_sans_constructeur(Champ_front_base);
public:
  Champ_front_base();
  inline virtual void completer() { }
  virtual int initialiser(double temps, const Champ_Inc_base& inco);
  virtual void associer_fr_dis_base(const Frontiere_dis_base&);
  using Champ_Proto::valeurs;
  inline virtual DoubleTab& valeurs() override;
  inline virtual const DoubleTab& valeurs() const override;
  virtual bool has_valeurs_au_temps(double temps) const { return true; }
  virtual DoubleTab& valeurs_au_temps(double temps)=0;
  virtual const DoubleTab& valeurs_au_temps(double temps) const = 0;
  virtual const Frontiere_dis_base& frontiere_dis() const;
  virtual Frontiere_dis_base& frontiere_dis();
  virtual const Domaine_dis_base& domaine_dis() const;
  virtual Champ_front_base& affecter_(const Champ_front_base& ch) = 0;
  virtual void fixer_nb_valeurs_temporelles(int nb_cases);
  virtual void mettre_a_jour(double temps);
  virtual void calculer_coeffs_echange(double temps);
  virtual void valeurs_face(int, DoubleVect&) const;
  virtual inline void verifier(const Cond_lim_base& la_cl) const;
  virtual double get_temps_defaut() const { return temps_defaut; }
  virtual void set_temps_defaut(double temps) { temps_defaut = temps; }
  virtual void changer_temps_futur(double temps, int i);
  virtual int avancer(double temps);
  virtual int reculer(double temps);
  virtual bool instationnaire() const { return instationnaire_; }
  virtual void set_instationnaire(bool flag) { instationnaire_ = flag; Gpoint_ = valeurs(); Gpoint_=0; } // Dimensionne Gpoint_
  virtual inline void set_derivee_en_temps(DoubleTab& Gpoint) { Gpoint_ = Gpoint; }
  virtual inline const DoubleTab& derivee_en_temps() const { return Gpoint_; }
  virtual void calculer_derivee_en_temps(double t1, double t2);

protected:
  double temps_defaut ; // The default time used when the parameter
  // is not specified. In particular, this is
  // the time used by operators and solvers.
  OBS_PTR(Frontiere_dis_base) la_frontiere_dis;
  Roue_ptr les_valeurs; // The field values
  DoubleTab Gpoint_; // Time derivative of the boundary condition values
private:
  bool instationnaire_ = false; // By default, stationary field
};


/*! @brief Returns the discretized boundary associated with the field.
 *
 * (const version)
 *
 * @return (Frontiere_dis_base&) the discretized boundary associated with the field
 */
inline const Frontiere_dis_base& Champ_front_base::frontiere_dis() const
{
  return la_frontiere_dis.valeur();
}


/*! @brief Returns the discretized boundary associated with the field.
 *
 * @return (Frontiere_dis_base&) the discretized boundary associated with the field
 */
inline Frontiere_dis_base& Champ_front_base::frontiere_dis()
{
  return la_frontiere_dis.valeur();
}


/*! @brief Returns the array of field values.
 *
 * @return (DoubleTab&) the array of field values at the default time.
 */
inline DoubleTab& Champ_front_base::valeurs()
{
  return valeurs_au_temps(temps_defaut);
}


/*! @brief Returns the array of field values at the default time.
 *
 * (const version)
 *
 * @return (DoubleTab&) the array of field values
 */
inline const DoubleTab& Champ_front_base::valeurs() const
{
  return valeurs_au_temps(temps_defaut);
}

// must be overridden by champ_front classes that want to verify the arguments
inline void Champ_front_base::verifier(const Cond_lim_base& la_cl) const
{
  return;
}

#endif
