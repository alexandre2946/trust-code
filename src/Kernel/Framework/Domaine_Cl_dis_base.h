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

#ifndef Domaine_Cl_dis_base_included
#define Domaine_Cl_dis_base_included

#include <Conds_lim.h>
#include <TRUST_Ref.h>
#include <MorEqn.h>
#include <Domaine_forward.h>

/*! @brief class Domaine_Cl_dis_base Domaine_Cl_dis_base objects represent discretized boundary conditions
 *
 *      . Domaine_Cl_dis_base is a piece of equation so it
 *      inherits from MorEqn, it is through this inheritance that each
 *      Domaine_Cl_dis_base object contains a reference to the equation to which it
 *      refers. The discretized boundary conditions represent the
 *      boundary conditions of the discretized domain associated with the equation
 *      referenced by Domaine_Cl_dis_base.
 *      Domaine_Cl_dis_base has a member representing the boundary conditions.
 *
 * @sa MorEqn Conds_lim Cond_lim_base, Abstract class., Abstract methods:, void completer(const Domaine_dis_base& ), void imposer_cond_lim(Champ_Inc_base& )
 */
class Domaine_Cl_dis_base : public MorEqn, public Objet_U
{

  Declare_base(Domaine_Cl_dis_base);

public:

  // Method overridden from Objet_U:
  void nommer(const Nom& nom) override;
  const Cond_lim&  les_conditions_limites(int ) const;
  Cond_lim&        les_conditions_limites(int );
  Conds_lim&       les_conditions_limites();
  const Conds_lim& les_conditions_limites() const;
  int           nb_cond_lim() const;
  virtual void     mettre_a_jour(double temps);
  virtual void     mettre_a_jour_ss_pas_dt(double temps);
  void             resetTime(double time);
  void             completer();
  int           contient_Cl(const Nom&);
  Domaine_dis_base&        domaine_dis();
  const Domaine_dis_base&  domaine_dis() const;
  Domaine&            domaine();
  const Domaine&      domaine() const;

  virtual void associer(const Domaine_dis_base& ddb) { /* Does nothing by default */ }

  virtual int   calculer_coeffs_echange(double temps);
  // !SC: passing OWN_PTR(Champ_Inc_base) is no longer necessary since there is a ref now
  virtual void     imposer_cond_lim(Champ_Inc_base&,double ) = 0;
  int           nb_faces_Cl() const;

  virtual const Cond_lim& la_cl_de_la_face(int num_face) const;

  const Cond_lim_base& condition_limite_de_la_face_reelle(int face_globale, int& face_locale) const;
  const Cond_lim_base& condition_limite_de_la_face_virtuelle(int face_globale, int& face_locale) const;
  const Cond_lim_base& condition_limite_de_la_frontiere(Nom frontiere) const;
  Cond_lim_base& condition_limite_de_la_frontiere(Nom frontiere);
  void changer_temps_futur(double temps,int i);
  void calculer_derivee_en_temps(double t1, double t2);
  void set_temps_defaut(double temps);
  int avancer(double temps);
  int reculer(double temps);
  virtual int initialiser(double temps);
  inline const Nom& le_nom() const override;

  virtual void associer_inconnue(const Champ_Inc_base&);
  virtual const Champ_Inc_base& inconnue() const;
  virtual Champ_Inc_base& inconnue();

protected:

  Nom nom_;
  Conds_lim  les_conditions_limites_;
  OBS_PTR(Champ_Inc_base) mon_inconnue;
  virtual void completer(const Domaine_dis_base& ) = 0;
};

inline const Nom& Domaine_Cl_dis_base::le_nom() const
{
  return nom_;
}

#endif

