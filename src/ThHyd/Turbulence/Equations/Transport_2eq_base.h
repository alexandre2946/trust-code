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
#ifndef Transport_2eq_base_included
#define Transport_2eq_base_included

#include <Modele_turbulence_hyd_2_eq_base.h>
#include <Equation_base.h>
#include <Operateur_Conv.h>
#include <Operateur_Diff.h>
#include <TRUST_Ref.h>

class Modele_turbulence_hyd_2_eq_base;
class Champ_Inc_base;
class Milieu_base;

/*! @brief Classe Transport_2eq_base Classe de base pour les equations
 *
 *     de transport des modeles a deux equations type k-DISSIP
 *
 */
class Transport_2eq_base: public Equation_base
{
  Declare_base(Transport_2eq_base);

public:

  void set_param(Param&) const override;
  void associer(const Equation_base&);
  inline void associer_vitesse(const Champ_base&);
  void associer_milieu_base(const Milieu_base&) override;
  Milieu_base& milieu() override;
  const Milieu_base& milieu() const override;
  double calculer_pas_de_temps() const override;
  inline const Modele_turbulence_hyd_2_eq_base& modele_turbulence() const;
  inline Modele_turbulence_hyd_2_eq_base& modele_turbulence();
  virtual void associer_modele_turbulence(const Modele_turbulence_hyd_2_eq_base& )=0;

protected:
  void get_position_cells(Nom&, int&);
  void get_position_faces(Nom&, int&);

  Entree& lire_op_diff_turbulent(Entree& is);

  Operateur_Diff terme_diffusif;
  Operateur_Conv terme_convectif;


  bool disable_VEF_mean_value_corrections_ = false;

  OBS_PTR(Milieu_base) le_fluide;
  OBS_PTR(Champ_Inc_base) la_vitesse_transportante;
  OBS_PTR(Modele_turbulence_hyd_2_eq_base) mon_modele;
};

inline void Transport_2eq_base::associer_vitesse(const Champ_base& vit)
{
  la_vitesse_transportante = ref_cast(Champ_Inc_base, vit);
}

/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * (version const)
 *
 * @return (Modele_turbulence_hyd_base&) le modele de turbulence associe a l'equation
 */
inline const Modele_turbulence_hyd_2_eq_base& Transport_2eq_base::modele_turbulence() const
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}

/*! @brief Renvoie le modele de turbulence associe a l'equation.
 *
 * @return (Modele_turbulence_hyd_base&) le modele de turbulence associe a l'equation
 */
inline Modele_turbulence_hyd_2_eq_base& Transport_2eq_base::modele_turbulence()
{
  assert(mon_modele.non_nul());
  return mon_modele.valeur();
}

#endif
