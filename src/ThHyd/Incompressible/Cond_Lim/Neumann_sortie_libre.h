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

#ifndef Neumann_sortie_libre_included
#define Neumann_sortie_libre_included

#include <Neumann_val_ext.h>

/*! @brief Neumann_sortie_libre This class represents an open boundary without imposed velocity.
 *
 *     For Navier-Stokes equations, pressure must be imposed on such a boundary.
 *     To handle the hydraulics, the class Sortie_libre_pression_imposee is derived from Neumann_sortie_libre.
 *     Boundary conditions of type Neumann_sortie_libre or derived types result in zero diffusive fluxes.
 *     However, the treatment of convective fluxes requires knowing the convected field outside the boundary in case of fluid re-entry.
 *     This is why the class carries an OWN_PTR(Champ_front_base) (member le_champ_ext).
 *
 *     In the computation operators, boundary conditions of type Neumann_sortie_libre and derived types will be treated identically.
 *
 * @sa Neumann Sortie_libre_pression_imposee
 */
class Neumann_sortie_libre: public Neumann_val_ext
{
  Declare_instanciable(Neumann_sortie_libre);
public:
  const DoubleTab& tab_ext() const override;
  DoubleTab& tab_ext() override;

  double val_ext(int i) const override;
  double val_ext(int i, int j) const override;
  const DoubleTab& val_ext() const;
  int initialiser(double temps) override;
  void associer_fr_dis_base(const Frontiere_dis_base&) override;
  void verifie_ch_init_nb_comp() const override;

  void fixer_nb_valeurs_temporelles(int nb_cases) override;
  void mettre_a_jour(double temps) override;
  void set_temps_defaut(double temps) override;
  void changer_temps_futur(double temps, int i) override;
  int avancer(double temps) override;
  int reculer(double temps) override;

protected:
  OWN_PTR(Champ_front_base) le_champ_ext;
  mutable DoubleTab val_ext_; // Stores all boundary condition values on all faces of the boundary (no assumption of a uniform field). Useful for GPU.
};

#endif
