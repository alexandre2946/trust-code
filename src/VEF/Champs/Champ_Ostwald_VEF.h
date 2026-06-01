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

#ifndef Champ_Ostwald_VEF_included
#define Champ_Ostwald_VEF_included

#include <Champ_Ostwald.h>
#include <Champ_P1NC.h>
#include <TRUST_Ref.h>

class Navier_Stokes_std;
class Domaine_VEF;

/*! @brief class Champ_Ostwald_VEF
 *
 *  @brief Represents a field in VEF discretization that varies as a function
 *         of the consistency and the structure index.
 *         Field used for the Ostwald fluid with VEF discretization.
 *         References domaine_VEF to use the correct domain with domaine_dis_base,
 *                    Champ_P1NC to compute D:D,
 *                    Navier_Stokes_std to access the thermo-hydraulic equation
 *                    and use one of its unknowns: the velocity.
 *
 * @sa Champ_Ostwald
 */
class Champ_Ostwald_VEF : public Champ_Ostwald
{
  Declare_instanciable(Champ_Ostwald_VEF);
public :
  void associer_eqn(const Navier_Stokes_std& );
  void mettre_a_jour(double temps) override;
  void me_calculer(double tps) override;
  int initialiser(const double temps) override;
  void init_mu(DoubleTab& );
  void calculer_dscald(DoubleTab&);

  inline const Champ_P1NC& mon_champs() const { return mon_champ_.valeur(); }
  inline void associer_champ(const Champ_P1NC& un_champ) { mon_champ_ = un_champ; }

protected :
  void calculer_mu(DoubleTab& );
  OBS_PTR(Champ_P1NC) mon_champ_;
  OBS_PTR(Navier_Stokes_std) eq_hydraulique;
};

#endif /* Champ_Ostwald_VEF_included */
