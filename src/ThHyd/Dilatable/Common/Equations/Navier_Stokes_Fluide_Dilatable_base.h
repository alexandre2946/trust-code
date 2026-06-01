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

#ifndef Navier_Stokes_Fluide_Dilatable_base_included
#define Navier_Stokes_Fluide_Dilatable_base_included

#include <Navier_Stokes_Fluide_Dilatable_Proto.h>
#include <Source_Masse_Fluide_Dilatable_base.h>
#include <TRUSTTabs_forward.h>
#include <Navier_Stokes_std.h>



class Matrice_Morse;

/*! @brief @brief Base class carrying the terms of the momentum equation for a fluid without turbulence modelling under the dilatable fluid assumption.
 *
 * The following form of the Navier-Stokes equations is used:
 *   DU/dt = div(viscous term) - gradP/rho + sources/rho
 *   div U = W
 * where DU/dt is the material derivative of velocity and rho is the density.
 * The implementation allows individual terms (viscous, convective, source) to be neglected.
 * The unknown is the velocity field.
 *
 * @sa Navier_Stokes_std
 */

class Navier_Stokes_Fluide_Dilatable_base : public Navier_Stokes_std, public Navier_Stokes_Fluide_Dilatable_Proto
{
  Declare_base(Navier_Stokes_Fluide_Dilatable_base);

public :
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  int preparer_calcul() override;
  void set_param(Param& param) const override;
  void discretiser() override;
  const Champ_Don_base& diffusivite_pour_transport() const override;
  const Champ_base& diffusivite_pour_pas_de_temps() const override;
  const Champ_base& vitesse_pour_transport() const override;

  // Virtual methods
  DoubleTab& derivee_en_temps_inco(DoubleTab& ) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  void completer() override;
  void assembler( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) override ;
  void assembler_avec_inertie( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) override ;
  void assembler_blocs_avec_inertie(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) override;
  int impr(Sortie& os) const override;
  bool initTimeStep(double dt) override;

  // Inline methods
  inline void mettre_a_jour(double temps) override
  {
    Navier_Stokes_std::mettre_a_jour(temps);
    if (source_masse_)
      source_masse_->mettre_a_jour(temps);
  }

  inline const Champ_Inc_base& rho_la_vitesse() const override { return rho_la_vitesse_; }

  inline bool has_source_masse() const { return bool(source_masse_); }
  inline const Source_Masse_Fluide_Dilatable_base& source_masse() const
  {
    assert(source_masse_);
    return source_masse_.valeur();
  }

protected:
  OWN_PTR(Source_Masse_Fluide_Dilatable_base) source_masse_;
};

#endif /* Navier_Stokes_Fluide_Dilatable_base_included */
