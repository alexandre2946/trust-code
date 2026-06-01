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

#ifndef Convection_Diffusion_Turbulent_included
#define Convection_Diffusion_Turbulent_included

#include <Modele_turbulence_scal_base.h>
#include <Equation_base.h>
#include <TRUST_Deriv.h>

class Operateur_Diff;

/*! @brief Mixin class for turbulent convection-diffusion of one or more scalar quantities.
 *
 *      Handles temperature, a single concentration, or multiple concentrations in a fluid
 *      satisfying the incompressibility condition div U = 0, with turbulence modelling.
 *      This class is not, on its own, part of the equation hierarchy (Objet_U/Equation_base);
 *      it is meant to be combined via multiple inheritance to produce turbulent equation classes
 *      (see e.g. Convection_Diffusion_Temperature_Turbulent).
 *
 * @sa Convection_Diffusion_std  Mod_turb_scal, Convection_Diffusion_Temperature_Turbulent, Convection_Diffusion_Concentration_Turbulent
 */
class Convection_Diffusion_Turbulent
{
public :
  Entree& lire_modele(Entree&, const Equation_base& );

  // Now accessed via equation().get_modele(TURBULENCE) to retrieve
  // the turbulence model and then the turbulent diffusion values array
  void completer();
  virtual bool initTimeStep(double dt);
  int preparer_calcul();
  virtual std::vector<YAML_data> data_a_sauvegarder() const;
  virtual int sauvegarder(Sortie&) const;
  virtual int reprendre(Entree&);
  virtual void mettre_a_jour(double);
  virtual ~Convection_Diffusion_Turbulent() {}

protected:
  Entree& lire_op_diff_turbulent(Entree&, const Equation_base&, Operateur_Diff&);
  OWN_PTR(Modele_turbulence_scal_base) le_modele_turbulence;
};

#endif /* Convection_Diffusion_Turbulent_included */
