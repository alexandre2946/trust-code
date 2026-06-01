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

#include <Modele_turbulence_scal_Schmidt.h>
#include <Modifier_pour_fluide_dilatable.h>
#include <Convection_Diffusion_std.h>
#include <Param.h>

Implemente_instanciable(Modele_turbulence_scal_Schmidt, "Modele_turbulence_scal_Schmidt", Modele_turbulence_scal_diffturb_base);
// XD schmidt modele_turbulence_scal_base schmidt INHERITS_BRACE The Schmidt model. For the scalar equations, only the
// XD_CONT model based on Reynolds analogy is available. If K_Epsilon was selected in the hydraulic equation, Schmidt
// XD_CONT must be selected for the convection-diffusion temperature equation coupled to the hydraulic equation and
// XD_CONT Schmidt for the concentration equations.

Sortie& Modele_turbulence_scal_Schmidt::printOn(Sortie& s) const { return Modele_turbulence_scal_diffturb_base::printOn(s); }

Entree& Modele_turbulence_scal_Schmidt::readOn(Entree& is)
{
  Modele_turbulence_scal_diffturb_base::readOn(is);
  Cerr << "The value of the turbulent Schmidt number is " << LeScturb_ << finl;
  return is;
}

void Modele_turbulence_scal_Schmidt::set_param(Param& param) const
{
  param.ajouter("ScTurb", &LeScturb_); // XD_ADD_P floattant
  // XD_CONT Keyword to modify the constant (Sct) of Schmlidt model : Dt=Nut/Sct Default value is 0.7.
  Modele_turbulence_scal_base::set_param(param);
}

/*! @brief Returns 1 if the keyword passed as parameter is a field name of the object.
 *
 * @param mot The keyword to compare against known field names.
 * @return 0 if the keyword is not a field name, 1 otherwise.
 */
int Modele_turbulence_scal_Schmidt::comprend_champ(const Motcle& mot) const
{
  if (mot == Motcle("diffusion_turbulente"))
    return 1;
  else
    return 0;
}

/*! @brief Returns 1 if a functional field (Champ_Fonc) with the specified name is owned by the turbulence model, 0 otherwise.
 *
 * @param mot Name of a functional field of the turbulence model.
 * @param ch_ref Reference to the found field (if found).
 * @return 1 if a functional field with the given name was found, 0 otherwise.
 */
int Modele_turbulence_scal_Schmidt::a_pour_Champ_Fonc(const Motcle& mot,
                                                      OBS_PTR(Champ_base) &ch_ref) const
{
  if (mot == Motcle("diffusion_turbulente"))
    {
      ch_ref = diffusivite_turbulente_.valeur();
      return 1;
    }
  return 0;
}

/*! @brief Computes the turbulent coefficient used in the equation and the wall law.
 *
 * @param Unused time parameter.
 */
void Modele_turbulence_scal_Schmidt::mettre_a_jour(double)
{
  calculer_diffusion_turbulente();
  const Milieu_base& mil = equation().probleme().milieu();
  if (loi_paroi_non_nulle())
    loipar_->calculer_scal(diffusivite_turbulente_);

  DoubleTab& lambda_t = conductivite_turbulente_->valeurs();
  lambda_t = diffusivite_turbulente_->valeurs();
  if (equation().probleme().is_dilatable())
    multiplier_par_rho_si_dilatable(lambda_t, mil);
  conductivite_turbulente_->valeurs().echange_espace_virtuel();
  diffusivite_turbulente_->valeurs().echange_espace_virtuel();
}

/*! @brief Computes the turbulent diffusion.
 *
 * turbulent_diffusion = turbulent_viscosity / turbulent_Schmidt_number
 *
 * @return The newly computed turbulent diffusion field.
 * @throws If diffusivite_turbulente and viscosite_turbulente fields do not have the same number of nodal values.
 */
Champ_Fonc_base& Modele_turbulence_scal_Schmidt::calculer_diffusion_turbulente()
{
  DoubleTab& alpha_t = diffusivite_turbulente_->valeurs();
  const DoubleTab& nu_t = la_viscosite_turbulente_->valeurs();
  double temps = la_viscosite_turbulente_->temps();
  int n = alpha_t.size();
  if (nu_t.size() != n)
    {
      Cerr << "The DoubleTab arrays of diffusivite_turbulente and viscosite_turbulente fields" << finl;
      Cerr << "must have the same number of nodal values" << finl;
      exit();
    }

  for (int i = 0; i < n; i++)
    alpha_t[i] = nu_t[i] / LeScturb_;
  diffusivite_turbulente_->changer_temps(temps);
  if (equation().probleme().is_dilatable())
    diviser_par_rho_si_dilatable(diffusivite_turbulente_->valeurs(), equation().probleme().milieu());
  return diffusivite_turbulente_;
}
