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

#include <Convection_Diffusion_Temperature_Turbulent.h>
#include <Fluide_Incompressible.h>
#include <Param.h>

Implemente_instanciable(Convection_Diffusion_Temperature_Turbulent, "Convection_Diffusion_Temperature_Turbulent", Convection_Diffusion_Temperature);
// XD convection_diffusion_temperature_turbulent eqn_base convection_diffusion_temperature_turbulent INHERITS_BRACE
// XD_CONT Energy equation (temperature diffusion convection) as well as the associated turbulence model equations.
// XD attr modele_turbulence modele_turbulence_scal_base modele_turbulence OPT Turbulence model for the energy equation.

Sortie& Convection_Diffusion_Temperature_Turbulent::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Convection_Diffusion_Temperature_Turbulent::readOn(Entree& is) { return Convection_Diffusion_Temperature::readOn(is); }

void Convection_Diffusion_Temperature_Turbulent::set_param(Param& param) const
{
  Convection_Diffusion_Temperature::set_param(param);
  param.ajouter_non_std("modele_turbulence", (this), Param::REQUIRED);
}

int Convection_Diffusion_Temperature_Turbulent::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      lire_op_diff_turbulent(is, *this, terme_diffusif);
      // the field for dt_stab is the same as that of the operator
      terme_diffusif.associer_diffusivite_pour_pas_de_temps(diffusivite_pour_pas_de_temps());
      return 1;
    }
  else if (mot == "modele_turbulence")
    {
      lire_modele(is, *this);
      RefObjU le_modele;
      le_modele = le_modele_turbulence.valeur();
      liste_modeles_.add_if_not(le_modele);
      return 1;
    }
  else
    return Convection_Diffusion_Temperature::lire_motcle_non_standard(mot, is);
}

/*! @brief for PDI IO: retrieve name, type and dimensions of the fields to save/restore
 *
 */
std::vector<YAML_data> Convection_Diffusion_Temperature_Turbulent::data_a_sauvegarder() const
{
  std::vector<YAML_data> data = Convection_Diffusion_Temperature::data_a_sauvegarder();
  std::vector<YAML_data> turb = Convection_Diffusion_Turbulent::data_a_sauvegarder();
  data.insert(data.end(), turb.begin(), turb.end());
  return data;
}

/*! @brief Saves to an output stream via a double call to: Convection_Diffusion_Temperature::sauvegarder(Sortie&);
 *
 *       Convection_Diffusion_Turbulent::sauvegarder(Sortie&);
 *
 * @param os Output stream.
 * @return Always returns 1.
 */
int Convection_Diffusion_Temperature_Turbulent::sauvegarder(Sortie& os) const
{
  int bytes = 0;
  bytes += Convection_Diffusion_Temperature::sauvegarder(os);
  bytes += Convection_Diffusion_Turbulent::sauvegarder(os);
  return bytes;
}

/*! @brief Restores from an input stream via a double call to: Convection_Diffusion_Temperature::reprendre(Entree&);
 *
 *       Convection_Diffusion_Turbulent::reprendre(Entree&);
 *
 * @param is Input stream.
 * @return Always returns 1.
 */
int Convection_Diffusion_Temperature_Turbulent::reprendre(Entree& is)
{
  Convection_Diffusion_Temperature::reprendre(is);
  Convection_Diffusion_Turbulent::reprendre(is);
  return 1;
}

/*! @brief Double call to: Convection_Diffusion_Turbulent::completer()
 *
 *      and Convection_Diffusion_Temperature::completer()
 *
 */
void Convection_Diffusion_Temperature_Turbulent::completer()
{
  Convection_Diffusion_Turbulent::completer();
  Convection_Diffusion_Temperature::completer();
}

/*! @brief Time update of the equation via a double call to: Convection_Diffusion_Temperature::mettre_a_jour(double);
 *
 *       and Convection_Diffusion_Turbulent::mettre_a_jour(double);
 *
 * @param temps Current time.
 */
void Convection_Diffusion_Temperature_Turbulent::mettre_a_jour(double temps)
{
  Convection_Diffusion_Temperature::mettre_a_jour(temps);
  Convection_Diffusion_Turbulent::mettre_a_jour(temps);
}

void Convection_Diffusion_Temperature_Turbulent::creer_champ(const Motcle& motlu)
{
  Convection_Diffusion_Temperature::creer_champ(motlu);

  if (le_modele_turbulence)
    le_modele_turbulence->creer_champ(motlu);
}

bool Convection_Diffusion_Temperature_Turbulent::has_champ(const Motcle& nom, OBS_PTR(Champ_base)& ref_champ) const
{
  if (Convection_Diffusion_Temperature::has_champ(nom, ref_champ))
    return true;

  if (le_modele_turbulence)
    if (le_modele_turbulence->has_champ(nom, ref_champ))
      return true;

  return false; /* nothing found */
}

bool Convection_Diffusion_Temperature_Turbulent::has_champ(const Motcle& nom) const
{
  if (Convection_Diffusion_Temperature::has_champ(nom))
    return true;

  if (le_modele_turbulence)
    if (le_modele_turbulence->has_champ(nom))
      return true;

  return false; /* nothing found */
}

const Champ_base& Convection_Diffusion_Temperature_Turbulent::get_champ(const Motcle& nom) const
{
  OBS_PTR(Champ_base) ref_champ;

  if (Convection_Diffusion_Temperature::has_champ(nom, ref_champ))
    return ref_champ;

  if (le_modele_turbulence)
    if (le_modele_turbulence->has_champ(nom, ref_champ))
      return ref_champ;

  throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));
}

void Convection_Diffusion_Temperature_Turbulent::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Convection_Diffusion_Temperature::get_noms_champs_postraitables(nom, opt);
  if (le_modele_turbulence)
    le_modele_turbulence->get_noms_champs_postraitables(nom, opt);
}

/*! @brief Double call to: Convection_Diffusion_Turbulent::preparer_calcul()
 *
 *       and Convection_Diffusion_Temperature::preparer_calcul()
 *
 * @return Always returns 1.
 */
int Convection_Diffusion_Temperature_Turbulent::preparer_calcul()
{
  Convection_Diffusion_Temperature::preparer_calcul();
  Convection_Diffusion_Turbulent::preparer_calcul();
  return 1;
}

bool Convection_Diffusion_Temperature_Turbulent::initTimeStep(double dt)
{
  bool ok = Convection_Diffusion_Turbulent::initTimeStep(dt);
  ok = ok && Convection_Diffusion_Temperature::initTimeStep(dt);
  return ok;
}

void Convection_Diffusion_Temperature_Turbulent::imprimer(Sortie& os) const
{
  Convection_Diffusion_Temperature::imprimer(os);
  le_modele_turbulence->imprimer(os);
}

const RefObjU& Convection_Diffusion_Temperature_Turbulent::get_modele(Type_modele type) const
{
  for (const auto &itr : liste_modeles_)
    {
      const RefObjU& mod = itr;
      if (mod)
        if ((sub_type(Modele_turbulence_scal_base, mod.valeur())) && (type == TURBULENCE))
          return mod;
    }
  return Equation_base::get_modele(type);
}
