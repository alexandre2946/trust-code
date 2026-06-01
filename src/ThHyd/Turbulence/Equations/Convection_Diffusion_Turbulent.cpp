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

#include <Convection_Diffusion_Turbulent.h>
#include <Modele_turbulence_scal_base.h>
#include <Operateur_Diff_base.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Operateur_Diff.h>
#include <Probleme_base.h>

#include <Domaine.h>
#include <Avanc.h>

Entree& Convection_Diffusion_Turbulent::lire_modele(Entree& is, const Equation_base& eqn)
{
  Cerr << "Reading and typing of the turbulence model : ";

  Motcle typ;
  is >> typ;

  Motcle nom1("Modele_turbulence_scal_");
  nom1 += typ;

  if (typ == "sous_maille_dyn")
    {
      nom1 += "_";
      Nom disc = eqn.discretisation().que_suis_je();
      if (disc == "VEFPreP1B") disc = "VEF";
      nom1 += disc;
    }
  Cerr << nom1 << finl;

  le_modele_turbulence.typer(nom1);
  le_modele_turbulence->associer_eqn(eqn);
  le_modele_turbulence->associer(eqn.domaine_dis(), eqn.domaine_Cl_dis());

  is >> le_modele_turbulence.valeur(); // On lit :-)

  le_modele_turbulence->discretiser();

  return is;
}

Entree& Convection_Diffusion_Turbulent::lire_op_diff_turbulent(Entree& is, const Equation_base& eqn, Operateur_Diff& terme_diffusif)
{
  Motcle accfermee = "}";
  Nom type = "Op_Dift_";

  Motcle motbidon;
  is >> motbidon;
  assert(motbidon == "{");
  is >> motbidon;
  if (motbidon == "negligeable")
    {
      type += motbidon;
      is >> motbidon;

      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(eqn);
      Cerr << terme_diffusif->que_suis_je() << finl;
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());

      assert(motbidon == accfermee);
    }
  else if (motbidon == "stab")
    {
      Nom disc = eqn.discretisation().que_suis_je();
      if (disc == "VEFPreP1B")
        disc = "VEF";

      Nom nb_inc;
      nb_inc = "_";
      if ((terme_diffusif.diffusivite().nb_comp() != 1) && ((disc == "VDF")))
        nb_inc = "_Multi_inco_";

      Nom type_inco = eqn.inconnue().que_suis_je();

      type += disc;
      type += nb_inc;
      type += (type_inco.suffix("Champ_"));
      if (Objet_U::axi)
        type += "_Axi";
      type += "_stab";

      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(eqn);
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
      is >> terme_diffusif.valeur();
    }
  else
    {
      Nom disc = eqn.discretisation().que_suis_je();
      // the C_D_Turb_T operators are common to VEF and VEFP1B discretizations
      if (disc == "VEFPreP1B")
        disc = "VEF";
      type += disc;

      Nom nb_inc;
      nb_inc = "_";
      if ((terme_diffusif.diffusivite().nb_comp() != 1) && ((disc == "VDF")))
        nb_inc = "_Multi_inco_";

      type += nb_inc;

      Nom type_inco = eqn.inconnue().que_suis_je();
      if (type_inco == "Champ_Q1_EF")
        type += "Q1";
      else
        type += (type_inco.suffix("Champ_"));

      if (Objet_U::axi)
        type += "_Axi";

      if (motbidon != accfermee)
        type += Nom("_") + motbidon, is >> motbidon;
      assert(motbidon == accfermee);

      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(eqn);
      Cerr << terme_diffusif->que_suis_je() << finl;
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
    }
  return is;
}

/*! @brief Completes the turbulence model.
 *
 */
void Convection_Diffusion_Turbulent::completer()
{
  le_modele_turbulence->completer();
  le_modele_turbulence->loi_paroi().completer();
}

bool Convection_Diffusion_Turbulent::initTimeStep(double dt)
{
  return le_modele_turbulence->initTimeStep(dt);
}

/*! @brief Prepares the computation.
 *
 * @return Always returns 1.
 */
int Convection_Diffusion_Turbulent::preparer_calcul()
{
  le_modele_turbulence->preparer_calcul();
  return 1;
}

/*! @brief for PDI IO: retrieve name, type and dimensions of the fields to save/restore
 *
 */
std::vector<YAML_data> Convection_Diffusion_Turbulent::data_a_sauvegarder() const
{
  std::vector<YAML_data> data = le_modele_turbulence->data_a_sauvegarder();
  return data;
}

/*! @brief Simple call to Modele_turbulence_scal_base::sauvegarder(Sortie&) on the turbulence member.
 *
 * Saves the turbulence model to an output stream.
 *
 * @param os Output stream.
 * @return Return code propagated from the turbulence model.
 */
int Convection_Diffusion_Turbulent::sauvegarder(Sortie& os) const
{
  return le_modele_turbulence->sauvegarder(os);
}

/*! @brief Restores from a checkpoint via an input stream.
 *
 * @param is Input stream.
 * @return Always returns 1.
 * @throws If end of file is reached during the restore.
 */
int Convection_Diffusion_Turbulent::reprendre(Entree& is)
{
  le_modele_turbulence->reprendre(is);
  return 1;
}

/*! @brief Time update of the turbulence model.
 *
 * @param temps Current time.
 */
void Convection_Diffusion_Turbulent::mettre_a_jour(double temps)
{
  le_modele_turbulence->mettre_a_jour(temps);
}
