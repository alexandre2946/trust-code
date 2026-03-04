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
#include <Transport_2eq_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <communications.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Domaine_VF.h>
#include <Param.h>
#include <Debog.h>

Implemente_base(Transport_2eq_base, "Transport_2eq_base", Equation_base);


/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_2eq_base::printOn(Sortie& is) const
{
  return is << que_suis_je() << "\n";
}

/*! @brief Simple appel a Equation_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_2eq_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  return is;
}

void Transport_2eq_base::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("diffusion", (this));
  param.ajouter_non_std("convection", (this));
  param.ajouter_condition("is_read_diffusion", "The diffusion operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_condition("is_read_convection", "The convection operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_flag("disable_VEF_mean_value_corrections", &disable_VEF_mean_value_corrections_); // X_D_ADD_P flag Flag to prevent usage of mean value of neighbors for corrections applied to k-eps and k-omega models. This only applies to VEF (mean values not used in VDF). See methods 'Transport_K_Eps_base::controler_K_Eps' and 'Transport_K_Omega_base::controler_K_Omega' for more details.
}

// Lecture et typage de l'operateur diffusion turbulente.
// Attention : il faut avoir fait "terme_diffusif.associer_diffusivite" avant d'enter ici.
Entree& Transport_2eq_base::lire_op_diff_turbulent(Entree& is)
{
  Motcle accouverte = "{", accfermee = "}";
  Nom type = "Op_Diff_K_";

  Nom eps_or_omega = que_suis_je();
  if (eps_or_omega.contient("Omega")) type += "Omega_";
  else type += "Eps_";

  Nom qc = modele_turbulence().equation().que_suis_je();
  Cerr << ">>>>>>> Nom eq = " << qc << finl;
  if (qc=="Navier_Stokes_QC" || qc == "Navier_Stokes_Turbulent_QC")
    {
      type+="QC_";
    }

  Nom discr = discretisation().que_suis_je();
  // les operateurs de diffusion sont communs aux discretisations VEF et VEFP1B
  if (discr == "VEFPreP1B") discr = "VEF";

  type += discr;

  Nom nb_inc;
  if (terme_diffusif.diffusivite().nb_comp() == 1)
    nb_inc = "_";
  else
    nb_inc = "_Multi_inco_";
  type += nb_inc;

  Nom type_inco = inconnue().que_suis_je();
  if (type_inco == "Champ_Q1_EF") type += "Q1";
  else type += (type_inco.suffix("Champ_"));

  if (axi) type += "_Axi";

  Motcle motbidon;
  is >> motbidon;
  if (motbidon != accouverte)
    {
      Cerr << "A { was expected while reading the turbulent diffusive term" << finl;
      Process::exit();
    }
  is >> motbidon;
  if (motbidon == "negligeable")
    {
      type = "Op_Dift_negligeable";
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
      is >> motbidon;
      //on lit la fin de diffusion { }
      if (motbidon != accfermee)
        Cerr << " On ne peut plus entrer d option apres negligeable " << finl;
    }
  else if (motbidon == accfermee)
    {
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      Cerr << terme_diffusif->que_suis_je() << finl;
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
    }
  else
    {
      type += motbidon;
      is >> motbidon;
      if (motbidon != accfermee)
        Cerr << " No option are now readable for the turbulent diffusive term" << finl;

      if (discr == "VEF")
        type += "_P1NC";
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      Cerr << terme_diffusif->que_suis_je() << finl;
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
    }
  return is;
}

//For VEF-like scheme
void Transport_2eq_base::get_position_faces(Nom& position, int& n)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());

  position = "x=";
  position += (Nom)domaine_vf.xv(n, 0);
  position += " y=";
  position += (Nom)domaine_vf.xv(n, 1);
  if (dimension == 3)
    {
      position += " z=";
      position += (Nom)domaine_vf.xv(n, 2);
    }
}

// For VDF-like scheme
void Transport_2eq_base::get_position_cells(Nom& position, int& n)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());

  position = "x=";
  position += (Nom)domaine_vf.xp(n, 0);
  position += " y=";
  position += (Nom)domaine_vf.xp(n, 1);
  if (dimension == 3)
    {
      position += " z=";
      position += (Nom)domaine_vf.xp(n, 2);
    }
}

/*! @brief Associe un milieu physique a l'equation.
 *
 * @param (Milieu_base& un_milieu) le milieu physique a associer a l'equation
 */
void Transport_2eq_base::associer_milieu_base(const Milieu_base& un_milieu)
{
  le_fluide =  un_milieu;
}

/*! @brief Renvoie le milieu (fluide) associe a l'equation.
 *
 * @return (Milieu_base&) le milieu (fluide) associe a l'equation
 */
Milieu_base& Transport_2eq_base::milieu()
{
  if(!le_fluide.non_nul())
    {
      Cerr << "No fluid has been associated to the two equations Transport"
           << que_suis_je() << " equation." << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

/*! @brief Renvoie le milieu (fluide) associe a l'equation.
 *
 * (version const)
 *
 * @return (Milieu_base&) le milieu (fluide) associe a l'equation
 */
const Milieu_base& Transport_2eq_base::milieu() const
{
  if(!le_fluide.non_nul())
    {
      Cerr << "No fluid has been associated to the Transport K_Epsilon"
           << que_suis_je() << " equation." << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

void Transport_2eq_base::associer(const Equation_base& eqn_hydr)
{
  Equation_base::associer_pb_base(eqn_hydr.probleme());
  Equation_base::associer_sch_tps_base(eqn_hydr.schema_temps());
  Equation_base::associer_domaine_dis(eqn_hydr.domaine_dis());
}

double Transport_2eq_base::calculer_pas_de_temps() const
{
  // on prend le pas de temps de l'eq de NS.
  return probleme().equation(0).calculer_pas_de_temps();
}
