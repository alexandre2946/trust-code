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

#include <Convection_Diffusion_Concentration.h>
#include <Frontiere_dis_base.h>
#include <Navier_Stokes_std.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Constituant.h>
#include <EChaine.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Convection_Diffusion_Concentration,"Convection_Diffusion_Concentration",Convection_Diffusion_std);
// XD convection_diffusion_concentration eqn_base convection_diffusion_concentration INHERITS_BRACE Constituent
// XD_CONT transport vectorial equation (concentration diffusion convection).

Convection_Diffusion_Concentration::Convection_Diffusion_Concentration():nb_constituants_(-1), masse_molaire_(-1.) { }

Sortie& Convection_Diffusion_Concentration::printOn(Sortie& is) const { return Convection_Diffusion_std::printOn(is); }

/*! @brief Verifies that the equation has an associated concentration and constituent, then calls Convection_Diffusion_std::readOn(Entree&).
 *
 * @param is an input stream
 * @return the modified input stream
 */
Entree& Convection_Diffusion_Concentration::readOn(Entree& is)
{
  assert(la_concentration);
  Convection_Diffusion_std::readOn(is);
  if (terme_convectif.op_non_nul())
    {
      Nom nom="Convection_";
      nom+=inconnue().le_nom(); // append the unknown name to plan for passive scalar equations
      terme_convectif.set_fichier(nom);
      terme_convectif.set_description((Nom)"Convective mass transfer rate=Integral(-C*u*ndS)[m"+(Nom)(dimension+bidim_axi)+".Mol.s-1]");
    }
  else
    {
      EChaine ech("{ negligeable }");
      ech >> terme_convectif;
    }
  if (terme_diffusif.op_non_nul())
    {
      Nom nom="Diffusion_";
      nom+=inconnue().le_nom();
      terme_diffusif.set_fichier(nom);
      terme_diffusif.set_description((Nom)"Diffusion mass transfer rate=Integral(alpha*grad(C)*ndS) [m"+(Nom)(dimension+bidim_axi)+".Mol.s-1]");
    }
  else
    {
      EChaine ech("{ negligeable }");
      ech >> terme_diffusif;
    }
  return is;
}

// returns the molar mass; it must have been read beforehand
const double& Convection_Diffusion_Concentration::masse_molaire() const
{
  assert(masse_molaire_>0);
  return masse_molaire_;
}

void Convection_Diffusion_Concentration::set_param(Param& param) const
{
  Convection_Diffusion_std::set_param(param);
  param.ajouter_non_std("nom_inconnue",(this)); // XD_ADD_P chaine
  // XD_CONT Keyword Nom_inconnue will rename the unknown of this equation with the given name. In the postprocessing
  // XD_CONT part, the concentration field will be accessible with this name. This is usefull if you want to track more
  // XD_CONT than one concentration (otherwise, only the concentration field in the first concentration equation can be
  // XD_CONT accessed).
  param.ajouter_non_std("alias",(this)); // XD_ADD_P chaine
  // XD_CONT not_set
  param.ajouter("masse_molaire",&masse_molaire_); // XD_ADD_P floattant
  // XD_CONT not_set
  param.ajouter_non_std("is_multi_scalar|is_multi_scalar_diffusion", (this)); // XD_ADD_P rien
  // XD_CONT Flag to activate the multi_scalar diffusion operator
}

int Convection_Diffusion_Concentration::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="nom_inconnue")
    {
      Motcle nom; // Question: should it be uppercased?
      is >> nom;
      Cerr << "The unknow of a Convection_Diffusion_Concentration equation is renamed"
           << "\n Old name : " << inconnue().le_nom()
           << "\n New name : " << nom << finl;
      inconnue().nommer(nom);
      champs_compris_.ajoute_champ(la_concentration);
      return 1;
    }
  else if (mot=="alias")
    {
      Motcle nom; // Question: veut-on le mettre en majuscules ?
      is >> nom;
      Cerr << "nom_inconnue: On renomme l'equation et son inconnue"
           << "\n Ancien nom : " << inconnue().le_nom()
           << "\n Nouveau nom : " << nom << finl;
      inconnue().nommer(nom);
      champs_compris_.ajoute_champ(la_concentration);
      return 1;
    }
  else if (mot=="is_multi_scalar" || mot=="is_multi_scalar_diffusion")
    {
      diffusion_multi_scalaire_ = true;
      return 1;
    }
  else
    return Convection_Diffusion_std::lire_motcle_non_standard(mot,is);
}

/*! @brief Associates a physical medium to the equation; the medium is cast to Constituant and associated with the equation.
 *
 * @param un_milieu the physical medium to associate
 * @throws diffusivity of the constituent in the fluid not defined
 */
void Convection_Diffusion_Concentration::associer_milieu_base(const Milieu_base& un_milieu)
{
  const Constituant& un_constituant = ref_cast(Constituant,un_milieu);
  associer_constituant(un_constituant);
}

const Champ_Don_base& Convection_Diffusion_Concentration::diffusivite_pour_transport() const
{
  return constituant().diffusivite_constituant();
}

/*! @brief Discretizes the equation.
 *
 */
void Convection_Diffusion_Concentration::discretiser()
{
  assert(le_constituant);
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  Cerr << "Transport concentration(s) equation discretization " << finl;
  nb_constituants_ = constituant().nb_constituants();
  dis.concentration(schema_temps(), domaine_dis(), la_concentration, nb_constituants_);
  champs_compris_.ajoute_champ(la_concentration);
  Equation_base::discretiser();
}

/*! @brief Returns the physical medium of the equation.
 *
 * (a Constituant upcast to Milieu_base)
 *     (const version)
 *
 * @return the Constituant upcast to Milieu_base
 */
const Milieu_base& Convection_Diffusion_Concentration::milieu() const
{
  return constituant();
}


/*! @brief Returns the physical medium of the equation.
 *
 * (a Constituant upcast to Milieu_base)
 *
 * @return the Constituant upcast to Milieu_base
 */
Milieu_base& Convection_Diffusion_Concentration::milieu()
{
  return constituant();
}


/*! @brief Returns the constituent (if one has been associated).
 *
 * (const version)
 *
 * @return the constituent associated with the equation
 */
const Constituant& Convection_Diffusion_Concentration::constituant() const
{
  if(!le_constituant)
    {
      Cerr << "You forgot to associate the constituent to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_constituant.valeur();
}


/*! @brief Returns the constituent (if one has been associated).
 *
 * @return the constituent associated with the equation
 */
Constituant& Convection_Diffusion_Concentration::constituant()
{
  if(!le_constituant)
    {
      Cerr << "No constituant has been associated "
           << "with a Convection_Diffusion_Concentration equation." << finl;
      exit();
    }
  return le_constituant.valeur();
}

int Convection_Diffusion_Concentration::preparer_calcul()
{
  Equation_base::preparer_calcul();
  double temps=schema_temps().temps_courant();
  constituant().initialiser(temps);
  return 1;
}

void Convection_Diffusion_Concentration::mettre_a_jour(double temps)
{
  Equation_base::mettre_a_jour(temps);
  constituant().mettre_a_jour(temps);
}

/*! @brief Prints the boundary fluxes to an output stream.
 *
 * Calls Equation_base::impr(Sortie&)
 *
 * @param os an output stream
 * @return propagated return code
 */
int Convection_Diffusion_Concentration::impr(Sortie& os) const
{
  return Equation_base::impr(os);
}

/*! @brief Returns 1 if the specified keyword designates an unknown field type of the equation.
 *
 *     Returns 1 if mot = "concentration"
 *     Returns 0 otherwise.
 *     If the method returns 1, ch_ref holds a reference to the field of the specified type.
 *
 * @param mot the field type keyword to look up
 * @param ch_ref reference to the field of the specified type
 * @return 1 if the field was found, 0 otherwise
 */
inline int string2int(char* digit, int& result)
{
  result = 0;

  //--- Convert each digit char and add into result.
  while (*digit >= '0' && *digit <='9')
    {
      result = (result * 10) + (*digit - '0');
      digit++;
    }

  //--- Check that there were no non-digits at end.
  if (*digit != 0)
    {
      return 0;
    }

  return 1;
}



/*! @brief Returns the name of the application domain of the equation.
 *
 * Here "Concentration".
 *
 * @return the name of the application domain of the equation
 */
const Motcle& Convection_Diffusion_Concentration::domaine_application() const
{
  static Motcle domaine = "Concentration";
  return domaine;
}

