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

#include <Frontiere_dis_base.h>
#include <Discret_Thermique.h>
#include <Probleme_base.h>
#include <Conduction.h>
#include <Solide.h>
#include <Param.h>

Implemente_instanciable(Conduction,"Conduction",Equation_base);
// XD Conduction eqn_base Conduction INHERITS_BRACE Heat equation.

Sortie& Conduction::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

Entree& Conduction::readOn(Entree& is )
{
  Equation_base::readOn(is);
  terme_diffusif.set_fichier("Diffusion_chaleur");
  terme_diffusif.set_description((Nom)"Conduction heat transfer rate=Integral(lambda*grad(T)*ndS) [W] if SI units used");
  return is;
}

void Conduction::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("diffusion",(this));
  param.ajouter_condition("is_read_diffusion","The diffusion operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_non_std("Traitement_particulier",(this));
}

//  Overrides the base Conduction::lire_motcle_non_standard()!
int Conduction::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  //Motcle motlu;
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;

      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      is >> terme_diffusif;
      // the dt_stab field is the same as the operator field. We always use diffusivite()
      // as in the Conduction case
      terme_diffusif.associer_diffusivite_pour_pas_de_temps(milieu().diffusivite());
      solveur_masse->set_name_of_coefficient_temporel("rho_cp_comme_T");

      return 1;
    }
  else if (mot=="Traitement_particulier")
    {
      Cerr << "Reading and typing of the Traitement_particulier object : " << finl;
      Nom type="Traitement_particulier_";
      Motcle motbidon;
      is >> motbidon;
      if (motbidon == "{")
        {
          Motcle le_cas;
          is >> le_cas;
          if (le_cas == "}")
            le_cas ="";
          else
            {
              type+= le_cas;
              type+= "_";
            }
          Nom discr=discretisation().que_suis_je();
          if (discr == "VEFPreP1B")
            discr = "VEF";
          type+=discr;
          Cerr << type << finl;
          le_traitement_particulier.typer(type);
          le_traitement_particulier->associer_eqn(*this);
          le_traitement_particulier->lire(is);
          le_traitement_particulier->preparer_calcul_particulier();
          return 1;
        }
      else
        {
          Cerr << "Error while reading Traitement_particulier" << finl;
          Cerr << "A { is expected." << finl;
          exit();
          return -1;
        }
    }
  else
    return Equation_base::lire_motcle_non_standard(mot,is);
}

// returns the *conductivity* rather than the diffusivity as in Conduction
const Champ_Don_base& Conduction::diffusivite_pour_transport() const
{
  return milieu().conductivite();
}

const Champ_base& Conduction::diffusivite_pour_pas_de_temps() const
{
  return terme_diffusif.diffusivite();
}

/*! @brief Associates a physical medium with the equation; the medium is cast to Solide.
 *
 * @param le_milieu Reference to the physical medium (cast to Solide internally).
 */
void Conduction::associer_milieu_base(const Milieu_base& le_milieu)
{
  associer_solide(ref_cast(Solide,le_milieu));
}

/*! @brief Associates the solid medium with the equation.
 *
 * @param un_solide The solid medium to associate with the equation.
 */
void Conduction::associer_solide(const Solide& un_solide)
{
  le_solide = un_solide;
}

/*! @brief Returns the number of operators in the equation. For the standard conduction equation this is always 1.
 *
 * @return Number of operators in the equation.
 */
int Conduction::nombre_d_operateurs() const
{
  return 1;
}

/*! @brief Returns the operator at the given index in the equation (const version).
 *
 * Returns terme_diffusif when i=0; exits when i>0.
 *
 * @param i Index of the operator to return.
 * @return The operator at the given index (only terme_diffusif is available).
 * @throws The standard conduction equation contains only one operator.
 */
const Operateur& Conduction::operateur(int i) const
{
  if (i == 0)
    return terme_diffusif;
  else
    {
      Cerr << "Conduction::operateur("<<i<<") !! " << finl;
      Cerr << "Equation Conduction has only one operator." << finl;
      exit();
    }
  // Required to satisfy compilers!!
  return terme_diffusif;
}

/*! @brief Returns the operator at the given index in the equation.
 *
 * Returns terme_diffusif when i=0; exits when i>0.
 *
 * @param i Index of the operator to return.
 * @return The operator at the given index (only terme_diffusif is available).
 * @throws The standard conduction equation contains only one operator.
 */
Operateur& Conduction::operateur(int i)
{
  if (i == 0)
    return terme_diffusif;
  else
    {
      Cerr << "Conduction::operateur("<<i<<") !! " << finl;
      Cerr << "Equation Conduction has only one operator." << finl;
      exit();
    }
  // Required to satisfy compilers!!
  return terme_diffusif;
}

/*! @brief Discretises the equation.
 *
 */
void Conduction::discretiser()
{
  const Discret_Thermique& dis=ref_cast(Discret_Thermique, discretisation());
  Cerr << "Conduction equation discretization" << finl;
  dis.temperature(schema_temps(), domaine_dis(), la_temperature);
  champs_compris_.ajoute_champ(la_temperature);
  Equation_base::discretiser();
}


/*! @brief Returns the physical medium associated with the equation (const version).
 *
 * Here Solide is upcast to Milieu_base.
 *
 * @return The physical medium associated with the equation (Solide upcast to Milieu_base).
 */
const Milieu_base& Conduction::milieu() const
{
  return solide();
}

/*! @brief Returns the physical medium associated with the equation.
 *
 * Here Solide is upcast to Milieu_base.
 *
 * @return The physical medium associated with the equation (Solide upcast to Milieu_base).
 */
Milieu_base& Conduction::milieu()
{
  return solide();
}

/*! @brief Returns the solid medium associated with the equation (const version).
 *
 * @return The solid medium associated with the equation.
 * @throws No solid medium has been associated with the equation.
 */
const Solide& Conduction::solide() const
{
  if(!le_solide)
    {
      Cerr << "You forgot to associate the solid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_solide.valeur();
}

void Conduction::creer_champ(const Motcle& motlu)
{
  if (motlu == "temperature_paroi" || motlu == "wall_temperature")
    {
      if (!temperature_paroi_)
        {
          const Discret_Thermique& dis = ref_cast(Discret_Thermique, discretisation());
          dis.t_paroi(domaine_dis(), domaine_Cl_dis(), la_temperature, temperature_paroi_);
          champs_compris_.ajoute_champ(temperature_paroi_);
        }
    }

  Equation_base::creer_champ(motlu);
}

bool Conduction::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  if (nom == "temperature_paroi" || nom == "wall_temperature")
    {
      ref_champ = Conduction::get_champ(nom);
      return true;
    }

  if (Equation_base::has_champ(nom, ref_champ))
    return true;

  if (le_traitement_particulier)
    if (le_traitement_particulier->has_champ(nom, ref_champ))
      return true;

  return false; /* nothing found */
}

bool Conduction::has_champ(const Motcle& nom) const
{
  if (nom == "temperature_paroi" || nom == "wall_temperature")
    return true;

  if (Equation_base::has_champ(nom))
    return true;

  if (le_traitement_particulier)
    if (le_traitement_particulier->has_champ(nom))
      return true;

  return false; /* nothing found */
}

const Champ_base& Conduction::get_champ(const Motcle& nom) const
{
  OBS_PTR(Champ_base) ref_champ;

  if (nom == "temperature_paroi" || nom == "wall_temperature")
    {
      double temps_init = schema_temps().temps_init();
      Champ_Fonc_base& ch_tp = ref_cast_non_const(Champ_Fonc_base, temperature_paroi_.valeur());
      if (((ch_tp.temps() != la_temperature->temps()) || (ch_tp.temps() == temps_init)) && ((la_temperature->mon_equation_non_nul())))
        ch_tp.mettre_a_jour(la_temperature->temps());
      return champs_compris_.get_champ(nom);
    }

  if (Equation_base::has_champ(nom, ref_champ))
    return ref_champ;

  if (le_traitement_particulier)
    if (le_traitement_particulier->has_champ(nom, ref_champ))
      return ref_champ;

  throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));
}

void Conduction::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Equation_base::get_noms_champs_postraitables(nom, opt);

  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("TEMPERATURE_PAROI");
  noms_compris.add("WALL_TEMPERATURE");

  if (opt == DESCRIPTION)
    Cerr << "Conduction : " << noms_compris << finl;
  else
    nom.add(noms_compris);

  if (le_traitement_particulier)
    le_traitement_particulier->get_noms_champs_postraitables(nom, opt);
}

/*! @brief Returns the solid medium associated with the equation.
 *
 * @return The solid medium associated with the equation.
 * @throws No solid medium has been associated with the equation.
 */
Solide& Conduction::solide()
{
  if(!le_solide)
    {
      Cerr << "A solide medium has not been associated to a Conduction equation"<<finl;
      exit();
    }
  return le_solide.valeur();
}

/*! @brief Prints the diffusive term to an output stream.
 *
 * @param os An output stream.
 * @return Always returns 1.
 */
int Conduction::impr(Sortie& os) const
{
  return Equation_base::impr(os);
}

/*! @brief Returns the application domain name of the equation.
 *
 * Here "Thermique".
 *
 * @return The application domain name of the equation.
 */
const Motcle& Conduction::domaine_application() const
{
  static Motcle domaine = "Thermique";
  return domaine;
}

void Conduction::mettre_a_jour(double temps)
{
  Equation_base::mettre_a_jour(temps);

  if (le_traitement_particulier)
    le_traitement_particulier->post_traitement_particulier();
}
