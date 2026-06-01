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

#include <Raccord_distant_homogene.h>
#include <Frontiere_dis_base.h>
#include <Champ_front_calc.h>
#include <Probleme_Couple.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Domaine_VF.h>

Implemente_instanciable_sans_constructeur(Champ_front_calc,"Champ_front_calc",Ch_front_var_instationnaire_dep);
// XD champ_front_calc front_field_base champ_front_calc NO_BRACE This keyword is used on a boundary to get a field from
// XD_CONT another boundary. The local and remote boundaries should have the same mesh. If not, the
// XD_CONT Champ_front_recyclage keyword could be used instead. It is used in the condition block at the limits of
// XD_CONT equation which itself refers to a problem called pb1. We are working under the supposition that pb1 is
// XD_CONT coupled to another problem.
// XD attr problem_name ref_Pb_base problem_name REQ Name of the other problem to which pb1 is coupled.
// XD attr bord chaine bord REQ Name of the side which is the boundary between the 2 domains in the domain object
// XD_CONT description associated with the problem_name object.
// XD attr field_name chaine field_name REQ Name of the field containing the value that the user wishes to use at the
// XD_CONT boundary. The field_name object must be recognized by the problem_name object.


Champ_front_calc::Champ_front_calc() { set_distant(1); }

Sortie& Champ_front_calc::printOn(Sortie& os) const { return os; }

/*! @brief Read the name of a problem, the name of a boundary and the name of an unknown field from an input stream.
 *
 *     Then creates the corresponding boundary field.
 *     Format:
 *       Champ_front_calc problem_name boundary_name field_name
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 */
Entree& Champ_front_calc::readOn(Entree& is)
{
  is >> nom_autre_pb_ >> nom_autre_bord_ >> nom_inco_;
  via_readon_ = true;
  (nom_inco_ == "VITESSE") ? fixer_nb_comp(dimension) : fixer_nb_comp(1);
  return is;
}

/*! @brief Create the Champ_front_calc object representing the trace of an unknown field on a boundary from names:
 *
 *          - of the problem carrying the unknown
 *          - of the concerned boundary
 *          - of the unknown
 *
 * @param (Nom& nom_pb) the name of the problem to which the unknown whose trace we want belongs
 * @param (Nom& nom) the name of the boundary on which the trace of the unknown is taken
 * @param (Motcle& nom_inco) the name of the unknown whose trace we want
 * @throws no problem with the specified name
 * @throws the problem does not have a field with the specified name
 */
void Champ_front_calc::creer(const Nom& nom_pb, const Nom& nom_bord, const Motcle& nom_inco)
{
  nom_autre_pb_ = nom_pb;
  nom_autre_bord_ = nom_bord;
  OBS_PTR(Probleme_base) autre_pb;
  Objet_U& ob = Interprete::objet(nom_autre_pb_);
  if (sub_type(Probleme_base, ob))
    {
      autre_pb = ref_cast(Probleme_base, ob);
    }
  else
    {
      Cerr << "We did not find problem with name " << nom_pb << finl;
      exit();
    }
  OBS_PTR(Champ_base) rch;
  rch = autre_pb->get_champ(nom_inco);
  if (sub_type(Champ_Inc_base, rch.valeur()))
    {
      l_inconnue = ref_cast(Champ_Inc_base, rch.valeur());
      fixer_nb_comp(rch->nb_comp());
    }
  else
    {
      Cerr << autre_pb->le_nom() << " did not have unknown field with name " << nom_inco << finl;
      exit();
    }
}

int Champ_front_calc::initialiser(double temps, const Champ_Inc_base& inco)
{
  // First thing to do: here it is fine! we create the field
  if (via_readon_) creer(nom_autre_pb_, nom_autre_bord_, nom_inco_);

  Ch_front_var_instationnaire_dep::initialiser(temps, inco);

  // Check/initialize Raccord boundaries in parallel:
  if (nproc() > 1)
    {
      const Domaine_dis_base& domaine_dis_opposee = front_dis().domaine_dis();
      const Domaine_dis_base& domaine_dis_locale = frontiere_dis().domaine_dis();
      const Frontiere& frontiere_opposee = front_dis().frontiere();
      const Frontiere& frontiere_locale = frontiere_dis().frontiere();
      if (distant_ && !sub_type(Raccord_distant_homogene, frontiere_opposee))
        {
          const Nom& nom_domaine_oppose = domaine_dis_opposee.domaine().le_nom();
          Cerr << "Error, the boundary " << frontiere_opposee.le_nom() << " should be a Raccord." << finl;
          Cerr << "Add in your data file between the definition and the partition of the domain " << nom_domaine_oppose << " :" << finl;
          Cerr << "Modif_bord_to_raccord " << nom_domaine_oppose << " " << frontiere_opposee.le_nom() << finl;
          exit();
        }
      if (distant_ == 1)
        {
          Raccord_distant_homogene& raccord_distant = ref_cast_non_const(Raccord_distant_homogene, frontiere_opposee);
          if (!raccord_distant.est_initialise())
            raccord_distant.initialise(frontiere_locale, domaine_dis_locale, domaine_dis_opposee);
        }
    }
  return 1;
}

/*! @brief Associate the unknown field to the object
 *
 * @param (Champ_Inc_base& inco) the unknown field whose trace will be taken
 */
void Champ_front_calc::associer_ch_inc_base(const Champ_Inc_base& inco)
{
  l_inconnue = inco;
}


/*! @brief Not implemented
 *
 * @param (Champ_front_base& ch)
 * @return (Champ_front_base&)
 */
Champ_front_base& Champ_front_calc::affecter_(const Champ_front_base& ch)
{
  return *this;
}

/*! @brief Time update of the field. We simply take the trace of the unknown field at the
 *
 *     time step at which it is located.
 *     WEC: verify that we take the unknown at the correct time!
 *
 */

// Precondition: the boundary name must be different from "??"
// Parameter: double
//    Meaning:
//    Default values:
//    Constraints:
//    Access:
// Return:
//    Meaning:
//    Constraints:
// Exception:
// Side effects:
// Postcondition:
void Champ_front_calc::mettre_a_jour(double temps)
{
  assert (nom_autre_bord_ != "??") ;
  DoubleTab& tab=valeurs_au_temps(temps);
  const Frontiere_dis_base& frontiere_dis_opposee = domaine_dis().frontiere_dis(nom_bord_oppose());
  l_inconnue->trace(frontiere_dis_opposee,tab,temps,distant_ /* distant */);
}

/*! @brief Returns the associated unknown field
 *
 * @return (Champ_Inc_base&) the associated unknown field
 */
const Champ_Inc_base& Champ_front_calc::inconnue() const
{
  return l_inconnue.valeur();
}

/*! @brief Returns the name of the boundary on which the trace is computed.
 *
 * @return (Nom&) the name of the boundary on which the trace is computed
 */
const Nom& Champ_front_calc::nom_bord_oppose() const
{
  return nom_autre_bord_;
}

/*! @brief Returns the equation associated with the unknown whose trace is taken.
 *
 * @return (Equation_base&) the equation associated with the unknown whose trace is taken
 */
const Equation_base& Champ_front_calc::equation() const
{
  if (!l_inconnue)
    {
      Cerr << "\nError in Champ_front_calc::equation() : not able to return the equation !" << finl;
      Process::exit();
    }
  return inconnue().equation();
}

/*! @brief Returns the medium associated with the equation carrying the unknown field whose trace is taken.
 *
 * @return (Milieu_base&) the medium associated with the equation carrying the unknown field whose trace is taken
 */
const Milieu_base& Champ_front_calc::milieu() const
{
  return equation().milieu();
}

/*! @brief Returns the discretized domain associated with the equation carrying the unknown field whose trace is taken.
 *
 * @return (Domaine_dis_base&) the discretized domain associated with the equation carrying the unknown field whose trace is taken
 */
const Domaine_dis_base& Champ_front_calc::domaine_dis() const
{
  return inconnue().domaine_dis_base();
}

/*! @brief Returns the domain of discretized boundary conditions carried by the equation carrying the unknown field
 *
 *     whose trace is taken
 *
 * @return (Domaine_Cl_dis_base&) the domain of discretized boundary conditions carried by the equation carrying the unknown field whose trace is taken
 */
const Domaine_Cl_dis_base& Champ_front_calc::domaine_Cl_dis() const
{
  return equation().domaine_Cl_dis();
}

/*! @brief Returns the discretized boundary corresponding to the domain on which the trace is taken.
 *
 * @return (Frontiere_dis_base&) discretized boundary corresponding to the domain on which the trace is taken
 * @throws boundary with the specified name not found
 */
const Frontiere_dis_base& Champ_front_calc::front_dis() const
{
  return domaine_dis().frontiere_dis(nom_autre_bord_);
}


void Champ_front_calc::verifier(const Cond_lim_base& la_cl) const
{
}
