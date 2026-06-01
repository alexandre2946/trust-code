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

#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Operateur_base.h>
#include <Probleme_base.h>
#include <TRUSTTrav.h>
#include <Operateur.h>
#include <Perf_counters.h>

Sortie& Operateur::ecrire(Sortie& os) const
{
  os << " { " << typ << " } " << finl;
  return os;
}


/*! @brief Reads an operator from an input stream.
 *
 * Types the operator and associates its equation with it.
 *     Format:
 *       {
 *        [A KEYWORD REPRESENTING A TYPE]
 *       }
 *
 * @param (Entree& is) the input stream from which to read the operator
 * @return (Entree&) the modified input stream
 * @throws opening brace expected
 * @throws closing brace expected
 */
Entree& Operateur::lire(Entree& is)
{
  typ="non defini";
  Motcle motlu;
  is >> motlu ;
  assert(motlu=="{");
  while(motlu!="}")
    {
      is>>motlu;
      if (motlu=="}")
        {
          if(typ=="non defini")
            {
              typ = "";
              typer();
              l_op_base().associer_eqn(equation());
              l_op_base().set_nb_ss_pas_de_temps(1);
              l_op_base().set_decal_temps(0);
            }
        }
      else if (motlu=="nb_sous_pas_de_temps")
        {
          if(typ=="non defini")
            {
              typ = "";
              typer();
              l_op_base().associer_eqn(equation());
            }
          int nb_ss_pas_de_temps;
          is >> nb_ss_pas_de_temps;
          l_op_base().set_nb_ss_pas_de_temps(nb_ss_pas_de_temps);
          l_op_base().set_decal_temps(0);
        }
      else if (motlu=="implicite")
        {
          // XD op_implicite objet_lecture nul NO_BRACE not_set
          // XD attr implicite chaine(into=["implicite"]) implicite REQ not_set
          // XD attr mot chaine(into=["solveur"]) mot REQ not_set
          // XD attr solveur solveur_sys_base solveur REQ not_set
          if(typ=="non defini")
            {
              typ = "";
              typer();
              l_op_base().associer_eqn(equation());
            }
          l_op_base().set_decal_temps(1);
          l_op_base().set_nb_ss_pas_de_temps(1);
          is >> motlu;
          if(motlu != "solveur")
            {
              Cerr << "We expected the keyword \"solveur\" instead of : "
                   << motlu << finl;;
              Process::exit();
            }
          l_op_base().lire_solveur(is);
        }
      else
        {
          if(typ!="non defini")
            {
              Cerr << "We must choose the type of operator in beginning "
                   <<"of reading block (before implicit in particular) "
                   << finl;
              Process::exit();
            }
          typ=motlu;
          typer();
          l_op_base().associer_eqn(equation());
          l_op_base().set_decal_temps(0);
          l_op_base().set_nb_ss_pas_de_temps(1);
          is >> l_op_base();
        }
    }
  return is;
}

/*! @brief Returns the field representing the unknown of the equation to which the operator belongs.
 *
 * @return (Champ_Inc_base&) the unknown field of the associated equation
 */
const Champ_Inc_base& Operateur::mon_inconnue() const
{
  return le_champ_inco.valeur();
}

/*! @brief Returns the discretization of the equation to which the operator belongs.
 *
 * @return (Discretisation_base&) the discretization of the associated equation
 */
const Discretisation_base& Operateur::discretisation() const
{
  return mon_equation->discretisation();
}

/*! @brief Updates the references of the objects associated with the operator.
 *
 * Operateur::le_champ_inco, Operateur::champ_inco
 *     Calls Operateur_base::completer()
 *
 */
void Operateur::completer()
{
  if (!le_champ_inco)
    le_champ_inco=mon_equation->inconnue();

  l_op_base().completer();
}

void Operateur::associer_champ(const Champ_Inc_base& ch, const std::string& nom_ch)
{
  le_champ_inco = ch;
  nom_inco_ = nom_ch;
  l_op_base().associer_champ(ch, nom_ch);
}

/*! @brief Performs a time update of the operator.
 *
 * Calls Operateur_base::mettre_a_jour(double)
 *
 * @param (double temps) the time step for update
 */
void Operateur::mettre_a_jour(double temps)
{
  l_op_base().mettre_a_jour(temps);
}
/*! @brief Computes the next time step.
 *
 */
double Operateur::calculer_pas_de_temps() const
{
  // If the operator's equation is not solved, we do not compute its stability time step
  if (equation().equation_non_resolue())
    return DMAXFLOAT;
  statistics().begin_count(STD_COUNTERS::compute_dt,statistics().get_last_opened_counter_level()+1);
  double dt_stab = l_op_base().calculer_dt_stab();
  statistics().end_count(STD_COUNTERS::compute_dt);
  // Check that the operator did perform an mp_min:
  assert(dt_stab==Process::mp_min(dt_stab));
  return dt_stab;
}
/*! @brief Calculate the next local time steps
 *
 */
void Operateur::calculer_pas_de_temps_locaux(DoubleTab& dt_locaux) const
{
  l_op_base().calculer_dt_local(dt_locaux);
}
/*! @brief Asks the equation whether printing is needed. Returns 1 for YES, 0 otherwise.
 *
 * @return (int) 1 if printing is needed, 0 otherwise
 */
int  Operateur::limpr() const
{
  return mon_equation->limpr();
}

/*! @brief Prints the operator to an output stream, if necessary.
 *
 * (see Schema_Temp_base::limpr())
 *
 * @param (Sortie& os) the output stream for printing
 */
void  Operateur::imprimer(Sortie& os) const
{
  if(limpr()) l_op_base().impr(os);
}


/*! @brief Prints the operator to an output stream unconditionally.
 *
 * @param (Sortie& os) the output stream for printing
 * @return (int) return code of Operateur_base::impr(Sortie&)
 */
int Operateur::impr(Sortie& os) const
{
  if (l_op_base().has_impr_file()) return l_op_base().impr(os);
  return 1;
}

/*! @brief Computes and adds the contribution of the operator to the right-hand side of the equation.
 *
 *     Calls Operateur::ajouter(const DoubleTab&, DoubleTab& )
 *
 * @param (Champ_Inc_base& ch) the unknown field on which the operator acts
 * @param[in,out] (DoubleTab& resu) the array storing the right-hand side values to which the operator contribution is added
 * @return (DoubleTab&) the right-hand side to which the operator contribution has been added
 */
DoubleTab& Operateur::ajouter(const Champ_Inc_base& ch, DoubleTab& resu) const
{
  int i ;
  int nstep=l_op_base().get_nb_ss_pas_de_temps();
  double dt=equation().schema_temps().pas_de_temps();
  dt/=nstep;
  if(nstep==1)
    return ajouter(ch.valeurs(), resu);
  DoubleTrav derivee(resu);
  DoubleTrav inco(ch.valeurs());
  inco=ch.valeurs();
  const Solveur_Masse_base& solveur_masse=equation().solv_masse();
  double dt_inv=1./(double(nstep));
  for (i=0; i<nstep; i++)
    {
      calculer(inco, derivee);
      derivee.echange_espace_virtuel();
      resu.ajoute(dt_inv, derivee) ;
      solveur_masse.appliquer(derivee);
      inco.ajoute_sans_ech_esp_virt(dt, derivee, VECT_ALL_ITEMS) ;
    }
  return resu;
}

/*! @brief Returns the (name of the) type of operator to create.
 *
 * @return (Nom&) the name of the type of operator to create
 */
const Nom& Operateur::type() const
{
  return typ;
}


/*! @brief Computes the contribution of the operator and returns the array of values.
 *
 * @param (Champ_Inc_base& ch) the unknown field on which the operator acts
 * @param (DoubleTab& resu) the array storing the values resulting from applying the operator to the unknown field.
 * @return (DoubleTab&) the result of applying the operator to the unknown field
 */
DoubleTab& Operateur::calculer(const Champ_Inc_base& ch,DoubleTab& resu) const
{
  return calculer(ch.valeurs(), resu);
}

/*! @brief Adds the contribution of the operator to the array passed as parameter.
 *
 *     Calls Operateur::ajouter(const Champ_Inc_base&, DoubleTab& )
 *
 * @param (DoubleTab& resu) the array storing the right-hand side values to which the operator contribution is added
 * @return (DoubleTab&) the right-hand side to which the operator contribution has been added
 */
DoubleTab& Operateur::ajouter(DoubleTab& resu) const
{
  return ajouter(le_champ_inco->valeurs(), resu);
}

/*! @brief Applies the operator to the unknown field and returns the result.
 *
 * Calls Operateur::calculer(const Champ_Inc_base&, DoubleTab& );
 *
 * @param (DoubleTab& resu) the array storing the values resulting from applying the operator to the unknown field.
 * @return (DoubleTab&) the result of applying the operator to the unknown field
 */
DoubleTab& Operateur::calculer(DoubleTab& resu) const
{
  resu=0.;
  return ajouter(le_champ_inco->valeurs(), resu);
}

void Operateur::set_fichier(const Nom& nom)
{
  l_op_base().set_fichier(nom);
}

void Operateur::set_description(const Nom& nom)
{
  l_op_base().set_description(nom);
}

void Operateur::ajouter_contribution_explicite_au_second_membre (const Champ_Inc_base& linconnue, DoubleTab& derivee) const
{
  l_op_base().ajouter_contribution_explicite_au_second_membre (linconnue, derivee);
}
