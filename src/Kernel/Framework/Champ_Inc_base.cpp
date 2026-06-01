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

#include <EcritureLectureSpecial.h>
#include <Scalaire_impose_paroi.h>
#include <Domaine_Cl_dis_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <Champ_Inc_P1_base.h>
#include <Neumann_val_ext.h>
#include <MD_Vector_tools.h>
#include <Champ_Inc_base.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <TRUST_2_PDI.h>
#include <Domaine_VF.h>
#include <YAML_data.h>
#include <Dirichlet.h>
#include <Domaine.h>

Implemente_base_sans_constructeur(Champ_Inc_base,"Champ_Inc_base",Champ_base);

Sortie& Champ_Inc_base::printOn(Sortie& os) const { return Champ_base::printOn(os); }
Entree& Champ_Inc_base::readOn(Entree& is) { return Champ_base::readOn(is); }

/*! @brief Sets the number of temporal values to keep.
 *
 * (a different number depending on the time scheme used)
 *     Calls Roue::fixer_nb_cases(int)
 *
 * @param (int i) the number of temporal values to keep
 * @return (int) the number of temporal values to keep
 */
int Champ_Inc_base::fixer_nb_valeurs_temporelles(int i)
{
  return les_valeurs->fixer_nb_cases(i);
}

/*! @brief Returns the number of temporal values currently kept.
 *
 * This value is stored by the Roue of Champ_Inc_base
 *
 * @return (int) the number of temporal values currently kept
 */
int Champ_Inc_base::nb_valeurs_temporelles() const
{
  return les_valeurs->nb_cases();
}

/*! @brief Reads the field values from an input stream.
 *
 * Reading format:
 *       int [THE NUMBER OF VALUES TO READ]
 *       [READ THE NUMBER OF VALUES WANTED]
 *
 * @param (Entree& is) the input stream
 * @return (int) returns 1 if the reading is correct
 * @throws the number of values to read is incorrect
 */
int Champ_Inc_base::lire_donnees(Entree& is)
{
  int n;
  is >> n;
  if (n != les_valeurs->valeurs().size())
    {
      Cerr << " the file does not contain the correct number of values to fill the field" << finl;
      Process::exit();
    }
  DoubleVect& tab = les_valeurs->valeurs();
  for (int i = 0; i < n; i++)
    is >> tab[i];
  return 1;
}

int Champ_Inc_base::fixer_nb_valeurs_nodales(int nb_noeuds)
{
  Cerr << "Internal error in Champ_Inc_base::fixer_nb_valeurs_nodales: method has not been implemented for class " << que_suis_je() << finl;
  Process::exit();
  return nb_noeuds;
}

void Champ_Inc_base::creer_tableau_distribue(const MD_Vector& md, RESIZE_OPTIONS opt)
{
  const int n = nb_valeurs_temporelles();
  for (int i = 0; i < n; i++)
    {
      DoubleTab& tab = futur(i);
      // Note B.M: This test is not symmetric with Champ_Fonc_base => inconsistency in nb_dim
      // for "multi-scalar" fields with one component.
      if (tab.size_array() == 0 && (!tab.get_md_vector()))
        {
          // Note B.M.: the fixer_nb_valeurs_nodales methods are called haphazardly.
          // Do nothing if the array already has the correct structure
          tab.resize(0, nb_compo_);
        }
      if (!(tab.get_md_vector() == md))
        {
          if (tab.get_md_vector())
            {
              Cerr << "Internal error in Champ_Inc_base::creer_tableau_distribue:\n" << " array has alreary a (wrong) parallel descriptor" << finl;
              Process::exit();
            }
          MD_Vector_tools::creer_tableau_distribue(md, tab, opt);
        }
    }
}

/*! @brief Returns the number of "real" geometric positions of the degrees of freedom, or -1 if not applicable (fields with multiple localisations)
 *
 */
int Champ_Inc_base::nb_valeurs_nodales() const
{
  int n;
  const DoubleTab& v = valeurs();
  if (v.size_reelle_ok())
    n = v.dimension(0);
  else
    n = -1;
  return n;
}

/*! @brief Returns the field values at time instant.
 *
 * @param (double temps) the time at which we want the values of the field
 * @return (DoubleTab&) the field values at time instant
 */
// WEC : Warning in the case of Pb_Couple we used the fact that this function returned the present when it didn't find
// a time greater than all available times!!!
// The behavior is now more explicit : a WARNING is displayed as soon as the present is returned instead of the requested time.
DoubleTab& Champ_Inc_base::valeurs(double tps)
{
  if (temps() == tps)
    return valeurs();
  else
    {
      Roue& la_roue = les_valeurs.valeur();
      if (temps() < tps)
        {
          for (int i = 0; i < nb_valeurs_temporelles(); i++)
            {
              if (la_roue.futur(i).temps() == tps)
                return la_roue.futur(i).valeurs();
              else if (la_roue.futur(i).temps() < temps())
                break;
            }
        }
      else if (temps() > tps)
        {
          for (int i = 0; i < nb_valeurs_temporelles(); i++)
            {
              if (la_roue.passe(i).temps() == tps)
                return la_roue.passe(i).valeurs();
              else if (la_roue.passe(i).temps() > temps())
                break;
            }
        }
    }
  Cerr << "ERROR : in Champ_Inc_base::valeurs(double), time " << tps << " not found, returns the present?" << finl;
  Cerr << "Contact TRUST support." << finl;
  Process::exit();
  return valeurs();
}

/*! @brief Returns the field values at time instant.
 *
 * @param (double temps) the time at which we want the values of the field
 * @return (DoubleTab&) the field values at time instant
 */
const DoubleTab& Champ_Inc_base::valeurs(double tps) const
// See above!
{
  if (temps() == tps)
    return valeurs();
  else
    {
      const Roue& la_roue = les_valeurs.valeur();

      if (temps() < tps)
        {
          // Future?
          for (int i = 0; i < nb_valeurs_temporelles(); i++)
            {
              if (la_roue.futur(i).temps() == tps)
                return la_roue.futur(i).valeurs();
              else if (la_roue.futur(i).temps() < temps())
                break;
            }
        }
      else if (temps() > tps)
        {
          // Past?
          for (int i = 1; i < nb_valeurs_temporelles(); i++)
            {
              if (la_roue.passe(i).temps() == tps)
                return la_roue.passe(i).valeurs();
              else if (la_roue.passe(i).temps() > temps())
                break;
            }
        }
    }
  Cerr << "ERROR : in Champ_Inc_base::valeurs(double), time " << tps << " not found, returns the present?" << finl;
  Cerr << "Contact TRUST support." << finl;
  Process::exit();
  return valeurs();
}

/*! @brief Advances the current pointer by i time steps, in the list of kept temporal values.
 *
 * @param (int i) the number of time steps to advance
 * @return (Champ_Inc_base&) returns *this, the field at the desired time step
 */
Champ_Inc_base& Champ_Inc_base::avancer(int i)
{
  while (i--)
    les_valeurs->avancer(les_valeurs);
  temps_ = les_valeurs->temps();
  return *this;
}

/*! @brief Rewinds the current pointer by i time steps, in the list of kept temporal values.
 *
 * @param (int i) the number of time steps to rewind
 * @return (Champ_Inc_base&) returns *this, the field at the desired time step
 */
Champ_Inc_base& Champ_Inc_base::reculer(int i)
{
  while (i--)
    les_valeurs->reculer(les_valeurs);
  temps_ = les_valeurs->temps();
  return *this;
}

/*! @brief Performs a time update of the unknown field.
 *
 * WEC : Now if we call it 2 times in a row with the same
 *     argument, the 2nd one does nothing.
 *
 * @param (double temps) the new time
 */
void Champ_Inc_base::mettre_a_jour(double un_temps)
{
  // Field with multiple temporal values:
  // Advance to the correct temporal value.
  if (les_valeurs->nb_cases() > 1)
    {
      for (int i = 0; i < les_valeurs->nb_cases(); i++)
        {
          if (les_valeurs[i].temps() == un_temps)
            {
              avancer(i);
              temps_ = un_temps;
              //Useless:
              //valeurs().echange_espace_virtuel();
              if (fonc_calc_)
                fonc_calc_(obj_calc_.valeur(), valeurs(), val_bord_, deriv_);
              /* first calculation of a Champ_Fonc_Calc -> copy the calculated values to all cases */
              if (fonc_calc_ && !fonc_calc_init_)
                for (int j = 1; j < les_valeurs->nb_cases(); j++, fonc_calc_init_ = 1)
                  les_valeurs[j].valeurs() = valeurs();
              return;
            }
        }
      Cerr << "In Champ_Inc_base::mettre_a_jour(double), " << finl;
      Cerr << "time " << un_temps << " not found in field " << le_nom() << finl;
      Cerr << "The times available are :" << finl;
      for (int i = 0; i < les_valeurs->nb_cases(); i++)
        Cerr << "  " << les_valeurs[i].temps() << finl;
      Process::exit();
    }
  // Field with a single temporal value:
  // Change the associated time.
  else
    {
      changer_temps(un_temps);
      if (fonc_calc_)
        fonc_calc_(obj_calc_.valeur(), valeurs(), val_bord_, deriv_);
      //Useless:
      //valeurs().echange_espace_virtuel();
    }
}

/*! @brief Sets the time of the i-th future field.
 *
 * @param (double t, int i) the new time
 * @return (double) the new time
 */
double Champ_Inc_base::changer_temps_futur(double t, int i)
{
  Roue& la_roue = les_valeurs.valeur();
  la_roue.futur(i).changer_temps(t);
  return t;
}

/*! @brief Sets the time of the i-th past field.
 *
 * @param (double t, int i) the new time
 * @return (double) the new time
 */
double Champ_Inc_base::changer_temps_passe(double t, int i)
{
  Roue& la_roue = les_valeurs.valeur();
  la_roue.passe(i).changer_temps(t);
  return t;
}

/*! @brief Returns the time of the i-th future field.
 *
 * @param (int i) the time
 * @return (double) the time
 */
double Champ_Inc_base::recuperer_temps_futur(int i) const
{
  const Roue& la_roue = les_valeurs.valeur();
  return la_roue.futur(i).temps();
}

/*! @brief Returns the time of the i-th past field.
 *
 * @param (int i) the time
 * @return (double) the time
 */
double Champ_Inc_base::recuperer_temps_passe(int i) const
{
  const Roue& la_roue = les_valeurs.valeur();
  return la_roue.passe(i).temps();
}
/*! @brief for PDI IO: retrieve the name of the HDF5 dataset in which the field will be saved or be restored from
 */
Nom Champ_Inc_base::get_PDI_dname() const
{
  Nom name = PDI_dname_;
  if(name == "??")
    {
      if(mon_equation_non_nul())
        name = equation().probleme().le_nom() + "_" + le_nom();
      else
        {
          // Sometimes (with Champ_fonc_reprise), the equation has not been associated yet so we can't prefix the name of the dataset with the name of the problem
          Cerr << "Champ_Inc_base::get_PDI_dname equation has not been associated yet. Please set the dataset name with Champ_Inc_base::set_PDI_dname()." << finl;
          Process::exit();
        }
    }
  return (Motcle)name;
}

/*! @brief for PDI IO: retrieve name, type and dimensions of the field to save/restore.
 */
std::vector<YAML_data> Champ_Inc_base::data_a_sauvegarder() const
{
  const Nom& name = get_PDI_dname();
  int nb_dim = valeurs().nb_dim();
  YAML_data d(name.getString(), "double", nb_dim);
  d.set_save_field_type(PDI_save_type_);
  std::vector<YAML_data> data;
  data.push_back(d);
  return data;
}

/*! @brief Saves the unknown field to an output stream.
 *
 *  Writes an identifier, the field values, and the date (the time at the moment of saving).
 *
 * @param (Sortie& fich) an output stream
 * @return (int) returns the size of array
 */
int Champ_Inc_base::sauvegarder(Sortie& fich) const
{
  // in special write mode only the master writes the header
  int a_faire, special;
  EcritureLectureSpecial::is_ecriture_special(special, a_faire);

  if (a_faire)
    {
      Nom mon_ident(nom_);
      mon_ident += que_suis_je();
      mon_ident += equation().probleme().domaine().le_nom();
      mon_ident += Nom(temps_, "%e");
      fich << mon_ident << finl;
      fich << que_suis_je() << finl;
      fich << temps_ << finl;
    }
  int bytes = 0;
  if (special)
    bytes = EcritureLectureSpecial::ecriture_special(*this, fich);
  else if (TRUST_2_PDI::is_PDI_checkpoint())
    {
      bytes = 8 * valeurs().size_array();

      // Sharing the dimensions of the unknown field with PDI
      TRUST_2_PDI pdi_interface;
      const Nom& name = get_PDI_dname();
      pdi_interface.share_TRUSTTab_dimensions(valeurs(), name, 1 /*write mode*/);

      if(PDI_save_type_)
        pdi_interface.share_type(name, que_suis_je());

      // Sharing the unknown field with PDI
      if( valeurs().dimension_tot(0) )
        pdi_interface.TRUST_start_sharing(name.getString(), valeurs().addr());
      else
        {
          // if the dimension is null in a direction - might happen in parallel - sharing an empty array
          ArrOfDouble garbage( valeurs().nb_dim() );
          pdi_interface.TRUST_start_sharing(name.getString(), garbage.addr());
        }
    }
  else
    {
      bytes = 8 * valeurs().size_array();
      valeurs().ecrit(fich);
    }

  if (a_faire)
    {
      // fich << flush ; Does not flush in binary mode!
      fich.flush();
    }
  Cerr << "Backup of the field " << nom_ << " performed on time : " << Nom(temps_, "%e") << finl;
  if (!est_egal(temps_, equation().probleme().schema_temps().temps_courant()))
    {
      Cerr.precision(12);
      Cerr << "Problem in Champ_Inc_base::sauvegarder, temps_=" << temps_ << " temps_courant()=" << equation().probleme().schema_temps().temps_courant() << finl;
      Process::exit();
    }
  // Return the number of bytes written
  return bytes;
}

/*! @brief Reads an unknown field from an input stream for a restart.
 *
 * @param (Entree& fich) an input stream
 * @return (int) always returns 1
 */
int Champ_Inc_base::reprendre(Entree& fich)
{
  double un_temps;
  int special = EcritureLectureSpecial::is_lecture_special();
  if (nom_ != Nom("anonyme")) // reading for restart
    {
      Cerr << "Resume of the field " << nom_ << finl;
      if(TRUST_2_PDI::is_PDI_restart())
        {
          TRUST_2_PDI pdi_interface;
          const Nom& name = get_PDI_dname();
          pdi_interface.share_TRUSTTab_dimensions(valeurs(), name, 0 /*read mode*/);
          if( valeurs().dimension_tot(0) )
            pdi_interface.read(name.getChar(), valeurs().addr());
          else
            {
              ArrOfDouble garbage( valeurs().nb_dim() );
              pdi_interface.read(name.getChar(), garbage.addr());
            }
        }
      else
        {
          int nb_val_nodales_old = nb_valeurs_nodales();
          fich >> un_temps;
          if (special)
            EcritureLectureSpecial::lecture_special(*this, fich);
          else
            valeurs().lit(fich);

          if (nb_val_nodales_old != nb_valeurs_nodales())
            {
              Cerr << finl << "Problem in the resumption " << finl;
              Cerr << "The field wich is read, does not have same number of nodal values" << finl;
              Cerr << "that the field created by the discretization " << finl;
              Process::exit();
            }
        }
      Cerr << " performed." << finl;
    }
  else // reading to skip the block
    {
      if(TRUST_2_PDI::is_PDI_restart())
        {
          Cerr << finl << "Problem in the resumption " << finl;
          Cerr << "PDI format does not require to navigate through file..." << finl;
          Process::exit();
        }
      BigDoubleTab tempo;
      fich >> un_temps;
      tempo.jump(fich);
    }
  return 1;
}

/*! @brief Computes the values of the unknown field at the specified positions.
 *
 * @param (DoubleTab& positions) the positions where the unknown field must be computed
 * @param (DoubleTab& valeurs) the array of unknown field values at the desired positions
 * @return (DoubleTab&) the array of unknown field values at the desired positions
 */
DoubleTab& Champ_Inc_base::valeur_aux(const DoubleTab& positions, DoubleTab& tab_valeurs) const
{
  const Domaine& domaine = domaine_dis_base().domaine();
  IntTrav les_polys;
  domaine.chercher_elements(positions, les_polys);

  return valeur_aux_elems(positions, les_polys, tab_valeurs);
}

/*! @brief Computes the values of the unknown field at the specified positions, for a given component of the field.
 *
 * @param (DoubleTab& positions) the positions where the unknown field must be computed
 * @param (DoubleTab& les_valeurs) the array of unknown field values at the desired positions
 * @param (int) the index of the field component to compute
 * @return (DoubleVect&) the array of values of the specified field component at the desired positions
 */
DoubleVect& Champ_Inc_base::valeur_aux_compo(const DoubleTab& positions, DoubleVect& tab_valeurs, int ncomp) const
{
  const Domaine& domaine = domaine_dis_base().domaine();
  IntTrav les_polys;
  domaine.chercher_elements(positions, les_polys);
  return valeur_aux_elems_compo(positions, les_polys, tab_valeurs, ncomp);
}

/*! @brief Computes the value of the unknown field at the specified position.
 *
 * @param (DoubleVect& position) the position at which the field is to be computed
 * @param (DoubleVect& les_valeurs) the value of the unknown field at the specified position
 * @return (DoubleVect&) the value of the unknown field at the specified position
 */
DoubleVect& Champ_Inc_base::valeur_a(const DoubleVect& position, DoubleVect& tab_valeurs) const
{
  const Domaine& domaine = domaine_dis_base().domaine();
  IntVect le_poly(1);
  domaine.chercher_elements(position, le_poly);
  return valeur_a_elem(position, tab_valeurs, le_poly(0));
}

/*! @brief Assignment of a generic OWN_PTR(Champ_base) (Champ_base) to an unknown field.
 *
 * @param (Champ_base& ch) the field on the right side of the assignment
 * @return (Champ_base&) the result of the assignment (*this)
 */
Champ_base& Champ_Inc_base::affecter_(const Champ_base& ch)
{
  DoubleTab noeuds;
  remplir_coord_noeuds(noeuds);

  if (valeurs().size_reelle_ok())
    {
      // Modif B.M. to avoid interpolation on virtual cells
      const int n = valeurs().dimension(0);
      DoubleTab pos, val;
      pos.ref_tab(noeuds, 0, n);
      val.ref_tab(valeurs(), 0, n);
      ch.valeur_aux(pos, val);
      //copy to all cases
      valeurs().echange_espace_virtuel();
      for (int i = 1; i < les_valeurs->nb_cases(); i++)
        les_valeurs[i].valeurs() = valeurs();
    }
  else
    {
      Cerr << "Champ_Inc_base::affecter_ not coded if size_reelle_ok()==0" << finl;
      Process::exit();
    }
  return *this;
}

//-Periodic BC case: ensures that the values on opposing periodic faces
// are identical. To do so, we take the half sum of the two values.
//The method must be overridden for fields discretized at faces.
void Champ_Inc_base::verifie_valeurs_cl()
{
}

/*! @brief Assignment of a component of a generic OWN_PTR(Champ_base) (Champ_base) to a component of an unknown field
 *
 * @param (Champ_base& ch) the right side of the assignment
 * @param (int compo) the index of the component to assign
 * @return (Champ_base&) the result of the assignment (with upcast)
 */
Champ_base& Champ_Inc_base::affecter_compo(const Champ_base& ch, int compo)
{
  DoubleTab noeuds;
  IntVect polys;
  if (!remplir_coord_noeuds_et_polys_compo(noeuds, polys, compo))
    {
      remplir_coord_noeuds_compo(noeuds, compo);
      ch.valeur_aux_compo(noeuds, valeurs(), compo);
    }
  else
    ch.valeur_aux_elems_compo(noeuds, polys, valeurs(), compo);
  return *this;
}

/*! @brief See Champ_base. Special case (unfortunately) of Champ_P0_VDF:
 *
 *     If the boundary is a connector, the result is computed on the associated connector. In this case, the DoubleTab x must be
 *     dimensioned on the associated connector.
 *
 */
DoubleTab& Champ_Inc_base::trace(const Frontiere_dis_base&, DoubleTab& x, double tps, int distant) const
{
  Cerr << que_suis_je() << "did not overloaded Champ_Inc_base::trace" << finl;
  return x;
}

/*! @brief DOES NOTHING. Method to override
 *
 * @param (DoubleTab&)
 * @param (IntVect&)
 * @return (int) always returns 0
 */
int Champ_Inc_base::remplir_coord_noeuds_et_polys(DoubleTab&, IntVect&) const
{
  return 0;
}

/*! @brief Simple call to Champ_Inc_base::remplir_coord_noeuds(DoubleTab&)
 *
 * @param (DoubleTab& coord) coordinates of the nodes to modify
 * @param (int) the index of the component to modify
 * @return (DoubleTab&)
 */
DoubleTab& Champ_Inc_base::remplir_coord_noeuds_compo(DoubleTab& coord, int) const
{
  return remplir_coord_noeuds(coord);
}

/*! @brief Simple call to: Champ_Inc_base::remplir_coord_noeuds_et_polys(DoubleTab&,IntVect& poly)
 *
 * @param (DoubleTab& coord)
 * @param (IntVect& poly)
 * @param (int)
 * @return (int) propagated return code
 */
int Champ_Inc_base::remplir_coord_noeuds_et_polys_compo(DoubleTab& coord, IntVect& poly, int) const
{
  return remplir_coord_noeuds_et_polys(coord, poly);
}

const Domaine& Champ_Inc_base::domaine() const
{
  return domaine_dis_base().domaine();
}

int Champ_Inc_base::imprime(Sortie& os, int ncomp) const
{
  Cerr << que_suis_je() << "::imprime not coded." << finl;
  Process::exit();
  return 1;
}

double Champ_Inc_base::integrale_espace(int ncomp) const
{
  Cerr << que_suis_je() << "::integrale_espace not coded." << finl;
  Process::exit();
  return 0.;
}

/*! @brief Sets the time of the field.
 *
 * @param (double t) the new time
 * @return (double) the new time
 */
double Champ_Inc_base::changer_temps(const double t)
{
  les_valeurs->changer_temps(t);
  return temps_ = t;
}

/*!
 * See comments in Probleme_base_interface_proto::resetTime_impl().
 * Here we force a new time value.
 */
void Champ_Inc_base::resetTime(double time)
{
  changer_temps(time);
}

/*! @brief Associates the field with the equation of which it represents an unknown.
 *
 * Simple call to MorEqn::associer_eqn(const Equation_base&)
 *
 * @param (Equation_base& eqn) the equation to which the field must be associated
 */
void Champ_Inc_base::associer_eqn(const Equation_base& eqn)
{
  MorEqn::associer_eqn(eqn);
}

void Champ_Inc_base::associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl)
{
  mon_dom_cl_dis = zcl;
}

void Champ_Inc_base::associer_domaine_dis_base(const Domaine_dis_base& z_dis)
{
  le_dom_VF = ref_cast(Domaine_VF, z_dis);
}

const Domaine_Cl_dis_base& Champ_Inc_base::domaine_Cl_dis() const
{
  if (!mon_dom_cl_dis)
    return equation().domaine_Cl_dis();
  else
    return mon_dom_cl_dis.valeur();
}

Domaine_Cl_dis_base& Champ_Inc_base::domaine_Cl_dis()
{
  if (!mon_dom_cl_dis)
    return equation().domaine_Cl_dis();
  else
    return mon_dom_cl_dis.valeur();
}

void Champ_Inc_base::init_champ_calcule(const Objet_U& obj, fonc_calc_t fonc)
{
  obj_calc_ = obj, fonc_calc_ = fonc, fonc_calc_init_ = 0;
  resize_val_bord();
}

void Champ_Inc_base::resize_val_bord()
{
  val_bord_.resize(ref_cast(Domaine_VF, domaine_dis_base()).xv_bord().dimension_tot(0), valeurs().line_size());
}

DoubleTab Champ_Inc_base::valeur_aux_bords() const
{
  if (fonc_calc_ || bord_fluide_multiphase_)
    {
      DoubleTab result;
      result.ref(val_bord_);
      return result;
    }
  //otherwise, compute from the BCs
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_dis_base());
  const IntTab& f_e = domaine.face_voisins(), &f_s = domaine.face_sommets();
  DoubleTrav result(domaine.xv_bord().dimension_tot(0), valeurs().line_size());

  const Conds_lim& cls = domaine_Cl_dis().les_conditions_limites();
  int j, k, f, fb, s, n, N = result.line_size(), is_p = (le_nom().debute_par("pression") || le_nom().debute_par("pressure")), n_som;
  for (const auto& itr : cls)
    {
      const Front_VF& fr = ref_cast(Front_VF, itr->frontiere_dis());
      //imposed boundary value, except if it is a wall (in which case the BC may have fewer components than the field -> Energie_Multiphase)
      if (is_p ? sub_type(Neumann, itr.valeur()) : (sub_type(Dirichlet, itr.valeur()) && !sub_type(Scalaire_impose_paroi, itr.valeur())))
        for (j = 0; j < fr.nb_faces_tot(); j++)
          for (f = fr.num_face(j), fb = domaine.fbord(f), n = 0; n < N; n++)
            result(fb, n) = is_p ? ref_cast(Neumann, itr.valeur()).flux_impose(j, n) : ref_cast(Dirichlet, itr.valeur()).val_imp(j, n);
      else if (sub_type(Neumann_val_ext, itr.valeur())) //externally imposed value
        for (j = 0; j < fr.nb_faces_tot(); j++)
          for (f = fr.num_face(j), fb = domaine.fbord(f), n = 0; n < N; n++)
            result(fb, n) = ref_cast(Neumann_val_ext, itr.valeur()).val_ext(j, n);
      else if (sub_type(Champ_Inc_P0_base, *this))
        for (j = 0; j < fr.nb_faces_tot(); j++) //P0 field: we can take the value at the element
          for (f = fr.num_face(j), fb = domaine.fbord(f), n = 0; n < N; n++)
            result(fb, n) = valeurs()(f_e(f, f_e(f, 0) == -1), n);
      else if (sub_type(Champ_Inc_P1_base, *this))
        for (j = 0; j < fr.nb_faces_tot(); j++) //P1 field: average of vertex values
          {
            f = fr.num_face(j), fb = domaine.fbord(f);
            for (n_som = 0; n_som < f_s.dimension(1) && f_s(f, n_som) >= 0;)
              n_som++;
            for (n = 0; n < N; n++)
              result(fb, n) = 0;
            for (k = 0; k < n_som; k++)
              for (s = f_s(f, k), n = 0; n < N; n++)
                result(fb, n) += valeurs()(s, n) / n_som;
          }
      else if (que_suis_je() == "Champ_P1NC")
        for (j = 0; j < fr.nb_faces_tot(); j++)
          for (f = fr.num_face(j), fb = domaine.fbord(f), n = 0; n < N; n++)
            result(fb, n) = valeurs()(f, n);
      else
        Process::exit("Champ_Inc_base::valeur_aux_bords() : must code something!");
    }
  return result;
}
