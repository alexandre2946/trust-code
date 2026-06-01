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

#include <Source_rayo_semi_transp_base.h>
#include <Cond_lim_rayo_semi_transp.h>
#include <Pb_rayo_semi_transp.h>
#include <Discretisation_base.h>
#include <Flux_radiatif_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <Symetrie.h>
#include <Domaine.h>

Implemente_instanciable(Pb_rayo_semi_transp, "Pb_rayo_semi_transp", Probleme_base);

Entree& Pb_rayo_semi_transp::readOn(Entree& is) { return Probleme_base::readOn(is); }

Sortie& Pb_rayo_semi_transp::printOn(Sortie& os) const { return os; }

bool Pb_rayo_semi_transp::initTimeStep(double dt)
{
  return eq_rayo().initTimeStep(dt);
}

bool Pb_rayo_semi_transp::iterateTimeStep(bool& converged)
{
  converged = true;
  return eq_rayo().resoudre();
}

void Pb_rayo_semi_transp::validateTimeStep()
{
  double temps = pb_fluide_->presentTime();
  eq_rayo().mettre_a_jour(temps);
  calculer_flux_radiatif();
  les_postraitements_.mettre_a_jour(temps);
  statistics().begin_count(STD_COUNTERS::update_variables, statistics().get_last_opened_counter_level() + 1);
  schema_temps().mettre_a_jour();
  statistics().end_count(STD_COUNTERS::update_variables);
}

void Pb_rayo_semi_transp::associer_sch_tps_base(const Schema_Temps_base& un_schema_en_temps)
{
  le_schema_en_temps_ = un_schema_en_temps;
  le_schema_en_temps_->associer_pb(*this);
}

void Pb_rayo_semi_transp::get_noms_champs_postraitables(Noms& noms,Option opt) const
{
  for (int i=0; i<nombre_d_equations(); i++)
    equation(i).get_noms_champs_postraitables(noms,opt);
}

void Pb_rayo_semi_transp::discretise_longueur_rayo()
{
  // Associate the fluid + various operations
  if (sub_type(Fluide_base, pb_fluide_->milieu()))
    {
      Fluide_base& fluide = ref_cast(Fluide_base, pb_fluide_->milieu());
      eq_rayo_.associer_fluide(fluide);

      if (fluide.is_rayo_semi_transp())
        {
          Champ_Don_base& coeff_abs = fluide.kappa();

          if (sub_type(Champ_Uniforme, coeff_abs))
            {
              // Type the OWN_PTR(Champ_Don_base) longueur_rayo_ as a Champ_Uniforme
              fluide.typer_longeur_rayo("Champ_Uniforme");
              Champ_Don_base& l_rayo = fluide.longueur_rayo();
              Champ_Uniforme& ch_l_rayo = ref_cast(Champ_Uniforme, l_rayo);
              ch_l_rayo.nommer("longueur_de_rayonnement");
              ch_l_rayo.fixer_nb_comp(1);
              // The number of nodal values is fixed to 1 here because this is a
              // Champ_Uniforme; in other cases it must be set equal to the number of elements or faces depending on the field location
              ch_l_rayo.fixer_nb_valeurs_nodales(1);
              ch_l_rayo.fixer_unite("m");
              ch_l_rayo.changer_temps(0);
            }
          else
            {
              Cerr << "The absorption coefficient is not a uniform OWN_PTR(Champ_base) but a " << coeff_abs.que_suis_je() << ". modify the method " << finl;
              Cerr << "Pb_Couple_rayo_semi_transp::discretiser to be able to handle this type of Champ_Don" << finl;
            }
          fluide.initialiser(pb_fluide_->schema_temps().temps_courant());
        }
      else
        {
          Cerr << "Error 0 in Pb_Couple_rayo_semi_transp::discretiser you have probably not specified all the" << finl;
          Cerr << "physical parameters of your incompressible fluid to handle a semi-transparent radiation problem" << finl;
          Process::exit();
        }
    }
  else
    {
      Cerr << "Error in Pb_rayo_semi_transp::readOn: the semi-transparent radiation problem can only be used" << finl;
      Cerr << "with a Fluide_base and not " << pb_fluide_->milieu().que_suis_je() << finl;
      Process::exit();
    }
}

void Pb_rayo_semi_transp::preparer_calcul()
{
  int contient_source_rayo_semi_transp = 0;

  for (int j = 0; j < pb_fluide_->nombre_d_equations(); j++)
    {

      // Associate the problem to the radiation boundary conditions.
      Domaine_Cl_dis_base& la_zcl = pb_fluide_->equation(j).domaine_Cl_dis();
      for (int num_cl = 0; num_cl < la_zcl.nb_cond_lim(); num_cl++)
        {
          Cond_lim_base& la_cl = la_zcl.les_conditions_limites(num_cl).valeur();
          Cond_lim_rayo_semi_transp *la_cl_rayo_semi_transp;

          if (la_cl.is_bc_rayo_semi_transp(la_cl_rayo_semi_transp))
            {
              la_cl_rayo_semi_transp->associer_pb_rayo_semi_transp(*this);
              la_cl_rayo_semi_transp->recherche_emissivite_et_A();

              // In the case of a contact exchange, the opposite BC must also be completed
              la_cl_rayo_semi_transp->completer_Cl_opposee_si_contact();
            }
        }

      // Associate the problem to the radiation source term of the temperature equation
      Sources& les_sources = pb_fluide_->equation(j).sources();
      for (int num_source = 0; num_source < les_sources.size(); num_source++)
        if ((sub_type(Source_rayo_semi_transp_base, les_sources[num_source].valeur())))
          contient_source_rayo_semi_transp = 1;
    }

  if (contient_source_rayo_semi_transp == 0)
    {
      Cerr << "Warning, you have not defined a semi-transparent radiation source term" << finl;
      Cerr << "remember to add the source term Source_rayo_semi_transp in the list of source terms " << finl;
      Cerr << "of the energy equation" << finl;
      Process::exit();
    }
  eq_rayo().completer();
}

void Pb_rayo_semi_transp::typer_lire_milieu(Entree& is)
{
  // Only discretise the equations
  discretiser_equations();
}

// Update the radiative flux for all boundaries of the problem
void Pb_rayo_semi_transp::calculer_flux_radiatif()
{
  Conds_lim& les_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites();

  for (int num_cl_rayo = 0; num_cl_rayo < les_cl_rayo.size(); num_cl_rayo++)
    {
      Cond_lim& la_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites(num_cl_rayo);
      if (sub_type(Flux_radiatif_base, la_cl_rayo.valeur()))
        {
          Flux_radiatif_base& la_cl_rayon = ref_cast(Flux_radiatif_base, la_cl_rayo.valeur());
          Equation_base& eq_temp = pb_fluide_->equation(1);
          la_cl_rayon.calculer_flux_radiatif(eq_temp);
        }
      else if (sub_type(Symetrie, la_cl_rayo.valeur()))
        {
          /* Do nothing */
        }
      else
        {
          Cerr << "Error: the boundary conditions of the radiation equation" << finl;
          Cerr << "must necessarily be of radiation type" << finl;
          Process::exit();
        }
    }
}

const Champ_front_base& Pb_rayo_semi_transp::flux_radiatif(const Nom& nom_bord) const
{
  // Loop over the boundaries to find the one whose name is nom_bord
  const Conds_lim& les_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites();

  for (int num_cl_rayo = 0; num_cl_rayo < les_cl_rayo.size(); num_cl_rayo++)
    {
      const Cond_lim& la_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites(num_cl_rayo);

      if (la_cl_rayo->frontiere_dis().le_nom() == nom_bord)
        {
          if (sub_type(Flux_radiatif_base, la_cl_rayo.valeur()))
            {
              Flux_radiatif_base& la_cl_rayon = ref_cast_non_const(Flux_radiatif_base, la_cl_rayo.valeur());
              return la_cl_rayon.flux_radiatif();
            }
          else
            {
              Cerr << "Error: the boundary conditions of the radiation equation" << finl;
              Cerr << "must necessarily be of radiation type" << finl;
              Process::exit();

            }
        }
    }
  Cerr << "Error: Pb_rayo_semi_transp::flux_radiatif" << finl;
  Cerr << "there is no boundary condition named " << nom_bord << finl;
  Process::exit();
  //for the compilers
  const Cond_lim& la_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites(0);
  Flux_radiatif_base& la_cl_rayon = ref_cast_non_const(Flux_radiatif_base, la_cl_rayo.valeur());
  return la_cl_rayon.flux_radiatif();
}
