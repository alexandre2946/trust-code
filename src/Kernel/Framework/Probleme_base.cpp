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
#include <Domaine_Cl_dis_base.h>
#include <Loi_Fermeture_base.h>
#include <Domaine_dis_base.h>
#include <Probleme_base.h>
#include <Synonyme_info.h>
#include <Postraitement.h>
#include <Perf_counters.h>
#include <EChaine.h>
#include <Perf_counters.h>
#include <Debog.h>

Implemente_base_sans_destructeur(Probleme_base,"Probleme_base",Probleme_U);

// XD pb_gen_base objet_u pb_gen_base INHERITS_BRACE Basic class for problems.

// XD Pb_base pb_gen_base Pb_base BRACE Resolution of equations on a domain. A problem is defined by creating an object
// XD_CONT and assigning the problem type that the user wishes to resolve. To enter values for the problem objects
// XD_CONT created, the Lire (Read) interpretor is used with a data block.
// XD attr milieu milieu_base milieu OPT The medium associated with the problem.
// XD attr constituant constituant constituant OPT Constituent.
// XD attr postraitement|Post_processing corps_postraitement postraitement OPT One post-processing (without name).
// XD attr postraitements|Post_processings postraitements postraitements OPT List of Postraitement objects (with name).
// XD attr liste_de_postraitements liste_post_ok liste_de_postraitements OPT This
// XD attr liste_postraitements liste_post liste_postraitements OPT This block defines the output files to be written
// XD_CONT during the computation. The output format is lata in order to use OpenDX to draw the results. This block can
// XD_CONT be divided in one or several sub-blocks that can be written at different frequencies and in different
// XD_CONT directories. Attention. The directory lata used in this example should be created before running the
// XD_CONT computation or the lata files will be lost.
// XD attr sauvegarde format_file_base sauvegarde OPT Keyword used when calculation results are to be backed up. When a
// XD_CONT coupling is performed, the backup-recovery file name must be well specified for each problem. In this case,
// XD_CONT you must save to different files and correctly specify these files when resuming the calculation.
// XD attr sauvegarde_simple format_file_base sauvegarde_simple OPT The same keyword than Sauvegarde except, the last
// XD_CONT time step only is saved.
// XD attr reprise format_file_base reprise OPT Keyword to resume a calculation based on the name_file file (see the
// XD_CONT class format_file). If format_reprise is xyz, the name_file file should be the .xyz file created by the
// XD_CONT previous calculation. With this file, it is possible to resume a parallel calculation on P processors,
// XD_CONT whereas the previous calculation has been run on N (N<>P) processors. Should the calculation be resumed,
// XD_CONT values for the tinit (see schema_temps_base) time fields are taken from the name_file file. If there is no
// XD_CONT backup corresponding to this time in the name_file, TRUST exits in error.
// XD attr resume_last_time format_file_base resume_last_time OPT Keyword to resume a calculation based on the name_file
// XD_CONT file, resume the calculation at the last time found in the file (tinit is set to last time of saved files).
// XD attr transparent_medium_radiation_model transparent_medium_radiation_model modele_rayonnement_milieu_transparent OPT Read a transparent medium radiation model associated to a fluid problem

// XD ref domaine domaine
// XD ref scheme schema_temps_base
// XD ref loi1 loi_fermeture_base
// XD ref loi2 loi_fermeture_base
// XD ref loi3 loi_fermeture_base
// XD ref loi4 loi_fermeture_base
// XD ref loi5 loi_fermeture_base

// XD problem_read_generic Pb_base problem_read_generic INHERITS_BRACE The probleme_read_generic differs rom the rest of
// XD_CONT the TRUST code : The problem does not state the number of equations that are enclosed in the problem. As the
// XD_CONT list of equations to be solved in the generic read problem is declared in the data file and not pre-defined
// XD_CONT in the structure of the problem, each equation has to be distinctively associated with the problem with the
// XD_CONT Associate keyword.
// XD ref eqn1 eqn_base
// XD ref eqn2 eqn_base
// XD ref eqn3 eqn_base
// XD ref eqn4 eqn_base
// XD ref eqn5 eqn_base
// XD ref eqn6 eqn_base
// XD ref eqn7 eqn_base
// XD ref eqn8 eqn_base
// XD ref eqn9 eqn_base
// XD ref eqn10 eqn_base


// XD constituant milieu_base constituant INHERITS_BRACE Constituent.
// XD attr coefficient_diffusion field_base coefficient_diffusion REQ Constituent diffusion coefficient value (m2.s-1).
// XD_CONT If a multi-constituent problem is being processed, the diffusivite will be a vectorial and each components
// XD_CONT will be the diffusion of the constituent.
// XD attr is_multi_scalar rien is_multi_scalar_diffusion OPT Flag to activate the multi_scalar diffusion operator

// XD format_file_base objet_lecture nul NO_BRACE Format of the file
// XD attr checkpoint_fname chaine checkpoint_fname REQ Name of file.

// XD binaire format_file_base binaire INHERITS_BRACE Format of the file - binary version
// XD formatte format_file_base formatte INHERITS_BRACE Format of the file - formatte version
// XD xyz format_file_base xyz INHERITS_BRACE Format of the file - xyz version
// XD single_hdf format_file_base single_hdf INHERITS_BRACE Format of the file - single_hdf version
// XD pdi format_file_base pdi INHERITS_BRACE Format of the file - pdi version

// XD pdi_expert format_file_base pdi_expert BRACE Format of the file - PDI expert version
// XD attr yaml_fname chaine yaml_fname REQ YAML file name

// Global variables to initialize est_le_premier_postraitement_pour_nom_fic
// and est_le_dernier_postraitement_pour_nom_fic in a single pass.
LIST(Nom) glob_noms_fichiers;
LIST(OBS_PTR(Postraitement)) glob_derniers_posts;

Sortie& Probleme_base::printOn(Sortie& os) const
{
  for (int i = 0; i < nombre_d_equations(); i++)
    os << equation(i).que_suis_je() << " " << equation(i) << finl;
  os << les_postraitements_;
  os << le_domaine_dis_.valeur();
  return os;
}

/*! @brief Reading of a problem in an input stream, and opening of the save stream.
 *
 *     Format:
 *      {
 *      medium_name block to read a medium
 *      equation_name block to read an equation
 *      Postraitement block to read postraitement
 *      reprise | sauvegarde | sauvegarde_simple
 *      formatte | binaire | pdi
 *      file_name
 *      }
 *  The sauvegarde_simple option allows saving the problem in the chosen file
 *  by overwriting each time the previous saves : this allows saving disk space.
 *
 * @param (Entree& is) input stream
 * @return (Entree&) the modified input stream
 * @throws no opening brace at the beginning of the format
 * @throws keyword "Postraitement" is not there
 * @throws save format must be "binary" or "formatted"
 * @throws no closing brace at the end of the data file
 */
Entree& Probleme_base::readOn(Entree& is)
{
  Cerr << "Reading of the problem " << le_nom() << finl;
  Motcle motlu;
  is >> motlu;
  if (motlu != "{")
    Process::exit("We expected { to start to read the problem !!! \n");

  /* 1 : solved_equations + milieu : NEW SYNTAX */
  lire_solved_equations(is);
  typer_lire_milieu(is);

  /* 2 : We read the equations */
  lire_equations(is, motlu); //"motlu" contains the first word after reading the equations

  /* 3 : Post-processing */
  // If the postraitement contains the word, we read another one...
  while (les_postraitements_.lire_postraitements(is, motlu, *this))
    {
      is >> motlu;
    }

  /* 4 : We complete ... */
  completer();
  Cerr << "Step verification of data being read in progress..." << finl;

  /* 5 : We verify ... */
  verifier();
  Cerr << "The read data are coherent" << finl;

  /* 6 : save/restart management ... */
  save_restart_.lire_sauvegarde_reprise(is, motlu);

  return is ;
}

void Probleme_base::typer_lire_milieu(Entree& is)
{
  // NOTA BENE :
  // Normally we have one medium per problem, except if the problem contains a concentration equation
  // In that case, we have an additional medium : constituent (don't ask why 2 media ... because I don't know !). to see if we can do better ...
  int nb_milieu = 1;

  // We search if it's a problem with concentration => with constituent
  const std::string conc = "Concentration", scal_pass = "Scalaires_Passifs", nom_pb = que_suis_je().getString();
  if (nom_pb.find(conc) != std::string::npos) nb_milieu = 2; // problem contains concentration !

  le_milieu_.resize(nb_milieu);

  for (int i = 0; i < nb_milieu; i++)
    {
      le_milieu_[i].typer_lire_simple(is, "Typing the medium ..."); // We start with reading the medium
      associer_milieu_base(le_milieu_[i].valeur()); // We associate it with each equation (virtual method for each problem ...)
    }

  // Medium(s) read ... Lets go ! We discretize the equations
  discretiser_equations();

  // uplifting of the unknown to the medium
  for (int i = 0; i < nombre_d_equations(); i++) equation(i).associer_milieu_equation();

  const bool is_constituant = nb_milieu == 1 ? false : true, is_scal_pass = (nom_pb.find(scal_pass) != std::string::npos) ? true : false;
  const int ns_ou_cond_eq = 0;
  int conc_eq = 1; // not const !

  // We discretize the medium/media ... And the concentration equation !
  if (is_constituant)
    {
      assert (nombre_d_equations() > 1);
      if (nombre_d_equations() == 2)
        for (int i = 0; i < nombre_d_equations(); i++) equation(i).milieu().discretiser((*this), la_discretisation_.valeur());
      else if (nombre_d_equations() == 3)
        {
          // We have 2 Cases :
          // - Pb_Thermohydraulique_Concentration (NS, Thermal, Conc)
          // - Pb_Hydraulique_Concentration_Scalaires_Passifs (NS, Conc + Equations_Scalaires_Passifs (the list) !)
          conc_eq = is_scal_pass ? 1 /* conc_eq */ : 2;
          equation(ns_ou_cond_eq).milieu().discretiser((*this), la_discretisation_.valeur()); // NS
          equation(conc_eq).milieu().discretiser((*this), la_discretisation_.valeur()); // Conc
        }
      else
        {
          // Rare case : Pb_Thermohydraulique_Concentration_Scalaires_Passifs
          // Here we have NS, Thermal, Conc + Equations_Scalaires_Passifs (the list) !
          assert (nombre_d_equations() == 4);
          conc_eq = 2;
          equation(ns_ou_cond_eq).milieu().discretiser((*this), la_discretisation_.valeur()); // NS
          equation(conc_eq).milieu().discretiser((*this), la_discretisation_.valeur()); // Conc
        }
    }
  else /* We discretize the medium of eq 1 and that's all :-) :-) */
    equation(ns_ou_cond_eq).milieu().discretiser((*this), la_discretisation_.valeur());
}

/*! @brief Reading of the equations of the problem.
 *
 */
Entree& Probleme_base::lire_equations(Entree& is, Motcle& mot)
{
  const int nb_eq = nombre_d_equations();
  is >> mot;
  if (nb_eq == 0)
    return is;

  bool already_read = true;

  if (mot == "correlations")
    {
      Cerr << "Reading of the correlations ..." << finl;
      lire_correlations(is);
      already_read = false;
    }

  if (mot == "Modele_rayonnement_milieu_transparent" || mot == "Transparent_medium_radiation_model")
    {
      Cerr << "Reading of the radiation model ..." << finl;
      lire_radiation_models(is, mot);
      already_read = false;
    }

  Cerr << "Reading of the equations ..." << finl;

  for (int i = 0; i < nb_eq; i++)
    {
      if (!already_read)
        is >> mot;

      is >> getset_equation_by_name(mot);
      already_read = false;
    }

  read_optional_equations(is, mot);

  return is;
}

Entree& Probleme_base::lire_radiation_models(Entree& is, Motcle& mot)
{
  Cerr << "The use of a transparent radiation model is not authorized for your problem " << que_suis_je() << " !!!" << finl;
  Cerr << "This model should only be used with a fluid problem. Update your data file." << finl;
  Process::exit();
  return is;
}

/*! @brief Associates the problem with all its equations.
 *
 */
void Probleme_base::associer()
{
  save_restart_.assoscier_pb_base(*this);

  for (int i = 0; i < nombre_d_equations(); i++)
    equation(i).associer_pb_base(*this);
}

/*! @brief override Objet_U::associer_(Objet_U& ob) Associates different objects with the problem by checking
 *
 *      the type of the object to be associated at execution.
 *      We can thus associate: a time scheme, a calculation domain.
 *      Uses the routine of the Type_Info class (Utilities)
 *
 * @param (Objet_U& ob) the object to associate
 * @return (int) 1 if association succeeded 0 otherwise 2 if the medium is already associated with another problem
 * @throws Unknown Objet_U type (not provided)
 */
int Probleme_base::associer_(Objet_U& ob)
{
  if (sub_type(Schema_Temps_base, ob))
    {
      associer_sch_tps_base(ref_cast(Schema_Temps_base, ob));
      return 1;
    }
  if (sub_type(Domaine, ob))
    {
      associer_domaine(ref_cast(Domaine, ob));
      return 1;
    }
  if (sub_type(Loi_Fermeture_base,ob))
    {
      Loi_Fermeture_base& loi=ref_cast(Loi_Fermeture_base,ob);
      liste_loi_fermeture_.add(loi);
      loi.associer_pb_base(*this);
      return 1;
    }
  if (sub_type(Domaine_32_64<trustIdType>, ob))
    {
      Cerr << "ERROR: You are trying to associate a 64-bit Domain to a Problem!" << finl;
      Cerr <<"  -> Keyword 'Domain_64' can *not* be used for a sequential run!" << finl;
      Cerr<< "  -> Keyword 'Domain_64' can only be used for initial partitioning. It must be changed into 'Domain' in the PAR_xxx dataset when running the parallel computation itself!" << finl;
      Process::exit();
    }
  if (sub_type(Milieu_base, ob))
    {
      Cerr << "YOU ARE USING AN OLD SYNTAX IN YOUR DATA FILE AND THIS IS NO MORE SUPPORTED !" << finl;
      Cerr << "STARTING FROM TRUST-v1.9.3 : THE MEDIUM SHOULD BE READ INSIDE THE PROBLEM AND NOT VIA ASSOSCIATION ... " << finl;
      Cerr << "HAVE A LOOK TO ANY TRUST TEST CASE TO SEE HOW IT SHOULD BE DONE ($TRUST_ROOT/tests/) ... " << finl;
      Cerr << "OR RUN -convert_data OPTION OF YOUR APPLICATION SCRIPT, FOR TRUST FOR EXAMPLE:" << finl;
      Cerr << "   trust -convert_data " << Objet_U::nom_du_cas() << ".data" << finl;
      Process::exit();
    }
  return 0;
}

/*! @brief Completes the equations associated with the problem.
 *
 * Filling of references, delegated to the equations.
 *
 */
void Probleme_base::completer()
{
  // Cerr << "Probleme_base::completer()" << finl;
  for (int i = 0; i < nombre_d_equations(); i++)
    equation(i).completer();

  for (auto& itr : liste_loi_fermeture_)
    itr->completer();

  les_postraitements_.completer();

  for (auto &corr : correlations_)
    corr.second->completer();
}

/*! @brief Verifies that the object is complete, coherent, .
 *
 * .. NOT DEVELOPED ALWAYS RETURNS 1
 *
 * @return (int) 1 if the object is correct
 */
int Probleme_base::verifier()
{
  return 1;
}

/*! @brief Associates a domain with the problem.
 *
 * Takes un_domaine as support.
 *      calls Domaine_dis::associer_dom(const Domaine& )
 *
 * @param (Domaine& un_domaine) the domain
 */
void Probleme_base::associer_domaine(const Domaine& un_domaine)
{
  le_domaine_ = un_domaine;
  statistics().record_nb_elem(un_domaine.nb_elem());
}

void Probleme_base::discretiser_equations()
{
  Cerr << "Discretization of the equations of problem " << que_suis_je() << " ..." <<  finl;
  for (int i = 0; i < nombre_d_equations(); i++)
    {
      equation(i).associer_domaine_dis(domaine_dis());
      equation(i).discretiser();
    }
}

/*! @brief Assigns a discretization to the problem Discretizes the Domain associated with the problem with the discretization
 *
 *      Associates the first Domain with the problem's equations
 *      Discretizes the equations associated with the problem
 *
 * @param (Discretisation_base& discretisation) a discretization for the problem
 */
void Probleme_base::discretiser(Discretisation_base& une_discretisation)
{
  associer();
  la_discretisation_ = une_discretisation;
  Cerr << "Discretization of the domain associated with the problem " << le_nom() << finl;

  if (!le_domaine_)
    Process::exit("ERROR: Discretize - You're trying to discretize a problem without having associated a Domain to it!!! Fix your dataset.");

  // Initialization of the renum_som_perio array
  le_domaine_->init_renum_perio();

  une_discretisation.associer_domaine(le_domaine_.valeur());
  le_domaine_dis_ = une_discretisation.discretiser();
  // Can not do this before, since the Domaine_dis is not typed yet:
  le_domaine_dis_->associer_domaine(le_domaine_);

  for (auto& itr : liste_loi_fermeture_)
    itr->discretiser(une_discretisation);
}

/*! @brief Flags the first and last post-processing for each file And initializes the post-processing
 *
 */
void Probleme_base::init_postraitements()
{
  for (auto& itr : les_postraitements_) // For each post-processing
    {
      OWN_PTR(Postraitement_base) &der_post = itr;

      // If it is of type Postraitement, initialize premier/dernier _pour_nom_fich
      if (sub_type(Postraitement, der_post.valeur()))
        {

          Postraitement& post = ref_cast(Postraitement, der_post.valeur());

          Nom nom_fichier = Sortie_Fichier_base::root;
          nom_fichier+=post.nom_fich();
          int rg = glob_noms_fichiers.rang(nom_fichier);
          if (rg == -1)   // This is the first time we encounter this name
            {
              glob_noms_fichiers.add(nom_fichier);
              glob_derniers_posts.add(post);
              post.est_le_premier_postraitement_pour_nom_fich() = 1;
            }
          else   // We have already seen this name
            {
              post.est_le_premier_postraitement_pour_nom_fich() = 0;
              Postraitement& autre_post = glob_derniers_posts[rg];
              autre_post.est_le_dernier_postraitement_pour_nom_fich() = 0;
              glob_derniers_posts[rg] = post;

              // Verify that the post-processing intervals are the same
              // for everything writing to the same file.
              if (post.champs_demande() && autre_post.dt_post() != post.dt_post())
                {
                  Cerr << "Error, the values of dt_post (" << autre_post.dt_post() << " and " << post.dt_post() << ") of two postprocessing blocks writing in the same file" << nom_fichier
                       << " are different!" << finl;
                  Cerr << "Specify the same dt_post, or use two different files for postprocessing." << finl;
                  Cerr << "INFO: For Probleme_couple, if you do not specify different filenames, postpros in subproblems must also have the same dt_post." << finl;
                  exit();
                }
            }
          post.est_le_dernier_postraitement_pour_nom_fich() = 1;
        }
    }
  les_postraitements_.init();
}

int Probleme_base::expression_predefini(const Motcle& motlu, Nom& expression)
{
  expression = "";
  return 0;
}

/*! @brief Writing of the problem to file for restart.
 *
 * Writes the name of the problem and saves the equations.
 *
 * @param (Sortie& os) output stream for backup
 * @return (int) always returns 1
 */
int Probleme_base::sauvegarder(Sortie& os) const
{
  Debog::set_nom_pb_actuel(le_nom());
  schema_temps().sauvegarder(os);
  Cerr << "Backup of problem " << le_nom() << finl;
  int bytes=0;
  for(int i=0; i<nombre_d_equations(); i++)
    {
      bytes += equation(i).sauvegarder(os);
      assert(bytes % 4 == 0);
    }
  bytes += domaine().save_additional_state(os, *this);  // save additional state (e.g. moving mesh data)
  bytes += les_postraitements_.sauvegarder(os);
  assert(bytes % 4 == 0); // To detect a sauvegarder() method which returns 1 instead of the number of bytes saved.
  return bytes;
}

/*! @brief Reading of an input stream (file) for restart after a backup with Probleme_base::sauvegarder(Sortie& os).
 *
 * @param (Entree& is) the input stream from which we read the restart
 * @return (int) always returns 1
 */
int Probleme_base::reprendre(Entree& is)
{
  statistics().begin_count(STD_COUNTERS::restart,statistics().get_last_opened_counter_level()+1);
  Debog::set_nom_pb_actuel(le_nom());
  schema_temps().reprendre(is);
  Cerr << "Resuming the problem " << le_nom() << finl;
  for(int i=0; i<nombre_d_equations(); i++)
    equation(i).reprendre(is);
  domaine().restore_additional_state(is, *this); //restore additional state (e.g. moving mesh data)
  les_postraitements_.reprendre(is);
  Cerr << "End of resuming the problem " << le_nom() << " after " << statistics().get_time_since_last_open(STD_COUNTERS::restart) << " s" << finl;
  statistics().end_count(STD_COUNTERS::restart);
  return 1;
}

/*! @brief Asks the time scheme if a print is needed
 *
 * @return (int) 1 a print is needed, 0 it's not.
 */
int Probleme_base::limpr() const
{
  return schema_temps().limpr();
}

/*! @brief Asks the time scheme if a backup is needed
 *
 * @return (int) 1 a backup is needed, 0 it's not.
 */
int Probleme_base::lsauv() const
{
  return schema_temps().lsauv();
}

/*! @brief Prints the equations associated with the problem if the associated time scheme indicates that it is necessary.
 *
 * @param (Sortie& os) output stream
 */
void Probleme_base::imprimer(Sortie& os) const
{
  for(int i=0; i<nombre_d_equations(); i++)
    equation(i).imprimer(os);
}

/*! @brief Associates a time scheme with the problem.
 *
 * Then associates the time scheme with all
 *     the equations of the problem.
 *
 * @param (Schema_Temps_base& un_schema_en_temps) the time scheme to associate
 */
void Probleme_base::associer_sch_tps_base(const Schema_Temps_base& un_schema_en_temps)
{
  if (le_schema_en_temps_)
    {
      Cerr << finl;
      Cerr<<"Error: Problem "<<le_nom()<<" was already associated with the scheme "<< le_schema_en_temps_->le_nom()<<" and we try to associate it with "<<un_schema_en_temps.le_nom() << "." <<finl;
      exit();
    }
  le_schema_en_temps_=un_schema_en_temps;
  le_schema_en_temps_->associer_pb(*this);
  for(int i=0; i<nombre_d_equations(); i++)
    equation(i).associer_sch_tps_base(un_schema_en_temps);
}

/*! @brief Returns the time scheme associated with the problem.
 *
 * (if it is not null) (const version)
 *
 * @return (Schema_Temps_base&) the time scheme associated with the problem
 * @throws the time scheme is not associated with the problem, the reference is null
 */
const Schema_Temps_base& Probleme_base::schema_temps() const
{
  if(!le_schema_en_temps_)
    {
      Cerr << le_nom() << " has not been associated to a time scheme !" << finl;
      exit();
    }
  return le_schema_en_temps_.valeur();
}


/*! @brief Returns the time scheme associated with the problem.
 *
 * (if it is not null)
 *
 * @return (Schema_Temps_base&) the time scheme associated with the problem
 * @throws if the time scheme is not associated with the problem, the reference is null
 */
Schema_Temps_base& Probleme_base::schema_temps()
{
  if(!le_schema_en_temps_)
    {
      Cerr << le_nom() << " has not been associated to a time scheme !" << finl;
      exit();
    }
  return le_schema_en_temps_.valeur();
}


/*! @brief Returns the domain associated with the problem.
 *
 * (const version)
 *
 * @return (Domaine&) a domain
 */
const Domaine& Probleme_base::domaine() const
{
  return le_domaine_.valeur();
}

/*! @brief Returns the domain associated with the problem.
 *
 * @return (Domaine&) a domain
 */
Domaine& Probleme_base::domaine()
{
  return le_domaine_.valeur();
}

/*! @brief Returns the discretized domain associated with the problem (const version).
 *
 * @return (Domaine_dis_base&) a discretized domain
 */
const Domaine_dis_base& Probleme_base::domaine_dis() const
{
  return le_domaine_dis_.valeur();
}

/*! @brief Returns the discretized domain associated with the problem.
 *
 * @return (Domaine_dis_base&) a discretized domain
 */
Domaine_dis_base& Probleme_base::domaine_dis()
{
  return le_domaine_dis_.valeur();
}

/*! @brief Associates a physical medium with the problem equations.
 *
 * Choice of the physical medium.
 *
 * @param (Milieu_base& mil) the medium to associate (Solide, Fluide Incompressible ...)
 */
void Probleme_base::associer_milieu_base(const Milieu_base& mil)
{
  for(int i=0; i<nombre_d_equations(); i++)
    equation(i).associer_milieu_base(mil);
}

/*! @brief Returns the physical medium associated with the problem (const version).
 *
 * Returns the medium associated with the first equation.
 *
 * @return (Milieu_base&) a physical medium
 */
const Milieu_base& Probleme_base::milieu() const
{
  return equation(0).milieu();
}

/*! @brief Returns the physical medium associated with the problem.
 *
 * Returns the medium associated with the first equation.
 *
 * @return (Milieu_base&) a physical medium
 */
Milieu_base& Probleme_base::milieu()
{
  return equation(0).milieu();
}

/*! @brief Returns the equation whose name is specified (const version).
 *
 * Equations are indexed by their associated name.
 *     Searches through all equations of the problem for the one
 *     carrying the specified name.
 *
 * @param (Nom& type) the name of the equation to return
 * @return (Equation_base&) an equation
 * @throws if no equation with the specified name exists
 */
const Equation_base& Probleme_base::equation(const Nom& type) const
{
  Motcle Type(type), Type_eqn;
  for (int i = 0; i < nombre_d_equations(); i++)
    {
      Type_eqn = equation(i).que_suis_je();
      if (Type_eqn == Type)
        return equation(i);

      // test if synonym ...
      const Synonyme_info *syn_info = Synonyme_info::synonyme_info_from_name(type);
      if (syn_info != 0)
        if (Motcle(syn_info->org_name_()) == Type_eqn)
          return equation(i);
    }
  Cerr << que_suis_je() << " does not contain any equation/medium of type: " << type << finl;
  Cerr << "Here is the list of possible equations for a " << que_suis_je() << " problem: " << finl;
  for (int i = 0; i < nombre_d_equations(); i++)
    Cerr << "\t- " << equation(i).que_suis_je() << finl;
  Process::exit();
  // For the compiler;
  return equation(0);
}

/*! @brief (B. Math): Virtual method added for problems having several equations
 *   of the same type (Probleme_FT_Disc_gen). In that case, the name of the equation
 *   is not its type...
 *
 */
const Equation_base& Probleme_base::get_equation_by_name(const Nom& un_nom) const
{
  return equation(un_nom);
}

/*! @brief (B. Math): Virtual method added for problems having several equations
 *   of the same type (Probleme_FT_Disc_gen). In that case, the name of the equation
 *   is not its type...
 *   Non-const version. This method is notably called when reading the problem.
 *
 */
Equation_base& Probleme_base::getset_equation_by_name(const Nom& un_nom)
{
  return equation(un_nom);
}

/*! @brief Returns the equation whose name is specified.
 *
 * Equations are indexed by their associated name.
 *     Searches through all equations of the problem for the one
 *     carrying the specified name.
 *
 * @param (Nom& type) the name of the equation to return
 * @return (Equation_base&) an equation
 * @throws if no equation with the specified name exists
 */
Equation_base& Probleme_base::equation(const Nom& type)
{
  Motcle Type(type), Type_eqn;
  for (int i = 0; i < nombre_d_equations(); i++)
    {
      Type_eqn = equation(i).que_suis_je();
      if (Type_eqn == Type)
        return equation(i);

      // test if synonym ...
      const Synonyme_info *syn_info = Synonyme_info::synonyme_info_from_name(type);
      if (syn_info != 0)
        if (Motcle(syn_info->org_name_()) == Type_eqn)
          return equation(i);
    }
  Cerr << que_suis_je() << " does not contain any equation/medium of type: " << type << finl;
  Cerr << "Here is the list of possible equations for a " << que_suis_je() << " problem: " << finl;
  for (int i = 0; i < nombre_d_equations(); i++)
    Cerr << "\t- " << equation(i).que_suis_je() << finl;
  Process::exit();
  // For the compiler;
  return equation(0);
}

void Probleme_base::creer_champ(const Motcle& motlu)
{
  domaine().creer_champ(motlu, *this);
  domaine_dis().creer_champ(motlu, *this);
  milieu().creer_champ(motlu);

  for (int i = 0; i < nombre_d_equations(); i++)
    equation(i).creer_champ(motlu);

  for (auto &itr : liste_loi_fermeture_)
    itr->creer_champ(motlu);
}

bool Probleme_base::has_champ(const Motcle& un_nom, OBS_PTR(Champ_base) &ref_champ) const
{
  if (domaine().has_champ(un_nom, ref_champ))
    return true;

  if (domaine_dis().has_champ(un_nom, ref_champ))
    return true;

  for (int i = 0; i < nombre_d_equations(); i++)
    {
      if (equation(i).has_champ(un_nom, ref_champ))
        return true;

      if (equation(i).milieu().has_champ(un_nom, ref_champ))
        return true;
    }

  for (const auto &corr : correlations_)
    if (corr.second->has_champ(un_nom, ref_champ))
      return true;

  for (const auto &itr : liste_loi_fermeture_)
    if (itr->has_champ(un_nom, ref_champ))
      return true;

  return false; /* nothing found */
}

bool Probleme_base::has_champ(const Motcle& un_nom) const
{
  if (domaine().has_champ(un_nom))
    return true;

  if (domaine_dis().has_champ(un_nom))
    return true;

  for (int i = 0; i < nombre_d_equations(); i++)
    {
      if (equation(i).has_champ(un_nom))
        return true;

      if (equation(i).milieu().has_champ(un_nom))
        return true;
    }

  for (const auto &corr : correlations_)
    if (corr.second->has_champ(un_nom))
      return true;

  for (const auto &itr : liste_loi_fermeture_)
    if (itr->has_champ(un_nom))
      return true;

  return false; /* nothing found */
}

const Champ_base& Probleme_base::get_champ(const Motcle& un_nom) const
{
  OBS_PTR(Champ_base) ref_champ;

  if (domaine().has_champ(un_nom, ref_champ))
    return ref_champ;

  if (domaine_dis().has_champ(un_nom, ref_champ))
    return ref_champ;

  for (int i = 0; i < nombre_d_equations(); i++)
    {
      if (equation(i).has_champ(un_nom, ref_champ))
        return ref_champ;

      if (equation(i).milieu().has_champ(un_nom, ref_champ))
        return ref_champ;
    }

  for (const auto &corr : correlations_)
    if (corr.second->has_champ(un_nom, ref_champ))
      return ref_champ;

  for (const auto &itr : liste_loi_fermeture_)
    if (itr->has_champ(un_nom, ref_champ))
      return ref_champ;

  Cerr << "The field of name " << un_nom << " do not correspond to a field understood by the problem." << finl;
  Cerr << "It may be a field dedicated only to post-process and defined in the Definition_champs set." << finl;
  Cerr << "1) If you have request the post-processing of " << un_nom << " in the Champs set" << finl;
  Cerr << "please remove the localisation elem or som that you may have specified." << finl;
  Cerr << "2) If you have used " << un_nom << " in Definition_champs, please use 'sources_reference { " << un_nom << " }'" << finl;
  Cerr << "instead of 'source refchamp { pb_champ " << le_nom() << " " << un_nom << " }'" << finl;
  Cerr << "3) Check reference manual." << finl;
  Cerr << "4) Contact TRUST support." << finl;
  Process::exit();

  throw std::runtime_error(std::string("Field ") + un_nom.getString() + std::string(" not found !"));
}

void Probleme_base::get_noms_champs_postraitables(Noms& noms,Option opt) const
{
  domaine().get_noms_champs_postraitables(noms, opt);
  domaine_dis().get_noms_champs_postraitables(noms, opt);
  milieu().get_noms_champs_postraitables(noms,opt);
  int nb_eq = nombre_d_equations();
  for (int i=0; i<nb_eq; i++)
    equation(i).get_noms_champs_postraitables(noms,opt);

  for (const auto& itr : liste_loi_fermeture_)
    {
      const Loi_Fermeture_base& loi=itr.valeur();
      loi.get_noms_champs_postraitables(noms,opt);
    }
}

int Probleme_base::comprend_champ_post(const Motcle& un_nom) const
{
  if (un_nom == "TEMPERATURE_PHYSIQUE") return 0;

  for (const auto &itr : postraitements())
    {
      const Postraitement& post = ref_cast(Postraitement, itr.valeur());
      if (post.comprend_champ_post(un_nom))
        return 1;
    }
  return 0;
}

bool Probleme_base::has_champ_post(const Motcle& un_nom) const
{
  for (const auto &itr : postraitements())
    if (sub_type(Postraitement, itr.valeur()))
      {
        const Postraitement& post = ref_cast(Postraitement, itr.valeur());
        if (post.has_champ_post(un_nom))
          return true;
      }

  return false; /* nothing found */
}

const Champ_Generique_base& Probleme_base::get_champ_post(const Motcle& un_nom) const
{
  for (const auto &itr : postraitements())
    if (sub_type(Postraitement, itr.valeur()))
      {
        const Postraitement& post = ref_cast(Postraitement, itr.valeur());
        if (post.has_champ_post(un_nom))
          return post.get_champ_post(un_nom);
      }

  Cerr << "The field named " << un_nom << " do not correspond to a field understood by the problem." << finl;
  Cerr << "Check the name of the field indicated into the postprocessing block of the data file " << finl;
  Cerr << "or in the list of post-processed fields above (in the block 'Reading of fields to be postprocessed')." << finl;
  Process::exit();

  throw std::runtime_error(std::string("Field ") + un_nom.getString() + std::string(" not found !"));
}

int Probleme_base::a_pour_IntVect(const Motcle&, OBS_PTR(IntVect)& ) const
{
  return 0;
}

/*! @brief Performs a time update of the problem.
 *
 * Performs the update on all equations of the problem.
 *
 * @param (double temps) the time step for the update
 */
void Probleme_base::mettre_a_jour(double temps)
{
  // Update the name of the problem being debugged
  Debog::set_nom_pb_actuel(le_nom());

  // Update the equations:
//statistics().begin_count(STD_COUNTERS::update_variables);
  for(int i=0; i<nombre_d_equations(); i++)
    equation(i).mettre_a_jour(temps);

  // Update the media:
  milieu().mettre_a_jour(temps);

  // Update the conserved fields in the equations (must be done after the media):
  for(int i=0; i<nombre_d_equations(); i++)
    equation(i).mettre_a_jour_champs_conserves(temps);

  // Update the post-processing:
  les_postraitements_.mettre_a_jour(temps);

  // Set conditions to update the domain in ALE:
  domaine().setUpdateTheGrid(true);


  for (auto& itr : liste_loi_fermeture_)
    {
      Loi_Fermeture_base& loi=itr.valeur();
      loi.mettre_a_jour(temps);
    }
  for (auto &corr : correlations_)
    corr.second->mettre_a_jour(temps);
  //statistics().end_count(STD_COUNTERS::update_variables);
}

/*! @brief Prepares the computation: initializes the medium parameters and prepares the computation of each equation.
 *
 */
void Probleme_base::preparer_calcul()
{
  const double temps = schema_temps().temps_courant();
  // Modification of the Qdm array held by domaine_dis() in the case
  // where there are periodic boundary conditions.
  // Note: if one of the equations has a periodic boundary condition
  //       then the others must necessarily have it too.
  equation(0).domaine_dis().modifier_pour_Cl(equation(0).domaine_Cl_dis().les_conditions_limites());
  milieu().initialiser(temps);
  for (int i = 0; i < nombre_d_equations(); i++)
    equation(i).preparer_calcul();
  milieu().preparer_calcul();
  for (int i = 0; i < nombre_d_equations(); i++) /* we can now fill the conserved fields */
    equation(i).mettre_a_jour_champs_conserves(temps);

  save_restart_.preparer_calcul();

  for (auto& itr : liste_loi_fermeture_)
    {
      Loi_Fermeture_base& loi = itr.valeur();
      loi.preparer_calcul();
    }

  if (correlations_.size() > 0)
    mettre_a_jour(temps);
}


/*! @brief Computes the value of the next time step for the problem.
 *
 * Computes the minimum of the time steps of the equations associated
 *     with the problem.
 *
 * @return (double) the maximum allowed time step for this problem
 */
double Probleme_base::calculer_pas_de_temps() const
{
  Debog::set_nom_pb_actuel(le_nom());
  double dt=schema_temps().pas_temps_max();
  for(int i=0; i<nombre_d_equations(); i++)
    dt=std::min(dt,equation(i).calculer_pas_de_temps());
  return dt;
}

void Probleme_base::lire_postraitement_interfaces(Entree& is)
{
  Cerr<<"The postprocessing of interfaces is only possible for"<<finl;
  Cerr<<"a problem type Pb_Front_Tracking, not a "<<que_suis_je()<<finl;
  exit();

}
void Probleme_base::postraiter_interfaces(const Nom& nomfich, Sortie& s, const Nom& format, double temps)
{
  Cerr<<que_suis_je()<<" must overloaded :postraiter_interfaces"<<finl;
  // exit();
}

bool Probleme_base::is_dilatable() const
{
  return milieu().is_dilatable();
}

/*! @brief Verifies that the necessary disk space exists.
 *
 */
void Probleme_base::allocation() const
{
  save_restart_.allocation();
}

/*! @brief If force=1, performs post-processing regardless of the post-processing frequencies.
 *
 *     The post-processing is updated and any treatments on
 *     probes, fields and statistics are performed.
 *   If force=0, respects the requested output frequencies.
 *
 */
int Probleme_base::postraiter(int force)
{
  statistics().begin_count(STD_COUNTERS::postreatment,statistics().get_last_opened_counter_level()+1);
  Schema_Temps_base& sch = schema_temps();
  Debog::set_nom_pb_actuel(le_nom());
  if (sch.nb_pas_dt() != 0)
    imprimer(Cout);
  if (force)
    {
      if (Process::nproc()>=100) Cerr << "[Post] Probleme_base::postraiter... " << finl;
      //Post-processable sources (Terme_Source_Acceleration) are not updated
      //for the final time and are not part of the champs_crees_ of the post-processing
      //which are updated by les_postraitements.mettre_a_jour.
      //We therefore update them here

      const int nb_pas_dt_max = sch.nb_pas_dt_max();
      bool& indice_nb_pas_dt = sch.set_indice_nb_pas_dt_max_atteint();
      bool& indice_tps_final = sch.set_indice_tps_final_atteint();
      const double t_init = sch.temps_init();
      const double t_max = sch.temps_max();

      //Test to avoid repeating post-processing at the initial time
      //This avoids a crash in the case of post-processing in meshtv format

      if (!(indice_nb_pas_dt && nb_pas_dt_max == 0) && !(indice_tps_final && est_egal(t_init, t_max)))
        {
          for (int i = 0; i < nombre_d_equations(); i++)
            equation(i).sources().mettre_a_jour(schema_temps().temps_courant());

          les_postraitements_.mettre_a_jour(schema_temps().temps_courant());
          les_postraitements_.postraiter();
          if (nb_pas_dt_max == 0)
            indice_nb_pas_dt = true;
          if (est_egal(t_init, t_max))
            indice_tps_final = true;
        }
      if (Process::nproc()>=100) Cerr << "[Post] Done in " << statistics().get_time_since_last_open(STD_COUNTERS::postreatment) << " s. If too slow consider using CGNS format." << finl;
    }
  else
    les_postraitements_.traiter_postraitement();

  statistics().end_count(STD_COUNTERS::postreatment);
  //Start specific postraitements for mobile domain (like ALE)
  if(!save_restart_.is_restart_in_progress() && le_domaine_dis_)
    {
      //no projection during the iteration of resumption of computation
      double temps = le_schema_en_temps_->temps_courant();
      le_domaine_dis_->domaine().update_after_post(temps);
    }

  save_restart_.set_restart_in_progress(false); //reset to false in order to make the following projections
  // end specific postraitements for mobile domain (like ALE)

  return 1;
}

/*! @brief Writes to file for restart (backup).
 *
 */
void Probleme_base::sauver() const
{
  statistics().begin_count(STD_COUNTERS::backup_file,statistics().get_last_opened_counter_level()+1);
  int bytes = save_restart_.sauver();
  Debog::set_nom_pb_actuel(le_nom());
  Cout << "[IO] " << statistics().get_time_since_last_open(STD_COUNTERS::backup_file) << " s to write save file." << finl;
  statistics().end_count(STD_COUNTERS::backup_file,1,bytes);
}

/*! @brief Finalizes post-processing and saves the problem to a file.
 *
 * Closes the file associated with post-processing. (Postraitement::finir())
 *
 */
void Probleme_base::finir()
{
  Debog::set_nom_pb_actuel(le_nom());
  schema_temps().finir(); // Close the .dt_ev file
  les_postraitements_.finir(); // Close post-processing files
  for (auto os : get_set_out_files())
    if (os->is_open())
      os->close(); // Close .out files
  // Clear the following global variables (useful for resetTime in the same directory)
  glob_noms_fichiers.vide();
  glob_derniers_posts.vide();

  if (schema_temps().temps_sauv() > 0.0)
    sauver();
  save_restart_.finir();
}

void Probleme_base::resetTime(double time)
{
  static const std::string param_name = "SORTIE_ROOT_DIRECTORY";
  std::string new_root_dir = (str_params_.count(param_name) == 0) ? "" : getOutputStringValue(param_name);

  resetTimeWithDir_impl(*this, time, new_root_dir);
}

/*! @brief Searches for parametric fields and for each, moves to the next parameter.
 *
 */
std::string Probleme_base::newCompute()
{
  std::string dirname="";
  // Loop over fields of boundary conditions:
  for (int i = 0; i < nombre_d_equations(); i++)
    {
      const Equation_base& eq = equation(i);
      const Conds_lim& condsLim = eq.domaine_Cl_dis().les_conditions_limites();
      for (auto const &condLim : condsLim)
        {
          const Cond_lim_base& la_cl_base = condLim.valeur();
          if (sub_type(Champ_front_Parametrique, la_cl_base.champ_front()))
            {
              const Champ_front_Parametrique& champ_front = ref_cast(Champ_front_Parametrique, la_cl_base.champ_front());
              dirname = champ_front.newCompute();
            }
        }
    }
  // Loop over source fields:
  if(Champ_Parametrique::enabled)
    {
      for (int i = 0; i < nombre_d_equations(); i++)
        {
          const Equation_base& eq = equation(i);
          const Sources& sources = eq.sources();
          for (auto &source: sources)
            {
              for (auto const &champ_don: source->champs_don())
                {
                  if (sub_type(Champ_Parametrique, champ_don.valeur()))
                    {
                      const Champ_Parametrique& champ = ref_cast(Champ_Parametrique, champ_don.valeur());
                      dirname = champ.newCompute();
                    }
                }
            }
        }
      // Loop over medium fields:
      for (auto const &champ_don: milieu().champs_don())
        {
          if (sub_type(Champ_Parametrique, champ_don.valeur()))
            {
              const Champ_Parametrique& champ = ref_cast(Champ_Parametrique, champ_don.valeur());
              dirname = champ.newCompute();
            }
        }
    }
  return dirname;
}

Entree& Probleme_base::read_optional_equations(Entree& is, Motcle& mot)
{
  /* reading of optional equations */
  Noms noms_eq, noms_eq_maj; //names of all possible equations!
  Type_info::les_sous_types(Nom("Equation_base"), noms_eq);
  for (auto& itr : noms_eq) noms_eq_maj.add(Motcle(itr)); //ha ha ha
  for (is >> mot; noms_eq_maj.rang(mot) >= 0; is >> mot)
    {
      eq_opt_.add(OWN_PTR(Equation_base)()); //another optional equation
      eq_opt_.dernier().typer(mot); //we give it the correct type
      Equation_base& eq = eq_opt_.dernier().valeur();
      //same associations as for the other equations: problem, medium, time scheme
      eq.associer_pb_base(*this);
      eq.associer_milieu_base(milieu());
      eq.associer_sch_tps_base(le_schema_en_temps_);
      eq.associer_domaine_dis(domaine_dis());
      eq.discretiser(); //must be done before reading the equation
      is >> eq; //and off we go!
      eq.associer_milieu_equation(); //propagate back to the medium
    }
  return is;
}

Entree& Probleme_base::lire_correlations(Entree& is)
{
  Motcle mot;
  is >> mot;
  if (mot != "{")
    {
      Cerr << "correlations : { expected instead of " << mot << finl;
      Process::exit();
    }
  for (is >> mot; mot != "}"; is >> mot)
    if (correlations_.count(mot.getString())) Process::exit(que_suis_je() + " : a correlation already exists for " + mot + " !");
    else
      Correlation_base::typer_lire_correlation(correlations_[mot.getString()], *this, mot, is);

  return is;
}

void Probleme_base::getOutputPointValues(const Nom& name,
                                         const std::vector<double>& x,
                                         const std::vector<double>& y,
                                         const std::vector<double>& z,
                                         std::vector<double>& vals, int compo)
{
  if (Process::is_parallel())
    Process::exit("Probleme_base::getOutputPointValues not implemented in // !! \n");

  assert (compo > -1);

  const int size_x = static_cast<int>(x.size());
  const int size_y = static_cast<int>(y.size());
  const int size_z = static_cast<int>(z.size());
  const int size_vals = static_cast<int>(vals.size());

  if (size_x != size_y)
    Process::exit("Error in Probleme_base::getOutputPointValues => vectors x and y must have same dimensions !!!");

  if (Objet_U::dimension > 2 && (size_x != size_z))
    Process::exit("Error in Probleme_base::getOutputPointValues => vectors x, y and z must have same dimensions !!!");

  if (size_vals == 0)
    vals.resize(size_x);
  else
    {
      if (size_vals != size_x)
        Process::exit("Error in Probleme_base::getOutputPointValues => vectors x, y/z and vals must have same dimensions !!!");
    }

  DoubleTrav les_positions; // TODO FIXME : attribute ?
  les_positions.resize(size_x, Objet_U::dimension);

  for (int i = 0; i < size_x; i++)
    {
      les_positions(i, 0) = x[i];
      les_positions(i, 1) = y[i];
      if (Objet_U::dimension > 2) les_positions(i, 2) = z[i];
    }

  IntVect elem; // TODO FIXME : attribute ?

  if(elem.size() != size_x)
    elem.resize(size_x);

  le_domaine_->chercher_elements(les_positions, elem, 1);

  // Check if some probes are outside the domain:
  ArrOfDouble tmp(size_x);
  for (int i = 0; i < size_x; i++)
    tmp[i] = elem[i];
  mp_max_for_each_item(tmp);
  for (int i = 0; i < size_x; i++)
    if (tmp[i] == -1)
      {
        Cerr << "Error in Probleme_base::getOutputPointValues => The point number " << i + 1 << " is outside the computational domain !!! " << finl;
        Process::exit();
      }

  DoubleTrav valeurs_locales;
  valeurs_locales.resize(size_x, 1);

  // TODO FIXME : remaining issue of som/grav ... to refactor with Sonde::initialiser()
  if (has_champ(Motcle(name)))
    {
      OBS_PTR(Champ_base) champ_ref = get_champ(Motcle(name));

      if (champ_ref->nb_comp() == 1)
        {
          assert(champ_ref->valeurs().line_size() == 1);
          champ_ref->valeur_aux_elems(les_positions, elem, valeurs_locales);
        }
      else
        {
          assert(compo < champ_ref->nb_comp());
          champ_ref->valeur_aux_elems_compo(les_positions, elem, valeurs_locales, compo);
        }
    }
  else /* from post like ICoCo*/
    {
      OBS_PTR(Champ_Generique_base) ref_ch = findOutputField(name);

      if (!ref_ch)
        {
          Cerr << "Error in Probleme_base::getOutputPointValues => No output fields of name " << name << finl;
          Process::exit();
        }

      OWN_PTR(Champ_base) espace_stockage;
      const Champ_base& ma_source = ref_ch->get_champ(espace_stockage);

      if (ma_source.nb_comp() == 1)
        {
          assert(ma_source.valeurs().line_size() == 1);
          ma_source.valeur_aux_elems(les_positions, elem, valeurs_locales);
        }
      else
        {
          assert(compo < ma_source.nb_comp());
          ma_source.valeur_aux_elems_compo(les_positions, elem, valeurs_locales, compo);
        }
    }

  for (int i = 0; i < size_x; i++)
    vals[i] = valeurs_locales(i, 0);
}
