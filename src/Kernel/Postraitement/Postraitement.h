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

#ifndef Postraitement_included
#define Postraitement_included

#include <Operateurs_Statistique_tps.h>
#include <Champ_Gen_de_Champs_Gen.h>
#include <Liste_Champ_Generique.h>
#include <Postraitement_base.h>
#include <Schema_Temps_base.h>
#include <Format_Post_base.h>
#include <Probleme_base.h>
#include <Champs_Fonc.h>
#include <Sondes_Int.h>
#include <TRUST_List.h>
#include <Interprete.h>
#include <TRUST_Ref.h>
#include <Parser_U.h>
#include <Sondes.h>

/*! @brief class Postraitement. The class holds -a list of generic fields champs_post_complet_ containing
 *
 *                                  all generic fields for post-processing
 *                            -a list of identifiers noms_champs_a_post_ containing the identifiers
 *                             of the fields to post-process
 *
 *      -Reading of generic fields declared in the "definition_champs" block (added to champs_post_complet_)
 *      and adding to this list the fields created by macro (triggered by old syntax in "champs" and "statistiques")
 *         Building the identifier list noms_champs_a_post_ during the reading operations
 *      -Updating generic fields, in practice updating the statistical operators for fields that carry one
 *      -post-processing performed in postraiter_champs()
 *         Iterating over identifiers
 *         The generic field corresponding to an identifier is retrieved by:             get_champ_post("identifiant")
 *         For this generic field: field values to write are computed by:                  get_champ(espace_stockage)
 *                                  supplementary info is retrieved by: get_property(), get_time() ...
 *         Writing the computed values:                                                postraiter(...)
 *
 *
 *   The macros used for creating generic fields are detailed in the .cpp of the class.
 *   See creer_champ_post() and creer_champ_post_stat()
 *
 *
 * @sa Hierarchy of generic fields (Champ_Generique_base), Syntax to follow in the data file, Postraitement {, Sondes, {, ..., }, Definition_champs, {, //Specification of generic fields, }, Champs, {, -Creation of generic field by macro if using old syntax (e.g. vitesse elem) and adding, identifier to noms_champs_a_post_ list, -Adding to the noms_champs_a_post_ list for a field declared in the "definition_champs" block, }, Statistiques, {, t_deb val_1 t_fin_val_2, -Creation of generic field by macro if using old syntax (e.g. Moyenne vitesse elem) and adding, identifier to noms_champs_a_post_ list, }, }
 */
class Postraitement : public Postraitement_base
{
  Declare_instanciable_sans_constructeur(Postraitement);
public:
  //
  // Overriding methods:
  //
  void associer_nom_et_pb_base(const Nom&, const Probleme_base&) override;
  void postraiter(int forcer) override;
  void mettre_a_jour(double temps) override;
  void finir() override;
  std::vector<YAML_data> data_a_sauvegarder() const override;
  int sauvegarder(Sortie& os) const override;
  int reprendre(Entree& is) override;
  void completer() override;
  void completer_sondes() override;
  void init() override;
  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void resetTime(double t, const std::string dirname) override;

  // specifique cgns : avoid duplicated file names
  void modify_cgns_basenames_and_reinit(const int, const int);

  //
  // Specific methods:
  //
  Postraitement();

  inline const Sondes& les_sondes() const { return les_sondes_; }
  inline Sondes& les_sondes() { return les_sondes_; }
  inline Probleme_base& probleme() { return mon_probleme.valeur(); }
  inline const Probleme_base& probleme() const { return mon_probleme.valeur(); }

  int postraiter_sondes();
  int traiter_sondes();
  virtual int postraiter_champs();
  // Called by postraiter_champs - only deal with the writing of the field values, not the geometrical parts:
  virtual void postprocess_field_values();

  int traiter_champs();
  virtual int lire_champs_a_postraiter(Entree& is, bool expect_acco);                //Optionally triggers creation of generic fields by macro
  //and builds the noms_champs_a_post_ list of post-processed fields
  int lire_champs_stat_a_postraiter(Entree&, bool expect_acco);        //same for statistics
  int lire_champs_operateurs(Entree& is);                //Reads a generic field, named and completed
  void complete_champ(Champ_Generique_base& champ,const Motcle& motlu);
  int postraiter_tableaux();
  int traiter_tableaux();
  int lire_tableaux_a_postraiter(Entree& );
  inline int lpost(double, double) const;
  inline int lpost_champ(double) const;
  inline int lpost_stat(double) const;
  inline int ind_post(int nb_pas_dt) const { return (nb_pas_dt%nb_pas_dt_post_==0) ? 1 : 0; }
  int nb_pas_dt_post() const { return nb_pas_dt_post_; }

  inline double dt_post() const { return dt_post_; }
  inline const Nom& format() const { return format_; }
  inline Nom nom_fich() const { return nom_fich_; }
  static inline LIST(Nom)& noms_fichiers_sondes() { return noms_fichiers_sondes_; }
  inline int& est_le_premier_postraitement_pour_nom_fich() { return est_le_premier_postraitement_pour_nom_fich_; }
  inline int& est_le_dernier_postraitement_pour_nom_fich() { return est_le_dernier_postraitement_pour_nom_fich_; }
  inline Operateurs_Statistique_tps& les_statistiques() { return les_statistiques_; }
  inline int sondes_demande() { return sondes_demande_; }
  inline int champs_demande() { return champs_demande_; }
  inline int stat_demande() const { return stat_demande_; }
  inline int stat_demande_definition_champs() const { return stat_demande_definition_champs_; }
  inline int tableaux_demande() { return tableaux_demande_; }
  inline bool besoin_postraiter_champs() { return (champs_demande_) || (stat_demande_) || (stat_demande_definition_champs_); }
  inline LIST(Nom)& noms_champs_a_post() { return noms_champs_a_post_; }
  inline const Liste_Champ_Generique& champs_post_complet() const { return champs_post_complet_; }

  //We distinguish post-processing of an array and of a tensor
  int postraiter(const Domaine& dom,const Noms& unites,const Noms& noms_compo,const int ncomp,
                 const double temps,
                 Nom nom_post,const Nom& localisation,const Nom& nature,const DoubleTab& valeurs,int tenseur);

  int postraiter_tableau(const Domaine& dom,const Noms& unites,const Noms& noms_compo,const int ncomp,
                         const double temps,
                         Nom nom_post,const Nom& localisation,const Nom& nature,const DoubleTab& valeurs);

  int postraiter_tenseur(const Domaine& dom,const Noms& unites,const Noms& noms_compo,const int ncomp,
                         const double temps,
                         Nom nom_post,const Nom& localisation,const Nom& nature,const DoubleTab& valeurs);


  virtual const Champ_Generique_base& get_champ_post(const Motcle& nom) const;
  virtual bool has_champ_post(const Motcle& nom) const;

  Nom set_expression_champ(const Motcle& motlu1,const Motcle& motlu2,
                           const Motcle& motlu3,const Motcle& motlu4,
                           const int trouve);

  //Macro methods to generate the creation of:
  //-Champ_Generique_Interpolation
  void creer_champ_post(const Motcle& motlu1,const Motcle& motlu2,Entree& s);
  //-Champ_Generique_Interpolation_Statistiques
  void creer_champ_post_stat(const Motcle& motlu1,const Motcle& motlu2,const Motcle& motlu3,const Motcle& motlu4,const double t_deb, const
                             double t_fin,Entree& s);
  //-Champ_Generique_Morceau_Equation
  void creer_champ_post_moreqn(const Motcle& type,const Motcle& option,const int num_eq,const int num_morceau,const int compo,Entree& s);

  //Macro method for the case of med fields
  void creer_champ_post_med(const Motcle& motlu1,const Motcle& motlu2,Entree& s);

  //Method comprend_champ_post() which indicates if the identifier corresponds to the name of a Champ_Generique_base
  //or to one of its components - sources are tested recursively
  int comprend_champ_post(const Motcle& identifiant) const;

  //temporaire a reviser
  void verifie_nom_et_sources(const Champ_Generique_base& champ);
  static Nom get_nom_localisation(const Entity& loc);

  int champ_fonc(Motcle& nom_champ, OBS_PTR(Champ_base)& mon_champ, OBS_PTR(Operateur_Statistique_tps_base)&
                 operateur_statistique) const;

  inline int& compteur_champ_stat();
  inline const double& tstat_deb() const;
  inline const double& tstat_fin() const;
  int cherche_stat_dans_les_sources(const Champ_Gen_de_Champs_Gen& ch, Motcle nom);

  /*! Calls by postraiter_champs() and allows a derived class to write extra meshes if needed
   * @return -1 if nothing more was written, 1 otherwise.
   */
  virtual int write_extra_mesh() { return -1; }
  const OBS_PTR(Domaine)& domaine() { return le_domaine_; }
  int DeprecatedKeepDuplicatedProbes=0; // Ancien format des sondes dans les .son qui autorise les sondes dupliquees

protected:

  int est_le_premier_postraitement_pour_nom_fich_, est_le_dernier_postraitement_pour_nom_fich_;
  double dt_post_;          ///< output of data (fields, stats, int_array) every dt_post (a time interval)
  int nb_pas_dt_post_;       ///< output of data (fields, stats, int_array) every dt_post (a period in number of iterations)
  Parser_U fdt_post_;

  Sondes les_sondes_;           // Probes to process
  Sondes_Int les_sondes_int_;   // Probes for integer arrays
  Operateurs_Statistique_tps les_statistiques_; // List of statistical operators to process

  LIST(Nom) noms_champs_a_post_;                 //contains the identifiers of fields to post-process
  Liste_Champ_Generique champs_post_complet_;   //contains all generic fields dedicated to post-processing

  //attributes for save-restart
  //redundant with attributes of Champ_Generique_Statistiques
  //The entire save-restart process should be managed by the statistical field
  int nb_champs_stat_;
  double tstat_deb_, tstat_fin_, tstat_dernier_calcul_;

  //Option to be managed by the statistical field
  int lserie_;
  double dt_integr_serie_;

  LIST(OBS_PTR(IntVect)) tableaux_a_postraiter_; // Liste de references a des tableaux a post-traiter
  LIST(Nom) noms_tableaux_;

  OWN_PTR(Format_Post_base) format_post_;
  std::vector<std::string> locs_required_;
  void add_locs_required_if_not(const Motcle& );

  static LIST(Nom) noms_fichiers_sondes_;
  bool sondes_demande_, champs_demande_, stat_demande_, stat_demande_definition_champs_, tableaux_demande_;
  int binaire_;
  Nom nom_fich_, format_, option_para_;
  Nom suffix_for_reset_; // Suffix appended to post base name when the method resetTime() was invoked - default to "_AFTER_RESET"
  double temps_, dernier_temps_; // time of the previous call to postraiter()
  OBS_PTR(Domaine) le_domaine_;
  OBS_PTR(Domaine_dis_base) domaine_dis_pour_faces_;
};


inline int Postraitement::lpost(double temps_courant, double dt_post) const
{
  double epsilon = 1.e-8;
  if (dt_post<=temps_courant - dernier_temps_)
    return 1;
  else
    {
      // See Schema_Temps_base::limpr for information about epsilon and modf
      double i, j;
      modf(temps_courant/dt_post + epsilon, &i);
      modf(dernier_temps_/dt_post + epsilon, &j);
      return ( i>j );
    }
}

/*! @brief Post-processing test taking into account the time evolution of the field.
 *
 *     Returns TRUE if the field requires post-processing, given
 *     the current time and the time step provided.
 *
 * @param (double temps_courant) current time
 * @param (double dt) the time step just completed
 * @return (int) boolean value, TRUE if the time step and current time provided indicate that post-processing is needed, FALSE otherwise.
 */
inline int& Postraitement::compteur_champ_stat() { return nb_champs_stat_; }
inline const double& Postraitement::tstat_deb() const { return tstat_deb_; }
inline const double& Postraitement::tstat_fin() const { return tstat_fin_; }

#endif /* Postraitement_included */
