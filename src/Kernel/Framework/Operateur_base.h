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

#ifndef Operateur_base_included
#define Operateur_base_included


#include <Equation_base.h>

#include <SolveurSys.h>
#include <TRUST_Ref.h>

#include <Matrice.h>

class Frontiere_dis_base;
class Matrice_Morse;
class EcrFicPartage;
class Conds_lim;
class SFichier;

/*! @brief class Operateur_base This class is the base of the hierarchy of objects representing an
 *
 *      operator used in TRUST equations. Its members are
 *      the attributes and methods common to all classes that
 *      represent an operator. An Operateur is a piece of an equation,
 *      which is why it derives from MorEqn, giving it a
 *      reference to the equation to which it is attached.
 *      Example of operator classes: Op_Diff_K_Eps_negligeable,
 *                                   Operateur_Div_base
 *
 * @sa MorEqn Operateur Equation_base
 *
 * Abstract class. A number of methods MUST be overridden in derived classes.
 */
class Operateur_base : public Objet_U, public MorEqn, public Champs_compris_interface
{
  Declare_base(Operateur_base);
public:
  virtual DoubleTab& ajouter(const DoubleTab&, DoubleTab&) const;
  virtual DoubleTab& calculer(const DoubleTab&, DoubleTab&) const;
  virtual void associer_champ(const Champ_Inc_base&, const std::string& nom_ch);
  virtual void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base& inco) =0;
  virtual void associer_domaine_cl_dis(const Domaine_Cl_dis_base&);
  virtual void dimensionner(Matrice_Morse&) const /* =0 */;
  virtual void dimensionner_bloc_vitesse(Matrice_Morse& matrice) const;
  virtual void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const /* =0 */;
  virtual void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const /* =0 */;
  /* allows filling dependencies on other variables */
  virtual void contribuer_bloc_vitesse(const DoubleTab&, Matrice_Morse&) const;
  virtual void contribuer_au_second_membre(DoubleTab&) const /* =0 */;
  void tester_contribuer_a_avec(const DoubleTab&, const Matrice_Morse&);

  /* interface {dimensionner,ajouter}_blocs -> cf Equation_base.h */
  virtual int has_interface_blocs() const { return 0; }
  virtual void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = { }) const;
  virtual void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = { }) const;

  virtual void dimensionner_termes_croises(Matrice_Morse&, const Probleme_base& autre_pb, int nl, int nc) const;
  virtual void ajouter_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, DoubleTab& resu) const;
  virtual void contribuer_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, Matrice_Morse& matrice) const;

  virtual double calculer_dt_stab() const;
  virtual void calculer_dt_local(DoubleTab&) const; //Local time step calculation
  virtual void completer();
  virtual void mettre_a_jour(double temps);
  virtual void abortTimeStep();
  virtual void resetTime(double time);
  virtual int impr(Sortie& os) const;
  inline void associer_eqn(const Equation_base&);
  inline int get_decal_temps() const;
  inline int set_decal_temps(int);
  inline int get_nb_ss_pas_de_temps() const;
  inline int set_nb_ss_pas_de_temps(int);
  inline const Matrice& get_matrice() const;
  inline Matrice& set_matrice();
  inline const SolveurSys& get_solveur() const;
  inline SolveurSys& set_solveur();
  inline Entree& lire_solveur(Entree&);
  virtual int systeme_invariant() const;
  virtual void ajouter_contribution_explicite_au_second_membre(const Champ_Inc_base& inconnue, DoubleTab& derivee) const;
  const Champ_Inc_base& mon_inconnue() const { return le_champ_inco.valeur(); }
  bool has_champ_inco() const { return bool(le_champ_inco); }
  const std::string& nom_inconnue() const
  {
    assert (has_champ_inco());
    return nom_inco_;
  }

  void ouvrir_fichier(SFichier& os, const Nom&, const int flag = 1) const;
  void ouvrir_fichier_partage(EcrFicPartage&, const Nom&, const int flag = 1) const;
  void set_fichier(const Nom&);
  inline const Nom fichier() const { return out_; }
  inline void set_description(const Nom& nom) { description_ = nom; }
  inline const Nom description() const { return description_; }
  inline DoubleTab& flux_bords() { return flux_bords_; }
  inline DoubleTab& flux_bords() const { return flux_bords_; }

  // Methods of the post-processable fields interface
  /////////////////////////////////////////////////////
  void creer_champ(const Motcle& motlu) override { }
  const Champ_base& get_champ(const Motcle& nom) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom, Option opt = NONE) const override;
  /////////////////////////////////////////////////////
  void calculer_pour_post(Champ_base& espace_stockage, const Nom& option, int comp) const override;
  Motcle get_localisation_pour_post(const Nom& option) const override;

  // Adding two methods for flux computation
  virtual void ajouter_flux(const DoubleTab& inconnue, DoubleTab& contribution) const;
  virtual void calculer_flux(const DoubleTab& inconnue, DoubleTab& flux) const;

  // Adding a preparer_calcul() method called during the equation's preparer_calcul()
  // This makes it easier to implement diffusion operators depending on
  // whether diffusivity varies or not.
  // The default implementation in Operateur_base.cpp does nothing
  virtual void preparer_calcul();
  int col_width_; // minimal size of a column for .out files (based on cl name length)
  bool has_impr_file() const { return out_ != "??"; }

protected:
  int decal_temps;
  int nb_ss_pas_de_temps;
  SolveurSys solveur;
  Matrice matrice_;
  Nom out_;                                 // Name of the .out file for printing
  Nom description_;
  mutable DoubleTab flux_bords_;         // Array containing the fluxes on the boundaries of the operator

  Champs_compris champs_compris_;
  OBS_PTR(Champ_Inc_base) le_champ_inco;
  std::string nom_inco_;
};


/*! @brief Associates an equation with the operator.
 *
 * Simple call to MorEqn::associer_eqn(const Equation_base&)
 *
 * @param (Equation_base& eqn) the equation to which the operator must be associated
 */
inline void Operateur_base::associer_eqn(const Equation_base& eqn)
{
  MorEqn::associer_eqn(eqn);
}
int Operateur_base::get_decal_temps() const
{
  return decal_temps;
}
int Operateur_base::set_decal_temps(int i)
{
  return decal_temps = i;
}
Entree& Operateur_base::lire_solveur(Entree& is)
{
  return is >> solveur;
}
const Matrice& Operateur_base::get_matrice() const
{
  return matrice_;
}
Matrice& Operateur_base::set_matrice()
{
  return matrice_;
}
int Operateur_base::get_nb_ss_pas_de_temps() const
{
  return nb_ss_pas_de_temps;
}
int Operateur_base::set_nb_ss_pas_de_temps(int i)
{
  return nb_ss_pas_de_temps = i;
}
const SolveurSys& Operateur_base::get_solveur() const
{
  return solveur;
}
SolveurSys& Operateur_base::set_solveur()
{
  return solveur;
}

#endif

