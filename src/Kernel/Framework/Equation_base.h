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

#ifndef Equation_base_included
#define Equation_base_included

#include <Ecrire_fichier_xyz_valeur.h>
#include <Parametre_equation_base.h>
#include <Domaine_Cl_dis_base.h>
#include <Discretisation_base.h>
#include <Solveur_Masse_base.h>
#include <Matrice_Morse_Diag.h>
#include <MD_Vector_tools.h>
#include <Interface_blocs.h>
#include <Value_Input_Int.h>
#include <TRUSTTab_parts.h>
#include <Champ_Inc_base.h>
#include <Matrice_Morse.h>
#include <Ecrire_YAML.h>
#include <Champs_Fonc.h>
#include <TRUST_Ref.h>
#include <TRUSTList.h>
#include <TRUSTTrav.h>
#include <Parser_U.h>
#include <Sources.h>
#include <vector>

class Schema_Temps_base;
class Cond_lim_base;
class Milieu_base;
class Operateur;
class Motcle;
class Param;

enum Type_modele { TURBULENCE };

/*! @brief class Equation_base The role of an equation is the calculation of one or more fields. This class is the base of the equations hierarchy.
 *
 *      Its members are the attributes and methods common to all classes that represent equations.
 *      An equation is modeled in the following way:
 *
 *          M * dU_h/dt + Sum_i(Op_i(U_h)) = Sum(Sources);
 *
 *      M is the mass matrix represented by a "Solveur_Masse" object
 *      U_h is the unknown represented by a "Champ_Inc" object
 *      Op_i is the i-th operator of the equation represented by an "Operateur" object
 *      Sources are the source terms (possibly non-existent) of the equation represented by "Source" objects.
 *      An equation is linked to a problem by a reference contained in the member OBS_PTR(Probleme_base) mon_probleme.
 *
 *      Abstract class from which all equations must derive.
 *      Abstract methods:
 *        int nombre_d_operateurs() const
 *        const Operateur& operateur(int) const
 *        Operateur& operateur(int)
 *        const Champ_Inc_base& inconnue() const
 *        Champ_Inc_base& inconnue()
 *        void associer_milieu_base(const Milieu_base&)
 *        const Milieu_base& milieu() const
 *        Milieu_base& milieu()
 *        Entree& lire(const Motcle&, Entree&) [protected]
 *
 */
class Equation_base : public Champs_compris_interface, public Objet_U
{
  Declare_base(Equation_base);

public :
  // Overridden method from Objet_U:
  void nommer(const Nom& nom) override;
  // MODIF ELI LAUCOIN (22/11/2007) : adding a forward and a backward method
  virtual void avancer(int i=1);
  virtual void reculer(int i=1);
  // FIN MODIF ELI LAUCOIN (22/11/2007)

  virtual int nombre_d_operateurs() const =0;
  virtual int nombre_d_operateurs_tot() const;
  virtual const Operateur& operateur(int) const =0;
  virtual Operateur& operateur(int) =0;
  virtual const Operateur& operateur_fonctionnel(int) const;
  virtual Operateur& operateur_fonctionnel(int);
  virtual const Champ_Inc_base& inconnue() const =0;
  virtual Champ_Inc_base& inconnue() =0;
  virtual void associer_milieu_base(const Milieu_base&)=0;
  virtual const Milieu_base& milieu() const =0;
  virtual Milieu_base& milieu() =0;

  virtual std::vector<YAML_data> data_a_sauvegarder() const;
  int sauvegarder(Sortie&) const override;
  int reprendre(Entree&) override;
  Nom create_polymacfamily_syno(const Nom& field_tag) const;
  // if some equations need to save some parts of their data in a different backup file, we need to override these 2 methods below
  // (useful if some backup formats are not available for every equations)
  virtual void init_save_file() { }
  virtual void close_save_file() { }

  int limpr() const;
  virtual void imprimer(Sortie& os) const;
  virtual int impr(Sortie& os) const;
  virtual void associer_milieu_equation();

  virtual DoubleTab& derivee_en_temps_inco(DoubleTab& );
  virtual DoubleTab& derivee_en_temps_inco_transport(DoubleTab& derivee) { return derivee_en_temps_inco(derivee); }
  virtual DoubleTab& corriger_derivee_expl(DoubleTab& );
  virtual DoubleTab& corriger_derivee_impl(DoubleTab& );
  virtual void mettre_a_jour(double temps);
  virtual void abortTimeStep();
  virtual void resetTime(double time);
  virtual void valider_iteration();
  virtual int preparer_calcul();
  virtual bool initTimeStep(double dt);
  virtual bool updateGivenFields();
  virtual void discretiser();
  virtual void associer_pb_base(const Probleme_base&);
  virtual void completer();
  virtual double calculer_pas_de_temps() const;
  void calculer_pas_de_temps_locaux(DoubleTab&) const;  //Computation of local time: Vect of size number of faces of the domain
  Sources& sources();
  const Sources& sources() const;
  inline Solveur_Masse_base& solv_masse();
  inline const Solveur_Masse_base& solv_masse() const;
  Probleme_base& probleme();
  const Probleme_base& probleme() const;
  Schema_Temps_base& schema_temps();
  const Schema_Temps_base& schema_temps() const;
  virtual void associer_sch_tps_base(const Schema_Temps_base&);
  virtual void associer_domaine_dis(const Domaine_dis_base&);

  const Discretisation_base& discretisation() const;

  virtual inline Domaine_Cl_dis_base& domaine_Cl_dis();
  virtual inline const Domaine_Cl_dis_base& domaine_Cl_dis() const;
  Domaine_dis_base& domaine_dis();
  const Domaine_dis_base& domaine_dis() const;
  //
  inline const Nom& le_nom() const override;
  inline DoubleVect& get_residu() { return residu_; }
  inline DoubleVect& residu_initial() { return residu_initial_; }
  void initialise_residu(int=0);
  virtual void imprime_residu(SFichier&);
  virtual Nom expression_residu();

  // methods for implicit scheme
  virtual void dimensionner_matrice(Matrice_Morse& mat_morse); //memorizes the matrix stencil after the 1st call
  virtual void dimensionner_matrice_sans_mem(Matrice_Morse& mat_morse); //internal method called by the above

  // adds contributions from operators and sources
  virtual void assembler( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) ;
  // modifies the matrix and the right-hand side according to boundary conditions
  virtual void modifier_pour_Cl( Matrice_Morse& mat_morse,DoubleTab& secmem) const;
  // assembles, adds inertia, and modifies for boundary conditions.
  virtual void assembler_avec_inertie( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem) ;
  virtual void dimensionner_termes_croises(Matrice_Morse& matrice, const Probleme_base& autre_pb, int nl, int nc);
  virtual void ajouter_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, DoubleTab& resu) const;
  virtual void contribuer_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, Matrice_Morse& matrice) const;

  /*
    interface {dimensionner/ajouter/assembler}_blocs
    specificities: - has_interface_blocs() returns 1 if all terms of the equation support this interface
                   - dimensionner_blocs() not memoized (to be managed by the caller) / callable on non-empty matrices
                   - assembler_blocs() uses the values of unknowns/fields at the current time (like inconnue().valeurs())
                   - assembler_blocs_*() reasons in increments: M.dInco = S -> beware of solver thresholds
                   - certain variables (semi_impl set) can be treated as "semi-implicit"
                     (predicted values are used, no derivatives provided)
  */
  virtual int  has_interface_blocs() const;
  virtual double get_time_factor() const { return 1.; }
  virtual void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const;
  virtual void assembler_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const;
  virtual void assembler_blocs_avec_inertie(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {});

  /* auxiliary methods of the _blocs interface: field conserved by the equation and its values on Dirichlet or Neumann_val_ext boundary conditions
     by default, champ_conserve = temporal_coefficient * unknown
     this field is mutable so that the time scheme can update it
  */
  //the field: as many spatial/temporal values as the unknown
  Champ_Inc_base& champ_conserve() const { return champ_conserve_.valeur(); }
  int has_champ_conserve() const { return bool(champ_conserve_); }

  void init_champ_conserve() const; //to be called in completer() of operators/sources that will need champ_conserve_
  /* default computation function for champ_conserve */
  static void calculer_champ_conserve(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
  /* returns the name of the conserved field and the function to compute it -> to be overridden  */
  virtual std::pair<std::string, fonc_calc_t> get_fonc_champ_conserve() const
  {
    return { inconnue().le_nom().getString(), calculer_champ_conserve};
  }

  //by default the conserved field
  virtual Champ_Inc_base& champ_convecte() const { return champ_conserve_.valeur(); }
  virtual int has_champ_convecte() const { return bool(champ_conserve_); }
  virtual void init_champ_convecte() const { init_champ_conserve(); }
  //update of champ_conserve / champ_convecte: called by Probleme_base::mettre_a_jour() after updating the medium
  //if reset = 1, forces computation of all temporal values (not just the current value)
  virtual void mettre_a_jour_champs_conserves(double temps, int reset = 0);

  //Methods of the post-processable fields interface
  /////////////////////////////////////////////////////
  void creer_champ(const Motcle& motlu) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  /////////////////////////////////////////////////////

  virtual const Motcle& domaine_application() const;
  virtual void verifie_ch_init_nb_comp(const Champ_Inc_base& ch_ref, const int nb_comp) const;
  virtual void verifie_ch_init_nb_comp_cl(const Champ_Inc_base& ch_ref, const int nb_comp, const Cond_lim_base& cl) const
  {
    verifie_ch_init_nb_comp(ch_ref, nb_comp);
  }

  DoubleTab& derivee_en_temps_conv(DoubleTab& , const DoubleTab& );
  // Diffusion implicit scheme
  Matrice_Morse_Diag diag_;
  void Gradient_conjugue_diff_impl(DoubleTrav& secmem, DoubleTab& solution)
  {
    return Gradient_conjugue_diff_impl(secmem,solution,0,NULL_);
  }
  void Gradient_conjugue_diff_impl(DoubleTrav& secmem, DoubleTab& solution, const DoubleTab& terme_mul)
  {
    return Gradient_conjugue_diff_impl(secmem,solution,terme_mul.dimension_tot(0),terme_mul);
  }
  inline OWN_PTR(Parametre_equation_base)& parametre_equation() { return parametre_equation_ ; }
  inline const OWN_PTR(Parametre_equation_base)& parametre_equation() const { return parametre_equation_ ; }
  virtual const RefObjU& get_modele(Type_modele type) const;
  virtual int equation_non_resolue() const;
  int disable_equation_residual() const { return disable_equation_residual_; };

  //for multi-step time schemes
  inline virtual const Champ_Inc_base& derivee_en_temps() const { return derivee_en_temps_; }
  inline virtual Champ_Inc_base& derivee_en_temps() { return derivee_en_temps_; }
  void set_calculate_time_derivative(int i) { calculate_time_derivative_=i; }
  int calculate_time_derivative() const { return calculate_time_derivative_; }

  void set_residuals(const DoubleTab& residual);
  virtual bool positive_unkown() { return false; }

  inline void add_champs_compris(const Champ_base& ch) { champs_compris_.ajoute_champ(ch); };

  // set to true if operator is multiscalar (mixes components together). Only coded for VDF-Elem at present !
  inline void set_diffusion_multi_scalaire(bool flg = true)
  {
    if (flg) assert (discretisation().is_vdf());
    diffusion_multi_scalaire_ = flg;
  }
  inline const bool& diffusion_multi_scalaire() const { return diffusion_multi_scalaire_; }

  public_for_cuda
  void Gradient_conjugue_diff_impl(DoubleTrav& secmem, DoubleTab& solution, int size_terme_mul, const DoubleTab& term_mul);

protected :

  Nom nom_;
  OWN_PTR(Solveur_Masse_base) solveur_masse;
  Sources les_sources;
  OBS_PTR(Schema_Temps_base) le_schema_en_temps;
  OBS_PTR(Domaine_dis_base) le_dom_dis;
  OWN_PTR(Domaine_Cl_dis_base) le_dom_Cl_dis;
  OBS_PTR(Probleme_base) mon_probleme;
  virtual void set_param(Param& titi) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  virtual Entree& lire_sources(Entree&);
  virtual Entree& lire_cond_init(Entree&);
  virtual Entree& lire_cl(Entree&);
  virtual int verif_Cl() const;
  mutable DoubleList dt_op_bak;
  //Method lire with a specific signature to cause a compilation failure
  //if the old lire method is present
  //virtual Entree& lire(const Motcle&, Entree&)
  virtual void lire() { exit(); }

  int sys_invariant_;
  int implicite_;
  bool has_time_factor_; // Parameter set to 1 if convection has a prefactor (eg rhoCp in energy)
  OWN_PTR(Parametre_equation_base) parametre_equation_;

  LIST(RefObjU) liste_modeles_; //The first element of the list is the null model
  Champs_compris champs_compris_;
  Champs_Fonc list_champ_combi;

  //memoization of the matrix for PolyMAC_HFV
  mutable Matrice_Morse matrice_stockee;
  mutable int matrice_init;

  //for the assembler_blocs interface
  mutable OWN_PTR(Champ_Inc_base) champ_conserve_;
  mutable OWN_PTR(Champ_Inc_base) champ_convecte_;

  // For multistep methods, store previous dI/dt(n), dI/dt(n-1),...
  OWN_PTR(Champ_Inc_base) derivee_en_temps_;
  int calculate_time_derivative_;

  // for positivization of the term at the end of an iteration if necessary
  // returns 1 for a positive field, 0 for a negative field

  bool diffusion_multi_scalaire_ = false;

private :
  virtual void derivee_en_temps_inco_sources(DoubleTrav& ) { /* Don nothing */ }
  virtual void verify_scheme() { /* Don nothing */ }

  Ecrire_fichier_xyz_valeur xyz_field_values_file_;

  //!SC: moved to protected (override of get_champ in Equation_Diphasique_base)
//  Champs_Fonc list_champ_combi;
  DoubleVect residu_;
  DoubleVect residu_initial_;
  // returns the FIELD (not the norm) of the residuals for each unknown of the problem
  OWN_PTR(Champ_Fonc_base)  field_residu_;

  mutable DoubleTab NULL_;
  int disable_equation_residual_ = 0;
  mutable Parser_U equation_non_resolue_;
  Value_Input_Int eq_non_resolue_input_;
};


/*! @brief Returns the name of the equation.
 *
 * @return (Nom&) the name of the equation
 */
inline const Nom& Equation_base::le_nom() const
{
  return nom_;
}

/*! @brief Returns the discretized boundary condition domain associated with the equation.
 *
 * @return (Domaine_Cl_dis_base&) discretized boundary condition domain
 */
inline Domaine_Cl_dis_base& Equation_base::domaine_Cl_dis()
{
  assert(le_dom_Cl_dis);
  return le_dom_Cl_dis.valeur();
}

/*! @brief Returns the discretized boundary condition domain associated with the equation (const version).
 *
 * @return (Domaine_Cl_dis_base&) discretized boundary condition domain
 */
inline const Domaine_Cl_dis_base& Equation_base::domaine_Cl_dis() const
{
  assert(le_dom_Cl_dis);
  return le_dom_Cl_dis.valeur();
}

/*! @brief Returns the mass solver associated with the equation.
 *
 * @return (Solveur_Masse_base&) the mass solver associated with the equation
 */
inline Solveur_Masse_base& Equation_base::solv_masse()
{
  return solveur_masse;
}

/*! @brief Returns the mass solver associated with the equation (const version).
 *
 * @return (Solveur_Masse_base&) the mass solver associated with the equation
 */
inline const Solveur_Masse_base& Equation_base::solv_masse() const
{
  return solveur_masse;
}

#endif /* Equation_base_included */
