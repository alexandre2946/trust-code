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

#ifndef Source_PDF_base_included
#define Source_PDF_base_included

#include <Interpolation_IBM_base.h>
#include <Prepro_IBM_base.h>
#include <Source_dep_inco_base.h>
#include <TRUST_Deriv.h>
#include <PDF_model.h>
#include <TRUSTTab.h>

class Probleme_base;
class Champ_Don_base;

/*! @brief class Source_PDF_base Base class for the source terms for the penalisation of the momentum in the Immersed Boundary Method (IBM)
 *
 *
 *
 * @sa Source_dep_inco_base, and Source_PDF
 */

class Source_PDF_base : public Source_dep_inco_base
{

  Declare_base(Source_PDF_base);

public:
  inline const int& getInterpolationBool() const { return interpolation_bool_; }
  inline const Interpolation_IBM_base& getInterpolationLu() const { return interpolation_lue_; }
  inline const Prepro_IBM_base& getpreproLu() const { return prepro_lu_; }
  void associer_pb(const Probleme_base& ) override;
  DoubleTab& ajouter(DoubleTab&) const override;
  DoubleTab& ajouter_(const DoubleTab&, DoubleTab&) const override;
  virtual DoubleTab& ajouter_(const DoubleTab&, DoubleTab&, const int) const ;
  DoubleTab& calculer(DoubleTab&) const override;
  DoubleTab& calculer(DoubleTab&, const int) const;
  DoubleTab& calculer_pdf(DoubleTab& ) const;
  void mettre_a_jour(double ) override;
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  double fonct_coeff(const double, const double, const double) const;
  double fonct_regul_PDF(const int, const double);
  virtual DoubleVect diag_coeff_elem(ArrOfDouble&, const DoubleTab&, int) const ;
  virtual DoubleTab compute_coeff_elem() const ;
  virtual DoubleTab compute_coeff_matrice() const;
  virtual void multiply_coeff_volume(DoubleTab&) const ;
  virtual void correct_variable(const DoubleTab&,DoubleTab&) const;
  void calculer_variable_imposee();
  void updateChampRho();
  inline const Champ_Don_base& get_champ_rotation() const { return champ_rotation_; }
  inline const Champ_Don_base& get_champ_aire() const { return champ_aire_; }
  inline const PDF_model& get_modele() const { return modele_lu_; }
  int impr(Sortie&) const override;
  void ouvrir_fichier(SFichier&, const Nom&, const int) const override;
  inline const DoubleTab& get_sec_mem_pdf() const { return sec_mem_pdf; }
  inline void set_sec_mem_pdf(DoubleTab& it) { sec_mem_pdf = it; }
  inline const DoubleTab& get_source_pdf() const { return source_term_PDF; }
  virtual inline void set_source_pdf(DoubleTab& it) { source_term_PDF = it; }
  DoubleVect& compute_source_term_PDF(int, DoubleTab&, DoubleVect&);
  DoubleVect& compute_source_term_PDF(DoubleVect&);
  virtual void volume_source_term_PDF(DoubleTab&);
  void update_elem_IBM(DoubleTab&, double, double);
  void update_pseudo_level_set_IBM(DoubleTab&, double);
  double aire_geometrique_IBM(DoubleTab&, int);
  void get_fields_from_prepro(Prepro_IBM_base&);
  void set_fields_to_prepro(Prepro_IBM_base&);

  inline const bool& get_matrice_pression_variable_bool_() const { return  matrice_pression_variable_bool_; }
  inline void set_matrice_pression_variable_bool_(bool& flag) { matrice_pression_variable_bool_ = flag; }
  virtual void correct_pressure(const DoubleTab&,DoubleTab&,const DoubleTab&) const;
  virtual void correct_incr_pressure(const DoubleTab&,DoubleTab&) const;
  virtual void filtre_CLD(DoubleTab&) const;
  inline void set_dt_computation_pdf(double dt) { dt_computation_pdf_ =  dt; }
  double calcul_dt_pdf();
  inline void set_temps_computation_pdf(double temps_pdf) { temps_computation_pdf_ =  temps_pdf; }
  inline double get_temps_pdf() {return temps_computation_pdf_;}
  void set_variable_imposee();
  void compute_NeighNode_IBM_elem(DoubleTab&, IntLists&, bool all_elem_vois = false);
  inline bool get_imm_wall_law() {return imm_wall_law_;}

  // Methodes de l interface des champs postraitables
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  const Champ_base& get_champ(const Motcle&) const override;
  void creer_champ(const Motcle& motlu) override;

protected:
  OWN_PTR(Champ_Don_base) champ_nodal_;
  virtual void compute_variable_imposee_projete(const DoubleTab&, const DoubleTab&, double, double);
  ArrOfDouble get_tuvw_local() const;
  void associer_domaines(const Domaine_dis_base&, const Domaine_Cl_dis_base&) override;
  virtual void compute_indicateur_nodal_champ_aire();
  void rotate_imposed_velocity(DoubleTab&);
  DoubleTab compute_pond(const DoubleTab& rho_m, const DoubleTab& aire, const DoubleVect& volume_thilde, int& nb_som_elem, int& nb_elems) const;
  virtual void calculer_variable_imposee_elem_fluid();
  virtual void calculer_variable_imposee_hybrid();
  virtual void calculer_variable_imposee_mean_grad();
  virtual void calculer_vitesse_imposee_power_law_tbl();
  virtual void calculer_vitesse_imposee_power_law_tbl_u_star();
  virtual void calculer_temperature_imposee_wall_law();
  int type_variable_imposee_ = -1;

  OWN_PTR(Champ_Don_base) champ_rotation_lu_, champ_rotation_, champ_aire_lu_, champ_aire_, champ_rho_;
  OWN_PTR(Champ_Don_base) champ_barycentre_lu_, champ_barycentre_;
  OWN_PTR(Prepro_IBM_base) prepro_lu_;
  void verify_results_prepro();

  bool transpose_rotation_ = false;
  DoubleTab indicateur_nodal_champ_aire_;
  IntTab indicateur_dead_cell_;
  DoubleTab variable_imposee_;
  PDF_model modele_lu_;
  double temps_relax_ = -100.;
  double echelle_relax_ = -100.;
  bool penalized_ = false;
  DoubleTab sec_mem_pdf; // part of the source term computed with the imposed variable
  double dt_computation_pdf_ = -123.; //time step used for computation of sec_mem_pdf
  double temps_computation_pdf_ = -123.; //time of the computation of sec_mem_pdf
  DoubleTab source_term_PDF; // PDF source term

  mutable OWN_PTR(Champ_Fonc_base)  champ_source_term_PDF_; //!< Champ pour postraitement
  mutable OWN_PTR(Champ_Fonc_base)  champ_barycentre_IBM_; //!< Champ pour postraitement
  mutable OWN_PTR(Champ_Fonc_base)  champ_aire_IBM_; //!< Champ pour postraitement
  mutable OWN_PTR(Champ_Fonc_base)  champ_normal_IBM_; //!< Champ pour postraitement
  mutable OWN_PTR(Champ_Fonc_base)  champ_vitesse_shape_IBM_; //!< Champ pour postraitement
  mutable OWN_PTR(Champ_Fonc_base)  champ_pseudo_level_set_IBM_; //!< Champ pour postraitement

  bool aire_from_prepro_ = false;
  bool rotation_from_prepro_ = false;
  bool barycentre_from_prepro_ = false;

  bool matrice_pression_variable_bool_ = false;
  bool imm_wall_law_ = false;

  // FOR THE INTERPOLATION
  OWN_PTR(Interpolation_IBM_base) interpolation_lue_;
  int interpolation_bool_ = 0;

  inline void comp_diff_L2_max(double dim0, double dim1, const DoubleTab& Array_src, const DoubleTab& Array_prep)
  {
    for (int k = 0; k < dim1; k++)
      {
        double err_L2 =0.;
        double err_max = 0.;
        double L2 = 0.;
        for (int dof = 0; dof <dim0; dof++)
          {
            double diff = abs(Array_prep(dof,k) - Array_src(dof,k));
            err_L2  += diff * diff ;
            L2 += Array_src(dof,k) * Array_src(dof,k) ;
            if (diff > err_max) err_max = diff;
          }
        if (L2 > 1.0e-12) err_L2 = sqrt(err_L2 / L2);

        Cerr<<"    composant # "<<k<<" => "<<err_L2<<" "<<err_max;
      }
    Cerr<<finl;
  }
};

#endif /* Source_PDF_base */
