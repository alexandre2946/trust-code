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

#ifndef Schema_Temps_base_included
#define Schema_Temps_base_included

#include <Interface_blocs.h>
#include <TRUST_Deriv.h>
#include <TRUST_Ref.h>
#include <TRUSTTab.h>
#include <Parser_U.h>
#include <SFichier.h>
#include <math.h>

class Probleme_base;
class Equation_base;
class Matrice_Base;
class SFichier;
class Motcle;
class Param;

/*! @brief class Schema_Temps_base
 *
 * This class represents a time scheme, i.e. a particular
 * resolution algorithm that will be associated with a
 * Probleme_base (a simple problem and not a coupling).
 * Schema_Temps_base is the abstract class at the base of
 * the time scheme hierarchy.
 *
 * We denote n as the current time, and n+1 as the time at the end of the time step.
 * A time scheme computes u(n+1) knowing u up to u(n).
 * It uses u(n), but may also need past values of u,
 * such as u(n-1), ...
 * It may also during computation use values of u at intermediate
 * times between n and n+1, e.g. n+1/2.
 * nb_valeurs_temporelles counts all allocated values:
 * n, n+1, the retained past values and the intermediate values
 * between n and n+1.
 * nb_valeurs_futures counts n+1 and the intermediate values between n and n+1.
 * It is therefore the number of notches the wheels turn when advancing one time step.
 * temps_futur(i) returns the i-th future time value.
 * Finally temps_defaut is the time that fields must return
 * when valeurs() is called - particularly in operators.
 * Currently only respected by Champ_Front of boundary conditions.
 *
 * @sa Equation_base Probleme_base Algo_MG_base
 *
 * Abstract class from which all time schemes must derive.
 *
 * Abstract methods:
 *   int faire_un_pas_de_temps_eqn_base(Equation_base&)
 */

class Schema_Temps_base : public Objet_U
{
  Declare_base(Schema_Temps_base);

public :
  virtual void initialize();
  virtual void terminate() ;
  virtual double computeTimeStep(bool& stop) const;
  virtual bool initTimeStep(double dt);
  virtual void validateTimeStep();
  virtual bool isStationary() const;
  virtual void abortTimeStep();
  virtual void resetTime(double time);

  virtual bool iterateTimeStep(bool& converged);
  int limpr() const;

  ////////////////////////////////
  //                            //
  // Characteristics of the     //
  // time scheme                //
  //                            //
  ////////////////////////////////

  virtual int nb_valeurs_temporelles() const =0;
  virtual int nb_valeurs_futures() const =0;
  virtual double temps_futur(int i) const =0;
  virtual double temps_defaut() const =0;

  /////////////////////////////////////////
  //                                     //
  // End of time scheme characteristics  //
  //                                     //
  /////////////////////////////////////////

  inline void nommer(const Nom&) override;
  inline const Nom& le_nom() const override;
  virtual int faire_un_pas_de_temps_eqn_base(Equation_base&) =0;

  virtual void set_param(Param& titi) const override;
  int sauvegarder(Sortie& ) const override;
  int reprendre(Entree& ) override ;
  virtual int mettre_a_jour();
  virtual void mettre_a_jour_dt(double toto)  { abort();  }
  virtual void mettre_a_jour_dt_stab();

  inline double pas_de_temps() const;
  inline const DoubleTab& pas_de_temps_locaux() const;

  virtual bool corriger_dt_calcule(double& dt) const;
  virtual void imprimer(Sortie& os) const;
  virtual int impr(Sortie& os) const;
  void imprimer(Sortie& os,Probleme_base& pb) const;
  void imprimer(Sortie& os,const Probleme_base& pb) const;
  virtual int impr(Sortie& os,Probleme_base& pb) const;
  virtual int impr(Sortie& os,const Probleme_base& pb) const;
  void imprimer_temps_courant(SFichier&) const;
  inline double pas_temps_min() const;
  inline double& pas_temps_min();
  inline double pas_temps_max() const;
  inline double& pas_temps_max();
  inline int nb_impr() const;
  inline double temps_courant() const;
  inline double temps_precedent() const;
  inline double temps_calcul() const;
  virtual inline void changer_temps_courant(const double );
  void update_critere_statio(const DoubleTab& tab_critere, Equation_base& equation);
  inline double facteur_securite_pas() const;
  inline double& facteur_securite_pas();
  inline void notify_failed_timestep();
  virtual int stop() const;
  int lsauv() const;
  inline int temps_final_atteint() const;
  inline int nb_pas_dt_max_atteint() const;
  inline int temps_cpu_max_atteint() const;
  inline int stationnaire_atteint() const
  {
    assert(stationnaire_atteint_ != -1);
    return stationnaire_atteint_;
  }
  inline int stationnaire_atteint_safe() const { return stationnaire_atteint_; }
  int stop_lu() const;
  inline int diffusion_implicite() const;
  inline double seuil_diffusion_implicite() const { return seuil_diff_impl_; }
  inline int impr_diffusion_implicite() const { return impr_diff_impl_; }
  inline int impr_extremums() const { return impr_extremums_; }
  inline int niter_max_diffusion_implicite() const { return niter_max_diff_impl_;  }
  inline int no_conv_subiteration_diffusion_implicite() const { return no_conv_subiteration_diff_impl_;  }
  inline int no_error_if_not_converged_diffusion_implicite() const { return no_error_if_not_converged_diff_impl_; }
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  virtual Entree& lire_nb_pas_dt_max(Entree&);
  virtual Entree& lire_periode_sauvegarde_securite_en_heures(Entree&);
  virtual Entree& lire_temps_cpu_max(Entree&);
  virtual Entree& lire_residuals(Entree&);

  virtual void completer() =0;

  inline double temps_init() const ;
  inline double temps_max() const ;
  inline double temps_sauv() const ;
  inline int nb_sauv_max() const;
  inline double temps_impr() const ;
  inline int precision_impr() const { return precision_impr_; }
  inline int wcol() const
  {
    // minimum column width for .out files
    // precision_impr_ + 9 because: -1.000e+150 we add the length of "-1." and "e+150" plus one space
    return precision_impr_ + 9;
  }
  inline int gnuplot_header() const { return gnuplot_header_; }
  inline double seuil_statio() const ;
  inline int nb_pas_dt_max() const ;
  inline int nb_pas_dt() const ;
  inline double mode_dt_start() const { return mode_dt_start_; }
  inline int indice_tps_final_atteint() const { return ind_tps_final_atteint; }
  inline int indice_nb_pas_dt_max_atteint() const { return ind_nb_pas_dt_max_atteint; }
  inline int lu() const { return lu_; }
  inline int file_allocation() const { return file_allocation_; }

  inline double& set_temps_init() { return tinit_; }
  inline double& set_temps_max() { return tmax_; }
  inline double& set_temps_courant() { return temps_courant_; }
  inline double& set_temps_precedent() { return temps_precedent_; }
  inline int& set_nb_pas_dt() { return nb_pas_dt_; }
  inline int& set_nb_pas_dt_max() { return nb_pas_dt_max_; }
  inline double& set_dt_min() { return dt_min_; }
  inline double& set_dt_max()
  {
    dt_max_str_ = Nom(); //deactivates the function dt_max = f(t)
    return dt_max_;
  }
  inline double& set_dt_sauv() { return dt_sauv_; }
  inline double& set_dt_impr() { return dt_impr_; }
  inline int& set_precision_impr() { return precision_impr_; }
  inline double& set_dt() { return dt_; }
  inline double& set_facsec() { return facsec_; }
  inline double& set_seuil_statio() { return seuil_statio_; }
  inline int& set_stationnaire_atteint()
  {
    if (stationnaire_atteint_ == -1)
      stationnaire_atteint_ = 1;
    return stationnaire_atteint_;
  }
  inline void set_stationnaires_atteints(bool flag) { stationnaires_atteints_=flag; }
  inline int& set_diffusion_implicite() { return ind_diff_impl_; }
  inline double& set_seuil_diffusion_implicite() { return seuil_diff_impl_; }
  inline int& set_niter_max_diffusion_implicite() { return niter_max_diff_impl_; }
  inline double& set_mode_dt_start() { return mode_dt_start_; }
  inline bool& set_indice_tps_final_atteint() { return ind_tps_final_atteint; }
  inline bool& set_indice_nb_pas_dt_max_atteint() { return ind_nb_pas_dt_max_atteint; }
  inline bool& set_lu() { return lu_; }
  inline const double& residu() const { return residu_ ; }
  inline double& residu() { return residu_ ; }
  inline const Nom& norm_residu() const { return norm_residu_ ; }

  inline int& schema_impr() { return schema_impr_ ; }
  inline const int& schema_impr() const { return schema_impr_ ; }

  virtual void associer_pb(const Probleme_base&);
  Probleme_base& pb_base();
  const Probleme_base& pb_base() const;

  // for multi-step time schemes
  inline virtual void modifier_second_membre(const Equation_base& eqn, DoubleTab& secmem) { };

  // for implicit schemes: add inertia to matrix and right-hand side
  virtual void ajouter_inertie(Matrice_Base& mat_morse,DoubleTab& secmem,const Equation_base& eqn) const;

  //interface ajouter_blocs
  virtual void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const Equation_base& eqn, const tabs_t& semi_impl = {}) const;

  // Flag to disable the writing of the .progress file
  inline bool disable_progress() const
  {
    return disable_progress_ ;
  }
  void write_dt_ev(bool init);
  void write_progress(bool init);
  // Flag to disable the writing of the .dt_ev file
  inline bool disable_dt_ev() const
  {
    return disable_dt_ev_ ;
  }
  void finir() const;

protected :
  OBS_PTR(Probleme_base) mon_probleme;
  Nom nom_;
  double dt_ = 0.0;                         ///< Computation time step
  DoubleTab dt_locaux_;                     ///< Local time steps: Vector of size nb faces of the mesh

  double temps_courant_ = -100.;
  double temps_precedent_ = -100.;
  double dt_failed_ = -100.;            ///< Value of a time step if it failed
  double dt_gf_ = DMAXFLOAT;
  double tinit_ = -DMAXFLOAT;
  double tmax_ = 1.e30;
  double tcpumax_ = 1.e30;
  int nb_pas_dt_ = 0;
  int nb_pas_dt_max_ = std::numeric_limits<int>::max();
  mutable int nb_impr_ = 0;
  double dt_min_ = 1.e-16;               ///< Minimum time step set by the user
  mutable double dt_max_ = 1.e30;        ///< Maximum time step set by the user
  Nom dt_max_str_;                       ///< setting of dt_max as a function of time
  mutable Parser_U dt_max_fn_;           ///< Associated Parser_U
  double dt_stab_=-100.;                 ///<  Stability time step
  mutable double facsec_ = 1.;
  double seuil_statio_ = 1.e-12;
  int seuil_statio_relatif_deconseille_ = 0; ///< Flag to specify whether seuil_statio_ is an absolute (default) or relative value
  Nom norm_residu_;
  double dt_sauv_ = 1.e30;
  mutable int nb_sauv_ = 0;              ///< how many checkpoints have we performed so far?
  int nb_sauv_max_ = 10;                 ///< Max number of checkpoints that will be performed (useful for PDI backup file)
  double limite_cpu_sans_sauvegarde_ = 23 * 3600;  ///< Default 23 hours;
  double periode_cpu_sans_sauvegarde_ = 23 * 3600;  ///< Default 23 hours;
  double temps_cpu_ecoule_ = 0;
  double dt_impr_ = 1.e30;            ///< Output time interval
  int precision_impr_ = 8;            ///< Number of significant digits for output
  double mode_dt_start_ = -2;         ///< Mode for computing the initial time step - contains a double if dt_init option is used
  double residu_ = 0;
  double residu_old_slope_ = -1000;
  double cumul_slope_ = 1e-20;
  int gnuplot_header_ = 0;

  bool adapt_dt_tmax_ = false;
  bool ind_tps_final_atteint = false;
  bool ind_nb_pas_dt_max_atteint = false;
  bool ind_temps_cpu_max_atteint = false;
  bool lu_ = false;
  int ind_diff_impl_ = 0;
  double seuil_diff_impl_ = 1.e-6;   ///<  Threshold for implicit treatment of diffusion by CG
  int impr_diff_impl_ = 0;
  int impr_extremums_ = 0;
  int niter_max_diff_impl_ = 1000;      ///< Maximum iterations for CG diffusion implicitation - Above 1000 iterations, diffusion implicit algorithm may be diverging
  int no_conv_subiteration_diff_impl_ = 0;
  int no_error_if_not_converged_diff_impl_ = 0;
  int schema_impr_ = -1;                  // 1 if the scheme is allowed to print to .out and dt_ev files
  int file_allocation_ = 0;                // 1 = disk space allocation (default), 0 otherwise
  int max_length_cl_ = -10;
private:
  int stationnaire_atteint_ = 0;          ///< Stationary reached by the problem using this scheme
  int stationnaires_atteints_ = 0;       ///< Stationary reached by the calculation (means all the problems reach stationary)
  SFichier progress_;
  bool disable_progress_ = false;              ///< Flag to disable the writing of the .progress file
  bool disable_dt_ev_ = false;                 ///< Flag to disable the writing of the .dt_ev file
};

/*! @brief overrides Objet_U::nommer(const Nom&) Gives a name to the time scheme
 *
 * @param (Nom& name) the name to give to the time scheme
 */
inline void Schema_Temps_base::nommer(const Nom& name)
{
  nom_=name;
}

/*! @brief overrides Objet_U::le_nom() Returns the name of the time scheme
 *
 * @return (Nom&) the name of the time scheme
 */
inline const Nom& Schema_Temps_base::le_nom() const
{
  return nom_;
}

inline void Schema_Temps_base::notify_failed_timestep()
{
  dt_failed_ = dt_;
}

/*! @brief Returns a reference to the maximum number of time steps
 *
 * @return (double&)
 */
inline int Schema_Temps_base::nb_pas_dt_max() const
{
  return nb_pas_dt_max_;
}

/*! @brief Returns a reference to the stationarity threshold
 *
 * @return (double&)
 */
inline double Schema_Temps_base::seuil_statio() const
{
  return seuil_statio_;
}

/*! @brief Returns a reference to the output time interval
 *
 * @return (double&)
 */
inline double Schema_Temps_base::temps_impr() const
{
  return dt_impr_;
}

/*! @brief Returns a reference to the checkpoint time interval
 *
 * @return (double&)
 */
inline double Schema_Temps_base::temps_sauv() const
{
  return dt_sauv_;
}

/*! @brief Returns the maximum number of checkpoints (estimate)
 *
 * @return (int)
 */
inline int Schema_Temps_base::nb_sauv_max() const
{
  return nb_sauv_max_;
}

/*! @brief Returns a reference to the maximum time
 *
 * @return (double&)
 */
inline double Schema_Temps_base::temps_max() const
{
  return tmax_;
}

/*! @brief Returns the current time step (delta_t).
 *
 * @return (double) the current time step
 */
inline double Schema_Temps_base::pas_de_temps() const
{
  return dt_;
}
inline const DoubleTab& Schema_Temps_base::pas_de_temps_locaux() const
{
  return dt_locaux_;
}
/*! @brief Returns the minimum time step.
 *
 * (const version)
 *
 * @return (double) the minimum time step of the time scheme
 */
inline double Schema_Temps_base::pas_temps_min() const
{
  return dt_min_;
}

/*! @brief Returns a reference to the minimum time step.
 *
 * @return (double&) the minimum time step of the time scheme
 */
inline double& Schema_Temps_base::pas_temps_min()
{
  return dt_min_;
}

/*! @brief Returns the maximum time step.
 *
 * (const version)
 *
 * @return (double) the maximum time step of the time scheme
 */
inline double Schema_Temps_base::pas_temps_max() const
{
  return dt_max_;
}

/*! @brief Returns a reference to the maximum time step.
 *
 * @return (double&) the maximum time step of the time scheme
 */
inline double& Schema_Temps_base::pas_temps_max()
{
  return dt_max_;
}

/*! @brief Returns the current time.
 *
 * @return (double) the current time of the time scheme
 */
inline double Schema_Temps_base::temps_courant() const
{
  return temps_courant_;
}

/*! @brief Returns the previous time.
 *
 * @return (double) the previous time of the time scheme
 */
inline double Schema_Temps_base::temps_precedent() const
{
  return temps_precedent_;
}

/*! @brief Returns the elapsed computation time i.
 *
 * e. (current time - initial time).
 *
 * @return (double) the elapsed computation time
 */
inline double Schema_Temps_base::temps_calcul() const
{
  return temps_courant_ - tinit_;
}

/*! @brief Returns the initial time.
 *
 * @return (double) the initial time
 */
inline double Schema_Temps_base::temps_init() const
{
  return tinit_;
}

/*! @brief Returns the number of time steps performed.
 *
 * @return (int) the number of time steps performed
 */
inline int Schema_Temps_base::nb_pas_dt() const
{
  return nb_pas_dt_;
}

/*! @brief Returns the number of outputs performed.
 *
 * @return (int) the number of outputs performed
 */
inline int Schema_Temps_base::nb_impr() const
{
  return nb_impr_;
}

/*! @brief Changes the current time.
 *
 * @param (double& t) the new current time
 */
inline void Schema_Temps_base::changer_temps_courant(const double t)
{
  temps_courant_ = t;
}

/*! @brief Returns the safety factor or multiplier of delta_t.
 *
 * This factor is used during the correction/verification of
 *     the time step. See Schema_Temps_base::corriger_dt_calcule(double&)
 *     (const version)
 *
 * @return (double) the safety factor of the time scheme
 */
inline double Schema_Temps_base::facteur_securite_pas() const
{
  return facsec_;
}

/*! @brief Returns a reference to the safety factor or multiplier of delta_t.
 *
 * This factor is used during the correction/verification of
 *     the time step. See Schema_Temps_base::corriger_dt_calcule(double&)
 *
 * @return (double&) the safety factor of the time scheme
 */
inline double& Schema_Temps_base::facteur_securite_pas()
{
  return facsec_;
}

/*! @brief Returns 1 if the final time has been reached (or exceeded).
 *
 * Returns 1 if temps_courant_ >= tmax
 *     Returns 0 otherwise
 *
 * @return (int) 1 if the final time is reached, 0 otherwise
 * @throws final time reached
 */
inline int Schema_Temps_base::temps_final_atteint() const
{
  return ind_tps_final_atteint;
}

/*! @brief Returns 1 if (the number of time steps >= maximum number of time steps).
 *
 * Returns 0 otherwise
 *
 * @return (int) 1 if the maximum number of time steps is exceeded, 0 otherwise
 * @throws maximum number of time steps reached
 */
inline int Schema_Temps_base::nb_pas_dt_max_atteint() const
{
  return ind_nb_pas_dt_max_atteint;
}

inline int Schema_Temps_base::temps_cpu_max_atteint() const
{
  return ind_temps_cpu_max_atteint;
}

/*! @brief Returns 1 if the time scheme has been read with diffusion_implicite.
 *
 * @return (int) 1 if the time scheme has been read, 0 otherwise.
 */
inline int Schema_Temps_base::diffusion_implicite() const
{
  return ind_diff_impl_;
}

#endif /* Schema_Temps_base_included */
