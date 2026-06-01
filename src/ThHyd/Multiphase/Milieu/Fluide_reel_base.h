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

#ifndef Fluide_reel_base_included
#define Fluide_reel_base_included

#include <Fluide_base.h>
#include <functional>
#include <span.hpp>
#include <Param.h>
#include <vector>
#include <array>
#include <map>

enum class Loi_en_T;
enum class Loi_en_h;

using MLoiSpanD = std::map<Loi_en_T, tcb::span<double>>;
using MLoiSpanD_h = std::map<Loi_en_h, tcb::span<double>>;
using MSpanD = std::map<std::string, tcb::span<double>>;
using VectorD = std::vector<double>;
using ArrayD = std::array<double,1>;
using SpanD = tcb::span<double>;


/*! @brief Represents a real fluid and its physical properties:
 *
 *         - kinematic viscosity, (mu)
 *         - dynamic viscosity,   (nu)
 *         - density,             (rho)
 *         - diffusivity,         (alpha)
 *         - thermal conductivity,(lambda)
 *         - heat capacity,       (Cp)
 *         - thermal expansion coefficient (beta_co)
 *
 * @sa Milieu_base
 */
class Fluide_reel_base: public Fluide_base
{
  Declare_base_sans_constructeur(Fluide_reel_base);
public :
  Fluide_reel_base()
  {
    converter_H_to_T_.set_instance(*this);
    converter_T_to_H_.set_instance(*this);
  }

  bool initTimeStep(double dt) override;
  int initialiser(const double temps) override;
  int check_unknown_range() const override; // verifies that each unknown "inco" lies between val_min[inco] and val_max[inco]
  int is_incompressible() const override { return (P_ref_ >= 0 && T_ref_ >= 0) || (P_ref_ >= 0 && h_ref_ >= 0); }
  void abortTimeStep() override;
  void mettre_a_jour(double temps) override;
  void preparer_calcul() override;
  void set_param(Param& param) const override;
  void discretiser(const Probleme_base& pb, const  Discretisation_base& dis) override;
  void creer_champs_non_lus() override { /* everything is done in discretiser */ }

  // range[inco] = { min, max}: by default, nothing to check
  virtual std::map<std::string, std::array<double, 2>> unknown_range() const { return {}; }
  virtual std::map<std::string, std::array<double, 2>> unknown_range_h() const { return {}; }

  // Methods used only in Pb_Euler
  inline virtual double calculer_vitesse_son(const double& rho, const double& p) const
  {
    Process::exit("Fluide_reel_base::calculer_vitesse_son is not implemented for your fluid !!! To call only for Pb_Euler also ... \n");
    return 0.;
  }
  inline virtual double calculer_pression(const double& rho, const double& rhou, const double& rhoE) const
  {
    Process::exit("Fluide_reel_base::calculer_vitesse_son is not implemented for your fluid !!! To call only for Pb_Euler also ... \n");
    return 0.;
  }
  inline virtual double init_energie_tot(const double& rho, const double& norm_U, const double& u) const
  {
    Process::exit("Fluide_reel_base::calculer_vitesse_son is not implemented for your fluid !!! To call only for Pb_Euler also ... \n");
    return 0.;
  }

protected :
  double T_ref_ = -1., P_ref_ = -1., h_ref_ = -1., t_init_ = -1.;
  int first_maj_ = 1;
  bool res_en_T_ = true; // by default resolution in T

  void calculate_fluid_properties_incompressible();
  void calculate_fluid_properties();

  void calculate_fluid_properties_enthalpie_incompressible();
  void calculate_fluid_properties_enthalpie();

  /*
   * *******************
   * For compressible
   * *******************
   */

  /* Laws in T */
  // density
  virtual void rho_(const SpanD T, const SpanD P, SpanD R, int ncomp = 1, int id = 0) const = 0;
  virtual void dP_rho_(const SpanD T, const SpanD P, SpanD dP_R, int ncomp = 1, int id = 0) const = 0;
  virtual void dT_rho_(const SpanD T, const SpanD P, SpanD dT_R, int ncomp = 1, int id = 0) const = 0;

  // enthalpy
  virtual void h_(const SpanD T, const SpanD P, SpanD H, int ncomp = 1, int id = 0) const = 0;
  virtual void dP_h_(const SpanD T, const SpanD P, SpanD dP_H, int ncomp = 1, int id = 0) const = 0;
  virtual void dT_h_(const SpanD T, const SpanD P, SpanD dT_H, int ncomp = 1, int id = 0) const = 0;

  // "weak" field laws -> no derivatives
  virtual void cp_(const SpanD T, const SpanD P, SpanD CP, int ncomp = 1, int id = 0) const = 0;
  virtual void beta_(const SpanD T, const SpanD P, SpanD B, int ncomp = 1, int id = 0) const = 0;
  virtual void mu_(const SpanD T, const SpanD P, SpanD M, int ncomp = 1, int id = 0) const = 0;
  virtual void lambda_(const SpanD T, const SpanD P, SpanD L, int ncomp = 1, int id = 0) const = 0;

  // application-specific methods to improve performance: used in Pb_Multiphase (for now!)
  virtual void compute_CPMLB_pb_multiphase_(const MSpanD , MLoiSpanD, int ncomp = 1, int id = 0) const;
  virtual void compute_all_pb_multiphase_(const MSpanD , MLoiSpanD, MLoiSpanD , int ncomp = 1, int id = 0) const;

  // Methods that can be called if point-to-point calculation is required
  double _rho_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::rho_>(T,P); }
  double _dP_rho_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::dP_rho_>(T,P); }
  double _dT_rho_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::dT_rho_>(T,P); }

  double _h_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::h_>(T,P); }
  double _dP_h_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::dP_h_>(T,P); }
  double _dT_h_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::dT_h_>(T,P); }

  double _cp_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::cp_>(T,P); }
  double _beta_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::beta_>(T,P); }
  double _mu_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::mu_>(T,P); }
  double _lambda_(const double T, const double P) const { return double_to_span<&Fluide_reel_base::lambda_>(T,P); }

  /* Laws in h */
  // density
  virtual void rho_h_(const SpanD h, const SpanD P, SpanD R, int ncomp = 1, int id = 0) const = 0;
  virtual void dP_rho_h_(const SpanD h, const SpanD P, SpanD dP_R, int ncomp = 1, int id = 0) const = 0;
  virtual void dh_rho_h_(const SpanD h, const SpanD P, SpanD dT_R, int ncomp = 1, int id = 0) const = 0;

  // temperature (from h)
  virtual void T_(const SpanD h, const SpanD P, SpanD H, int ncomp = 1, int id = 0) const = 0;
  virtual void dP_T_(const SpanD h, const SpanD P, SpanD dP_H, int ncomp = 1, int id = 0) const = 0;
  virtual void dh_T_(const SpanD h, const SpanD P, SpanD dT_H, int ncomp = 1, int id = 0) const = 0;

  // "weak" field laws -> no derivatives
  virtual void cp_h_(const SpanD h, const SpanD P, SpanD CP, int ncomp = 1, int id = 0) const = 0;
  virtual void beta_h_(const SpanD h, const SpanD P, SpanD B, int ncomp = 1, int id = 0) const = 0;
  virtual void mu_h_(const SpanD h, const SpanD P, SpanD M, int ncomp = 1, int id = 0) const = 0;
  virtual void lambda_h_(const SpanD h, const SpanD P, SpanD L, int ncomp = 1, int id = 0) const = 0;

  // application-specific methods to improve performance: used in Pb_Multiphase (for now!)
  virtual void compute_CPMLB_pb_multiphase_h_(const MSpanD , MLoiSpanD_h, int ncomp = 1, int id = 0) const;
  virtual void compute_all_pb_multiphase_h_(const MSpanD , MLoiSpanD_h, MLoiSpanD_h , int ncomp = 1, int id = 0) const;

  // Methods that can be called if point-to-point calculation is required
  double _rho_h_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::rho_h_>(h,P); }
  double _dP_rho_h_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::dP_rho_h_>(h,P); }
  double _dh_rho_h_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::dh_rho_h_>(h,P); }

  double _T_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::T_>(h,P); }
  double _dP_T_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::dP_T_>(h,P); }
  double _dh_T_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::dh_T_>(h,P); }

  double _cp_h_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::cp_h_>(h,P); }
  double _beta_h_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::beta_h_>(h,P); }
  double _mu_h_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::mu_h_>(h,P); }
  double _lambda_h_(const double h, const double P) const { return double_to_span<&Fluide_reel_base::lambda_h_>(h,P); }

private:
  typedef void(Fluide_reel_base::*function_span_generic)(const SpanD , const SpanD , SpanD , int , int ) const;

  template <function_span_generic FUNC>
  void double_to_span(const double T_ou_h, const double P, SpanD res) const
  {
    ArrayD Tt = {T_ou_h}, Pp = {P}, res_ = {0.};
    (this->*FUNC)(SpanD(Tt), SpanD(Pp), SpanD(res_),1,0); // fill res_
    for (auto& val : res) val = res_[0]; // fill res
  }

  template <function_span_generic FUNC>
  double double_to_span(const double T_ou_h, const double P) const
  {
    ArrayD Tt = {T_ou_h}, Pp = {P}, res_ = {0.};
    (this->*FUNC)(SpanD(Tt), SpanD(Pp), SpanD(res_),1,0); // fill res_
    return res_[0];
  }

  /*
   * *********************
   * For incompressible
   * *********************
   */

  /* Laws in T */
  void _rho_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::rho_>(T,P,res); }
  void _dP_rho_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dP_rho_>(T,P,res); }
  void _dT_rho_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dT_rho_>(T,P,res); }

  void _h_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::h_>(T,P,res); }
  void _dP_h_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dP_h_>(T,P,res); }
  void _dT_h_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dT_h_>(T,P,res); }

  void _cp_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::cp_>(T,P,res); }
  void _beta_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::beta_>(T,P,res); }
  void _mu_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::mu_>(T,P,res); }
  void _lambda_(const double T, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::lambda_>(T,P,res); }

  void _compute_CPMLB_pb_multiphase_(MLoiSpanD ) const;
  void _compute_all_pb_multiphase_(MLoiSpanD , MLoiSpanD ) const;

public:
  /*
   * Elie Saikali: internal struct to convert derivatives from h to T (for Pb_Multiphase).
   *
   * We seek dX/dP|T and dX/dT|P
   *
   * We know that dX = dX/dP|T dP + dX/dT|P dT = dX/dP|h dP + dX/dh|P dh
   *
   * Expanding:
   *
   *    dX/dP|h dP + dX/dh|P dh = dX/dP|h dP + dX/dh|P {  dh/dP|T dP + dh/dT|P dT }
   *                            = { dX/dP|h + dX/dh|P * dh/dP|T } dP + { dX/dh|P * dh/dT|P } dT
   *                            = dX/dP|T dP + dX/dT|P dT
   *
   * Therefore,
   *
   *    dX/dP|T = dX/dP|h + dX/dh|P * dh/dP|T
   *    dX/dT|P = dX/dh|P * dh/dT|P
   *
   *    /         \   /             \   /         \
   *    | dX/dP|T |   | 1  dh/dP|T  |   | dX/dP|h |
   *    |         | = |             | * |         |
   *    | dX/dT|P |   | 0  dh/dT|P  |   | dX/dh|P |
   *    \         /   \             /   \         /
   */


  struct H_to_T
  {
    void set_instance(const Fluide_reel_base& fld) { z_fld_ = fld; }

    void dX_dP_T(const SpanD dX_dP_h, const SpanD dX_dh_P, SpanD dX_dP);
    void dX_dT_P(const SpanD dX_dP_h, const SpanD dX_dh_P, SpanD dX_dT);
  private:
    OBS_PTR(Fluide_reel_base) z_fld_;
  };

  struct T_to_H
  {
    void set_instance(const Fluide_reel_base& fld) { z_fld_ = fld; }

    void dX_dP_h(const SpanD dX_dP_T, const SpanD dX_dT_P, SpanD dX_dP);
    void dX_dh_P(const SpanD dX_dP_T, const SpanD dX_dT_P, SpanD dX_dh);
  private:
    OBS_PTR(Fluide_reel_base) z_fld_;
  };

  H_to_T converter_H_to_T_;
  T_to_H converter_T_to_H_;

  /* Laws in h */
  void _rho_h_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::rho_h_>(h,P,res); }
  void _dP_rho_h_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dP_rho_h_>(h,P,res); }
  void _dh_rho_h_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dh_rho_h_>(h,P,res); }

  void _T_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::T_>(h,P,res); }
  void _dP_T_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dP_T_>(h,P,res); }
  void _dh_T_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::dh_T_>(h,P,res); }

  void _cp_h_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::cp_h_>(h,P,res); }
  void _beta_h_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::beta_h_>(h,P,res); }
  void _mu_h_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::mu_h_>(h,P,res); }
  void _lambda_h_(const double h, const double P, SpanD res) const { double_to_span<&Fluide_reel_base::lambda_h_>(h,P,res); }

  void _compute_CPMLB_pb_multiphase_h_(MLoiSpanD_h ) const;
  void _compute_all_pb_multiphase_h_(MLoiSpanD_h , MLoiSpanD_h ) const;
};

#endif /* Fluide_reel_base_included */
