
#ifndef Milieu_composite_Euler_included
#define Milieu_composite_Euler_included

#include <Saturation_base.h>
#include <Interface_base.h>
#include <Interface_Baer_Nunziato.h>
#include <TRUST_Deriv.h>
#include <Fluide_base.h>
#include <vector>
#include <set>

class Milieu_composite_Euler: public Fluide_base
{

  Declare_instanciable(Milieu_composite_Euler);
public :
  int check_unknown_range() const override; // ok ( mais comprendre celle de chaque phase
  int initialiser(const double temps) override;

  bool initTimeStep(double dt) override;
  bool has_saturation(int k, int l) const;
  bool has_interface(int k, int l) const;
  inline bool has_saturation() const { return has_saturation_; }
  inline bool has_interface() const { return has_interface_; }

  void discretiser(const Probleme_base& pb, const  Discretisation_base& dis) override;
  void abortTimeStep() override;
  void preparer_calcul() override;
  void mettre_a_jour(double temps) override;
  void associer_equation(const Equation_base* eqn) const override;

  Interface_base& get_interface(int k, int l) const;
  Saturation_base& get_saturation(int k, int l) const;

  const Fluide_base& get_fluid(const int i) const;

  Fluide_base& get_fluid(const int i);

  inline const Noms& noms_phases() const { return noms_phases_; }

  inline bool are_fluid_properties_initialised() const { return fluid_properties_initialised_; }

  const Interface_base& interface_phase() const { return inter_lu_.valeur(); }
  Interface_base& interface_phase() { return inter_lu_.valeur(); }

  void init_energie_tot(DoubleTab& energie_tot_jdd) const;
  void calculer_pression(DoubleTab& pression_old) const ;
  void calculer_vitesse_son(DoubleTab& c) const;

protected :
  OWN_PTR(Champ_Don_base) rho_m_, h_m_;
  Noms noms_phases_;
  double t_init_ = -1.;
  bool has_saturation_ = false, has_interface_ = false;
  bool res_en_T_ = true; // par defaut resolution en T
  bool fluid_properties_initialised_ = false;
  std::vector<std::vector<Interface_base *>> tab_interface_;
  std::vector<OWN_PTR(Fluide_base)> fluides_;
  OWN_PTR(Interface_base) sat_lu_, inter_lu_;
  //OWN_PTR(Interface_Baer_Nunziato) inter_BN_;

  std::pair<std::string, int> check_fluid_name(const Nom& name);
  void mettre_a_jour_tabs_Euler();
  //void mettre_a_jour_tabs();
  static void calculer_masse_volumique(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
  static void calculer_energie_interne(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
  static void calculer_enthalpie(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
  static void calculer_temperature_multiphase(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv);
};

#endif /* Milieu_composite_Euler_included */
