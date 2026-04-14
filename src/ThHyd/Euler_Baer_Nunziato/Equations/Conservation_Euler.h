#ifndef Conservation_Euler_included
#define Conservation_Euler_included

#include <Convection_Diffusion_std.h>
//#include <Operateur_Evanescence.h>
//#include <Operateur_Grad.h>
#include <TRUST_Ref.h>
#include <Domaine_Coloc.h>
#include <Operateur_NConserv.h>

class Fluide_base;

class Conservation_Euler : public Convection_Diffusion_std
{
  Declare_instanciable(Conservation_Euler);
public :
  int nombre_d_operateurs() const override { Process::exit();  return 1; } //ok
  const Operateur& operateur(int) const override { Process::exit();  return terme_convectif;} //ok
  Operateur& operateur(int) override { Process::exit();  return terme_convectif;} //ok
  void completer() override {Equation_base::completer();} ;
  inline const Champ_Inc_base& inconnue() const override { return l_inco_ch_; };
  inline Champ_Inc_base& inconnue() override { return l_inco_ch_; };
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  double calculer_pas_de_temps() const override {return 1e10;} ;

  void associer_milieu_base(const Milieu_base& ) override;
  void associer_fluide(const Fluide_base& );
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  const Fluide_base& fluide() const;
  Fluide_base& fluide();

  virtual inline DoubleTab flux(const int& f, const int& elem)  const {return 0; Process::exit();}
  virtual inline double flux_bord(const double& inco_bord, const double& vit_n_bord, const double& p_bord) const {return 0; Process::exit();};
  virtual inline double termes_NonConservatif(const double& alpha_bord, const double& vitesse_n_inter, const double& p_bord) const {return 0; Process::exit();};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // TODO
  void dimensionner_matrice_sans_mem(Matrice_Morse& matrice) override {  Process::exit();} ;
  int has_interface_blocs() const override {Process::exit(); return -1;};
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override {Process::exit();};
  void assembler_blocs_avec_inertie(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) override {Process::exit();};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

protected :
  Operateur_NConserv terme_nconserv_;
  OWN_PTR(Champ_Inc_base) l_inco_ch_;
  OBS_PTR(Fluide_base) le_fluide_;
};

#endif /* Conservation_Euler_included */
