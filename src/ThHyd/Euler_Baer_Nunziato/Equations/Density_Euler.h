#ifndef Density_Euler_included
#define Density_Euler_included

#include <Conservation_Euler.h>

class Density_Euler : public Conservation_Euler
{
  Declare_instanciable(Density_Euler);
public :
  void discretiser() override; //OK
  Entree& lire_cond_init(Entree& is) override;
  inline DoubleTab flux(const int& f, const int& e ) const override;
  inline double flux_bord(const double& alpha_rho_bord, const double& vit_n_bord, const double& p_bord ) const override;
  int verif_Cl() const override {return 1;};
  const Champ_Inc_base& densite() const { return densite_.valeur();};
  Champ_Inc_base& densite() { return densite_.valeur();}
  void mettre_a_jour_champs_conserves(double temps, int reset) override;
  void init_alpha_rho();
  void set_param(Param& param) override;
  int nombre_d_operateurs() const override { return 1;}
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
protected:
  OWN_PTR(Champ_Inc_base) densite_;
};
#endif /* Density_Euler_included */
