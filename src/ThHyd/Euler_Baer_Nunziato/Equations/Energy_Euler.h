#ifndef Energy_Euler_included
#define Energy_Euler_included

#include <Conservation_Euler.h>
#include <Milieu_composite_Euler.h>
class Energy_Euler : public Conservation_Euler
{
  Declare_instanciable(Energy_Euler);
public :
  void discretiser() override; //ok
  int verif_Cl() const override {return 1;}; // TODO
  void set_param(Param& param) override;
  int nombre_d_operateurs() const override { return 2; }
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  Entree& lire_cond_init(Entree& is) override;
  inline double flux_bord(const double& rhoE_bord, const double& vit_n_bord, const double& p_bord) const override;
  inline DoubleTab flux(const int& f, const int& left_or_right) const override;
  inline double termes_NonConservatif(const double& alpha_bord, const double& vitesse_normale_interieur, const double& p_inter) const override
  {return -alpha_bord * vitesse_normale_interieur * p_inter;};

  inline void init_energie_tot()
  {
    const Milieu_composite_Euler& mil=ref_cast(Milieu_composite_Euler,milieu());
    mil.init_energie_tot(inconnue().valeurs());
  };

};
#endif /* Energy_Euler_included */
