#ifndef Fraction_Euler_included
#define Fraction_Euler_included

#include <Conservation_Euler.h>

class Fraction_Euler : public Conservation_Euler
{
  Declare_instanciable(Fraction_Euler);
public :
  void discretiser() override; //ok
  int verif_Cl() const override {return 1;}; // TODO
  void set_param(Param& param) override;
  int nombre_d_operateurs() const override { return 1; }
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;
  inline double termes_NonConservatif(const double& alpha_bord, const double& vitesse_normale_interieur, const double& p_inter) const override
  {return  alpha_bord * vitesse_normale_interieur;};
};

#endif /* Fraction_Euler_included */
