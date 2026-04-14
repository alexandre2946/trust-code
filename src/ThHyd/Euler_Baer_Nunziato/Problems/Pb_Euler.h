

#ifndef Pb_Euler_included
#define Pb_Euler_included

#include <Density_Euler.h>
#include <Momentum_Euler.h>
#include <Energy_Euler.h>
#include <Fraction_Euler.h>
#include <Correlation_base.h>
#include <Pb_Fluide_base.h>
#include <Nom.h>
#include <TRUST_Deriv.h>
#include <TRUST_List.h>
#include <Interprete.h>
#include <Verif_Cl.h>


class Pb_Euler : public Pb_Fluide_base
{
  Declare_instanciable(Pb_Euler);
public:
  void associer_milieu_base(const Milieu_base& ) override; //ok
  void typer_lire_milieu(Entree& is) override; // ok
  int nombre_d_equations() const override {return 4;}; // ok
  int nb_phases() const {return noms_phases_.size();} ; //ok
  // pas d heritage mais ...... ok
  virtual Momentum_Euler& equation_qdm() { return eq_qdm_; }
  virtual const Momentum_Euler& equation_qdm() const { return eq_qdm_; } //ok
  virtual Density_Euler& equation_masse() { return eq_masse_; }//ok
  virtual const Density_Euler& equation_masse() const { return eq_masse_; }//ok
  virtual Energy_Euler& equation_energie() { return eq_energie_; }//ok
  virtual const Energy_Euler& equation_energie() const { return eq_energie_; }//ok
  virtual Fraction_Euler& equation_fraction() { return eq_fraction_; }//ok
  virtual const Fraction_Euler& equation_fraction() const { return eq_fraction_; }//ok
  const Nom& nom_phase(int i) const { return noms_phases_[i]; }  //ok
  const Noms& noms_phases() const { return noms_phases_; } //ok
  Entree& lire_equations(Entree& is, Motcle& dernier_mot) override; // ok
  const Equation_base& equation(int) const override ; // ok
  Equation_base& equation(int) override; //
  void preparer_calcul() override; //ok
  void mettre_a_jour(double temps) override;
//////////////////////////////////////////////////////////////////////////////////////////////////
  int verifier() override {return 1;}; // A modifier !!!!!!
  inline virtual bool resolution_en_T() const { return true; } // a anelever
/////////////////////////////////////////////////////////////////////////////////////////////////////

protected:
  Noms noms_phases_;
  Density_Euler eq_masse_;
  Momentum_Euler eq_qdm_;
  Energy_Euler eq_energie_;
  Fraction_Euler eq_fraction_;
};

#endif /* Pb_Euler_included */
