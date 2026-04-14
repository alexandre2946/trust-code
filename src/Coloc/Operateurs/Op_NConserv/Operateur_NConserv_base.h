
#ifndef Operateur_NConserv_base_included
#define Operateur_NConserv_base_included

#include <Operateur_base.h>
#include <TRUST_Ref.h>

class Champ_base;
class Operateur_NConserv_base  : public Operateur_base
{
  Declare_base(Operateur_NConserv_base);
public :
  inline DoubleTab& calculer(const DoubleTab& inco, DoubleTab& resu) const override;
  void dimensionner(Matrice_Morse& ) const override {} ; //a enelver TODO
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override {};

  int has_interface_blocs() const override {  return 1; }

//protected:
//  mutable SFichier Flux_NConserv, Flux_NConserv_moment, Flux_NConserv_sum;


};

inline DoubleTab& Operateur_NConserv_base::calculer(const DoubleTab& inco, DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(inco, resu);
}

#endif
