#ifndef Op_NConserv_HLL_Coloc_Elem_included
#define Op_NConserv_HLL_Coloc_Elem_included

#include <Op_NConserv_Coloc_base.h>

class Op_NConserv_HLL_Coloc_Elem : public Op_NConserv_Coloc_base_Elem
{
  Declare_instanciable( Op_NConserv_HLL_Coloc_Elem ) ;
public:
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;
  void calculer_terme_NC(DoubleTab& num_flux_left, DoubleTab& num_flux_right, const int& f) const;

  void calculer_terme_NC_fraction(DoubleTab& num_flux_left, DoubleTab& num_flux_right, const int& f) const;
  void calculer_terme_NC_energie(DoubleTab& num_flux_left, DoubleTab& num_flux_right, const int& f) const;
  void Abgral_scheme(DoubleTab& , DoubleTab& ) const override;
};

#endif /*Op_NConserv_HLL_Coloc_Elem_included*/

