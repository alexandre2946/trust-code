#ifndef Op_NConserv_HLL_Coloc_Vect_included
#define Op_NConserv_HLL_Coloc_Vect_included

#include <Op_NConserv_Coloc_base.h>

class Op_NConserv_HLL_Coloc_Vect : public Op_NConserv_Coloc_base_Vect
{
  Declare_instanciable( Op_NConserv_HLL_Coloc_Vect ) ;
public:
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;
  void Abgral_scheme(DoubleTab& num_flux_left, DoubleTab& num_flux_right) const override;
};

#endif /*Op_NConserv_HLL_Coloc_Vect_included*/

