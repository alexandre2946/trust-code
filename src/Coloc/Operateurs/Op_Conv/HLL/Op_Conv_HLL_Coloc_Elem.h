#ifndef Op_Conv_HLL_Coloc_Elem_included
#define Op_Conv_HLL_Coloc_Elem_included

#include <Op_Conv_Coloc_base.h>

class Op_Conv_HLL_Coloc_Elem : public Op_Conv_Coloc_base_Elem
{
  Declare_instanciable( Op_Conv_HLL_Coloc_Elem ) ;
public:
  inline void scheme(DoubleTab& num_flux, const int& f) const override ;
};

#endif /*Op_Conv_HLL_Coloc_Elem_included*/

