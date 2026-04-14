#ifndef Op_Conv_Rusanov_Coloc_Vect_included
#define Op_Conv_Rusanov_Coloc_Vect_included


#include <Op_Conv_Coloc_base.h>

class Op_Conv_Rusanov_Coloc_Vect : public Op_Conv_Coloc_base_Vect
{
  Declare_instanciable( Op_Conv_Rusanov_Coloc_Vect ) ;
public:
  inline void scheme(DoubleTab&, const int&) const override;
};

#endif /*Op_Conv_Rusanov_Coloc_Vect_included*/

