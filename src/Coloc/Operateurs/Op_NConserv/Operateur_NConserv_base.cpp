
#include <Operateur_NConserv_base.h>
#include <TRUSTTrav.h>

Implemente_base(Operateur_NConserv_base,"Operateur_NConserv_base",Operateur_base);


Sortie& Operateur_NConserv_base::printOn(Sortie& os) const
{
  return os;
}

Entree& Operateur_NConserv_base::readOn(Entree& is)
{
  return is;
}




