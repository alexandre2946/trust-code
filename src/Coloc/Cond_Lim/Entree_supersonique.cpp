#include <Entree_supersonique.h>

Implemente_instanciable(Entree_supersonique,"Entree_supersonique",Entree_fluide_alpha_impose);

Sortie& Entree_supersonique::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

Entree& Entree_supersonique::readOn(Entree& s )
{
  Entree_fluide_alpha_impose::readOn(s);
  return s ;
}

