#include <Sortie_supersonique.h>

Implemente_instanciable(Sortie_supersonique,"Sortie_supersonique",Symetrie);

Sortie& Sortie_supersonique::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

Entree& Sortie_supersonique::readOn(Entree& s )
{
  le_champ_front.typer("Champ_front_uniforme");
  le_champ_front->fixer_nb_comp(0);
  return s ;
}

