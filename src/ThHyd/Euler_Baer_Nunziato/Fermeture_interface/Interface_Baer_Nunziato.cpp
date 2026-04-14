#include <Interface_Baer_Nunziato.h>
#include <Pb_Euler.h>

Implemente_instanciable(Interface_Baer_Nunziato, "Interface_Baer_Nunziato", Interface_base);
// XD Interface_Baer_Nunziato Interface_base Interface_Baer_Nunziato -1 Liquid-gas interface with a constant surface tension sigma

Sortie& Interface_Baer_Nunziato::printOn(Sortie& os) const { return os; }
Entree& Interface_Baer_Nunziato::readOn(Entree& is)
{

  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Interface_Baer_Nunziato::set_param(Param& param)
{
  param.ajouter("vitesse", &id_vitesse_interface_ , Param::REQUIRED);
  param.ajouter("pression", &id_pression_interface_, Param::REQUIRED);
}

