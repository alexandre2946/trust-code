#include <Op_NConserv_negligeable.h>
#include <Champ_base.h>

Implemente_instanciable(Op_NConserv_negligeable,"Op_NConserv_negligeable",Operateur_NConserv_base);

Sortie& Op_NConserv_negligeable::printOn(Sortie& os) const { return os; }

Entree& Op_NConserv_negligeable::readOn(Entree& is) { return is; }

/*! @brief Associe la vitesse a l'operateur.
 *
 * @param (Champ_Inc_base& ch) le champ inconnue representant la vitesse
 */
//void Op_NConserv_negligeable::associer_vitesse(const Champ_base& ch)
//{
//  la_vitesse = ch;
//}

/*! @brief Renvoie le champ inconnue representant la vitesse
 *
 * @return (Champ_Inc_base&) le champ inconnue representant la vitesse
 */
const Champ_base& Op_NConserv_negligeable::vitesse() const
{
  return la_vitesse.valeur();
}

void Op_NConserv_negligeable::ajouter_flux(const DoubleTab& inconnue, DoubleTab& contribution) const
{
  // nothing to do
}

void Op_NConserv_negligeable::calculer_flux(const DoubleTab& inconnue, DoubleTab& flux) const
{
  flux = 0.0;
}
