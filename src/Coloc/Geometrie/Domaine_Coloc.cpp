
#include <Domaine_Coloc.h>

Implemente_instanciable(Domaine_Coloc, "Domaine_Coloc", Domaine_PolyMAC_P0);

Sortie& Domaine_Coloc::printOn(Sortie& os) const { return Domaine_PolyMAC_P0::printOn(os); }

Entree& Domaine_Coloc::readOn(Entree& is) { return Domaine_PolyMAC_P0::readOn(is); }

