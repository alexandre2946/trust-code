#include <Op_NConserv_Coloc_base.h>
#include <Domaine_Coloc.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_P0_base.h>
#include <Conservation_Euler.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
#include <Milieu_composite_Euler.h>
#include <Fluide_reel_base.h>


Implemente_base(Op_NConserv_Coloc_base,"Op_NConserv_Coloc_base",Operateur_NConserv_base);
Implemente_instanciable(Op_NConserv_Coloc_base_Elem,"Op_NConserv_Coloc_base_Elem",Op_NConserv_Coloc_base);
Implemente_instanciable(Op_NConserv_Coloc_base_Vect,"Op_NConserv_Coloc_base_Vect",Op_NConserv_Coloc_base);

Sortie& Op_NConserv_Coloc_base::printOn(Sortie& os) const { return Operateur_NConserv_base::printOn(os); }
Entree& Op_NConserv_Coloc_base::readOn(Entree& is) { Operateur_NConserv_base::readOn(is);  return is; }
Sortie& Op_NConserv_Coloc_base_Elem::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_Coloc_base_Elem::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is);  return is; }
Sortie& Op_NConserv_Coloc_base_Vect::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_Coloc_base_Vect::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is);  return is; }
void Op_NConserv_Coloc_base::completer()
{
  Operateur_base::completer();
  assert(le_dom_poly_.non_nul());
}

void Op_NConserv_Coloc_base::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& inc)
{
  le_dom_poly_ = ref_cast(Domaine_Coloc, domaine_dis);
  la_zcl_poly_ = ref_cast(Domaine_Cl_Coloc, zcl);
  le_champ_inco = ref_cast(Champ_Inc_base,inc);
}

void Op_NConserv_Coloc_base::associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl)
{
  la_zcl_poly_ = ref_cast(Domaine_Cl_Coloc, zcl);
}




