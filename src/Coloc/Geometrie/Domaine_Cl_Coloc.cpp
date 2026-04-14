

#include <Dirichlet_entree_fluide_leaves.h>
#include <Champ_front_softanalytique.h>
#include <Dirichlet_paroi_defilante.h>
#include <Champ_Face_PolyMAC_P0P1NC.h>
#include <Dirichlet_paroi_fixe.h>
#include <Discretisation_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Face_PolyMAC.h>
#include <Dirichlet_homogene.h>
#include <Champ_Inc_P0_base.h>
#include <Domaine_PolyMAC.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <Periodique.h>
#include <Symetrie.h>
#include <Debog.h>

Implemente_instanciable(Domaine_Cl_Coloc, "Domaine_Cl_Coloc", Domaine_Cl_PolyMAC);

Sortie& Domaine_Cl_Coloc::printOn(Sortie& os) const { return os; }

Entree& Domaine_Cl_Coloc::readOn(Entree& is) { return Domaine_Cl_PolyMAC::readOn(is); }



void Domaine_Cl_Coloc::imposer_cond_lim(Champ_Inc_base& ch, double temps)
{

  DoubleTab& ch_tab = ch.valeurs(temps);

  if ( (ch.nature_du_champ() == scalaire) || (sub_type(Champ_Inc_P0_base, ch))) {}

  else
    {
      Cerr << "Le type de OWN_PTR(Champ_Inc_base) " << ch.que_suis_je() << " n'est pas prevu en Coloc family " << finl;
      Process::exit();
    }

  ch_tab.echange_espace_virtuel();
  Debog::verifier("Domaine_Cl_Coloc::imposer_cond_lim ch_tab", ch_tab);
}


