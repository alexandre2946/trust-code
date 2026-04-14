
#ifndef Domaine_Cl_Coloc_included
#define Domaine_Cl_Coloc_included

//  Cette classe porte les tableaux qui servent a mettre en oeuvre
//  les condition aux limites dans la formulation Coloc
#include <Domaine_Cl_PolyMAC.h>



class Domaine_VF;

class Domaine_Cl_Coloc : public Domaine_Cl_PolyMAC
{

  Declare_instanciable(Domaine_Cl_Coloc);

public :

  //void completer(const Domaine_dis_base& ) override;
  //int initialiser(double temps) override;
  void imposer_cond_lim(Champ_Inc_base&, double) override;

  //int nb_faces_sortie_libre() const;

  //Domaine_VF& domaine_vf();
  //const Domaine_VF& domaine_vf() const;

  //int nb_bord_periodicite() const;
//protected:

  //int modif_perio_fait_ = 0;
};

#endif /* Domaine_Cl_Coloc_included */
