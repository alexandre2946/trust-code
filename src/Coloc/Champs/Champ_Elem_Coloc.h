#ifndef Champ_Elem_Coloc_included
#define Champ_Elem_Coloc_included


#include <Champ_Inc_P0_base.h>
#include <Domaine_Coloc.h>

//#include <Operateur.h>
//#include <Op_Diff_PolyMAC_base.h>

class Champ_Elem_Coloc: public Champ_Inc_P0_base
{
  Declare_instanciable(Champ_Elem_Coloc);
public:
  const Domaine_Coloc& domaine_Coloc() const; //ok
  //Champ_base& affecter_(const Champ_base& ch) override {return Champ_Inc_P0_base::affecter_(ch);} //ok
  int fixer_nb_valeurs_nodales(int n) override; //ok
  int nb_valeurs_nodales() const override ; //ok

  // DoubleTab& valeur_aux_faces(DoubleTab& dst) const override; // mettre exit , parce que pas besoin
  //int imprime(Sortie& os, int ncomp) const override; // aux sommets mais y a imprimer_P0 // je peux l enlever pour coloc
  int reprendre(Entree& fich) override;
};


#endif /* Champ_Elem_Coloc_included */


//TODO

