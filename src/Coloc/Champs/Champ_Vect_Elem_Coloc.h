#ifndef Champ_Vect_Elem_Coloc_included
#define Champ_Vect_Elem_Coloc_included

#include <Champ_Elem_Coloc.h>
//#include <Domaine_Coloc.h>
//#include <TRUST_Ref.h>

class Champ_Vect_Elem_Coloc: public Champ_Elem_Coloc
{
  Declare_instanciable(Champ_Vect_Elem_Coloc);

protected:
  void init_fcl() const override;
//  const Domaine_Coloc& domaine_Coloc() const; //ok
//  int nb_valeurs_nodales() const override ; //ok
//  int fixer_nb_valeurs_nodales(int n) override; //ok
//
//
//  Champ_base& affecter_(const Champ_base& ch) override {return Champ_Inc_P0_base::affecter_(ch);}
//  // DoubleTab& valeur_aux_faces(DoubleTab& dst) const override; // a optm
//  //int imprime(Sortie& os, int ncomp) const override; //aux sommets mais y a imprimer_P0 // je peux l enlever pour coloc
//  int reprendre(Entree& fich) override; // non

};

#endif /* Champ_Vect_Elem_Coloc_included */

