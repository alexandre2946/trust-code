

#ifndef Coloc_discretisation_included
#define Coloc_discretisation_included

#include <PolyMAC_P0_discretisation.h>

class Coloc_discretisation: public PolyMAC_P0_discretisation
{
  Declare_instanciable(Coloc_discretisation);
public:
  using PolyMAC_P0_discretisation::discretiser_champ;
  bool is_coloc() const override { return true; }
  bool is_polymac_p0() const override { return false; }
  Nom domaine_cl_dis_type() const override { return "Domaine_Cl_Coloc"; }

  void discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& nom, const Noms& unite, int nb_comp, int nb_pas_dt, double temps, OWN_PTR(Champ_Inc_base)& champ,
                         const Nom& sous_type = NOM_VIDE) const override;
  Nom get_name_of_type_for(const Nom& class_operateur, const Nom& type_operateur, const Equation_base& eqn, const OBS_PTR(Champ_base) &champ_sup) const override;
  //void pression(const Schema_Temps_base& sch, Domaine_dis_base& z, OWN_PTR(Champ_Inc_base) &ch) const override;

private:
  void discretiser_champ_fonc_don(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& nom, const Noms& unite, int nb_comp, double temps, Objet_U& champ) const override ;
};

#endif /* Coloc_discretisation_included */


//TODO

// pression n est pas un inc ; enlever de discretiser_champ

