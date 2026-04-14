#ifndef Op_NConserv_negligeable_included
#define Op_NConserv_negligeable_included

#include <Operateur_negligeable.h>
#include <Operateur_NConserv_base.h>
#include <TRUST_Ref.h>

class Champ_base;

/*! @brief Classe Op_NConserv_negligeable Cette classe represente un opperateur de convection negligeable.
 *
 *     Lorsqu'un operateur de ce type est utilise dans une equation
 *     cela revient a negliger le terme de convection.
 *     Les methodes de modification et de participation a un calcul de
 *     l'operateur sont en fait des appels aux meme methodes de
 *     Operateur_negligeable qui ne font rien.
 *
 * @sa Operateur_negligeable Operateur_Conv_base
 */
class Op_NConserv_negligeable: public Operateur_negligeable,
  public Operateur_NConserv_base
{
  Declare_instanciable(Op_NConserv_negligeable);

public :

  inline void contribuer_au_second_membre(DoubleTab& ) const override { }
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override { }
  inline void associer_domaine_cl_dis(const Domaine_Cl_dis_base&) override { }
  /* interface {dimensionner,ajouter}_blocs -> ne font rien */
  int  has_interface_blocs() const override { return 1; }
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override { }
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override { }

  inline void mettre_a_jour(double) override;
  //void associer_vitesse(const Champ_base& ) override ;
  const Champ_base& vitesse() const;

  void ajouter_flux(const DoubleTab& inconnue, DoubleTab& contribution) const override;
  void calculer_flux(const DoubleTab& inconnue, DoubleTab& flux) const override;

  void check_multiphase_compatibility() const override { };
  //void set_incompressible(const int) override { };

protected :

  OBS_PTR(Champ_base) la_vitesse;
  inline void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base& ) override ;
};

/*! @brief Mise a jour en temps d'un operateur negligeable: NE FAIT RIEN Simple appel a Operateur_negligeable::mettre_a_jour(double)
 *
 * @param (double temps)
 */
inline void Op_NConserv_negligeable::mettre_a_jour(double temps)
{
  Operateur_negligeable::mettre_a_jour(temps);
}


/*! @brief Associe divers objets a un operateurs negligeable: NE FAIT RIEN Simple appel a Operateur_negligeable::associer(const Domaine_dis_base&,
 *
 *                                                      const Domaine_Cl_dis_base&,
 *                                                      const Champ_Inc_base&)
 *
 * @param (Domaine_dis_base& z)
 * @param (Domaine_Cl_dis_base& zcl)
 * @param (Champ_Inc_base& ch)
 */
inline void Op_NConserv_negligeable::associer(const Domaine_dis_base& z,
                                              const Domaine_Cl_dis_base& zcl,
                                              const Champ_Inc_base& ch)
{
  Operateur_negligeable::associer(z, zcl, ch);
}
#endif
