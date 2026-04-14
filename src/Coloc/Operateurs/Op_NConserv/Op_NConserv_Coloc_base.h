#ifndef Op_NConserv_Coloc_base_included
#define Op_NConserv_Coloc_base_included

#include <Operateur_NConserv_base.h>
#include <Domaine_Coloc.h>
#include <TRUST_Ref.h>
#include <SFichier.h>

class Domaine_Cl_PolyMAC;
class Domaine_PolyMAC;

class Op_NConserv_Coloc_base : public Operateur_NConserv_base
{
  Declare_base(Op_NConserv_Coloc_base) ;
public:
  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;
  void associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl) override;
  void mettre_a_jour(double temps) override {};
  int has_interface_blocs() const override {  return 1; }
  void completer() override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override {Process::exit();};
  virtual void Abgral_scheme(DoubleTab& num_flux_left, DoubleTab& num_flux_right) const {Process::exit();};
protected:
  OBS_PTR(Domaine_PolyMAC) le_dom_poly_;
  OBS_PTR(Domaine_Cl_PolyMAC) la_zcl_poly_;

};


class Op_NConserv_Coloc_base_Elem : public Op_NConserv_Coloc_base
{
  Declare_instanciable(Op_NConserv_Coloc_base_Elem) ;
public:
  //void Riemann_solver(DoubleTab& num_flux) const override;
};


class Op_NConserv_Coloc_base_Vect : public Op_NConserv_Coloc_base
{
  Declare_instanciable(Op_NConserv_Coloc_base_Vect) ;
public:
  //void Riemann_solver(DoubleTab& num_flux) const override;
};



#endif /*Op_NConserv_Coloc_base_included*/

