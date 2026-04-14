#include <Op_Conv_HLL_Coloc_Elem.h>

#include <Domaine_Coloc.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_P0_base.h>
#include <Conservation_Euler.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
#include <Milieu_composite_Euler.h>
#include <Fluide_reel_base.h>


Implemente_instanciable(Op_Conv_HLL_Coloc_Elem,"Op_Conv_HLL_Coloc_Elem",Op_Conv_Coloc_base_Elem);

Sortie& Op_Conv_HLL_Coloc_Elem::printOn(Sortie& os) const { return Op_Conv_Coloc_base::printOn(os); }
Entree& Op_Conv_HLL_Coloc_Elem::readOn(Entree& is) {  Op_Conv_Coloc_base::readOn(is); return is;}



inline void Op_Conv_HLL_Coloc_Elem::scheme(DoubleTab& num_flux, const int& f) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const Conservation_Euler& eq = ref_cast(Conservation_Euler,equation());
  const IntTab& f_e = domaine.face_voisins();
  const DoubleTab& vit_n= ref_cast(Momentum_Euler, equation().probleme().equation(0)).vitesse_normale();
  const int nb_phases = ref_cast(Pb_Euler,eq.probleme()).nb_phases();
  const DoubleTab& c = ref_cast(Momentum_Euler,equation().probleme().equation(0)).vitesse_son();
  const DoubleTab& w = le_champ_inco->valeurs();
  const int el = f_e(f,0), er = f_e(f,1);
  const DoubleTab flux_l = eq.flux(f,0), flux_r = eq.flux(f,1);

  double Sp=0, Sm=0;

  for (int n = 0; n < nb_phases; n++)
    {
      double un_l = vit_n(f,n), un_r = vit_n(f, n + nb_phases), c_l = c(el,n), c_r = c(er,n);
      double Sm_k = std::min(un_l - c_l, un_r - c_r);
      double Sp_k = std::max(un_l + c_l, un_r + c_r);
      Sm = std::min(Sm,std::min(0.0, Sm_k));
      Sp = std::max(Sp,std::max(0.0, Sp_k));

    }

  for (int n = 0; n < nb_phases; n++)
    {
      // double un_l = vit_n(f,n), un_r = vit_n(f, n + nb_phases), c_l = c(el,n), c_r = c(er,n);
//      double Sm = std::min(un_l - c_l, un_r - c_r);
//      double Sp = std::max(un_l + c_l, un_r + c_r);
//      Sm = std::min(0.0, Sm);
//      Sp = std::max(0.0, Sp);



      num_flux(f, n) = (Sp * flux_l(n) - Sm * flux_r(n) + Sp * Sm * (w(er, n) - w(el, n)));
      num_flux(f, n) /= (Sp - Sm);

    }
}
