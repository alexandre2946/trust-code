#include <Op_Conv_Rusanov_Coloc_Elem.h>
#include <Domaine_Coloc.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_P0_base.h>
#include <Conservation_Euler.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
#include <Milieu_composite_Euler.h>
#include <Fluide_reel_base.h>


Implemente_instanciable(Op_Conv_Rusanov_Coloc_Elem,"Op_Conv_Rusanov_Coloc_Elem",Op_Conv_Coloc_base_Elem);

Sortie& Op_Conv_Rusanov_Coloc_Elem::printOn(Sortie& os) const { return Op_Conv_Coloc_base::printOn(os); }
Entree& Op_Conv_Rusanov_Coloc_Elem::readOn(Entree& is) {  Op_Conv_Coloc_base::readOn(is); return is;}


inline void Op_Conv_Rusanov_Coloc_Elem::scheme(DoubleTab& num_flux, const int& f) const
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

  for (int n =0 ; n< nb_phases; n++)
    {
      double un_l = vit_n(f,n), un_r = vit_n(f, n + nb_phases), c_l = c(el,n), c_r = c(er,n);
      double s = std::max(fabs(un_l-c_l),fabs(un_l+c_l));
      s = std::max(s,std::max(fabs(un_r-c_r),fabs(un_r+c_r)));
      num_flux(f,n) = 0.5 * ( flux_l(n) + flux_r(n) ) - s *0.5 * (w(er,n)-w(el,n));
    }
}
