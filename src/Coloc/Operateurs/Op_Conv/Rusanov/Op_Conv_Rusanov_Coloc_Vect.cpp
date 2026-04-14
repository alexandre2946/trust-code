#include <Op_Conv_Rusanov_Coloc_Vect.h>
#include <Domaine_Coloc.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_P0_base.h>
#include <Conservation_Euler.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
#include <Milieu_composite_Euler.h>
#include <Fluide_reel_base.h>


Implemente_instanciable(Op_Conv_Rusanov_Coloc_Vect,"Op_Conv_Rusanov_Coloc_Vect",Op_Conv_Coloc_base_Vect);

Sortie& Op_Conv_Rusanov_Coloc_Vect::printOn(Sortie& os) const { return Op_Conv_Coloc_base::printOn(os); }
Entree& Op_Conv_Rusanov_Coloc_Vect::readOn(Entree& is) { Op_Conv_Coloc_base::readOn(is); return is;}

inline void Op_Conv_Rusanov_Coloc_Vect::scheme(DoubleTab& num_flux, const int& f) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const DoubleTab& w = le_champ_inco->valeurs();
  const IntTab& f_e = domaine.face_voisins();
  const Momentum_Euler& eq = ref_cast(Momentum_Euler,equation());
  const DoubleTab& vit_n= eq.vitesse_normale();
  const DoubleTab& c = eq.vitesse_son();
  const int nb_phase = ref_cast(Pb_Euler,eq.probleme()).nb_phases();
  const int el = f_e(f,0), er = f_e(f,1);
  const DoubleTab flux_l = eq.flux_(f,0), flux_r = eq.flux_(f,1);

  for (int n =0 ; n< nb_phase; n++)
    {
      double un_l = vit_n(f,n), un_r = vit_n(f, n + nb_phase), c_l = c(el,n), c_r = c(er,n);
      double s = std::max(fabs(un_l-c_l),fabs(un_l+c_l));
      s = std::max(s,std::max(fabs(un_r-c_r),fabs(un_r+c_r)));
      for (int d = 0; d< dimension; d++)
        num_flux(f, n + nb_phase * d ) = 0.5* (flux_l(n + nb_phase * d) + flux_r(n + nb_phase*d)) - s *0.5 * (w(er,n + nb_phase * d)-w(el,n + nb_phase * d));
    }
}


