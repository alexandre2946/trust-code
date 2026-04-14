#include <Op_NConserv_HLL_Coloc_Elem.h>
#include <Domaine_Coloc.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_P0_base.h>
#include <Conservation_Euler.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
#include <Milieu_composite_Euler.h>
#include <Fluide_reel_base.h>
#include <Interface_Baer_Nunziato.h>
#include <Sortie_supersonique.h>
#include <Entree_supersonique.h>
#include <Neumann_paroi_flux_nul.h>
Implemente_instanciable(Op_NConserv_HLL_Coloc_Elem,"Op_NConserv_HLL_Coloc_Elem",Op_NConserv_Coloc_base_Elem);

Sortie& Op_NConserv_HLL_Coloc_Elem::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_HLL_Coloc_Elem::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is); return is;}

void Op_NConserv_HLL_Coloc_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();
  assert(secmem.line_size()==2);
  const int N = domaine.nb_faces_tot();
  DoubleTab num_flux_left(N);
  DoubleTab num_flux_right(N);

  Abgral_scheme(num_flux_left,num_flux_right);

  for (int f = 0; f < N; f++)
    {
      for (int i = 0; i < 2; i++)
        {
          int e = f_e(f, i);
          if ( e >= 0 && e < domaine.nb_elem())
            {
              secmem(e,0) -= (i ? num_flux_right(f) : num_flux_left(f)) * fs(f);
              secmem(e,1) += (i ? num_flux_right(f) : num_flux_left(f)) * fs(f); //
            }
        }
    }
}

void Op_NConserv_HLL_Coloc_Elem::Abgral_scheme(DoubleTab& num_flux_left,DoubleTab& num_flux_right) const
{

  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const IntTab& f_e = domaine.face_voisins();
  const IntTab&  fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();

  const Pb_Euler& pb = ref_cast(Pb_Euler,equation().probleme());
  const Conservation_Euler& eq = ref_cast(Conservation_Euler,equation());

  const Interface_Baer_Nunziato& interface = ref_cast(Interface_Baer_Nunziato,
                                                      ref_cast(Milieu_composite_Euler,pb.milieu()).interface_phase());
  const int n = interface.id_phase_vitesse_inter();
  const int m = interface.id_phase_pression_inter();

  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  const DoubleTab& vit_n= pb.equation_qdm().vitesse_normale();
  //const DoubleTab& c= pb.equation_qdm().vitesse_son();
  const DoubleTab& p = pb.equation_qdm().pression().valeurs();
  //const int nb_phases = pb.nb_phases();

  //const bool equation_alpha= (sub_type(Fraction_Euler, equation()));

  const Conds_lim& cls = equation().domaine_Cl_dis().les_conditions_limites();

  const Conds_lim& cls_alpha = pb.equation_fraction().domaine_Cl_dis().les_conditions_limites();
  //const Conds_lim& cls_energie = pb.equation_energie().domaine_Cl_dis().les_conditions_limites();

  for (int f = 0; f < domaine.nb_faces_tot(); f++)
    {
      if (fcl(f, 0) == 0 ) calculer_terme_NC(num_flux_left,num_flux_right,f);

      else
        {

          const int& e = f_e(f,0) >= 0 ?  f_e(f,0) :  f_e(f,1); //pas besoin
          assert(f_e(f,0) >= 0 && vit_n (f,0)!=123.123); //pas besoin
          num_flux_right(f) = 123.123; //pas besoin

          DoubleTab normal(dimension);
          for (int d = 0; d< dimension; d++) normal(d) = domaine.face_normales(f,d) / domaine.face_surfaces(f);

          if (sub_type(Sortie_supersonique,cls[fcl(f, 1)].valeur()))
            {
              num_flux_left(f) = eq.termes_NonConservatif(alpha(e,0),vit_n (f,n),p(e, m));
            }

          else if (sub_type(Entree_supersonique,cls[fcl(f, 1)].valeur())) // Dirichlet  : 6
            {
              const double alpha_bord = ref_cast(Dirichlet, cls_alpha[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), 0);
              num_flux_left(f) = eq.termes_NonConservatif(alpha_bord,vit_n(f,n),p(e, m));
            }

          else if ( sub_type(Neumann_paroi_flux_nul,cls[fcl(f, 1)].valeur())) //Neumann_homogene : 5
            {
              num_flux_left(f) = eq.termes_NonConservatif(alpha(e,0),vit_n (f,n),p(e, m));
            }
          else
            {
              Cerr<<"The BC of type " << fcl(f, 0) << "for the equation " << eq.que_suis_je() << " is not available \n";
              Process::exit();
            }
        }
    }
}


void Op_NConserv_HLL_Coloc_Elem::calculer_terme_NC(DoubleTab& num_flux_left,DoubleTab& num_flux_right, const int& f) const
{
  if (sub_type(Fraction_Euler, equation()))  calculer_terme_NC_fraction(num_flux_left,num_flux_right, f);
  else if (sub_type(Energy_Euler, equation())) calculer_terme_NC_energie(num_flux_left,num_flux_right, f);
  else Process::exit();
}

void Op_NConserv_HLL_Coloc_Elem::calculer_terme_NC_fraction(DoubleTab& num_flux_left,DoubleTab& num_flux_right, const int& f)  const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const Pb_Euler& pb = ref_cast(Pb_Euler,equation().probleme());
  const Interface_Baer_Nunziato& interface = ref_cast(Interface_Baer_Nunziato,
                                                      ref_cast(Milieu_composite_Euler,pb.milieu()).interface_phase());
  const int n = interface.id_phase_vitesse_inter();
  const int m = interface.id_phase_pression_inter();
  const IntTab& f_e = domaine.face_voisins();
  //const IntTab&  fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();
  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  //const Conservation_Euler& eq = ref_cast(Conservation_Euler,equation());
  const DoubleTab& vit_n= pb.equation_qdm().vitesse_normale();
  const DoubleTab& c= pb.equation_qdm().vitesse_son();

  const int nb_phases = pb.nb_phases();

  const int el = f_e(f,0), er = f_e(f,1);

  int k = m;
  double un_l = vit_n(f,k), un_r = vit_n(f, k + nb_phases), c_l = c(el,k), c_r = c(er,k);
  double Sm1 = std::min(un_l - c_l, un_r - c_r);
  double Sp1 = std::max(un_l + c_l, un_r + c_r);

  k=n;
  un_l = vit_n(f,k), un_r = vit_n(f, k + nb_phases), c_l = c(el,k), c_r = c(er,k);
  double Sm2 = std::min(un_l - c_l, un_r - c_r);
  double Sp2 = std::max(un_l + c_l, un_r + c_r);

  double Sm = std::min(0.0, std::min(Sm1,Sm2));
  double Sp = std::max(0.0, std::max(Sp1,Sp2));


  num_flux_left(f) = (Sp*alpha(el,0)-Sm*alpha(er,0)) * un_l + Sp*Sm*(alpha(er,0)-alpha(el,0));
  num_flux_left(f) /= (Sp - Sm);

  k = m;
  un_r = - vit_n(f,k), un_l = - vit_n(f, k + nb_phases), c_l = c(er,k), c_r = c(el,k);
  Sm1 = std::min(un_l - c_l, un_r - c_r);
  Sp1 = std::max(un_l + c_l, un_r + c_r);
  k=n;
  un_r = - vit_n(f,k), un_l = - vit_n(f, k + nb_phases), c_l = c(er,k), c_r = c(el,k);
  Sm2 = std::min(un_l - c_l, un_r - c_r);
  Sp2 = std::max(un_l + c_l, un_r + c_r);
  Sm = std::min(0.0, std::min(Sm1,Sm2));
  Sp = std::max(0.0, std::max(Sp1,Sp2));

  num_flux_right(f) = (Sp*alpha(er,0)-Sm*alpha(el,0)) * un_l + Sp*Sm*(alpha(el,0)-alpha(er,0));
  num_flux_right(f) /= (Sp - Sm);
}

void Op_NConserv_HLL_Coloc_Elem::calculer_terme_NC_energie(DoubleTab& num_flux_left,DoubleTab& num_flux_right, const int& f) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const Pb_Euler& pb = ref_cast(Pb_Euler,equation().probleme());
  const Interface_Baer_Nunziato& interface = ref_cast(Interface_Baer_Nunziato,
                                                      ref_cast(Milieu_composite_Euler,pb.milieu()).interface_phase());
  const int n = interface.id_phase_vitesse_inter();
  const int m = interface.id_phase_pression_inter();

  const IntTab& f_e = domaine.face_voisins();
  //const IntTab&  fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();
  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  //const Conservation_Euler& eq = ref_cast(Conservation_Euler,equation());
  const DoubleTab& vit_n= pb.equation_qdm().vitesse_normale();
  const DoubleTab& c= pb.equation_qdm().vitesse_son();
  const DoubleTab& p = pb.equation_qdm().pression().valeurs();
  const int nb_phases = pb.nb_phases();

  const int el = f_e(f,0), er = f_e(f,1);

  int k = m;
  double un_l = vit_n(f,k), un_r = vit_n(f, k + nb_phases), c_l = c(el,k), c_r = c(er,k);
  double Sm1 = std::min(un_l - c_l, un_r - c_r);
  double Sp1 = std::max(un_l + c_l, un_r + c_r);
  k=n;
  un_l = vit_n(f,k), un_r = vit_n(f, k + nb_phases), c_l = c(el,k), c_r = c(er,k);
  double Sm2 = std::min(un_l - c_l, un_r - c_r);
  double Sp2 = std::max(un_l + c_l, un_r + c_r);
  double Sm = std::min(0.0, std::min(Sm1,Sm2));
  double Sp = std::max(0.0, std::max(Sp1,Sp2));


  num_flux_left(f) = (Sp*alpha(el,0)-Sm*alpha(er,0)) * un_l * p(el,m);
  num_flux_left(f) /= -(Sp - Sm);

  k = m;
  un_r = - vit_n(f,k), un_l = - vit_n(f, k + nb_phases), c_l = c(er,k), c_r = c(el,k);
  Sm1 = std::min(un_l - c_l, un_r - c_r);
  Sp1 = std::max(un_l + c_l, un_r + c_r);
  k=n;
  un_r = - vit_n(f,k), un_l = - vit_n(f, k + nb_phases), c_l = c(er,k), c_r = c(el,k);
  Sm2 = std::min(un_l - c_l, un_r - c_r);
  Sp2 = std::max(un_l + c_l, un_r + c_r);
  Sm = std::min(0.0, std::min(Sm1,Sm2));
  Sp = std::max(0.0, std::max(Sp1,Sp2));

  num_flux_right(f) = (Sp*alpha(er,0)-Sm*alpha(el,0)) * un_l * p(er,m);
  num_flux_right(f) /= -(Sp - Sm);
}

