#include <Op_NConserv_HLL_Coloc_Vect.h>
#include <Domaine_Coloc.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_P0_base.h>
#include <Conservation_Euler.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
#include <Milieu_composite_Euler.h>
#include <Fluide_reel_base.h>
#include <Sortie_supersonique.h>

Implemente_instanciable(Op_NConserv_HLL_Coloc_Vect,"Op_NConserv_HLL_Coloc_Vect",Op_NConserv_Coloc_base_Vect);

Sortie& Op_NConserv_HLL_Coloc_Vect::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_HLL_Coloc_Vect::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is); return is;}

void Op_NConserv_HLL_Coloc_Vect::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();

  assert(secmem.line_size()/dimension==2);
  const int N = domaine.nb_faces_tot();
  DoubleTab num_flux_left(N,dimension);
  DoubleTab num_flux_right(N,dimension);

  Abgral_scheme(num_flux_left,num_flux_right);

  for (int f = 0; f < N; f++)
    {
      for (int i = 0; i < 2; i++)
        {
          int e = f_e(f, i);
          if ( e >= 0 && e < domaine.nb_elem())
            {
              double val = (i ? num_flux_right(f,0) : num_flux_left(f,0)) * fs(f);
              secmem(e,0) -= val;
              secmem(e,1) += val;
              val = (i ? num_flux_right(f,1) : num_flux_left(f,1)) * fs(f);;
              secmem(e,2) -= val;
              secmem(e,3) += val;//
            }
        }
    }

}

void Op_NConserv_HLL_Coloc_Vect::Abgral_scheme(DoubleTab& num_flux_left, DoubleTab& num_flux_right) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc,le_dom_poly_.valeur());
  const IntTab& f_e = domaine.face_voisins();
  const IntTab&  fcl = ref_cast(Champ_Inc_P0_base, equation().inconnue()).fcl();

  const Momentum_Euler& eq = ref_cast(Momentum_Euler,equation());
  const DoubleTab& vit_n= eq.vitesse_normale();
  //const DoubleTab& vit= eq.vitesse().valeurs();
  const DoubleTab& p = eq.pression().valeurs();

  const Pb_Euler& pb = ref_cast(Pb_Euler,equation().probleme());
  const Interface_Baer_Nunziato& interface = ref_cast(Interface_Baer_Nunziato,
                                                      ref_cast(Milieu_composite_Euler,pb.milieu()).interface_phase());
  const int n = interface.id_phase_vitesse_inter();
  const int m = interface.id_phase_pression_inter();
  const int nb_phase = pb.nb_phases();
  //const DoubleTab& rho = pb.equation_masse().densite().valeurs();
  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  const DoubleTab& c= pb.equation_qdm().vitesse_son();
  for (int f = 0; f < domaine.nb_faces_tot(); f++)
    {
      const int el = f_e(f,0), er = f_e(f,1);
      if (fcl(f, 0) == 0 )
        {
          // XXX
//          int k = n; // ce point n est pas claire : on utilise Sp et Sm de la phase pour les 2 phases
//          k=0;
//          double un_l = vit_n(f,k), un_r = vit_n(f, k + nb_phase), c_l = c(el,k), c_r = c(er,k);
//          double Sm = std::min(un_l - c_l, un_r - c_r);
//          double Sp = std::max(un_l + c_l, un_r + c_r);
//          Sm = std::min(0.0, Sm);
//          Sp = std::max(0.0, Sp);

          int k = m;
          double un_l = vit_n(f,k), un_r = vit_n(f, k + nb_phase), c_l = c(el,k), c_r = c(er,k);
          double Sm1 = std::min(un_l - c_l, un_r - c_r);
          double Sp1 = std::max(un_l + c_l, un_r + c_r);

          k=n;
          un_l = vit_n(f,k), un_r = vit_n(f, k + nb_phase), c_l = c(el,k), c_r = c(er,k);
          double Sm2 = std::min(un_l - c_l, un_r - c_r);
          double Sp2 = std::max(un_l + c_l, un_r + c_r);

          double Sm = std::min(0.0, std::min(Sm1,Sm2));
          double Sp = std::max(0.0, std::max(Sp1,Sp2));
          for (int d = 0; d< dimension; d++)
            {
              double n_d = domaine.face_normales(f,d)/domaine.face_surfaces(f);
              num_flux_left(f,d) = (Sp*alpha(el,0)-Sm*alpha(er,0)) * p(el,m) * n_d;
              num_flux_left(f,d) /= -(Sp - Sm);

            }

//          un_r = - vit_n(f,k), un_l = - vit_n(f, k + nb_phase), c_l = c(er,k), c_r = c(el,k);
//
//          Sm = std::min(un_l - c_l, un_r - c_r);
//          Sp = std::max(un_l + c_l, un_r + c_r);
//          Sm = std::min(0.0, Sm);
//          Sp = std::max(0.0, Sp);
          k = m;
          un_r = - vit_n(f,k), un_l = - vit_n(f, k + nb_phase), c_l = c(er,k), c_r = c(el,k);
          Sm1 = std::min(un_l - c_l, un_r - c_r);
          Sp1 = std::max(un_l + c_l, un_r + c_r);
          k=n;
          un_r = - vit_n(f,k), un_l = - vit_n(f, k + nb_phase), c_l = c(er,k), c_r = c(el,k);
          Sm2 = std::min(un_l - c_l, un_r - c_r);
          Sp2 = std::max(un_l + c_l, un_r + c_r);
          Sm = std::min(0.0, std::min(Sm1,Sm2));
          Sp = std::max(0.0, std::max(Sp1,Sp2));
          for (int d = 0; d< dimension; d++)
            {
              double n_d = - domaine.face_normales(f,d)/domaine.face_surfaces(f);
              num_flux_right(f,d) = (Sp*alpha(er,0)-Sm*alpha(el,0)) * p(er,m) * n_d;
              num_flux_right(f,d) /= -(Sp - Sm);
            }
        }

      else
        {
          //tableaux de correspondance lies aux CLs : fcl(f, .) = { type de CL, num de la CL, indice de la face dans la CL }
          //types de CL : 0 -> pas de CL
          //              1 -> Neumann
          //              2 -> Navier ou symetrie
          //              3 -> Dirichlet ou Neumann_homogene
          //              4 -> Dirichlet_homogene
          //              5 -> Periodique

          const int& e = f_e(f,0) >= 0 ?  f_e(f,0) :  f_e(f,1); //pas besoin
          assert(f_e(f,0) >= 0 && vit_n (f,0)!=123.123); //pas besoin


          const Conds_lim& cls_qdm = pb.equation_qdm().domaine_Cl_dis().les_conditions_limites();
          const Conds_lim& cls_alpha = pb.equation_fraction().domaine_Cl_dis().les_conditions_limites();
          //const Conds_lim& cls_p = pb.equation_energie().domaine_Cl_dis().les_conditions_limites();

          DoubleTab normal(dimension);
          for (int d = 0; d< dimension; d++) normal(d) = domaine.face_normales(f,d) / domaine.face_surfaces(f);

          if (sub_type(Sortie_supersonique,cls_qdm[fcl(f, 1)].valeur()))
            {

              const double& alpha_bord = alpha(e,0);
              for (int d = 0; d< dimension; d++)
                {
                  num_flux_left(f,d) = -alpha_bord*p(e,m)*normal(d);
                  num_flux_right(f,d) =123.123;
                }
            }
          else if (sub_type(Dirichlet,cls_qdm[fcl(f, 1)].valeur()))
            {
              //const double p_bord = ref_cast(Dirichlet, cls_p[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), m);
              const double alpha_bord = ref_cast(Dirichlet, cls_alpha[fcl(f, 1)].valeur()).val_imp(fcl(f, 2), 0);


              for (int d = 0; d< dimension; d++)
                {
                  num_flux_left(f,d) = -alpha_bord*p(e,m)*normal(d);
                  num_flux_right(f,d) =123.123;
                }
            }


          else if ( sub_type(Symetrie,cls_qdm[fcl(f, 1)].valeur()) && !sub_type(Sortie_supersonique,cls_qdm[fcl(f, 1)].valeur()) )
            {
              //Slip wall : u_n=-u_n


              const double alpha_bord = alpha(e,0);

              for (int d = 0; d< dimension; d++)
                {
                  num_flux_left(f,d) = -alpha_bord*p(e,m)*normal(d);
                  num_flux_right(f,d) =123.123;
                }
            }

          else
            {
              Cerr << " La CL de type " << fcl(f, 0) <<" pour l'equation "<< eq.que_suis_je() << " n est pas diponible .....\n";
              Process::exit();
            }

//      else
//        {
//          const int& e = f_e(f,0) >= 0 ?  f_e(f,0) :  f_e(f,1); //pas besoin
//
//          //tableaux de correspondance lies aux CLs : fcl(f, .) = { type de CL, num de la CL, indice de la face dans la CL }
//          //types de CL : 0 -> pas de CL
//          //              1 -> Neumann
//          //              2 -> Navier ou symetrie
//          //              3 -> Dirichlet ou Neumann_homogene
//          //              4 -> Dirichlet_homogene
//          //              5 -> Periodique
//
//
//          if (fcl(f, 0) == 2)
//            {
//              const double inco_bord =  alpha(e,0) ;
//              const double p_bord =  p(e, m);//
//              for (int d = 0; d< dimension; d++)
//                {
//                  double n_d = domaine.face_normales(f,d) / domaine.face_surfaces(f);
//                  num_flux_left(f,d) = -inco_bord*p_bord*n_d;
//                  num_flux_right(f,d) =123.123;
//                }
//            }
//
//          else if (fcl(f, 0) == 4) //Paroi_fixe dans .data
//            {
//              const double inco_bord =  alpha(e,0) ;
//              const double p_bord =  p(e, m);//
//              for (int d = 0; d< dimension; d++)
//                {
//                  double n_d = domaine.face_normales(f,d) / domaine.face_surfaces(f);
//                  num_flux_left(f,d) = -inco_bord*p_bord*n_d;
//                  num_flux_right(f,d) =123.123;
//                }
//            }


        }
    }
}
