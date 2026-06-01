/****************************************************************************
* Copyright (c) 2024, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#include <Terme_Source_Canal_RANS_LES_VEF_Face.h>
#include <Champ_Uniforme.h>
#include <Domaine_VEF.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
#include <EFichier.h>
#include <Champ_P1NC.h>
#include <SFichier.h>

Implemente_instanciable_sans_destructeur(Terme_Source_Canal_RANS_LES_VEF_Face,"Source_Canal_RANS_LES_VEF_P1NC",Source_base);

Terme_Source_Canal_RANS_LES_VEF_Face::~Terme_Source_Canal_RANS_LES_VEF_Face()
{
  //The destructor is called at initialization, so
  //the field saving is done outside initialization

  if(utemp.size()!=0)
    {
      SFichier fic ("utemp.dat");
      fic << utemp;
    }
}

Sortie& Terme_Source_Canal_RANS_LES_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Terme_Source_Canal_RANS_LES_VEF_Face::readOn(Entree& is )
{
  int compteur = 0;
  Motcle mot_lu;
  Motcle acc_ouverte("{");
  Motcle acc_fermee("}");

  // 0 => spatial average
  // 1 => sliding temporal average (alpha average)
  // 2 => temporal average
  // 3 => temporal average starting at t=0

  Motcles les_mots(11);
  {
    les_mots[0]="alpha_tau";   //source term relaxation coefficient
    les_mots[1]="Ly";          //height of the plane channel (useful for spatial averaging)
    les_mots[2]="t_moy_start"; //averaging start time
    les_mots[3]="f_start";     //time from which the source term starts to be activated (weighted)
    les_mots[4]="f_tot";       //time from which the source term is fully applied
    les_mots[5]="type_moyenne";// 0: spatial average 1: sliding temporal average 2: temporal average
    les_mots[6]="dir";         // forcing direction of the flow
    les_mots[7]="u_tau";       // theoretical friction velocity
    les_mots[8]="nu";          // molecular viscosity
    les_mots[9]="rayon";       // cylinder radius
    les_mots[10]="u_target";   // key for target field: 1=storage 2=reading
  }
  is >> mot_lu;
  if(mot_lu != acc_ouverte)
    {
      Cerr << "Expected { instead of " << mot_lu
           << " while reading the parameters of the wall law " << finl;
    }
  is >> mot_lu;
  while(mot_lu != acc_fermee)
    {
      int rang=les_mots.search(mot_lu);
      switch(rang)
        {
        case 0 :
          is >> alpha_tau;
          Cerr << "alpha_tau = " << alpha_tau << finl;
          compteur++;
          break;
        case 1  :
          is >> Ly;
          Cerr << "Ly = "<< Ly << finl;
          compteur++;
          break;
        case 2  :
          is >> t_moy_start;
          Cerr << " t_moy_start = "<< t_moy_start << finl;
          compteur++;
          break;
        case 3  :
          is >> f_start;
          Cerr << " f_start = "<< f_start << finl;
          compteur++;
          break;
        case 4  :
          is >> f_tot;
          Cerr << " f_tot = "<< f_tot << finl;
          compteur++;
          break;
        case 5  :
          is >> moyenne;
          Cerr << "type_moyenne : " << moyenne << finl;
          if(moyenne==0)
            {
              Cerr << " Relevance of spatial averaging in VEF? -> STOP (not implemented)" << finl;
              exit();
            }
          if(moyenne==1)
            {
              Cerr << " YOUNES, YOUNES !!! !"  << finl;
              Cerr << " YOUNES EST PARTI !!!"  << finl;
              exit();
            }
          compteur++;
          break;
        case 6  :
          is >> dir_nom;
          if     ( (dir_nom=="-X") || (dir_nom=="+X") ) dir=0;
          else if( (dir_nom=="-Y") || (dir_nom=="+Y") ) dir=1;
          else if( (dir_nom=="-Z") || (dir_nom=="+Z") ) dir=2;
          else
            {
              Cerr << " ERROR when reading direction for forcing term " << finl;
              Cerr << " 'dir' parameter has to be selected among following list : +X, -X, +Y, -Y, +Z or -Z " << finl;
              exit();
            }
          Cerr << " forcing along axis: " << dir_nom << " -> component: " << dir <<finl;
          compteur++;
          break;
        case 7  :
          is >> u_tau;
          Cerr << "friction velocity: " << u_tau << finl;
          compteur++;
          break;
        case 8  :
          is >> nu;
          Cerr << "molecular viscosity: " << nu << finl;
          compteur++;
          break;
        case 9  :
          is >> rayon;
          Cerr << "rayon : " << rayon << finl;
          compteur++;
          break;
        case 10  :
          is >> u_target;
          Cerr << " u_target: " << u_target << finl;
          compteur++;
          break;
        default :
          {
            Cerr << mot_lu << " is not a keyword understood by Source_Canal_RANS_LES_VEF_P1NC" << finl;
            Cerr << "Understood keywords are: " << les_mots << finl;
            exit();
          }
        }
      is >> mot_lu;
    }
  Cerr << compteur << " keywords were read in the readOn of Source_Canal_RANS_LES_VEF_P1NC" << finl;

  init();

  return is;

}

void Terme_Source_Canal_RANS_LES_VEF_Face::associer_domaines(const Domaine_dis_base& domaine_dis,
                                                             const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_VEF = ref_cast(Domaine_VEF, domaine_dis);
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis);
}

void Terme_Source_Canal_RANS_LES_VEF_Face::associer_pb(const Probleme_base& pb)
{
}

void Terme_Source_Canal_RANS_LES_VEF_Face::init()
{
  const Domaine_dis_base& zdisbase=mon_equation->inconnue().domaine_dis_base();
  const Domaine_VEF& domaine_VEF=ref_cast(Domaine_VEF, zdisbase);
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  int nb_faces = domaine_VEF.nb_faces();
  const DoubleTab& xv = domaine_VEF.xv();
  const double tps = mon_equation->schema_temps().temps_courant();

  double x1,x2;
  int comp1=0;
  int comp2=0;


  U_RANS.resize(nb_faces,3);
  U_RANS=0.;

  utemp.resize(nb_faces,3);
  utemp = 0.;

  utemp_sum.resize(nb_faces,3);
  utemp_sum = 0.;

  if (tps > t_moy_start )
    {
      EFichier fic1 ("utemp.dat");

      fic1 >> utemp;

      //      int trash;
      //fic1 >> trash;

      //for(int num_face = 0 ; num_face<nb_faces ; num_face++)
      //        {
      //          fic1 >> utemp(num_face,0) >> utemp(num_face,1) >> utemp(num_face,2);
      //        }
    }


  if(u_target==1) // Write target field
    {
      SFichier fic ("utarget.dat");

      for(int num_face=0 ; num_face<nb_faces ; num_face++)
        {
          fic << vitesse(num_face,0) << vitesse(num_face,1) << vitesse(num_face,2);
        }
    }
  if(u_target==2) // Read target field
    {
      EFichier fic ("vitesse_RANS.dat");
      SFichier ficv ("verif_target.dat");

      for(int num_face=0 ; num_face<nb_faces ; num_face++)
        {
          fic >> U_RANS(num_face,0) >> U_RANS(num_face,1) >> U_RANS(num_face,2);
          ficv << U_RANS(num_face,0) << " " << U_RANS(num_face,1) << " " << U_RANS(num_face,2) << finl;
        }
    }
  else // target field: analytical function
    {
      if(dir==0)
        {
          comp1=1 ;
          comp2=2 ;
        }
      if(dir==1)
        {
          comp1=0 ;
          comp2=2 ;
        }
      if(dir==2)
        {
          comp1=0 ;
          comp2=1 ;
        }

      for(int num_face=0 ; num_face<nb_faces ; num_face++)
        {
          x1=xv(num_face,comp1);
          x2=xv(num_face,comp2);

          U_RANS(num_face,dir)= u_tau*(2.44*log((rayon-sqrt(x1*x1+x2*x2))*u_tau/nu)+5.32) ;

          //      FAT : HOT : U_RANS(num_face)= 0.1145*(2.44*log((0.077-sqrt(y*y+z*z))*0.1145/1.52e-7)+5.1) ;
          //      FAT : COLD: U_RANS(num_face)= 0.029*(2.44*log((0.077-sqrt(y*y+x*x))*0.029/6.45e-7)+5.1) ;
        }
      if( dir_nom=="-X" || dir_nom=="-Y" || dir_nom=="-Z") U_RANS*=-1.;
    }

}//end init


void Terme_Source_Canal_RANS_LES_VEF_Face::mettre_a_jour(double temps)
{
  //Cerr << "Inside mettre_a_jour" << finl;

  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  //  const DoubleTab& xv = domaine_VEF.xv();

  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  const double dt = mon_equation->schema_temps().pas_de_temps();
  const double dt_min = mon_equation->schema_temps().pas_temps_min();
  const double tps = mon_equation->schema_temps().temps_courant();
  int nb_faces = domaine_VEF.nb_faces();

  static int cpt=0;

  static double duree=0;

  if(moyenne==1)
    {

      Cerr << "coded in VDF, requires VEF implementation? " << finl;
      exit();

    }// end moyenne == 1
  else if(moyenne==2)
    {
      //******************************************************
      //*************** TEMPORAL AVERAGE *********************
      //******************************************************


      if(tps>=t_moy_start)
        {
          double DT_moy=(tps-dt)-t_moy_start;

          utemp*=DT_moy;

          for (int num_face=0; num_face<nb_faces; num_face++)
            for(int comp=0; comp<3; comp++)
              {
                utemp(num_face,comp) += dt*vitesse(num_face,comp);
              }

          utemp/=(DT_moy+dt);

          //           int face=123;
          //           SFichier fic_temp("utemp.dat", ios::app);
          //           SFichier fic_inst("uinst.dat", ios::app);
          //           fic_temp << tps << "  " <<  utemp(face,0) << "  "
          //                    << xv(face,0) << "  " << xv(face,1) << "  " <<  xv(face,2) << finl;
          //           fic_inst << tps << "  " << vitesse(face,0) << finl;
        }


      //***********************************************************
      //*************** END TEMPORAL AVERAGE  *****************
      //*******************************************************
    }
  else if(moyenne==3)
    {
      //********************************************************************************
      //*************** TIME AVERAGE (starting at t=0 or close...) *******************
      //********************************************************************************

      if(tps>0)
        {
          if(cpt==0)
            {
              for (int num_face=0; num_face<nb_faces; num_face++)
                for(int comp=0; comp<3; comp++)
                  {
                    duree=400*dt;
                    utemp_sum(num_face,comp) = U_RANS(num_face,comp)*duree;
                    duree-=dt;
                    utemp(num_face,comp) = utemp_sum(num_face,comp)/(tps-dt_min+duree);

                    utemp_sum(num_face,comp) += dt*vitesse(num_face,comp);
                    utemp(num_face,comp) = vitesse(num_face,comp);
                  }
              cpt++;

            }
          else
            {
              for (int num_face=0; num_face<nb_faces; num_face++)
                for(int comp=0; comp<3; comp++)
                  {
                    utemp_sum(num_face,comp) += dt*vitesse(num_face,comp);
                    utemp(num_face,comp) = utemp_sum(num_face,comp)/(tps-dt_min+dt+duree);
                  }
            }
        }

    }
  else
    {
      Cerr << "moyenne = " << moyenne << finl;
      Cerr << "Problem selecting the averaging type" << finl;
    }


}//end mettre_a_jour

DoubleTab& Terme_Source_Canal_RANS_LES_VEF_Face::ajouter(DoubleTab& resu) const
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  int nb_faces = domaine_VEF.nb_faces();
  const DoubleVect& porosite_surf = equation().milieu().porosite_face();
  const DoubleVect& volumes_entrelaces = domaine_VEF.volumes_entrelaces();
  const double tps = mon_equation->schema_temps().temps_courant();
  const double dt = mon_equation->schema_temps().pas_de_temps();
  const double dt_min = mon_equation->schema_temps().pas_temps_min();

  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();

  //   DoubleTab resu_old(U_RANS);
  //   DoubleTab utemp_bar(utemp);

  //   if (sub_type(Champ_P1NC,mon_equation->inconnue()))
  //   {
  //    const Champ_P1NC& Ch = ref_cast(Champ_P1NC,mon_equation->inconnue());
  //    Ch.filtrer_L2(U_RANS_bar);
  //    Ch.filtrer_L2(utemp_bar);
  //   }

  double vol = 0.;
  SFichier fic_f("forcing_term.dat", ios::app);
  SFichier fic_finst("finst.dat", ios::app);
  double mbf0 = 0.; // mean body force
  double mbf1 = 0.; // mean body force
  double mbf2 = 0.; // mean body force
  double mbf = 0.; //body force
  //static int cpt = 0;
  double coef=1.;
  double scaling=0.;

  if((tps>f_start)||((moyenne==3)&&(tps>dt_min)&&(tps>f_start)))
    {
      //cpt++;
      if (tps < f_tot) coef*=(tps-f_start)/(f_tot-f_start);

      for(int num_face = 0 ; num_face<nb_faces ; num_face++)
        {
          vol = volumes_entrelaces(num_face)*porosite_surf(num_face);

          for(int comp=0; comp<3; comp++)
            {
              scaling=U_RANS(num_face,comp)/((U_RANS(num_face,0)*U_RANS(num_face,0)
                                              +U_RANS(num_face,1)*U_RANS(num_face,1)
                                              +U_RANS(num_face,2)*U_RANS(num_face,2))+1e-10);
              //                    if(cpt > 150)
              //                      resu(num_face,comp)=resu_old(num_face,comp);
              //                    else
              //                      {
              resu(num_face,comp) = scaling*coef*((U_RANS(num_face,comp)-utemp(num_face,comp))/(alpha_tau*dt))*vol;
              //f(num_face,comp) = /*coef**/((U_RANS_bar(num_face,comp)-utemp_bar(num_face,comp))/(alpha_tau*dt))*vol;
              //                        resu_old(num_face,comp)=resu(num_face,comp);
              //                      }
            }

          mbf0 += resu(num_face,0);
          mbf1 += resu(num_face,1);
          mbf2 += resu(num_face,2);
        }

      mbf0/=nb_faces;
      mbf1/=nb_faces;
      mbf2/=nb_faces;

      int num_face=123;

      scaling=U_RANS(num_face,0)/((U_RANS(num_face,0)*U_RANS(num_face,0)
                                   +U_RANS(num_face,1)*U_RANS(num_face,1)
                                   +U_RANS(num_face,2)*U_RANS(num_face,2))+1e-10);

      mbf = scaling*((U_RANS(num_face,0)-utemp(num_face,0))/(alpha_tau*dt))*vol;

      fic_finst << tps << " " << mbf << " " <<  U_RANS(num_face,0)
                << " " << utemp(num_face,0) << "  " << vitesse(num_face,0)
                << " " << scaling << finl;

      fic_f << tps << " " << mbf0 << " " << mbf1 <<" " << mbf2 << finl;

      if (sub_type(Champ_P1NC,mon_equation->inconnue()))
        {
          const Champ_P1NC& Ch = ref_cast(Champ_P1NC,mon_equation->inconnue());
          Ch.filtrer_resu(resu);
        }


    }//end if f_start

  return resu;
}

DoubleTab& Terme_Source_Canal_RANS_LES_VEF_Face::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

