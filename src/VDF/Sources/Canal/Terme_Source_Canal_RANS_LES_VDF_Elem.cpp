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

#include <Terme_Source_Canal_RANS_LES_VDF_Elem.h>
#include <Champ_Uniforme.h>
#include <Domaine_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Pb_Hydraulique.h>
#include <Pb_Thermohydraulique.h>
#include <EFichier.h>
#include <Interprete.h>
#include <SFichier.h>

Implemente_instanciable_sans_destructeur(Terme_Source_Canal_RANS_LES_VDF_Elem,"Canal_RANS_LES_VDF_P0_VDF",Source_base);



Terme_Source_Canal_RANS_LES_VDF_Elem::~Terme_Source_Canal_RANS_LES_VDF_Elem()
{
  //The destructor is called at initialization so
  //the field saving is done outside initialization

  if(umoy.size()!=0)
    {
      SFichier vit_sauv ("Ttemp_sum.dat");
      vit_sauv << utemp_sum;
      SFichier vit_sauv2 ("Tmoy.dat");
      vit_sauv2 << umoy;
    }
}

//// printOn
//

Sortie& Terme_Source_Canal_RANS_LES_VDF_Elem::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Terme_Source_Canal_RANS_LES_VDF_Elem::readOn(Entree& is )
{
  int compteur = 0;
  Motcle mot_lu;
  Motcle acc_ouverte("{");
  Motcle acc_fermee("}");
  nom_pb_rans="non_couple";

  // 0 => moyenne spatiale
  // 1 => moyenne temporelle glissante (moyenne en alpha)
  // 2 => moyenne temporelle

  Motcles les_mots(6);
  {
    les_mots[0]="alpha_tau"; //relaxation coefficient of the source term
    les_mots[1]="Ly"; //height of the plane channel (useful for spatial averaging)
    les_mots[2]="f_start"; //time from which the source term is activated
    les_mots[3]="t_av"; //averaging time (sliding temporal average)
    les_mots[4]="type_moyenne"; //
    les_mots[5]="nom_pb_rans"; //
  }
  is >> mot_lu;
  if(mot_lu != acc_ouverte)
    {
      Cerr << "Expected { instead of " << mot_lu
           << " while reading the wall law parameters" << finl;
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
          is >> f_start;
          Cerr << "f_start = "<< f_start << finl;
          compteur++;
          break;
        case 3  :
          is >> t_av;
          Cerr << "t_av = " << t_av << finl;
          Cerr << "averaging_type: " << moyenne << finl;
          if((moyenne==1)&&(t_av<=0))
            {
              Cerr << "The time period for the sliding temporal average is not specified!"
                   << finl;
              exit();
            }
          compteur++;
          break;
        case 4  :
          is >> moyenne;
          Cerr << "averaging_type: " << moyenne << finl;
          compteur++;
          break;
        case 5  :
          is >> nom_pb_rans;

          compteur++;
          break;
        default :
          {
            Cerr << mot_lu << " is not a recognized keyword" << finl;
            Cerr << "The recognized keywords are: " << les_mots << finl;
            exit();
          }
        }
      is >> mot_lu;
    }
  Cerr << "nom_pb_rans = " << nom_pb_rans << finl;
  Cerr << compteur << " keywords were read in the readOn of the thermal forcing term" << finl;

  init();

  return is;

}

void Terme_Source_Canal_RANS_LES_VDF_Elem::associer_domaines(const Domaine_dis_base& domaine_dis,
                                                             const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_VDF = ref_cast(Domaine_VDF, domaine_dis);
  le_dom_Cl_VDF = ref_cast(Domaine_Cl_VDF, domaine_Cl_dis);
}

void Terme_Source_Canal_RANS_LES_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  ;
}

void Terme_Source_Canal_RANS_LES_VDF_Elem::init()
{
  const Domaine_dis_base& zdisbase=mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF& domaine_VDF=ref_cast(Domaine_VDF, zdisbase);
  const double tps = mon_equation->schema_temps().temps_courant();

  int nb_elems = domaine_VDF.nb_elem();


  tau.resize(nb_elems);
  U_RANS.resize(nb_elems);
  if(nom_pb_rans == "non_couple")
    {
      SFichier fic_verif("Tverif.RANS");
      EFichier fic_vit("temperature_RANS.dat");

      for(int num_elem=0 ; num_elem<nb_elems ; num_elem++)
        {
          int elem;
          fic_vit >> elem ;
          fic_vit >> U_RANS(elem) ;
          fic_verif << elem << " " << U_RANS(elem) << finl;
        }
    }

  for(int num_elem=0 ; num_elem<nb_elems ; num_elem++)
    {
      tau(num_elem) = alpha_tau;
    }

  utemp_gliss.resize(nb_elems);
  utemp.resize(nb_elems);
  utemp_sum.resize(nb_elems);

  utemp_gliss = 0.;
  utemp = 0.;
  utemp_sum = 0.;

  umoy.resize(nb_elems);

  if (tps > f_start)
    {
      EFichier vit_umoy ("Ttemp_sum.dat");
      EFichier vit_umoy2 ("Tmoy.dat");
      SFichier vit_reprise ("TLES.reprise");
      SFichier vit_reprise2 ("TLES2.reprise");

      int trash;
      vit_umoy >> trash;
      Cerr << "trash = " << trash << finl;
      vit_umoy2 >> trash;
      Cerr << "trash = " << trash << finl;

      for(int num_elem = 0 ; num_elem<nb_elems ; num_elem++)
        {
          vit_umoy >> utemp_sum(num_elem);
          vit_umoy2 >> umoy(num_elem);

          vit_reprise << num_elem << " " << utemp_sum(num_elem) << finl;
          vit_reprise2 << num_elem << " " << umoy(num_elem) << finl;
        }
    }
}//fin init


///////////// Modif elem a partir d'ici///////////////

void Terme_Source_Canal_RANS_LES_VDF_Elem::mettre_a_jour(double temps)
{
  //Cerr << "Je suis dans le mettre_a_jour" << finl;

  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();

  //velocity=temperature

  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  const double dt = mon_equation->schema_temps().pas_de_temps();
  const double tps = mon_equation->schema_temps().temps_courant();
  int nb_elems = domaine_VDF.nb_elem();
  const double dt_min = mon_equation->schema_temps().pas_temps_min();

  int cptbis=0;

  //****************************************************
  //******* Update of the target velocity (RANS) ***********
  //**************************************************
  if(nom_pb_rans != "non_couple")
    {
      OBS_PTR(Probleme_base) pb_rans;
      Objet_U& obj=Interprete::objet(nom_pb_rans);

      if( sub_type(Probleme_base, obj) )
        {
          pb_rans = ref_cast(Probleme_base, obj);
        }

      U_RANS = pb_rans->equation(1).inconnue().valeurs();

    }
  //***************************************************

  if(moyenne==2)
    {
      //******************************************************
      //*************** TEMPORAL AVERAGE *******************
      //******************************************************

      //Compute a first meaningful temporal average

      if((tps>(f_start-t_av))&&(tps<f_start))
        {
          if (cptbis==0)
            {
              for (int num_elem=0; num_elem<nb_elems; num_elem++)
                {
                  utemp_sum(num_elem) = vitesse(num_elem)*(tps-(f_start-t_av));
                  umoy(num_elem) = vitesse(num_elem);
                }
              cptbis = 3;
            }
          else
            {
              for (int num_elem=0; num_elem<nb_elems; num_elem++)
                {
                  utemp_sum(num_elem) += vitesse(num_elem)*dt;
                  umoy(num_elem) = utemp_gliss(num_elem)/(tps-(f_start-t_av));
                }
            }

        }


      if(tps>=f_start)
        {
          for (int num_elem=0; num_elem<nb_elems; num_elem++)
            {
              utemp_sum(num_elem) += dt*vitesse(num_elem);
              umoy(num_elem)=utemp_sum(num_elem)/(tps-(f_start-t_av));
            }
          cptbis=4;
        }


      //***********************************************************
      //*************** END TEMPORAL AVERAGE  ***************
      //*******************************************************
    }

  else if(moyenne==3)
    {
      if (tps>0)
        {
          if (cptbis==0)
            {
              for (int num_elem=0; num_elem<nb_elems; num_elem++)
                {
                  utemp_sum(num_elem) = vitesse(num_elem)*dt;
                  umoy(num_elem) = vitesse(num_elem);
                }
              cptbis ++;
            }
          else
            {
              for (int num_elem=0; num_elem<nb_elems; num_elem++)
                {
                  utemp_sum(num_elem) += vitesse(num_elem)*dt;
                  umoy(num_elem) = utemp_sum(num_elem)/(tps-dt_min+dt);
                }
            }
        }

    }
  else
    {
      Cerr << "moyenne = " << moyenne << finl;
      Cerr << "Problem with the choice of averaging type" << finl;
    }
  compteur_reprise++;

}//fin mettre_a_jour


void Terme_Source_Canal_RANS_LES_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  int nb_elems = domaine_VDF.nb_elem();
  const DoubleVect& volume = domaine_VDF.volumes();
  const double tps = mon_equation->schema_temps().temps_courant();
  const double dt = mon_equation->schema_temps().pas_de_temps();
  const double dt_min = mon_equation->schema_temps().pas_temps_min();

  //velocity=temperature
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();

  double vol=0.;
  SFichier fic_f("f_temp.dat", ios::app);
  SFichier fic_finst("f_temp_inst.dat", ios::app);

  double mbf = 0.; //maximum body force
  double mbf2 = 0.; //maximum body force

  static int cpt2=0;

  // Compute the norm of velocities at element centers

  if(((tps>f_start)&&(compteur_reprise > 1))||((moyenne==3)&&(tps>dt_min)))
    {
      for(int num_elem = 0 ; num_elem<nb_elems ; num_elem++)
        {
          vol = volume(num_elem);

          secmem(num_elem) += ((U_RANS(num_elem)-umoy(num_elem))/(tau(num_elem)*dt))*vol;

          mbf2 +=((U_RANS(num_elem)-umoy(num_elem))/(tau(num_elem)*dt))*vol;

          cpt2++;
        }

      mbf2 = mbf2/cpt2;

      int num_elem=12;

      mbf = (U_RANS(num_elem)-umoy(num_elem))
            /(tau(num_elem)*dt)*vol;

      fic_finst << tps << " " << mbf << " " <<  U_RANS(num_elem)
                << " " << umoy(num_elem) << "  " << vitesse(num_elem) << finl;

      fic_f << tps << " " << mbf2 << " "
            << U_RANS(num_elem) << " " << umoy(num_elem) << " "
            << utemp_sum(num_elem) << " " << dt <<finl;
    }//fin if f_start

}

DoubleTab& Terme_Source_Canal_RANS_LES_VDF_Elem::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}
