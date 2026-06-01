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

#include <Traitement_particulier_NS_Profils_thermo_VDF.h>
#include <Domaine_VDF.h>
#include <Probleme_base.h>
#include <LecFicDistribueBin.h>
#include <EcrFicCollecteBin.h>
#include <Navier_Stokes_std.h>
#include <Convection_Diffusion_Temperature.h>
#include <Convection_Diffusion_Chaleur_QC.h>
#include <communications.h>

Implemente_instanciable_sans_constructeur(Traitement_particulier_NS_Profils_thermo_VDF,"Traitement_particulier_NS_Profils_thermo_VDF",Traitement_particulier_NS_Profils_VDF);
// XD profils_thermo traitement_particulier_base profils_thermo NO_BRACE non documente
// XD attr bloc bloc_lecture bloc REQ not_set

Traitement_particulier_NS_Profils_thermo_VDF::Traitement_particulier_NS_Profils_thermo_VDF():oui_repr_stats_thermo(0)
{
}

Sortie& Traitement_particulier_NS_Profils_thermo_VDF::printOn(Sortie& is) const
{
  return is;
}

Entree& Traitement_particulier_NS_Profils_thermo_VDF::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_Profils_thermo_VDF::lire(Entree& is)
{
  return Traitement_particulier_NS_Profils_VDF::lire(is);
}

Entree& Traitement_particulier_NS_Profils_thermo_VDF::lire(const Motcle& motlu, Entree& is)
{
  // Par defaut on ne fait pas de statistiques sur la temperature.

  if (motlu=="stats_thermo")
    {
      //Statistiques sur le champ de temperature
      oui_stats_thermo = 1;
    }
  else
    {
      if(motlu=="reprise_thermo")
        {
          // Reprise des statistiques pour le champs de temperature
          oui_repr_stats_thermo = 1; // On veut reprendre les statistiques sur la temperature
          is  >> fich_repr_stats_thermo ;  // Indication du nom du fichier de reprise des stats pour le champs de temperature
          Cerr << "On reprend les statistiques sur le champs thermique" << finl;
        }
      else
        {
          if(motlu=="tmin_tmax")
            {
              // On demande les max et min de temperature dans un fichier
              tmin_tmax = 1;
              Cerr << "User asked to monitor min and max of temperature" << finl;
              Cerr << "in channel flow computation." << finl;
            }
          else
            {
              return Traitement_particulier_NS_Profils_VDF::lire(motlu,is);
            }
        }
    }
  return is;
}

void Traitement_particulier_NS_Profils_thermo_VDF::associer_eqn(const Equation_base& eq_ns )
{
  Traitement_particulier_NS_Profils_VDF::associer_eqn(eq_ns);

  const Probleme_base& pb = mon_equation->probleme();
  int flag=0;
  for(int i=0; i<pb.nombre_d_equations(); i++)
    {
      if(sub_type(Convection_Diffusion_Temperature,pb.equation(i)))
        {
          mon_equation_NRJ = ref_cast(Convection_Diffusion_Temperature,pb.equation(i));
          flag=1;
        }
      else if(sub_type(Convection_Diffusion_Chaleur_QC,pb.equation(i)))
        {
          mon_equation_NRJ = ref_cast(Convection_Diffusion_Chaleur_QC,pb.equation(i));
          flag=1;
        }

    }
  if (flag==0)
    {
      Cerr << "Error : Equation of NRJ was not found..." << finl;
      Cerr << "User can not ask for statistics on temperature" << finl;
      Cerr << "if we do not solve a heat equation." << finl;
      Cerr << "Try to remove the '_thermo' part at key-word 'Profils'" << finl;
      exit();
    }
}

// #################### Postraitement Statistique du Champ de Temperature #########

void Traitement_particulier_NS_Profils_thermo_VDF::post_traitement_particulier()
{
  Traitement_particulier_NS_Profils_VDF::post_traitement_particulier();

  DoubleTab Tmoy_m(n_probes,Nap);
  Tmoy_m=0.;
  DoubleTab Trms_m(n_probes,Nap);
  Trms_m=0.;
  DoubleTab upTp_m(n_probes,Nap);
  upTp_m=0.;
  DoubleTab vpTp_m(n_probes,Nap);
  vpTp_m=0.;
  DoubleTab wpTp_m(n_probes,Nap);
  wpTp_m=0.;

  DoubleTab Tmoy_p(n_probes,Nap);
  Tmoy_p=0.;
  DoubleTab Trms_p(n_probes,Nap);
  Trms_p=0.;
  DoubleTab upTp_p(n_probes,Nap);
  upTp_p=0.;
  DoubleTab vpTp_p(n_probes,Nap);
  vpTp_p=0.;
  DoubleTab wpTp_p(n_probes,Nap);
  wpTp_p=0.;

  double tps = mon_equation->inconnue().temps();


  // Statistics on the thermal field are only computed simultaneously with statistics on the dynamic field.
  if ((oui_u_inst != 0)&&(oui_stats_thermo != 0))
    {
      calculer_moyennes_spatiales_thermo(Tmoy_m,Trms_m,upTp_m,vpTp_m,wpTp_m,corresp_uv_m,compt_uv_m,Nuv,xUVm);
      calculer_moyennes_spatiales_thermo(Tmoy_p,Trms_p,upTp_p,vpTp_p,wpTp_p,corresp_uv_p,compt_uv_p,Nuv,xUVp);

      static double temps_dern_post_inst = -100.;
      if (std::fabs(tps-temps_dern_post_inst)>=dt_post_inst)
        {
          ecriture_fichier_moy_spat_thermo(Tmoy_m,Trms_m,upTp_m,vpTp_m,wpTp_m,Tmoy_p,Trms_p,upTp_p,vpTp_p,wpTp_p,Yuv_m,Nuv,delta_UVm, delta_UVp);

          temps_dern_post_inst = tps;
        }
    }


  // Temporal averages:

  // Statistics on the thermal field are only computed simultaneously with statistics on the dynamic field.
  if ((oui_u_inst != 0)&&(oui_stats_thermo != 0)&&(oui_stat != 0))
    {
      double tpsbis = mon_equation->inconnue().temps();
      if ((tpsbis>=temps_deb)&&(tpsbis<=temps_fin))
        {
          static int init_stat_temps = 0;
          if((init_stat_temps==0)&&( oui_repr != 1))  // if this is not a restart: otherwise values are read from the file
            {
              double dt_v = mon_equation->schema_temps().pas_de_temps();
              temps_deb = tpsbis-dt_v;
              init_stat_temps++;
            }

          calculer_integrales_temporelles(Tmoy_temp,Tmoy_m,Tmoy_p,delta_UVm,delta_UVp);
          calculer_integrales_temporelles(Trms_temp,Trms_m,Trms_p,delta_UVm,delta_UVp);
          calculer_integrales_temporelles(upTp_temp,upTp_m,upTp_p,delta_UVm,delta_UVp);
          calculer_integrales_temporelles(vpTp_temp,vpTp_m,vpTp_p,delta_UVm,delta_UVp);
          calculer_integrales_temporelles(wpTp_temp,wpTp_m,wpTp_p,delta_UVm,delta_UVp);

          static double temps_dern_post_stat = -100.;
          if (std::fabs(tpsbis-temps_dern_post_stat)>=dt_post_stat)
            {
              double dt = tpsbis-temps_deb;
              ecriture_fichier_moy_temp_thermo(Tmoy_temp,Trms_temp,upTp_temp,vpTp_temp,wpTp_temp,Yuv_m,dt,Nuv);
              temps_dern_post_stat = tpsbis;
            }
        }
    }
}





// #################### Compute Averages #################################

void Traitement_particulier_NS_Profils_thermo_VDF::calculer_moyennes_spatiales_thermo(DoubleTab& tmoy, DoubleTab& trms, DoubleTab& uptp, DoubleTab& vptp, DoubleTab& wptp, const IntTab& corresp, const IntTab& compt, const IntVect& NN, const DoubleTab& xUV)
{

  const Domaine_dis_base& zdisbase=mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF& domaine_VDF=ref_cast(Domaine_VDF, zdisbase);
  const IntTab& elem_faces = domaine_VDF.elem_faces();

  // To find the coordinates of the max and min temperature points:
  const DoubleTab& xp = domaine_VDF.xp();

  // The number of elements in the VDF domain.
  int nb_elems = domaine_VDF.domaine().nb_elem();

  // Access the values of the 3 velocity components from mon_equation.
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();

  // Access the temperature values from mon_equation_NRJ.
  const DoubleTab& Temp = mon_equation_NRJ->inconnue().valeurs();



  // These 3 vectors hold the values of the 3 velocity components U, V and W averaged at element centers over the entire space.
  // They correspond to <U>, <V> and <W>.
  DoubleTab u_moy_cent(n_probes,Nap);
  u_moy_cent= 0.;
  DoubleTab v_moy_cent(n_probes,Nap);
  v_moy_cent= 0.;
  DoubleTab w_moy_cent(n_probes,Nap);
  w_moy_cent= 0.;

  // t2moy corresponds to
  DoubleTab t2moy(n_probes,Nap);
  t2moy     = 0.;

  // vitu, vitv and vitw are the values of the 3 velocity components taken at element centers,
  // i.e. averaged.
  double vitu,vitv,vitw;

  // For checking temperature bounds in plane channel flow
  double Tmin=1000000000.,Tmax=0.;
  int elem_min=0, elem_max=0;

  // face_ui_j is used to get the face number carrying the velocity component ui,
  // with 0 and 1 being the two opposing faces used to average at the element center.
  // This approach is necessary among other reasons because we iterate over elements.
  int face_u_0,face_u_1,face_v_0,face_v_1,face_w_0,face_w_1;

  int i,j,deja_fait=0;
  int num_elem;

  // tmoy is the mean of T : <T>(y,t)
  tmoy = 0.;
  // trms is the standard deviation of temperature : sqrt(<Tp^2>)(y,t)=<T^2>-<T>^2
  trms = 0.;
  // uptp : <U><T>-<UT>
  uptp = 0.;
  // vptp : <V><T>-<VT>
  vptp = 0.;
  // wptp : <W><T>-<WT>
  wptp = 0.;

  for(i=0; i<n_probes; i++)
    {
      // We iterate over all elements to compute all averages at element centers.
      for (num_elem=0; num_elem<nb_elems; num_elem++)
        {
          // <T>
          if(xp(num_elem,dir_profil)==xUV(i))
            {
              tmoy(i,corresp(i,num_elem)) += Temp[num_elem];
              face_u_0 = elem_faces(num_elem,0);
              face_u_1 = elem_faces(num_elem,dimension);
              vitu = 0.5*(vitesse[face_u_0]+vitesse[face_u_1]);

              face_v_0 = elem_faces(num_elem,1);
              face_v_1 = elem_faces(num_elem,dimension+1);
              vitv = 0.5*(vitesse[face_v_0]+vitesse[face_v_1]);

              face_w_0 = elem_faces(num_elem,2);
              face_w_1 = elem_faces(num_elem,dimension+2);
              vitw = 0.5*(vitesse[face_w_0]+vitesse[face_w_1]);
              // <U>
              u_moy_cent(i,corresp(i,num_elem)) += vitu;
              // <V>
              v_moy_cent(i,corresp(i,num_elem)) += vitv;
              // <W>
              w_moy_cent(i,corresp(i,num_elem)) += vitw;
              // Trms : here we compute <T^2> but then subtract the <T>^2 part
              trms(i,corresp(i,num_elem)) += Temp[num_elem]*Temp[num_elem];
              // Compute <uT>, <vT>, and <wT>
              uptp(i,corresp(i,num_elem)) += vitu*Temp[num_elem];
              vptp(i,corresp(i,num_elem)) += vitv*Temp[num_elem];
              wptp(i,corresp(i,num_elem)) += vitw*Temp[num_elem];
            }

          // To find the max and min temperatures in the flow.
          if(Temp[num_elem]<Tmin)
            {
              Tmin=Temp[num_elem];
              elem_min=num_elem;
            }
          if(Temp[num_elem]>Tmax)
            {
              Tmax=Temp[num_elem];
              elem_max=num_elem;
            }

        }//END loop over elements

      if((tmin_tmax==1)&&(deja_fait==0))
        {
          SFichier fic1("T_min_max.dat",ios::app);
          double tps = mon_equation->inconnue().temps();
          fic1 << tps << " " << Tmin << " " << xp(elem_min,0) << " " << xp(elem_min,1);
          fic1 << " " << Tmax << " " << xp(elem_max,0) << " " << xp(elem_max,1) <<finl;
          fic1.flush();
          fic1.close();
          deja_fait=1;
        }
    } // All profiles have been iterated over

  // NN is the number of distinct Y values, i.e. the number of points in the final curve after averaging.

  // compt[j] is the number of elements that were used to compute a given mean temperature value,
  // i.e. the number of elements sharing the same Y coordinate.
  // FOR PARALLEL !!
  IntTab compt_p(compt);
  envoyer(compt_p,Process::me(),0,Process::me());

  DoubleTab tmoy_p(tmoy);
  envoyer(tmoy_p,Process::me(),0,Process::me());

  DoubleTab trms_p(trms);
  envoyer(trms_p,Process::me(),0,Process::me());

  DoubleTab uptp_p(uptp);
  envoyer(uptp_p,Process::me(),0,Process::me());

  DoubleTab vptp_p(vptp);
  envoyer(vptp_p,Process::me(),0,Process::me());

  DoubleTab wptp_p(wptp);
  envoyer(wptp_p,Process::me(),0,Process::me());

  DoubleTab u_moy_cent_p(u_moy_cent);
  envoyer(u_moy_cent_p,Process::me(),0,Process::me());

  DoubleTab v_moy_cent_p(v_moy_cent);
  envoyer(v_moy_cent_p,Process::me(),0,Process::me());

  DoubleTab w_moy_cent_p(w_moy_cent);
  envoyer(w_moy_cent_p,Process::me(),0,Process::me());

  if(je_suis_maitre())
    {
      IntTab compt_tot(compt);
      DoubleTab tmoy_tot(tmoy);
      DoubleTab trms_tot(trms);
      DoubleTab uptp_tot(uptp);
      DoubleTab vptp_tot(vptp);
      DoubleTab wptp_tot(wptp);
      DoubleTab u_moy_cent_tot(u_moy_cent);
      DoubleTab v_moy_cent_tot(v_moy_cent);
      DoubleTab w_moy_cent_tot(w_moy_cent);

      compt_tot=0;
      tmoy_tot=0.;
      trms_tot=0.;
      uptp_tot=0.;
      vptp_tot=0.;
      wptp_tot=0.;
      u_moy_cent_tot=0.;
      v_moy_cent_tot=0.;
      w_moy_cent_tot=0.;

      for(int p=0; p<Process::nproc(); p++)
        {
          recevoir(compt_p,p,0,p);
          compt_tot+=compt_p;

          recevoir(tmoy_p,p,0,p);
          tmoy_tot+=tmoy_p;

          recevoir(trms_p,p,0,p);
          trms_tot+=trms_p;

          recevoir(uptp_p,p,0,p);
          uptp_tot+=uptp_p;

          recevoir(vptp_p,p,0,p);
          vptp_tot+=vptp_p;

          recevoir(wptp_p,p,0,p);
          wptp_tot+=wptp_p;

          recevoir(u_moy_cent_p,p,0,p);
          u_moy_cent_tot+=u_moy_cent_p;

          recevoir(v_moy_cent_p,p,0,p);
          v_moy_cent_tot+=v_moy_cent_p;

          recevoir(w_moy_cent_p,p,0,p);
          w_moy_cent_tot+=w_moy_cent_p;


        }

      for(i=0; i<n_probes; i++)
        {
          for (j=0; j<NN(i); j++)
            {
              tmoy(i,j) = tmoy_tot(i,j) / compt_tot(i,j);
              trms(i,j) = trms_tot(i,j) / compt_tot(i,j);
              uptp(i,j) = uptp_tot(i,j) / compt_tot(i,j);
              vptp(i,j) = vptp_tot(i,j) / compt_tot(i,j);
              wptp(i,j) = wptp_tot(i,j) / compt_tot(i,j);
              u_moy_cent(i,j) = u_moy_cent_tot(i,j) / compt_tot(i,j);
              v_moy_cent(i,j) = v_moy_cent_tot(i,j) / compt_tot(i,j);
              w_moy_cent(i,j) = w_moy_cent_tot(i,j) / compt_tot(i,j);

              trms(i,j) -= tmoy(i,j)*tmoy(i,j);
              uptp(i,j) -= u_moy_cent(i,j)*tmoy(i,j);
              vptp(i,j) -= v_moy_cent(i,j)*tmoy(i,j);
              wptp(i,j) -= w_moy_cent(i,j)*tmoy(i,j);
            }
        }// End loop over profiles

    }// END Parallel section

}



// #################### Compute Temporal Integral ###############################

void Traitement_particulier_NS_Profils_thermo_VDF::calculer_integrales_temporelles(DoubleTab& moy_temp, const DoubleTab& moy_spat_m, const DoubleTab& moy_spat_p, const DoubleVect& delta_m, const DoubleVect& delta_p)
{
  int i,j;
  double dt_v = mon_equation->schema_temps().pas_de_temps();
  DoubleTab moy(moy_temp);
  for(i=0; i<n_probes; i++)
    {
      for(j=0; j<Nap; j++)
        {
          moy(i,j)=moy_spat_m(i,j)*delta_m(i)/(delta_m(i)+delta_p(i));
          moy(i,j)+=moy_spat_p(i,j)*delta_p(i)/(delta_m(i)+delta_p(i));
        }
    }
  if(je_suis_maitre())
    ////moy_temp.ajoute(dt_v,moy);
    moy_temp.ajoute_sans_ech_esp_virt(dt_v,moy);

}





// #################### Writing Spatial Averages to File ##################


void Traitement_particulier_NS_Profils_thermo_VDF::ecriture_fichier_moy_spat_thermo(const DoubleTab& Tmoy_m, const DoubleTab& Trms_m, const DoubleTab& upTp_m, const DoubleTab& vpTp_m, const DoubleTab& wpTp_m, const DoubleTab& Tmoy_p, const DoubleTab& Trms_p, const DoubleTab& upTp_p, const DoubleTab& vpTp_p, const DoubleTab& wpTp_p, const DoubleTab& Y, const IntVect& NN, const DoubleVect& delta_m,  const DoubleVect& delta_p)
{

  int i,j;
  double tps = mon_equation->inconnue().temps();
  Nom temps = Nom(tps);

  for(i=0; i<n_probes; i++)
    {
      Nom nom_fic = "./Space_Avg/Avg_temp_";
      switch(dir_profil)
        {
        case 0:
          {
            nom_fic +="X_";
            break;
          }
        case 1:
          {
            nom_fic +="Y_";
            break;
          }
        case 2:
          {
            nom_fic +="Z_";
            break;
          }
        }
      Nom pos  = Nom(positions(i));
      nom_fic += pos;
      nom_fic += "_t_";
      nom_fic += temps;
      nom_fic += ".dat";

      if(je_suis_maitre())
        {
          SFichier fic (nom_fic);
          fic << "# Space Averaged Statistics" << finl ;
          fic << "# Y       <T>       Trms       <upTp>       -<vpTp>       -<wpTp>" << finl ;

          DoubleTab Tmoy(Tmoy_m);
          DoubleTab Trms(Trms_m);
          DoubleTab upTp(upTp_m);
          DoubleTab vpTp(vpTp_m);
          DoubleTab wpTp(wpTp_m);

          for(j=0; j<NN(i); j++)
            {
              Tmoy(i,j)=Tmoy_m(i,j)*delta_m(i)/(delta_m(i)+delta_p(i));
              Tmoy(i,j)+=Tmoy_p(i,j)*delta_p(i)/(delta_m(i)+delta_p(i));
              Trms(i,j)=Trms_m(i,j)*delta_m(i)/(delta_m(i)+delta_p(i));
              Trms(i,j)+=Trms_p(i,j)*delta_p(i)/(delta_m(i)+delta_p(i));
              upTp(i,j)=upTp_m(i,j)*delta_m(i)/(delta_m(i)+delta_p(i));
              upTp(i,j)+=upTp_p(i,j)*delta_p(i)/(delta_m(i)+delta_p(i));
              vpTp(i,j)=vpTp_m(i,j)*delta_m(i)/(delta_m(i)+delta_p(i));
              vpTp(i,j)+=vpTp_p(i,j)*delta_p(i)/(delta_m(i)+delta_p(i));
              wpTp(i,j)=wpTp_m(i,j)*delta_m(i)/(delta_m(i)+delta_p(i));
              wpTp(i,j)+=wpTp_p(i,j)*delta_p(i)/(delta_m(i)+delta_p(i));
            }

          for (j=0; j<NN(i); j++)
            fic << Y(i,j) << " " << Tmoy(i,j) << " " << sqrt(std::max(Trms(i,j),0.0)) << " " << upTp(i,j) << " " << -vpTp(i,j) << " " << -wpTp(i,j) << finl;
          fic.flush();
          fic.close();
        } //End of "if I am master"
    }//END loop over probes
}





// #################### Writing Temporal Averages to File ################

void Traitement_particulier_NS_Profils_thermo_VDF::ecriture_fichier_moy_temp_thermo(const DoubleTab& Tmoy, const DoubleTab& Trms, const DoubleTab& upTp, const DoubleTab& vpTp, const DoubleTab& wpTp, const DoubleTab& Y, const double dt, const IntVect& NN)
{

  int i,j;
  for(i=0; i<n_probes; i++)
    {
      Nom nom_fic = "./Time_Avg/Avg_time_temp_";
      double tps = mon_equation->inconnue().temps();

      switch(dir_profil)
        {
        case 0:
          {
            nom_fic +="X_";
            break;
          }
        case 1:
          {
            nom_fic +="Y_";
            break;
          }
        case 2:
          {
            nom_fic +="Z_";
            break;
          }
        }
      Nom pos  = Nom(positions(i));
      nom_fic += pos;
      nom_fic +="_t_";

      Nom temps = Nom(tps);
      nom_fic+= temps;

      nom_fic+=".dat";

      if(je_suis_maitre())
        {
          SFichier fic (nom_fic);
          fic << "# Time Averaged Statistics" << finl ;
          fic << "# Y       <T>       Trms       <upTp>       -<vpTp>       -<wpTp>" << finl ;


          for (j=0; j<NN(i)  ; j++)
            fic << Y(i,j) << " " << Tmoy(i,j)/dt << " " << sqrt(std::max(Trms(i,j)/dt,0.0)) << " " << upTp(i,j)/dt << " " << -vpTp(i,j)/dt << " " << -wpTp(i,j)/dt << finl;
          fic.flush();
          fic.close();
        }
    }//END loop over probes
}





// #################### Saving temporal statistics ###################

void Traitement_particulier_NS_Profils_thermo_VDF::sauver_stat() const
{

  Traitement_particulier_NS_Profils::sauver_stat();
  double tps = mon_equation->inconnue().temps();

  if (  (oui_stat == 1)&&(tps>=temps_deb)&&(tps<=temps_fin) )
    {
      Cerr << "In Traitement_particulier_NS_Profils_thermo_VDF::sauver_stat" << finl;
      // Save the sum (without dividing by dt)
      int i,j;
      Nom temps = Nom(tps);
      Nom fich_sauv_temp ="temperature_field_avg_time_";
      fich_sauv_temp+=temps;
      fich_sauv_temp+=".sauv";

      // Save u_moy!!
      EcrFicCollecteBin fict (fich_sauv_temp);
      fict << temps << finl;
      for(i=0; i<n_probes; i++)
        {
          fict << Nuv(i) << finl;
          for (j=0; j<Nuv(i); j++)
            {
              fict << Tmoy_temp(i,j) << " " << finl;
              fict << Trms_temp(i,j) << " " << finl;
              fict << upTp_temp(i,j) << " " << finl;
              fict << vpTp_temp(i,j) << " " << finl;
              fict << wpTp_temp(i,j) << " " << finl;
            }
        }
      fict << temps_deb << finl;
      fict.flush();
      fict.close();
      Cerr << "Temperature statistics saved at t=" << tps << finl;
    }
}





void Traitement_particulier_NS_Profils_thermo_VDF::reprendre_stat()
{
  Traitement_particulier_NS_Profils::reprendre_stat();
  double tps = mon_equation->inconnue().temps();
  int i,j;
  double ti;
  // double ti2;
  Nom temps;
  if (oui_stat == 1)
    {
      ifstream fict(fich_repr_stats_thermo);
      if (fict)
        {
          LecFicDistribueBin fict2 (fich_repr_stats_thermo);
          fict2 >> temps;
          Tmoy_temp.resize(n_probes,Nap);
          Trms_temp.resize(n_probes,Nap);
          upTp_temp.resize(n_probes,Nap);
          vpTp_temp.resize(n_probes,Nap);
          wpTp_temp.resize(n_probes,Nap);

          for(i=0; i<n_probes; i++)
            {
              fict2 >> Nuv(i);
              for (j=0; j<Nuv(i); j++)
                {
                  fict2 >> Tmoy_temp(i,j) ;
                  fict2 >> Trms_temp(i,j);
                  fict2 >> upTp_temp(i,j) ;
                  fict2 >> vpTp_temp(i,j) ;
                  fict2 >> wpTp_temp(i,j) ;
                }
            }
          fict2 >> ti ;
          Cerr << "ti=" << ti << finl;
          Nom chc_ti = Nom(ti), chc_temps_deb = Nom(temps_deb);
          Cerr << "chc_ti=" << chc_ti << "   chc_temps_deb=" << chc_temps_deb << finl;
          if (chc_ti!=chc_temps_deb)
            {
              if (ti > temps_deb)
                {
                  Cerr << "Pb de reprise des stats du champ thermique :" << finl;
                  Cerr << "Le temps de debut des stats demande " << temps_deb  << finl;
                  Cerr << "est inferieur a celui de debut des stats sauvees !" << ti << finl;
                  exit();
                }
              else
                {
                  if (temps_deb>=tps)
                    {
                      Cerr << "On recommence le calcul des stats thermiques a partir de t=" << temps_deb << "s." << finl;
                      Tmoy_temp = 0.;
                      Trms_temp = 0.;
                      upTp_temp = 0.;
                      vpTp_temp = 0.;
                      wpTp_temp = 0.;
                      oui_repr_stats_thermo = 0;
                    }
                  else
                    {
                      Cerr << "On a deja depasse le temps auquel vous voulez commencer les stats thermiques !" << finl;
                      exit();
                    }
                }
            }
          else
            {
              Cerr << "On continue le calcul des stats thermiques, debute a t=" <<  ti << "s." << finl;
              temps_deb = ti;
            }
        }
      else
        {
          if (tps<=temps_deb)
            {
              Cerr << "On n a pas encore debute le calcul des stats thermiques." << finl;
              oui_repr_stats_thermo = 0;
            }
          else
            {
              Cerr << "Il faut donner le fichier pour reprendre les stats thermiques !" << finl;
              exit();
            }
        }
    }
  else
    {
      //      Cerr << "Pas de calcul des stats commence." << finl;
      //      Cerr << "Pas de reprise des statistiques thermiques." << finl;
    }
}





void  Traitement_particulier_NS_Profils_thermo_VDF::init_calcul_moyenne()
{
  Traitement_particulier_NS_Profils_VDF::init_calcul_moyenne();
}


void  Traitement_particulier_NS_Profils_thermo_VDF::preparer_calcul_particulier()
{
  if ((oui_u_inst != 0)||(oui_profil_nu_t != 0))
    // Call only the method in NS_Profils_VDF to initialize
    // the averaging computation, since only the correspondence tables
    // and the sizing of Yuv and Nuv quantities are done there.

    Traitement_particulier_NS_Profils_VDF::init_calcul_moyenne();

  if (oui_stat != 0)
    // However, to start computing the stats a small processing step is needed
    // for the thermo since the arrays for the temporal averages
    // specific to thermo, as well as those of the classical hydraulics in NS_Profils, must be sized.

    init_calcul_stats();
}




void Traitement_particulier_NS_Profils_thermo_VDF::init_calcul_stats()
{
  Traitement_particulier_NS_Profils::init_calcul_stats();

  if (oui_repr_stats_thermo!=1)
    {
      // The user does not request a restart, meaning there will be no reading
      // from a file and the arrays must therefore be sized manually to start the statistics.
      // These arrays are declared in Traitement_particulier_NS_Profils_thermo_VDF.h

      Tmoy_temp.resize(n_probes,Nap);
      Tmoy_temp=0;
      Trms_temp.resize(n_probes,Nap);
      Trms_temp=0;
      upTp_temp.resize(n_probes,Nap);
      upTp_temp=0;
      vpTp_temp.resize(n_probes,Nap);
      vpTp_temp=0;
      wpTp_temp.resize(n_probes,Nap);
      wpTp_temp=0;
    }
}
