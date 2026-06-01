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

#include <Traitement_particulier_NS_canal_VDF.h>
#include <Domaine_VDF.h>

#include <Fluide_base.h>
#include <Navier_Stokes_std.h>
#include <Modele_turbulence_hyd_base.h>

Implemente_instanciable(Traitement_particulier_NS_canal_VDF,"Traitement_particulier_NS_canal_VDF",Traitement_particulier_NS_canal);


/*! @brief
 *
 * @param is output stream
 * @return modified output stream
 */
Sortie& Traitement_particulier_NS_canal_VDF::printOn(Sortie& is) const
{
  return is;
}


/*! @brief
 *
 * @param is input stream
 * @return modified input stream
 */
Entree& Traitement_particulier_NS_canal_VDF::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_canal_VDF::lire(Entree& is)
{
  return Traitement_particulier_NS_canal::lire(is);
}

void Traitement_particulier_NS_canal_VDF::remplir_Y(DoubleVect& tab_Y,  DoubleVect& tab_compt, int& nNy) const
{
  // Initialize the various member parameters of the class
  // needed for computing the various spatial averages
  // Initialization of: Y, compt

  const Domaine_dis_base& zdisbase = mon_equation->inconnue().domaine_dis_base();
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, zdisbase);
  const DoubleTab& xp = domaine_VF.xp();
  int nb_elems = domaine_VF.domaine().nb_elem();
  int num_elem,j,indic,trouve;
  double y;

  j=0;
  indic = 0;

  tab_Y.resize(1);
  tab_compt.resize(1);
  tab_Y = -100.;
  tab_compt = 0;

  //Fill the Y array
  ////////////////////////////////////////////////////////

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      y = xp(num_elem,1);
      trouve = 0;

      for (j=0; j<indic+1; j++)
        {
          if(est_egal(y,tab_Y[j]))
            {
              tab_compt[j] ++;
              j=indic+1;
              trouve = 1;
              break;
            }
        }
      if (trouve==0)
        {
          tab_Y[indic]=y;
          tab_compt[indic] ++;
          indic++;

          tab_Y.resize(indic+1);
          tab_Y(indic)=-100.;
          tab_compt.resize(indic+1);
        }
    }

  nNy = indic;

  tab_Y.resize(nNy);
  tab_compt.resize(nNy);
}

//Addition F.A 15/02/11: reorganizing operations that were repeated
// to make full use of the memory space (loops),
// by creating a large array (~7M per proc but not exchanged).
// The array has the following structure: the row is the element index,
// index of the element above, index of the element below, position in the Y vector.
// i.e. an array of nelem x 3.
// After remplir_Y, this function is called to perform the various calculations,
// using this array as a function argument.

void Traitement_particulier_NS_canal_VDF::remplir_Tab_recap(IntTab& Tab_rec) const
{
  const Domaine_dis_base& zdisbase=mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF& domaine_VDF=ref_cast(Domaine_VDF, zdisbase);
  const DoubleTab& xp = domaine_VDF.xp();
  const IntTab& elem_faces = domaine_VDF.elem_faces();

  int face; //face receiver
  int elem_test,elem_test2; // element for testing fictitious elements
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot(); // total number of elements (real + fictitious)
  int nb_elems = domaine_VDF.domaine().nb_elem();

  IntTab trouve(1);// array of elements already processed
  double y=0;
  int i,num_elem; // counters
  int q=1; //Cursor for the upper arrays
  trouve[0]=0;


  Tab_rec.resize(nb_elems,3); // Size the array.

  for (num_elem=nb_elems; num_elem<nb_elem_tot; num_elem++) // loop over fictitious elements
    {
      face = elem_faces(num_elem,1+dimension);
      elem_test=domaine_VDF.elem_voisin(num_elem,face,0);
      face = elem_faces(num_elem,1);
      elem_test2=domaine_VDF.elem_voisin(num_elem,face,1);

      if ((elem_test>0) && (elem_test<nb_elems)) // if the element above is a real element
        {
          trouve[q-1]=elem_test;
          q =q +1;
          trouve.resize(q);

          Tab_rec(elem_test,0)=num_elem; // assign the same value to both upper and lower slots
          Tab_rec(elem_test,1)=num_elem; //so the function computing values sees a normal element

          y=xp(elem_test,1);
          for (i=0; i<Ny; i++)
            {
              if(est_egal(y,Y[i]))
                break;
            }

          Tab_rec(elem_test,2)=i; // store i to avoid repeating the loop at each time step
        }
      else if ((elem_test2<nb_elems)&&(elem_test2>0)) //otherwise if the element below is a real element
        {
          trouve[q-1]=elem_test2;
          q =q +1;
          trouve.resize(q);

          Tab_rec(elem_test2,0)=num_elem; // assign the same value to both upper and lower slots
          Tab_rec(elem_test2,1)=num_elem; //so the function computing values sees a normal element

          y=xp(elem_test2,1);
          for (i=0; i<Ny; i++)
            {
              if(est_egal(y,Y[i]))
                break;
            }
          Tab_rec(elem_test2,2)=i; // store i to avoid repeating the loop at each time step

        }
      // otherwise do nothing
    }

  Cerr << "Traitement particulier canal: there is a possible improvement to make for boundary faces!! " << finl;
  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      q=0;// reuse counter q (no longer needed) to check whether an equivalent was found
      for(i=0; i<(trouve.size()-1); i++) // trouve is one slot too large, but instead of resizing it we use the criterion size-1
        if(num_elem==trouve[i])
          {
            q = 0;  // fix q so it cannot pass the next test. // correction: set q=0 because it was a false problem
            break;
          }
      // in reality lambda diverges at the interface
      if(q==0) //
        {
          face=elem_faces(num_elem,1); //lower face
          elem_test=domaine_VDF.elem_voisin(num_elem,face,1);


          if (elem_test+1)
            {
              Tab_rec(num_elem,1)=elem_test; // wrong if elem_test=-1, otherwise fills with the element below
            }
          else
            {
              Tab_rec(num_elem,1)=domaine_VDF.elem_voisin(num_elem,elem_faces(num_elem,1+dimension),0); // treat it as a virtual element
            }

          face= elem_faces(num_elem,1+dimension); //upper face
          elem_test=domaine_VDF.elem_voisin(num_elem,face,0);

          if (elem_test+1)
            {
              Tab_rec(num_elem,0)=elem_test; // wrong if elem_test=-1, otherwise fills with the element above
            }
          else
            {
              Tab_rec(num_elem,0)=domaine_VDF.elem_voisin(num_elem,elem_faces(num_elem,1),1); // treat it as a virtual element
            }

          y = xp(num_elem,1);
          for (i=0; i<Ny; i++)
            if(est_egal(y,Y[i])) break;

          Tab_rec(num_elem,2)=i;

        }

    }
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_vitesse_rho_mu(DoubleTab& val_moy) const
{
  const Domaine_dis_base& zdisbase=mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF& domaine_VDF=ref_cast(Domaine_VDF, zdisbase);
  //  const DoubleTab& xp = domaine_VDF.xp();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  double u,v,wl;
  int nb_elems = domaine_VDF.domaine().nb_elem();
  int num_elem,i;
  int face_x_0,face_y_0,face_y_1,face_z_0;

  const Fluide_base& le_fluide = ref_cast(Fluide_base,mon_equation->milieu());
  const DoubleTab& visco_dyn = le_fluide.viscosite_dynamique().valeurs();
  const DoubleTab& tab_rho_elem = le_fluide.masse_volumique().valeurs();
  int taille_mu=visco_dyn.dimension(0);
  int taille_rho=tab_rho_elem.dimension(0);

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      //y=xp(num_elem,1);

      face_x_0 = elem_faces(num_elem,0);
      //      face_x_1 = elem_faces(num_elem,dimension);
      face_y_0 = elem_faces(num_elem,1);
      face_y_1 = elem_faces(num_elem,1+dimension);

      // PQ : 12/10 : to avoid averaging the velocity u and w locally,
      //              we choose to "shift" them
      //              from the face center to the element center
      //              based on the assumption that the flow is homogeneous
      //                    in the xz planes.
      //                    For v, the average is necessary to return to the element center.

      //  u = .5*(vitesse[face_x_0]+vitesse[face_x_1]);
      u = vitesse[face_x_0];
      v = .5*(vitesse[face_y_0]+vitesse[face_y_1]);

      i= Tab_recap(num_elem,2);

      val_moy(i,0) += u;
      val_moy(i,1) += v;
      val_moy(i,3) += u*u;
      val_moy(i,4) += v*v;
      val_moy(i,6) += u*v;

      if(dimension==2)   val_moy(i,9) += sqrt(u*u);      //tangential velocity for friction computation

      if(dimension==3)
        {
          face_z_0 = elem_faces(num_elem,2);
          //  face_z_1 = elem_faces(num_elem,2+dimension);

          //     w = .5*(vitesse[face_z_0]+vitesse[face_z_1]);
          wl = vitesse[face_z_0];

          val_moy(i,2) += wl;
          val_moy(i,5) += wl*wl;
          val_moy(i,7) += u*wl;
          val_moy(i,8) += v*wl;
          val_moy(i,9) += sqrt(u*u+wl*wl);      //tangential velocity for friction computation
        }

      if (taille_rho==1)  val_moy(i,10) += tab_rho_elem(0,0);
      else                val_moy(i,10) += tab_rho_elem[num_elem];


      if (taille_mu==1)   val_moy(i,11) += visco_dyn(0,0);
      else                val_moy(i,11) += visco_dyn[num_elem];
    }
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_nut(DoubleTab& val_moy) const
{
  const Domaine_dis_base& zdisbase=mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF& domaine_VDF=ref_cast(Domaine_VDF, zdisbase);
  //const DoubleTab& xp = domaine_VDF.xp();
  const RefObjU& modele_turbulence = mon_equation->get_modele(TURBULENCE);
  const Modele_turbulence_hyd_base& mod_turb = ref_cast(Modele_turbulence_hyd_base,modele_turbulence.valeur());
  const DoubleTab& nu_t = mod_turb.viscosite_turbulente().valeurs();

  int nb_elems = domaine_VDF.domaine().nb_elem();
  int num_elem,i;
  // double y;

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      //y=xp(num_elem,1);

      i= Tab_recap(num_elem,2);

      val_moy(i,12) += nu_t[num_elem];
    }
}

void Traitement_particulier_NS_canal_VDF::calculer_moyenne_spatiale_Temp(DoubleTab& val_moy) const
{
  const Domaine_dis_base& zdisbase=mon_equation->inconnue().domaine_dis_base();
  const Domaine_VDF& domaine_VDF=ref_cast(Domaine_VDF, zdisbase);
  //const DoubleTab& xp = domaine_VDF.xp();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  const DoubleTab& temperature = Temp->valeurs();
  const DoubleTab& vitesse = mon_equation->inconnue().valeurs();
  double u,v,wl,T;
  int nb_elems = domaine_VDF.domaine().nb_elem();
  int num_elem,i;
  int face_x_0,face_x_1,face_y_0,face_y_1,face_z_0,face_z_1;

  for (num_elem=0; num_elem<nb_elems; num_elem++)
    {
      //y=xp(num_elem,1);

      T = temperature[num_elem];

      face_x_0 = elem_faces(num_elem,0);
      face_x_1 = elem_faces(num_elem,dimension);
      face_y_0 = elem_faces(num_elem,1);
      face_y_1 = elem_faces(num_elem,1+dimension);

      u = .5*(vitesse[face_x_0]+vitesse[face_x_1]);
      v = .5*(vitesse[face_y_0]+vitesse[face_y_1]);


      i=Tab_recap(num_elem,2);
      // for (i=0; i<Ny; i++)     if(est_egal(y,Y[i])) break;

      val_moy(i,13) += T;
      val_moy(i,14) += T*T;
      val_moy(i,15) += u*T;
      val_moy(i,16) += v*T;

      if(dimension==3)
        {
          face_z_0 = elem_faces(num_elem,2);
          face_z_1 = elem_faces(num_elem,2+dimension);

          wl = .5*(vitesse[face_z_0]+vitesse[face_z_1]);

          val_moy(i,17) += wl*T;
        }
    }
}

