/****************************************************************************
* Copyright (c) 2026, CEA
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

#include <Op_Conv_Vort_VEF_Face.h>
#include <Champ_P1NC.h>
#include <Periodique.h>

Implemente_instanciable(Op_Conv_Vort_VEF_Face,"Op_Conv_Conserve_Ec_VEF_P1NC",Op_Conv_VEF_base);


//// printOn
//

Sortie& Op_Conv_Vort_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Op_Conv_Vort_VEF_Face::readOn(Entree& s )
{
  return s ;
}


DoubleTab& Op_Conv_Vort_VEF_Face::ajouter(const DoubleTab& transporte,
                                          DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  const Champ_Inc_base& la_vitesse=vitesse_.valeur();
  const Champ_P1NC& vit = ref_cast(Champ_P1NC,vitesse_.valeur());

  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const auto& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();
  const int nb_faces = domaine_VEF.nb_faces();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem = domaine_VEF.nb_elem();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();


  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();

  int nfac = domaine.nb_faces_elem();

  const DoubleVect& volumes = domaine_VEF.volumes();
  int comp0;
  double flux;//,flux_int;
  int num_face;
  int elem0,elem1;
  double vol0,vol1;
  double inter,a0,a1,a2,f_int;

  IntVect face(nfac);
  DoubleVect cc(dimension);
  DoubleTab psc(nfac);

  int num_int;
  IntTab autre_num_face(dimension-1);
  IntTab autre_num_face_loc(dimension-1);
  int poly,face_adj,fa7,i,j,n_bord;
  int rang ;
  int num10, num20;
  int nu1, nu2;
  int num_calc;


  // For the convection treatment, standard polyhedra (not "seeing" boundary conditions)
  // are distinguished from non-standard polyhedra (having at least one boundary face).
  // A standard polyhedron has n facets on which the convection scheme is applied.
  // For a non-standard polyhedron with Dirichlet boundary conditions, part of its
  // facets are carried by the boundary faces.
  // In short, for a polyhedron the convection treatment depends on the type
  // (triangle, tetrahedron ...) and the number of Dirichlet faces.

  int ncomp_ch_transporte;
  if (transporte.nb_dim() == 1)
    ncomp_ch_transporte=1;
  else
    ncomp_ch_transporte= transporte.dimension(1);

  //  Cerr << "ncomp_ch_transporte=" << ncomp_ch_transporte << finl;

  // Reset the array used for
  // computing the stability time step
  fluent_ = 0;

  // WARNING: issue with determining the flux (fluent)
  // ********  set to 1 for now!!!
  //  fluent_ = 1.;


  // Special treatment for periodic faces

  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            nb_faces_perio++;
        }
    }

  DoubleTab tab;
  if (ncomp_ch_transporte == 1)
    tab.resize(nb_faces_perio);
  else
    tab.resize(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          //          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (ncomp_ch_transporte == 1)
                tab(nb_faces_perio) = resu(num_face);
              else
                for (int comp=0; comp<ncomp_ch_transporte; comp++)
                  tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }

  //  Cerr << "tab=" << tab << finl;

  // Non-standard polyhedra are grouped into 2 sets in Domaine_VEF:
  //  - boundary and joint polyhedra
  //  - boundary and non-joint polyhedra
  // Polyhedra are processed in the order in which they appear in the domain

  // loop over polyhedra


  // Loop to add the part: Gradient(U^2/2)
  // ******* loop over elements
  // 06/01/2000 Boundary conditions not yet handled here (except periodic)
  const IntTab& KEL=domaine_VEF.type_elem().KEL();
  for (poly=0; poly<nb_elem_tot; poly++)
    {

      rang = rang_elem_non_std(poly);

      // compute the face indices of the polyhedron
      for (face_adj=0; face_adj<nfac; face_adj++)
        face[face_adj]= elem_faces(poly,face_adj);

      // Find the global indices of all faces
      for (fa7=0; fa7<nfa7; fa7++)
        {
          nu1=-1;
          nu2=-1;
          num10 = face[KEL(0,fa7)];
          num20 = face[KEL(1,fa7)];
          // The facet is surrounded by faces num1 and num2
          //        Cerr << "num1=" << num1 << "  num2=" << num2 << finl;

          //           i=0;
          //           j=0;
          //           while(i<nfac)
          //             {
          //               num_int = face[i];
          //               if ((num_int!= num1)&&(num_int!= num2))
          //                 {
          //                   autre_num_face(j)=num_int;
          //                   j++;
          //                 }
          //               i++;
          //             }

          // Find the indices of the other faces (local and global)

          i=0;
          j=0;
          //                k=0;
          while(i<nfac)
            {
              num_int = face[i];
              if (num_int == num10)
                {
                  nu1=i;
                }
              else if (num_int == num20)
                {
                  nu2=i;
                }
              else
                {
                  autre_num_face_loc(j)=i;
                  autre_num_face(j)=num_int;
                  j++;
                }
              i++;
            }

          //           Cerr << "num1=" << num1 << "  num2=" << num2 << "  autre_num_face(0)=" << autre_num_face(0) << finl;
          //           if (dimension==3)
          //             Cerr << "autre_num_face(1)=" << autre_num_face(1) << finl;

          if (rang==-1)
            {
              for (i=0; i<dimension; i++)
                cc[i] = facette_normales(poly,fa7,i);
            }
          else
            for (i=0; i<dimension; i++)
              cc[i] = normales_facettes_Cl(rang,fa7,i);

          // Compute the scalar products u(xi).n.S  // >>> fluent computation!!
          for (i=0; i<nfac; i ++)
            {
              psc[i] = 0.;
              for (j=0; j<dimension; j++)
                {
                  psc[i]+= la_vitesse.valeurs()(face[i],j)*cc[j];
                }
            }

          // Compute the flux
          // Loop over components: uu+vv+(ww)
          flux = 0.;
          if (dimension == 2)
            {
              f_int =  2.*((psc[nu1]+psc[nu2])- psc[autre_num_face_loc(0)])/3.;
            }
          else
            {
              // (dimension == 3)
              assert(dimension == 3);
              {
                f_int = 3.*(psc[nu1]+psc[nu2]);
                f_int -= (psc[autre_num_face_loc(0)]+psc[autre_num_face_loc(1)]);
                f_int /= 4.;
              }
            }
          if (f_int >= 0.)
            num_calc = num10;
          else
            num_calc = num20;

          flux = 0.;
          for (comp0=0; comp0<dimension; comp0++)
            flux += la_vitesse.valeurs()(num_calc,comp0)*la_vitesse.valeurs()(num_calc,comp0);

          for (comp0=0; comp0<dimension; comp0++)
            {
              resu(num10, comp0) -= 0.5*flux*cc[comp0];
              resu(num20, comp0) += 0.5*flux*cc[comp0];
            }

          // *** ??? : flux (fluent) evaluation
          if (f_int>0.)
            {
              // fluent_[num2] += std::fabs(f_int);
              fluent_[num20] = ( fluent_[num20] > std::fabs(f_int))? fluent_[num20] : std::fabs(f_int);
            }
          else
            {
              fluent_[num10] = ( fluent_[num10] > std::fabs(f_int))? fluent_[num10] : std::fabs(f_int);
              // fluent_[num1] += std::fabs(f_int);
            }

        }
    }

  // FIN DE LA BOUCLE SUR LES ELEMENTS
  ////////// Apply compensation here, because the next loop is over faces.
  ////////// Doing it at the end would count the face contribution twice.
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Dimensioning the array of convective fluxes at the boundary
  // of the computational domain
  DoubleTab& flux_b = flux_bords_;
  flux_b.resize(domaine_VEF.nb_faces_bord(),ncomp_ch_transporte);
  flux_b = 0.;

  // Loop over the boundaries to process the boundary conditions

  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          IntVect fait(le_bord.nb_faces());
          fait = 0;
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (fait[num_face-num1] == 0)
                {
                  voisine = la_cl_perio.face_associee(num_face-num1) + num1;

                  if (ncomp_ch_transporte == 1)
                    {
                      diff1 = resu(num_face)-tab(nb_faces_perio);
                      diff2 = resu(voisine)-tab(nb_faces_perio+voisine-num_face);
                      resu(voisine)  += diff1;
                      resu(num_face) += diff2;
                      flux_b(voisine,0) += diff1;
                      flux_b(num_face,0) += diff2;
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  += diff1;
                        resu(num_face,comp) += diff2;
                        flux_b(voisine,comp) += diff1;
                        flux_b(num_face,comp) += diff2;
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }

  /////////////////////////////////////////////////////
  // Loop to add the vorticity part
  // ****** Loop over faces

  // Compute the vorticity
  DoubleTab vorticite;
  if (dimension == 2)
    vorticite.resize(nb_elem);
  else if (dimension == 3)
    vorticite.resize(nb_elem,dimension);

  vit.cal_rot_ordre1(vorticite);

  //  Cerr << "vorticite=" << vorticite << finl;

  for (num_face=0; num_face<nb_faces; num_face++)
    {
      vol0=-1;
      vol1=-1;
      elem0 = face_voisins(num_face,0);
      elem1 = face_voisins(num_face,1);

      if (elem0 != -1)
        vol0 = volumes(elem0);

      if (elem1 != -1)
        vol1 = volumes(elem1);

      //      Cerr << "vol0=" << vol0 << "  vol1=" << vol1 << finl;

      if (dimension == 2)
        {
          //           for (comp=0;comp<dimension;comp++)
          //             {
          assert(vol0>0);
          assert(vol1>0);
          inter  = vorticite[elem0]*vol0/3.+vorticite[elem1]*vol1/3.;

          resu(num_face,0) -= -inter*la_vitesse.valeurs()(num_face,1);
          resu(num_face,1) -= inter*la_vitesse.valeurs()(num_face,0);

          // minus sign because we are in the right-hand side

          // *** PBL : flux (fluent) evaluation
          //               if(psc >= 0)
          //                 fluent_[num2] += psc;
          //               else
          //                 fluent_[num1] -= psc;
        }
      else if (dimension == 3)
        {
          assert(vol0>0);
          assert(vol1>0);
          // vect(a) = vorticite*Vol
          a0 = vorticite(elem0,0)*vol0/4. + vorticite(elem1,0)*vol1/4.;
          a1 = vorticite(elem0,1)*vol0/4. + vorticite(elem1,1)*vol1/4.;
          a2 = vorticite(elem0,2)*vol0/4. + vorticite(elem1,2)*vol1/4.;

          resu(num_face,0) -= a1*la_vitesse.valeurs()(num_face,2)-a2*la_vitesse.valeurs()(num_face,1) ;
          resu(num_face,1) -= a2*la_vitesse.valeurs()(num_face,0)-a0*la_vitesse.valeurs()(num_face,2) ;
          resu(num_face,2) -= a0*la_vitesse.valeurs()(num_face,1)-a1*la_vitesse.valeurs()(num_face,0) ;

          // minus sign because we are in the right-hand side

          // *** PBL : flux (fluent) evaluation
          //               if(psc >= 0)
          //                 fluent_[num2] += psc;
          //               else
          //                 fluent_[num1] -= psc;
        }
    }

  //******* VERIF PERIO
  Cerr << "DEBUT VERIF PERIO" << finl;
  //  Cerr << "nb_front_Cl=" << domaine_VEF.nb_front_Cl() << finl;
  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          //          Cerr << "num1=" << num1 << "  num2=" << num2 << finl;
          for (num_face=num1; num_face<num2; num_face++)
            {
              voisine = la_cl_perio.face_associee(num_face-num1) + num1;
              for (int ii=0; ii<dimension; ii++)
                {
                  if ( resu(num_face,ii)!=resu(voisine,ii) )
                    {
                      Cerr << "Pbl de periodicite a la face" << num_face << finl;
                      Cerr << "diff = " << resu(num_face,ii)-resu(voisine,ii)  << finl;
                    }
                }
            }
        }
    }
  Cerr << "FIN VERIF PERIO" << finl;
  //******* FIN VERIF PERIO

  modifier_flux(*this);
  return resu;
}
