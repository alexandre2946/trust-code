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

#include <Op_Conv_Centre_EF_VEF_Face.h>
#include <Neumann_sortie_libre.h>
#include <Periodique.h>

Implemente_instanciable(Op_Conv_Centre_EF_VEF_Face,"Op_Conv_Centre_EF_VEF_P1NC",Op_Conv_VEF_base);
// XD bloc_ef objet_lecture nul NO_BRACE not_set
// XD attr mot1 chaine(into=["transportant_bar","transporte_bar","filtrer_resu","antisym"]) mot1 REQ not_set
// XD attr val1 entier(into=[0,1]) val1 REQ not_set
// XD attr mot2 chaine(into=["transportant_bar","transporte_bar","filtrer_resu","antisym"]) mot2 REQ not_set
// XD attr val2 entier(into=[0,1]) val2 REQ not_set
// XD attr mot3 chaine(into=["transportant_bar","transporte_bar","filtrer_resu","antisym"]) mot3 REQ not_set
// XD attr val3 entier(into=[0,1]) val3 REQ not_set
// XD attr mot4 chaine(into=["transportant_bar","transporte_bar","filtrer_resu","antisym"]) mot4 REQ not_set
// XD attr val4 entier(into=[0,1]) val4 REQ not_set
// XD convection_ef convection_deriv ef NO_BRACE For VEF calculations, a centred convective scheme based on Finite
// XD_CONT Elements formulation can be called through the following data:NL2 NL2 Convection { EF transportant_bar val
// XD_CONT transporte_bar val antisym val filtrer_resu val }NL2 NL2 This scheme is 2nd order accuracy (and get better
// XD_CONT the property of kinetic energy conservation). Due to possible problems of instabilities phenomena, this
// XD_CONT scheme has to be coupled with stabilisation process (see Source_Qdm_lambdaup).These two last data are
// XD_CONT equivalent from a theoretical point of view in variationnal writing to : div(( u. grad ub , vb) - (u. grad
// XD_CONT vb, ub)), where vb corresponds to the filtered reference test functions.NL2 NL2 Remark:NL2 This class
// XD_CONT requires to define a filtering operator : see solveur_bar
// XD attr mot1 chaine(into=["defaut_bar"]) mot1 OPT equivalent to transportant_bar 0 transporte_bar 1 filtrer_resu 1
// XD_CONT antisym 1
// XD attr bloc_ef bloc_ef bloc_ef OPT not_set

Sortie& Op_Conv_Centre_EF_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Entree& Op_Conv_Centre_EF_VEF_Face::readOn(Entree& s )
{
  return s ;
}

// WARNING!!!!!!! modifications regarding fluent_ and autre_num_face_loc are made only in 3D!!!!
// C.A. 30/06/99

////////////////////////////////////////////////////////////////////
//
//                      Implementation of functions
//
//                   of class Op_Conv_Centre_EF_VEF_Face
//
////////////////////////////////////////////////////////////////////


DoubleTab& Op_Conv_Centre_EF_VEF_Face::ajouter(const DoubleTab& transporte,
                                               DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  const Champ_Inc_base& la_vitesse=vitesse_.valeur();

  const DoubleVect& porosite_face = equation().milieu().porosite_face();

  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& face_normales = domaine_VEF.face_normales();
  const auto& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();


  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();

  int nfac = domaine.nb_faces_elem();
  //int nsom = domaine.nb_som_elem();


  // For the convection treatment, standard polyhedra (not "seeing" boundary conditions)
  // are distinguished from non-standard polyhedra (having at least one boundary face).
  // A standard polyhedron has n facets on which the convection scheme is applied.
  // For a non-standard polyhedron with Dirichlet boundary conditions, part of its
  // facets are carried by the boundary faces.
  // In short, for a polyhedron the convection treatment depends on the type
  // (triangle, tetrahedron ...) and the number of Dirichlet faces.

  double flux;
  int poly,face_adj,fa7,i,j,comp0,n_bord;
  int num_face, rang ,itypcl;
  int num10, num20;

  int ncomp_ch_transporte;
  if (transporte.nb_dim() == 1)
    ncomp_ch_transporte=1;
  else
    ncomp_ch_transporte= transporte.dimension(1);

  IntVect face(nfac);
  DoubleVect cc(dimension);

  //////////////////////////////
  DoubleTab psc(nfac);
  int num_int;
  IntTab autre_num_face(dimension-1);
  IntTab autre_num_face_loc(dimension-1);
  double coef1=0.,coef2=0.,coef3=0.;
  //  double psc1;
  int nu1,nu2;
  //  int k;
  double f_int;

  //DoubleVect vs(dimension);
  //DoubleVect vc(dimension);
  //DoubleTab vsom(nsom,dimension);

  // Reset to zero the array used for
  // the stability time step computation
  fluent_ = 0;

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


  // Non-standard polyhedra are arranged in 2 groups in Domaine_VEF:
  //  - boundary and joint polyhedra
  //  - boundary and non-joint polyhedra
  // Polyhedra are processed following the order in which they appear
  // in the domain

  // loop over polyhedra
  const IntTab& KEL=domaine_VEF.type_elem().KEL();
  for (poly=0; poly<nb_elem_tot; poly++)
    {

      rang = rang_elem_non_std(poly);
      if (rang==-1)
        itypcl=0;
      else
        itypcl=domaine_Cl_VEF.type_elem_Cl(rang);

      // compute the face indices of the polyhedron
      for (face_adj=0; face_adj<nfac; face_adj++)
        face[face_adj]= elem_faces(poly,face_adj);

      // Find the local indices of all faces
      for (fa7=0; fa7<nfa7; fa7++)
        {
          nu1=-1;
          nu2=-1;
          num10 = face[KEL(0,fa7)];
          num20 = face[KEL(1,fa7)];
          // The facet is surrounded by faces num1 and num2
          //        Cerr << "num1=" << num1 << "  num2=" << num2 << finl;

          // Find the indices of the other faces

          i=0;
          j=0;
          //                k=0;
          while(i<nfac)
            {
              num_int = face[i];
              if (num_int == num10)
                {
                  nu1=i;
                  //                Cerr << "nu1 (in loops)=" << nu1 << finl;
                }
              else if (num_int == num20)
                {
                  nu2=i;
                  //                  Cerr << "nu2 (in loops)=" << nu2 << finl;
                }
              else
                {
                  autre_num_face_loc(j)=i;
                  autre_num_face(j)=num_int;
                  //                  Cerr << "autre_num_face (in loops)=" << autre_num_face(j) << finl;
                  j++;
                  //                          k++;
                }
              i++;
            }

          if (rang==-1)
            {
              for (i=0; i<dimension; i++)
                cc[i] = facette_normales(poly,fa7,i);
            }
          else
            for (i=0; i<dimension; i++)
              cc[i] = normales_facettes_Cl(rang,fa7,i);

          // Compute the dot products u(xi).n.S
          for (i=0; i<nfac; i ++)
            {
              psc[i] = 0.;
              for (j=0; j<dimension; j++)
                {
                  //                   Cerr << "cc[j]=" << cc[j] << finl;
                  //                   Cerr << "la_vitesse.valeurs()(face[i],j)=" << la_vitesse.valeurs()(face[i],j) << finl;

                  psc[i]+= la_vitesse.valeurs()(face[i],j)*cc[j]*porosite_face(face[i]);
                }
              //            Cerr << "psc[" << i << "]=" <<  psc[i] << finl;
            }

          //                assert(ncomp_ch_transporte==dimension);
          // This scheme is valid (at least I think...) only when ch_transporte = ch_tranportant = velocity??


          if (dimension == 2)
            {
              switch(itypcl)
                {
                case 0:
                  {
                    coef1  = 13.0*( psc[nu1]+psc[nu2] ) ;
                    coef1 -=  8.0*psc[autre_num_face_loc(0)];

                    coef2  =  8.0*( psc[nu1]+psc[nu2] ) ;
                    coef2 -=  7.0*psc[autre_num_face_loc(0)];

                    if (ncomp_ch_transporte == 1)
                      {
                        flux  = (transporte(num10)+transporte(num20))*coef1;
                        flux -=  transporte(autre_num_face(0))*coef2;
                        flux /= 27.;
                        resu(num10) -= flux;
                      }
                    else
                      {
                        for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                          {
                            flux  = (transporte(num10,comp0)+transporte(num20,comp0))*coef1;
                            flux -=  transporte(autre_num_face(0),comp0)*coef2;
                            flux /= 27.;
                            resu(num10, comp0) -= flux;
                            resu(num20, comp0) += flux;
                          }
                      }
                    // For the stability time step computation
                    fluent_[num10] += 2.*((psc[nu1]+psc[nu2])- psc[autre_num_face_loc(0)])/3.;
                    fluent_[num10] -= 2.*((psc[nu1]+psc[nu2])- psc[autre_num_face_loc(0)])/3.;
                    break;
                  }
                default :
                  {
                    int numfa7;
                    //  !!!!!! only the case where the transported field is a vector has been handled!!!
                    if ((itypcl==1)||(itypcl==2)||(itypcl==4))     // 1 Dirichlet face!!
                      {
                        switch(itypcl)
                          {
                          case 1 :
                            {
                              numfa7 = 2;
                              break;
                            }
                          case 2 :
                            {
                              numfa7 = 1;
                              break;
                            }
                          case 4 :
                            {
                              numfa7 = 0;
                              break;
                            }
                          default :
                            {
                              numfa7=-1;
                              Cerr << "This should not be possible!!!" << finl;
                              exit();
                              break;
                            }
                          }

                        //                    if (fa7 == itypcl)
                        if (fa7 == numfa7)  // We are on the fa7 not coinciding with the Dirichlet face
                          {
                            coef1  = 2.*( psc[nu1]+psc[nu2] ) ;
                            coef1 -= psc[numfa7];

                            coef2  =   psc[nu1]+psc[nu2] ;
                            coef2 -=  2.*psc[numfa7];

                            for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                              {
                                flux  = (transporte(num10,comp0)+transporte(num20,comp0))*coef1;
                                flux -=  transporte(numfa7,comp0)*coef2;
                                flux /= 6.;
                                resu(num10, comp0) -= flux;
                                resu(num20, comp0) += flux;
                              }
                            // For the stability time step computation
                            fluent_[num10] += 0.5*(psc[nu1]+psc[nu2]);
                            fluent_[num20] -= 0.5*(psc[nu1]+psc[nu2]);
                            //                            fluent_[num1] = ( fluent_[num1] > f_int) ? fluent_[num1] : f_int ;
                          }
                        else
                          {
                            // For coinciding fa7
                            if (fa7 == nu2)
                              {
                                coef1  = 2.*( psc[nu1]-psc[nu2] ) ;
                                coef1 += 3.*psc[numfa7];

                                coef2  =  3.*(psc[nu1]-psc[nu2]) ;
                                coef2 +=  6.*psc[numfa7];


                                for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                                  {
                                    flux  = (transporte(face[nu1],comp0)-transporte(face[nu2],comp0))*coef1;
                                    flux +=  transporte(face[numfa7],comp0)*coef2;
                                    flux /= 6.;
                                    resu(face[nu1],comp0) -= flux;  // is this the right sign??????????????????????
                                    // Cerr << "num1=" << num1 << "  num2=" << num2 << finl;
                                  }
                                // For the stability time step computation
                                fluent_[face[nu1]] -=  0.5*(psc[nu1]-psc[nu2])+psc[numfa7]; // sign????????

                                //////////////////////////////////////////////
                              }
                            else
                              {
                                // fa7 == nu1
                                coef1  = 2.*( psc[nu2]-psc[nu1] ) ;
                                coef1 += 3.*psc[numfa7];

                                coef2  =  3.*(psc[nu2]-psc[nu1]) ;
                                coef2 +=  6.*psc[numfa7];  // is this the right numbering?


                                for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                                  {
                                    flux  = (transporte(face[nu2],comp0)-transporte(face[nu1],comp0))*coef1;
                                    flux +=  transporte(face[numfa7],comp0)*coef2;
                                    flux /= 6.;
                                    resu(face[nu1],comp0) -= flux;  // is this the right sign?????????????????
                                    // Cerr << "num1=" << num1 << "  num2=" << num2 << finl;
                                  }
                                // For the stability time step computation
                                fluent_[face[nu1]] -=  0.5*(psc[nu2]-psc[nu1])+psc[numfa7];   // sign?????????????????

                                //////////////////////////////////////////////
                              }
                          }
                      }
                    else
                      {
                        // 2 Dirichlet faces!!!!!
                        // same expression for the 3 fa7 (coinciding with boundary faces)
                        // Assuming num fa7 = num face

                        // Retrieve the indices of the other faces
                        switch(fa7)
                          {
                          case 0 :
                            {
                              nu1 = 1;
                              nu2 = 2;
                              break;
                            }
                          case 1 :
                            {
                              nu1 = 0;
                              nu2 = 2;
                              break;
                            }
                          case 2 :
                            {
                              nu1 = 1;
                              nu2 = 0;
                              break;
                            }
                          default :
                            {
                              nu1=-1;
                              nu2=-1;
                              Cerr << "Stopping everything, this should not be possible!!!!" << finl;
                              Cerr << "otherwise it means I have misunderstood something!!!" << finl;
                              exit();
                            }
                          }
                        // Cerr << "fa7=" << fa7 << "   nu1=" << nu1 << "   nu2=" << nu2 << finl;
                        coef1  = psc[fa7];
                        coef2  = (psc[nu1]-psc[nu2])/3.;
                        for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                          {
                            flux  = transporte(face[fa7],comp0)*coef1;
                            flux += (transporte(face[nu1],comp0)-transporte(face[nu2],comp0))*coef2;
                            flux /= 6.;
                            resu(face[nu1],comp0) -= flux;  // is this the right sign?
                          }
                        // For the stability time step computation
                        fluent_[face[fa7]] -=  psc[fa7];
                        //                             f_int =  psc[fa7];
                        //                             fluent_[face[fa7]] = ( fluent_[face[fa7]] > f_int) ? fluent_[face[fa7]]  : f_int ;
                        //////////////////////////////////////////////
                      }
                    break;
                  }
                }
            }          // END of if(dimension == 2)
          else if (dimension == 3)
            {
              switch(itypcl)
                {
                case 0:
                  {
                    coef1  =  19.0*( psc[nu1]+psc[nu2] ) ;
                    coef1 -=   7.0*(psc[autre_num_face_loc(0)]+psc[autre_num_face_loc(1)]);
                    //                      Cerr << "coef1=" << coef1 << finl;

                    coef2  =   7.0*( psc[nu1]+psc[nu2] ) ;
                    coef2 -=  15.0*psc[autre_num_face_loc(0)];
                    coef2 +=   9.0*psc[autre_num_face_loc(1)];
                    //                      Cerr << "coef2=" << coef2 << finl;

                    coef3  =   7.0*( psc[nu1]+psc[nu2] ) ;
                    coef3 +=   9.0*psc[autre_num_face_loc(0)];
                    coef3 -=  15.0*psc[autre_num_face_loc(1)];
                    //                      Cerr << "coef3=" << coef3 << finl;

                    for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                      {
                        flux  = (transporte(num10,comp0)+transporte(num20,comp0))*coef1;
                        flux -=  transporte(autre_num_face(0),comp0)*coef2;
                        flux -=  transporte(autre_num_face(1),comp0)*coef3;
                        flux /= 32.;
                        resu(num10, comp0) -= flux;
                        resu(num20, comp0) += flux;
                        // Cerr << "transporte(num1,comp)=" << transporte(num1,comp) << finl;
                        //                           Cerr << "transporte(num2,comp)=" << transporte(num2,comp) << finl;
                        //                           Cerr << "transporte(autre_num_face(0),comp)=" << transporte(autre_num_face(0),comp) << finl;
                        //                           Cerr << "transporte(autre_num_face(1),comp)=" << transporte(autre_num_face(1),comp) << finl;
                        //                           Cerr << "flux=" << flux << finl;
                        //************Compute fluent for the time step
                        // f_int = 0.5*cc[comp]*(la_vitesse.valeurs()(num1,comp)+la_vitesse.valeurs()(num2,comp));

                        //  if (f_int>=0.)
                        //                             fluent_[num2] += f_int ;
                        //                           else
                        //                             fluent_[num1] -= f_int ;
                      }
                    f_int = 3.*(psc[nu1]+psc[nu2]);
                    f_int -= (psc[autre_num_face_loc(0)]+psc[autre_num_face_loc(1)]);
                    f_int /= 4.;
                    if (f_int>=0.)
                      {
                        //fluent_[num2] += f_int ;
                        fluent_[num20] = ( fluent_[num20] > std::fabs(f_int))? fluent_[num20] : std::fabs(f_int);
                      }
                    else
                      {
                        //fluent_[num1] -= f_int ;
                        fluent_[num10] = ( fluent_[num10] > std::fabs(f_int))? fluent_[num10] : std::fabs(f_int);
                      }



                    break;
                  }
                default:
                  {
                    Cerr << "Dans le default!!" << finl;
                    int nu3=-1;
                    if (fa7 == itypcl)
                      {
                        i=0;
                        while(i<nfac) // search for the 4th point!!
                          {
                            nu3 = i;
                            if (nu3 == fa7)
                              i++;
                            else if (nu3 == nu1)
                              i++;
                            else if (nu3 == nu2)
                              i++;
                            else
                              i=nfac;
                          }
                        coef1  =  6.*( psc[nu1]+psc[nu2] ) ;
                        coef1 -=  3.*psc[nu3]+psc[fa7];

                        coef2  =   3.0*( psc[nu1]+psc[nu2] ) ;
                        coef2 -=   7.0*psc[nu3];
                        coef2 +=   4.0*psc[fa7];

                        coef3  =   psc[nu1]+psc[nu2];
                        coef3 +=   psc[nu3];
                        coef3 -=   7.*psc[fa7];

                        for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                          {
                            flux  = (transporte(num10,comp0)+transporte(num20,comp0))*coef1;
                            flux -=  transporte(face[nu3],comp0)*coef2;
                            flux -=  transporte(face[fa7],comp0)*coef3;
                            flux /= 12.;
                            resu(num10, comp0) -= flux;
                            resu(num20, comp0) += flux;
                            // For the stability time step computation
                            // SIGNE????????????
                            //                                  f_int = std::fabs(flux/transporte(num1,comp));
                            f_int = 3.*(psc[nu1]+psc[nu2]);
                            f_int -= (psc[autre_num_face_loc(0)]+psc[autre_num_face_loc(1)]);
                            f_int /= 4.;
                            // fluent_[num1] -= f_int ;
                            //                               fluent_[num2] += f_int ;
                            if (f_int >=0.)
                              fluent_[num20] += f_int ;
                            else
                              fluent_[num10] -= f_int ;

                            //                              fluent_[num1] = ( fluent_[num1] > f_int) ? fluent_[num1]  : f_int ;
                            //////////////////////////////////////////////
                            //                          Cerr << "flux=" << flux << finl;
                          }
                      }
                    else
                      {
                        // This appears to be unfinished???!!!!!!!!!!
                      }
                    break;
                  }
                }
            }
        }
    }
  //  Cerr << "vitesse=" << la_vitesse.valeurs() << finl;


  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;
  double pscav;

  // Sizing the array of convective fluxes at the domain boundary
  DoubleTab& flux_b = flux_bords_;
  flux_b.resize(domaine_VEF.nb_faces_bord(),ncomp_ch_transporte);
  flux_b = 0.;

  // Loop over boundaries to handle boundary conditions
  // a convection term is taken into account only for
  // Neumann_sortie_libre boundary conditions

  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {

      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          ////////////WARNING!!!!!!!!!!! This corresponds to the old schema and not to EF!!
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              pscav =0;
              for (i=0; i<dimension; i++)
                pscav += la_vitesse.valeurs()(num_face,i)*face_normales(num_face,i)*porosite_face(num_face);
              if (pscav>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= pscav*transporte(num_face);
                    flux_b(num_face,0) -= pscav*transporte(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= pscav*transporte(num_face,i);
                      flux_b(num_face,i) -= pscav*transporte(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= pscav*la_sortie_libre.val_ext(num_face-num1);
                      flux_b(num_face,0) -= pscav*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= pscav*la_sortie_libre.val_ext(num_face-num1,i);
                        flux_b(num_face,i) -= pscav*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                  fluent_[num_face] -= pscav;
                }
            }
          //        Cerr << "For now Neumann_sortie_libre not supported!!!" << finl;
        }
      else if (sub_type(Periodique,la_cl.valeur()))
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
                      flux_b(num_face,0) += diff1;
                      // For the stability time step computation
                      // NOTHING in periodic ??? (faces already processed before??????)
                      //////////////////////////////////////////////
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        //                       Cerr << "num_face=" << num_face << "  diff1=" << diff1 << finl;
                        //                       Cerr << "voisine=" << voisine << "  diff2=" << diff2 << finl;
                        resu(voisine,comp)  += diff1;
                        resu(num_face,comp) += diff2;
                        flux_b(voisine,comp) += diff1;
                        flux_b(num_face,comp) += diff1;
                        // For the stability time step computation
                        // NOTHING in periodic ??? (faces already processed before??????)
                        //////////////////////////////////////////////
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }
  modifier_flux(*this);
  return resu;

}
