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

#include <Op_Conv_kschemas_VEF.h>
#include <Periodique.h>
#include <Neumann_sortie_libre.h>

Implemente_base(Op_Conv_kschemas_VEF,"Op_Conv_kschemas_VEF_P1NC",Op_Conv_VEF_base);


Sortie& Op_Conv_kschemas_VEF::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Entree& Op_Conv_kschemas_VEF::readOn(Entree& s )
{
  return s ;
}

void Op_Conv_kschemas_VEF::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_cl_dis, const Champ_Inc_base& ch)
{
  // CCa le 28/05/99 Le schema Kquick ne marche pas en paralle !!
  if (Process::is_parallel())
    {
      Cerr << "ATTENTION le kquick ne marche pas en parallele !!!" << finl;
      exit();
    }

  Op_Conv_VEF_base::associer(domaine_dis, domaine_cl_dis, ch);
}


//////////////////////////////////////////////////////////////
//   Functions for the k-schemes.
////////////////////////////////////////////////////////////////

// convkschemas : utility function for convection

void convkschemas(const double K, const int ncomp, int dimension, const int poly ,
                  const int poly1, const int poly2,const int jel0,
                  const int jel1,const double psc ,const DoubleTab& tab1 ,
                  DoubleVect& fluent, DoubleVect& flux,
                  const DoubleVect& rx0, const DoubleTab& gradient_elem )
{

  int comp,amont,i,elem1,elem2;
  double CF,UTC,deltat0,deltat1,deltat;
  DoubleVect rx(dimension);
  deltat0 = 0.;
  deltat1 = 0.;

  ////////////////////////////////////////////////////////////////////////
  // Test on boundaries
  ////////////////////////////////////////////////////////////////////////
  if ((poly1==-1) || (poly2==-1))
    {
      if (psc >= 0)
        {
          amont = jel0;
          fluent[jel1] += psc;
        }
      else
        {
          amont = jel1;
          fluent[jel0] -= psc;
        }

      for (comp=0; comp<ncomp; comp++)
        flux(comp) = tab1(amont,comp);

    }
  else
    {
      if (psc >= 0)
        {
          amont = jel0;
          rx = rx0;
          elem1 = poly;
          elem2 = poly1;
          fluent(jel1)  += psc;
        }
      else
        {
          amont = jel1;
          rx = rx0;
          rx *= -1.;
          elem1 = poly;
          elem2 = poly2;
          fluent(jel0)  -= psc;
        }

      for (comp=0; comp<ncomp; comp++)
        {
          deltat0 = deltat1 = 0.0;
          flux(comp) = tab1(amont,comp);
          //Cerr << " flux(" << comp << ") ie phiamont= " <<  flux(comp) << finl;

          for (i=0; i<dimension; i++)
            {
              deltat0 += gradient_elem(elem1,comp,i)*rx(i);
              deltat1 += gradient_elem(elem2,comp,i)*rx(i);
            }

          if (K == 0.5)
            {
              deltat = deltat0 + deltat1;

              if  (std::fabs(deltat) <= 1.e-5)
                {
                  CF = 0.125;
                }
              else
                {
                  UTC = deltat1 / deltat;

                  if ( (UTC <= -1.) || (UTC >= 1.5) )      CF = 0.125;
                  else if ((UTC > -1.) && (UTC <= 0.))     CF = 0.5 + 0.375*UTC;
                  else if ((UTC > 0.) && (UTC <= 0.25))    CF = 0.5 - 0.625*sqrt(UTC);
                  else if ((UTC > 0.25) && (UTC < 1.5 ))   CF = 0.25* std::fabs(UTC - 1.);
                  else
                    {
                      CF=0.;
                      Process::exit();
                    }
                }
              // Compute the flux
              flux(comp) += (0.5 - CF)*deltat0 + CF*deltat1 ;
              //Cerr << " flux(" << comp << ")= " <<  flux(comp) << finl;
            }
          else
            {
              // Compute the flux
              flux(comp) += 0.25*((1.+K)*deltat0 + (1.-K)*deltat1) ;
            }
        }

    }
}

//////////////////////////////////////////////////////////////////////////
// Procedure AJOUTER
/////////////////////////////////////////////////////////////////////////

DoubleTab& Op_Conv_kschemas_VEF::ajouter(const DoubleTab& transporte,
                                         DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  const DoubleVect& porosite_face = equation().milieu().porosite_face();
  const Champ_Inc_base& la_vitesse=vitesse_.valeur();
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& face_normales = domaine_VEF.face_normales();
  const auto& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();
  const int nb_faces = domaine_VEF.nb_faces();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem = domaine_VEF.nb_elem();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleVect& volumes = domaine_VEF.volumes();
  const DoubleTab& xv = domaine_VEF.xv();
  const DoubleTab& xg = domaine_VEF.xp();
  const DoubleTab& coord = domaine.coord_sommets();
  int premiere_face_int = domaine_VEF.premiere_face_int();
  const IntTab& les_Polys = domaine.les_elems();

  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();

  int nfac = domaine.nb_faces_elem();
  int nsom = domaine.nb_som_elem();
  int nb_som_facette = domaine.type_elem()->nb_som_face();

  // For the convection treatment, standard polyhedra (not "seeing" boundary conditions)
  // are distinguished from non-standard polyhedra (having at least one boundary face).
  // A standard polyhedron has n facets on which the convection scheme is applied.
  // For a non-standard polyhedron with Dirichlet boundary conditions, part of its
  // facets are carried by the boundary faces.
  // In short, for a polyhedron the convection treatment depends on the type
  // (triangle, tetrahedron ...) and the number of Dirichlet faces.

  double psc;
  int poly,poly1,poly2,face_adj,fa7,i,j,n_bord;
  int num_face, rang ,itypcl;
  int num10,num20,num3,num_som;

  // MODIF SB su 10/09/03
  // For the following 3 elements, there are as many vertices as faces
  // making up the geometric element.
  // Problem with hexahedra: 8 vertices and 6 faces, so using the array
  // face[i] no longer works.
  // The method retained to avoid computing velocity at vertices without
  // the shape functions is therefore not usable for hexahedra,
  // where the Face=>vertices array exists but not its inverse.
  // Too costly; porosity extension to hexahedra is not done for now.

  int istetra=0;
  const Elem_VEF_base& type_elemvef= domaine_VEF.type_elem();
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF")) istetra=1;

  const int ncomp_ch_transporte= transporte.line_size();
  int fac,elem1,elem2,comp0;
  int nb_faces_ = domaine_VEF.nb_faces();
  IntVect face(nfac);

  DoubleVect flux(ncomp_ch_transporte);
  DoubleVect fluxsom(ncomp_ch_transporte);
  DoubleVect fluxg(ncomp_ch_transporte);

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

  DoubleTab tab(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              for (int comp=0; comp<ncomp_ch_transporte; comp++)
                tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }


  ///////////////////////////////////////////////////////////////////////////////////////////////
  //                        <
  // gradient computation;  < [ Ujp*np/vol(j) ]
  //                         j
  ////////////////////////////////////////////////////////////////////////////////////////////////
  DoubleTab gradient_elem(0, ncomp_ch_transporte, dimension);
  domaine_VEF.domaine().creer_tableau_elements(gradient_elem);
  // Loop over faces


  for (fac=0; fac< premiere_face_int; fac++)
    {
      elem1=face_voisins(fac,0);
      if(ncomp_ch_transporte==1)
        for (i=0; i<dimension; i++)
          {
            gradient_elem(elem1, 0, i) +=
              face_normales(fac,i)*transporte(fac);
          }
      else
        for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
          for (i=0; i<dimension; i++)
            gradient_elem(elem1, comp0, i) +=
              face_normales(fac,i)*transporte(fac,comp0);
      // dUcomp/dXi
    } // end of for faces

  for (; fac<nb_faces_; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      if(ncomp_ch_transporte==1)
        for (i=0; i<dimension; i++)
          {
            gradient_elem(elem1, 0, i) +=
              face_normales(fac,i)*transporte(fac);
            gradient_elem(elem2, 0, i) -=
              face_normales(fac,i)*transporte(fac);
          }
      else
        for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
          for (i=0; i<dimension; i++)
            {
              gradient_elem(elem1, comp0, i) +=
                face_normales(fac,i)*transporte(fac,comp0);
              gradient_elem(elem2, comp0, i) -=
                face_normales(fac,i)*transporte(fac,comp0);
            }
      // dUcomp/dXi
    } // end of for faces

  for (int elem=0; elem<nb_elem; elem++)
    for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
      for (i=0; i<dimension; i++)
        gradient_elem(elem,comp0,i) /= volumes(elem);

  gradient_elem.echange_espace_virtuel();

  ////////////////////////////////////////////////////////////////////////////////////
  // We have gradient_elem per element
  ////////////////////////////////////////////////////////////////////////////////////

  DoubleVect vs(dimension);
  DoubleVect vc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleVect cc(dimension);
  double xm;

  // Reset the array used for
  // computing the stability time step
  fluent_ = 0;

  // Non-standard polyhedra are grouped into 2 sets in Domaine_VEF:
  //  - boundary and joint polyhedra
  //  - boundary and non-joint polyhedra
  // Polyhedra are processed in the order in which they appear in the domain

  //////////////////////////////////////////////////////////////////////////////////////
  // loop over polyhedra
  //////////////////////////////////////////////////////////////////////////////////////
  const IntTab& KEL=domaine_VEF.type_elem().KEL();
  for (poly=0; poly<nb_elem; poly++)
    {
      rang = rang_elem_non_std(poly);
      if (rang==-1)
        itypcl=0;
      else
        itypcl=domaine_Cl_VEF.type_elem_Cl(rang);

      // compute the face indices of the polyhedron
      for (face_adj=0; face_adj<nfac; face_adj++)
        face(face_adj)= elem_faces(poly,face_adj);

      int scom;
      DoubleVect rx0(dimension);

      // compute velocity at the vertices of the polyhedra
      for (j=0; j<dimension; j++)
        {
          vs(j) = la_vitesse.valeurs()(face(0),j)*porosite_face(face(0));
          for (i=1; i<nfac; i++)
            vs(j)+= la_vitesse.valeurs()(face(i),j)*porosite_face(face(i));
        }

      //int ncomp;
      if (istetra==1)
        {
          for (j=0; j<nsom; j++)
            {
              for (int ncomp=0; ncomp<Objet_U::dimension; ncomp++)
                vsom(j,ncomp) =vs[ncomp] - Objet_U::dimension*la_vitesse.valeurs()(face[j],ncomp)*porosite_face(face[j]);
            }
        }
      else
        {
          // to be valid with hexahedra
          // Use the shape functions implemented in class Champs_P1_impl or Champs_Q1_impl
          //int ncomp;
          for (j=0; j<nsom; j++)
            {
              num_som = domaine.sommet_elem(poly,j);
              for (int ncomp=0; ncomp<dimension; ncomp++)
                {
                  vsom(j,ncomp) = la_vitesse.valeur_a_sommet_compo(num_som,poly,ncomp);
                }
            }
        }
      // calcul de la vitesse au centre de gravite
      domaine_VEF.type_elem().calcul_vc(face,vc,vs,vsom,vitesse(),itypcl,porosite_face);

      // Boucle sur les facettes du polyedre non standard:
      for (fa7=0; fa7<nfa7; fa7++)
        {
          //Cerr << "la facette etudiee est " << fa7 << finl;
          // fa7 separe num1 et num2. num3 est la troisieme face (2D).

          num10 = face(KEL(0,fa7));
          num20 = face(KEL(1,fa7));
          num3 = face(KEL(2,fa7));

          // Determination des elements voisins aux faces num1 et num2

          poly1 = face_voisins(num10,0);
          if (poly1==poly)
            poly1 = face_voisins(num10,1);

          poly2 = face_voisins(num20,0);
          if (poly2==poly)
            poly2 = face_voisins(num20,1);

          scom = les_Polys(poly,KEL(2,fa7));

          // compute rx0, distance between the midpoints of 'num i'

          for (i=0; i<dimension; i++)
            rx0(i) = xv(num20,i)-xv(num10,i);

          // facet normals

          if (rang==-1)
            for (i=0; i<dimension; i++)
              cc[i] = facette_normales(poly, fa7, i);
          else
            for (i=0; i<dimension; i++)
              cc[i] = normales_facettes_Cl(rang,fa7,i);

          /////////////////////////////////////////////////////////////////////////
          // Process the point for which velocity = 0.5*(vertex_velocity + midpoint_velocity)
          /////////////////////////////////////////////////////////////////////////

          for (i=0; i<nb_som_facette-1; i++)
            {
              //////////////////////////////////////////////////////////////////////////
              //Determination of PhiIJ at the midpoint between the vertex and the center of num3
              /////////////////////////////////////////////////////////////////////////

              psc = 0;
              for (j=0; j<dimension; j++)
                psc+=((vsom(KEL(i+2,fa7),j) + la_vitesse.valeurs()(num3,j) * porosite_face(num3)))*cc[j];
              psc *=0.5;
              convkschemas(K,ncomp_ch_transporte,dimension,poly,poly1,poly2,num10,num20,psc,transporte,
                           fluent_,flux,rx0,gradient_elem);

              ////////////////////////////////////////////////////////////////////////////////////////////////////////
              //Gradient limiter. Computed at the same time as the flux.
              // gradient(K0) = teta*gradient(K0)+ (1-teta)*gradient(K1 or K2)
              ////////////////////////////////////////////////////////////////////////////////////////////////////////

              double teta = 0.5;

              /////////////////////////////////////////////////////////////////////////
              // Process the vertices that are also vertices of the polyhedron
              /////////////////////////////////////////////////////////////////////////

              // XXX XXX XXX : Attention : we can not factorize more... the code is not the same
              if (ncomp_ch_transporte == 1)
                {
                  for (j=0; j<dimension; j++)
                    {
                      xm = 0.5 *(coord(scom,j)+xv(num3,j));
                      if (psc >= 0)
                        {
                          if (poly1==-1)
                            fluxsom(0) += gradient_elem(poly,0,j)*(coord(scom,j)-xm);
                          else
                            fluxsom(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly1,0,j))*(coord(scom,j)-xm);
                        }
                      else
                        {
                          if (poly2==-1)
                            fluxsom(0) += gradient_elem(poly,0,j)*(coord(scom,j)-xm);
                          else
                            fluxsom(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly2,0,j))*(coord(scom,j)-xm);
                        }
                    }
                  fluxsom(0) += flux(0);
                  fluxsom(0) *= psc;
                }
              else
                {
                  for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                    {
                      fluxsom(comp0) = flux(comp0);
                      for (j=0; j<dimension; j++)
                        {
                          xm = 0.5 *(coord(scom,j)+xv(num3,j));
                          if (psc >= 0)
                            {
                              if (poly1==-1)
                                fluxsom(comp0) += gradient_elem(poly,comp0,j)*(coord(scom,j)-xm);
                              else
                                fluxsom(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly1,comp0,j))*(coord(scom,j)-xm);
                            }
                          else
                            {
                              if (poly2==-1)
                                fluxsom(comp0) += gradient_elem(poly,comp0,j)*(coord(scom,j)-xm);
                              else
                                fluxsom(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly2,comp0,j))*(coord(scom,j)-xm);
                            }
                        }
                      fluxsom(comp0) *= psc;
                    }
                }

              ////////////////////////////////////////////////////////////////////////////
              // process the centre of gravity
              ////////////////////////////////////////////////////////////////////////////

              // XXX XXX XXX : Attention : we can not factorize more... the code is not the same
              if (ncomp_ch_transporte == 1)
                {
                  for (j=0; j<dimension; j++)
                    {
                      xm = 0.5 *(coord(scom,j)+xv(num3,j));
                      if (psc >= 0)
                        {
                          if (poly1==-1)
                            fluxg(0) += gradient_elem(poly,0,j)*(xg(poly,j)-xm);
                          else
                            fluxg(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly1,0,j))*(xg(poly,j)-xm);
                        }

                      else
                        {
                          if (poly2==-1)
                            fluxg(0) += gradient_elem(poly,0,j)*(xg(poly,j)-xm);
                          else
                            fluxg(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly2,0,j))*(xg(poly,j)-xm);
                        }
                    }

                  fluxg(0) += flux(0);
                  fluxg(0) *= psc;
                }
              else
                {
                  for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                    {
                      fluxg(comp0) = flux(comp0);
                      for (j=0; j<dimension; j++)
                        {
                          xm = 0.5 *(coord(scom,j)+xv(num3,j));
                          if (psc >= 0)
                            {
                              if (poly1==-1)
                                fluxg(comp0) += gradient_elem(poly,comp0,j)*(xg(poly,j)-xm);
                              else
                                fluxg(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly1,comp0,j))*(xg(poly,j)-xm);
                            }

                          else
                            {
                              if (poly2==-1)
                                fluxg(comp0) += gradient_elem(poly,comp0,j)*(xg(poly,j)-xm);
                              else
                                fluxg(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly2,comp0,j))*(xg(poly,j)-xm);
                            }
                        }
                      fluxg(comp0) *= psc;
                    }
                }
              //////////////////////////////////////////////////////////////////////////////
              // Integration of u.n.flux
              /////////////////////////////////////////////////////////////////////////////
              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                {
                  resu(num10,comp0) -= ( 0.5*(fluxsom(comp0)+fluxg(comp0)) );
                  resu(num20,comp0) += ( 0.5*(fluxsom(comp0)+fluxg(comp0)) );
                }
            }
        }
    } // end of the loop

  // Process joint elements of thickness 1
  for (poly=0; poly<nb_elem_tot; poly++)
    {
      // Check if a face of the polyhedron is a joint face
      for (face_adj=0; face_adj<nfac; face_adj++)
        if(face_adj<nb_faces) break;
      if(face_adj<nfac)
        {
          rang = rang_elem_non_std(poly);
          if (rang==-1)
            itypcl=0;
          else
            itypcl=domaine_Cl_VEF.type_elem_Cl(rang);

          // calcul des numeros des faces du polyedre
          for (face_adj=0; face_adj<nfac; face_adj++)
            {
              face(face_adj)= elem_faces(poly,face_adj);
              //Cerr << "les faces de l'elements sont : " << face(face_adj) << finl;
            }

          int scom;
          DoubleVect rx0(dimension);

          // calcul de la vitesse aux sommets des polyedres
          for (j=0; j<dimension; j++)
            {
              vs(j) = la_vitesse.valeurs()(face(0),j)*porosite_face(face(0));
              for (i=1; i<nfac; i++)
                vs(j)+= la_vitesse.valeurs()(face(i),j)*porosite_face(face(j));
            }
          int ncomp;
          for (j=0; j<nsom; j++)
            {
              num_som = domaine.sommet_elem(poly,j);
              for (ncomp=0; ncomp<dimension; ncomp++)
                vsom(j,ncomp) = la_vitesse.valeur_a_sommet_compo(num_som,poly,ncomp);
            }
          // calcul de la vitesse au centre de gravite

          domaine_VEF.type_elem().calcul_vc(face,vc,vs,vsom,vitesse(),itypcl,porosite_face);


          // Boucle sur les facettes du polyedre non standard:

          for (fa7=0; fa7<nfa7; fa7++)
            {
              //Cerr << "la facette etudiee est " << fa7 << finl;
              // fa7 separe num1 et num2. num3 est la troisieme face (2D).

              num10 = face(KEL(0,fa7));
              num20 = face(KEL(1,fa7));
              num3 = face(KEL(2,fa7));

              // Determination des elements voisins aux faces num1 et num2

              poly1 = face_voisins(num10,0);
              if (poly1==poly)
                {
                  poly1 = face_voisins(num10,1);
                }

              poly2 = face_voisins(num20,0);
              if (poly2==poly)
                {
                  poly2 = face_voisins(num20,1);
                }

              scom = les_Polys(poly,KEL(2,fa7));

              // compute rx0, distance between the midpoints of 'num i'

              for (i=0; i<dimension; i++)
                rx0(i) = xv(num20,i)-xv(num10,i);

              // facet normals

              if (rang==-1)
                for (i=0; i<dimension; i++)
                  cc[i] = facette_normales(poly, fa7, i);
              else
                for (i=0; i<dimension; i++)
                  cc[i] = normales_facettes_Cl(rang,fa7,i);

              /////////////////////////////////////////////////////////////////////////
              // Process the point where velocity = 0.5(vitsommet + vitmilieu)
              /////////////////////////////////////////////////////////////////////////

              for (i=0; i<nb_som_facette-1; i++)
                {
                  //////////////////////////////////////////////////////////////////////////
                  //Determine PhiIJ at the midpoint between the vertex and the midpoint of num3
                  /////////////////////////////////////////////////////////////////////////

                  psc = 0;
                  for (j=0; j<dimension; j++)
                    psc+=((vsom(KEL(i+2,fa7),j) + la_vitesse.valeurs()(num3,j) * porosite_face(num3)))*cc[j];
                  psc *=0.5;
                  convkschemas(K,ncomp_ch_transporte,dimension,poly,poly1,poly2,num10,num20,psc,transporte,
                               fluent_,flux,rx0,gradient_elem);


                  ////////////////////////////////////////////////////////////////////////////////////////////////////////
                  //Gradient limiter. Computed at the same time as the flux.
                  // gradient(K0) = teta*gradient(K0)+ (1-teta)*gradient(K1 or K2)
                  ////////////////////////////////////////////////////////////////////////////////////////////////////////

                  double teta = 0.5;

                  /////////////////////////////////////////////////////////////////////////
                  // Process the vertices that are also vertices of the polyhedron
                  /////////////////////////////////////////////////////////////////////////

                  // XXX XXX XXX : Attention : we can not factorize more... the code is not the same
                  if (ncomp_ch_transporte == 1)
                    {
                      for (j=0; j<dimension; j++)
                        {
                          xm = 0.5 *(coord(scom,j)+xv(num3,j));
                          if (psc >= 0)
                            {
                              if (poly1==-1)
                                fluxsom(0) += gradient_elem(poly,0,j)*(coord(scom,j)-xm);
                              else
                                fluxsom(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly1,0,j))*(coord(scom,j)-xm);
                            }
                          else
                            {
                              if (poly2==-1)
                                fluxsom(0) += gradient_elem(poly,0,j)*(coord(scom,j)-xm);
                              else
                                fluxsom(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly2,0,j))*(coord(scom,j)-xm);
                            }
                        }
                      fluxsom(0) += flux(0);
                      fluxsom(0) *= psc;
                    }
                  else
                    {
                      for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                        {
                          fluxsom(comp0) = flux(comp0);
                          for (j=0; j<dimension; j++)
                            {
                              xm = 0.5 *(coord(scom,j)+xv(num3,j));
                              if (psc >= 0)
                                {
                                  if (poly1==-1)
                                    fluxsom(comp0) += gradient_elem(poly,comp0,j)*(coord(scom,j)-xm);
                                  else
                                    fluxsom(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly1,comp0,j))*(coord(scom,j)-xm);
                                }
                              else
                                {
                                  if (poly2==-1)
                                    fluxsom(comp0) += gradient_elem(poly,comp0,j)*(coord(scom,j)-xm);
                                  else
                                    fluxsom(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly2,comp0,j))*(coord(scom,j)-xm);
                                }
                            }
                          fluxsom(comp0) *= psc;
                        }
                    }

                  ////////////////////////////////////////////////////////////////////////////
                  // Process the center of gravity
                  ////////////////////////////////////////////////////////////////////////////

                  // XXX XXX XXX : Attention : we can not factorize more... the code is not the same
                  if (ncomp_ch_transporte == 1)
                    {
                      for (j=0; j<dimension; j++)
                        {
                          xm = 0.5 *(coord(scom,j)+xv(num3,j));
                          if (psc >= 0)
                            {
                              if (poly1==-1)
                                fluxg(0) += gradient_elem(poly,0,j)*(xg(poly,j)-xm);
                              else
                                fluxg(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly1,0,j))*(xg(poly,j)-xm);
                            }
                          else
                            {
                              if (poly2==-1)
                                fluxg(0) += gradient_elem(poly,0,j)*(xg(poly,j)-xm);
                              else
                                fluxg(0) += (teta*gradient_elem(poly,0,j) + (1.- teta)*gradient_elem(poly2,0,j))*(xg(poly,j)-xm);
                            }
                        }
                      fluxg(0) += flux(0);
                      fluxg(0) *= psc;
                    }
                  else
                    {
                      for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                        {
                          fluxg(comp0) = flux(comp0);
                          for (j=0; j<dimension; j++)
                            {
                              xm = 0.5 *(coord(scom,j)+xv(num3,j));
                              if (psc >= 0)
                                {
                                  if (poly1==-1)
                                    fluxg(comp0) += gradient_elem(poly,comp0,j)*(xg(poly,j)-xm);
                                  else
                                    fluxg(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly1,comp0,j))*(xg(poly,j)-xm);
                                }
                              else
                                {
                                  if (poly2==-1)
                                    fluxg(comp0) += gradient_elem(poly,comp0,j)*(xg(poly,j)-xm);
                                  else
                                    fluxg(comp0) += (teta*gradient_elem(poly,comp0,j) + (1.-teta)*gradient_elem(poly2,comp0,j))*(xg(poly,j)-xm);
                                }
                            }
                          fluxg(comp0) *= psc;
                        }
                    }

                  //////////////////////////////////////////////////////////////////////////////
                  // Integration of u.n.flux
                  /////////////////////////////////////////////////////////////////////////////
                  for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                    {
                      resu(num10,comp0) -= ( 0.5*(fluxsom(comp0)+fluxg(comp0)) );
                      resu(num20,comp0) += ( 0.5*(fluxsom(comp0)+fluxg(comp0)) );
                    }
                }
            }
        }
    } // end of the loop
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Dimensioning the array of convective fluxes at the boundary
  // of the computational domain
  DoubleTab& flux_b = flux_bords_;
  flux_b.resize(domaine_VEF.nb_faces_bord(),ncomp_ch_transporte);
  flux_b = 0.;

  // Loop over the boundaries to process the boundary conditions
  // a convection term is included for Neumann_sortie_libre boundary conditions only

  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += la_vitesse.valeurs()(num_face,i)*face_normales(num_face,i)*porosite_face(num_face);
              if (psc>0)
                for (i=0; i<ncomp_ch_transporte; i++)
                  {
                    resu(num_face,i) -= psc*transporte(num_face,i);
                    flux_b(num_face,i) -= psc*transporte(num_face,i);
                  }
              else
                {
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                      flux_b(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                    }
                  fluent_(num_face) -= psc;
                }
            }
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
  modifier_flux(*this);
  return resu;
}
