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

#include <Op_Conv_Muscl_old_VEF_Face.h>
#include <Champ_P1NC.h>
#include <Neumann_sortie_libre.h>
#include <Periodique.h>

Implemente_instanciable_sans_constructeur(Op_Conv_Muscl_old_VEF_Face,"Op_Conv_Muscl_old_VEF_P1NC",Op_Conv_VEF_base);
// XD convection_muscl_old convection_deriv muscl_old NO_BRACE Only for VEF discretization.

Sortie& Op_Conv_Muscl_old_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Entree& Op_Conv_Muscl_old_VEF_Face::readOn(Entree& s )
{
  return s ;
}

//
//   Member functions of class Op_Conv_Muscl_old_VEF_Face
//



//////////////////////////////////////////////////////////////
//   MUSCL functions
////////////////////////////////////////////////////////////////

#define sgn(x) (x>0) ? 1:-1
inline double minmod(double grad1, double grad2, double gradc)
{
  int s1=sgn(grad1);
  int s2=sgn(grad2);
  int sc=sgn(gradc);
  double gradlim;
  if ((s1==s2) && (s2==sc))
    {
      gradlim=std::min(std::fabs(grad1), std::fabs(grad2));
      gradlim=std::min(std::fabs(gradlim), std::fabs(gradc));
      return sc*gradlim;
    }
  else
    return 0;
}



////////////////////////////////////////////////////////////////////
//
//                      Implementation of member functions
//
//                   of class Op_Conv_Muscl_old_VEF_Face
//
////////////////////////////////////////////////////////////////////


DoubleTab& Op_Conv_Muscl_old_VEF_Face::ajouter(const DoubleTab& transporte,
                                               DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_vef.valeur());
  const Champ_Inc_base& la_vitesse=vitesse_.valeur();
  const DoubleVect& porosite_face = equation().milieu().porosite_face();
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& face_normales = domaine_VEF.face_normales();
  const auto& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();

  const int nfa7 = domaine_VEF.type_elem().nb_facette();

  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleVect& volumes = domaine_VEF.volumes();

  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();
  int premiere_face_int = domaine_VEF.premiere_face_int();
  int nfac = domaine.nb_faces_elem();
  int nsom = domaine.nb_som_elem();
  int nb_som_facette = domaine.type_elem()->nb_som_face();
  double inverse_nb_som_facette=1./nb_som_facette;
  DoubleTab& vecteur_face_facette = ref_cast_non_const(Domaine_VEF,domaine_VEF).vecteur_face_facette();

  // For convection treatment, we distinguish standard polyhedra
  // which do not "see" boundary conditions, from non-standard
  // polyhedra which have at least one boundary face.
  // A standard polyhedron has n facets on which the convection
  // scheme is applied.
  // For a non-standard polyhedron with Dirichlet boundary conditions,
  // some facets are carried by boundary faces.
  // In short, for a polyhedron, convection treatment depends
  // on the type (triangle, tetrahedron ...) and the number of Dirichlet faces.

  const Elem_VEF_base& type_elemvef= domaine_VEF.type_elem();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang ,itypcl;
  int num10,num20,num_som;
  int ncomp_ch_transporte;
  if (transporte.nb_dim() == 1)
    ncomp_ch_transporte=1;
  else
    ncomp_ch_transporte= transporte.dimension(1);

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

  int fac=0,elem1,elem2,comp;
  int nb_faces_ = domaine_VEF.nb_faces();
  IntVect face(nfac);

  DoubleTab gradient_elem(nb_elem_tot,ncomp_ch_transporte,dimension);
  // (du/dx du/dy dv/dx dv/dy) per element
  DoubleTab gradient(0, ncomp_ch_transporte, dimension);
  domaine_VEF.creer_tableau_faces(gradient);

  // (du/dx du/dy dv/dx dv/dy) per face
  // gradient_elem=0.; already done by the constructor

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  Champ_P1NC::calcul_gradient(transporte,gradient_elem,domaine_Cl_VEF);

  // We have gradient_elem per element

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Minmod limitation
  //
  // Loop over faces
  //
  for (fac=0; fac< premiere_face_int; fac++)
    {
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (i=0; i<dimension; i++)
          gradient(fac, comp, i)
          /* upwind: */= 0;

    } // end of the face loop

  for (; fac<nb_faces_; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      double vol1=volumes(elem1);
      double vol2=volumes(elem2);
      double inverse_voltot=1./(vol1+vol2);
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (i=0; i<dimension; i++)
          {
            double grad1=gradient_elem(elem1, comp, i);
            double grad2=gradient_elem(elem2, comp, i);
            double gradc=(vol1*grad1 + vol2*grad2)*inverse_voltot;
            gradient(fac, comp, i) = minmod(grad1, grad2, gradc);
          }
    } // end of the face loop

  gradient.echange_espace_virtuel();


  DoubleVect vs(dimension);
  DoubleVect vc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleVect cc(dimension);


  const IntTab& KEL=type_elemvef.KEL();

  // Reset to zero the array used for
  // the computation of the stability time step
  fluent_ = 0;

  // Non-standard polyhedra are sorted in 2 groups in Domaine_VEF:
  //  - boundary and joint polyhedra
  //  - boundary and non-joint polyhedra
  // Polyhedra are processed in the order they appear
  // in the domain

  // loop over polyhedra
  const DoubleTab& vitesse_face=la_vitesse.valeurs();

  for (poly=0; poly<nb_elem_tot; poly++)
    {

      rang = rang_elem_non_std(poly);
      if (rang==-1)
        itypcl=0;
      else
        itypcl=domaine_Cl_VEF.type_elem_Cl(rang);

      // compute the face indices of the polyhedron
      for (face_adj=0; face_adj<nfac; face_adj++)
        face(face_adj)= elem_faces(poly,face_adj);

      for (j=0; j<dimension; j++)
        {
          vs(j) = vitesse_face(face(0),j)*porosite_face(face(0));
          for (i=1; i<nfac; i++)
            vs(j)+= vitesse_face(face(i),j)*porosite_face(face(i));
        }
      // compute velocity at the vertices of the polyhedra
      // Shape functions implemented in Champs_P1_impl or Champs_Q1_impl will be used
      if (istetra==1)
        {
          for (i=0; i<nsom; i++)
            for (j=0; j<dimension; j++)
              vsom(i,j) = vs[j] - dimension*vitesse_face(face[i],j)*porosite_face(face[i]);
        }
      else
        {
          // to be valid with hexahedra (slower to compute...)
          int ncomp;
          for (j=0; j<nsom; j++)
            {
              num_som = domaine.sommet_elem(poly,j);
              for (ncomp=0; ncomp<dimension; ncomp++)
                vsom(j,ncomp) = la_vitesse.valeur_a_sommet_compo(num_som,poly,ncomp);
            }
        }


      // compute vc (at the intersection of the 3 facets)
      type_elemvef.calcul_vc(face,vc,vs,vsom,vitesse(),itypcl,porosite_face);

      // Loop over facets of the non-standard polyhedron:
      for (fa7=0; fa7<nfa7; fa7++)
        {
          num10 = face(KEL(0,fa7));
          num20 = face(KEL(1,fa7));
          // facet normals
          if (rang==-1)
            for (i=0; i<dimension; i++)
              cc[i] = facette_normales(poly, fa7, i);
          else
            for (i=0; i<dimension; i++)
              cc[i] = normales_facettes_Cl(rang,fa7,i);

          // Apply the convection scheme to each vertex of the facet
          for (i=0; i<nb_som_facette; i++)
            {
              psc =0;
              if (i==nb_som_facette-1)
                {
                  // Treat the vertex coincident with the centroid of the polyhedron
                  for (j=0; j<dimension; j++)
                    psc += vc[j]*cc[j];
                }
              else
                {
                  // Treat the vertex or vertices that are also vertices of the polyhedron
                  for (j=0; j<dimension; j++)
                    psc+= vsom(KEL(i+2,fa7),j)*cc[j];
                }
              psc *= inverse_nb_som_facette;

              // Compute convmuscl (previously in a function that was not always inlined by the compiler)
              int comp2,amont,ii;
              double flux;
              if (psc >= 0)
                {
                  amont = num10;
                  fluent_(num20)  += psc;
                }
              else
                {
                  amont = num20;
                  fluent_(num10)  -= psc;
                }
              if (ncomp_ch_transporte == 1)
                {
                  flux = transporte(amont);
                  if (psc >= 0)
                    for (ii=0; ii<dimension; ii++)
                      flux += gradient(amont,0,ii)*vecteur_face_facette(poly,fa7,ii,0);
                  else
                    for (ii=0; ii<dimension; ii++)
                      flux += gradient(amont,0,ii)*vecteur_face_facette(poly,fa7,ii,1);

                  flux*=psc;
                  resu(num10) -= flux;
                  resu(num20) += flux;
                }
              else
                for (comp2=0; comp2<ncomp_ch_transporte; comp2++)
                  {
                    flux = transporte(amont,comp2);
                    if (psc >= 0)
                      for (ii=0; ii<dimension; ii++)
                        flux += gradient(amont,comp2,ii)*vecteur_face_facette(poly,fa7,ii,0);
                    else
                      for (ii=0; ii<dimension; ii++)
                        flux += gradient(amont,comp2,ii)*vecteur_face_facette(poly,fa7,ii,1);
                    flux *= psc;
                    resu(num10,comp2) -= flux;
                    resu(num20,comp2) += flux;
                  }

            } // end of the loop over facet vertices
        } // end of the loop over facets
    } // end of the loop




  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Loop over the boundaries to process the boundary conditions
  // a convection term is taken into account for
  // Neumann_sortie_libre boundary conditions only

  // Dimensioning the array of convective fluxes at the boundary of the computational
  // domain
  DoubleTab& flux_b = flux_bords_;
  flux_b.resize(domaine_VEF.nb_faces_bord(),ncomp_ch_transporte);
  flux_b = 0.;

  for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += la_vitesse.valeurs()(num_face,i)*face_normales(num_face,i)*porosite_face(num_face);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*transporte(num_face);
                    flux_b(num_face,0) -= psc*transporte(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*transporte(num_face,i);
                      flux_b(num_face,i) -= psc*transporte(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                      flux_b(num_face,0) -= psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
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
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
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
                    for (int comp2=0; comp2<ncomp_ch_transporte; comp2++)
                      {
                        diff1 = resu(num_face,comp2)-tab(nb_faces_perio,comp2);
                        diff2 = resu(voisine,comp2)-tab(nb_faces_perio+voisine-num_face,comp2);
                        resu(voisine,comp2)  += diff1;
                        resu(num_face,comp2) += diff2;
                        flux_b(voisine,comp2) += diff1;
                        flux_b(num_face,comp2) += diff2;
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
