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

#include <Op_Conv_Centre_old_VEF_Face.h>
#include <Periodique.h>
#include <Neumann_sortie_libre.h>

Implemente_instanciable(Op_Conv_Centre_old_VEF_Face,"Op_Conv_Centre_old_VEF_P1NC",Op_Conv_VEF_base);
// XD convection_centre_old convection_deriv centre_old NO_BRACE Only for VEF discretization.

//// printOn
//

Sortie& Op_Conv_Centre_old_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Op_Conv_Centre_old_VEF_Face::readOn(Entree& s )
{
  return s ;
}

//
//   Functions of class Op_Conv_Centre_old_VEF_Face
//
void Op_Conv_Centre_old_VEF_Face::associer(const Domaine_dis_base& domaine_dis,
                                           const Domaine_Cl_dis_base& domaine_cl_dis,
                                           const Champ_Inc_base& ch_transporte)
{
  const Domaine_VEF& zvef = ref_cast(Domaine_VEF,domaine_dis);
  const Domaine_Cl_VEF& zclvef = ref_cast(Domaine_Cl_VEF,domaine_cl_dis);
  const Champ_Inc_base& le_ch_transporte = ref_cast(Champ_Inc_base,ch_transporte);

  le_dom_vef = zvef;
  la_zcl_vef = zclvef;
  champ_transporte = le_ch_transporte;

  fluent_.reset();
  le_dom_vef->creer_tableau_faces(fluent_);
}

DoubleTab& Op_Conv_Centre_old_VEF_Face::ajouter(const DoubleTab& transporte,
                                                DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  //  const Champ_Inc_base& le_transporte = champ_transporte.valeur();
  const Champ_Inc_base& la_vitesse =vitesse_.valeur();

  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& face_normales = domaine_VEF.face_normales();
  const auto& facette_normales = domaine_VEF.facette_normales();
  //  const DoubleVect& volumes_entrelaces = domaine_VEF.volumes_entrelaces();
  const Domaine& domaine = domaine_VEF.domaine();
  //  const int nb_faces = domaine_VEF.nb_faces();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  //  const int nb_elem = domaine_VEF.nb_elem();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();


  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();
  //  const DoubleVect& volumes_entrelaces_Cl = domaine_Cl_VEF.volumes_entrelaces_Cl();

  const DoubleVect& porosite_face = equation().milieu().porosite_face();

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
  double flux;
  int poly,face_adj,fa7,i,j,comp0,n_bord;
  int num_face, rang ,itypcl;
  int num10, num20, num_som;

  int ncomp_ch_transporte;
  if (transporte.nb_dim() == 1)
    ncomp_ch_transporte=1;
  else
    ncomp_ch_transporte= transporte.dimension(1);

  // MODIF SB on 10/09/03
  // For the following 3 elements, there are as many vertices as faces
  // making up the geometric element.
  // Problem with hexahedra: 8 vertices and 6 faces, so use of the array
  // face[i] no longer works.
  // The chosen method to avoid computing the velocity at vertices without
  // shape functions is therefore not applicable;
  // for the hexa we have no access to the face.
  // The Face=>vertices array exists but not the inverse.
  // Too costly and for now porosities are not extended to hexahedra.

  int istetra=0;
  const Elem_VEF_base& type_elemvef= domaine_VEF.type_elem();
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  IntVect face(nfac);
  DoubleVect vs(dimension);
  DoubleVect vc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleVect cc(dimension);

  // declaration for the transported field

  DoubleVect ts(ncomp_ch_transporte);
  DoubleVect tc(ncomp_ch_transporte);
  DoubleTab tsom(nsom,ncomp_ch_transporte);


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

      for (j=0; j<dimension; j++)
        {
          vs[j] = la_vitesse.valeurs()(face[0],j)*porosite_face(face[0]);
          for (i=1; i<nfac; i++)
            vs[j]+= la_vitesse.valeurs()(face[i],j)*porosite_face(face[i]);
        }
      // int ncomp;
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

      // compute vc
      domaine_VEF.type_elem().calcul_vc(face,vc,vs,vsom,vitesse(),
                                        itypcl,porosite_face);

      // compute the transported field at the polyhedron vertices, tsom
      if(ncomp_ch_transporte == 1)
        {
          ts[0]=transporte(face[0]);
          for (i=1; i<nfac; i++)
            ts[0]+= transporte(face[i]);

          for (i=0; i<nsom; i++)
            tsom(i,0) = ts[0] - dimension*transporte(face[i],0);
        }
      else
        {
          for (j=0; j<ncomp_ch_transporte; j++)
            {
              ts[j] = transporte(face[0],j);
              for (i=1; i<nfac; i++)
                ts[j]+= transporte(face[i],j);
            }
          for (i=0; i<nsom; i++)
            for (j=0; j<ncomp_ch_transporte; j++)
              tsom(i,j) = ts[j] - dimension*transporte(face[i],j);
        }

      // compute the transported field at the centre of gravity, tc

      for (j=0; j<ncomp_ch_transporte; j++)
        tc[j] = ts[j]/nfac;


      // Loop over the facets of the non-standard polyhedron:

      for (fa7=0; fa7<nfa7; fa7++)
        {
          num10 = face[KEL(0,fa7)];
          num20 = face[KEL(1,fa7)];
          if (rang==-1)
            for (i=0; i<dimension; i++)
              cc[i] = facette_normales(poly,fa7,i);
          else
            for (i=0; i<dimension; i++)
              cc[i] = normales_facettes_Cl(rang,fa7,i);

          // Apply the convection scheme at each vertex of the facet

          // Treat the vertex/vertices that are also vertices of the polyhedron

          int isom;
          for (i=0; i<nb_som_facette-1; i++)
            {
              isom = KEL(i+2,fa7);
              psc =0;
              for (j=0; j<dimension; j++)
                psc+= vsom(isom,j)*cc[j];
              psc /= nb_som_facette;

              if(psc >= 0)
                fluent_[num20] += psc;
              else
                fluent_[num10] -= psc;

              // write the flux
              if (ncomp_ch_transporte == 1)
                {
                  flux = tsom(isom,0)*psc;
                  resu(num10) -= flux;
                  resu(num20) += flux;
                }
              else
                for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                  {
                    flux = tsom(isom,comp0)*psc;
                    resu(num10, comp0) -= flux;
                    resu(num20, comp0) += flux;
                  }
            }


          // Treat the vertex coinciding with the centre of gravity of the polyhedron

          psc=0;
          for (j=0; j<dimension; j++)
            psc += vc[j]*cc[j];
          psc /= nb_som_facette;

          if(psc >= 0)
            fluent_[num20] += psc;
          else
            fluent_[num10] -= psc;

          // write the flux

          if (ncomp_ch_transporte == 1)
            {
              flux = tc[0]*psc;
              resu(num10) -= flux;
              resu(num20) += flux;
            }
          else
            for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
              {
                flux = tc[comp0]*psc;
                resu(num10, comp0) -= flux;
                resu(num20, comp0) += flux;
              }
        }

    } // end of loop

  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

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
                  fluent_[num_face] -= psc;
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
  modifier_flux(*this);
  return resu;

}
