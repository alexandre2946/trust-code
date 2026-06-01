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

#include <Domaine_Cl_EF.h>
#include <Domaine_EF.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>
#include <Symetrie.h>
#include <Neumann.h>
#include <Neumann_homogene.h>
#include <Periodique.h>
//#include <Champ_P0_EF.h>
#include <Champ_P1_EF.h>
#include <Domaine.h>
#include <Tri_EF.h>
#include <Tetra_EF.h>
#include <Quadri_EF.h>
#include <Equation_base.h>
#include <Champ_front_txyz.h>
#include <Champ_front_softanalytique.h>
#include <Static_Int_Lists.h>
#include <Probleme_base.h>
#include <Discretisation_base.h>
#include <Matrice_Morse.h>
#include <Dirichlet_paroi_fixe_iso_Genepi2.h>
#include <Champ_Don_base.h>

Implemente_instanciable(Domaine_Cl_EF,"Domaine_Cl_EF",Domaine_Cl_dis_base);

Sortie& Domaine_Cl_EF::printOn(Sortie& os ) const
{
  return os;
}

Entree& Domaine_Cl_EF::readOn(Entree& is )
{
  return Domaine_Cl_dis_base::readOn(is) ;
}

/*! @brief Fill the internal arrays.
 *
 */
void Domaine_Cl_EF::completer(const Domaine_dis_base& un_domaine_dis)
{
  if (sub_type(Domaine_EF,un_domaine_dis))
    {
      const Domaine_EF& le_dom_EF = ref_cast(Domaine_EF, un_domaine_dis);
      remplir_type_elem_Cl(le_dom_EF);
    }
  else
    {
      Cerr << "Domaine_Cl_EF::completer() prend comme argument un Domaine_EF " << finl;
      exit();
    }
}
void construit_connectivite_sommet(int type_cl,Static_Int_Lists& som_face_bord,const Conds_lim& les_conditions_limites_,const Domaine_EF& domaine_EF)
{
  if (type_cl!=1) Process::exit();
  int compt=0;
  int nb_faces_tot=domaine_EF.nb_faces_tot();
  ArrOfInt face_bords(nb_faces_tot);

  for(int i=0; i<les_conditions_limites_.size(); i++)
    {
      const Cond_lim_base& la_cl=les_conditions_limites_[i].valeur();
      const Front_VF& le_bord=  domaine_EF.front_VF(i);
      int num2 =  le_bord.nb_faces_tot();

      if ( (sub_type(Symetrie,la_cl)))
        {
          for (int ind_face=0; ind_face<num2; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              face_bords[compt++]=face;
            }
        }
    }
  face_bords.resize_array(compt);
  ArrOfInt is_sommet_sur_bord(domaine_EF.nb_som_tot());
  const IntTab& face_sommets=domaine_EF.face_sommets();
  int nb_som_face=face_sommets.dimension(1);
  for (int fac=0; fac<compt; fac++)
    {
      int face=face_bords[fac];
      for (int som=0; som<nb_som_face; som++)
        {
          int sommet=face_sommets(face,som);
          is_sommet_sur_bord[sommet]++;
        }
    }

  // Use is_sommet_sur_bord to set the size of the lists stored in som_face_bord.
  // - som_face_bord contains nb_som_tot lists
  // - each list is sized by the number of boundary faces (bearing a Dirichlet BC)
  //   connected to the considered vertex
  // - a given list for a vertex stores the indices of the faces connected to it

  som_face_bord.set_list_sizes(is_sommet_sur_bord);
  is_sommet_sur_bord=0;
  for (int fac=0; fac<compt; fac++)
    {
      int face=face_bords[fac];
      for (int som=0; som<nb_som_face; som++)
        {
          int sommet=face_sommets(face,som);
          int n=(is_sommet_sur_bord[sommet])++;
          som_face_bord.set_value(sommet,n,face);

        }
    }

}

static void construire_normale_locale_face(const DoubleTab& face_normales,
                                           const IntTab& faces_sommets,
                                           const DoubleTab& coord_sommets,
                                           int face,
                                           int dimension,
                                           int nb_som_face,
                                           bool is_bidim_axi,
                                           ArrOfDouble& normale_locale)
{
  for (int d = 0; d < dimension; d++)
    normale_locale[d] = face_normales(face, d);

  if (!(is_bidim_axi && norme_array(normale_locale) < 1e-12))
    return;

  for (int d = 0; d < dimension; d++)
    normale_locale[d] = 0.;

  const int som0 = faces_sommets(face, 0);
  const int som1 = faces_sommets(face, 1);
  const double dx = coord_sommets(som1, 0) - coord_sommets(som0, 0);
  const double dy = coord_sommets(som1, 1) - coord_sommets(som0, 1);
  normale_locale[0] = -dy;
  normale_locale[1] = dx;
}

/*! @brief Called by remplir_volumes_entrelaces_Cl(): fills type_elem_Cl_.
 *
 */
void Domaine_Cl_EF::remplir_type_elem_Cl(const Domaine_EF& le_dom_EF)
{
  const Domaine& z = le_dom_EF.domaine();

  const IntTab& faces_sommets=le_dom_EF.face_sommets();
  int nb_som_face=faces_sommets.dimension(1);
  const DoubleTab& coord_sommets = z.coord_sommets();
  int nb_som_tot=z.nb_som_tot();
  type_sommet_.resize_array(z.nb_som_tot());
  type_sommet_=-1;
  IntTab titi(nb_som_tot);

  for(int i=0; i<les_conditions_limites_.size(); i++)
    {
      Cond_lim_base& la_cl=les_conditions_limites_[i].valeur();
      const Front_VF& le_bord= le_dom_EF.front_VF(i);
      int num2 =  le_bord.nb_faces_tot();

      if ( (sub_type(Dirichlet,la_cl))|| (sub_type(Dirichlet_homogene,la_cl)) )
        {
          for (int ind_face=0; ind_face<num2; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              for (int s=0; s<nb_som_face; s++)
                {
                  int som=faces_sommets(face,s);
// If we have Dirichlet_paroi_fixe_iso_Genepi2, the 0 is not taken into account for the average
                  if (!sub_type(Dirichlet_paroi_fixe_iso_Genepi2,la_cl))
                    titi(som)++;
                  if ((type_sommet_[som]!=1)&& (type_sommet_[som]!=3))
                    type_sommet_[som]=2;
                  else
                    type_sommet_[som]=3;
                }
            }
        }
      else if ( (sub_type(Symetrie,la_cl)))
        {
          for (int ind_face=0; ind_face<num2; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              for (int s=0; s<nb_som_face; s++)
                {
                  int som=faces_sommets(face,s);
                  if (type_sommet_[som]<=1)
                    type_sommet_[som]=1;
                  else
                    type_sommet_[som]=3;
                }
            }
        }
      else  if ( (sub_type(Neumann,la_cl))|| (sub_type(Neumann_homogene,la_cl)) )
        {
          for (int ind_face=0; ind_face<num2; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              for (int s=0; s<nb_som_face; s++)
                {
                  int som=faces_sommets(face,s);
                  if (type_sommet_[som]<0)
                    type_sommet_[som]=0;
                }
            }
        }
      else
        {
          Cerr<<__FILE__<<":" <<(int)__LINE__<<" non code pour cette cl "<<la_cl.que_suis_je()<<finl;
        }

    }
  // if we have a Dirichlet, store 2*nb_participant + 1 if the vertex also belongs to a symmetry boundary
  for (int som=0; som<nb_som_tot; som++)
    {
      if (type_sommet_[som]>1)
        {
          //assert(titi(som)>0);
          type_sommet_[som]+=2*(titi(som));
        }
    }
  // Build the vertex -> symmetry boundary face connectivity
  if (equation().inconnue().nature_du_champ()==vectoriel)
    {
      equation().probleme().discretisation().discretiser_champ("VITESSE",le_dom_EF,"normales_nodales","1",dimension,0.,normales_symetrie_);
      equation().probleme().discretisation().discretiser_champ("CHAMP_SOMMETS",le_dom_EF,"normales_nodales_bis","1",dimension,0., normales_symetrie_bis_);
      Static_Int_Lists sommet_face_symetrie;
      int type_cl=1;
      construit_connectivite_sommet(type_cl,sommet_face_symetrie,les_conditions_limites_,le_dom_EF);
      // sommet_face_symetrie contains the number of symmetry faces associated with each vertex
      const DoubleTab& face_normales = le_dom_EF.face_normales();
      ArrOfDouble n(dimension),t1(dimension),t2(dimension),normale_locale(dimension);
      for (int som=0; som<nb_som_tot; som++)
        {
          int nbf= sommet_face_symetrie.get_list_size(som);
          if (nbf>0)
            //if ( type_sommet_(som)==1)
            {
              n=0;
              t1=0;
              t2=0;
              // determine the vertex normal
              //int nbf= sommet_face_symetrie.get_list_size(som);
              for (int f=0; f<nbf; f++)
                {
                  int face=sommet_face_symetrie(som,f);
                  construire_normale_locale_face(face_normales, faces_sommets, coord_sommets, face, dimension, nb_som_face, bidim_axi, normale_locale);
                  for (int d=0; d<dimension; d++)
                    n[d]+=normale_locale[d];
                }
              n/=nbf;

              double norm_n=norme_array(n);
              n/=norm_n;
              for (int d=0; d<dimension; d++)
                normales_symetrie_->valeurs()(som,d)=n[d];
              //	    Cerr<<som<<" must cancel a first direction "<<n(0) << " " <<n(1)<<" "<<n(dimension==3?2:1)<<finl;

              for (int f=0; f<nbf; f++)
                {
                  int face=sommet_face_symetrie(som,f);
                  construire_normale_locale_face(face_normales, faces_sommets, coord_sommets, face, dimension, nb_som_face, bidim_axi, normale_locale);
                  double prod=0;
                  for (int d=0; d<dimension; d++)
                    prod+=normale_locale[d]*n[d];

                  double s=0;
                  // double v = 0;

                  for (int d=0; d<dimension; d++)
                    {
                      t1[d]=normale_locale[d]-n[d]*prod;
                      s+=normale_locale[d]*normale_locale[d];
                      // v+=t1[d]*n[d];
                    }

                  //Cerr<<" vv"<< v<<" "<<norme_array(t1)<<" "<<norm_n<<finl;
                  if (norme_array(t1)>(1e-4*sqrt(s)))
                    {

                      // ease debugging
                      if (std::fabs(min_array(t1))>max_array(t1))
                        t1*=-1;
                      t1/=norme_array(t1);


                      //	    Cerr<<som<<" must cancel a second direction "<<t1(0) << " " <<t1(1)<<" "<<t1(dimension==3?2:1)<<finl;
                      f=nbf;
                      for (int d=0; d<dimension; d++)
                        normales_symetrie_bis_->valeurs()(som,d)=t1[d];
                      //assert(v==0);
                    }

                }
              for (int f=0; f<nbf; f++)
                {
                  int face=sommet_face_symetrie(som,f);
                  construire_normale_locale_face(face_normales, faces_sommets, coord_sommets, face, dimension, nb_som_face, bidim_axi, normale_locale);
                  double prod=0,prod1=0,s=0;
                  for (int d=0; d<dimension; d++)
                    {
                      prod+=normale_locale[d]*n[d];

                      prod1+=normale_locale[d]*t1[d];
                      s+=normale_locale[d]*normale_locale[d];

                    }
                  for (int d=0; d<dimension; d++)
                    t2[d]=normale_locale[d]-n[d]*prod-t1[d]*prod1;

                  if (norme_array(t2)>(1e-4*sqrt(s)))
                    {
                      // ease debugging
                      if (std::fabs(min_array(t2))>max_array(t2))
                        t2*=-1;
                      t2/=norme_array(t2);
                      Cerr<<face<<" "<<nbf<<" sommet "<<som<<" "<<norme_array(t2)/s<<" on doit annuler une troiseme direction"<<t2[0] << " " <<t2[1]<<" "<<t2[2]<<finl;
                      Cerr<<som<<" "<<t1[0] << " " <<t1[1]<<" "<<t1[2]<<finl;
                      Cerr<<som<<" "<<n[0] << " " <<n[1]<<" "<<n[2]<<finl;
                      f=nbf;
                      if (!normales_symetrie_ter_)
                        equation().probleme().discretisation().discretiser_champ("CHAMP_SOMMETS",le_dom_EF,"normales_nodales_bis","1",dimension,0., normales_symetrie_ter_);
                      for (int d=0; d<dimension; d++)
                        normales_symetrie_ter_->valeurs()(som,d)=t2[d];
                      //exit();
                    }
                }
            }
        }
      normales_symetrie_->valeurs().echange_espace_virtuel();
      normales_symetrie_bis_->valeurs().echange_espace_virtuel();
      //exit();
    }
}
/*! @brief Imposes symmetry conditions, i.e. cancels the field components along the normal(s).
 *
 * If tous_les_sommets_sym = 1, even vertices that also belong to a Dirichlet boundary are zeroed.
 * @param values Field values array to modify.
 * @param tous_les_sommets_sym If 1, apply symmetry to all symmetry vertices including Dirichlet ones.
 */
void Domaine_Cl_EF::imposer_symetrie(DoubleTab& values,int tous_les_sommets_sym) const
{
  // return;

  assert(values.dimension(1)==dimension);
  const Domaine& z = domaine_dis().domaine();
  int nb_som_tot=z.nb_som_tot();
  assert(values.dimension_tot(0)==nb_som_tot);

  const DoubleTab& n =normales_symetrie_->valeurs();
  const DoubleTab& n_bis =normales_symetrie_bis_->valeurs();
  int dirmax=2;
  if (normales_symetrie_ter_) dirmax=3;
  for (int som=0; som<nb_som_tot; som++)
    if (( type_sommet_[som]==1)|| ( tous_les_sommets_sym&&(type_sommet_[som]%2==1)))
      {
        for (int dir=0; dir<dirmax; dir++)
          {
            const DoubleTab& nn=(dir==0?n:(dir==1?n_bis:normales_symetrie_ter_->valeurs()));
            double prod=0;
            for (int d=0; d<dimension; d++)
              prod+=values(som,d)*nn(som,d);
            for (int d=0; d<dimension; d++)
              values(som,d)-=prod*nn(som,d);
          }
      }
}

void Domaine_Cl_EF::imposer_symetrie_partiellement(DoubleTab& values,const Noms& a_exclure) const
{
  const IntTab& faces_sommets=domaine_EF().face_sommets();
  int nb_som_face=faces_sommets.dimension(1);


  const DoubleTab& n =normales_symetrie_->valeurs();
  const DoubleTab& n_bis =normales_symetrie_bis_->valeurs();
  int dirmax=2;
  if (normales_symetrie_ter_) dirmax=3;
  int nbcond=nb_cond_lim();

  ArrOfInt type_sommet_bis(type_sommet_);
  for (int n_bord=0; n_bord<nbcond; n_bord++)
    {

      const Cond_lim_base& la_cl = les_conditions_limites(n_bord).valeur();
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num2 =  le_bord.nb_faces_tot();
      const Nom& nom_bord=le_bord.le_nom();

      if (sub_type(Symetrie,la_cl)&&(a_exclure.rang(nom_bord)>-1))
        {
          Cerr<<__FILE__<<(int)__LINE__<<" on  impose pas symetrie sur "<<nom_bord<<finl;
          for (int ind_face=0; ind_face<num2; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              for (int s=0; s<nb_som_face; s++)
                {
                  int som=faces_sommets(face,s);
                  type_sommet_bis[som]=3;
                }
            }
        }
    }

  for (int n_bord=0; n_bord<nbcond; n_bord++)
    {

      const Cond_lim_base& la_cl = les_conditions_limites(n_bord).valeur();
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num2 =  le_bord.nb_faces_tot();
      const Nom& nom_bord=le_bord.le_nom();

      if (sub_type(Symetrie,la_cl)&&(a_exclure.rang(nom_bord)<0))
        {
          Cerr<<__FILE__<<(int)__LINE__<<" on impose symetrie sur "<<nom_bord<<finl;
          for (int ind_face=0; ind_face<num2; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              for (int s=0; s<nb_som_face; s++)
                {
                  int som=faces_sommets(face,s);
                  if ( type_sommet_bis[som]==1)  //|| ( tous_les_sommets_sym&&(type_sommet_(som)%2==1)))
                    {
                      for (int dir=0; dir<dirmax; dir++)
                        {
                          const DoubleTab& nn=(dir==0?n:(dir==1?n_bis:normales_symetrie_ter_->valeurs()));
                          double prod=0;
                          for (int d=0; d<dimension; d++)
                            prod+=values(som,d)*nn(som,d);
                          for (int d=0; d<dimension; d++)
                            values(som,d)-=prod*nn(som,d);
                        }
                    }
                }
            }
        }
    }
}
void Domaine_Cl_EF::modifie_gradient(ArrOfDouble& grad_mod, const ArrOfDouble& grad, int som) const
{
  if (type_sommet_[som]!=1) return;
  assert(grad_mod.size_array()==dimension);


  const DoubleTab& n =normales_symetrie_->valeurs();
  const DoubleTab& n_bis =normales_symetrie_bis_->valeurs();

  assert ( type_sommet_[som]>=1);
  int dirmax=2;
  if (normales_symetrie_ter_) dirmax=3;
  for (int dir=0; dir<dirmax; dir++)
    {
      const DoubleTab& nn=(dir==0?n:(dir==1?n_bis:normales_symetrie_ter_->valeurs()));
      double prod=0;
      for (int d=0; d<dimension; d++)
        prod+=grad[d]*nn(som,d);
      for (int d=0; d<dimension; d++)
        grad_mod[d]+=prod*nn(som,d);
    }

}


/*! @brief Transforms la_matrice and secmem to produce a secmem normal to boundaries,
 *  plus the matrix needed to ensure that the solution is correct.
 * @param la_matrice Morse matrix to modify.
 * @param secmem Right-hand side vector to modify.
 */
void  Domaine_Cl_EF::imposer_symetrie_matrice_secmem(Matrice_Morse& la_matrice, DoubleTab& secmem) const
{
  // return;
  assert(secmem.dimension(1)==dimension);
  const Domaine& z = domaine_dis().domaine();
  int nb_som=z.nb_som();
  assert(secmem.dimension(0)==nb_som);
  int nb_comp=secmem.dimension(1);
  const DoubleTab& n =normales_symetrie_->valeurs();
  const DoubleTab& n_bis =normales_symetrie_bis_->valeurs();
  ArrOfDouble normale(dimension);

  const auto& tab1 = la_matrice.get_tab1();
  const auto& tab2 = la_matrice.get_tab2();

  const DoubleTab& champ_inconnue = equation().inconnue().valeurs();
  int dirmax=2;
  if (normales_symetrie_ter_) dirmax=3;
  for (int som=0; som<nb_som; som++)
    if ( type_sommet_[som]==1)
      {
        for (int dir=0; dir<dirmax; dir++)
          {
            const DoubleTab& nn=(dir==0?n:(dir==1?n_bis:normales_symetrie_ter_->valeurs()));
            for (int d=0; d<dimension; d++) normale[d]=nn(som,d);
            // First recompute secmem = secmem - A*present in order to be able to modify A (and project at the same time)
            auto nb_coeff_ligne=tab1[som*nb_comp+1] - tab1[som*nb_comp];
            for (int k=0; k<nb_coeff_ligne; k++)
              {
                for (int comp=0; comp<nb_comp; comp++)
                  {
                    int j=tab2[tab1[som*nb_comp+comp]-1+k]-1;
                    //assert(j!=(som*nb_comp+comp));
                    //if ((j>=(som*nb_comp))&&(j<(som*nb_comp+nb_comp)))

                    const    double coef_ij=la_matrice(som*nb_comp+comp,j);
                    int som2=j/nb_comp;
                    int comp2=j-som2*nb_comp;
                    secmem(som,comp)-=coef_ij*champ_inconnue(som2,comp2);
                  }
              }
            double somme_b=0;

            for (int comp=0; comp<nb_comp; comp++)
              somme_b+=secmem(som,comp)*normale[comp];
            //Cerr<<som<<" sommet " <<somme_b<<" "<<secmem(som,0)<<" "<<secmem(som,1)<<finl;
            // subtract secmem.n * n
            for (int comp=0; comp<nb_comp; comp++)
              secmem(som,comp)-=somme_b*normale[comp];

            if (1)
              {
                // restore the same diagonal everywhere, using the average
                double ref=0;
                for (int comp=0; comp<nb_comp; comp++)
                  {

                    int j0=som*nb_comp+comp;
                    ref+=la_matrice.coef(j0,j0);
                  }
                ref/=nb_comp;

                for (int comp=0; comp<nb_comp; comp++)
                  {
                    int j0=som*nb_comp+comp;
                    double rap=ref/la_matrice.coef(j0,j0);
                    //Cerr<<dir<<" "<<som <<" "<<comp<<" rpp "<<rap <<finl;
                    for (int k=0; k<nb_coeff_ligne; k++)
                      {

                        int j=tab2[tab1[j0]-1+k]-1;
                        la_matrice(j0,j)*=rap;
                      }
                    assert(est_egal(la_matrice(j0,j0),ref));
                  }
              }
            if (1)
              {
                // cancel all off-diagonal block coefficients
                // on rows i for which normale(d) != 0
                // i.e. fabs(normale(d)) > tol
                const double tol = 1e-12;
                for (int k=0; k<nb_coeff_ligne; k++)
                  {

                    for (int comp=0; comp<nb_comp; comp++)
                      if (std::fabs(normale[comp])>tol)
                        {
                          int j=tab2[tab1[som*nb_comp+comp]-1+k]-1;
                          if (j!=(som*nb_comp+comp))
                            if ((j>=(som*nb_comp))&&(j<(som*nb_comp+nb_comp)))
                              {
                                la_matrice(som*nb_comp+comp,j)=0;
                              }
                        }
                  }
              }
            {
              // for off-diagonal blocks, ensure that Aij.ni = 0

              ArrOfDouble somme((int)nb_coeff_ligne);
              for (int k=0; k<nb_coeff_ligne; k++)
                {

                  int j=tab2[tab1[som*nb_comp]-1+k]-1;
                  for (int comp=0; comp<nb_comp; comp++)
                    somme[k]+=la_matrice(som*nb_comp+comp,j)*normale[comp];
                }
              // subtract somme * ni
              for (int k=0; k<nb_coeff_ligne; k++)
                {

                  int j=tab2[tab1[som*nb_comp]-1+k]-1;
                  for (int comp=0; comp<nb_comp; comp++)
                    if ((j<(som*nb_comp))||(j>=(som*nb_comp+nb_comp)))
                      la_matrice(som*nb_comp+comp,j)-=(somme[k])*normale[comp];
                }
            }
            // Finally recompute secmem = secmem + A*champ_inconnue (A has been heavily modified)
            for (int k=0; k<nb_coeff_ligne; k++)
              {
                for (int comp=0; comp<nb_comp; comp++)
                  {
                    int j=tab2[tab1[som*nb_comp+comp]-1+k]-1;
                    int som2=j/nb_comp;
                    int comp2=j-som2*nb_comp;

                    const    double coef_ij=la_matrice(som*nb_comp+comp,j);
                    secmem(som,comp)+=coef_ij*champ_inconnue(som2,comp2);

                  }
              }
            {
              // verification
              double somme_b2=0;

              for (int comp=0; comp<nb_comp; comp++)
                somme_b2+=secmem(som,comp)*normale[comp];
              //Cerr<<" lllllllll "<<somme_b2<<" "<<tt<<finl;
              if (std::fabs(somme_b2) >= 1e-8)
                Cerr << "Domaine_Cl_EF::imposer_symetrie_matrice_secmem: secmem.n != 0 ("
                     << somme_b2 << ") au sommet " << som << ", projection appliquee." << finl;
              // subtract secmem.n * n
              for (int comp=0; comp<nb_comp; comp++)
                secmem(som,comp)-=somme_b2*normale[comp];

            }
          }
      }
  //  exit();
}

/*! @brief Imposes boundary conditions at time "temps" of the Champ_Inc.
 * @param ch Instationary field to apply boundary conditions to.
 * @param temps Current time value.
 */
void Domaine_Cl_EF::imposer_cond_lim(Champ_Inc_base& ch, double temps)
{
  DoubleTab& ch_tab = ch.valeurs(temps);
  int nb_comp = ch.nb_comp();
  const Domaine_EF& domaineEF =  domaine_EF(); //ref_cast(Domaine_EF,ch.equation().domaine_dis());
  const IntTab& faces_sommets=domaineEF.face_sommets();
  int nb_som_face=faces_sommets.dimension(1);
  const DoubleTab& coords= domaineEF.domaine().coord_sommets();

  // first zero the field on Dirichlet boundaries
  // then add 1/nb_cl * val_imp

  for (int n_bord=0; n_bord<nb_cond_lim(); n_bord++)
    {

      const Cond_lim_base& la_cl = les_conditions_limites(n_bord).valeur();
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num2 =  le_bord.nb_faces_tot();

      if (sub_type(Dirichlet,la_cl)||(sub_type(Dirichlet_homogene,la_cl)))
        {
          //const Dirichlet& la_cl_diri = ref_cast(Dirichlet,la_cl);

          for (int ind_face=0; ind_face<num2; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              for (int s=0; s<nb_som_face; s++)
                {
                  int som=faces_sommets(face,s);
                  assert(type_sommet_[som]>=2);


                  if (nb_comp == 1)
                    ch_tab[som]=0;
                  else
                    for (int ncomp=0; ncomp<nb_comp; ncomp++)
                      ch_tab(som,ncomp)=0;
                }
            }
        }
    }
  for (int n_bord=0; n_bord<nb_cond_lim(); n_bord++)
    {

      const Cond_lim_base& la_cl = les_conditions_limites(n_bord).valeur();
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num2 =  le_bord.nb_faces_tot();

      if (sub_type(Dirichlet,la_cl))
        {
          const Dirichlet& la_cl_diri = ref_cast(Dirichlet,la_cl);
          if (sub_type(Champ_front_softanalytique,la_cl_diri.champ_front()))
            {
              Cerr<<" You must use Champ_front_fonc_txyz instead of "<<la_cl_diri.champ_front().que_suis_je()<<finl;
              exit();
            }
          int avec_valeur_aux_sommets=0;
          if (sub_type(Champ_front_var_instationnaire,la_cl_diri.champ_front()))
            {
              const Champ_front_var_instationnaire& ch_txyz=ref_cast(Champ_front_var_instationnaire,la_cl_diri.champ_front());
              avec_valeur_aux_sommets=ch_txyz.valeur_au_temps_et_au_point_disponible();
            }
          // For Dirichlet faces, impose the unknown at the vertex
          if (avec_valeur_aux_sommets)
            {
              const Champ_front_var_instationnaire& ch_txyz=ref_cast(Champ_front_var_instationnaire,la_cl_diri.champ_front());
              for (int ind_face=0; ind_face<num2; ind_face++)
                {
                  int face=le_bord.num_face(ind_face);
                  for (int s=0; s<nb_som_face; s++)
                    {
                      int som=faces_sommets(face,s);
                      assert(type_sommet_[som]>=4);
                      double coef=1./(type_sommet_[som]/2-1);
                      //Cerr<<"iciPB "<<coef<<finl;
                      double x,y,z=0;
                      x=coords(som,0);
                      y=coords(som,1);
                      if (dimension==3)
                        z=coords(som,2);
                      if (nb_comp == 1)
                        ch_tab[som]+=coef*ch_txyz.valeur_au_temps_et_au_point(temps,som,x,y,z,0);
                      else
                        for (int ncomp=0; ncomp<nb_comp; ncomp++)
                          ch_tab(som,ncomp)+=coef*ch_txyz.valeur_au_temps_et_au_point(temps,som,x,y,z,ncomp);
                    }
                }
            }
          else
            for (int ind_face=0; ind_face<num2; ind_face++)
              {
                int face=le_bord.num_face(ind_face);
                for (int s=0; s<nb_som_face; s++)
                  {
                    int som=faces_sommets(face,s);
                    assert(type_sommet_[som]>=4);
                    double coef=1./(type_sommet_[som]/2-1);
                    //Cerr<<"iciPB "<<coef<<finl;
                    if (nb_comp == 1)
                      ch_tab[som]+=coef*la_cl_diri.val_imp_au_temps(temps,ind_face);
                    else
                      for (int ncomp=0; ncomp<nb_comp; ncomp++)
                        ch_tab(som,ncomp)+=coef*la_cl_diri.val_imp_au_temps(temps,ind_face,ncomp);
                  }
              }

        }
      /*
      else if (sub_type(Dirichlet_homogene,la_cl))
      {
        // For Dirichlet faces, the unknown is imposed at the vertex
        for (int ind_face=0; ind_face<num2; ind_face++)
      {
      int face=le_bord.num_face(ind_face);
      for (int s=0;s<nb_som_face;s++)
        {
          int som=faces_sommets(face,s);
          assert(type_sommet_(som)>=2);
          if (nb_comp == 1)
      ch_tab[som] = 0;
          else
      for (int ncomp=0; ncomp<nb_comp; ncomp++)
        ch_tab(som,ncomp) =0;
        }
      }
      }
      */
      // provisional: to be done only once
      else if ( (sub_type(Symetrie,la_cl) ) &&
                (ch.nature_du_champ()==vectoriel) )
        {
          imposer_symetrie(ch_tab);
        }
    }



}

int Domaine_Cl_EF::nb_faces_sortie_libre() const
{
  exit();
  /*
  int compteur=0;
  for(int cl=0; cl<les_conditions_limites_.size(); cl++)
    {
      if(sub_type(Neumann_sortie_libre, les_conditions_limites_[cl].valeur()))
  {
    const Front_VF& le_bord=ref_cast(Front_VF,les_conditions_limites_[cl]->frontiere_dis());
    compteur+=le_bord.nb_faces();
  }
    }
  return compteur;
  */
  return -1;
}

int Domaine_Cl_EF::nb_bord_periodicite() const
{
  int compteur = 0;
  for (const auto &itr : les_conditions_limites_)
    {
      if (sub_type(Periodique, itr.valeur()))
        compteur++;
    }
  return compteur;
}


int Domaine_Cl_EF::initialiser(double temps)
{
  Domaine_Cl_dis_base::initialiser(temps);

  if (nb_bord_periodicite()>0)
    {
      Cerr<<" Periodicity is not implemented !!!"<<finl;
      abort();
    }
  return 1;
}

Domaine_EF& Domaine_Cl_EF::domaine_EF()
{
  return ref_cast(Domaine_EF, domaine_dis());
}

const Domaine_EF& Domaine_Cl_EF::domaine_EF() const
{
  return ref_cast(Domaine_EF, domaine_dis());
}
