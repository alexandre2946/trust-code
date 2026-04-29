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

#include <Source_PDF_VEF.h>
#include <Domaine.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Schema_Temps_base.h>
#include <Discretisation_base.h>
#include <Matrice_Morse.h>
#include <Matrice_Bloc.h>
#include <TRUSTTrav.h>
#include <Interpolation_IBM_elem_fluid.h>
#include <Interpolation_IBM_mean_gradient.h>
#include <Interpolation_IBM_hybrid.h>
#include <Interpolation_IBM_power_law_tbl.h>
#include <Interpolation_IBM_power_law_tbl_u_star.h>
#include <Interpolation_IBM_thermal_wall_law.h>
#include <Dirichlet.h>
#include <SFichier.h>
#include <Champ_P1NC.h>
#include <Op_Diff_VEF_base.h>

Implemente_instanciable(Source_PDF_VEF,"Source_PDF_VEF_P1NC",Source_PDF_base);

/*##################################################################################################
####################################################################################################
################################# READON AND PRINTON ###############################################
####################################################################################################
##################################################################################################*/



Entree& Source_PDF_VEF::readOn(Entree& s)
{
  Source_PDF_base::readOn(s);
  return s;
}

Sortie& Source_PDF_VEF::printOn(Sortie& s ) const
{
  Source_PDF_base::printOn(s);
  return s;
}

/*##################################################################################################
####################################################################################################
################################# PROBLEM ASSOCIATION AND ZONES ####################################
####################################################################################################
##################################################################################################*/

void Source_PDF_VEF::associer_pb(const Probleme_base& pb)
{
  Source_PDF_base::associer_pb(pb);

  int nb_som=le_dom_VEF->nb_som();
  tab_u_star_ibm_.resize(nb_som);
  tab_y_plus_ibm_.resize(nb_som);
}

void Source_PDF_VEF::compute_indicateur_nodal_champ_aire()
{
  Source_PDF_base::compute_indicateur_nodal_champ_aire();
}

void Source_PDF_VEF::associer_domaines(const Domaine_dis_base& domaine_dis,
                                       const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_VEF = ref_cast(Domaine_VEF, domaine_dis);
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis);
}

void Source_PDF_VEF::multiply_coeff_volume(DoubleTab& coeff) const
{
  const DoubleVect& vol_diam=ref_cast(Domaine_VEF, le_dom_VEF.valeur()).volumes_entrelaces();
  int n = vol_diam.size_totale();
  int nc = coeff.dimension_tot(0) ;
  if (n != nc)
    {
      Cerr<<" dimensions differentes n nc = "<<n<<" "<<nc<<finl;
      exit();
    }
  int ndim = coeff.dimension(1) ;
  for (int i=0; i<n; i++)
    {
      for (int comp=0; comp<ndim; comp++) coeff(i,comp) = coeff(i,comp)*vol_diam(i);
    }
}

/*##################################################################################################
####################################################################################################
################################# AJOUTER_ #########################################################
####################################################################################################
##################################################################################################*/

/*
This function redirects toward the ajouter_ which correspond to the model chosen.
*/

DoubleTab& Source_PDF_VEF::ajouter_(const DoubleTab& variable, DoubleTab& resu, const int i_traitement_special) const
{
  // calcul de : coefficient rho/dt * coeff_relax * (volume_thilde / 8) * variable
  /* i_traitement_special = 0 => traitement classique; coefficient (1/eta) */
  /* i_traitement_special = 1 => coefficient (1/eta -> 1) */
  /* i_traitement_special = 2 => coefficient (1/eta -> 1 + 1/eta) */

  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& elems= domaine_VEF.elem_faces() ;
  int nb_face_elem=domaine_VEF.domaine().nb_faces_elem();
  int nb_elems=domaine_VEF.domaine().nb_elem_tot();
  const DoubleVect& volume=domaine_VEF.volumes_entrelaces();
  int dim_esp = Objet_U::dimension ;
  int dim_var = variable.dimension(1);
  // Cerr << "erreur ici" << finl;
  if (dim_var != 1 && dim_var != dim_esp)
    {
      Cerr<<"ajouter_: dimension variable differente 1 ou "<<dim_esp<<"! "<<finl;
      exit();
    }
  ArrOfDouble tuvw(dim_var);
  const DoubleTab& rotation=champ_rotation_->valeurs();
  const DoubleTab& aire=champ_aire_->valeurs();
  //OWN_PTR(Champ_Don_base) rho_test;
  //champ_rho_->affecter(equation().probleme().get_champ("masse_volumique"));
  const DoubleTab& rho_m=champ_rho_->valeurs();
  DoubleTab pond;

  // return = rho/dt * coeff_relax * volume
  int coef1 = nb_face_elem;
  pond = compute_pond(rho_m, aire, volume, coef1, nb_elems);

  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      if (aire(num_elem)>0.)
        {
          if (i_traitement_special == 0)
            {
              for (int k=0; k<dim_var; k++) tuvw[k] =  1.0 / modele_lu_.eta_;
            }
          else if (i_traitement_special == 1) //terme temps en rho v
            {
              for (int k=0; k<dim_var; k++) tuvw[k] =  1.0;
            }
          else if (i_traitement_special == 101) //terme temps en v
            {
              for (int k=0; k<dim_var; k++) tuvw[k] =  1.0 / rho_m(num_elem);
            }
          else if (i_traitement_special == 2)
            {
              for (int k=0; k<dim_var; k++) tuvw[k] =  1.0 + 1.0 / modele_lu_.eta_;
            }
          else if (i_traitement_special == 102)
            {
              for (int k=0; k<dim_var; k++) tuvw[k] =  1.0 / rho_m(num_elem) + 1.0 / modele_lu_.eta_;
            }
          else
            {
              Cerr << "Source_PDF_VEF::ajouter_ : i_traitement_special should be 0, 1, 2, 101 or 102" << finl;
              exit();
            }
          for (int comp=0; comp<dim_var; comp++)
            {
              double tijvj=0;
              if (dim_var == 1)
                {
                  tijvj = tuvw[0];
                }
              else
                // dim_var = dim_esp
                {
                  for (int c=0; c<dim_var; c++)
                    {
                      double coeff_im=0;
                      for(int k=0; k<dim_esp; k++)
                        {
                          coeff_im+=rotation(num_elem,3*k+comp)*tuvw[k]*rotation(num_elem,3*k+c);
                        }
                      tijvj+=coeff_im ;
                    }
                }
              // Cerr<< "Source_PDF_VEF::ajouter_ tijvj=" <<tijvj<<" pond("<<num_elem<<") = "<<pond(num_elem)<< finl;
              // cas generique : 1/eta rho/dt * coeff_relax * volume * variable
              for (int i=0; i<nb_face_elem; i++)
                {
                  resu(elems(num_elem,i),comp)-=pond(num_elem)*coef1*tijvj*variable(elems(num_elem,i),comp);
                  // Cerr<< "Source_PDF_VEF::ajouter_ noeud: "<<i <<" "<<pond(num_elem)*coef1*tijvj*variable(elems(num_elem,i),comp)<<" variable = "<<variable(elems(num_elem,i),comp)<< finl;
                }
            }
        }
    }
  return resu;
}

DoubleTab& Source_PDF_VEF::ajouter_(const DoubleTab& variable, DoubleTab& resu) const
{
  const int i_traitement_special = 0 ;
  ajouter_(variable, resu, i_traitement_special) ;
  return resu;
}

/*##################################################################################################
####################################################################################################
################################# CONTRIBUER_A_AVEC ################################################
####################################################################################################
##################################################################################################*/

/*
This function redirects toward the contribuer_avec_ which correspond to the model chosen.
*/

void  Source_PDF_VEF::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  // calcul : 1/eta rho/dt * coeff_relax * volume
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& elems= domaine_VEF.elem_faces();
  int nb_face_elem=domaine_VEF.domaine().nb_faces_elem();
  int nb_elems=domaine_VEF.domaine().nb_elem_tot();
  const DoubleVect& volume=domaine_VEF.volumes_entrelaces();
  const DoubleTab& variable=equation().inconnue().valeurs();
  int dim_var = variable.dimension(1);
  int dim_esp = Objet_U::dimension ;
  if (dim_var != 1 && dim_var != dim_esp)
    {
      Cerr<<"contribuer_a_avec: dimension variable differente 1 ou "<<dim_esp<<"! "<<finl;
      exit();
    }
  ArrOfDouble tuvw(dim_var);
  const DoubleTab& rotation=champ_rotation_->valeurs();
  const DoubleTab& aire=champ_aire_->valeurs();
  //champ_rho_->affecter(equation().probleme().get_champ("masse_volumique"));
  const DoubleTab& rho_m=champ_rho_->valeurs();
  DoubleTab pond ;

  // return = rho/dt * coeff_relax * volume pour elements interceptes
  int coef1 = nb_face_elem;
  pond = compute_pond(rho_m, aire, volume, coef1, nb_elems);

  for (int num_elem=0; num_elem<nb_elems; num_elem++) //elements loop
    {
      if (aire(num_elem)>0.)
        {
          for (int k=0; k<dim_var; k++) tuvw[k] = 1.0 / modele_lu_.eta_;
          // Cerr<< "Source_PDF_VEF::contribuer_a_avec: tuvw=" <<tuvw<<" pond("<<num_elem<<") = "<<pond(num_elem)<< finl;

          for (int comp=0; comp<dim_var; comp++)
            {
              for (int comp2=0; comp2<dim_var; comp2++)
                {
                  for (int i=0; i<nb_face_elem; i++)
                    {
                      for (int j=0; j<nb_face_elem; j++)
                        {
                          int s1=elems(num_elem,i);

                          double coeff_im=0;
                          if (dim_var == 1)
                            {
                              coeff_im = tuvw[0];
                            }
                          else
                            // dim_var = dim_esp
                            {
                              for (int k=0; k<dim_esp; k++)
                                {
                                  coeff_im+=rotation(num_elem,3*k+comp)*tuvw[k]*rotation(num_elem,3*k+comp2);
                                }
                            }
                          int c1=s1*dim_var+comp;
                          // Cerr << c1 << finl;

                          // 1.0/eta * rho/dt * coeff_relax * volume
                          matrice.coef(c1,c1)+=pond(num_elem)*coeff_im;
                          // Cerr<< "Source_PDF_VEF::contribuer_a_avec: contribution matrice c1c1 : "<<c1<<" = "<<pond(num_elem)*coeff_im<< finl;
                        }
                    }
                }
            }
        }
    }
  // verif_ajouter_contrib(variable, matrice) ;
}

void  Source_PDF_VEF::verif_ajouter_contrib(const DoubleTab& variable, Matrice_Morse& matrice) const
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& elems= domaine_VEF.elem_faces();
  int nb_face_elem=domaine_VEF.domaine().nb_faces_elem();
  int nb_elems=domaine_VEF.domaine().nb_elem_tot();
  const DoubleTab& aire=champ_aire_->valeurs();

  DoubleTrav force(variable);
  int dim_var = variable.dimension(1);
  ajouter_(variable,force);
  DoubleTrav force2(variable);
  matrice.ajouter_multvect_(variable, force2) ;
  DoubleTab diforce(force);
  diforce += force2 ;
  double difmax = 0.;
  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      if (aire(num_elem)>0.)
        {
          for (int i=0; i<nb_face_elem; i++)
            {
              int s1=elems(num_elem,i);
              double difmod2 = 0.;
              for (int k=0; k<dim_var; k++) difmod2 += diforce(s1,k)*diforce(s1,k);
              if (difmod2 > 1.0e-5)
                {
                  if (difmod2 > difmax) difmax = difmod2 ;
                  // for (int k=0; k<dim_var; k++) Cerr << " ajouter multvect diforce x = "<< force(s1,k) << " "<< force2(s1,k) << " " << diforce(s1,k)<<endl;
                }
            }
        }
    }
  if (difmax > 0.)
    {
      Cerr<< "Source_PDF_VEF: Max norme caree diff. force absolue = "<<difmax<<endl;
    }
}

/*##################################################################################################
####################################################################################################
################################# calculer_variable_imposee ########################################
####################################################################################################
##################################################################################################*/

void Source_PDF_VEF::calculer_variable_imposee_elem_fluid()
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  int nb_faces=domaine_VEF.nb_faces();
  int dim_esp = Objet_U::dimension;
  const Domaine& dom = domaine_VEF.domaine();
  const IntTab& faces_sommets=domaine_VEF.face_sommets();
  int nb_som_face=faces_sommets.dimension(1);
  int nb_som = domaine_VEF.nb_som();
  DoubleTab& variable_imposee_mod = modele_lu_.variable_imposee_->valeurs();
  DoubleTab& variable_imposee_calculee = variable_imposee_;
  Interpolation_IBM_elem_fluid& interp = ref_cast(Interpolation_IBM_elem_fluid,interpolation_lue_.valeur());
  DoubleTab& fluid_points = interp.fluid_points_->valeurs();
  DoubleTab& solid_points = interp.solid_points_->valeurs();
  DoubleTab& fluid_elems = interp.fluid_elems_->valeurs();
  Champ_P1NC& champ_variable_inconnue = ref_cast(Champ_P1NC,equation().inconnue());
  DoubleTab& val_variable_inconnue = champ_variable_inconnue.valeurs();
  int nb_comp = val_variable_inconnue.dimension(1);
  DoubleTab variable_imposee_sommet(nb_som,nb_comp);
  IntTab sommet_interp(1, nb_som);
  DoubleTab xf(1, dim_esp);
  DoubleTab vf(1, nb_comp);
  ArrOfInt cells(1);
  double eps = 1e-12;
  variable_imposee_calculee = 0.0;
  variable_imposee_sommet = 0.0;
  sommet_interp = 0;
  for (int i = 0; i < nb_som; i++)
    {
      if ((fluid_elems(i) >= 0.0))
        {
          sommet_interp(0, i) = 1;
          double d1 = 0.0;
          double d2 = 0.0;
          for(int j = 0; j < dim_esp; j++)
            {
              double xj = dom.coord(i,j);
              double xjf = fluid_points(i,j);
              xf(0, j) = xjf;
              double xjs = solid_points(i,j);
              d1 += (xj-xjs)*(xj-xjs);
              d2 += (xjf-xj)*(xjf-xj);
            }
          d1 = sqrt(d1);
          d2 = sqrt(d2);
          double inv_d ;
          if ( (d1 + d2) > eps ) inv_d = 1.0 / (d1 + d2);
          else inv_d = 0. ;
          cells[0] = int(fluid_elems(i));
          champ_variable_inconnue.valeur_a_elem(xf, vf, cells[0]);

          for (int j = 0; j < nb_comp; j++)
            {
              double vjs = variable_imposee_mod(i,j);
              variable_imposee_sommet(i,j) = vjs + (vf(0, j)-vjs)*inv_d*d1;
            }
        }
      else
        {
          for(int j = 0; j < nb_comp; j++)
            {
              variable_imposee_sommet(i,j) = variable_imposee_mod(i,j);
            }
        }
    }
  for (int f = 0; f < nb_faces; f++)
    {
      int compteur = 0;
      int som_f = 0;
      while ((som_f < nb_som_face) and (sommet_interp(0, faces_sommets(f,som_f))==1))
        {
          compteur++;
          som_f++;
        }
      if (compteur == nb_som_face)
        {
          for (int som_=0; som_ < nb_som_face; som_++)
            {
              int i = faces_sommets(f,som_);
              for(int j = 0; j < nb_comp; j++)
                {
                  variable_imposee_calculee(f,j) += variable_imposee_sommet(i,j)/nb_som_face;
                }
            }
        }
    }
}

void Source_PDF_VEF::calculer_variable_imposee_mean_grad()
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  int nb_faces=domaine_VEF.nb_faces();
  int dim_esp = Objet_U::dimension;
  const Domaine& dom = domaine_VEF.domaine();
  const IntTab& faces_sommets=domaine_VEF.face_sommets();
  int nb_som_face=faces_sommets.dimension(1);
  int nb_som = domaine_VEF.nb_som();
  DoubleTab& variable_imposee_mod = modele_lu_.variable_imposee_->valeurs();
  DoubleTab& variable_imposee_calculee = variable_imposee_;
  Interpolation_IBM_mean_gradient& interp = ref_cast(Interpolation_IBM_mean_gradient,interpolation_lue_.valeur());
  DoubleTab& solid_points = interp.solid_points_->valeurs();
  DoubleTab& is_dirichlet = interp.is_dirichlet_->valeurs();
  double eps = 1e-12;
  DoubleTab& variable_inconnue = equation().inconnue().valeurs();
  int nb_comp = variable_inconnue.dimension(1);
  variable_imposee_mod.echange_espace_virtuel();
  DoubleTab variable_imposee_sommet(nb_som,nb_comp);
  IntTab sommet_interp(1, nb_som);
  variable_imposee_calculee = 0.0;
  variable_imposee_sommet = 0.0;
  sommet_interp = 0;
  for (int i = 0; i < nb_som; i++)
    {
      if (is_dirichlet(i) > 0.0)
        {
          sommet_interp(0, i) = 1;
          double d2 = 0.0;
          for(int j = 0; j < dim_esp; j++)
            {
              double x = dom.coord(i,j);
              double xp = solid_points(i,j);
              d2 += (x - xp)*(x - xp);
            }
          for(int j = 0; j <nb_comp ; j++) variable_imposee_sommet(i,j) = variable_imposee_mod(i,j);
          d2 = sqrt(d2);
          if (d2 > eps)
            {
              IntList& voisins = interp.getSommetsVoisinsOf(i);
              if (!(voisins.est_vide()))
                {
                  ArrOfDouble mean_grad(nb_comp);
                  mean_grad = 0.0;
                  int taille = voisins.size();
                  int nb_contrib = 0;
                  for (int k = 0; k < taille; k++)
                    {
                      int num_som = voisins[k];
                      double d1 = 0.0;
                      for (int j = 0; j < dim_esp; j++)
                        {
                          double xf = dom.coord(num_som,j);
                          double xpf = solid_points(num_som,j);
                          d1 += (xf - xpf)*(xf - xpf);
                        }
                      d1 = sqrt(d1);
                      if (d1 > eps)
                        {
                          for (int j = 0; j < nb_comp; j++)
                            {
                              double vpf = variable_imposee_mod(num_som,j);
                              double vf = variable_inconnue(num_som,j);
                              mean_grad[j] += (vf - vpf)/d1;
                            }
                          nb_contrib += 1;
                        }
                    }
                  for (int j = 0; j < nb_comp; j++)
                    {
                      mean_grad[j] /= nb_contrib;
                      variable_imposee_sommet(i,j) += mean_grad[j]*d2;
                    }
                }
            }
        }
      else
        {
          for(int j = 0; j < nb_comp; j++)
            {
              variable_imposee_sommet(i,j) = variable_imposee_mod(i,j);
            }
        }
    }
  for (int f = 0; f < nb_faces; f++)
    {
      int compteur = 0;
      int som_f = 0;
      while ((som_f < nb_som_face) and (sommet_interp(0, faces_sommets(f,som_f))==1))
        {
          compteur++;
          som_f++;
        }
      if (compteur == nb_som_face)
        {
          for (int som_=0; som_ < nb_som_face; som_++)
            {
              int i = faces_sommets(f,som_);
              for(int j = 0; j < nb_comp; j++)
                {
                  variable_imposee_calculee(f,j) += variable_imposee_sommet(i,j)/nb_som_face;
                }
            }
        }
    }
}

void Source_PDF_VEF::calculer_variable_imposee_hybrid()
{
  Cerr << "on passe ici" << endl;
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  int nb_faces=domaine_VEF.nb_faces();
  int dim_esp = Objet_U::dimension;
  const Domaine& dom = domaine_VEF.domaine();
  const IntTab& faces_sommets=domaine_VEF.face_sommets();
  int nb_som_face=faces_sommets.dimension(1);
  int nb_som = domaine_VEF.nb_som();
  DoubleTab& variable_imposee_mod = modele_lu_.variable_imposee_->valeurs();
  DoubleTab& variable_imposee_calculee = variable_imposee_;
  Interpolation_IBM_hybrid& interp = ref_cast(Interpolation_IBM_hybrid,interpolation_lue_.valeur());
  DoubleTab& solid_points = interp.solid_points_->valeurs();
  DoubleTab& fluid_elems = interp.fluid_elems_->valeurs();
  DoubleTab& fluid_points = interp.fluid_points_->valeurs();
  double eps = 1e-12;
  Champ_P1NC& champ_variable_inconnue = ref_cast(Champ_P1NC,equation().inconnue());
  DoubleTab& variable_inconnue = champ_variable_inconnue.valeurs();
  int nb_comp = variable_inconnue.dimension(1);
  DoubleTab xf(1, dim_esp);
  DoubleTab vf(1, nb_comp);
  ArrOfInt cells(1);
  DoubleTab variable_imposee_sommet(nb_som,nb_comp);
  IntTab sommet_interp(1, nb_som);
  variable_imposee_calculee = 0.0;
  variable_imposee_sommet = 0.0;
  sommet_interp = 0;
  Cerr << "on passe ici" << endl;
  for (int i = 0; i < nb_som; i++)
    {
      Cerr << "on passe la" << endl;
      if (fluid_elems(i) >= 0.0)
        {
          sommet_interp(0, i) = 1;
          double d1 = 0.0;
          double d2 = 0.0;
          for(int j = 0; j < dim_esp; j++)
            {
              double xj = dom.coord(i,j);
              double xjf = fluid_points(i,j);
              xf(0, j) = xjf;
              double xjs = solid_points(i,j);
              d1 += (xj-xjs)*(xj-xjs);
              d2 += (xjf-xj)*(xjf-xj);
            }
          d1 = sqrt(d1);
          d2 = sqrt(d2);
          double inv_d = 1.0 / (d1 + d2);
          cells[0] = int(fluid_elems(i));
          champ_variable_inconnue.valeur_a_elem(xf, vf, cells[0]);
          for (int j = 0; j < nb_comp; j++)
            {
              double vjs = variable_imposee_mod(i,j);
              variable_imposee_sommet(i,j) = vjs + (vf(0, j)-vjs)*inv_d*d1;
            }
        }
      else if (fluid_elems(i) > -1.5 && fluid_elems(i) < -0.5)
        {
          sommet_interp(0, i) = 1;
          double d2 = 0.0;
          for(int j = 0; j < dim_esp; j++)
            {
              double x = dom.coord(i,j);
              double xp = solid_points(i,j);
              d2 += (x - xp)*(x - xp);
            }
          d2 = sqrt(d2);
          for(int j = 0; j <nb_comp ; j++) variable_imposee_sommet(i,j) = variable_imposee_mod(i,j);
          if (d2 > eps)
            {
              IntList& voisins = interp.getSommetsVoisinsOf(i);
              if (!(voisins.est_vide()))
                {
                  ArrOfDouble mean_grad(nb_comp);
                  mean_grad = 0.0;
                  int taille = voisins.size();
                  int nb_contrib = 0;
                  for (int k = 0; k < taille; k++)
                    {
                      int num_som = voisins[k];
                      double d1 = 0.0;
                      for (int j = 0; j < dim_esp; j++)
                        {
                          double xjf = dom.coord(num_som,j);
                          double xpf = solid_points(num_som,j);
                          d1 += (xjf - xpf)*(xjf - xpf);
                        }
                      d1 = sqrt(d1);
                      if (d1 > eps)
                        {
                          for (int j = 0; j < nb_comp; j++)
                            {
                              double vpf = variable_imposee_mod(num_som,j);
                              double vjf = variable_inconnue(num_som,j);
                              mean_grad[j] += (vjf - vpf)/d1;
                            }
                          nb_contrib += 1;
                        }
                    }
                  for (int j = 0; j < nb_comp; j++)
                    {
                      mean_grad[j] /= nb_contrib;
                      variable_imposee_sommet(i,j) += mean_grad[j]*d2;
                    }
                }
            }
        }
      else
        {
          for(int j = 0; j < nb_comp; j++)
            {
              variable_imposee_sommet(i,j) = variable_imposee_mod(i,j);
            }
        }
    }
  for (int f = 0; f < nb_faces; f++)
    {
      int compteur = 0;
      int som_f = 0;
      while ((som_f < nb_som_face) and (sommet_interp(0, faces_sommets(f,som_f))==1))
        {
          compteur++;
          som_f++;
        }
      if (compteur == nb_som_face)
        {
          for (int som_=0; som_ < nb_som_face; som_++)
            {
              int i = faces_sommets(f,som_);
              for(int j = 0; j < nb_comp; j++)
                {
                  variable_imposee_calculee(f,j) += variable_imposee_sommet(i,j)/nb_som_face;
                }
            }
        }
    }
}

/*##################################################################################################
####################################################################################################
################################# calculer_vitesse_imposee #########################################
####################################################################################################
##################################################################################################*/

void Source_PDF_VEF::calculer_vitesse_imposee_power_law_tbl()
{
  assert(Objet_U::dimension==3);
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  int nb_faces=domaine_VEF.nb_faces();
  int nb_comp = Objet_U::dimension;
  const Domaine& dom = domaine_VEF.domaine();
  const IntTab& faces_sommets=domaine_VEF.face_sommets();
  int nb_som_face=faces_sommets.dimension(1);
  int nb_som=domaine_VEF.nb_som();
  DoubleTab& vitesse_imposee_mod = modele_lu_.variable_imposee_->valeurs();
  DoubleTab& vitesse_imposee_calculee = variable_imposee_;
  Interpolation_IBM_power_law_tbl& interp = ref_cast(Interpolation_IBM_power_law_tbl,interpolation_lue_.valeur());
  int form_WJSP = interp.get_formulation_WJSP();
  DoubleTab& fluid_points = interp.fluid_points_->valeurs();
  DoubleTab& solid_points = interp.solid_points_->valeurs();
  DoubleTab& fluid_elems = interp.fluid_elems_->valeurs();
  double A_pwl = interp.get_A_pwl(form_WJSP);
  double B_pwl = interp.get_B_pwl();
  double y_c_p_pwl = 0;
  double C_pwl_WJSP = 0;
  double D_pwl_WJSP = 0;
  double p_pwl_WJSP = 0;
  double y_c1_p_pwl_WJSP = 0;
  double y_c2_p_pwl_WJSP = 0;
  if (!form_WJSP)
    {
      y_c_p_pwl = interp.get_y_c_p_pwl();
    }
  else
    {
      C_pwl_WJSP = interp.get_C_pwl_WJSP();
      D_pwl_WJSP = interp.get_D_pwl_WJSP();
      p_pwl_WJSP = interp.get_p_pwl_WJSP();
      y_c1_p_pwl_WJSP = interp.get_y_c1_p_pwl_WJSP();
      y_c2_p_pwl_WJSP = interp.get_y_c2_p_pwl_WJSP();
    }
  int impr_yplus = interp.get_impr() ;
  Champ_P1NC& champ_vitesse_inconnue = ref_cast(Champ_P1NC,equation().inconnue());
  DoubleTab& val_vitesse_inconnue = champ_vitesse_inconnue.valeurs();
  if (nb_comp != val_vitesse_inconnue.dimension(1))
    {
      Cerr<<"Source_PDF_VEF::calculer_vitesse_imposee_power_law_tbl: variable must be a vector of dimension "<<Objet_U::dimension<<" ! "<<finl;
      abort();
    }
  double form_lin_pwl = interp.get_formulation_linear_pwl();
  DoubleTab xf(1, nb_comp);
  DoubleTab vf(1, nb_comp);
  ArrOfInt cells(1);
  DoubleTab vitesse_imposee_sommet(nb_som,nb_comp);
  IntTab sommet_interp(1, nb_som);
  vitesse_imposee_calculee = 0.0;
  vitesse_imposee_sommet = 0.0;
  sommet_interp = 0;
  // operateur(0) : diffusivite
  if (equation().nombre_d_operateurs()<1)
    {
      Cerr << "Source_PDF_VEF : nombre_d_operateurs = "<<equation().nombre_d_operateurs()<<" < 1"<<finl;
      exit();
    }

  int flag = 0;
  const Op_Diff_VEF_base& opdiffu = ref_cast(Op_Diff_VEF_base,equation().operateur(0).l_op_base());
  flag = opdiffu.diffusivite().valeurs().dimension(0)>1 ? 1 : 0;

  double eps = 1e-12;
  double d1_min = 1.0e+10;
  double d1_max = 0.0;
  double d1_mean = 0.0;
  double yplus_min = 1.0e+10;
  double yplus_max = 0.0;
  double yplus_mean = 0.0;
  double u_tau_min = 1.0e+10;
  double u_tau_max = 0.0;
  double u_tau_mean = 0.0;
  double yplus_ref_min = 1.0e+10;
  double yplus_ref_max = 0.0;
  double yplus_ref_mean = 0.0;
  double u_tau_ref_min = 1.0e+10;
  double u_tau_ref_max = 0.0;
  double u_tau_ref_mean = 0.0;
  double h_yplus_min = 1.0e+10;
  double h_yplus_max = 0.0;
  double h_yplus_mean = 0.0;
  int yplus_count=0 ;
  double y_plus=0;

  int N = interp.get_N_histo();
  DoubleTab tab_h(1, nb_som);
  DoubleTab abs_h(1, N+1);
  DoubleTab compteur_h(1, N);
  compteur_h = 0.;
  tab_h = 0.;

  tab_u_star_ibm_ = 0.;
  tab_y_plus_ibm_ = 0.;

  for (int i = 0; i < nb_som; i++)
    {
      for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = vitesse_imposee_mod(i,j);
      int itisok = 1;
      if ((fluid_elems(i) >= 0.0))// and (is_dirichlet(i) == 1))
        {
          sommet_interp(0, i) = 1;
          double d1 = 0.0;
          double d2 = 0.0;
          DoubleTab normale(1, nb_comp);
          double norme_de_la_normale= 0.0;
          for(int j = 0; j < nb_comp; j++)
            {
              double xj = dom.coord(i,j);
              double xjf = fluid_points(i,j);
              xf(0, j) = xjf;
              double xjs = solid_points(i,j);
              d1 += (xj-xjs)*(xj-xjs);
              d2 += (xjf-xj)*(xjf-xj);
              normale(0,j) = xjf - xjs;
              norme_de_la_normale += normale(0,j)*normale(0,j);
            }
          d1 = sqrt(d1);
          d2 = sqrt(d2);
          double y_ref= d1+d2;

          // traitement des exceptions
          if (d1 < eps) itisok = 0; //le point P est sur la frontiere immergee : affectation
          norme_de_la_normale = sqrt(norme_de_la_normale);
          if ( norme_de_la_normale > eps )
            for(int j = 0; j < nb_comp; j++) normale(0,j) /= norme_de_la_normale; // la on met la norme unite tout en gardant la direction, on normalise quoi
          else
            {
              for(int j = 0; j < nb_comp; j++) normale(0,j) = 0.; // le point fluide est sur la frontiere immergee : affectation
              itisok = 0;
            }
          if ( y_ref < eps ) // le point fluide est sur la frontiere immergee : affectation
            {
              y_ref = 1.;
              itisok = 0;
            }
          // calcul vitesse power law
          if (itisok)
            {
              cells(0) = int(fluid_elems(i));
              champ_vitesse_inconnue.valeur_a_elem(xf, vf, cells[0]); // vf la vitesse totale interpolee au pt fluide
              double Vn = 0.;
              for(int j = 0; j < nb_comp; j++) Vn += vf(0, j) * normale(0,j);
              DoubleTab v_ref_t(1, nb_comp);

              for(int j = 0; j < nb_comp; j++) v_ref_t(0, j) =vf(0, j) - Vn*normale(0,j);

              double norme_v_ref_t = 0.;
              for(int j = 0; j < nb_comp; j++) norme_v_ref_t +=  v_ref_t(0, j)*v_ref_t(0,j)  ;
              norme_v_ref_t = sqrt ( norme_v_ref_t );
              if (norme_v_ref_t < 1.0e-6) norme_v_ref_t = 1.e-6;

              double nu = (flag ? opdiffu.diffusivite().valeurs()(cells) : opdiffu.diffusivite().valeurs()(0,0));
              // On calcule U_tau_ref pour pouvoir calculer les quantites adimensionees ( y+ ...)
              // Hypothese : sous-couche puissance
              double u_tau_ref = pow ( norme_v_ref_t , (1/(1+B_pwl)) ) * pow ( A_pwl, (-1/(1+B_pwl)) )  * pow ( y_ref , (-B_pwl/(1+B_pwl)) ) * pow ( nu,(B_pwl/(1+B_pwl)) ) ;
              // Cerr<<"u_tau_ref y_ref   nu ="<<u_tau_ref<<" "<<y_ref<<" "<<nu<<finl;;
              double y_ref_p = y_ref * u_tau_ref  / nu;  //la on a enfin tout ce qu'il faut pour etablir la loi polynomiale pour y r+
              double test_ref;

              if (!form_WJSP)
                {
                  if ( y_ref_p > y_c_p_pwl)  // a partir de la commence l'expression de la loi de paroi polynomiale turbulente
                    {
                      if (!form_lin_pwl)
                        {
                          for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  v_ref_t(0,j) * pow ( d1 / y_ref , B_pwl )  ;
                        }
                      else // Formulation lineaire de la loi polynomilale ------------------------------
                        {
                          for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  v_ref_t(0,j) *(1. - B_pwl*(1-d1/y_ref)) ;
                        }
                      test_ref = 1.;
                      // Cerr << "zone log/sous-couche inertielle en coherence avec le calcul de u tau" << finl;
                    }
                  else
                    {
                      // En fait, on est en sous-couche visqueuse
                      if (!form_lin_pwl)
                        {
                          for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  v_ref_t(0,j) * ( d1 / y_ref )   ;
                        }
                      else
                        {
                          for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = v_ref_t(0,j) * (1. - pow(u_tau_ref, 2.)*(y_ref-d1)/(nu *norme_v_ref_t));
                        }
                      test_ref = -1. ;
                      u_tau_ref = pow ( (nu * norme_v_ref_t / y_ref) , 0.5 ) ; // on recalcule u_tau et y_plus en lineaire au pt reference
                      y_ref_p = y_ref * u_tau_ref / nu;
                      if ( y_ref_p > y_c_p_pwl )
                        {
                          // Incoherence : on n utilise pas de lois de paroi
                          itisok = 0;
                        }
                      // Cerr << "zone lineaire/sous-couche visqueuse" << finl;
                    }                   // Fin LdPturb

                  double norme_v_imp_c = 0;
                  for(int j=0 ; j < nb_comp; j++) norme_v_imp_c += vitesse_imposee_sommet(i ,j) * vitesse_imposee_sommet(i ,j);
                  norme_v_imp_c = sqrt( norme_v_imp_c );

                  double u_tau = pow ( norme_v_imp_c , (1/(1+B_pwl)) ) * pow ( A_pwl, (-1/(1+B_pwl)) )  * pow ( d1 , (-B_pwl/(1+B_pwl)) ) * pow ( nu,(B_pwl/(1+B_pwl)) ) ; // Hypothese : sous-couche puissance
                  y_plus = u_tau * d1 / nu;
                  double test_;
                  if (y_plus >y_c_p_pwl)
                    {
                      test_ = 1.;
                      // Cerr<<"u_tau ="<<u_tau<<" y+ = "<<y_plus<<" y+ref = "<<y_ref_p<<finl;
                    }
                  else
                    {
                      test_ = -1.;
                      u_tau = pow ( (nu * norme_v_imp_c/ d1) , 0.5 ) ; // on recalcule u_tau et y_plus en lineaire au pt force
                      y_plus = u_tau * d1 / nu;
                      if ( y_plus > y_c_p_pwl )
                        {
                          // Incoherence : on n utilise pas de lois de paroi
                          Cerr<<"INCOHERENCE: u_tau ="<<u_tau<<" y+ = "<<y_plus<<" y+ref = "<<y_ref_p<<finl;
                          itisok = 0;
                        }
                    }
                  // Cerr<<"u_tau d1   nu ="<<u_tau<<" "<<d1<<" "<<nu<<finl;;

                  // ici on effectue le test pour savoir si le noeud de frontiere et le point fluide se trouvent dans la meme zone: si oui ok, si non on impose la vitesse

                  // Traitement des exceptions
                  if ( (test_ * test_ref < 0) || (itisok == 0) )  //si non on affecte la vitesse paroi
                    {
                      for(int j = 0; j < nb_comp; j++)  vitesse_imposee_sommet(i,j) = vitesse_imposee_mod(i,j);

                      // Cerr << "erreur" << finl;
                      itisok = 0;
                    }
                  else
                    {
                      tab_u_star_ibm_(i) = u_tau;
                      tab_y_plus_ibm_(i) = y_plus;
                    }
                  if (impr_yplus && itisok && (indicateur_nodal_champ_aire_(i)==1.))
                    {
                      if (d1 > d1_max) d1_max = d1;
                      if (d1 < d1_min) d1_min = d1;
                      d1_mean +=  d1;
                      if (y_plus > yplus_max) yplus_max = y_plus;
                      if (y_plus < yplus_min) yplus_min = y_plus;
                      yplus_mean +=  y_plus;
                      if (u_tau > u_tau_max) u_tau_max = u_tau;
                      if (u_tau < u_tau_min) u_tau_min = u_tau;
                      u_tau_mean += u_tau;
                      if (y_ref_p > yplus_ref_max) yplus_ref_max = y_ref_p;
                      if (y_ref_p < yplus_ref_min) yplus_ref_min = y_ref_p;
                      yplus_ref_mean += y_ref_p;
                      if (u_tau_ref > u_tau_ref_max) u_tau_ref_max = u_tau_ref;
                      if (u_tau_ref < u_tau_ref_min) u_tau_ref_min = u_tau_ref;
                      u_tau_ref_mean += u_tau_ref;
                      // Distribution des y+ au voisinage des parois
                      if (y_plus > h_yplus_max) h_yplus_max = y_plus;
                      if (y_plus < h_yplus_min) h_yplus_min = y_plus;
                      h_yplus_mean += y_plus;
                      tab_h(0,yplus_count) = y_plus;
                      yplus_count += 1;
                    }
                }
              else
                {
                  //double u_tau;
                  if ( y_ref_p > y_c2_p_pwl_WJSP)  // a partir de la commence l'expression de la loi de paroi polynomiale turbulente
                    {
                      y_plus = u_tau_ref * d1/nu;
                      if (y_plus > y_c2_p_pwl_WJSP)
                        {
                          if (!form_lin_pwl)
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  v_ref_t(0,j) * pow ( d1 / y_ref , B_pwl )  ;
                            }
                          else // Formulation lineaire de la loi polynomilale ------------------------------
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  v_ref_t(0,j) * (1. - B_pwl*(1-d1/y_ref)) ;
                            }
                        }
                      else if (y_plus < y_c2_p_pwl_WJSP and y_plus > y_c1_p_pwl_WJSP)
                        {
                          if (!form_lin_pwl)
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = (C_pwl_WJSP*pow(y_plus,2*p_pwl_WJSP-1) *u_tau_ref + D_pwl_WJSP*pow(y_plus,p_pwl_WJSP-1) *u_tau_ref)* v_ref_t(0,j)/norme_v_ref_t;
                            }
                          else // Formulation lineaire de la loi polynomilale ------------------------------
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = 0 ;
                            }
                        }
                      else
                        {
                          if (!form_lin_pwl)
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  y_plus * u_tau_ref * v_ref_t(0,j)/norme_v_ref_t;
                            }
                          else // Formulation lineaire de la loi polynomilale ------------------------------
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  v_ref_t(0,j) * (1. - B_pwl*(1-d1/y_ref)) ;
                            }
                        }
                      // Cerr << "zone log/sous-couche inertielle en coherence avec le calcul de u tau" << finl;
                    }
                  else // hypothèse buffer layer
                    {
                      u_tau_ref = (nu/y_ref) * pow(-(D_pwl_WJSP/(2*C_pwl_WJSP)) * (1 + pow(1 + 4 * C_pwl_WJSP/pow(D_pwl_WJSP,2) * norme_v_ref_t * y_ref /nu ,1/2)),1/p_pwl_WJSP) ;
                      y_ref_p = y_ref * u_tau_ref  / nu;
                      if (y_ref_p < y_c2_p_pwl_WJSP and y_ref_p > y_c1_p_pwl_WJSP)
                        {
                          y_plus = u_tau_ref * d1/nu;
                          if (y_plus < y_c2_p_pwl_WJSP and y_plus > y_c1_p_pwl_WJSP)
                            {
                              if (!form_lin_pwl)
                                {
                                  double phi = -pow(1 + 4 * C_pwl_WJSP/pow(D_pwl_WJSP,2) *  norme_v_ref_t * y_ref /nu,1/2) - 1;
                                  double u_norm = pow(D_pwl_WJSP,2) * nu / (4 * C_pwl_WJSP * y_ref) * (pow ( d1 / y_ref , 2*p_pwl_WJSP - 1 ) * pow(phi,2) + 2 * pow ( d1 / y_ref , p_pwl_WJSP - 1 ) * phi);
                                  for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  u_norm * v_ref_t(0,j)/norme_v_ref_t;
                                }
                              else // Formulation lineaire de la loi polynomilale ------------------------------
                                {
                                  for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  0 ; //à faire
                                }
                            }
                          else
                            {
                              if (!form_lin_pwl)
                                {
                                  for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  y_plus * u_tau_ref * v_ref_t(0,j)/norme_v_ref_t  ;
                                }
                              else // Formulation lineaire de la loi polynomilale ------------------------------
                                {
                                  for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  0 ;
                                }
                            }
                        }
                      else
                        {
                          // En fait, on est en sous-couche visqueuse
                          if (!form_lin_pwl)
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) =  v_ref_t(0,j) * ( d1 / y_ref )   ;
                            }
                          else
                            {
                              for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = v_ref_t(0,j) * (1. - pow(u_tau_ref, 2.)*(y_ref-d1)/(nu *norme_v_ref_t));
                            }
                          u_tau_ref = pow ( (nu * norme_v_ref_t / y_ref) , 0.5 ) ; // on recalcule u_tau et y_plus en lineaire au pt reference
                          y_ref_p = y_ref * u_tau_ref / nu;
                          if ( y_ref_p > y_c1_p_pwl_WJSP)
                            {
                              // Incoherence : on n utilise pas de lois de paroi
                              itisok = 0;
                              Cerr << "Incohérence" << endl;
                            }
                        }
                    }
                  // Cerr<<"u_tau d1   nu ="<<u_tau<<" "<<d1<<" "<<nu<<finl;;

                  // ici on effectue le test pour savoir si le noeud de frontiere et le point fluide se trouvent dans la meme zone: si oui ok, si non on impose la vitesse

                  // Traitement des exceptions
                  if (itisok == 0)  //si non on affecte la vitesse paroi
                    {
                      for(int j = 0; j < nb_comp; j++)  vitesse_imposee_sommet(i,j) = vitesse_imposee_mod(i,j);

                      // Cerr << "erreur" << finl;
                      itisok = 0;
                    }
                  else
                    {
                      tab_u_star_ibm_(i) = u_tau_ref;
                      tab_y_plus_ibm_(i) = y_plus;
                    }
                  if (impr_yplus && itisok && (indicateur_nodal_champ_aire_(i)==1.))
                    {
                      if (d1 > d1_max) d1_max = d1;
                      if (d1 < d1_min) d1_min = d1;
                      d1_mean +=  d1;
                      if (y_plus > yplus_max) yplus_max = y_plus;
                      if (y_plus < yplus_min) yplus_min = y_plus;
                      yplus_mean +=  y_plus;
                      if (u_tau_ref > u_tau_max) u_tau_max = u_tau_ref;
                      if (u_tau_ref < u_tau_min) u_tau_min = u_tau_ref;
                      u_tau_mean += u_tau_ref;
                      if (y_ref_p > yplus_ref_max) yplus_ref_max = y_ref_p;
                      if (y_ref_p < yplus_ref_min) yplus_ref_min = y_ref_p;
                      yplus_ref_mean += y_ref_p;
                      if (u_tau_ref > u_tau_ref_max) u_tau_ref_max = u_tau_ref;
                      if (u_tau_ref < u_tau_ref_min) u_tau_ref_min = u_tau_ref;
                      u_tau_ref_mean += u_tau_ref;
                      // Distribution des y+ au voisinage des parois
                      if (y_plus > h_yplus_max) h_yplus_max = y_plus;
                      if (y_plus < h_yplus_min) h_yplus_min = y_plus;
                      h_yplus_mean += y_plus;
                      tab_h(0,yplus_count) = y_plus;
                      yplus_count += 1;
                    }
                }
            }
        }
    }
  get_champ("u_star_ibm");
  get_champ("y_plus_ibm");
  for (int f = 0; f < nb_faces; f++)
    {
      int compteur = 0;
      int som_f = 0;
      while ((som_f < nb_som_face) and (sommet_interp(0, faces_sommets(f,som_f))==1))
        {
          compteur++;
          som_f++;
        }
      if (compteur == nb_som_face)
        {
          for (int som_=0; som_ < nb_som_face; som_++)
            {
              int i = faces_sommets(f,som_);
              for(int j = 0; j < nb_comp; j++)
                {
                  vitesse_imposee_calculee(f,j) += vitesse_imposee_sommet(i,j)/nb_som_face;
                }
            }
        }
    }

  if (impr_yplus && (yplus_count >= 1) )
    {
      d1_mean /= yplus_count;
      yplus_mean /= yplus_count;
      yplus_ref_mean /= yplus_count;
      u_tau_mean /= yplus_count;
      u_tau_ref_mean /= yplus_count;
      Cerr<<"min mean max y  = "<<d1_min<<" "<<d1_mean<<" "<<d1_max<<finl;
      Cerr<<"min mean max y+ = "<<yplus_min<<" "<<yplus_mean<<" "<<yplus_max<<finl;
      Cerr<<"min mean u_tau = "<<u_tau_min<<" "<<u_tau_mean<<" "<<u_tau_max<<finl;
      Cerr<<"min mean max y+_ref = "<<yplus_ref_min<<" "<<yplus_ref_mean<<" "<<yplus_ref_max<<finl;
      Cerr<<"min mean u_tau_ref = "<<u_tau_ref_min<<" "<<u_tau_ref_mean<<" "<<u_tau_ref_max<<finl;

      h_yplus_mean /= yplus_count;
      double range = ( h_yplus_max - h_yplus_min)/N;
      for (int k = 0; k < N+1; k++) abs_h(0,k) = (h_yplus_min + k*range);
      for (int i=0; i < yplus_count; i++)
        {
          for (int j = 0; j < N; j++)
            {
              if ( abs_h(0,j) < tab_h(0,i) &&  tab_h(0,i) <= abs_h(0,j+1) ) compteur_h(0,j) +=1;
            }
        }
      Cerr<<"histogramme y+ compteur =";
      for (int j = 0; j < N; j++) Cerr<<" "<<compteur_h(0,j)<<" ";
      Cerr<<finl;
      Cerr<<"histogramme y+ abscisse =";
      for (int j = 0; j < N+1; j++) Cerr<<" "<<abs_h(0,j)<<" ";
      Cerr<<finl;
      Cerr<<"histogramme y+ : min = "<<h_yplus_min<<" - mean = "<<h_yplus_mean<<" - max =  "<<h_yplus_max<<finl;

      Cerr<<"Moyenne sur "<<yplus_count<<" points"<<finl;
    }
}

void Source_PDF_VEF::calculer_vitesse_imposee_power_law_tbl_u_star()
{
  assert(Objet_U::dimension==3);
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  int nb_faces=domaine_VEF.nb_faces();
  int nb_comp = Objet_U::dimension;
  const Domaine& dom = domaine_VEF.domaine();
  const IntTab& faces_sommets=domaine_VEF.face_sommets();
  int nb_som_face=faces_sommets.dimension(1);
  int nb_som = domaine_VEF.nb_som();
  DoubleTab& vitesse_imposee_mod = modele_lu_.variable_imposee_->valeurs();
  DoubleTab& vitesse_imposee_calculee = variable_imposee_;
  Interpolation_IBM_power_law_tbl_u_star& interp = ref_cast(Interpolation_IBM_power_law_tbl_u_star,interpolation_lue_.valeur());
  DoubleTab& is_dirichlet = interp.is_dirichlet_->valeurs();
  DoubleTab& solid_points = interp.solid_points_->valeurs();
  DoubleTab& solid_elems = interp.solid_elems_->valeurs();
  double A_pwl = interp.get_A_pwl(0);
  double B_pwl = interp.get_B_pwl();
  double y_c_p_pwl = interp.get_y_c_p_pwl();
  int impr_yplus = interp.get_impr() ;
  DoubleTab& vitesse_inconnue = equation().inconnue().valeurs();
  DoubleTab vitesse_imposee_sommet(nb_som,nb_comp);
  IntTab sommet_interp(1, nb_som);
  vitesse_imposee_calculee = 0.0;
  vitesse_imposee_sommet = 0.0;
  sommet_interp = 0;
  if (nb_comp != vitesse_inconnue.dimension(1))
    {
      Cerr<<"Source_PDF_EF::calculer_vitesse_imposee_power_law_tbl_u_star: variable must be a vector of dimension "<<Objet_U::dimension<<" ! "<<finl;
      abort();
    }

  ArrOfInt cells(1);

  // operateur(0) : diffusivite
  if (equation().nombre_d_operateurs()<1)
    {
      Cerr << "Source_PDF_VEF : nombre_d_operateurs = "<<equation().nombre_d_operateurs()<<" < 1"<<finl;
      exit();
    }

  int flag = 0;
  const Op_Diff_VEF_base& opdiffu = ref_cast(Op_Diff_VEF_base,equation().operateur(0).l_op_base());
  flag = opdiffu.diffusivite().valeurs().dimension(0)>1 ? 1 : 0;

  double eps = 1e-12;
  double d1_min = 1.0e+10;
  double d1_max = 0.0;
  double d1_mean = 0.0;
  double yplus_min = 1.0e+10;
  double yplus_max = 0.0;
  double yplus_mean = 0.0;
  double u_tau_min = 1.0e+10;
  double u_tau_max = 0.0;
  double u_tau_mean = 0.0;
  double yplus_ref_min = 1.0e+10;
  double yplus_ref_max = 0.0;
  double yplus_ref_mean = 0.0;
  double u_tau_ref_min = 1.0e+10;
  double u_tau_ref_max = 0.0;
  double u_tau_ref_mean = 0.0;
  double h_yplus_min = 1.0e+10;
  double h_yplus_max = 0.0;
  double h_yplus_mean = 0.0;
  int yplus_count=0 ;
  int yplus_ref_count=0 ;

  int N = interp.get_N_histo();
  DoubleTab tab_h(1, nb_som);
  DoubleTab abs_h(1, N+1);
  DoubleTab compteur_h(1, N);
  compteur_h = 0.;
  tab_h = 0.;

  vitesse_imposee_mod.echange_espace_virtuel();
  tab_u_star_ibm_ = 0.;
  tab_y_plus_ibm_ = 0.;
  // DoubleTrav nb_vois(is_dirichlet);

  for (int i = 0; i < nb_som; i++)
    {
      double x = 0.;
      if (is_dirichlet(i) > 0.0)
        {
          sommet_interp(0, i) = 1;
          int itisok_P = 1;
          double d1P = 0.0;
          DoubleTab normaleP(1, nb_comp);
          double norme_de_la_normaleP= 0.0;
          //  normal et distance normale au point P
          for(int j = 0; j < nb_comp; j++)
            {
              vitesse_imposee_sommet(i,j) = vitesse_imposee_mod(i,j);
              x = dom.coord(i,j);
              double xp = solid_points(i,j);
              normaleP(0,j) = x - xp;
              norme_de_la_normaleP += normaleP(0,j)*normaleP(0,j);
            }
          norme_de_la_normaleP = sqrt(norme_de_la_normaleP);
          d1P = norme_de_la_normaleP;
          if (d1P > eps)
            for(int j = 0; j < nb_comp; j++) normaleP(0,j) /= norme_de_la_normaleP;
          else normaleP = 0.;

          cells(0) = int(solid_elems(i)); //pour definir la viscosite laminaire

          double nu = (flag ? opdiffu.diffusivite().valeurs()(cells) : opdiffu.diffusivite().valeurs()(0,0));

          double u_tau = 0.;
          double y_plus = 0.;
          if (d1P > eps)
            //le point P n'est pas sur la frontiere immergee
            {
              IntList& voisins = interp.getSommetsVoisinsOf(i);
              if (!(voisins.est_vide()))
                // il y a une liste de dof voisins
                {
                  double u_tau_ref;
                  double y_plus_ref;
                  DoubleTab mean_u_tau_neighbour(1, nb_comp); // moyenne de u_tau vectoriel sur les dof voisins
                  mean_u_tau_neighbour = 0.0;
                  double mean_y_plus_neighbour = 0.; // moyenne de y_plus sur les dof voisins
                  int taille = voisins.size();
                  // nb_vois(i) = taille * 1.0 ;
                  double pond_tot = 0.;
                  double ponderation_k;
                  for (int k = 0; k < taille; k++)
                    // boucle sur les dof voisins
                    {
                      int num_som = voisins[k];
                      int itisok = 1;
                      double d1 = 0.0;
                      DoubleTab a(1, nb_comp);
                      DoubleTab b(1, nb_comp);
                      double norme_a = 0.0;
                      double norme_b = 0.0;

                      for (int j = 0; j < nb_comp; j++)
                        {
                          double xf = dom.coord(num_som,j);
                          double xpf = solid_points(num_som,j);
                          d1 += (xf - xpf)*(xf - xpf);
                          a(0,j) = (xf - x);
                          b(0,j) = (xf - x)*normaleP(0,j);
                          norme_a += a(0,j)*a(0,j);
                          norme_b += b(0,j)*b(0,j);
                        }
                      d1 = sqrt(d1);
                      if ( d1 < eps ) itisok = 0;

                      if ( (norme_a - norme_b) > eps )
                        ponderation_k = 1/(sqrt(norme_a - norme_b));
                      else
                        ponderation_k = 1./eps;

                      u_tau_ref = 0.;
                      y_plus_ref = 0.;
                      if (itisok)
                        {
                          // normale et tangente vitesse vis a vis de normalP
                          double Vn = 0.;
                          for(int j = 0; j < nb_comp; j++) Vn += vitesse_inconnue(num_som,j) * normaleP(0,j);
                          DoubleTab vtf_k(1, nb_comp);
                          for(int j = 0; j < nb_comp; j++) vtf_k(0, j) = vitesse_inconnue(num_som, j) - Vn * normaleP(0,j);
                          double norme_vtf_k = 0.;
                          for(int j = 0; j < nb_comp; j++) norme_vtf_k += vtf_k(0, j)*vtf_k(0,j);
                          norme_vtf_k = sqrt(norme_vtf_k);
                          u_tau_ref = pow ( norme_vtf_k , (1/(1+B_pwl)) ) * pow ( A_pwl, (-1/(1+B_pwl)) ) * pow ( d1 , (-B_pwl/(1+B_pwl)) ) * pow ( nu,(B_pwl/(1+B_pwl)) ) ;  // vitesse de frottement power law au point fluide k
                          y_plus_ref = d1 * u_tau_ref / nu;
                          if (y_plus_ref < y_c_p_pwl ) // En fait, on est en lois lineaire !
                            {
                              u_tau_ref = pow ( (nu * norme_vtf_k / d1) , 0.5 ) ; // on recalcule u_tau_ref et y_plus_ref en lineaire au point fluide k
                              y_plus_ref = d1 * u_tau_ref / nu;
                              if ( y_plus_ref > y_c_p_pwl )
                                {
                                  // Incoherence : on n utilise pas de lois de paroi
                                  itisok = 0;
                                }
                            }

                          // On moyenne les vecteurs vitesses de frottement tangents (vis a vis de P) par l'inverse des distances
                          if (norme_vtf_k < eps) itisok = 0;
                          if (itisok)
                            {
                              for(int j = 0; j < nb_comp; j++) mean_u_tau_neighbour(0,j) += (u_tau_ref * vtf_k(0, j) / norme_vtf_k) * ponderation_k;
                              mean_y_plus_neighbour += y_plus_ref * ponderation_k;
                              pond_tot += ponderation_k;
                            }
                        }

                      if (impr_yplus && itisok)
                        {
                          if (y_plus_ref > yplus_ref_max) yplus_ref_max = y_plus_ref;
                          if (y_plus_ref < yplus_ref_min) yplus_ref_min = y_plus_ref;
                          yplus_ref_mean += y_plus_ref;
                          if (u_tau_ref > u_tau_ref_max) u_tau_ref_max = u_tau_ref;
                          if (u_tau_ref < u_tau_ref_min) u_tau_ref_min = u_tau_ref;
                          u_tau_ref_mean += u_tau_ref;
                          yplus_ref_count += 1;
                        }

                    }

                  // On calcul u_tau moyen, y_plus et la vitesse tangente au point force P
                  if (pond_tot > 0.)
                    {
                      mean_u_tau_neighbour /= pond_tot; // u_tau moyen vectoriel
                      mean_y_plus_neighbour /= pond_tot; // y_plus moyen
                      u_tau = 0. ;
                      for(int j = 0; j < nb_comp; j++) u_tau += mean_u_tau_neighbour(0,j) * mean_u_tau_neighbour(0,j);
                      u_tau = sqrt(u_tau); // u_tau = norme de u_tau moyen vectoriel
                      y_plus = u_tau * d1P / nu;

                      if (u_tau > eps) mean_u_tau_neighbour /= u_tau;
                      else mean_u_tau_neighbour = 0. ;// tangente en P
                      if (y_plus > y_c_p_pwl ) // loi polynomiale
                        {
                          if (mean_y_plus_neighbour > y_c_p_pwl ) // coherence loi polynomiale
                            {
                              double vit_coeff = A_pwl * pow (d1P, B_pwl) / pow (nu, B_pwl);
                              for (int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = pow (u_tau, (1+B_pwl)) * vit_coeff * mean_u_tau_neighbour(0,j);
                            }
                          else itisok_P = 0; // non coherence loi polynomiale
                        }
                      else // loi lineaire
                        {
                          if (mean_y_plus_neighbour < y_c_p_pwl ) // coherence loi lineaire
                            {
                              for (int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = (d1P * u_tau * u_tau / nu) * mean_u_tau_neighbour(0,j);
                            }
                          else itisok_P = 0; // non coherence loi lineaire
                        }
                    }
                  // il n y a pas de contribution des dof voisins
                  else itisok_P = 0;
                }
              else
                // il n y a pas de liste des dof voisins
                itisok_P = 0;
            }
          //le point P est sur la frontiere immergee : affectation
          else itisok_P = 0;

          // pour post-pro
          if(itisok_P)
            {
              // pour post-pro & nécessité pour T wall law
              tab_u_star_ibm_(i) = u_tau;
              tab_y_plus_ibm_(i) = y_plus;
            }
          if (impr_yplus && itisok_P)
            {
              if (d1P > d1_max) d1_max = d1P;
              if (d1P < d1_min) d1_min = d1P;
              d1_mean +=  d1P;
              if (y_plus > yplus_max) yplus_max = y_plus;
              if (y_plus < yplus_min) yplus_min = y_plus;
              yplus_mean +=  y_plus;
              if (u_tau > u_tau_max) u_tau_max = u_tau;
              if (u_tau < u_tau_min) u_tau_min = u_tau;
              u_tau_mean += u_tau;
              // Distribution des y+ au voisinage des parois
              if (y_plus > h_yplus_max) h_yplus_max = y_plus;
              if (y_plus < h_yplus_min) h_yplus_min = y_plus;
              h_yplus_mean += y_plus;
              tab_h(0,yplus_count) = y_plus;
              yplus_count += 1;
            }
        }
      else
        {
          for(int j = 0; j < nb_comp; j++) vitesse_imposee_sommet(i,j) = vitesse_imposee_mod(i,j);
        }
    }
  //vitesse_imposee_calculee.echange_espace_virtuel();
  for (int f = 0; f < nb_faces; f++)
    {
      int compteur = 0;
      int som_f = 0;
      while ((som_f < nb_som_face) and (sommet_interp(0, faces_sommets(f,som_f))==1))
        {
          compteur++;
          som_f++;
        }
      if (compteur == nb_som_face)
        {
          for (int som_=0; som_ < nb_som_face; som_++)
            {
              int i = faces_sommets(f,som_);
              for(int j = 0; j < nb_comp; j++)
                {
                  vitesse_imposee_sommet(f,j) += vitesse_imposee_sommet(i,j)/nb_som_face;
                }
            }
        }
    }

  if (impr_yplus && (yplus_count >= 1) && (yplus_ref_count >= 1) )
    {
      d1_mean /= yplus_count;
      yplus_mean /= yplus_count;
      yplus_ref_mean /= yplus_ref_count;
      u_tau_mean /= yplus_count;
      u_tau_ref_mean /= yplus_ref_count;
      Cerr<<"min mean max y  = "<<d1_min<<" "<<d1_mean<<" "<<d1_max<<finl;
      Cerr<<"min mean max y+ = "<<yplus_min<<" "<<yplus_mean<<" "<<yplus_max<<finl;
      Cerr<<"min mean u_tau = "<<u_tau_min<<" "<<u_tau_mean<<" "<<u_tau_max<<finl;
      Cerr<<"min mean max y+_ref = "<<yplus_ref_min<<" "<<yplus_ref_mean<<" "<<yplus_ref_max<<finl;
      Cerr<<"min mean u_tau_ref = "<<u_tau_ref_min<<" "<<u_tau_ref_mean<<" "<<u_tau_ref_max<<finl;

      h_yplus_mean /= yplus_count;
      double range = ( h_yplus_max - h_yplus_min)/N;
      for (int k = 0; k < N+1; k++) abs_h(0,k) = (h_yplus_min + k*range);
      for (int i=0; i < yplus_count; i++)
        {
          for (int j = 0; j < N; j++)
            {
              if ( abs_h(0,j) < tab_h(0,i) &&  tab_h(0,i) <= abs_h(0,j+1) ) compteur_h(0,j) +=1;
            }
        }
      Cerr<<"histogramme y+ compteur =";
      for (int j = 0; j < N; j++) Cerr<<" "<<compteur_h(0,j)<<" ";
      Cerr<<finl;
      Cerr<<"histogramme y+ abscisse =";
      for (int j = 0; j < N+1; j++) Cerr<<" "<<abs_h(0,j)<<" ";
      Cerr<<finl;
      Cerr<<"histogramme y+ : min = "<<h_yplus_min<<" - mean = "<<h_yplus_mean<<" - max =  "<<h_yplus_max<<finl;

      Cerr<<"Moyenne sur "<<yplus_count<<" points forces et "<<yplus_ref_count<<" points references"<<finl;
    }

}


/*##################################################################################################
####################################################################################################
################################# calculer_temperature_imposee #########################################
####################################################################################################
##################################################################################################*/


void Source_PDF_VEF::calculer_temperature_imposee_wall_law()
{
  int dim_esp = Objet_U::dimension;
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const Domaine& dom = domaine_VEF.domaine();
// Variables du domaine
  const IntTab& faces_sommets=domaine_VEF.face_sommets();
  int nb_faces=domaine_VEF.nb_faces();
  int nb_som_face=faces_sommets.dimension(1);
  int nb_som = domaine_VEF.nb_som();
// Tableaux valeurs imposees et calc
  DoubleTab& variable_imposee_mod = modele_lu_.variable_imposee_->valeurs();
  DoubleTab& temperature_imposee_calculee = variable_imposee_;
// Tableaux informations de la paroi immergée
  Interpolation_IBM_thermal_wall_law& interp = ref_cast(Interpolation_IBM_thermal_wall_law,interpolation_lue_.valeur());
  DoubleTab& solid_points = interp.solid_points_->valeurs();
  DoubleTab& fluid_points = interp.fluid_points_->valeurs();
  DoubleTab& fluid_elems = interp.fluid_elems_->valeurs();
// Discrétisation de la grandeur calculée
  Champ_P1NC& champ_variable_inconnue = ref_cast(Champ_P1NC,equation().inconnue());
  DoubleTab& val_variable_inconnue = champ_variable_inconnue.valeurs();
// Récupération des données d'une loi de paroi en vitesse
  const Champ_base& champ_ustar_som = equation().probleme().get_champ("u_star_ibm");
  const DoubleTab& val_ustar = champ_ustar_som.valeurs();
  const Champ_base& champ_yplus = equation().probleme().get_champ("y_plus_ibm");
  const DoubleTab& val_yplus = champ_yplus.valeurs();
// Composante température
  int nb_comp = val_variable_inconnue.dimension(1);
// Variables pour calculs
  DoubleTab temperature_imposee_sommet(nb_som,nb_comp); // Tableau des valeurs aux sommets calc.
  IntTab sommet_interp(nb_som); // Tableau disant si un noeud est forcé
  DoubleTab x(1, dim_esp);
  DoubleTab xf(1, dim_esp);
  DoubleTab Tref(1, nb_comp);
  ArrOfInt cells(1);
  double thetap = 0.0;
  double thetapref = 0.0;
// Choix de modélisation
  int form_Tp = interp.get_formulation_Tp(); // Choix de la loi de paroi utilisée
  int boundary_type = interp.get_boundary_type(); // Choix du type de condition aux limites
  double T_inlet = interp.get_T_inlet();
  double Pr = interp.get_Prandlt_mol();
// Initialisation tableaux
  temperature_imposee_calculee = 0.0;
  temperature_imposee_sommet = 0.0;
  sommet_interp = 0;
// Précision
  double eps = 1e-6;

// Début du calcul de T_imp aux noeuds
  for (int i = 0; i < nb_som; i++) // On parcoure l'ensemble des noeuds du domaine
    {
      if ((fluid_elems(i) >= 0.0)) // Noeuds forcés uniquement
        {
          int itisok = 1;
          sommet_interp(i) = 1; //On note que le sommet i est forcé
          temperature_imposee_sommet(i,0) = variable_imposee_mod(i,0);
          // Obtention du y+ et u* au noeud considéré
          double utau = val_ustar(i);
          double yplus = val_yplus(i);
          // Calcul si utau != 0 => loi classique
          if (utau != 0)
            {
              // Obtention des y et yref, ainsi que la normale à la paroi
              double d1 = 0.0;
              double d2 = 0.0;
              DoubleTab normale(1, nb_comp);
              double norme_de_la_normale= 0.0;
              for(int j = 0; j < nb_comp; j++)
                {
                  double xj = dom.coord(i,j);
                  double xjf = fluid_points(i,j);
                  xf(0, j) = xjf;
                  double xjs = solid_points(i,j);
                  d1 += (xj-xjs)*(xj-xjs);
                  d2 += (xjf-xj)*(xjf-xj);
                  normale(0,j) = xjf - xjs;
                  norme_de_la_normale += normale(0,j)*normale(0,j);
                }
              d1 = sqrt(d1);
              d2 = sqrt(d2);
              double y_ref= d1+d2;
              // traitement des exceptions
              if (d1 < eps) itisok = 0; // le point P est sur la frontiere immergee : affectation
              else
                {
                  norme_de_la_normale = sqrt(norme_de_la_normale);
                  for(int j = 0; j < nb_comp; j++) normale(0,j) /= norme_de_la_normale;
                }
              // Obtention des caractéristiques du fluide utiles au calcul
              double nu = d1 * utau / yplus;
              double rho = equation().milieu().masse_volumique().valeurs()(0,0);
              double Cp = equation().milieu().capacite_calorifique().valeurs()(0, 0);
              //double alpha = equation().milieu().diffusivite().valeurs()(0, 0);
              // -------
              cells[0] = int(fluid_elems(i)); // N° élément fluide
              champ_variable_inconnue.valeur_a_elem(xf, Tref, cells[0]); // Détermination de la température au point fluide de référence
              double yplusref = y_ref * utau / nu;
              if (!form_Tp) //Loi de paroi de Kader
                {
                  thetap = interp.Kader(yplus,Pr); // Calcul du theta plus avec Kader pour le point forcé (connaissance de y+)
                  thetapref = interp.Kader(yplusref,Pr); // Calcul du theta plus avec Kader pour le point fluide
                }
              else
                {
                  Cerr << "Loi non implémentée, choix d'une autre loi obligatoire" << endl;
                  exit();
                }
              if ((itisok == 1) and (!boundary_type))
                {
                  double Tw = variable_imposee_mod(i,0);
                  if (Tw > T_inlet)
                    {
                      temperature_imposee_sommet(i,0) = Tw - abs((thetap / thetapref) * (Tref(0,0) - Tw));
                    }
                  else
                    {
                      temperature_imposee_sommet(i,0) = Tw + abs((thetap / thetapref) * (Tref(0,0) - Tw));
                    }
                }
              else if ((itisok == 1) and (boundary_type == 1))
                {
                  double qw = variable_imposee_mod(i,0);
                  double T_tau = qw/(rho * Cp * utau);
                  temperature_imposee_sommet(i,0) = Tref(0,0) + T_tau * (thetapref - thetap);
                }
              else
                {
                  Cerr << "Erreur dans le choix du type de condition aux limites" << endl;
                  exit();
                }
            }
        }
    }

// Calcul de la température pour les faces forcées
  for (int f = 0; f < nb_faces; f++)
    {
      int compteur = 0;
      int som_f = 0;
      while ((som_f < nb_som_face) and (sommet_interp(faces_sommets(f,som_f))==1))  //On récupère les faces dont l'ensemble des noeuds sont forcées
        {
          compteur++;
          som_f++;
        }
      if (compteur == nb_som_face)
        {
          for (int som_faces =0; som_faces < nb_som_face; som_faces++)
            {
              int i = faces_sommets(f,som_faces);
              temperature_imposee_calculee(f,0) += temperature_imposee_sommet(i,0)/nb_som_face; //Calcul de la température à la face par moyenne des valeurs aux noeuds
            }
        }
    }
}

/*void Source_PDF_VEF::calculer_temperature_imposee_mean_wall_law()
{
//à faire
  int dim_esp = Objet_U::dimension;
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const Domaine& dom = domaine_VEF.domaine();
// Variables du domaine
  const IntTab& faces_sommets=domaine_VEF.face_sommets();
  int nb_faces=domaine_VEF.nb_faces();
  int nb_som_face=faces_sommets.dimension(1);
  int nb_som = domaine_VEF.nb_som();
// Tableaux valeurs imposees et calc
  DoubleTab& variable_imposee_mod = modele_lu_.variable_imposee_->valeurs();
  DoubleTab& temperature_imposee_calculee = variable_imposee_;
// Tableaux informations de la paroi immergée
  Interpolation_IBM_thermal_wall_law& interp = ref_cast(Interpolation_IBM_thermal_wall_law,interpolation_lue_.valeur());
  DoubleTab& solid_points = interp.solid_points_->valeurs();
  DoubleTab& fluid_points = interp.fluid_points_->valeurs();
  DoubleTab& fluid_elems = interp.fluid_elems_->valeurs();
// Discrétisation de la grandeur calculée
  Champ_P1NC& champ_variable_inconnue = ref_cast(Champ_P1NC,equation().inconnue());
  DoubleTab& val_variable_inconnue = champ_variable_inconnue.valeurs();
// Récupération des données d'une loi de paroi en vitesse
  const Champ_base& champ_ustar_som = equation().probleme().get_champ("u_star_ibm");
  const DoubleTab& val_ustar = champ_ustar_som.valeurs();
  const Champ_base& champ_yplus = equation().probleme().get_champ("y_plus_ibm");
  const DoubleTab& val_yplus = champ_yplus.valeurs();
// Composante température
  int nb_comp = val_variable_inconnue.dimension(1);
// Variables pour calculs
  DoubleTab temperature_imposee_sommet(nb_som,nb_comp); // Tableau des valeurs aux sommets calc.
  IntTab sommet_interp(nb_som); // Tableau disant si un noeud est forcé
  DoubleTab x(1, dim_esp);
  DoubleTab xf(1, dim_esp);
  DoubleTab Tref(1, nb_comp);
  ArrOfInt cells(1);
  double thetap = 0.0;
  double thetapref = 0.0;
// Choix de modélisation
  int form_Tp = interp.get_formulation_Tp(); // Choix de la loi de paroi utilisée
  int boundary_type = interp.get_boundary_type(); // Choix du type de condition aux limites
  double T_inlet = interp.get_T_inlet();
  double Pr = interp.get_Prandlt_mol();
// Initialisation tableaux
  temperature_imposee_calculee = 0.0;
  temperature_imposee_sommet = 0.0;
  sommet_interp = 0;
// Précision
  double eps = 1e-6;

// Début du calcul de T_imp aux noeuds
  for (int i = 0; i < nb_som; i++) // On parcoure l'ensemble des noeuds du domaine
    {
      if ((fluid_elems(i) >= 0.0)) // Noeuds forcés uniquement
        {
          int itisok = 1;
          sommet_interp(i) = 1; //On note que le sommet i est forcé
          temperature_imposee_sommet(i,0) = variable_imposee_mod(i,0);
          // Obtention du y+ et u* au noeud considéré
          double utau = val_ustar(i);
          double yplus = val_yplus(i);
          // Calcul si utau != 0 => loi classique
          if (utau != 0)
            {
              // Obtention des y et yref, ainsi que la normale à la paroi
              double d1 = 0.0;
              double d2 = 0.0;
              DoubleTab normale(1, nb_comp);
              double norme_de_la_normale= 0.0;
              for(int j = 0; j < nb_comp; j++)
                {
                  double xj = dom.coord(i,j);
                  double xjf = fluid_points(i,j);
                  xf(0, j) = xjf;
                  double xjs = solid_points(i,j);
                  d1 += (xj-xjs)*(xj-xjs);
                  d2 += (xjf-xj)*(xjf-xj);
                  normale(0,j) = xjf - xjs;
                  norme_de_la_normale += normale(0,j)*normale(0,j);
                }
              d1 = sqrt(d1);
              d2 = sqrt(d2);
              double y_ref= d1+d2;
              // traitement des exceptions
              if (d1 < eps) itisok = 0; // le point P est sur la frontiere immergee : affectation
              else
                {
                  norme_de_la_normale = sqrt(norme_de_la_normale);
                  for(int j = 0; j < nb_comp; j++) normale(0,j) /= norme_de_la_normale;
                }
              // Obtention des caractéristiques du fluide utiles au calcul
              double nu = d1 * utau / yplus;
              double rho = equation().milieu().masse_volumique().valeurs()(0,0);
              double Cp = equation().milieu().capacite_calorifique().valeurs()(0, 0);
              //double alpha = equation().milieu().diffusivite().valeurs()(0, 0);
              double Pr = 1; //nu/alpha;
              // -------
              cells[0] = int(fluid_elems(i)); // N° élément fluide
              champ_variable_inconnue.valeur_a_elem(xf, Tref, cells[0]); // Détermination de la température au point fluide de référence
              double yplusref = y_ref * utau / nu;
              if (!form_Tp) //Loi de paroi de Kader
                {
                  thetap = interp.Kader(yplus,Pr); // Calcul du theta plus avec Kader pour le point forcé (connaissance de y+)
                  thetapref = interp.Kader(yplusref,Pr); // Calcul du theta plus avec Kader pour le point fluide
                }
              else
                {
                  Cerr << "Loi non implémentée, choix d'une autre loi obligatoire" << endl;
                  exit();
                }
              if ((itisok == 1) and (!boundary_type))
                {
                  double Tw = variable_imposee_mod(i,0);
                  if (Tw > T_inlet)
                    {
                      temperature_imposee_sommet(i,0) = Tw - abs((thetap / thetapref) * (Tref(0,0) - Tw));
                    }
                  else
                    {
                      temperature_imposee_sommet(i,0) = Tw + abs((thetap / thetapref) * (Tref(0,0) - Tw));
                    }
                }
              else if ((itisok == 1) and (boundary_type == 1))
                {
                  double qw = variable_imposee_mod(i,0);
                  double T_tau = qw/(rho * Cp * utau);
                  temperature_imposee_sommet(i,0) = Tref(0,0) + T_tau * (thetapref - thetap);
                }
              else
                {
                  Cerr << "Erreur dans le choix du type de condition aux limites" << endl;
                  exit();
                }
            }
        }
    }

// Calcul de la température pour les faces forcées
  for (int f = 0; f < nb_faces; f++)
    {
      int compteur = 0;
      int som_f = 0;
      while ((som_f < nb_som_face) and (sommet_interp(faces_sommets(f,som_f))==1))  //On récupère les faces dont l'ensemble des noeuds sont forcées
        {
          compteur++;
          som_f++;
        }
      if (compteur == nb_som_face)
        {
          for (int som_faces =0; som_faces < nb_som_face; som_faces++)
            {
              int i = faces_sommets(f,som_faces);
              temperature_imposee_calculee(f,0) += temperature_imposee_sommet(i,0)/nb_som_face; //Calcul de la température à la face par moyenne des valeurs aux noeuds
            }
        }
    }
}*/

/*##################################################################################################
####################################################################################################
################################# Pression #########################################################
####################################################################################################
##################################################################################################*/

void Source_PDF_VEF::correct_incr_pressure(const DoubleTab& coeff_node, DoubleTab& correction_en_pression) const
{
  const DoubleTab& aire=champ_aire_->valeurs();
  //  int nb_elem = correction_en_pression.size();
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& elems= domaine_VEF.elem_faces();
  int nb_face_elem=domaine_VEF.domaine().nb_faces_elem();
  int nb_elems=domaine_VEF.domaine().nb_elem_tot();
  int ncomp = coeff_node.dimension(1) ;
  DoubleTrav coef_elem(nb_elems);

  // Filtre par face
  DoubleTrav flag_cl(coeff_node);
  filtre_CLD(flag_cl);
  // int nb_face=domaine_VEF.domaine().nb_faces();
  // for(int face=0; face<nb_face; face++)Cerr<<"face flag_c = "<<face<<" "<<flag_cl(face,0)<<" "<<flag_cl(face,1)<<" "<<flag_cl(face,2)<<finl;

  for(int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      if (aire(num_elem)>0.)
        {
          correction_en_pression(num_elem) = 0.;
        }
      else
        {
          coef_elem(num_elem) = 0. ;
          for (int i=0; i<nb_face_elem; i++)
            {
              int sii=elems(num_elem,i);
              for (int comp=0; comp<ncomp; comp++)
                {
                  if ((coeff_node(sii,comp)!= 1.0) || (flag_cl(sii,comp) != 1.0)) coef_elem(num_elem) += 1. ;
                }
            }
          if (coef_elem(num_elem) == (ncomp*nb_face_elem)) correction_en_pression(num_elem) = 0.;
        }
    }
}

void Source_PDF_VEF::correct_pressure(const DoubleTab& coeff_node, DoubleTab& pression, const DoubleTab& correction_en_pression) const
{
  const DoubleTab& aire=champ_aire_->valeurs();
  //int nb = pression.size();
  //Cerr << pression << finl;
  //Cerr << nb << finl;
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& elems= domaine_VEF.elem_faces();
  int nb_face_elem=domaine_VEF.domaine().nb_faces_elem();
  int nb_elem= domaine_VEF.domaine().nb_elem_tot();
  int ncomp = coeff_node.dimension(1) ;
  DoubleTrav coef_elem(nb_elem);

  DoubleTrav flag_cl(coeff_node);
  filtre_CLD(flag_cl);

  for(int num_elem=0; num_elem<nb_elem; num_elem++)
    {
      if (aire(num_elem)<=0.)
        {
          coef_elem(num_elem) = 0. ;
          for (int i=0; i<nb_face_elem; i++)
            {
              int sii=elems(num_elem,i);
              for (int comp=0; comp<ncomp; comp++)
                {
                  if ((coeff_node(sii,comp)!= 1.0) || (flag_cl(sii,comp) != 1.0)) coef_elem(num_elem) += 1. ;
                }
            }
          if (coef_elem(num_elem) != (ncomp*nb_face_elem)) pression(num_elem)+=correction_en_pression(num_elem);
        }
    }
}

void Source_PDF_VEF::filtre_CLD(DoubleTab& flag_cl) const
{
  // Filtre et dof par face
  const Domaine_Cl_VEF& le_dom_cl = le_dom_Cl_VEF.valeur();
  int nb_cond_lim = le_dom_cl.nb_cond_lim();
  int nb_comp = flag_cl.dimension(1);
  flag_cl = 1.;
  for (int ij=0; ij<nb_cond_lim; ij++)
    {
      const Cond_lim_base& la_cl_base = le_dom_cl.les_conditions_limites(ij).valeur();
      const Front_VF& le_bord = ref_cast(Front_VF,le_dom_cl.les_conditions_limites(ij)->frontiere_dis());
      int nfaces = le_bord.nb_faces_tot();
      if (sub_type(Dirichlet,la_cl_base))
        {
          for (int ind_face=0; ind_face < nfaces; ind_face++)
            {
              int face=le_bord.num_face(ind_face);
              for (int comp=0; comp<nb_comp; comp++) flag_cl(face,comp) = 0. ;
            }
        }
    }
  //vector; NS Attention seulement defini en EF
  // if (nb_comp == Objet_U::dimension) le_dom_cl.imposer_symetrie(flag_cl,0);
  flag_cl.echange_espace_virtuel();
  return;
}

/*##################################################################################################
####################################################################################################
################################# Impression #######################################################
####################################################################################################
##################################################################################################*/

int Source_PDF_VEF::impr(Sortie& os) const
{
  if (out_=="??")
    {
      Cerr << __FILE__ << (int)__LINE__ << " ERROR / Source_PDF_VEF::impr" << finl;
      Cerr << "No balance printed for " << que_suis_je() << finl;
      Cerr << "Because output file name is not specified." << finl;
      return 0;
    }
  else
    {
      int nb_compo=bilan_.size();
      if (nb_compo==0)
        {
          Cerr << __FILE__ << (int)__LINE__ << " ERROR / Source_PDF_VEF::impr" << finl;
          Cerr << "No balance printed for " << que_suis_je() << finl;
          Cerr << "Because bilan_ array is not filled." << finl;
          return 0;
        }
      else
        {
          const Schema_Temps_base& sch=equation().probleme().schema_temps();
          double temps=sch.temps_courant();
          double pdtps = sch.pas_de_temps();
          if (temps == pdtps)
            {
              int flag = Process::je_suis_maitre();
              if(flag) ouvrir_fichier(Flux,"",flag);
              return 0;
            }
          const DoubleTab& variable=equation().inconnue().valeurs();
          Nom espace=" \t";

          int nb_comp = variable.dimension(1);
          assert(nb_comp == modele_lu_.dim_variable_);
          int pdf_dt_conv = modele_lu_.pdf_bilan();

          DoubleTrav secmem_conv(variable);
          // defaut pdf_dt_conv = 0
          secmem_conv = 0.;
          int i_traitement_special = 0;

          if (pdf_dt_conv == 1 || pdf_dt_conv == 2)
            // pdf_dt_conv = 1 ou 2
            {
              if (nb_comp == Objet_U::dimension)
                {
                  //vector; NS
                  Navier_Stokes_std& eq_NS = ref_cast_non_const(Navier_Stokes_std, equation());
                  int transport_rhou=0;
                  if (eq_NS.nombre_d_operateurs() > 1)
                    {
                      transport_rhou= (eq_NS.vitesse_pour_transport().le_nom()=="rho_u"?1:0);
                      if (pdf_dt_conv == 2) eq_NS.derivee_en_temps_conv(secmem_conv , variable);
                    }
                  i_traitement_special = (transport_rhou==1?2:102);
                }
              else if (nb_comp == 1)
                {
                  //scalar; temps en rho
                  Equation_base& eq_Ba = ref_cast_non_const(Equation_base, equation());
                  if (equation().nombre_d_operateurs() > 1)
                    {
                      if (pdf_dt_conv == 2) eq_Ba.derivee_en_temps_conv(secmem_conv , variable);
                    }
                  i_traitement_special = 2;
                }
              else
                {
                  Cerr << "Source_PDF_VEF: for scalar or vector only; dim = " <<nb_comp<< finl;
                  Process::exit();
                }
            }
          else if(pdf_dt_conv != 0 )
            {
              Cerr<<"Source_PDF_VEF: Modele pdf_bilan must be 0; 1 or 2 only"<<endl;
              Process::exit();
            }

          Source_PDF_base& my_src = ref_cast_non_const(Source_PDF_base,*this);
          my_src.compute_source_term_PDF(i_traitement_special, secmem_conv, bilan_);

          int flag = Process::je_suis_maitre();
          if(flag)
            {
              ouvrir_fichier(Flux,"",flag);
              if(nb_comp == Objet_U::dimension)
                {
                  // vector
                  assert(Objet_U::dimension==3);
                  Flux << temps << espace << bilan_(0) << espace << bilan_(1) << espace << bilan_(2) << finl;
                }
              else
                {
                  // scalar
                  Flux << temps << espace << bilan_(0) << finl;
                }
            }
        }
      return 1;
    }
}
//  void imprimer_ustar_yplus__mean_only(Sortie&, const Nom& ) const override;


/*##################################################################################################
####################################################################################################
################################# POSTRAITEMENT ####################################################
####################################################################################################
##################################################################################################*/

void Source_PDF_VEF::creer_champ(const Motcle& motlu)
{
  Source_PDF_base::creer_champ(motlu);

  if (motlu=="u_star_ibm" && !champ_u_star_ibm_.non_nul() && imm_wall_law_)
    {
      int nb_comp = 1;
      Noms noms(1);
      noms[0]="u_star_ibm";
      Noms unites(1);
      unites[0] = "m/s";
      double temps=0.;
      const Discretisation_base& discr = equation().probleme().discretisation();
      discr.discretiser_champ("champ_sommets",equation().domaine_dis(),scalaire,noms,unites,nb_comp,temps,champ_u_star_ibm_);
      champs_compris_.ajoute_champ(champ_u_star_ibm_);
    }
  else if (motlu=="y_plus_ibm" && !champ_y_plus_ibm_.non_nul() && imm_wall_law_)
    {
      int nb_comp = 1;
      Noms noms(1);
      noms[0]="y_plus_ibm";
      Noms unites(1);
      unites[0] = "-";
      double temps=0.;
      const Discretisation_base& discr = equation().probleme().discretisation();
      discr.discretiser_champ("champ_sommets",equation().domaine_dis(),scalaire,noms,unites,nb_comp,temps,champ_y_plus_ibm_);
      champs_compris_.ajoute_champ(champ_y_plus_ibm_);
    }
}

bool Source_PDF_VEF::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  if (Source_PDF_base::has_champ(nom)) return Source_PDF_base::has_champ(nom, ref_champ);

  if (nom == "u_star_ibm" && champ_u_star_ibm_.non_nul())
    {
      ref_champ = Source_PDF_VEF::get_champ(nom);
      return true;
    }
  else if (nom == "y_plus_ibm" && champ_y_plus_ibm_.non_nul())
    {
      ref_champ = Source_PDF_VEF::get_champ(nom);
      return true;
    }
  else
    return false; /* rien trouve */
}

bool Source_PDF_VEF::has_champ(const Motcle& nom) const
{
  if (Source_PDF_base::has_champ(nom)) return true;

  if (nom == "u_star_ibm" && champ_u_star_ibm_.non_nul())
    return true;
  else if (nom == "y_plus_ibm" && champ_y_plus_ibm_.non_nul())
    return true;
  else
    return false; /* rien trouve */
}

const Champ_base& Source_PDF_VEF::get_champ(const Motcle& nom) const
{
  if (Source_PDF_base::has_champ(nom)) return Source_PDF_base::get_champ(nom);
  if (nom=="u_star_ibm")
    {
      if (champ_u_star_ibm_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));
      // Initialisation a 0 du champ volumique u_star
      DoubleTab& valeurs = champ_u_star_ibm_->valeurs();
      valeurs=0;
      if (tab_u_star_ibm_.size_array()>0)
        {
          const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
          int nb_som=domaine_VEF.nb_som();
          for (int num_node=0; num_node<nb_som; num_node++)
            valeurs(num_node)=tab_u_star_ibm_(num_node);
        }
      valeurs.echange_espace_virtuel();
      champ_u_star_ibm_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else if (nom=="y_plus_ibm")
    {
      if (champ_y_plus_ibm_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ volumique u_star
      DoubleTab& valeurs = champ_y_plus_ibm_->valeurs();
      valeurs=0;
      if (tab_y_plus_ibm_.size_array()>0)
        {
          const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
          int nb_som=domaine_VEF.nb_som();
          for (int num_node=0; num_node<nb_som; num_node++)
            valeurs(num_node)=tab_y_plus_ibm_(num_node);
        }
      valeurs.echange_espace_virtuel();
      champ_y_plus_ibm_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else
    return champs_compris_.get_champ(nom);
}

void Source_PDF_VEF::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Source_PDF_base::get_noms_champs_postraitables(nom,opt);

  Noms noms_compris = champs_compris_.liste_noms_compris();
  if(imm_wall_law_)
    {
      noms_compris.add("u_star_ibm");
      noms_compris.add("y_plus_ibm");
    }
  if (opt==DESCRIPTION)
    Cerr<<" Source_PDF_VEF : "<< noms_compris <<finl;
  else
    nom.add(noms_compris);
}
