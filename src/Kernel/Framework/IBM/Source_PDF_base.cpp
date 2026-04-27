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

#include <Domaine_Cl_dis_base.h>
#include <Schema_Temps_base.h>
#include <Source_PDF_base.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Type_info.h>
#include <SFichier.h>
#include <Param.h>
#include <Interpolation_IBM_elem_fluid.h>
#include <Interpolation_IBM_mean_gradient.h>
#include <Interpolation_IBM_hybrid.h>
#include <Interpolation_IBM_power_law_tbl.h>
#include <Interpolation_IBM_power_law_tbl_u_star.h>
#include <Interpolation_IBM_thermal_wall_law.h>

Implemente_base(Source_PDF_base,"Source_PDF_base",Source_dep_inco_base);
// XD source_pdf_base Source_dep_inco_base source_pdf_base 1 Basic class of source_PDF terms introduced in the equation.

Entree& Source_PDF_base::readOn(Entree& s)
{
  Param param(que_suis_je());
  param.ajouter("prepro_ibm", &prepro_lu_,Param::OPTIONAL); // to realise the PDF IBM preprocessing
  param.ajouter("aire", &champ_aire_lu_,Param::OPTIONAL); // XD_ADD_P field_base volumic field: a boolean for the cell (0 or 1) indicating if the obstacle is in the cell
  param.ajouter_flag("get_aire_from_prepro", &aire_from_prepro_); // XD_ADD_P get aire IBM from prepro.
  param.ajouter("barycentre", &champ_barycentre_lu_,Param::OPTIONAL); // XD_ADD_P field_base volumic field with 3 components representing the face barycenters
  param.ajouter_flag("get_barycenter_from_prepro", &barycentre_from_prepro_); // XD_ADD_P get aire IBM from prepro.
  param.ajouter("rotation", &champ_rotation_lu_,Param::OPTIONAL); // XD_ADD_P field_base volumic field with 9 components representing the change of basis on cells (local to global). Used for rotating cases for example.
  param.ajouter_flag("get_rotation_from_prepro", &rotation_from_prepro_); // XD_ADD_P get aire IBM from prepro.
  param.ajouter_flag("transpose_rotation", &transpose_rotation_); // XD_ADD_P rien  whether to transpose the basis change matrix.
  param.ajouter("modele",&modele_lu_,Param::REQUIRED);   // XD_ADD_P bloc_pdf_model model used for the Penalized Direct Forcing
  temps_relax_ = modele_lu_.temps_relax_;
  echelle_relax_ =  modele_lu_.echelle_relax_;
  param.ajouter("interpolation",&interpolation_lue_,Param::OPTIONAL); // XD_ADD_P interpolation_ibm_base interpolation method

  param.lire_avec_accolades(s);
  if (interpolation_lue_)
    {
      if (!(interpolation_lue_->que_suis_je() == "Interpolation_IBM_aucune"))
        {
          interpolation_bool_ = true;
          if ((interpolation_lue_->que_suis_je() == "Interpolation_IBM_power_law_tbl") || (interpolation_lue_->que_suis_je() == "Interpolation_IBM_power_law_tbl_u_star") ) imm_wall_law_ = true;
          interpolation_lue_->set_source(*this);
        }
    }

  if (!equation().probleme().que_suis_je().contient("_IBM"))
    {
      Cerr << "Error in the equation " << equation().que_suis_je() << " !" << finl;
      Cerr << "You are using a " << que_suis_je() << " source term but your problem is not an IBM problem !" << finl;
      Cerr << "Fix your data set and use one of the following problems : " << finl;

      Noms noms_pb;
      Type_info::les_sous_types(Nom("Probleme_base"), noms_pb);

      for (auto &itr : noms_pb)
        if (itr.contient("_IBM"))
          Cerr << "  - " << itr << finl;

      Cerr << "  - Pb_Couple_Optimisation_IBM" << finl;

      Process::exit();
    }

  return s;
}

Sortie& Source_PDF_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

void Source_PDF_base::associer_pb(const Probleme_base& pb)
{
  if (prepro_lu_.non_nul())
    {
      prepro_lu_->associer_pb(pb);

      // Verification eventuelle
      if( prepro_lu_->verify_results_prepro() == 1) verify_results_prepro();
    }

  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();

  // Matrice Rotation
  int dim_esp = Objet_U::dimension;
  assert(dim_esp==3);
  int nb_comp=dim_esp*dim_esp;
  Noms nom_c(nb_comp);
  Noms unites(nb_comp);
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nom_c,unites,nb_comp,0.,champ_rotation_);
  if (prepro_lu_.non_nul() && rotation_from_prepro_)
    {
      champ_rotation_->valeurs() = prepro_lu_->get_champ_rotation();
    }
  else
    {
      if (champ_rotation_lu_.non_nul())
        {
          champ_rotation_->affecter(champ_rotation_lu_);
        }
      else
        {
          Cerr<<"Source_PDF_base::associer_pb: rotation term is missing "<<finl;
          exit();
        }
    }

  // Barycentre
  int nb_comp0=dim_esp;
  Noms nom_c0(nb_comp0);
  Noms unites0(nb_comp0);
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nom_c0,unites0,nb_comp0,0.,champ_barycentre_);
  if (prepro_lu_.non_nul() && barycentre_from_prepro_)
    {
      champ_barycentre_->valeurs() = prepro_lu_->get_champ_barycentre();
    }
  else
    {
      if (champ_barycentre_lu_.non_nul())
        {
          champ_barycentre_->affecter(champ_barycentre_lu_);
        }
      // else
      // barycentres non requis sans prepro_ibm
      //   {
      //     Cerr<<"Source_PDF_base::associer_pb: barycenter term is missing "<<finl;
      //     exit();
      //   }
    }

  // Champ pour terme source PDF
  if (equation().discretisation().is_ef())
    pb.discretisation().discretiser_champ("champ_sommets",le_dom_dis,"","",1,0., champ_nodal_);
  else
    pb.discretisation().discretiser_champ("champ_face",le_dom_dis,"","",1,0., champ_nodal_);

  // Aire
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis ,"aire","m-1",1,0., champ_aire_);
  if (prepro_lu_.non_nul() && aire_from_prepro_)
    {
      champ_aire_->valeurs() = prepro_lu_->get_champ_aire();
    }
  else
    {
      if (champ_aire_lu_.non_nul())
        {
          champ_aire_->affecter(champ_aire_lu_);
        }
      else
        {
          Cerr<<"Source_PDF_base::associer_pb: aire term is missing "<<finl;
          exit();
        }
    }

  // Vitesse Imposee Shape IB
  if (get_modele().get_PDF_mobile())
    {
      nb_comp=dim_esp;
      Noms nom_c11(nb_comp);
      Noms unites11(nb_comp);
      pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nom_c11,unites11,nb_comp,0.,modele_lu_.vitesse_shape_IBM_);
    }

  // Variable Imposee sur IB
  const DoubleTab& variable=equation().inconnue().valeurs();
  Motcle directive("temperature");
  if (equation().inconnue().nature_du_champ()==vectoriel) directive="vitesse";
  nb_comp=variable.dimension(1);
  assert(nb_comp==modele_lu_.dim_variable_);
  Nom nom_c1(equation().inconnue().le_nom());
  Nom unites1(equation().inconnue().unites()[0]);
  type_variable_imposee_ = modele_lu_.type_variable_imposee_;
  pb.discretisation().discretiser_champ(directive,le_dom_dis ,nom_c1,unites1,nb_comp,0.,modele_lu_.variable_imposee_);

  if (interpolation_bool_)
    {
      Domaine_dis_base& le_dom_dis_base = ref_cast_non_const(Domaine_dis_base,le_dom_dis);
      interpolation_lue_->discretise(pb.discretisation(),le_dom_dis_base);
    }

  // set variable imposee
  set_variable_imposee();

  // Rho
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis ,"rho","kg.m-3",1,0., champ_rho_);
  Source_PDF_base::updateChampRho();

  // Relaxation
  temps_relax_ = modele_lu_.temps_relax_;
  echelle_relax_ = modele_lu_.echelle_relax_;

  matrice_pression_variable_bool_ = false;
  if (temps_relax_ != 1.0e+12) matrice_pression_variable_bool_ = true;
}

void Source_PDF_base::set_variable_imposee()
{
  int dim_esp = Objet_U::dimension;
  const Probleme_base& pb = equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();

  // Transposition Matrice Rotation
  if (transpose_rotation_)
    {
      DoubleTab& val=champ_rotation_->valeurs();

      int nb_case=val.dimension_tot(0);
      for (int ele=0; ele<nb_case; ele++)
        for (int k=0; k<dim_esp; k++)
          for (int i=k+1; i<dim_esp; i++)
            {
              double tmp=val(ele,3*k+i);
              val(ele,3*k+i)=val(ele,3*i+k);
              val(ele,3*i+k)=tmp;
            }
    }

  // Variable Imposee sur IB
  const DoubleTab& variable=equation().inconnue().valeurs();
  int nb_comp=variable.dimension(1);

  if (interpolation_bool_) // interpolation
    {
      if (type_variable_imposee_ == 1) // fonction xyzt
        {
          if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_gradient_moyen" || interpolation_lue_->que_suis_je() == "Interpolation_IBM_power_law_tbl_u_star")
            {
              const Interpolation_IBM_mean_gradient& interp = ref_cast(Interpolation_IBM_mean_gradient,interpolation_lue_.valeur());
              this->compute_variable_imposee_projete(interp.solid_elems_->valeurs(), interp.solid_points_->valeurs(), -2.0, 1e-6);
            }
          else
            {
              const Interpolation_IBM_elem_fluid& interp = ref_cast(Interpolation_IBM_elem_fluid,interpolation_lue_.valeur());
              compute_variable_imposee_projete(interp.fluid_elems_->valeurs(), interp.solid_points_->valeurs(), -2.0, 1e-6);
            }
        }
      else // data
        {
          modele_lu_.variable_imposee_->affecter(modele_lu_.variable_imposee_lu_);
        }
    }
  else // pas d'interpolation
    {
      if (type_variable_imposee_ != 1) // data
        {
          modele_lu_.variable_imposee_->affecter(modele_lu_.variable_imposee_lu_);
        }
      else // fonction xyzt
        {
          const DoubleTab& coords = le_dom_dis.domaine().coord_sommets();
          Domaine_VF& le_dom_VF = ref_cast_non_const(Domaine_VF, pb.domaine_dis());
          modele_lu_.affecter_variable_imposee(le_dom_VF, coords);
        }
    }
  if (modele_lu_.local_ == 1)
    {
      if (nb_comp != dim_esp)
        {
          Cerr<<"Source_PDF_base::calcul_variable_imposee: dimension variable differente "<<dim_esp<<"! "<<finl;
          abort();
        }
      rotate_imposed_velocity(modele_lu_.variable_imposee_->valeurs());
    }

  variable_imposee_ = modele_lu_.variable_imposee_->valeurs();

  compute_indicateur_nodal_champ_aire();
}

void Source_PDF_base::set_fields_to_prepro(Prepro_IBM_base& un_prepro)
{
  un_prepro.set_champ_aire(champ_aire_->valeurs());
  un_prepro.set_champ_rotation(champ_rotation_->valeurs());
  un_prepro.set_champ_barycentre(champ_barycentre_->valeurs());
}

void Source_PDF_base::get_fields_from_prepro(Prepro_IBM_base& un_prepro)
{
  champ_aire_->valeurs() = un_prepro.get_champ_aire();
  champ_rotation_->valeurs() = un_prepro.get_champ_rotation();
  champ_barycentre_->valeurs() = un_prepro.get_champ_barycentre();
}

void Source_PDF_base::compute_NeighNode_IBM_elem(DoubleTab& aireArray, IntLists& elem_voisins, bool all_elem_vois)
{
  const Probleme_base& pb = equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine_VF& the_dom_VF = ref_cast(Domaine_VF,le_dom_dis);
  const Domaine& le_dom =  le_dom_dis.domaine();
  int nb_elem_tot = le_dom.nb_elem_tot();
  int nb_som_elem = le_dom.nb_som_elem();
  int nb_faces_elem = le_dom.nb_faces_elem();
  const IntTab& face_voisins = the_dom_VF.face_voisins();
  const IntTab& elems = le_dom.les_elems() ;
  const IntTab& elem_face = the_dom_VF.elem_faces();

  // On vide la Lists
  assert (elem_voisins.size() == nb_elem_tot);
  for (int num_elem = 0; num_elem < nb_elem_tot; num_elem++) elem_voisins[num_elem].vide();

  for (int num_elem = 0; num_elem < nb_elem_tot; num_elem++) // inclus les elem ghost
    {
      if (aireArray(num_elem) > 0.)
        {
          if (!elem_voisins[num_elem].est_vide()) // Liste non vide !
            {
              Cerr<<"IBM : Source_PDF_base::compute_NeighNode_IBM_elem:  Neighbour liste is not empty for element "<<num_elem<<" : Abort"<<finl;
              exit();
            }
          // Cerr << "******Elem = "<<num_elem<<finl;

          IntList node_of_elem;
          node_of_elem.vide();
          for (int nod=0; nod<nb_som_elem; nod++)
            {
              int nod_glob = elems(num_elem, nod);
              node_of_elem.add_if_not(nod_glob);
            }
          for (int fac=0; fac<nb_faces_elem; fac++)
            {
              int num_fac = elem_face(num_elem, fac);
              // Cerr << "face "<<fac<<" : num_face = "<<num_fac<<finl;
              for (int voisin=0; voisin<2; voisin++) // Deux voisins seulement
                {
                  int num_elem_v = face_voisins(num_fac,voisin);
                  if ( (num_elem_v!=-1) && (num_elem_v!=num_elem) )
                    {
                      if (aireArray(num_elem_v) > 0. || all_elem_vois) elem_voisins[num_elem].add_if_not(num_elem_v); // voisin de num_elem par une face
                      // Cerr << "voisin "<<voisin<<" : num_elem_v = "<<num_elem_v<<finl;
                      // Cerr<<"list elem_voisins => ";
                      // for (int v=0; v<elem_voisins[num_elem].size(); v++) Cerr<<(elem_voisins[num_elem])[v]<< " ";
                      // Cerr<<finl;

                      // element num_elem_v_v, voisin de num_elem_v par une face
                      for (int fac_v=0; fac_v<nb_faces_elem; fac_v++)
                        {
                          int num_fac_v = elem_face(num_elem_v, fac_v);
                          for (int voisin_v=0; voisin_v<2; voisin_v++)
                            {
                              int num_elem_v_v = face_voisins(num_fac_v,voisin_v);
                              if ( (num_elem_v_v!=-1) && (num_elem_v_v!=num_elem_v) && (num_elem_v_v!=num_elem) )
                                {
                                  if (aireArray(num_elem_v_v) > 0. || all_elem_vois)
                                    {
                                      int iok = 0 ;
                                      for (int nod_v_v=0; nod_v_v<nb_som_elem; nod_v_v++)
                                        {
                                          int num_nod_v_v = elems(num_elem_v_v, nod_v_v);
                                          if (node_of_elem.contient(num_nod_v_v)) iok = 1;
                                        }
                                      // Cerr << "num_elem_v_v = "<<num_elem_v_v<<" iok = "<<iok<<finl;
                                      if (iok == 1) elem_voisins[num_elem].add_if_not(num_elem_v_v); // voisin de num_elem par un noeud
                                    }
                                }
                            }
                        }
                    }
                }
            }
          // Cerr<<"=> "<<elem_voisins[num_elem].size()<<" voisins de "<<num_elem<<" : ";
          // for (int voisin=0; voisin<elem_voisins[num_elem].size(); voisin++) Cerr<<(elem_voisins[num_elem])[voisin]<< " ";
          // Cerr<<finl;
        }
    }
  return;
}

double Source_PDF_base::aire_geometrique_IBM(DoubleTab& rotation, int elem)
{
  int dim_esp = Objet_U::dimension;
  const Probleme_base& pb = equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine_VF& the_dom_VF = ref_cast(Domaine_VF,le_dom_dis);
  const Domaine& dom = le_dom_dis.domaine();
  const IntTab& elems = dom.les_elems() ;
  const DoubleTab coordsDom3D=dom.les_sommets();
  int nb_som_elem = dom.nb_som_elem();

  DoubleTrav hmax(dim_esp);
  for (int k=0; k < dim_esp; k++)
    {
      double xmin = 1.e30;
      double xmax = -1.e30;
      for (int id_loc=0; id_loc < nb_som_elem; id_loc++)
        {
          int id_glob = elems(elem,id_loc);
          xmin = std::min(xmin, coordsDom3D(id_glob, k));
          xmax = std::max(xmax, coordsDom3D(id_glob, k));
        }
      hmax(k) = xmax-xmin;
    }
  double hcharac = 0.;
  for (int k=0; k<dim_esp; k++) hcharac += hmax(k) * abs(rotation(elem,3*k+(dim_esp-1)));
  hcharac = abs(hcharac);
  // hcharac = (prepro_lu_->get_h_max_elem())(e);
  double aire = the_dom_VF.volumes(elem) / hcharac;

  // Cerr<< "=> element: "<<elem<<" ; hmax "<< hmax<<finl;
  // Cerr<< "normal "<<rotation(elem,3*0+(dim_esp-1))<<" "<<rotation(elem,3*1+(dim_esp-1))<<" "<<rotation(elem,3*2+(dim_esp-1))<<finl;
  // Cerr<< "hcharac "<< hcharac<<" ; aire = "<<aire<<finl;

  return aire;
}

void Source_PDF_base::update_pseudo_level_set_IBM(DoubleTab& vecteur_deplacement, double alpha)
{
  // Update a partir de la pseudo level set et vecteur deplacement vertex

  if (!(getInterpolationBool() == true))
    {
      Cerr << "Source_PDF_base::update_pseudo_level_set_IBM: No Interpolation. Aborting..." << finl;
      abort();
    }
  if (!(prepro_lu_.non_nul()))
    {
      Cerr << "Source_PDF_base::update_pseudo_level_set_IBM: No IBM prepro. Aborting..." << finl;
      abort();
    }

  int dim_esp = Objet_U::dimension;
  assert(dim_esp == 3);
  Interpolation_IBM_base& my_interp = ref_cast_non_const(Interpolation_IBM_base, getInterpolationLu());
  DoubleTab& level_set = ref_cast_non_const(DoubleTab, my_interp.get_pseudo_level_set());
  int nb_som = level_set.dimension(0);

  DoubleTab& bary = champ_barycentre_->valeurs();
  DoubleTab& rotation = champ_rotation_->valeurs();
  DoubleTab& aire = champ_aire_->valeurs();
  int nb_elem = aire.dimension(0);

  const DoubleTab& nor = my_interp.get_normal_proj_solid();
  assert (nor.dimension(0) == nb_som);
  assert (nor.dimension(1) == dim_esp);

  assert (vecteur_deplacement.dimension(0) == nb_som);
  assert (vecteur_deplacement.dimension(1) == dim_esp);

  const Probleme_base& pb = equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine& dom = le_dom_dis.domaine();
  const IntTab& elems= dom.les_elems() ;
  int nb_som_elem=dom.nb_som_elem();
  const DoubleTab coordsDom3D=dom.les_sommets();

  // Preparation de la level-set
  double nor_nor = 0.;
  for (int e=0; e<nb_elem; e++)
    {
      int nb_nor_nul = 0;
      for (int il=0; il<nb_som_elem; il++)
        {
          int i = elems(e,il);
          if (level_set(i) == -99999.) // level_set non initialise
            {
              nb_nor_nul += 1;
              continue;
            }
          nor_nor = 0.;
          for (int k=0; k<dim_esp; k++) nor_nor += nor(i, k)*nor(i, k);
          if (nor_nor <= 1.e-10) nb_nor_nul += 1; // vertex non calcule PDF
        }
      if ( nb_nor_nul == nb_som_elem )
        for (int il=0; il<nb_som_elem; il++) level_set(elems(e,il)) = -99999.;
    }

  // Gradient de phi par element (schema explicite) + sommet (robustesse)
  DoubleTrav grad_phi_e_0(nb_elem,dim_esp) ;
  DoubleTrav grad_phi_i_0(nb_som,dim_esp) ;
  IntTrav contrib_i(nb_som);
  int etarg = -1; // pour debug
  // grad_phi par element (normee)
  for (int e=0; e<nb_elem; e++)
    {
      int initialised_vertex_in_e = 1;
      if (e == etarg) Cerr<<">>> elem = "<<e<<" : grad phi"<<finl;
      // Ajout contributions dans l element
      for (int il=0; il<nb_som_elem; il++)
        {
          int i = elems(e,il);
          if (level_set(i) == -99999.) // level_set non initialisee => pas de contrib.
            {
              initialised_vertex_in_e = 0;
              continue;
            }
          int dir_grad_phi_0 = (level_set(i)<0.?-1:(level_set(i)==0.?0.:1)); // grad_phi suivant phi croissant (dir_grad_phi * nor)
          nor_nor = 0.;
          for (int k=0; k<dim_esp; k++) nor_nor += nor(i, k)*nor(i, k);
          if (e == etarg) Cerr<<"i = "<<i<<" norme carree  normal = "<<nor_nor<<" dir_grad_phi_0 = "<<dir_grad_phi_0<<finl;
          if (nor_nor > 1.e-10)
            {
              for (int k=0; k<dim_esp; k++) grad_phi_e_0(e,k) += dir_grad_phi_0*nor(i, k);
            }
        }
      // normalisation grad_phi_e_0 par element si non null
      // rem. : si level_set non initialisee pour tout noeud de l element => grad_phi_e_0 null
      // Traitement initialised_vertex_in_e = 0 et Assemblage sommet
      nor_nor = 0.;
      for (int k=0; k<dim_esp; k++) nor_nor += grad_phi_e_0(e,k)*grad_phi_e_0(e,k);
      if (nor_nor > 1.e-10)
        {
          for (int k=0; k<dim_esp; k++) grad_phi_e_0(e,k) /= sqrt(nor_nor);

          // Si initialised_vertex_in_e = 0
          //   => initialisation de ces vertex a partir de grad_phi_e_0 et des autres vertex
          if ( initialised_vertex_in_e == 0)
            {
              // boucle recherche vertex de reference
              double level_set_ref = -99999.;
              DoubleTrav XYZ_ref(dim_esp);
              for (int il=0; il<nb_som_elem; il++)
                {
                  int i = elems(e,il);
                  nor_nor = 0.;
                  for (int k=0; k<dim_esp; k++) nor_nor += nor(i, k)*nor(i, k);
                  if (level_set(i) != -99999. && nor_nor > 1.e-10)
                    {
                      level_set_ref = level_set(i) ;
                      for (int k=0; k<dim_esp; k++) XYZ_ref(k) = coordsDom3D(i,k);
                      break; // on a trouve
                    }
                }
              if (e == etarg) Cerr<<"Vertex non initialise : level_set_ref = "<<level_set_ref<<finl;
              // boucle definition level_set vertex non initialises
              for (int il=0; il<nb_som_elem; il++)
                {
                  int i = elems(e,il);
                  if (level_set(i) == -99999.) // on initialise
                    {
                      double prod_sca = 0.;
                      for (int k=0; k<dim_esp; k++) prod_sca += (coordsDom3D(i,k)-XYZ_ref(k))*grad_phi_e_0(e,k);
                      level_set(i) = level_set_ref + prod_sca;
                    }
                }
            }

          // Assemblage sommet
          for (int il=0; il<nb_som_elem; il++)
            {
              int i = elems(e,il);
              for (int k=0; k<dim_esp; k++) grad_phi_i_0(i,k) += grad_phi_e_0(e,k);
              contrib_i(i) += 1;
            }
        }

      if (e == etarg)
        {
          for (int il=0; il<nb_som_elem; il++)
            {
              int i = elems(e,il);
              Cerr<<"i = "<<i<<" ; old level set = "<<level_set(i)<<" ; nor = ";
              for (int k=0; k<dim_esp; k++) Cerr<<nor(i,k)<<" ";
              Cerr<<" ; deplacement = ";
              for (int k=0; k<dim_esp; k++) Cerr<<vecteur_deplacement(i,k)<<" ";
              Cerr<<finl;
            }
          Cerr<<"grad_phi_e_0 = ";
          for (int k=0; k<dim_esp; k++) Cerr<<grad_phi_e_0(e,k)<<" ";
          Cerr<<finl;
        }
    }

  // grad_phi par sommet (norme)
  for (int i=0; i<nb_som; i++)
    {
      if (contrib_i(i) != 0)
        {
          for (int k=0; k<dim_esp; k++) grad_phi_i_0(i,k) /=contrib_i(i);

          // normalisation
          nor_nor = 0.;
          for (int k=0; k<dim_esp; k++) nor_nor += grad_phi_i_0(i,k)*grad_phi_i_0(i,k);
          if (nor_nor > 1.e-10)
            for (int k=0; k<dim_esp; k++) grad_phi_i_0(i,k) /= sqrt(nor_nor);
        }
    }

  // Advecte phi par sommet pour tous les sommets + new_bary pour element traverse par level set zero (moyenne arith)
  // normale egale a grad_phi_e_0 pour ces elements
  IntTrav advect_faite(nb_som) ;
  aire = 0.; // Mise a zero champ aire IB
  bary = 0.; // Mise a zero champ barycentre IB
  DoubleTrav new_normal_e(bary);
  for (int e=0; e<nb_elem; e++)
    {
      if (e == etarg) Cerr<<">>> elem = "<<e<<" : Advecte phi"<<finl;
      double nor_nor_e_0 = 0.;
      for (int k=0; k<dim_esp; k++) nor_nor_e_0 += grad_phi_e_0(e,k)*grad_phi_e_0(e,k);
      if (nor_nor_e_0 <= 1.0e-10) continue; // element avec grad phi element non defini

      double contrib_e = 0.;
      for (int il=0; il<nb_som_elem; il++)
        {
          int i = elems(e,il);
          // Advection phi pour le sommet i de l element (si pas deja fait)
          if (advect_faite(i) == 0)
            {
              if (e == etarg)
                {
                  Cerr<<" ** node = "<<i<<finl;
                  Cerr<<" vecteur_deplacement = ";
                  for (int k=0; k<dim_esp; k++) Cerr<<vecteur_deplacement(i,k)<<" ";
                  Cerr<<finl;
                }

              if (level_set(i) != -99999.)
                {
                  double vdotnorm = 0.0 ;
                  for (int k=0; k<dim_esp; k++) vdotnorm += vecteur_deplacement(i,k)*grad_phi_i_0(i,k);
                  level_set(i) -= alpha*vdotnorm;
                }
              advect_faite(i) = (level_set(i)<0.?-1:1);


              if (e == etarg)
                {
                  Cerr<<" nor = ";
                  for (int k=0; k<dim_esp; k++) Cerr<<nor(i,k)<<" ";
                  Cerr<<finl;
                  Cerr<<" new level set : "<<level_set(i)<<" ; advect_faite = "<<advect_faite(i)<<finl;
                  Cerr<<" grad_phi_i_0 = ";
                  for (int k=0; k<dim_esp; k++) Cerr<<grad_phi_i_0(i,k)<<" ";
                  Cerr<<finl;
                }
            }

          // calcul de la projection solide nodale et du barycentre (moyenne arith.)
          for (int k=0; k<dim_esp; k++)
            bary(e,k) += coordsDom3D(i,k)-(level_set(i)*grad_phi_i_0(i,k)); // OP^1 = OX - phi^1 grad_phi^0
          contrib_e += 1.;
          if (e == etarg)
            {
              Cerr<<" ** node = "<<i;
              Cerr<<" ; new level set = "<<level_set(i);
              Cerr<<" ; coords = ";
              for (int k=0; k<dim_esp; k++) Cerr<<coordsDom3D(i,k)<<" ";
              Cerr<<finl;
              Cerr<<" Projection : ";
              for (int k=0; k<dim_esp; k++) Cerr<<coordsDom3D(i,k)-(level_set(i)*grad_phi_i_0(i,k))<<" ";
              Cerr<<finl;
            }
        }
      if (contrib_e != 0.)
        for (int k=0; k<dim_esp; k++) bary(e,k) /= contrib_e;

      // Si element traverse par level set zero + grad_phi_e element OK => contribution au champ aire
      int iok = 0;
      // element avec level set differents signe ?
      double level_ref = 0.;
      for (int il=0; il<nb_som_elem; il++)
        {
          int i = elems(e,il);
          if (level_set(i) != 0. && level_set(i) != -99999.)
            {
              level_ref = level_set(i); // reference 1er vertex avec phi<>0 et initialisee
              break;
            }
        }
      if (level_ref != 0.)
        {
          for (int il=0; il<nb_som_elem; il++)
            {
              int i = elems(e,il);
              if (level_set(i) != 0. && level_set(i) != -99999.)
                {
                  iok = (level_ref*level_set(i) < 0. ? 1:0) ;
                  if (iok == 1)  break;
                }
            }
        }
      // grad_phi_e element OK ?
      if (nor_nor_e_0 > 1.0e-10) iok += 1;
      if (iok == 2)
        {
          // Definition de l aire a faire une fois les rotations mises a jour
          aire(e) = 9999.;
          // Definition de la normale element
          for (int k=0; k<dim_esp; k++) new_normal_e(e, k) = grad_phi_e_0(e,k);
          if (e == etarg)
            {
              Cerr<<"aire = "<<aire(e)<<finl;
              Cerr<<"new_bary_e = ";
              for (int k=0; k<dim_esp; k++) Cerr<<bary(e,k)<<" ";
              Cerr<<finl;
              Cerr<<"new_normal_e = ";
              for (int k=0; k<dim_esp; k++) Cerr<<new_normal_e(e,k)<<" ";
              Cerr<<finl;
            }
        }
      else
        for (int k=0; k<dim_esp; k++) bary(e,k) =0.;

    }
  level_set.echange_espace_virtuel();

  // Modification champ rotation par element
  rotation = 0.; // Mise a zero des rotations
  rotation.echange_espace_virtuel();
  DoubleTrav t1EulerArr(bary);
  DoubleTrav t2EulerArr(bary);
  for (int e=0; e<nb_elem; e++)
    {
      if(aire(e) == 9999.)
        {
          prepro_lu_->computeLocalFrame(new_normal_e, t1EulerArr, t2EulerArr, e);
          prepro_lu_->computeMatRot(new_normal_e, t1EulerArr, t2EulerArr, e);
        }
    }
  rotation = prepro_lu_->get_champ_rotation();

  // Rotations à jour => modification possible du champ aire
  for (int e=0; e<nb_elem; e++)
    {
      if(aire(e) == 9999.)
        {
          aire(e) = aire_geometrique_IBM(rotation, e);
          if (e == etarg) Cerr<<"aire = "<<aire(e)<<finl;
        }
    }

  aire.echange_espace_virtuel();
  bary.echange_espace_virtuel();
  rotation.echange_espace_virtuel();

  // finalisation : mise a jour de l IB
  set_fields_to_prepro(prepro_lu_);
  prepro_lu_->compute_solid_fluid(1);

  my_interp.set_fields_from_prepro_to_interp(prepro_lu_);
  my_interp.calculer_normal_et_distance_proj_solid();
  my_interp.set_pseudo_level_set(level_set);
  // my_interp.definir_pseudo_level_set(); //re-initialisation phi

  set_variable_imposee();
  prepro_lu_-> Save_Med_File();

}

void Source_PDF_base::update_elem_IBM(DoubleTab& vecteur_deplacement, double alpha, double raid)
{
  // Update a partir des barycentres elementaires de Source_PDF_base

  int dim_esp = Objet_U::dimension;
  assert (dim_esp == 3);

  DoubleTab& aire = ref_cast_non_const(DoubleTab, champ_aire_->valeurs());
  DoubleTab& bary = ref_cast_non_const(DoubleTab, champ_barycentre_->valeurs());
  DoubleTab& rotation = ref_cast_non_const(DoubleTab, champ_rotation_->valeurs());
  const Probleme_base& pb = equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine& dom = le_dom_dis.domaine();
  int nb_elem_tot = dom.nb_elem_tot();
  const DoubleTab coordsDom3D=dom.les_sommets();

  int nb_elem = aire.dimension(0);
  assert (nb_elem == vecteur_deplacement.dimension(0));
  assert (dim_esp == vecteur_deplacement.dimension(1));
  IntTab indic_dead_cell(nb_elem);
  indic_dead_cell = 0;

  // calcul voisins de chaque element traverse
  IntLists elem_voisins(nb_elem_tot);
  compute_NeighNode_IBM_elem(aire, elem_voisins);

  // update du champ barycentre
  IntTrav toward_new_elem(nb_elem);

  // stpe1 : prediction
  DoubleTrav bary_zero(bary);
  bary_zero = bary;
  double norm_vd2;
  for (int e=0; e<nb_elem; e++)
    {
      norm_vd2 = 0.;
      for (int k=0; k<dim_esp; k++) norm_vd2 += vecteur_deplacement(e,k)*vecteur_deplacement(e,k);
      if (aire(e)>0.)
        {
          for (int k=0; k<dim_esp; k++) bary(e,k) = bary_zero(e,k) + alpha * vecteur_deplacement(e,k) ;
          if (elem_voisins[e].size() >= 2 && sqrt(norm_vd2) > 1.e-10)
            {
              for (int voi=0; voi<elem_voisins[e].size(); voi++)
                {
                  for (int k=0; k<dim_esp; k++) bary(e,k) += alpha * raid * bary_zero((elem_voisins[e])[voi],k);
                }
              for (int k=0; k<dim_esp; k++) bary(e,k) /= (1 + alpha * raid * elem_voisins[e].size());
            }
          toward_new_elem(e) = dom.chercher_elements(bary(e,0),bary(e,1),bary(e,2));
          if (aire(toward_new_elem(e)) <= 0. ) indic_dead_cell(toward_new_elem(e)) = 1;
        }
      else //pour eviter que null soit considere comme element numero o
        {
          toward_new_elem(e) =-2; // rien a deplacer
        }
    }
  indicateur_dead_cell_ = indic_dead_cell;

  // stpe2 : deplacement des champs vers des elements nouveaux
  DoubleTrav new_aire(aire);
  DoubleTrav new_bary(bary);
  DoubleTrav new_rotation(rotation);
  IntTrav contrib(nb_elem);

  for (int e=0; e<nb_elem; e++)
    {
      if (toward_new_elem(e) == -1 ) //element non trouve
        {
          Cerr<<"Source_PDF_base::update_elem_IBM: element displacement is not valid."<<finl;
          exit();
        }
      else if (toward_new_elem(e) == -2 ) // element non coupe par gamma
        {
        }
      else //deplacement des donnees (barycentres et matrices rotation)
        {
          for (int k=0; k<dim_esp; k++) new_bary(toward_new_elem(e),k) += bary(e,k);
          for (int k=0; k<rotation.dimension(1); k++) new_rotation(toward_new_elem(e),k) += rotation(e,k);
          new_aire(toward_new_elem(e)) += aire(e);
          contrib(toward_new_elem(e)) += 1;
        }
    }

  // Moyenne arithmetique pour barycentres et matrices rotation; moyenne arithmetique ou definition geometrique pour aire

  int aire_geometriq_ok = 1;

  for (int e=0; e<nb_elem; e++)
    {
      if(contrib(e)>0)
        {
          for (int k=0; k<dim_esp; k++) new_bary(e,k) /= contrib(e);
          for (int k=0; k<rotation.dimension(1); k++) new_rotation(e,k) /= contrib(e);

          if (aire_geometriq_ok)
            new_aire(e) = aire_geometrique_IBM(new_rotation, e);
          else new_aire(e) /= contrib(e);

        }
    }

  // Maj champs
  for (int e=0; e<nb_elem; e++)
    {
      aire(e) = new_aire(e);
      for (int k=0; k<dim_esp; k++) bary(e,k) = new_bary(e,k);
      for (int k=0; k<rotation.dimension(1); k++) rotation(e,k) = new_rotation(e,k);
    }

  // Verification de la conservation de la topologie
  // Traitement des elements elem ayant 0 ou 1 voisin
  // les voisins de l'antecedent de elem doivent etre voisins de elem

  // calcul voisins de chaque nouveau element traverse
  IntLists elem_voisins_new(nb_elem_tot);
  compute_NeighNode_IBM_elem(aire, elem_voisins_new);

  // recherche des elements deplaces ayant moins de 2 voisins
  for (int e=0; e<nb_elem_tot; e++)
    {
      if ( elem_voisins_new[e].size() < 2)
        {
          IntList antecedents;
          antecedents.vide();
          // list des elements (antecedents) s'etant deplace vers e
          for (int e_old=0; e_old<nb_elem_tot; e_old++)
            if (toward_new_elem(e_old) == e) antecedents.add_if_not(e_old);
//          if (antecedents.size() != 0) Cerr<<"element "<<e<<" => "<<finl;
          IntList to_link;
          to_link.vide();
          // boucle sur les antecedents de e
          for (int ant=0; ant<antecedents.size(); ant++)
            {
              // Cerr<<" voisins de antecedent = "<<antecedents[ant] <<" ";
              // boucle sur les anciens voisins de chaque antecedent de e
              int nb_vois_ant = elem_voisins[antecedents[ant]].size();
              for (int vois_ant=0; vois_ant<nb_vois_ant; vois_ant++)
                {
                  int old_vois = (elem_voisins[antecedents[ant]])[vois_ant];
                  // l'ancien voisin est il toujours voisin de e ?
                  int i_present = 0;
                  for (int vois_e=0; vois_e<elem_voisins_new[e].size(); vois_e++)
                    if ( (elem_voisins_new[e])[vois_e] == toward_new_elem(old_vois) ) i_present = 1;
                  if ( (!i_present) && (toward_new_elem(old_vois) != e) ) to_link.add_if_not(toward_new_elem(old_vois));
                  // Cerr<<old_vois<<" ";
                }
              // Cerr<<finl;

              // if (to_link.size() != 0)  Cerr<<"to_link ";
              DoubleTrav Xtest(dim_esp);
              for (int n=0; n <to_link.size(); n++)
                {
                  // Cerr<<to_link[n]<<" ";
                  // Xdeb (bary(e,0),bary(e,1),bary(e,2));
                  // Xfin (bary(to_link[n],0),bary(to_link[n],1),bary(to_link[n],2));
                  int decoup = 5;
                  for (int nb_step=0; nb_step<(decoup-1); nb_step++)
                    {
                      for (int k=0; k<dim_esp; k++) Xtest(k) = bary(e,k) + (bary(to_link[n], k) - bary(e,k))/decoup*(nb_step+1);
                      int elem_Xtest = dom.chercher_elements(Xtest(0),Xtest(1),Xtest(2));
                      if (elem_Xtest != -1)
                        {
                          if (aire(elem_Xtest) <= 0.)
                            {
                              for (int k=0; k<dim_esp; k++) bary(elem_Xtest,k) = Xtest(k);
                              for (int k=0; k<rotation.dimension(1); k++) rotation(elem_Xtest,k) = (rotation(e,k) + rotation(to_link[n],k)) / 2.;
                              aire(elem_Xtest) = aire_geometrique_IBM(rotation, elem_Xtest);
                            }
                        }
                    }
                }
              // Cerr<<finl;
            }
        }

    }

  // finalisation
  aire.echange_espace_virtuel();
  bary.echange_espace_virtuel();
  rotation.echange_espace_virtuel();


  if (prepro_lu_.non_nul())
    {
      set_fields_to_prepro(prepro_lu_);
      prepro_lu_->compute_solid_fluid(1);
      if(getInterpolationBool() == true)
        {
          Interpolation_IBM_base& my_interp = ref_cast_non_const(Interpolation_IBM_base, getInterpolationLu());
          my_interp.set_fields_from_prepro_to_interp(prepro_lu_);
        }
    }
  set_variable_imposee();
}

void Source_PDF_base::verify_results_prepro()
{
  if ( !((&prepro_lu_)->non_nul()) ) return;
  Cerr<<"(IBM) source PDF base: Prepro IBM verification if any ...."<<finl;

  const DoubleTab& aireArray = prepro_lu_->champ_aire_->valeurs();
  const DoubleTab& baryArray = prepro_lu_->champ_bary_->valeurs();
  const DoubleTab& normalArray = prepro_lu_->champ_normal_->valeurs();
  const DoubleTab& rotationArray = prepro_lu_->champ_rotation_->valeurs();
  const DoubleTab& isNodeDirichletArray = prepro_lu_->isNodeDirichlet_->valeurs();
  const DoubleTab& solideArray = prepro_lu_->solid_points_->valeurs();
  const DoubleTab& fluidArray = prepro_lu_->fluid_points_->valeurs();

  const DoubleTab& solid_elemsArray = prepro_lu_->solid_elems_->valeurs();
  const DoubleTab& fluid_elemsArray = prepro_lu_->fluid_elems_->valeurs();
  const DoubleTab& corresp_elemsArray = prepro_lu_->corresp_elems_->valeurs();

  int dim_esp = Objet_U::dimension;

  ///////////////////////////////// Champs de base /////////////////////////////////
  // Aire
  if (champ_aire_lu_.non_nul())
    {
      const DoubleTab& aireArray_src_pdf = champ_aire_lu_->valeurs();
      assert(aireArray_src_pdf.dimension(0)==aireArray.dimension(0));
      Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for aire :"<<finl;
      comp_diff_L2_max(aireArray.dimension(0), aireArray.dimension(1), aireArray_src_pdf, aireArray);
    }

  // Barycentre
  if (champ_barycentre_lu_.non_nul())
    {
      const DoubleTab& baryArray_src_pdf = champ_barycentre_lu_->valeurs();
      assert(baryArray_src_pdf.dimension(0)==baryArray.dimension(0));
      assert(baryArray_src_pdf.dimension(1)==dim_esp);
      Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for barycenter :"<<finl;
      comp_diff_L2_max(baryArray.dimension(0), dim_esp, baryArray_src_pdf, baryArray);
    }

  // Normale
  if (champ_rotation_lu_.non_nul())
    {
      const DoubleTab& rotArray_src_pdf = champ_rotation_lu_->valeurs();
      assert(rotArray_src_pdf.dimension(1)==dim_esp*dim_esp);
      assert(rotArray_src_pdf.dimension(0)==normalArray.dimension(0));
      assert(normalArray.dimension(1)==dim_esp);
      Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for normal :"<<finl;
      for (int k = 0; k < dim_esp; k++)
        {
          double err_L2_norm =0.;
          double err_max_norm = 0.;
          double L2_norm = 0.;
          for (int elem = 0; elem <normalArray.dimension(0); elem++)
            if (aireArray(elem) > 0.)
              {
                double diff_norm = abs(normalArray(elem,k) - rotArray_src_pdf(elem,3*k+(dim_esp-1)));
                err_L2_norm  += diff_norm * diff_norm ;
                L2_norm += rotArray_src_pdf(elem,3*k+(dim_esp-1)) * rotArray_src_pdf(elem,3*k+(dim_esp-1)) ;
                if (diff_norm > err_max_norm) err_max_norm = diff_norm;
              }
          if (L2_norm > 1.0e-12) err_L2_norm = sqrt(err_L2_norm / L2_norm);

          Cerr<<"    composant # "<<k<<" => "<<err_L2_norm<<" "<<err_max_norm;
        }
      Cerr<<finl;
    }

  // Rotation
  if (champ_rotation_lu_.non_nul())
    {
      const DoubleTab& rotArray_src_pdf = champ_rotation_lu_->valeurs();
      assert(rotArray_src_pdf.dimension(1)==dim_esp*dim_esp);
      assert(rotArray_src_pdf.dimension(0)==rotationArray.dimension(0));
      assert(rotationArray.dimension(1)==dim_esp*dim_esp);
      Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for rotation :"<<finl;
      comp_diff_L2_max(rotationArray.dimension(0), dim_esp*dim_esp, rotArray_src_pdf, rotationArray);
    }

  ///////////////////////////////// Interpolation /////////////////////////////////

  if(interpolation_bool_)
    {
      Cerr<<"Prepro_IBM::verify_results_prepro: Interpolation_IBM"<<finl;

      // is_dirichlet
      if (interpolation_lue_->is_dirichlet_lu_.non_nul())
        {
          const DoubleTab& isNodeDirichletArray_src_pdf = interpolation_lue_->is_dirichlet_lu_->valeurs();
          assert(isNodeDirichletArray_src_pdf.dimension(0)==isNodeDirichletArray.dimension(0));
          assert(isNodeDirichletArray_src_pdf.dimension(1)==isNodeDirichletArray.dimension(1));
          Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for isNodeDirichlet :"<<finl;
          comp_diff_L2_max(isNodeDirichletArray.dimension(0), isNodeDirichletArray.dimension(1), isNodeDirichletArray_src_pdf, isNodeDirichletArray);
        }

      // Projection solide
      if (interpolation_lue_->solid_points_lu_.non_nul())
        {
          const DoubleTab& solideArray_src_pdf = interpolation_lue_->solid_points_lu_->valeurs();
          assert(solideArray_src_pdf.dimension(0)==solideArray.dimension(0));
          assert(solideArray_src_pdf.dimension(1)==dim_esp);
          Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for solid projection :"<<finl;
          comp_diff_L2_max(solideArray.dimension(0), dim_esp, solideArray_src_pdf, solideArray);
        }

      // correspondance elems
      if (interpolation_lue_->corresp_elems_lu_.non_nul())
        {
          const DoubleTab& corresp_elems_src_pdf = interpolation_lue_->corresp_elems_lu_->valeurs();
          assert(corresp_elems_src_pdf.dimension(1)==1);
          assert(corresp_elems_src_pdf.dimension(0)==corresp_elemsArray.dimension(0));
          assert(corresp_elemsArray.dimension(1)==1);
          Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for correspondance elems :"<<finl;
          comp_diff_L2_max(corresp_elemsArray.dimension(0), 1, corresp_elems_src_pdf, corresp_elemsArray);
        }

      Cerr<<"interpolation_lue_->que_suis_je() = "<<interpolation_lue_->que_suis_je()<<finl;
      // Interpolation_IBM_mean_gradient
      if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_gradient_moyen")
        {
          Cerr<<"Prepro_IBM::verify_results_prepro: >> mean_gradient"<<finl;
          Interpolation_IBM_mean_gradient& interp = ref_cast(Interpolation_IBM_mean_gradient,interpolation_lue_.valeur());

          // Element projection solide
          if (interp.solid_elems_lu_.non_nul())
            {
              const DoubleTab& solid_elemsArray_src_pdf = interp.solid_elems_lu_->valeurs();
              assert(solid_elemsArray_src_pdf.dimension(0)==solid_elemsArray.dimension(0));
              assert(solid_elemsArray_src_pdf.dimension(1)==solid_elemsArray.dimension(1));
              Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for solid projection element :"<<finl;
              comp_diff_L2_max(solid_elemsArray.dimension(0), solid_elemsArray.dimension(1), solid_elemsArray_src_pdf, solid_elemsArray);
            }
        }

      // Interpolation_IBM_elem_fluid
      if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_element_fluide")
        {
          Cerr<<"Prepro_IBM::verify_results_prepro: >> elem_fluide"<<finl;
          Interpolation_IBM_elem_fluid& interp = ref_cast(Interpolation_IBM_elem_fluid,interpolation_lue_.valeur());

          // Projection fluide
          if (interp.fluid_points_lu_.non_nul())
            {
              const DoubleTab& fluidArray_src_pdf = interp.fluid_points_lu_->valeurs();
              assert(fluidArray_src_pdf.dimension(0)==fluidArray.dimension(0));
              assert(fluidArray_src_pdf.dimension(1)==dim_esp);
              Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for fluid projection :"<<finl;
              comp_diff_L2_max(fluidArray.dimension(0), dim_esp, fluidArray_src_pdf, fluidArray);
            }

          // Element projection fluide
          if (interp.fluid_elems_lu_.non_nul())
            {
              const DoubleTab& fluid_elemsArray_src_pdf = interp.fluid_elems_lu_->valeurs();
              assert(fluid_elemsArray_src_pdf.dimension(0)==fluid_elemsArray.dimension(0));
              assert(fluid_elemsArray_src_pdf.dimension(1)==fluid_elemsArray.dimension(1));
              Cerr<<"Prepro_IBM::verify_results_prepro: relative L2-norm and abs max-norm errors for fluid projection element :"<<finl;
              comp_diff_L2_max(fluid_elemsArray.dimension(0), fluid_elemsArray.dimension(1), fluid_elemsArray_src_pdf, fluid_elemsArray);
            }
        }

    }
}

void Source_PDF_base::filtre_CLD(DoubleTab& filtre) const
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::rotate_imposed_velocity(DoubleTab& vitesse_imposee)
{
  const Domaine_dis_base& Domaine_dis = equation().probleme().domaine_dis();
  int nb_elem_tot=Domaine_dis.domaine().nb_elem_tot();
  DoubleTab& rotation = champ_rotation_->valeurs();

  // VDF/VEF by default
  const Domaine_VF& the_dom = ref_cast(Domaine_VF,Domaine_dis);
  const IntTab& elems = (equation().discretisation().is_ef() ? Domaine_dis.domaine().les_elems() : the_dom.elem_faces());
  int nb_som_tot = the_dom.nb_faces_tot();
  int nb_dof_elem = Domaine_dis.domaine().nb_faces_elem();
  if (equation().discretisation().is_ef())
    {
      nb_dof_elem = Domaine_dis.domaine().nb_som_elem();
      nb_som_tot= Domaine_dis.domaine().nb_som_tot();
    }
  else if (!(equation().discretisation().is_vdf()) && !(equation().discretisation().is_vef()))
    {
      Cerr<<"Source_PDF_base::rotate_imposed_velocit: discretisation is not EF, VEF, VDF : "<<equation().discretisation()<<finl;
    }

  // COPIES
  DoubleTab vitesse = vitesse_imposee;
  assert(Objet_U::dimension==3);
  int dim_var=vitesse_imposee.dimension(1);
  assert(Objet_U::dimension==dim_var);

  //vitesse_imposee.echange_espace_virtuel();
  DoubleTrav marqueur=equation().probleme().get_champ("vitesse").valeurs();
  vitesse_imposee = 0.;
  marqueur = 0.;

  for (int num_elem=0; num_elem<nb_elem_tot; num_elem++)
    {
      double norm =  sqrt(rotation(num_elem,0)*rotation(num_elem,0)
                          +rotation(num_elem,1)*rotation(num_elem,1)
                          +rotation(num_elem,2)*rotation(num_elem,2)
                          +rotation(num_elem,3)*rotation(num_elem,3)
                          +rotation(num_elem,4)*rotation(num_elem,4)
                          +rotation(num_elem,5)*rotation(num_elem,5)
                          +rotation(num_elem,6)*rotation(num_elem,6)
                          +rotation(num_elem,7)*rotation(num_elem,7)
                          +rotation(num_elem,8)*rotation(num_elem,8));
      if (norm > 1e-6)
        {
          for (int i=0; i<nb_dof_elem; i++)
            {
              int s = elems(num_elem,i);
              for (int c=0; c<dim_var; c++)
                {
                  for (int k=0; k<dim_var; k++)
                    {
                      vitesse_imposee(s,c) += rotation(num_elem,3*c+k)*vitesse(s,k);
                      //Cerr << "k = " << k << ", c = " << c << ", 3*c+k = " << 3*c+k << ", rotation = " << rotation(num_elem,3*c+k) << finl;
                    }
                }
              marqueur(s,0) += 1.0;
            }
        }
    }
  for (int i=0; i<nb_som_tot; i++)
    {
      if(marqueur(i,0) > 0.)
        {
          for (int c=0; c<dim_var; c++)
            {
              vitesse_imposee(i,c) /= marqueur(i,0);
            }
        }
    }
}

void Source_PDF_base::compute_variable_imposee_projete(const DoubleTab& marqueur, const DoubleTab& points, double val, double eps)
{
  const Domaine_dis_base& le_dom_dis = equation().probleme().domaine_dis();

  int nb_som=le_dom_dis.domaine().nb_som();
  int dim = Objet_U::dimension;
  const DoubleTab& coords = le_dom_dis.domaine().coord_sommets();
  ArrOfDouble x(dim);
  int dim_var = equation().inconnue().valeurs().dimension(1);
  for (int i = 0; i < nb_som; i++)
    {
      if ((marqueur(i) >= val + eps) || (marqueur(i) <= val - eps))
        {
          for (int j = 0; j < dim; j++)
            {
              x[j] = points(i,j);
            }
          for (int j = 0; j < dim_var; j++)
            {
              modele_lu_.variable_imposee_->valeurs()(i,j) = modele_lu_.get_variable_imposee(x,j);
            }
        }
      else
        {
          for (int j = 0; j < dim; j++)
            {
              x[j] = coords(i,j);
            }
          for (int j = 0; j < dim_var; j++)
            {
              modele_lu_.variable_imposee_->valeurs()(i,j) = modele_lu_.get_variable_imposee(x,j);
            }
        }
    }
}

void Source_PDF_base::associer_domaines(const Domaine_dis_base& domaine_dis,
                                        const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::compute_indicateur_nodal_champ_aire()
{
  // indicateur_nodal_champ_aire_ = 1 si sommet appartient a un element dont l'aire <> 0
  // 0 sinon
  const DoubleTab& aire=champ_aire_->valeurs();
  const Domaine_dis_base& Domaine_dis = equation().probleme().domaine_dis();
  int nb_elems=Domaine_dis.domaine().nb_elem_tot();
  int nb_nodes=Domaine_dis.domaine().nb_som_tot();
  int nb_som_elem=Domaine_dis.domaine().nb_som_elem();
  const IntTab& elems= Domaine_dis.domaine().les_elems() ;

  DoubleTab indic(nb_nodes);
  indic = 0.;
  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      if (aire(num_elem)>0.)
        {
          for (int i=0; i<nb_som_elem; i++)
            {
              int s1=elems(num_elem,i);
              indic(s1) = 1.;
            }
        }
    }
  indicateur_nodal_champ_aire_ = indic;
}

double Source_PDF_base::fonct_regul_PDF(const int elem, const double dist)
{
  double coeff = 1.;

  int size_dead_cell_tab = indicateur_dead_cell_.size();
  if (size_dead_cell_tab == 0) return coeff;
  if (elem >= size_dead_cell_tab)
    {
      Cerr<<"Source_PDF_base::fonct_regul_PDF: element = "<<elem<<" >= size of tab indicateur_dead_cell = "<<size_dead_cell_tab<<finl;
      exit();
    }
  if (indicateur_dead_cell_(elem) != 1) return coeff;

  const Probleme_base& pb = equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine_VF& the_dom_VF = ref_cast(Domaine_VF,le_dom_dis);
  DoubleTab& rotation = ref_cast_non_const(DoubleTab, champ_rotation_->valeurs());

  double h_max_elem = the_dom_VF.volumes(elem) / aire_geometrique_IBM(rotation, elem);
  double coeff_regul_PDF = modele_lu_.regul_coeff_PDF_;
  double h_ref = coeff_regul_PDF * h_max_elem;
  coeff = max((h_ref - dist)/h_ref , 0.);

  return coeff;
}

double Source_PDF_base::fonct_coeff(const double rho_m, const double aire, const double dt) const
{
  // return = 0 si aire dans element <= 0
  // return = rho/dt * coeff_relax sinon
  double val_coeff = 0.;
  if (aire<=0.)
    {
      return val_coeff;
    }

  double inv_dt = 1./dt;
  // Cerr<<"dt pour Source_PDF_base::fonct_coeff : "<< dt <<finl;
  val_coeff = rho_m * inv_dt;
  const double tps_cour = equation().probleme().schema_temps().temps_courant();
  if (temps_relax_ != 1.0e+12)
    {
      // double bco = 100.*(1.-exp(-0.1/echelle_relax_)); // tps_cour_intersection = 0.1 * temps_relax_
      double deltrel = tps_cour - temps_relax_;
      double coeff_relax = (1.+tanh(deltrel/(echelle_relax_ * temps_relax_)))/2.;
      // double coeff_relax = min (1., std::max(1.-bco*deltrel*deltrel/temps_relax_/temps_relax_, exp(deltrel/(echelle_relax_ * temps_relax_))));
      // double coeff_relax = 1.0 - exp(-tps_cour / temps_relax_);
      // double coeff_relax = tanh(tps_cour / temps_relax_);
      val_coeff *= coeff_relax ;
      // Cerr<< "coeff. relax. pour PDF = "<<coeff_relax<<" val_coeff = "<<val_coeff <<endl;
    }
  return val_coeff;
}

// TENSOR CALCULATION
// Pour pouvoir eventuellement traiter differements les directions d'espace
ArrOfDouble Source_PDF_base::get_tuvw_local() const
{
  assert(Objet_U::dimension==3);
  ArrOfDouble tuvw(dimension) ;
  tuvw[0] = 1.0 / modele_lu_.eta_;
  tuvw[1] = 1.0 / modele_lu_.eta_;
  tuvw[2] = 1.0 / modele_lu_.eta_;
  return tuvw ;
}

DoubleVect Source_PDF_base::diag_coeff_elem(ArrOfDouble& variable_elem, const DoubleTab& rotation, int num_elem) const
{
  assert(Objet_U::dimension==3);
  ArrOfDouble tuvw(dimension);

  tuvw = get_tuvw_local();

  ArrOfDouble sum_dir_loc(dimension) ;
  for (int i=0; i<dimension; i++)
    {
      sum_dir_loc[i]=0.;
      for (int k=0; k<dimension; k++)
        sum_dir_loc[i]+=rotation(num_elem,3*i+k);
    }

  DoubleVect diag_coef(dimension);
  for (int comp=0; comp<dimension; comp++)
    {
      diag_coef(comp) = 0. ;
      for (int k=0; k<dimension; k++)
        diag_coef(comp)+=tuvw[k]*sum_dir_loc[k]*rotation(num_elem,3*k+comp);
    }

  return diag_coef ;
}

DoubleTab Source_PDF_base::compute_coeff_elem() const
{
  // coeff = 1 si aire dans element <= 0
  // coeff = 1 + (Ksi/eta) * coeff_relax sinon
  const Domaine_dis_base& domaine_dis = equation().probleme().domaine_dis();
  const IntTab& elems= domaine_dis.domaine().les_elems() ;
  int nb_som_elem=domaine_dis.domaine().nb_som_elem();
  int nb_elems=domaine_dis.domaine().nb_elem_tot();
  const DoubleTab& rotation=champ_rotation_->valeurs();
  const DoubleTab& aire = champ_aire_->valeurs();

  int dim_var = equation().inconnue().valeurs().dimension(1);
  ArrOfDouble variable_elem(dim_var);
  int dim_esp = Objet_U::dimension ;
  if (dim_var != 1 && dim_var != dim_esp)
    {
      Cerr<<"compute_coeff_elem: dimension variable differente 1 ou "<<dim_esp<<"! "<<finl;
      exit();
    }
  DoubleTab coeff(nb_elems, dim_var);

  double dt = (ref_cast_non_const(Source_PDF_base, *this)).calcul_dt_pdf();
  const DoubleTab& rho_m=champ_rho_->valeurs();
  if (equation().probleme().schema_temps().temps_courant()==0)
    {
      dt = 1.0;
    }
  // Cerr<<"dt pour compute_coeff_elem : "<< dt <<finl;

  double val_coeff ;
  DoubleVect val_diag_coef(dim_var);
  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      if (aire(num_elem)<=0.)
        {
          for (int comp=0; comp<dim_var; comp++) coeff(num_elem, comp) = 1. ;
        }
      else
        {
          const DoubleTab& variable=equation().inconnue().valeurs();
          variable_elem=0;
          for (int s=0; s<nb_som_elem; s++)
            {
              int som_glob=elems(num_elem,s);
              for (int comp=0; comp<dim_var; comp++)
                variable_elem[comp]+=variable(som_glob,comp);
            }
          // rho/dt * coeff_relax
          val_coeff = fonct_coeff(rho_m(num_elem), aire(num_elem), dt) ;
          // 1/eta
          if (dim_var == dim_esp)
            {
              val_diag_coef = diag_coeff_elem(variable_elem, rotation, num_elem) ;
            }
          else
            {
              val_diag_coef(0) = 1.0 / modele_lu_.eta_ ;
            }
          // coeff = 1 + (Ksi/eta) * coeff_relax
          for (int comp=0; comp<dim_esp; comp++)
            coeff(num_elem, comp) = 1.+ val_coeff * val_diag_coef(comp) * dt / rho_m(num_elem) ;
        }
    }

  return coeff;
}

DoubleTab Source_PDF_base::compute_coeff_matrice() const
{
  // coeff sommet = Sigma_elem coeff_elem / Nb contributions;  avec :
  // coeff_elem = Ksi/eta * coeff_relax si element dont l'aire <>= 0
  // coeff_elem = 0 sinon
  const Domaine_dis_base& Domaine_dis = equation().probleme().domaine_dis();
  int nb_elems=Domaine_dis.domaine().nb_elem_tot();

  // VDF/VEF by default
  const Domaine_VF& the_dom = ref_cast(Domaine_VF,Domaine_dis);
  const IntTab& elems = (equation().discretisation().is_ef() ? Domaine_dis.domaine().les_elems() : the_dom.elem_faces());
  int nb_dof_tot = the_dom.nb_faces_tot();
  int nb_dof_elem = Domaine_dis.domaine().nb_faces_elem();
  if (equation().discretisation().is_ef())
    {
      nb_dof_tot= Domaine_dis.domaine().nb_som_tot();
      nb_dof_elem = Domaine_dis.domaine().nb_som_elem();
    }
  else if (!(equation().discretisation().is_vdf()) && !(equation().discretisation().is_vef()))
    {
      Cerr<<"Source_PDF_base::compute_coeff_matrice: discretisation is not EF, VEF, VDF : "<<equation().discretisation()<<finl;
    }

  const DoubleTab& rotation=champ_rotation_->valeurs();
  const DoubleTab& aire = champ_aire_->valeurs();

  int dim_esp = Objet_U::dimension ;
  const DoubleTab& variable=equation().inconnue().valeurs();
  int dim_var = variable.dimension(1);
  if (dim_var != 1 && dim_var != dim_esp)
    {
      Cerr<<"compute_coeff_matrice: dimension variable differente 1 ou "<<dim_esp<<"! "<<finl;
      exit();
    }
  ArrOfDouble variable_elem(dim_var);
  DoubleTab coeff(variable);
  IntTab contrib(nb_dof_tot, dim_var);
  contrib = 0;
  coeff = 0.;

  double dt = (ref_cast_non_const(Source_PDF_base, *this)).calcul_dt_pdf();
  const DoubleTab& rho_m=champ_rho_->valeurs();
  // Cerr<<"dt pour compute_coeff_matrice : "<< dt <<finl;

  double val_coeff ;
  DoubleVect val_diag_coef(dim_var), coeff_el(dim_var) ;
  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      if (aire(num_elem)<=0.)
        {
          for (int comp=0; comp<dim_var; comp++) coeff_el(comp) = 0. ;
        }
      else
        {
          variable_elem=0.;
          for (int s=0; s<nb_dof_elem; s++)
            {
              int dof_glob=elems(num_elem,s);
              for (int comp=0; comp<dim_var; comp++)
                variable_elem[comp]+=variable(dof_glob,comp)/nb_dof_elem;
            }
          // rho/dt * coeff_relax
          val_coeff = fonct_coeff(rho_m(num_elem), aire(num_elem), dt) ;
          //  1/eta
          if (dim_var == dim_esp)
            {
              val_diag_coef = diag_coeff_elem(variable_elem, rotation, num_elem) ;
            }
          else
            {
              val_diag_coef(0) = 1.0 / modele_lu_.eta_ ;
            }
          // coeff = Ksi/eta * coeff_relax
          for (int comp=0; comp<dim_var; comp++)
            coeff_el(comp) = val_coeff * val_diag_coef(comp) * dt / rho_m(num_elem) ;
        }
      for (int i = 0; i<nb_dof_elem; i++)
        {
          int s = elems(num_elem,i);
          for (int comp=0; comp<dim_var; comp++)
            {
              if (coeff_el(comp) > 0.)
                {
                  coeff(s,comp)+= coeff_el(comp);
                  contrib(s, comp) += 1;
                }
              // Cerr<<"Source_PDF_EF: coeff_el ("<<comp<<"): "<<" "<<coeff_el(comp)<<finl;
            }
        }
    }

  double val = 0.;
  for (int s = 0; s<nb_dof_tot; s++)
    {
      for (int comp=0; comp<dim_var; comp++)
        {
          val = 1.;
          if (contrib(s,comp) != 0)
            {
              val += coeff(s,comp)/contrib(s,comp);
            }
          coeff(s,comp) = val;
          // Cerr<<"Source_PDF_EF: contrib coeff ("<<s<<","<<comp<<"): "<<contrib(s,comp)<<" "<<coeff(s,comp)<<finl;
        }
    }
  // Cerr << "Coeff matrice" << coeff << finl;
  return coeff;
}

void Source_PDF_base::multiply_coeff_volume(DoubleTab& coeff) const
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

DoubleTab Source_PDF_base::compute_pond(const DoubleTab& rho_m, const DoubleTab& aire, const DoubleVect& volume, int& coef1, int& nb_elems) const
{
// return = 0 si aire dans element <= 0
// return = rho/dt * coeff_relax * volume / (coef1*coef1)  sinon
  DoubleTab pond = rho_m;
  double inv_coef1 = 1. / coef1;

  double dt = (ref_cast_non_const(Source_PDF_base, *this)).calcul_dt_pdf();

  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    {
      pond(num_elem) = fonct_coeff(rho_m(num_elem), aire(num_elem), dt) ;
      pond(num_elem) *= volume(num_elem)*inv_coef1*inv_coef1;
    }

  return pond;
}

DoubleTab& Source_PDF_base::ajouter(DoubleTab& secmem) const
{
  if(has_interface_blocs())
    {
      ajouter_blocs({}, secmem);
      return secmem;
    }
  return ajouter_(variable_imposee_,secmem);
}

DoubleTab& Source_PDF_base::ajouter_(const DoubleTab& variable, DoubleTab& resu) const
{
  const int i_traitement_special = 0 ;
  ajouter_(variable, resu, i_traitement_special) ;
  return resu;
}

DoubleTab& Source_PDF_base::ajouter_(const DoubleTab& variable, DoubleTab& resu, const int i_traitement_special) const
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
  return resu;
}

void  Source_PDF_base::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

DoubleTab& Source_PDF_base::calculer(DoubleTab& resu) const
{
  resu = 0;
  const DoubleTab& variable=equation().inconnue().valeurs();
  return ajouter_(variable,resu);
}

DoubleVect& Source_PDF_base::compute_source_term_PDF(DoubleVect& bilan)
{
  int i_traitement_special = 0;
  const DoubleTab& variable=equation().inconnue().valeurs();
  DoubleTrav add_zero_term(variable);
  add_zero_term = 0.;
  return compute_source_term_PDF(i_traitement_special, add_zero_term, bilan);
}

DoubleVect& Source_PDF_base::compute_source_term_PDF(int i_traitement_special, DoubleTab& add_term, DoubleVect& bilan)
{
  assert(Objet_U::dimension <= 3);
  // double temps_corrige = temps_computation_pdf_ + dt_computation_pdf_;
  // Cerr<<"compute_source_term_PDF:: temps courant, temps computation_pdf et temps corrige = "<<equation().probleme().schema_temps().temps_courant() <<" "<<temps_computation_pdf_<<" "<<temps_corrige<<finl;
  DoubleTab& variable=equation().inconnue().valeurs();
  // Si temps courant n'a pas ete mis a jour, on prend temps futur pour inconnue
  if (equation().probleme().schema_temps().temps_courant() == temps_computation_pdf_) variable=equation().inconnue().futur();
  int nb_som= variable.dimension(0);
  int nb_comp = variable.dimension(1);
  assert(nb_comp == modele_lu_.dim_variable_);

  DoubleTrav resu(variable);
  calculer(resu, i_traitement_special);
  // On s'assure d'avoir le meme dt que pour sec_mem_pdf (calculer_pdf)
  double dt_courant = calcul_dt_pdf();
  resu *= dt_courant/dt_computation_pdf_;
  DoubleTrav flag_cl(resu);
  filtre_CLD(flag_cl);

  source_term_PDF.resize(sec_mem_pdf.dimension(0), sec_mem_pdf.dimension(1));
  source_term_PDF = 0.;

  DoubleVect source_term[3];
  bilan.resize(equation().inconnue().nb_comp());
  bilan = 0.0;
  Cerr<<"(IBM) Source_PDF_base: Bilan terme PDF = ";
  for (int i=0; i<nb_comp; i++)
    {
      source_term[i] = champ_nodal_->valeurs();
      source_term[i] = 0.;
      for (int j=0; j<nb_som; j++)
        {
          double  filter = (std::fabs(resu(j,i)) > 0?1.:0.);
          filter *= flag_cl(j,i);
          source_term[i](j) = (resu(j,i) - sec_mem_pdf(j,i) + add_term(j,i) )*filter;
          source_term_PDF(j,i) = source_term[i](j);
          // if(  (abs(resu(j,i)) > 1.0) )
          //   {
          //     Cerr<<"i,j = "<<i<<" "<<j<<" ; resu = "<<resu(j,i)<<" ; sec_mem_pdf = "<<sec_mem_pdf(j,i)<<" ; bilan = "<<source_term[i](j)<<endl;
          //     Cerr<<"added term (secmem_conv) = "<<add_term(j,i)<<endl;
          //   }
        }
      bilan(i) = mp_somme_vect(source_term[i]);
      Cerr<< bilan(i)<<" ";
    }
  source_term_PDF.echange_espace_virtuel();
  Cerr<<finl;
  return bilan;
}

void Source_PDF_base::volume_source_term_PDF(DoubleTab& volume_term_PDF)
{
}

DoubleTab& Source_PDF_base::calculer(DoubleTab& resu, const int i_traitement_special) const
{
  resu = 0;
  const DoubleTab& variable=equation().inconnue().valeurs();
  return ajouter_(variable,resu, i_traitement_special);
}

void Source_PDF_base::calculer_variable_imposee_elem_fluid()
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::calculer_variable_imposee_mean_grad()
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::calculer_variable_imposee_hybrid()
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::calculer_vitesse_imposee_power_law_tbl()
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::calculer_vitesse_imposee_power_law_tbl_u_star()
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::calculer_temperature_imposee_wall_law()
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::calculer_variable_imposee()
{
  if (interpolation_bool_)
    {
      if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_element_fluide")
        calculer_variable_imposee_elem_fluid();
      else if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_gradient_moyen")
        calculer_variable_imposee_mean_grad();
      else if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_hybride")
        calculer_variable_imposee_hybrid();
      else if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_power_law_tbl")
        calculer_vitesse_imposee_power_law_tbl();
      else if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_power_law_tbl_u_star")
        calculer_vitesse_imposee_power_law_tbl_u_star();
      else if (interpolation_lue_->que_suis_je() == "Interpolation_IBM_thermal_wall_law")
        calculer_temperature_imposee_wall_law();
    }
}

DoubleTab& Source_PDF_base::calculer_pdf(DoubleTab& resu) const
{
  resu = 0;
  ajouter_(variable_imposee_,resu);

  double my_dt = (ref_cast_non_const(Source_PDF_base, *this)).calcul_dt_pdf();
  (ref_cast_non_const(Source_PDF_base, *this)).set_dt_computation_pdf(my_dt);
  double my_temps = equation().probleme().schema_temps().temps_courant();
  (ref_cast_non_const(Source_PDF_base, *this)).set_temps_computation_pdf(my_temps);

  return resu;
}

double Source_PDF_base::calcul_dt_pdf()
{
  double dt_ref = equation().probleme().schema_temps().pas_de_temps();
  double dt_min = equation().probleme().schema_temps().pas_temps_min();
  double my_dt = std::max(dt_ref,dt_min);
  if (equation().probleme().schema_temps().temps_courant()==0) my_dt = 1.;
  return my_dt;
}

void Source_PDF_base::mettre_a_jour(double temps)
{
  //la_source->mettre_a_jour(temps);
}

void Source_PDF_base::correct_variable(const DoubleTab& coeff_node, DoubleTab& variable) const
{
  int nb_data_vit = variable.size();
  int nb_node_2 = coeff_node.dimension(0) ;
  int ncomp = coeff_node.dimension(1) ;
  if (nb_data_vit != nb_node_2 * ncomp)
    {
      Cerr<<"Source_PDF_base: numbers of nodes are different for variable and coeff_node."<<finl;
      abort();
    }
  const DoubleTab& vit_dir=variable_imposee_;
  for(int i=0; i<nb_node_2; i++)
    {
      for (int j=0; j<ncomp; j++)
        {
          if (coeff_node(i,j)!= 1.0)
            {
              variable(i,j)=vit_dir(i,j);
            }
          // Cerr<<"Source_PDF_base: coeff_node variable("<<i<<","<<j<<"): "<<coeff_node(i,j)<<" "<<variable(i,j)<<finl;
        }
    }
}

int Source_PDF_base::impr(Sortie& os) const
{
  Cerr << "Source_PDF_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
  return 0;
}

void Source_PDF_base::ouvrir_fichier(SFichier& os, const Nom& type, const int flag=1) const
{
  // flag nul on n'ouvre pas le fichier
  if (flag==0)
    return ;

  const Probleme_base& pb=equation().probleme();
  const Schema_Temps_base& sch=pb.schema_temps();
  const int precision=sch.precision_impr();
  Nom nomfichier(out_);
  if (type!="") nomfichier+=(Nom)"_"+type;
  nomfichier+=".out";

  // On cree le fichier a la premiere impression avec l'en tete
  if (sch.nb_impr()==1 && !pb.reprise_effectuee())
    {
      os.ouvrir(nomfichier);
      SFichier& fic=os;
      Nom espace="\t\t";
      fic << (Nom)"# Printing of the source term "+que_suis_je()+" of the equation "+equation().que_suis_je()+" of the problem "+equation().probleme().le_nom() << finl;
      fic << "# " << description() << finl;
      if(modele_lu_.dim_variable_ == Objet_U::dimension)
        {
          // vector
          assert(Objet_U::dimension==3);
          fic << "# Time" << espace << "Fx" << espace << "Fy" << espace << "Fz";
        }
      else if(modele_lu_.dim_variable_ == 1)
        {
          // scalar
          fic << "# Time" << espace << "Sum";
        }
      else
        {
          Cerr << "Source_PDF_base: for scalar or vector only; dim = " << modele_lu_.dim_variable_ << finl;
          Process::exit();
        }
      fic << finl;
    }
  // Sinon on l'ouvre
  else
    os.ouvrir(nomfichier,ios::app);

  os.precision(precision);
  os.setf(ios::scientific);
}

void Source_PDF_base::updateChampRho()
{
  if (equation().probleme().que_suis_je() == "Pb_Melange")
    champ_rho_->affecter(equation().probleme().get_champ("masse_volumique_melange"));
  else
    champ_rho_->affecter(equation().probleme().get_champ("masse_volumique"));
}

void Source_PDF_base::correct_incr_pressure(const DoubleTab& coeff_node, DoubleTab& correction_en_pression) const
{
  Cerr << "Source_PDF_NS_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::correct_pressure(const DoubleTab& coeff_node, DoubleTab& pression, const DoubleTab& correction_en_pression) const
{
  Cerr << "Source_PDF_NS_base: Not implemented for current discretisation. Aborting..." << finl;
  abort();
}

void Source_PDF_base::creer_champ(const Motcle& motlu)
{
  const Probleme_base& pb = equation().probleme();
  const DoubleTab& variable=equation().inconnue().valeurs();
  int nb_comp = variable.dimension(1);
  Motcle directive("temperature");
  if (equation().inconnue().nature_du_champ()==vectoriel) directive="vitesse";
  const Domaine_dis_base& le_dom_dis = equation().domaine_dis();

  double temps=0.;
  if (motlu=="source_term_pdf" && !champ_source_term_PDF_.non_nul())
    {
      if (nb_comp == 1)
        {
          Noms noms(1);
          noms[0]="source_term_PDF";
          Noms unites(1);
          unites[0] = "-";
          pb.discretisation().discretiser_champ(directive,le_dom_dis,scalaire,noms,unites,nb_comp,temps,champ_source_term_PDF_);
        }
      else if (nb_comp == Objet_U::dimension)
        {
          Noms noms(nb_comp);
          noms[0]="source_term_PDF";
          Noms unites(nb_comp);
          unites[0] = "-";
          pb.discretisation().discretiser_champ(directive,le_dom_dis,vectoriel,noms,unites,nb_comp,temps,champ_source_term_PDF_);
        }
      champs_compris_.ajoute_champ(champ_source_term_PDF_);
    }

  if (motlu=="barycentre_IBM" && !champ_barycentre_IBM_.non_nul())
    {
      nb_comp = Objet_U::dimension;
      Noms nomsb(nb_comp);
      nomsb[0]="barycentre_IBM";
      Noms unitesb(nb_comp);
      unitesb[0] = "m";
      pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nomsb,unitesb,nb_comp,0., champ_barycentre_IBM_);
      champs_compris_.ajoute_champ(champ_barycentre_IBM_);
    }

  if (motlu=="normal_IBM" && !champ_normal_IBM_.non_nul())
    {
      nb_comp = Objet_U::dimension;
      Noms nomsn(nb_comp);
      nomsn[0]="normal_IBM";
      Noms unitesn(nb_comp);
      unitesn[0] = "m";
      pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nomsn,unitesn,nb_comp,0., champ_normal_IBM_);
      champs_compris_.ajoute_champ(champ_normal_IBM_);
    }

  if (motlu=="aire_IBM" && !champ_aire_IBM_.non_nul())
    {
      pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,"aire_IBM","m^2",1,0., champ_aire_IBM_);
      champs_compris_.ajoute_champ(champ_aire_IBM_);
    }

  if (motlu=="velocity_shape_IBM" && !champ_vitesse_shape_IBM_.non_nul())
    {
      nb_comp = Objet_U::dimension;
      Noms nomsv(nb_comp);
      for (int i=0; i< nb_comp; i++) nomsv[i]="velocity_shape_IBM";
      Noms unitesv(nb_comp);
      for (int i=0; i< nb_comp; i++) unitesv[i] = "m/s";
      pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nomsv,unitesv,nb_comp,0., champ_vitesse_shape_IBM_);
      champs_compris_.ajoute_champ(champ_vitesse_shape_IBM_);
    }

  if (motlu=="pseudo_level_set_IBM" && !champ_pseudo_level_set_IBM_.non_nul())
    {
      nb_comp = 1;
      Noms nompl(nb_comp);
      for (int i=0; i< nb_comp; i++) nompl[i]="pseudo_level_set_IBM";
      Noms unitepl(nb_comp);
      for (int i=0; i< nb_comp; i++) unitepl[i] = "m";
      Motcle directiveps("temperature");
      pb.discretisation().discretiser_champ(directive,le_dom_dis,scalaire,nompl,unitepl,nb_comp,0., champ_pseudo_level_set_IBM_);
      champs_compris_.ajoute_champ(champ_pseudo_level_set_IBM_);
    }
}

void Source_PDF_base::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("source_term_pdf");
  noms_compris.add("barycentre_IBM");
  noms_compris.add("normal_IBM");
  noms_compris.add("aire_IBM");
  noms_compris.add("velocity_shape_IBM");
  noms_compris.add("pseudo_level_set_IBM");
  if (opt==DESCRIPTION)
    Cerr<<" Source_PDF_base : "<< noms_compris <<finl;
  else
    nom.add(noms_compris);
}

bool Source_PDF_base::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  if (nom == "source_term_pdf" && champ_source_term_PDF_.non_nul())
    {
      ref_champ = get_champ(nom);
      return true;
    }
  else if (nom == "barycentre_IBM" && champ_barycentre_IBM_.non_nul())
    {
      ref_champ = get_champ(nom);
      return true;
    }
  else if (nom == "normal_IBM" && champ_normal_IBM_.non_nul())
    {
      ref_champ = get_champ(nom);
      return true;
    }
  else if (nom == "aire_IBM" && champ_aire_IBM_.non_nul())
    {
      ref_champ = get_champ(nom);
      return true;
    }
  else if (nom == "velocity_shape_IBM" && champ_vitesse_shape_IBM_.non_nul())
    {
      ref_champ = get_champ(nom);
      return true;
    }
  else if (nom == "pseudo_level_set_IBM" && champ_pseudo_level_set_IBM_.non_nul())
    {
      ref_champ = get_champ(nom);
      return true;
    }
  else
    return false; /* rien trouve */
}

bool Source_PDF_base::has_champ(const Motcle& nom) const
{
  if (nom == "source_term_pdf" && champ_source_term_PDF_.non_nul())
    return true;
  else if (nom == "barycentre_IBM" && champ_barycentre_IBM_.non_nul())
    return true;
  else if (nom == "normal_IBM" && champ_normal_IBM_.non_nul())
    return true;
  else if (nom == "aire_IBM" && champ_aire_IBM_.non_nul())
    return true;
  else if (nom == "velocity_shape_IBM" && champ_vitesse_shape_IBM_.non_nul())
    return true;
  else if (nom == "pseudo_level_set_IBM" && champ_pseudo_level_set_IBM_.non_nul())
    return true;
  else
    return false; /* rien trouve */
}

const Champ_base& Source_PDF_base::get_champ(const Motcle& nom) const
{
  int dim_esp =  Objet_U::dimension;

  if (nom=="source_term_pdf")
    {
      if (champ_source_term_PDF_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ_source_term_PDF_
      DoubleTab& valeurs = champ_source_term_PDF_->valeurs();
      valeurs=0.;
      if (source_term_PDF.size_array()>0)
        {
          int nb_d0 = source_term_PDF.dimension(0);
          int nb_d1 = source_term_PDF.dimension(1);
          for (int num0=0; num0<nb_d0; num0++)
            for (int num1=0; num1<nb_d1; num1++)
              valeurs(num0,num1)=source_term_PDF(num0,num1);
        }
      valeurs.echange_espace_virtuel();
      champ_source_term_PDF_ ->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else if (nom=="barycentre_IBM")
    {
      if (champ_barycentre_IBM_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ_barycentre_IBM_
      DoubleTab& valeurs = champ_barycentre_IBM_->valeurs();
      valeurs=0.;
      if (champ_barycentre_.non_nul())
        {
          const DoubleTab& barycentre = champ_barycentre_->valeurs();
          int nb_d0 = barycentre.dimension(0);
          int nb_d1 = barycentre.dimension(1);
          for (int num0=0; num0<nb_d0; num0++)
            for (int num1=0; num1<nb_d1; num1++)
              valeurs(num0,num1)=barycentre(num0,num1);
        }
      valeurs.echange_espace_virtuel();
      champ_barycentre_IBM_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else if (nom=="normal_IBM")
    {
      if (champ_normal_IBM_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ_normal_IBM_
      DoubleTab& valeurs = champ_normal_IBM_->valeurs();
      valeurs=0.;
      if (champ_rotation_.non_nul())
        {
          const DoubleTab& rotation = champ_rotation_->valeurs();
          int nb_d0 = rotation.dimension(0);
          for (int num0=0; num0<nb_d0; num0++)
            for (int num1=0; num1<dim_esp; num1++)
              valeurs(num0,num1)=rotation(num0,3*num1+(dim_esp-1));
        }
      valeurs.echange_espace_virtuel();
      champ_normal_IBM_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else if (nom=="aire_IBM")
    {
      if (champ_aire_IBM_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ_aire_IBM_
      DoubleTab& valeurs = champ_aire_IBM_->valeurs();
      valeurs=0.;
      if (champ_aire_.non_nul())
        {
          const DoubleTab& aire = champ_aire_->valeurs();
          int nb_d0 = aire.dimension(0);
          int nb_d1 = aire.dimension(1);
          for (int num0=0; num0<nb_d0; num0++)
            for (int num1=0; num1<nb_d1; num1++)
              valeurs(num0,num1)=aire(num0,num1);
        }
      valeurs.echange_espace_virtuel();
      champ_aire_IBM_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else if (nom=="velocity_shape_IBM")
    {
      if (champ_vitesse_shape_IBM_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ_vitesse_shape_IBM_
      DoubleTab& valeurs = champ_vitesse_shape_IBM_->valeurs();
      valeurs=0.;
      if (modele_lu_.vitesse_shape_IBM_.non_nul())
        {
          const DoubleTab& vitesse = modele_lu_.get_vitesse_shape_IBM();
          int nb_d0 = vitesse.dimension(0);
          int nb_d1 = vitesse.dimension(1);
          for (int num0=0; num0<nb_d0; num0++)
            for (int num1=0; num1<nb_d1; num1++)
              valeurs(num0,num1)=vitesse(num0,num1);
        }
      valeurs.echange_espace_virtuel();
      champ_vitesse_shape_IBM_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else if (nom=="pseudo_level_set_IBM")
    {
      if (champ_pseudo_level_set_IBM_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ_pseudo_level_set_IBM_
      DoubleTab& valeurs = champ_pseudo_level_set_IBM_->valeurs();
      valeurs=0.;
      if (interpolation_bool_)
        {
          if (interpolation_lue_->pseudo_level_set_.non_nul())
            {
              Interpolation_IBM_base& my_interp = ref_cast_non_const(Interpolation_IBM_base, getInterpolationLu());
              const DoubleTab& pseudo_level_set = my_interp.get_pseudo_level_set();
              int nb_d0 = pseudo_level_set.dimension(0);
              for (int num0=0; num0<nb_d0; num0++) valeurs(num0)=pseudo_level_set(num0);
            }
        }
      valeurs.echange_espace_virtuel();
      champ_pseudo_level_set_IBM_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else
    return champs_compris_.get_champ(nom);
}
