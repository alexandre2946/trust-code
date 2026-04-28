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

#include <Interpolation_IBM_base.h>
#include <Source_PDF_base.h>

Implemente_base(Interpolation_IBM_base, "Interpolation_IBM_base", Objet_U);
// XD interpolation_ibm_base objet_u interpolation_ibm_base 0 Base class for all the interpolation methods available in the Immersed Boundary Method (IBM).

void Interpolation_IBM_base::discretise(const Discretisation_base& dis, Domaine_dis_base& le_dom_dis)
{
  int nb_comp = Objet_U::dimension;
  Noms units(nb_comp);
  Noms c_nam(nb_comp);

  dis.discretiser_champ("champ_sommets",le_dom_dis,vectoriel,c_nam,units,nb_comp,0.,solid_points_);
  if (solid_points_from_prepro_)
    {
      OBS_PTR(Prepro_IBM_base) my_prep = my_source_->getpreproLu();
      if ((&my_prep)->non_nul())
        {
          DoubleTab& the_values = ref_cast_non_const(DoubleTab, my_prep->get_champ_solid_points());
          solid_points_->valeurs() = the_values;
        }
    }
  else
    {
      if (solid_points_lu_.non_nul()) solid_points_->affecter(solid_points_lu_);
    }

  if (is_dirichlet_from_prepro_)
    {
      OBS_PTR(Prepro_IBM_base) my_prep = my_source_->getpreproLu();
      if ((&my_prep)->non_nul())
        {
          DoubleTab& the_values = ref_cast_non_const(DoubleTab, my_prep->get_isNodeDirichlet());
          dis.discretiser_champ("champ_sommets",le_dom_dis,"is_dirichlet","none",1,0., is_dirichlet_);
          is_dirichlet_->valeurs() = the_values;
        }
    }
  else
    {
      if(is_dirichlet_lu_.non_nul())
        {
          dis.discretiser_champ("champ_sommets",le_dom_dis,"is_dirichlet","none",1,0., is_dirichlet_);
          is_dirichlet_->affecter(is_dirichlet_lu_);
        }
      else
        {
          if (my_source_->equation().discretisation().is_vef() && my_source_->get_imm_wall_law())
            {
              Cerr<<"Interpolation_IBM_base::discretise field est_dirichlet not provided for imm. wall law with VEF discretization. Exit."<<finl;
              exit();
            }

        }
    }

  if (corresp_elems_from_prepro_)
    {
      OBS_PTR(Prepro_IBM_base) my_prep = my_source_->getpreproLu();
      if ((&my_prep)->non_nul())
        {
          DoubleTab& the_values = ref_cast_non_const(DoubleTab, my_prep->get_champ_corresp_elems());
          dis.discretiser_champ("champ_elem",le_dom_dis,"corresp_elems","none",1,0., corresp_elems_);
          corresp_elems_->valeurs() = the_values;
        }
    }
  else
    {
      if (corresp_elems_lu_.non_nul())
        {
          dis.discretiser_champ("champ_elem",le_dom_dis,"corresp_elems","none",1,0., corresp_elems_);
          corresp_elems_->affecter(corresp_elems_lu_);
        }
    }

  if (my_source_->get_modele().get_PDF_mobile()) discretise_PDF_mobile(dis, le_dom_dis);

  return;
}

void Interpolation_IBM_base::discretise_PDF_mobile(const Discretisation_base& dis, Domaine_dis_base& le_dom_dis)
{
  int nb_comp = Objet_U::dimension;

  dis.discretiser_champ("champ_sommets",le_dom_dis,"distance_projete_solide","",1,0., champ_dis_proj_solid_);
  DoubleTab& distPSArray = champ_dis_proj_solid_->valeurs();
  distPSArray = 0.;

  dis.discretiser_champ("champ_sommets",le_dom_dis,"normal_projete_solide","",nb_comp,0., champ_normal_proj_solid_);
  DoubleTab& normPSArray = champ_normal_proj_solid_->valeurs();
  normPSArray = 0.;

  calculer_normal_et_distance_proj_solid();

  if (my_source_->get_modele().get_use_pseudo_level_set_moving_PDF())
    {
      discretise_pseudo_level_set(dis, le_dom_dis);
      definir_pseudo_level_set();
    }

  return;
}

void Interpolation_IBM_base::discretise_pseudo_level_set(const Discretisation_base& dis, Domaine_dis_base& le_dom_dis)
{
  dis.discretiser_champ("champ_sommets",le_dom_dis,"pseudo_level_set","",1,0., pseudo_level_set_);
  DoubleTab& PSLevelStSArray = pseudo_level_set_->valeurs();
  PSLevelStSArray = 0.;
}

Sortie& Interpolation_IBM_base::printOn( Sortie& os ) const
{
  return os;
}

Entree& Interpolation_IBM_base::readOn( Entree& is )
{
  return is;
}

void Interpolation_IBM_base::set_param(Param& param) const
{
  param.ajouter_flag("impr",&impr_);  // XD_ADD_P flag To print IBM-related data
  param.ajouter("nb_histo_boxes_impr",&N_histo_,Param::OPTIONAL);  // XD_ADD_P entier number of histogram boxes for printed data
  param.ajouter("est_dirichlet",&is_dirichlet_lu_,Param::OPTIONAL);   // XD_ADD_P field_base Node field of booleans indicating whether the node belong to an element where the interface is
  param.ajouter_flag("get_solid_points_from_prepro", &solid_points_from_prepro_); // XD_ADD_P get IBM solid points from prepro.
  param.ajouter_flag("get_is_dirichlet_from_prepro", &is_dirichlet_from_prepro_); // XD_ADD_P get IBM is_dirichlet from prepro.
  param.ajouter_flag("get_corresp_elems_from_prepro", &corresp_elems_from_prepro_); // XD_ADD_P get IBM corresp_elems from prepro.
}

void Interpolation_IBM_base::definir_pseudo_level_set()
{
  DoubleTab& distance_signee = pseudo_level_set_->valeurs();

  int dim_esp = Objet_U::dimension;
  const DoubleTab& aire = my_source_->get_champ_aire().valeurs();
  int nb_elem = aire.dimension(0);
  const Probleme_base& pb = my_source_->equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine& dom = le_dom_dis.domaine();
  int nb_node_tot = dom.nb_som_tot();
  int nb_elem_tot = dom.nb_elem_tot();

  IntLists elem_voisins(nb_elem_tot);
  DoubleTab& aire_ncst = ref_cast_non_const(DoubleTab, aire);
  bool all_elem_vois = true;
  my_source_->compute_NeighNode_IBM_elem(aire_ncst, elem_voisins, all_elem_vois);

  DoubleTab& nor = champ_normal_proj_solid_->valeurs();
  DoubleTab& d1 = champ_dis_proj_solid_->valeurs();
  if (! nor.size_array() )
    {
      Cerr << "Interpolation_IBM_base::definir_pseudo_level_set : pas de champ_normal_proj_solid_. Exit"<< finl;
      exit();
    }
  if (! d1.size_array() )
    {
      Cerr << "Interpolation_IBM_base::definir_pseudo_level_set : pas de champ_dis_proj_solid_. Exit"<< finl;
      exit();
    }

  distance_signee=-99999.;
  IntTrav deja_fait_elem(nb_elem_tot);
  DoubleTrav nor_ref_i(nb_node_tot,dim_esp);
  IntList to_do_cut_elem; // Liste voisin cut cell a traiter
  to_do_cut_elem.vide(); // On vide la List
  int idebug_lset = 0;
  for (int e=0; e<nb_elem; e++)
    {
      if (aire(e)>0. && deja_fait_elem(e) == 0) // Determination d un element racine
        {
          if (idebug_lset) Cerr<<">>>>>>>>>>>>>> element racine = "<<e<<finl;
          define_pseudo_level_set_for_one_cut_cell(elem_voisins, e, deja_fait_elem, nor_ref_i, to_do_cut_elem, idebug_lset);

          // Traitement des elements dans to_do_cut_elem
          if (idebug_lset)
            {
              Cerr<<"to_do_cut_elem after define_pseudo_level_set_for_one_cut_cell : ";
              for (int ip=0; ip<to_do_cut_elem.size(); ip++)  Cerr<<to_do_cut_elem[ip]<<" ";
              Cerr<<finl;
            }
          while (to_do_cut_elem.size() != 0)
            {
              int elem_voi = to_do_cut_elem[0];
              if (idebug_lset) Cerr<<"cutted elem_voi = "<<elem_voi<<finl;
              define_pseudo_level_set_for_one_cut_cell(elem_voisins, elem_voi, deja_fait_elem, nor_ref_i, to_do_cut_elem, idebug_lset);
              to_do_cut_elem.suppr(elem_voi);
              if (idebug_lset)
                {
                  Cerr<<"to_do_cut_elem after define_pseudo_level_set_for_one_cut_cell : ";
                  for (int ip=0; ip<to_do_cut_elem.size(); ip++)  Cerr<<to_do_cut_elem[ip]<<" ";
                  Cerr<<finl;
                  Cerr<<"///////////////////////////////////////"<<finl;
                }
            }
        }
    }
}

void Interpolation_IBM_base::define_pseudo_level_set_for_one_cut_cell(IntLists& elem_voisins, int elem, IntTab& deja_fait_elem, DoubleTab& nor_ref_i, IntList& to_do_cut_elem, int idebug_lset)
{
  DoubleTab& distance_signee = pseudo_level_set_->valeurs();
  DoubleTab& nor = champ_normal_proj_solid_->valeurs();
  DoubleTab& d1 = champ_dis_proj_solid_->valeurs();
  int dim_esp = Objet_U::dimension;
  const Probleme_base& pb = my_source_->equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine& dom = le_dom_dis.domaine();
  const IntTab& elems= dom.les_elems() ;
  int nb_som_elem=dom.nb_som_elem();

  bool exist_ref_e = false;
  DoubleTrav nor_ref_e(dim_esp);
  // definition d une reference pour l element elem
  for (int ilr=0; ilr<nb_som_elem; ilr++)
    {
      int ir = elems(elem,ilr);
      // On cherche si il existe une reference pour un sommet
      if (distance_signee(ir) != -99999.)
        {
          for (int k=0; k<dim_esp; k++) nor_ref_e(k) = nor_ref_i(ir,k);
          exist_ref_e = true;
          if (idebug_lset) Cerr<<" node = "<<ir;
          break;
        }
    }
  if (idebug_lset) Cerr<<" > element = "<<elem<<" exist_ref_e = "<<int(exist_ref_e)<<finl;

  // Distance signee pour les sommets de l element elem
  for (int il=0; il<nb_som_elem; il++)
    {
      int i = elems(elem,il);
      if (distance_signee(i) == -99999.) // Sommet pas traite
        {
          if (!exist_ref_e) // Pas de reference pour l element
            {
              double norm = 0.;
              for (int k=0; k<dim_esp; k++) norm += nor(i,k) * nor(i,k);
              if (norm > 1.e-10)
                {
                  for (int k=0; k<dim_esp; k++) nor_ref_e(k) = nor(i,k);
                  exist_ref_e = true;
                  if (idebug_lset)
                    {
                      Cerr<<"i = "<<i<<" : exist_ref_e nor_ref_e = "<<int(exist_ref_e);
                      for (int k=0; k<dim_esp; k++) Cerr<<" "<<nor_ref_e(k);
                      Cerr<<finl;
                    }
                }
              else if (idebug_lset) Cerr<<"*** i = "<<i<<" norme  = "<<norm<<finl;
            }
          if (exist_ref_e)
            {
              double pscal = 0.;
              for (int k=0; k<dim_esp; k++) pscal += nor(i,k) * nor_ref_e(k);
              if (abs(pscal) > 1.e-10) distance_signee(i) = d1(i) * pscal / abs(pscal);
              for (int k=0; k<dim_esp; k++) nor_ref_i(i,k) = nor(i,k)*distance_signee(i)/abs(distance_signee(i));
              if (idebug_lset)
                {
                  Cerr<<"*** i = "<<i<<" nor = ";
                  for (int k=0; k<dim_esp; k++) Cerr<<" "<<nor(i,k);
                  Cerr<<" distance_signee = "<<distance_signee(i)<<finl;
                }
            }
          else
            {
              Cerr<<"Interpolation_IBM_base::definir_pseudo_level_set: pas de reference pour l element "<<elem<<finl;
              exit();
            }
        }
    }
  deja_fait_elem(elem) = 1;

  // definition distance signee pour le cluster d'elements non coupes voisins de e
  if (exist_ref_e) calcul_cluster_pseudo_level_set(elem_voisins, elem, deja_fait_elem, nor_ref_i, to_do_cut_elem, idebug_lset);
}

void Interpolation_IBM_base::calcul_cluster_pseudo_level_set(IntLists& elem_voisins, int elem, IntTab& deja_fait_elem, DoubleTab& nor_ref_i, IntList& to_do_cut_elem, int idebug_lset)
{
  int dim_esp = Objet_U::dimension;
  const DoubleTab& aire = my_source_->get_champ_aire().valeurs();
  const Probleme_base& pb = my_source_->equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine& dom = le_dom_dis.domaine();
  const IntTab& elems= dom.les_elems() ;
  int nb_som_elem=dom.nb_som_elem();

  DoubleTab& nor = champ_normal_proj_solid_->valeurs();
  DoubleTab& d1 = champ_dis_proj_solid_->valeurs();
  DoubleTab& distance_signee = pseudo_level_set_->valeurs();
  if (! nor.size_array() )
    {
      Cerr << "Interpolation_IBM_base::calcul_cluster_pseudo_level_set : pas de champ_normal_proj_solid_. Exit"<< finl;
      exit();
    }
  if (! d1.size_array() )
    {
      Cerr << "Interpolation_IBM_base::calcul_cluster_pseudo_level_set : pas de champ_dis_proj_solid_. Exit"<< finl;
      exit();
    }

  int nb_elem_voi = elem_voisins[elem].size();
  if (nb_elem_voi != 0)
    {
      for (int voi=0; voi<nb_elem_voi; voi++)
        {
          int elem_voi = (elem_voisins[elem])[voi];
          if (idebug_lset) Cerr<<">>>> element of cluster = "<<elem_voi<<" aire = "<<aire(elem_voi)<<" deja_fait_ele = "<<deja_fait_elem(elem_voi)<<finl;

          if (elem_voi != elem && aire(elem_voi) > 0. && deja_fait_elem(elem_voi) == 0) to_do_cut_elem.add_if_not(elem_voi);

          bool exist_ref_e = false;
          DoubleTrav nor_ref_e(dim_esp);
          if (elem_voi != elem && deja_fait_elem(elem_voi) == 0 && aire(elem_voi) <= 0.)
            {
              // definition d une reference pour l element e
              for (int ilr=0; ilr<nb_som_elem; ilr++)
                {
                  int ir = elems(elem_voi,ilr);
                  // On cherche si il existe une reference pour un sommet
                  if (distance_signee(ir) != -99999.)
                    {
                      for (int k=0; k<dim_esp; k++) nor_ref_e(k) = nor_ref_i(ir,k);
                      exist_ref_e = true;
                      if (idebug_lset) Cerr<<" node = "<<ir;
                      break;
                    }
                }
              if (idebug_lset)
                {
                  Cerr<<" > exist_ref_e = "<<int(exist_ref_e)<<" nor_ref_e =";
                  for (int k=0; k<dim_esp; k++) Cerr<<" "<<nor_ref_e(k);
                  Cerr<<finl;
                }

              if (!exist_ref_e) // Pas de reference pour l element
                {
                  Cerr<<"Interpolation_IBM_base::definir_pseudo_level_set: pas de reference pour l element "<<elem_voi<<finl;
                  exit();
                }

              for (int il_v=0; il_v<nb_som_elem; il_v++)
                {
                  int i_v = elems(elem_voi,il_v);
                  if (distance_signee(i_v) == -99999.)
                    {
                      double pscal = 0.;
                      for (int k=0; k<dim_esp; k++) pscal += nor(i_v,k) * nor_ref_e(k);
                      if (abs(pscal) > 1.e-10) distance_signee(i_v) = d1(i_v) * pscal / abs(pscal);
                      for (int k=0; k<dim_esp; k++) nor_ref_i(i_v,k) = nor(i_v,k)*distance_signee(i_v)/abs(distance_signee(i_v));
                      if (idebug_lset)
                        {
                          Cerr<<"*** i_v = "<<i_v<<" nor = ";
                          for (int k=0; k<dim_esp; k++) Cerr<<" "<<nor(i_v,k);
                          Cerr<<" distance_signee = "<<distance_signee(i_v)<<finl;
                        }
                    }
                }
              deja_fait_elem(elem_voi) = 1;
            }
        }
    }
}

void Interpolation_IBM_base::calculer_normal_et_distance_proj_solid()
{
  // 3 layers
  const DoubleTab& solid_points = solid_points_->valeurs();
  if (! solid_points.size_array() )
    {
      Cerr << "Interpolation_IBM_base::maj_normal_proj_solid : pas de projection solid. Exit"<< finl;
      exit();
    }
  double eps = 1e-12;
  int dim_esp = Objet_U::dimension;
  assert (solid_points.dimension(1) == dim_esp);
  const Probleme_base& pb = my_source_->equation().probleme();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  const Domaine& dom = le_dom_dis.domaine();
  const IntTab& elems= dom.les_elems() ;
  int nb_som_elem=dom.nb_som_elem();
  int nb_elem_tot = dom.nb_elem_tot();
  int nb_node_tot = dom.nb_som_tot();

  const DoubleTab& aire = my_source_->get_champ_aire().valeurs();
  int nb_elem = aire.dimension(0);

  IntLists elem_voisins(nb_elem_tot);
  DoubleTab& aire_ncst = ref_cast_non_const(DoubleTab, aire);
  bool all_elem_vois = true;
  my_source_->compute_NeighNode_IBM_elem(aire_ncst, elem_voisins, all_elem_vois);

  DoubleTab& nor = champ_normal_proj_solid_->valeurs();
  nor = 0.;
  DoubleTab& d1 = champ_dis_proj_solid_->valeurs();
  IntTrav deja_fait_vert(nb_node_tot);
  d1 = 0.;
  for (int e=0; e<nb_elem; e++)
    {
      if (aire(e)>0.)
        {
          // element coupe IB
          for (int il=0; il<nb_som_elem; il++)
            {
              int i = elems(e,il);
              if (deja_fait_vert(i) == 0)
                {
                  d1(i) = 0.0;
                  for (int k=0; k<dim_esp; k++)
                    {
                      double xk = dom.coord(i,k);
                      double xks = solid_points(i,k);
                      nor(i, k) = xk - xks;
                      d1(i) += (xk-xks)*(xk-xks);
                    }
                  d1(i) = sqrt(d1(i));
                  if (d1(i)  > eps )
                    for (int k=0; k<dim_esp; k++) nor(i, k) /= d1(i) ;
                  deja_fait_vert(i) = 1;
                }
            }
          //voisins de e
          int nb_elem_voi = elem_voisins[e].size();
          for (int voi=0; voi<nb_elem_voi; voi++)
            {
              int elem_voi = (elem_voisins[e])[voi];
              for (int il_v=0; il_v<nb_som_elem; il_v++)
                {
                  int i_v = elems(elem_voi,il_v);
                  if (deja_fait_vert(i_v) == 0)
                    {
                      d1(i_v) = 0.0;
                      for (int k=0; k<dim_esp; k++)
                        {
                          double xk = dom.coord(i_v,k);
                          double xks = solid_points(i_v,k);
                          nor(i_v, k) = xk - xks;
                          d1(i_v) += (xk-xks)*(xk-xks);
                        }
                      d1(i_v) = sqrt(d1(i_v));
                      if (d1(i_v)  > eps )
                        for (int k=0; k<dim_esp; k++) nor(i_v, k) /= d1(i_v) ;
                      deja_fait_vert(i_v) = 1;
                    }
                }
            }
        }
    }
}

void Interpolation_IBM_base::set_fields_from_prepro_to_interp(Prepro_IBM_base& un_prepro)
{
  solid_points_->valeurs() = un_prepro.get_champ_solid_points();
  corresp_elems_->valeurs() = un_prepro.get_champ_corresp_elems();
  if(is_dirichlet_.non_nul())
    {
      is_dirichlet_->valeurs() = un_prepro.get_isNodeDirichlet();
    }
  else
    {
      Cerr<<"Interpolation_IBM_base::set_fields_from_prepro_to_interp field est_dirichlet not provided. Exit."<<finl;
      exit();
    }

  if (my_source_->get_modele().get_PDF_mobile()) calculer_normal_et_distance_proj_solid();
}
