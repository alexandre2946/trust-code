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

#include <Champ_Generique_Reduction_0D.h>
#include <Discretisation_base.h>
#include <TRUSTTab_parts.h>
#include <communications.h>
#include <Probleme_base.h>
#include <Synonyme_info.h>

#include <Domaine_VF.h>
#include <Param.h>

Implemente_instanciable(Champ_Generique_Reduction_0D,"Reduction_0D",Champ_Gen_de_Champs_Gen);
Add_synonym(Champ_Generique_Reduction_0D,"Champ_Post_Reduction_0D");
// XD reduction_0d champ_post_de_champs_post reduction_0d INHERITS_BRACE To calculate the min, max, sum, average,
// XD_CONT weighted sum, weighted average, weighted sum by porosity, weighted average by porosity, euclidian norm,
// XD_CONT normalized euclidian norm, L1 norm, L2 norm of a field.

Sortie& Champ_Generique_Reduction_0D::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

//see Champ_Gen_de_Champs_Gen::readOn
Entree& Champ_Generique_Reduction_0D::readOn(Entree& s )
{
  LIST(Motcle) mot_compris;
  mot_compris.add("min");
  mot_compris.add("max");
  mot_compris.add("euclidian_norm"); // new name for norme_L2
  mot_compris.add("normalized_euclidian_norm"); // new name for normalized_norm_L2
  mot_compris.add("moyenne");
  mot_compris.add("somme");
  mot_compris.add("moyenne_ponderee");
  mot_compris.add("somme_ponderee");
  mot_compris.add("moyenne_ponderee_porosite");
  mot_compris.add("somme_ponderee_porosite");
  mot_compris.add("valeur_a_gauche");
  mot_compris.add("L2_norm");  // L2 norm
  mot_compris.add("L1_norm");  // L1 norm
  mot_compris.add("average");  // new name for moyenne
  mot_compris.add("sum");  // new name for somme
  mot_compris.add("weighted_average"); // new name for moyenne_ponderee
  mot_compris.add("weighted_sum"); // new name weighted_sum for somme_ponderee
  mot_compris.add("weighted_average_porosity");  // new name for moyenne_ponderee_porosite
  mot_compris.add("weighted_sum_porosity");  // new name for somme_ponderee_porosite
  mot_compris.add("left_value");  // new name for valeur_a_gauche

  Champ_Gen_de_Champs_Gen::readOn(s);
  if (mot_compris.rang(methode_)<0)
    {
      Cerr << "Method " << methode_ << " is an unknown option for methode keyword in "<< que_suis_je() << "." << finl;
      Cerr << "Choose from " << mot_compris << finl;
      Process::exit();
    }
  return s ;
}

//  methode : indicates the type of reduction to be performed
//              (min, max, moyenne, moyenne_ponderee_volume_elem, somme, somme_ponderee)
void Champ_Generique_Reduction_0D::set_param(Param& param) const
{
  Champ_Gen_de_Champs_Gen::set_param(param);
  param.ajouter("methode",&methode_,Param::REQUIRED); // XD_ADD_P chaine(into=["min","max","moyenne","average","moyenne_ponderee","weighted_average","somme","sum","somme_ponderee","weighted_sum","somme_ponderee_porosite","weighted_sum_porosity","euclidian_norm","normalized_euclidian_norm","L1_norm","L2_norm","valeur_a_gauche","left_value"])
  // XD_CONT name of the reduction method: NL2 - min for the minimum value, NL2 - max for the maximum value, NL2 -
  // XD_CONT average (or moyenne) for a mean, NL2 - weighted_average (or moyenne_ponderee) for a mean ponderated by
  // XD_CONT integration volumes, e.g: cell volumes for temperature and pressure in VDF, volumes around faces for
  // XD_CONT velocity and temperature in VEF, NL2 - sum (or somme) for the sum of all the values of the field, NL2 -
  // XD_CONT weighted_sum (or somme_ponderee) for a weighted sum (integral), NL2 - weighted_average_porosity (or
  // XD_CONT moyenne_ponderee_porosite) and weighted_sum_porosity (or somme_ponderee_porosite) for the mean and sum
  // XD_CONT weighted by the volumes of the elements, only for ELEM localisation, NL2 - euclidian_norm for the euclidian
  // XD_CONT norm, NL2 - normalized_euclidian_norm for the euclidian norm normalized, NL2 - L1_norm for norm L1, NL2 -
  // XD_CONT L2_norm for norm L2
}

void Champ_Generique_Reduction_0D::completer(const Postraitement_base& post)
{
  Champ_Gen_de_Champs_Gen::completer(post);
  //Champ source_espace_stockage;
  Motcle directive;
  directive = get_directive_pour_discr();

  if (directive=="pression")
    {
      const Noms& nom_champ = get_property("nom");
      const Noms& nom_source = get_source(0).get_property("nom");
      Cerr<<"Problem with the "<<que_suis_je()<<" post processing field named "<<nom_champ[0]<<"."<<finl;
      Cerr<<"The source field named "<<nom_source[0]<<" of this last "<<que_suis_je()<<" field "<<finl;
      Cerr<<"must be previously interpolated at the elem or som location."<<finl;
      Cerr<<"Please use instead the syntax : "<<finl;
      Cerr<<"..."<<que_suis_je()<<" { source Interpolation { localisation ... } ... }"<<finl;
      Cerr<<"or contact TRUST support."<<finl;
      exit();
    }
  if (methode_=="valeur_a_gauche" || methode_=="left_value")
    {
      numero_proc_=-1;
    }

  Journal()<<"METHODE "<<methode_<<finl;
  if ((methode_=="valeur_a_gauche" || methode_=="left_value")&&(numero_proc_==-1))
    {
      const Domaine_dis_base& domaine_dis = get_source(0).get_ref_domaine_dis_base();

      const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
      // leftmost position
      const DoubleTab& coords=zvf.domaine().les_sommets();
      const IntTab& conn=zvf.domaine().les_elems();
      double minp=mp_min_vect(coords);
      //Cerr<<" uu "<<minp<<finl;
      minp=Process::mp_min(minp);
      // Cerr<<" uu "<<minp<<finl;
      const Domaine& domaine =zvf.domaine();
      int elemin=-1;
      double dmin=DMAXFLOAT;
      int nb_elem=domaine.nb_elem();
      //ArrOfDouble xp(dimension);
      for (int ele=0; ele<nb_elem; ele++)
        {

          double d=0;
          for (int dir=0; dir<dimension; dir++)
            {
              double xp=0;
              int nb_som_elem=conn.dimension(1);
              for (int s=0; s<nb_som_elem; s++)
                {
                  int som=conn(ele,s);
                  xp+=coords(som,dir);
                }
              xp/=nb_som_elem;
              d+=(xp-minp)*(xp-minp);
            }
          if (d<dmin)
            {
              dmin=d;
              elemin=ele;
            }
        }
      double dming=mp_min(dmin);
      // we send to the master
      if (!je_suis_maitre())
        {
          envoyer(dmin,0,97);
          envoyer(elemin,0,97);
        }
      if (je_suis_maitre())
        {
          numero_proc_=me();
          numero_elem_=elemin;
          for (int p=1; p<nproc(); p++)
            {
              double dminloc;
              int eleminloc;
              recevoir(dminloc,p,97);
              recevoir(eleminloc,p,97);

              if (dminloc<dmin)
                {
                  dmin=dminloc;
                  numero_proc_=p;
                  numero_elem_=eleminloc;
                }
            }

          if (!est_egal(dmin,dming))
            {
              Cerr<<" iiiiiiiiiii "<<dmin<<" "<<dming<<finl;
              exit();
            }
        }
      envoyer_broadcast(numero_proc_,0);
      envoyer_broadcast(numero_elem_,0);

    }

}

const Champ_base& Champ_Generique_Reduction_0D::get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const
{

  OWN_PTR(Champ_base) source_espace_stockage_tmp;
  const Champ_base& source = get_source(0).get_champ_without_evaluation(source_espace_stockage_tmp);
  Nature_du_champ nature_source = source.nature_du_champ();
  int nb_comp = source.nb_comp();

  OWN_PTR(Champ_Fonc_base)  es_tmp;
  espace_stockage = creer_espace_stockage(nature_source,nb_comp,es_tmp);
  return espace_stockage;
}
/*! @brief 0D reduction of the source field (in the sense that we make it uniform) according to the method (min, max, moyenne, moyenne_ponderee_volume_elem, somme, somme_ponderee)
 *
 *  In the case where the field has multiple components, they are processed one by one
 *
 */
const Champ_base& Champ_Generique_Reduction_0D::get_champ(OWN_PTR(Champ_base)&) const
{
  const Champ_base& source = get_source(0).get_champ(source_espace_stockage_);
  const Domaine_dis_base& domaine_dis = get_source(0).get_ref_domaine_dis_base();
  Nature_du_champ nature_source = source.nature_du_champ();
  bool basis_function = source.is_basis_function();
  int order = source.order_field();
  int nb_comp = source.nb_vect_comp();

  // dimension() on the value array of PolyMAC_HFV fields returns -1 (multiple supports)
  // ToDo: completely rewrite this method (horrible, very poorly written) by delegating the min/max/sum/... methods for each OWN_PTR(Champ_base)!
  if (source.que_suis_je()=="Champ_Face_PolyMAC_HFV" || source.que_suis_je()=="Champ_Face_PolyMAC_MPFA")
    Process::exit("PolyMAC_HFV/PolyMAC_MPFA face field not supported yet for Reduction_0D");


  if (!espace_stockage_)
    creer_espace_stockage(nature_source,nb_comp,espace_stockage_);
  else
    espace_stockage_->changer_temps(get_time());

  int nb_dim = source.valeurs().nb_dim();
  // correction for 3D
  if (nb_dim==2)
    nb_dim= source.valeurs().dimension(1);

  ConstDoubleTab_parts valeurs_source_parts(source.valeurs()); // to ignore auxiliary variables
  const DoubleTab& valeurs_source = valeurs_source_parts[0];   // of PolyMAC_HFV (otherwise: min, moyenne WRONG)
  DoubleTab& espace_valeurs = espace_stockage_->valeurs();
  const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
  double val_extraite=-100.;

  if (source.is_basis_function() or source.is_quadrature())
    {
      if (nb_comp==1)
        {
          extraire(val_extraite,valeurs_source,basis_function,order);
          espace_valeurs = val_extraite;
        }
      else
        {
          assert(nb_comp==Objet_U::dimension);
          int size_vect = valeurs_source.dimension(0);
          DoubleTrav vect_source;
          vect_source.resize(size_vect);
          for (int comp=0; comp<nb_comp; comp++)
            {
              {
                ToDo_Kokkos("Code but check test!");
                CDoubleTabView valeurs = valeurs_source.view_ro();
                DoubleArrView vect = static_cast<ArrOfDouble&>(vect_source).view_wo();
                Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), size_vect, KOKKOS_LAMBDA(const int i)
                {
                  vect(i) = valeurs(i,comp);
                });
                end_gpu_timer(__KERNEL_NAME__);
              }
              extraire(val_extraite,vect_source,basis_function,order);
              {
                DoubleTabView valeurs = espace_valeurs.view_wo();
                Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), size_vect, KOKKOS_LAMBDA(const int i)
                {
                  valeurs(i,comp) = val_extraite;
                });
                end_gpu_timer(__KERNEL_NAME__);
              }
            }
        }
    }
  else if (nb_comp==1)
    {
      extraire(val_extraite,valeurs_source,basis_function);
      espace_valeurs = val_extraite;
    }
  else
    {
      for (int comp=0; comp<nb_comp; comp++)
        {
          int size_vect = valeurs_source.dimension(0);
          if (nb_dim!=nb_comp) //Cas des Champ_Face_VDF
            {
              size_vect=0;
              ToDo_Kokkos("critical, warning check you have a NR test case with .son !");
              for (int i=0; i<valeurs_source.dimension(0); i++)
                if (zvf.orientation(i)==comp)
                  ++size_vect;
            }

          DoubleTrav vect_source;
          //For the somme option, vect_source must have a parallel structure
          //to apply val_extraite = mp_prodscal(vect_source,un)
          //Its dimension is then set relative to the number of items in the source
          //ex: zvf.nb_faces() if loc==FACE
          if (methode_=="somme" || methode_=="moyenne" || methode_=="sum" || methode_=="average")
            {
              Entity loc = get_localisation();
              if (loc==Entity::ELEMENT)
                zvf.domaine().creer_tableau_elements(vect_source,RESIZE_OPTIONS::NOCOPY_NOINIT);
              else if (loc==Entity::NODE)
                zvf.domaine().creer_tableau_sommets(vect_source,RESIZE_OPTIONS::NOCOPY_NOINIT);
              else if (loc==Entity::FACE)
                zvf.creer_tableau_faces(vect_source,RESIZE_OPTIONS::NOCOPY_NOINIT);
              vect_source = 0.;
            }
          else
            //In the VDF case with nb_dim!=nb_comp
            //Its dimension is set relative to the nb_faces whose orientation is comp
            vect_source.resize(size_vect);

          // Filling
          if (nb_dim==nb_comp)
            {
              CDoubleTabView valeurs = valeurs_source.view_ro();
              DoubleArrView vect = static_cast<ArrOfDouble&>(vect_source).view_wo();
              Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), size_vect, KOKKOS_LAMBDA(const int i)
              {
                vect(i) = valeurs(i,comp);
              });
              end_gpu_timer(__KERNEL_NAME__);
            }
          else
            {
              ToDo_Kokkos("critical, warning check you have a NR test case with .son !");
              int k=0;
              for (int i=0; i<valeurs_source.dimension(0); i++)
                if (zvf.orientation(i)==comp)
                  {
                    if (methode_=="somme" || methode_=="moyenne" || methode_=="sum" || methode_=="average")
                      vect_source(i) = valeurs_source(i);
                    else
                      {
                        vect_source(k) = valeurs_source(i);
                        k++;
                      }
                  }
            }
          // Pass the component if necessary for Champ_face fields
          extraire(val_extraite,vect_source,basis_function,(nb_dim==nb_comp?-1:comp));

          if (nb_dim==nb_comp)
            {
              DoubleTabView valeurs = espace_valeurs.view_wo();
              Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), size_vect, KOKKOS_LAMBDA(const int i)
              {
                valeurs(i,comp) = val_extraite;
              });
              end_gpu_timer(__KERNEL_NAME__);
            }
          else
            {
              ToDo_Kokkos("critical, warning check you have a NR test case with .son !");
              for (int i=0; i<valeurs_source.dimension(0); i++)
                if (zvf.orientation(i)==comp)
                  espace_valeurs(i) = val_extraite;
            }
        }
    }
  espace_valeurs.echange_espace_virtuel();

  return espace_stockage_;
}

//Extracts the value from the vector val_source into val_extraite
void Champ_Generique_Reduction_0D::extraire(double& val_extraite,const DoubleVect& val_source, const bool basis_function, const int composante_VDF) const
{

  // Careful :: the composante_VDF for DG gives the order of the source basis_functions

  if (methode_=="min")
    {
      val_extraite = mp_min_vect(val_source);
    }
  else if (methode_=="max")
    {
      val_extraite = mp_max_vect(val_source);
    }
  /*
    else if (methode_=="moyenne") {
    // We do not use mp_moyenne_vect because of problems in VDF with a vector field to make a distributed val_source
  // To calculate the mean, we use somme(val_source)/somme(1)
  val_extraite = mp_moyenne_vect(val_source);
  } */
  else if (methode_=="euclidian_norm")
    {
      val_extraite = mp_norme_vect(val_source);
    }
  else if (methode_=="normalized_euclidian_norm")
    {
      DoubleVect val_un(val_source);
      val_un=1.;
      val_extraite = mp_norme_vect(val_source)/mp_norme_vect(val_un);
    }
  else if (methode_ =="L1_norm" || methode_ =="L2_norm")
    {
      // If we are:
      // - at ELEM -> we weight by element volumes,
      // - at FACE -> we weight by interleaved volumes (we do not account for extended volumes since they are not accessible),
      // - at NODE -> we weight by nodal control volumes [Vol(som)= Sum_over_elem_surrounding_som(Vol_elem/nb_som_per_elem)].
      const Domaine_dis_base& domaine_dis = get_ref_domaine_dis_base();
      const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
      double sum=0;
      const DoubleVect& volumes = zvf.volumes();
      //int volumes_size_tot = mp_sum(volumes.size_array());
      if (volumes.size_array()<zvf.nb_elem())
        {
          Cerr << "The mesh volumes of the domain " << zvf.domaine().le_nom() << " are not available yet." << finl;
          Cerr << "It is not implemented yet." << finl;
          exit();
        }

      // at ELEM
      if (get_localisation()==Entity::ELEMENT)
        {
          if (methode_ =="L1_norm")
            {
              sum = zvf.compute_L1_norm(val_source,basis_function,composante_VDF);
            }
          else if (methode_ =="L2_norm")
            {
              sum = zvf.compute_L2_norm(val_source,basis_function,composante_VDF);
            }
          else
            {
              Cerr << "Error in Champ_Generique_Reduction_0D::extraire" << finl;
              exit();
            }
        }

      // at FACE
      if (get_localisation()==Entity::FACE)
        {
          // Computation of control volumes at each face
          int nb_face = zvf.nb_faces();
          if (!volume_controle_.size())
            {
              volume_controle_.resize(nb_face);
              volume_controle_=0;
              int nb_faces_par_elem = zvf.elem_faces().dimension_tot(1);
              int nb_elem = zvf.nb_elem();
              ToDo_Kokkos("Code but check test!");
              for (int i=0; i<nb_elem; i++)
                for (int j=0; j<nb_faces_par_elem; j++)
                  {
                    int face=zvf.elem_faces(i,j);
                    volume_controle_(face)+=volumes(i)/nb_faces_par_elem;
                  }
            }
          if (composante_VDF>=0)
            {
              const IntVect& ori = zvf.orientation();
              int k=0;
              if (methode_ =="L1_norm")
                {
                  ToDo_Kokkos("Code but check test!");
                  for (int i=0; i<nb_face; i++)
                    {
                      if (ori(i)==composante_VDF)
                        {
                          sum+=std::fabs(val_source(k))*volume_controle_(i);
                          k++;
                        }
                    }
                }
              else if (methode_ =="L2_norm")
                {
                  ToDo_Kokkos("Code but check test!");
                  for (int i=0; i<nb_face; i++)
                    {
                      if (ori(i)==composante_VDF)
                        {
                          sum+=val_source(k)*val_source(k)*volume_controle_(i);
                          k++;
                        }
                    }
                }
              else
                {
                  Cerr << "Error in Champ_Generique_Reduction_0D::extraire" << finl;
                  exit();
                }
            }
          else
            {
              if (methode_ =="L1_norm")
                {
                  ToDo_Kokkos("Code but check test!");
                  for (int i=0; i<nb_face; i++)
                    {
                      sum+=std::fabs(val_source(i))*volume_controle_(i);
                    }
                }
              else if (methode_ =="L2_norm")
                {
                  ToDo_Kokkos("Code but check test!");
                  for (int i=0; i<nb_face; i++)
                    {
                      sum+=val_source(i)*val_source(i)*volume_controle_(i);
                    }
                }
              else
                {
                  Cerr << "Error in Champ_Generique_Reduction_0D::extraire" << finl;
                  exit();
                }
            }
        }

      // at NODE
      if (get_localisation()==Entity::NODE)
        {
          // Computation of control volumes at each node
          int nb_som = zvf.nb_som();
          if (!volume_controle_.size())
            {
              volume_controle_.resize(nb_som);
              volume_controle_=0;
              int nb_som_par_elem = zvf.domaine().les_elems().dimension_tot(1);
              int nb_elem = zvf.nb_elem();
              ToDo_Kokkos("Code but check test!");
              for (int i=0; i<nb_elem; i++)
                for (int j=0; j<nb_som_par_elem; j++)
                  {
                    int som=zvf.domaine().sommet_elem(i,j);
                    volume_controle_(som)+=volumes(i)/nb_som_par_elem;
                  }
            }
          if (methode_ =="L1_norm")
            {
              ToDo_Kokkos("Code but check test!");
              for (int i=0; i<nb_som; i++)
                {
                  sum+=std::fabs(val_source(i))*volume_controle_(i);
                }
            }
          else if (methode_ =="L2_norm")
            {
              ToDo_Kokkos("Code but check test!");
              for (int i=0; i<nb_som; i++)
                {
                  sum+=val_source(i)*val_source(i)*volume_controle_(i);
                }
            }
          else
            {
              Cerr << "Error in Champ_Generique_Reduction_0D::extraire" << finl;
              exit();
            }
        }
      val_extraite = mp_sum(sum);
      if (methode_ =="L2_norm")
        {
          val_extraite = sqrt(val_extraite);
        }
    }

  else if (methode_=="weighted_average" || methode_=="weighted_sum" || methode_=="moyenne_ponderee" || methode_=="somme_ponderee")
    {
      // If we are:
      // - at ELEM -> we weight by element volumes,
      // - at FACE -> we weight by interleaved volumes (we do not account for extended volumes since they are not accessible),
      // - at NODE -> we weight by nodal control volumes [Vol(som)= Sum_over_elem_surrounding_som(Vol_elem/nb_som_per_elem)].

      const Domaine_dis_base& domaine_dis = get_ref_domaine_dis_base();
      const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
      double sum=0;
      double volume=0;
      //int volumes_size_tot = mp_sum(volumes.size_array());
      if (zvf.volumes().size_array()<zvf.nb_elem())
        {
          Cerr << "The mesh volumes of the domain " << zvf.domaine().le_nom() << " are not available yet." << finl;
          Cerr << "It is not implemented yet." << finl;
          exit();
        }

      // at ELEM
      if (get_localisation()==Entity::ELEMENT)
        {
          zvf.compute_average(val_source, sum, volume, basis_function, composante_VDF);
        }

      // at FACE
      else if (get_localisation()==Entity::FACE)
        {
          // Computation of control volumes at each face
          int nb_face = zvf.nb_faces();
          if (!volume_controle_.size() || zvf.domaine().deformable())
            {
              volume_controle_.resize(nb_face);
              volume_controle_=0;
              int nb_faces_par_elem = zvf.elem_faces().dimension_tot(1);
              int nb_elem = zvf.nb_elem();
              CIntTabView elem_faces = zvf.elem_faces().view_ro();
              CDoubleArrView volumes = zvf.volumes().view_ro();
              DoubleArrView volume_controle = volume_controle_.view_rw();
              Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_2D({0, 0}, {nb_elem, nb_faces_par_elem}), KOKKOS_LAMBDA(const int i, const int j)
              {
                int face = elem_faces(i,j);
                Kokkos::atomic_add(&volume_controle(face), volumes(i)/nb_faces_par_elem);
              });
              end_gpu_timer(__KERNEL_NAME__);
            }
          if (composante_VDF>=0)
            {
              ToDo_Kokkos("Code but check test!");
              //const IntVect& ori = zvf.orientation();
              int k=0;
              for (int i=0; i<nb_face; i++)
                //if (ori(i)==composante_VDF)
                if (zvf.orientation(i)==composante_VDF)
                  {
                    sum+=val_source(k)*volume_controle_(i);
                    volume+=volume_controle_(i);
                    k++;
                  }
            }
          else
            {
              CDoubleArrView volume_controle = volume_controle_.view_ro();
              CDoubleArrView val = val_source.view_ro();
              Kokkos::parallel_reduce(start_gpu_timer(__KERNEL_NAME__), nb_face, KOKKOS_LAMBDA(const int i, double & sum_tmp, double & volume_tmp)
              {
                double vc = volume_controle(i);
                sum_tmp += val(i) * vc;
                volume_tmp += vc;
              }, sum, volume);
              end_gpu_timer(__KERNEL_NAME__);
            }
        }

      // at NODE
      else if (get_localisation()==Entity::NODE)
        {
          // Computation of control volumes at each node
          int nb_som = zvf.nb_som();
          if (!volume_controle_.size())
            {
              volume_controle_.resize(nb_som);
              volume_controle_=0;
              const DoubleVect& volumes = zvf.volumes();
              int nb_som_par_elem = zvf.domaine().les_elems().dimension_tot(1);
              int nb_elem = zvf.nb_elem();
              ToDo_Kokkos("Code but check test!");
              for (int i=0; i<nb_elem; i++)
                for (int j=0; j<nb_som_par_elem; j++)
                  {
                    int som=zvf.domaine().sommet_elem(i,j);
                    volume_controle_(som)+=volumes(i)/nb_som_par_elem;
                  }
            }
          ToDo_Kokkos("Code but check test!");
          for (int i=0; i<nb_som; i++)
            {
              sum+=val_source(i)*volume_controle_(i);
              volume+=volume_controle_(i);
            }
        }
      // Optimization: combine 2 mp_sum into 1 collective call
      mp_sum_for_each(sum, volume);
      val_extraite = sum;
      if (methode_=="moyenne_ponderee" || methode_=="weighted_average")
        val_extraite /= volume;
    }
  else if (methode_=="moyenne_ponderee_porosite" || methode_=="somme_ponderee_porosite" || methode_=="weighted_average_porosity" || methode_=="weighted_sum_porosity")
    {
      // - at ELEM -> we weight by element volumes *volumetric_porosity
      // if not at ELEM, error

      const Domaine_dis_base& domaine_dis = get_ref_domaine_dis_base();
      const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
      double sum=0;
      double volume=0;
      const DoubleVect& volumes = zvf.volumes();
      // int volumes_size_tot = mp_sum(volumes.size_array());
      if (volumes.size_array()<zvf.nb_elem())
        {
          Cerr << "The mesh volumes of the domain " << zvf.domaine().le_nom() << " are not available yet." << finl;
          Cerr << "It is not implemented yet." << finl;
          exit();
        }
      if (get_nb_sources()!=2)
        {
          Cerr<<" you must define the porosity "<<finl;
          exit();
        }

      // at ELEM
      if (get_localisation()==Entity::ELEMENT)
        {
          OWN_PTR(Champ_base) source_espace_stockage2;
          const Champ_base& source2 = get_source(1).get_champ(source_espace_stockage2);
          Motcle nom_source_1(get_source(1).get_nom_post());
          if (!(nom_source_1.debute_par("porosite_volumique")||(nom_source_1.debute_par("beta"))))
            {
              Cerr<<" Error with option pondere_porosite !!! "<< get_source(1).get_nom_post()<<finl;
              exit();

            }
          const DoubleVect& poro= source2.valeurs();
          assert(volumes.size_array()==poro.size_array());
          ToDo_Kokkos("Code but check test!");
          zvf.compute_average_porosity(val_source,poro,sum,volume,basis_function,composante_VDF);
        }
      else
        {
          Cerr<<que_suis_je()<<" not implemented for this localisation "<<finl;
          exit();

        }

      // Optimization: combine 2 mp_sum into 1 collective call
      mp_sum_for_each(sum, volume);
      val_extraite = sum;
      if (methode_=="moyenne_ponderee_porosite" || methode_=="weighted_average_porosity")
        val_extraite /= volume;
    }

  else if (methode_=="somme" || methode_=="moyenne" || methode_=="sum" || methode_=="average")
    {
      const Domaine_dis_base& domaine_dis = get_source(0).get_ref_domaine_dis_base();
      const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
      if (!un_.get_md_vector())
        {
          Entity loc = get_localisation();
          if (loc == Entity::ELEMENT)
            zvf.domaine().creer_tableau_elements(un_, RESIZE_OPTIONS::NOCOPY_NOINIT);
          else if (loc == Entity::NODE)
            zvf.domaine().creer_tableau_sommets(un_, RESIZE_OPTIONS::NOCOPY_NOINIT);
          else if (loc == Entity::FACE)
            zvf.creer_tableau_faces(un_, RESIZE_OPTIONS::NOCOPY_NOINIT);
          mapToDevice(un_);
        }
      un_ = 1.;
      if (methode_=="somme" || methode_=="sum")
        {
          val_extraite = mp_prodscal(val_source,un_);
          // Why not use val_extraite = mp_somme_vect(val_source); ?
        }
      else if (methode_=="moyenne" || methode_=="average")
        {
          Entity loc = get_localisation();
          if (loc==Entity::FACE && composante_VDF>=0)
            {
              // In the vector case in VDF, we must only count the
              // faces of the studied component:
              int nb_face = zvf.nb_faces();
              CIntArrView ori = zvf.orientation().view_ro();
              DoubleArrView un = un_.view_wo();
              Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_face, KOKKOS_LAMBDA(const int i)
              {
                if (ori(i)!=composante_VDF)
                  un(i) = 0;
              });
              end_gpu_timer(__KERNEL_NAME__);
            }
          val_extraite = mp_somme_vect(val_source) / mp_somme_vect(un_);
        }
      else
        {
          Cerr << "Error in Champ_Generique_Reduction_0D::extraire" << finl;
          exit();
        }
    }
  else if (methode_=="valeur_a_gauche" || methode_=="left_value")
    {

      assert(numero_proc_!=-1);
      if (me()==numero_proc_)
        {
          val_extraite=val_source(numero_elem_);
          envoyer(val_extraite,me(),-1,98);
        }
      else
        recevoir(val_extraite,numero_proc_,me(),98 );
    }
  else
    {
      Cerr << "Method " << methode_ << " is an unknown option for methode keyword." << finl;
      exit();
    }
}

const Noms Champ_Generique_Reduction_0D::get_property(const Motcle& query) const
{

  Motcles motcles(1);
  motcles[0] = "composantes";
  int rang = motcles.search(query);
  switch(rang)
    {

    case 0:
      {
        Noms source_compos = get_source(0).get_property("composantes");
        int nb_comp = source_compos.size();
        Noms compo(nb_comp);

        for (int i=0; i<nb_comp; i++)
          {
            Nom nume(i);
            compo[i] = nom_post_+nume;
          }

        return compo;
      }

    }
  return Champ_Gen_de_Champs_Gen::get_property(query);
}

//Name the field as a source by default
//"Reduction_0D_"+nom_champ_source
void Champ_Generique_Reduction_0D::nommer_source()
{
  if (nom_post_=="??")
    {
      Nom nom_post_source, nom_champ_source;
      const Noms nom = get_source(0).get_property("nom");
      nom_champ_source = nom[0];
      nom_post_source =  "Reduction_0D_";
      nom_post_source +=  nom_champ_source;
      nommer(nom_post_source);
    }
}

const Motcle Champ_Generique_Reduction_0D::get_directive_pour_discr() const
{
  OWN_PTR(Champ_base) espace_stockage_tmp;
  const Champ_base& source = get_source(0).get_champ_without_evaluation(espace_stockage_tmp);
  if (source.is_basis_function() or source.is_quadrature())
    return "champ_fonc_quad_dg";

  return get_source(0).get_directive_pour_discr();
}
