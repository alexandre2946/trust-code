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

#include <Pb_Couple_Optimisation_IBM.h>
#include <Terme_Derivee_Forme_base.h>
#include <Operateur_Diff_base.h>
#include <Equation_IBM_proto.h>
#include <Prepro_IBM_base.h>
#include <Probleme_base.h>
#include <Operateur.h>
#include <TRUSTTrav.h>

Implemente_instanciable(Pb_Couple_Optimisation_IBM,"Pb_Couple_Optimisation_IBM",Probleme_Couple);

Entree& Pb_Couple_Optimisation_IBM::readOn(Entree& is)
{
  int nb_pb_lu=0;
  Cerr << "<<<<<< Reading of Probleme_Couple_Optimisation_IBM " << le_nom() << " >>>>>>>>>>"<<finl;
  Motcle motlu;
  is >> motlu;
  if (motlu != Motcle("{"))
    {
      Cerr << "We expected { to start to read the Probleme_Couple" << finl;
      exit();
    }
  is >> motlu;
  while (motlu!=Motcle("}"))   // fin du readOn
    {
      if (motlu=="groupes")
        {
          if (nb_problemes())
            {
              Cerr << "We can associate problems to Probleme_Couple" << finl;
              Cerr << "* either by \"associer prob_couple pb\" (in which case they are all in the same group)" << finl;
              Cerr << "* either by the keyword \"groupes\" while reading the object Probleme_Couple" << finl;
              Cerr << "but not both!" << finl;
            }
          assert(nb_problemes()==0);

          LIST(LIST(Nom)) les_noms;
          is >> les_noms;

          groupes.resize_array(les_noms.size());
          for (int i=0; i<les_noms.size(); i++)
            {
              groupes[i]=les_noms[i].size();
              for (int j=0; j<les_noms[i].size(); j++)
                {
                  Nom nom_pb=les_noms[i][j];
                  Objet_U& ob=Interprete::objet(nom_pb);
                  Probleme_base& pb=ref_cast(Probleme_base,ob);
                  ajouter(pb);
                }
            }
        }
      else if (motlu=="state_problem")
        {
          Nom nom_pb;
          is >> nom_pb;
          Objet_U& ob=Interprete::objet(nom_pb);
          Probleme_base& pb_lu=ref_cast(Probleme_base,ob);
          pb_etat_opt_ = pb_lu;
          nb_pb_lu += 1;
          Cerr << "We define state problem as ... " << pb_etat_opt_->le_nom() << finl;
        }

      else if (motlu=="adjoint_problem")
        {
          Nom nom_pb;
          is >> nom_pb;
          Objet_U& ob=Interprete::objet(nom_pb);
          Probleme_base& pb_lu=ref_cast(Probleme_base,ob);
          pb_adjt_opt_ = pb_lu;
          nb_pb_lu += 1;
          Cerr << "We define adjoint problem as ... " << pb_adjt_opt_->le_nom() << finl;
        }

      else if (motlu=="projection_problem")
        {
          Nom nom_pb;
          is >> nom_pb;
          Objet_U& ob=Interprete::objet(nom_pb);
          Probleme_base& pb_lu=ref_cast(Probleme_base,ob);
          pb_projection_opt_ = pb_lu;
          nb_pb_lu += 1;
          Cerr << "We define projection problem as ... " << pb_projection_opt_->le_nom() << finl;
          if ( (pb_projection_opt_->que_suis_je() != "Pb_Conduction") && (pb_projection_opt_->que_suis_je() != "Pb_Conduction_IBM"))
            {
              Cerr << "Projection problem is not of type Pb_Conduction or Pb_Conduction_IBM : "<<pb_projection_opt_->que_suis_je()<<finl;
              abort();
            }
          if ( pb_projection_opt_->nombre_d_equations() != 1 )
            {
              Cerr<<"Pb_Couple_Optimisation_IBM : nb equation != 1 for projection problem."<<finl;
              abort();
            }

        }
      else if (motlu=="objective_function")
        {
          Motcle type;
          is >> type;
          fonction_cout_lu_.typer(type);
          Champ_Don_base& ch_fonction_cout_lu = ref_cast(Champ_Don_base,fonction_cout_lu_.valeur());
          is >> ch_fonction_cout_lu;
          const int nb_comp = ch_fonction_cout_lu.nb_comp();
          if (nb_comp != 1)
            {
              Cerr<<"Pb_Couple_Optimisation_IBM : dimension objective_function != 1 ."<<finl;
              abort();
            }
          fonction_cout_lu_->fixer_nb_comp(nb_comp);
          if (ch_fonction_cout_lu.le_nom()=="anonyme") ch_fonction_cout_lu.nommer("fonction_cout");
          for (int n = 0; n < nb_comp; n++) fonction_cout_lu_->fixer_nom_compo(n, ch_fonction_cout_lu.le_nom() + (nb_comp > 1 ? Nom(n) :""));
        }

      else if (motlu=="objective_function_visualization")
        {
          is >> visu_cout_;
        }

      else if (motlu=="IBM_moving_parameter")
        {
          is >> alpha_;
          Cerr << "IBM_moving_parameter = " << alpha_ << finl;
        }

      else if (motlu=="ponderation_shap_deriv_for_projection")
        {
          is >> pond_shap_deriv_for_proj_;
          Cerr << "ponderation_shap_deriv_for_projection = " << pond_shap_deriv_for_proj_ << finl;
        }

      else if (motlu=="PDF_regularization_for_shape_deriv")
        {
          is >> regul_PDF_shape_deriv_;
        }

      else if (motlu=="area_constraints_low_high")
        {
          is >> modif_aire_pc_low_;
          is >> modif_aire_pc_high_;
        }

      else
        {
          Cerr << "We expected state_problem, adjoint_problem, projection_problem, objective_function, objective_function_visualization, IBM_moving_parameter, ponderation_shap_deriv_for_projection, PDF_regularization_for_shape_deriv, deplacement_vector_visualization or area_constraints_low_high" << finl;
          exit();
        }

      is >> motlu;
    }

  if (nb_pb_lu != 3)
    {
      Cerr <<"We expected three problems: the state_problem, the adjoint_problem and the projection_problem" << finl;
      exit();
    }
  Cerr << "<<<<<<<<< End of reading of Probleme_Couple_Optimisation_IBM >>>>>>>>>>"<<finl;
  return is;
}

Sortie& Pb_Couple_Optimisation_IBM::printOn(Sortie& os) const
{
  return Probleme_Couple::printOn(os);
}

void Pb_Couple_Optimisation_IBM::initialize( )
{
  Probleme_Couple::initialize();

  if (!pb_etat_opt_)
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : no detected state problem."<<finl;
      Cerr<<"Case not supported"<<finl;
      Cerr<<"Aborting..."<<finl;
      abort();
    }

  int nb_source_pdf = 0;
  int nb_eq_IBM = 0;
  int mobile = 1;
  Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : search for the Source_PDF in the "<<pb_etat_opt_->nombre_d_equations()<<" equation(s) of the state problem "<<pb_etat_opt_->le_nom()<<finl;
  for (int i = 0; i <pb_etat_opt_->nombre_d_equations() ; i++)
    {
      Cerr<<"Equation : "<<pb_etat_opt_->equation(i).le_nom()<<finl;
      if(pb_etat_opt_->equation(i).que_suis_je().finit_par("IBM"))
        {
          nb_eq_IBM++;
          Equation_IBM_proto& eq_ibm = dynamic_cast<Equation_IBM_proto&>(pb_etat_opt_->equation(i));
          int id_source_pdf = eq_ibm.get_i_source_pdf();
          if (id_source_pdf == -1)
            {
              Cerr<<"Pb_Couple_Optimisation_IBM : Source_PDF not detected in IBM equation."<<finl;
              Cerr << "Either this source term is missing or not read again !!!"<< finl;
              Cerr<<"Aborting..."<<finl;
              abort();
            }
          Equation_base& eq_i = pb_etat_opt_->equation(i);
          Source_PDF_base& src = dynamic_cast<Source_PDF_base&>((eq_i.sources())[id_source_pdf].valeur());
          my_source_PDF_opt_ = src;
          PDF_model& pdf_mod_st = ref_cast_non_const(PDF_model, src.get_modele());
          if (pdf_mod_st.pdf_bilan() != 0)
            {
              Cerr<<"Model_PDF with bilan_PDF <> 0. Exit."<<finl;
              exit();
            }
          pdf_mod_st.set_PDF_mobile(mobile); // on declare l'IB mobile
          pdf_mod_st.discretiser_vitesse_shape_IBM(pb_etat_opt_); // on discretise vitesse_shape_IBM_ (PDF_mobile)
          Numero_eq_optimis_ = i;
          nb_source_pdf++;
        }
    }

  if (nb_eq_IBM==0)
    {
      Cerr<<"Pb_Couple_Optimisation_IBM : no IBM equation detected."<<finl;
      Cerr<<"The state probleme needs one IBM equation"<<finl;
      Cerr<<"Aborting..."<<finl;
      abort();
    }

  if (nb_source_pdf>1)
    {
      Cerr<<"Pb_Couple_Optimisation_IBM : more than one Source_PDF detected in state problem."<<finl;
      Cerr<<"Case not supported"<<finl;
      Cerr<<"Aborting..."<<finl;
      abort();
    }
  else if (nb_source_pdf>0)
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : one Source_PDF detected in state problem."<<finl;
    }
  else
    {
      Cerr<<"Pb_Couple_Optimisation_IBM : no Source_PDF detected."<<finl;
      Cerr<<"It needs one Source_PDF in state problem"<<finl;
      Cerr<<"Aborting..."<<finl;
      abort();
    }

  const int exist_interp = my_source_PDF_opt_->getInterpolationBool();
  if (exist_interp) my_interpolation_opt_ = my_source_PDF_opt_->getInterpolationLu();
  if (my_interpolation_opt_)
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : interpolation detected."<<finl;
      my_interpolation_opt_->discretise_PDF_mobile(pb_etat_opt_->equation(Numero_eq_optimis_).discretisation(), pb_etat_opt_->equation(Numero_eq_optimis_).domaine_dis());
    }
  my_prepro_opt_ = my_source_PDF_opt_->getpreproLu();
  if (my_prepro_opt_)
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : prepro_IBM detected for state equ."<<finl;
    }
  else
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : prepro_IBM is missing for state equ."<<finl;
      exit();
    }

  if( pb_etat_opt_->equation(Numero_eq_optimis_).que_suis_je() != pb_adjt_opt_->equation(Numero_eq_optimis_).que_suis_je() )
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : the state pb and the adjoint pb are of different type."<<finl;
      Cerr<<"      pb_etat : "<<pb_etat_opt_->equation(Numero_eq_optimis_).que_suis_je() <<finl;
      Cerr<<"      pb_adjt : "<<pb_adjt_opt_->equation(Numero_eq_optimis_).que_suis_je() <<finl;
      exit();
    }

  Equation_IBM_proto& eq_ibm_adj = dynamic_cast<Equation_IBM_proto&>(pb_adjt_opt_->equation(Numero_eq_optimis_));
  int id_source_pdf_adj = eq_ibm_adj.get_i_source_pdf();
  if (id_source_pdf_adj == -1)
    {
      Cerr<<"Pb_Couple_Optimisation_IBM : Source_PDF not detected in adjoint IBM equation."<<finl;
      Cerr << "Either this source term is missing or not read again !!!"<< finl;
      Cerr<<"Aborting..."<<finl;
      abort();
    }
  Source_PDF_base& src_adj = dynamic_cast<Source_PDF_base&>((pb_adjt_opt_->equation(Numero_eq_optimis_).sources())[id_source_pdf_adj].valeur());
  my_source_PDF_opt_adjt_ = src_adj;
  PDF_model& pdf_mod_adj = ref_cast_non_const(PDF_model, src_adj.get_modele());
  pdf_mod_adj.set_PDF_mobile(mobile);

  const int exist_interp_adjt = my_source_PDF_opt_adjt_->getInterpolationBool();
  if (exist_interp_adjt) my_interpolation_opt_adjt_ = my_source_PDF_opt_adjt_->getInterpolationLu();
  if (my_interpolation_opt_adjt_)
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : interpolation detected for adjoint problem."<<finl;
      my_interpolation_opt_adjt_->discretise_PDF_mobile(pb_adjt_opt_->equation(Numero_eq_optimis_).discretisation(), pb_adjt_opt_->equation(Numero_eq_optimis_).domaine_dis());
    }

  Equation_base& eqn_proj = pb_projection_opt_->equation(0);
  Sources& sources_proj = eqn_proj.sources();
  int size_s = sources_proj.size();
  if(size_s == 0) { Cerr<<"Pb_Couple_Optimisation_IBM::solveTimeStep : no source for projection problem."<<finl; abort();};

  Numero_src_deriv_form_ = -1;
  Nom type = "Derivee_Forme";
  if (eqn_proj.discretisation().is_vdf())
    type += "_VDF";
  else  if (eqn_proj.discretisation().is_vef())
    type += "_VEF";
  else  if (eqn_proj.discretisation().is_ef())
    type += "_EF";

  for (int nsrc=0; nsrc<size_s; nsrc++)
    {
      Source_base& ma_source = sources_proj(nsrc).valeur();
      if (ma_source.que_suis_je() == type) Numero_src_deriv_form_ = nsrc;
    }
  if (Numero_src_deriv_form_ == -1)
    {
      Cerr<<"Pb_Couple_Optimisation_IBM : no source for projection problem is of type derivee_forme."<<finl;
      abort();
    }

  const Domaine_dis_base& le_dom_dis = pb_projection_opt_->domaine_dis();

  // source_derivee_forme : scalaire par element
  pb_projection_opt_->discretisation().discretiser_champ("champ_elem",le_dom_dis,"source_derivee_forme","",1,0., source_derivee_forme_);
  DoubleTab& sourceArray = source_derivee_forme_->valeurs();
  sourceArray = 0.;

  // normal_derivee_forme : normal par element
  pb_projection_opt_->discretisation().discretiser_champ("champ_elem",le_dom_dis,"normal_derivee_forme","",3,0., normal_derivee_forme_);
  DoubleTab& normSourceArray = normal_derivee_forme_->valeurs();
  normSourceArray = 0.;

  // Fonction cout scalaire par element

  if (fonction_cout_lu_)
    {
      Champ_Don_base& ch_fonction_cout_lu = ref_cast(Champ_Don_base,fonction_cout_lu_.valeur());
      const int nb_comp = ch_fonction_cout_lu.nb_comp();
      pb_projection_opt_->discretisation().discretiser_champ("champ_elem",le_dom_dis,"fonction_cout","",nb_comp,0., fonction_cout_);
      for (int n = 0; n < nb_comp; n++) fonction_cout_->fixer_nom_compo(n, ch_fonction_cout_lu.le_nom() + (nb_comp > 1 ? Nom(n) :""));
      // PL: Il faut faire nommer_completer_champ_physique les 2 champs (plantage sinon pour une fonction de type Champ_fonc_tabule)
      eqn_proj.discretisation().nommer_completer_champ_physique(eqn_proj.domaine_dis(),ch_fonction_cout_lu.le_nom(),"",fonction_cout_lu_,eqn_proj.probleme());
      eqn_proj.discretisation().nommer_completer_champ_physique(eqn_proj.domaine_dis(),ch_fonction_cout_lu.le_nom(),"",fonction_cout_,eqn_proj.probleme());
      fonction_cout_->valeurs() = 0.;
      fonction_cout_->affecter(ch_fonction_cout_lu);
    }
  else
    {
      pb_projection_opt_->discretisation().discretiser_champ("champ_elem",le_dom_dis,"fonction_cout","",1,0., fonction_cout_);
      DoubleTab& coutArray = fonction_cout_->valeurs();
      coutArray = 0.;
    }


  Cerr<<"Shape optimization coupled problem :";
  for (int ipb=0; ipb < nb_problemes(); ipb++) Cerr<<" "<<probleme(ipb).le_nom();
  Cerr<<finl;

}

bool Pb_Couple_Optimisation_IBM::initTimeStep(double dt)
{
  const Domaine_dis_base& le_dom_dis = pb_etat_opt_->domaine_dis();
  int nb_elem = le_dom_dis.nb_elem();
  const Schema_Temps_base& sch=pb_etat_opt_->schema_temps();
  double temps=sch.temps_courant();
  double pdtps = sch.pas_de_temps();
  bool level_set = my_source_PDF_opt_->get_modele().get_use_pseudo_level_set_moving_PDF();

  const Domaine& dom = le_dom_dis.domaine();
  int nb_faces_elem = dom.nb_faces_elem();
  const Domaine_VF& the_dom_VF = ref_cast(Domaine_VF,le_dom_dis);
  const IntTab& elem_face = the_dom_VF.elem_faces();
  const IntTab& face_voisins = the_dom_VF.face_voisins();
  // VDF/VEF by default
  int nb_dof_elem = nb_faces_elem;
  // Formulation EF
  if (pb_projection_opt_->equation(0).discretisation().is_ef())
    {
      nb_dof_elem = dom.nb_som_elem();
    }
  else if (!(pb_projection_opt_->equation(0).discretisation().is_vdf()) && !(pb_projection_opt_->equation(0).discretisation().is_vef()))
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM::calcul_derivee_forme_IBM: discretisation is not EF, VEF or VDF. Aborting... "<<pb_projection_opt_->equation(0).discretisation()<<finl;
    }
  const IntTab& elems = (pb_projection_opt_->equation(0).discretisation().is_ef() ? dom.les_elems() : the_dom_VF.elem_faces());

  // Mise a jour de la fonction cout (quantité élémentaire)
  if (fonction_cout_lu_)
    {
      fonction_cout_lu_->mettre_a_jour(temps);
      fonction_cout_->affecter(fonction_cout_lu_.valeur());
    }
  if (visu_cout_) Save_Med_File_fonction_cout(fonction_cout_->valeurs());

  // Mise a jour de la frontiere IBM (quantité élémentaire): barycentre, normale et aire
  DoubleTab& normSourceArray = normal_derivee_forme_->valeurs(); //vecteur normalise par element donnant le sens
  const DoubleTab& deplacement_scalaire = pb_projection_opt_->equation(0).inconnue().valeurs();
  int dim0 = (level_set ? deplacement_scalaire.dimension(0) : normSourceArray.dimension(0));
  DoubleTab vecteur_deplacement(dim0, normSourceArray.dimension(1));
  if (level_set)  // <vecteur deplacement pour approche pseudo level set (par dof)
    {
      int relev_vois_e = 0; // copie ou non normSourceArray pour elements voisins par une face
      // on boucle sur les cellules ayant normSourceArray non nul
      vecteur_deplacement = 0.;
      DoubleTrav contrib_i(dim0);
      for (int e=0; e<normSourceArray.dimension(0); e++)
        {
          double norm = 0.;
          for (int k=0; k<normSourceArray.dimension(1); k++) norm += normSourceArray(e,k)*normSourceArray(e,k);
          if (norm > 1.e-10) // Existe normSourceArray non nul pour l element e
            {
              for (int il=0; il<nb_dof_elem; il++)
                {
                  int i = elems(e,il);
                  for (int k=0; k<vecteur_deplacement.dimension(1); k++) vecteur_deplacement(i,k) += normSourceArray(e,k);
                  contrib_i(i) += 1.;
                }
              if (relev_vois_e)
                {
                  // copie normSourceArray pour elements voisins par une face
                  for (int fac=0; fac<nb_faces_elem; fac++)
                    {
                      int num_fac = elem_face(e, fac);
                      for (int voisin=0; voisin<2; voisin++) // Deux voisins seulement
                        {
                          int elem_voi = face_voisins(num_fac,voisin);
                          if ( (elem_voi!=-1) && (elem_voi!=e) )
                            {
                              norm = 0.;
                              for (int k=0; k<normSourceArray.dimension(1); k++) norm += normSourceArray(elem_voi,k)*normSourceArray(elem_voi,k);
                              if (norm < 1.e-10) // N'existe pas normSourceArray chez le voisin elem_voi
                                {
                                  for (int il_v=0; il_v<nb_dof_elem; il_v++)
                                    {
                                      int i_v = elems(elem_voi,il_v);
                                      for (int k=0; k<vecteur_deplacement.dimension(1); k++) vecteur_deplacement(i_v,k) += normSourceArray(e,k);
                                      contrib_i(i_v) += 1.;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
      for (int i=0; i<dim0; i++)
        {
          if (contrib_i(i) >0.)
            {
              double nor_nor = 0.;
              for (int k=0; k<vecteur_deplacement.dimension(1); k++) nor_nor += (vecteur_deplacement(i,k)/contrib_i(i))*(vecteur_deplacement(i,k)/contrib_i(i));
              if (nor_nor > 1.e-10)
                for (int k=0; k<vecteur_deplacement.dimension(1); k++) vecteur_deplacement(i,k) /= sqrt(nor_nor);
            }
        }

    }
  else // <vecteur deplacement pour approche particulaire sous-contrainte (par elem)
    {
      for (int e=0; e<dim0; e++)
        {
          for (int k=0; k<vecteur_deplacement.dimension(1); k++) vecteur_deplacement(e,k) = normSourceArray(e,k);
        }
    }

  bool iShapDeplOK = false;
  double bilan = 0.;

  if (temps >= pdtps && alpha_ > 0.)
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : updating the shape following shape displacements"<<finl;
      const DoubleTab& aire = my_source_PDF_opt_->get_champ_aire().valeurs();
      bilan = 0.;
      for (int e=0; e<nb_elem; e++)
        if (aire(e)>0.) bilan += aire(e);
      Cerr<<"(IBM) Area balance before shape moving = "<<bilan<<finl;

      // Appel au prepro IBM pour le champ h_max_elem (limitation du deplacement scalaire a une maille voisine si alpha_ = 1)
      const DoubleTab& h_max_elem = (level_set ? my_prepro_opt_->get_h_max_node() : my_prepro_opt_->get_h_max_elem());

      // Definition du vecteur deplacement elementaire suivant le sens de normal_derivee_forme_ (sens du deplacement)
      double deplacement_scalaire_e;
      for (int e=0; e<vecteur_deplacement.dimension(0); e++)
        {
          if (level_set)  // deplacement pour approche pseudo level set
            {
              //on est aux dof de la projection (e = dof)
              deplacement_scalaire_e = deplacement_scalaire(e);
            }
          else // deplacement pour approche particulaire sous-contrainte
            {
              //on est aux elements (e = elem)
              deplacement_scalaire_e =0.;
              for (int il=0; il<nb_dof_elem; il++)
                {
                  int i = elems(e,il);
                  deplacement_scalaire_e += deplacement_scalaire(i) / nb_dof_elem;
                }
            }
          double abs_depl = abs(deplacement_scalaire_e);
          double sign_depl = 0.;
          if (abs_depl > 1.e-10) sign_depl = deplacement_scalaire_e / abs_depl;
          double min_en_depl = min(abs_depl * alpha_, h_max_elem(e)) * sign_depl;
          for (int k=0; k<vecteur_deplacement.dimension(1); k++) vecteur_deplacement(e,k) *= min_en_depl;
        }
      iShapDeplOK = true;
    }
  else vecteur_deplacement *= 0.;

  // Vitesse deplacement IB pour postraitement
  DoubleTrav vitesse_deplacement(normSourceArray);
  if (level_set)  // projection aux elements
    {
      const DoubleTab& aire = my_source_PDF_opt_->get_champ_aire().valeurs();
      for (int e=0; e<normSourceArray.dimension(0); e++)
        {
          for (int il=0; il<nb_dof_elem; il++)
            {
              int i = elems(e,il);
              for (int k=0; k<vecteur_deplacement.dimension(1); k++) vitesse_deplacement(e,k) += vecteur_deplacement(i,k)*(aire(e)>0.?1.:0.) / nb_dof_elem;
            }
        }
    }
  else
    {
      for (int e=0; e<vecteur_deplacement.dimension(0); e++)
        {
          for (int k=0; k<vecteur_deplacement.dimension(1); k++) vitesse_deplacement(e,k) = vecteur_deplacement(e,k) ;
        }
    }
  vitesse_deplacement /= pdtps;
  PDF_model& pdf_mod_etat = ref_cast_non_const(PDF_model, my_source_PDF_opt_->get_modele());
  pdf_mod_etat.set_vitesse_shape_IBM(vitesse_deplacement);

  // Deplacement vectoriel elementaire des barycentres des faces IBM et maj prepro_IBM et interpolation de l'eq. etat
  if (iShapDeplOK)
    {
      double raid = my_source_PDF_opt_->get_modele().raid();
      if (level_set)
        my_source_PDF_opt_->update_pseudo_level_set_IBM(vecteur_deplacement, 1.0);
      else
        my_source_PDF_opt_->update_elem_IBM(vecteur_deplacement, 1.0, raid);

      const DoubleTab& aire = my_source_PDF_opt_->get_champ_aire().valeurs();
      bilan = 0.;
      for (int e=0; e<nb_elem; e++) bilan += aire(e);
      if (!area_ref_set_)
        {
          area_ref_set_ = 1;
          area_ref_ = bilan;
          Cerr<<"(IBM) Reference area balance after first shape moving = "<<bilan<<finl;
        }
      else
        {
          Cerr<<"(IBM) Area balance after shape moving = "<<bilan<<finl;
          if ( (bilan <= area_ref_ * (modif_aire_pc_low_)) || (bilan >= area_ref_ * (modif_aire_pc_high_)) )
            {
              alpha_ =-1.0; // Bloquage de la modification de la geometrie
              Cerr<<"(IBM) Area constraint: shape moving is stopped. Percent =  "<<(bilan/area_ref_ )<<finl;
            }
        }


      // Appel au prepro IBM pour generer la nouvelle description discrete de la frontiere
      // Eq. adjointe
      my_source_PDF_opt_adjt_->get_fields_from_prepro(my_prepro_opt_);
      if(my_source_PDF_opt_adjt_->getInterpolationBool() == true)
        {
          Interpolation_IBM_base& my_interp_adjt = ref_cast_non_const(Interpolation_IBM_base, my_source_PDF_opt_adjt_->getInterpolationLu());
          my_interp_adjt.set_fields_from_prepro_to_interp(my_prepro_opt_);
        }
      my_source_PDF_opt_adjt_->set_variable_imposee();
    }

  bool ok =  Probleme_Couple::initTimeStep(dt);
  return ok;
}

bool Pb_Couple_Optimisation_IBM::solveTimeStep()
{

  bool ok;
  suppProblem(pb_projection_opt_);
  ok = Probleme_Couple::solveTimeStep();

  addProblem(pb_projection_opt_);
  suppProblem(pb_etat_opt_);
  suppProblem(pb_adjt_opt_);

  const Schema_Temps_base& sch=pb_etat_opt_->schema_temps();
  double temps=sch.temps_courant();
  double pdtps = sch.pas_de_temps();
  if (temps > pdtps)
    {
      // Calcul terme derivee de forme
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM : update of shape derivative term in projection eq." <<finl;
      DoubleVect my_bilan;
      my_source_PDF_opt_->compute_source_term_PDF(my_bilan);
      my_source_PDF_opt_adjt_->compute_source_term_PDF(my_bilan);
      calcul_derivee_forme_IBM();

      // derivee de forme dans le terme source de equation de projection
      Source_base& ma_source = (pb_projection_opt_->equation(0).sources())(Numero_src_deriv_form_).valeur();
      Terme_Derivee_Forme_base& ma_source_DF_eq = ref_cast(Terme_Derivee_Forme_base, ma_source);
      DoubleTab& sourceEqArray = ref_cast_non_const(DoubleTab, ma_source_DF_eq.get_source_derivee_forme());
      DoubleTab& sourceDFArray = source_derivee_forme_->valeurs();
      assert(sourceEqArray.dimension(0) == sourceDFArray.dimension(0));
      assert(sourceEqArray.dimension(1) == sourceDFArray.dimension(1));
      for (int i=0; i<sourceEqArray.dimension(0); i++)
        {
          for (int j=0; j<sourceEqArray.dimension(1); j++) sourceEqArray(i,j) = 0.-sourceDFArray(i,j)*pond_shap_deriv_for_proj_;
        }
      sourceEqArray.echange_espace_virtuel();
      ma_source_DF_eq.set_source_derivee_forme(sourceEqArray);
      ma_source_DF_eq.mettre_a_jour(pb_projection_opt_->equation(0).probleme().schema_temps().temps_courant());
    }

  ok = Probleme_Couple::solveTimeStep();

  addProblem(pb_etat_opt_);
  addProblem(pb_adjt_opt_);

  return ok;
}

void Pb_Couple_Optimisation_IBM::calcul_derivee_forme_IBM()
{

  // Calcul terme derivee de forme à partir du source_pdf de pb_etat + pb_adjoint
  DoubleTab& sourceArray = source_derivee_forme_->valeurs();
  DoubleTab& normSourceArray = normal_derivee_forme_->valeurs();
  DoubleTab& fonctCout = fonction_cout_->valeurs();
  DoubleTab coutArray(fonctCout);
  assert(fonctCout.dimension(1)==1);
  assert(fonctCout.dimension(1)==sourceArray.dimension(1));
  const Champ_Don_base& champ_aire = my_source_PDF_opt_->get_champ_aire();
  const DoubleTab& aire = champ_aire.valeurs();
  const DoubleTab& edp_val_i = my_source_PDF_opt_->get_source_pdf();
  const DoubleTab& adjoint = pb_adjt_opt_->equation(Numero_eq_optimis_).inconnue().valeurs();
  const int interp_adjoint = my_source_PDF_opt_adjt_->getInterpolationBool();

  // dimentional verifications (0)
  int nb_dof = edp_val_i.dimension(0);
  // dimentional verifications (1)
  int nb_comp = edp_val_i.dimension(1);
  int dim_esp = Objet_U::dimension;
  assert(nb_comp<=dim_esp);

  const Domaine_dis_base& le_dom_dis = pb_etat_opt_->domaine_dis();
  const Domaine& dom = le_dom_dis.domaine();
  const Domaine_VF& the_dom_VF = ref_cast(Domaine_VF,le_dom_dis);
  int nb_elem = le_dom_dis.nb_elem();
  int nb_elem_tot = le_dom_dis.nb_elem_tot();
  int nb_som_elem = dom.nb_som_elem();
  int nb_faces_elem = dom.nb_faces_elem();

  // VDF/VEF by default
  const IntTab& faces_som = the_dom_VF.face_sommets();
  int nb_som_face = the_dom_VF.nb_som_face();
  int nb_dof_elem = nb_faces_elem;
  int is_face_dis = 1;
  const IntTab& elems = (pb_adjt_opt_->equation(Numero_eq_optimis_).discretisation().is_ef() ? le_dom_dis.domaine().les_elems() : the_dom_VF.elem_faces());
  if (pb_adjt_opt_->equation(Numero_eq_optimis_).discretisation().is_ef())
    {
      // Formulation EF
      nb_dof_elem = nb_som_elem;
      is_face_dis = 0;
    }
  else if (!(pb_adjt_opt_->equation(Numero_eq_optimis_).discretisation().is_vdf()) && !(pb_adjt_opt_->equation(Numero_eq_optimis_).discretisation().is_vef()))
    {
      Cerr<<"(IBM) Pb_Couple_Optimisation_IBM::calcul_derivee_forme_IBM: discretisation is not EF, VEF or VDF. Aborting... "<<pb_adjt_opt_->equation(Numero_eq_optimis_).discretisation()<<finl;
    }

  DoubleTrav nor(nb_dof, dim_esp);
  DoubleTrav dist(nb_dof);
  my_prepro_opt_->calculer_normal_proj_solid(nor, dist);
  IntLists elem_voisins(nb_elem_tot);
  DoubleTab& aire_ncst = ref_cast_non_const(DoubleTab, aire);
  bool all_elem_vois = true;
  my_source_PDF_opt_->compute_NeighNode_IBM_elem(aire_ncst, elem_voisins, all_elem_vois);

  // Loop on elements for weight computation on dof
  DoubleTrav weight(nb_dof);
  weight = 0.;
  for (int e=0; e<nb_elem; e++)
    {
      if (aire(e)>0.)
        {
          for (int il=0; il<nb_dof_elem; il++)
            {
              int i = elems(e,il);
              weight(i) += 1.;
            }
        }
    }

  // Loop on elements for shape derivative computation
  sourceArray = 0.;
  normSourceArray = 0.;
  coutArray = 0.;
  double regualpha_i=1.;
  double eps = 1e-12;
  double d1;
  double PDF_contr_min = 1.e+30;
  double PDF_contr_max = -1.e+30;
  double cout_contr_min = 1.e+30;
  double cout_contr_max = -1.e+30;
  double lim_p = 0.5;
  double lim_n = 0. - lim_p;
  DoubleTrav nor_e(1, dim_esp);

  for (int e=0; e<nb_elem; e++)
    {
      if (aire(e)>0.)
        {
          nor_e = 0.;
          for (int il=0; il<nb_dof_elem; il++)
            {
              int i = elems(e,il);

              // scalar PDF_source contribution at dof
              double ps_edp_i = 0.;
              if (regul_PDF_shape_deriv_) regualpha_i= my_source_PDF_opt_->fonct_regul_PDF(e, dist(i));
              // Cerr<<"e, i, dist(i), regualpha_i = "<<e<<" "<<i<<" "<<dist(i)<<" "<<regualpha_i<<finl;
              for (int nc=0; nc<nb_comp; nc++)
                if(interp_adjoint) ps_edp_i += regualpha_i * edp_val_i(i,nc) * adjoint(i,nc) / weight(i);
                else ps_edp_i += regualpha_i * edp_val_i(i,nc) / weight(i);

              // sum of nodal normals weighted by nodal PDF_source
              if (!(is_face_dis))
                {
                  for (int k=0; k<dim_esp; k++) nor_e(0, k) += nor(i, k)*ps_edp_i;
                }
              else if (is_face_dis)
                {
                  for (int no=0; no<nb_som_face; no++)
                    {
                      int inode = faces_som(i, no);
                      for (int k=0; k<dim_esp; k++) nor_e(0, k) += nor(inode, k)*ps_edp_i;
                    }
                }
            }

          // arrith. mean of PDF weighted normal resulting in normalized vector folowing (-grad_gamma J) by element, from lower to higher PDF_source
          d1 = 0.0;
          for (int k=0; k<dim_esp; k++) d1 += nor_e(0, k)*nor_e(0, k);
          d1 = sqrt(d1);
          if (d1  > eps )
            for (int k=0; k<dim_esp; k++) normSourceArray(e, k) = nor_e(0, k) / d1 ;

          sourceArray(e,0) = std::min(d1, lim_p);
          if (sourceArray(e,0) < PDF_contr_min) PDF_contr_min = sourceArray(e,0);
          if (sourceArray(e,0) > PDF_contr_max) PDF_contr_max = sourceArray(e,0);

          // Elementary cost function contribution

          // boucle sur les voisin vois de e :
          int nb_elem_voi = elem_voisins[e].size();
          int ps_pos = 0; // nb element voisin ayant une contribution positive a prendre en compte
          int ps_neg = 0; // nb element voisin ayant une contribution  negative a prendre en compte
          double meanCout_pos = 0.;
          double meanCout_neg = 0.;
          if (nb_elem_voi != 0)
            {
              for (int voi=0; voi<nb_elem_voi; voi++)
                {
                  // Le voisin doit avoir aire() = 0
                  // Cerr<<"voisin = "<<(elem_voisins[e])[voi]<<" aire_voisin = "<<aire((elem_voisins[e])[voi])<<finl;
                  if (aire((elem_voisins[e])[voi]) <= 0.)
                    {
                      int elem_voi = (elem_voisins[e])[voi];
                      // pour les noeuds i de chaque voisin: ps nor(i, k) * normSourceArray(e,k) donne le signe
                      // de la contribution qui doit etre la meme pour tous les noeuds i de chaque voisin
                      double ps_voi_ref = 0.;
                      int iok = 0;
                      for (int il=0; il<nb_dof_elem; il++)
                        {
                          int i = elems(elem_voi,il);
                          double ps_voi_i =0.;
                          for (int k=0; k<dim_esp; k++) ps_voi_i += nor(i, k) * normSourceArray(e,k) ;
                          if (iok == 0 && ps_voi_i != 0.)
                            {
                              ps_voi_ref = ps_voi_i; // premier ps non nul
                              iok = 1;
                            }
                          if (ps_voi_ref*ps_voi_i < 0)
                            {
                              iok = 0;
                              break;
                            }
                          // Cerr<<"Element voisin = "<<elem_voi<<" , ps_voi_i = "<<ps_voi_i<<" , iok, ps_voi_ref = "<<iok<<" "<<ps_voi_ref<<finl;
                        }
                      // recuperer le signe du ps et faire une moyenne sur l'element voisin elem_voi
                      if (iok == 1 && ps_voi_ref > 0.)
                        {
                          ps_pos += 1;
                          meanCout_pos += fonctCout(elem_voi,0);
                        }
                      else if (iok == 1 && ps_voi_ref<0.)
                        {
                          ps_neg += 1;
                          meanCout_neg += fonctCout(elem_voi,0);
                        }
                    }
                }
              // Moyenne fct cout pour les contributions + et -
              if (ps_pos != 0 && ps_neg != 0)
                coutArray(e,0) = std::max(std::min(meanCout_pos/ps_pos - meanCout_neg/ps_neg, lim_p), lim_n);
              else coutArray(e,0) = 0.;
            }
          if (coutArray(e,0) < cout_contr_min) cout_contr_min = coutArray(e,0);
          if (coutArray(e,0) > cout_contr_max) cout_contr_max = coutArray(e,0);

          // Cerr<<"element "<<e<<" : normSourceArray = ";
          // for (int k=0; k<dim_esp; k++) Cerr<<normSourceArray(e, k) <<" ";
          // Cerr<<" ; sourceArray = "<<sourceArray(e,0);
          // Cerr<<" ; coutArray = "<<coutArray(e,0)<<" nb voisins ="<<nb_elem_voi<<" ; pos/neg = "<<ps_pos<<" "<<ps_neg<<finl;
        }
    }

  //Normalisation eventuelle contribution PDF et cout entre 0 et +1
  double bilan = 0.;
  for (int e=0; e<nb_elem; e++)
    {
      if (aire(e)>0.)
        {
          // double source = (sourceArray(e,0) - PDF_contr_min) /(PDF_contr_max - PDF_contr_min + 1.e-6) ;
          // double source = (log(sourceArray(e,0)) - log(PDF_contr_min)) /(log(PDF_contr_max) - log(PDF_contr_min) + 1.e-6);
          // double coutsd = (coutArray(e,0) - cout_contr_min) /(cout_contr_max - cout_contr_min + 1.e-6) ;
          double source = sourceArray(e,0) ;
          double coutsd = coutArray(e,0) ;
          sourceArray(e,0) = source - coutsd;
          bilan += sourceArray(e,0);
        }
    }

  sourceArray.echange_espace_virtuel();
  normSourceArray.echange_espace_virtuel();
  Cerr<<"(IBM) Pb_Couple_Optimisation_IBM::calcul_derivee_forme_IBM: Balance (// Balance) = ";
  Cerr<<bilan<<" "<<mp_sum(bilan)<<finl;
}

void Pb_Couple_Optimisation_IBM::Save_Med_File_fonction_cout(DoubleTab& fcout)
{
  const Domaine& le_dom = pb_projection_opt_.valeur().domaine();
  const Schema_Temps_base& sch=pb_etat_opt_->schema_temps();
  double temps=sch.temps_courant();
  double temps_ecoule=sch.temps_calcul();
  double dt_impr=sch.temps_impr();
  int step=int(temps_ecoule/dt_impr);

  Nom nom_fichier_med_Out_ = "objective_function.med";
  if (sch.nb_pas_dt() < 1)
    {
      ecr_med_fonction_cout_.set_file_name_and_dom(nom_fichier_med_Out_, le_dom);
      ecr_med_fonction_cout_.ecrire_domaine(false);
    }
  else if (step > nb_save_call_)
    {
      nb_save_call_ += 1;
      const Nom& type_elem = le_dom.type_elem()->que_suis_je();
      int nb_comp = fcout.dimension(1);
      Noms nom_cout(nb_comp);
      Noms unites_cout(nb_comp);
      for (int k=0; k<nb_comp; k++) unites_cout[k] = "-";
      for (int k=0; k<nb_comp; k++) nom_cout[k] = "-";
      ecr_med_fonction_cout_.ecrire_champ("CHAMPMAILLE", "objective_function", fcout, unites_cout, nom_cout, type_elem, temps);
    }
}
