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

#include <Champ_Generique_Transformation.h>
#include <Probleme_base.h>
#include <Domaine_VF.h>
#include <TRUSTTabs.h>
#include <ParserView.h>

#include <Discretisation_base.h>
#include <Interprete.h>
#include <Equation_base.h>
#include <Entree_complete.h>
#include <Milieu_base.h>
#include <Postraitement.h>
#include <Synonyme_info.h>
#include <Param.h>

#include <Linear_algebra_tools_impl.h>
#include <MD_Vector.h>
#include <MD_Vector_composite.h>
#include <Domaine_dis_base.h>

Implemente_instanciable(Champ_Generique_Transformation,"Transformation",Champ_Gen_de_Champs_Gen);
Add_synonym(Champ_Generique_Transformation,"Champ_Post_Transformation");

// XD transformation champ_post_de_champs_post champ_post_transformation BRACE To create a field with a transformation
// XD_CONT using source fields and x, y, z, t. If you use in your datafile source refChamp { Pb_champ pb pression }, the
// XD_CONT field pression may be used in the expression with the name pression_natif_dom; this latter is the same as
// XD_CONT pression. If you specify nom_source in refChamp bloc, you should use the alias given to pressure field. This
// XD_CONT is avail for all equations unknowns in transformation.

Sortie& Champ_Generique_Transformation::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Champ_Generique_Transformation::readOn(Entree& s )
{
  Champ_Gen_de_Champs_Gen::readOn(s);
  verifier_coherence_donnees();
  return s ;
}

//-methode : type of transformation to perform
// possible methods :
//        formule : function of generic source fields, x, y, z and t
//        vecteur : construction of a vector from the provided data
//        produit_scalaire : dot product of two specified vectors
//        norme : norm of a specified vector
//        composante : creation of a scalar field from a component of a vector field
//-expression : expression of the transformation (for the formule and vecteur methods)
//-numero : number of the component to extract (only for the composante method)
//-localisation : localization of the storage field
//-unite : to specify the unit of a field to improve the readability of post-processings
void Champ_Generique_Transformation::set_param(Param& param) const
{
  Champ_Gen_de_Champs_Gen::set_param(param);
  param.ajouter("methode",&methode_,Param::REQUIRED); // XD_ADD_P chaine(into=["produit_scalaire","norme","vecteur","formule","composante"])
  // XD_CONT methode 0 methode norme : will calculate the norm of a vector given by a source field NL2 methode
  // XD_CONT produit_scalaire : will calculate the dot product of two vectors given by two sources fields NL2 methode
  // XD_CONT composante numero integer : will create a field by extracting the integer component of a field given by a
  // XD_CONT source field NL2 methode formule expression 1 : will create a scalar field located to elements using
  // XD_CONT expressions with x,y,z,t parameters and field names given by a source field or several sources fields. NL2
  // XD_CONT methode vecteur expression N f1(x,y,z,t) fN(x,y,z,t) : will create a vector field located to elements by
  // XD_CONT defining its N components with N expressions with x,y,z,t parameters and field names given by a source
  // XD_CONT field or several sources fields.
  param.ajouter("unite", &unite_);  // XD_ADD_P chaine
  // XD_CONT will specify the field unit
  param.ajouter_non_std("expression",(this)); // XD_ADD_P listchaine
  // XD_CONT expression 1 see methodes formule and vecteur
  param.ajouter_non_std("numero",(this)); // XD_ADD_P entier
  // XD_CONT numero 1 see methode composante
  param.ajouter("localisation",&localisation_); // XD_ADD_P chaine
  // XD_CONT localisation 1 type_loc indicate where is done the interpolation (elem for element or som for node). The
  // XD_CONT optional keyword methode is limited to calculer_champ_post for the moment
}

int Champ_Generique_Transformation::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="expression")
    {
      int dim;
      is >> dim;
      les_fct.dimensionner(dim);
      if ((Motcle(methode_) == "formule") && (dim!=1) )
        {
          Cerr << "The integer N of the expression with FORMULE method must be equal to 1." <<finl;
          Cerr << "Use the VECTEUR expression method for N="<<dim<<finl;
          exit();
        }
      for (int i=0; i<dim; i++)
        is>>les_fct[i];
      return 1;
    }
  else if (mot=="numero")
    {
      int num_compo;
      is >> num_compo;
      nb_comp_ = 1;
      nature_ch = scalaire;
      /*     test component number: we cannot rely on the space dimension in the case of a multi_scalaire field
             Problem: it is nb_comp_ that will subsequently carry the component dimension and it is initialized to 1 here to avoid a bug
             The test is moved to Champ_Generique_Transformation::get_champ

             if ((num_compo<0) || (num_compo>dimension-1))
             {
             //   Cerr<<"Component number must be >=0 and < "<<dimension<<" for composante option."<<finl;
             Cerr<<"Component number must be >=0 for composante option."<<finl;
             exit();
             } */
      les_fct.dimensionner(num_compo+1);
      return 1;
    }
  else
    return Champ_Gen_de_Champs_Gen::lire_motcle_non_standard(mot,is);
}

void Champ_Generique_Transformation::verifier_coherence_donnees()
{
  if ((Motcle(methode_) == "produit_scalaire") || (Motcle(methode_) == "norme"))
    les_fct.dimensionner(1);
  else if ((Motcle(methode_) == "formule") || (Motcle(methode_) == "vecteur"))
    {
      if (les_fct.size()<1)
        {
          Cerr<<"The expression must be read for the methods fonction and vecteur."<<finl;
          exit();
        }
    }
  else
    {
      Motcles methods(6);
      methods[0] = "formule";
      methods[1] = "vecteur";
      methods[2] = "produit_scalaire";
      methods[3] = "norme";
      methods[4] = "composante";
      methods[5] = "composante_normale";
      if (methods.search(methode_) < 0)
        {
          Cerr<<"Error in  "<<que_suis_je()<<finl;
          Cerr<<"Read keyword : " << methode_ << "\n Allowed items are : "<<methods<<finl;
          exit();
        }
    }

  if ((Motcle(methode_) == "vecteur"))
    {
      // We try to set the number of components immediately with vecteur
      nb_comp_ = les_fct.size();

      // the verification of sources is done during completer for reference sources
    }
}
void Champ_Generique_Transformation::verifier_localisation()
{
  //TODO DG temporary ?
  const Domaine_dis_base& domaine_dis = get_source(0).get_ref_domaine_dis_base();
  if (domaine_dis.que_suis_je()=="Domaine_DG")
    {
      if (!(localisation_=="elem"))
        {
          Cerr << "Error in Champ_Generique_Transformation::verifier_localisation" << finl;
          Cerr << "with DG the only possible localisation for Transformation is \"elem\"" << finl;
          exit();
        }
    }

  Motcles localisations(4);
  localisations[0] = "elem";
  localisations[1] = "som";
  localisations[2] = "faces";
  localisations[3] = "elem_som";

  if (localisations.search(localisation_) < 0)
    {
      Cerr << "Error in "<<que_suis_je()<<finl;
      Cerr << "Read keyword : " << localisation_ << "\n Allowed items are " << localisations << finl;
      exit();
    }
}
void Champ_Generique_Transformation::completer(const Postraitement_base& post)
{

  bool sources_location_ok = true;
  int nb_source_fictive=0;
  int nb_sources = get_nb_sources();
  if (nb_sources==0)
    {
      nb_source_fictive = 1;
      OWN_PTR(Champ_Generique_base)& new_src = sources_.add(OWN_PTR(Champ_Generique_base)());
      Nom nom_probleme = ref_cast(Postraitement,post).probleme().le_nom();
      const Objet_U& ob = interprete().objet(nom_probleme);
      const Probleme_base& pb = ref_cast(Probleme_base,ob);
      const Nom& nom_champ = pb.equation(0).inconnue().le_nom();
      Entree s_bidon;
      Nom ajout("");
      ajout += " Champ_Post_refChamp { Pb_champ ";
      ajout += pb.le_nom();
      ajout += " ";
      ajout += nom_champ;
      ajout += " }";
      Entree_complete s_complete(ajout,s_bidon);
      s_complete>>new_src;
      nb_sources = get_nb_sources();
    }

  Champ_Gen_de_Champs_Gen::completer(post);

  if (nb_source_fictive==1)
    {
      if (get_nb_sources()>1)
        {
          assert( sources_.size( ) == 1 );
          // a source was added unnecessarily ...
          sources_.vide();
          nb_source_fictive = 0;
        }
    }

  nb_sources = get_nb_sources();

  //Determine the number of components to be assigned
  //to the storage space returned by this post-processing field by
  //searching for the source that has the largest number of components

  nb_comp_ = 1;
  nature_ch = scalaire;

  if (nb_source_fictive==0)
    {
      for (int i=0; i<nb_sources; i++)
        {
          OWN_PTR(Champ_base) source_espace_stockage_tmp;
          const Champ_base& source = get_source(i).get_champ_without_evaluation(source_espace_stockage_tmp);

          if ((Motcle(methode_) == "vecteur"))
            {

              if (source.is_vectorial())
                {
                  Cerr<<que_suis_je()<<" The source fields must be of scalar nature for option vecteur."<<finl;
                  exit();
                }
            }
          int nb_comp = source.nb_vect_comp();
          nb_comp_ = (nb_comp_<nb_comp)?  nb_comp:nb_comp_;

          if (source.is_vectorial())
            nature_ch = vectoriel;
          else if (source.nature_du_champ()==multi_scalaire)
            nature_ch = multi_scalaire;
        }
    }
  //If no source field was specified, we add one
  //to provide access to the discretized domain ...

  Noms sources_location;

  if( nb_source_fictive > 0 )
    {
      //special treatment when we have fictives sources
      //the user needs to provide a location
      if( localisation_ == "??" )
        {
          // Take the localization of the unknown


          Cerr << "Error in Champ_Generique_Transformation::completer "<<finl;
          Cerr << "No sources were specified. So you need to specify at least the location."<<finl;
          Cerr << "Aborting..."<<finl;
          Cerr << "ToDo for developer: guess location by taking the pb.equation(0).inconnue() localization." << finl;
          Process::abort( );
        }
      else
        {
          fictive_source_ = true;
        }
    }
  else
    {

      //real sources
      for( int i=0; i<nb_sources; i++)
        {
          //loop over sources to check that all sources have the same location
          bool use_directive = false;
          OWN_PTR(Champ_base) my_field;

          const Champ_base& source_i = get_source( i ).get_champ_without_evaluation( my_field );
          if( source_i.a_un_domaine_dis_base( ) )
            {
              const Domaine_dis_base& domaine_source_i = source_i.domaine_dis_base( );
              if( ! sub_type( Domaine_VF , domaine_source_i ) )
                {
                  Cerr << "Warning in Champ_Generique_Transformation::completer "<<finl;
                  Cerr << "The source number "<<i<<" has a domaine_dis_base but not of type Domaine_VF so use directive to determine the location"<<finl;
                  use_directive=true;
                }
              if( ! use_directive )
                {
                  const DoubleTab& values_source_i = source_i.valeurs( );
                  const Domaine_VF& zvf_source_i = ref_cast( Domaine_VF, domaine_source_i );


                  MD_Vector md;
                  md = values_source_i.get_md_vector( );

                  //composite case, in particular Champ_{P0,Face}_PolyMAC_HFV...
                  if (get_source(i).get_discretisation().is_poly_family() && sub_type( MD_Vector_composite, md.valeur( )))
                    {
                      const MD_Vector& md0 = ref_cast(MD_Vector_composite, md.valeur()).get_desc_part(0);
                      if (md0 == zvf_source_i.domaine( ).les_elems().get_md_vector( ))
                        sources_location.add( "elem" );
                      else if (md0 == zvf_source_i.face_sommets( ).get_md_vector( ))
                        sources_location.add("faces");
                    }
                  else if( ( sub_type( MD_Vector_composite, md.valeur( ) ) ) && (localisation_=="??"))
                    {
                      Cerr << "Error in Champ_Generique_Transformation::completer "<<finl;
                      Cerr << "The source number "<<i<<" is composite and the location was not provided. It is forbidden to apply a 'transformation' on a composite source without providing a location. "<<finl;
                      Cerr << "You must perform an interpolation before or specify a location."<<finl;
                      Cerr << "Aborting..."<<finl;
                      Process::abort( );
                      // if we want to deal with composite case we can do something like...
                      // MD_Vector_composite md_comp = ref_cast( MD_Vector_composite, values_source_i.get_md_vector( ).valeur( ) );
                      // nb_parts = md_comp.nb_parts( ); //to perform loop over parts and then get the considered part :
                      // md = md_comp.get_desc_part( part );
                      // when dealing with pression_pa it leads to have elem and som locations
                    }
                  else if ( md == zvf_source_i.face_sommets( ).get_md_vector( ) )
                    {
                      sources_location.add( "faces" );
                    }
                  else if( md == zvf_source_i.domaine( ).les_elems( ).get_md_vector( ) )
                    {
                      sources_location.add( "elem" );
                    }
                  else if( md == zvf_source_i.domaine( ).les_sommets( ).get_md_vector( ) )
                    {
                      sources_location.add( "som" );
                    }
                  else if( md == zvf_source_i.xa( ).get_md_vector( ) )
                    {
                      sources_location.add( "xa" );
                    }
                  else
                    {
                      Cerr << "Warning in Champ_Generique_Transformation::completer "<<finl;
                      Cerr << "Impossible to determine the localisation of the source number "<<i<<" using md_vectors so try using directive"<<finl;
                      use_directive = true;
                    }
                  // Cout << "sources location "<<i<<" : "<<sources_location( i )<<finl;
                }//end of not using directive case
            }//end of has domaine dis
          else
            {
              //in the case where there is no domaine dis base, we look at the directive
              use_directive=true;
            }

          if( use_directive )
            {
              Motcle directive = get_source( i ).get_directive_pour_discr( );
              if (directive=="champ_elem")
                sources_location.add( "elem" );
              else if (directive=="champ_sommets")
                sources_location.add( "som" );
              else if (directive=="champ_face")
                sources_location.add( "faces" );
              else if (directive=="pression")
                {
                  if  (localisation_=="??")
                    {
                      Cerr << "Error in Champ_Generique_Transformation::completer"<<finl;
                      Cerr << "The directive associated to the source number "<<i<<" is pressure and the location was not provided.";
                      Cerr <<" It is forbidden to apply a 'transformation' on a composite source without providing a location. "<<finl;
                      Cerr << "You must perform an interpolation before or specify a location." <<finl;
                      Cerr << "Aborting..."<<finl;
                      Process::abort( );
                    }
                  sources_location.add( "unknown" );
                }
              else if ( directive == "champ_uniforme" )
                sources_location.add( "elem" );
              else if ( directive == "champ_don" )
                sources_location.add( "elem" );
              else
                {
                  Cerr<<"Error in Champ_Generique_Transformation::completer"<<finl;
                  Cerr<<"The discretization of the first source of the generic field "<<nom_post_<<finl;
                  Cerr<<"does not correspond to any TRUST field discretization."<<finl;
                  Cerr<<directive<<finl;
                  Cerr<<"Aborting..."<<finl;
                  Process::abort( );
                }
            }
        }//loop over sources
      const int nb_locations = sources_location.size( );
      if( nb_sources != nb_locations )
        {
          Cerr << "Error in Champ_Generique_Transformation::completer "<<finl;
          Cerr << nb_sources<<" were detected but "<<nb_locations<<" were found"<<finl;
          Cerr << "Aborting..."<<finl;
          Process::abort( );
        }
      Nom reference_location = sources_location[0];
      for(int i=1; i<nb_locations; i++)
        {
          if( sources_location[i] != reference_location )
            {
              sources_location_ok = false;
              Cerr << "Warning in Champ_Generique_Transformation::completer "<<finl;
              Cerr << "All sources have location : "<<finl;
              Cerr <<sources_location;
              Cerr <<"But the specified location is "<<localisation_<<finl;
              //Cerr << "Aborting..."<<finl;
              //Process::abort( );
            }
        }
    } //end dealing with nb_source_fictive <= 0

  if( localisation_ == "??" )
    {
      if( sources_location_ok  )
        {
          Cerr <<"sources have same location : "<< sources_location[0]<<finl;
        }
      else
        {
          Cerr << "Warning sources have not the same location, taking "<<sources_location[0]<<" as location"<<finl;
        }
      localisation_ = sources_location[0] ;
    }
  else
    {
      if( ! fictive_source_ && localisation_ != sources_location[0] )
        {
          Cerr << "Warning first source is located to "<<sources_location[0]<<" but the user has specified "<< localisation_<<finl;
        }
    }

  verifier_localisation();
  preparer_macro();

  if (Motcle(methode_)=="formule")
    {
      // we re-query nb_sources because sources may have been added by completer
      nb_sources = get_nb_sources();
      fxyz.dimensionner(1);
      fxyz[0].setNbVar(4+nb_sources);
      fxyz[0].setString(les_fct[0]);
      fxyz[0].addVar("x");
      fxyz[0].addVar("y");
      fxyz[0].addVar("z");
      fxyz[0].addVar("t");
      for (int i=0; i<nb_sources; i++)
        {
          const Noms nom = get_source(i).get_property("nom");
          fxyz[0].addVar(nom[0]);
          Journal()<<" Name of source "<<nom[0]<<finl;
        }
      fxyz[0].parseString();
    }
}


void projette(DoubleTab& valeurs_espace,const DoubleTab& val_source,const Domaine_VF& zvf,const Motcle& loc,bool champ_normal_faces)
{
  const Nom& type_elem=zvf.domaine().type_elem()->que_suis_je();
  if ((champ_normal_faces)||(loc!="elem")|| (type_elem!="Quadrangle"))
    {
      Cerr<<"option composante_normale not coded in this case"<<finl;
      Process::exit();
    }
  const DoubleTab& coord=zvf.domaine().coord_sommets();
  const IntTab& elems=zvf.domaine().les_elems();

  int nb_v= valeurs_espace.dimension(0);
  assert(nb_v==zvf.nb_elem());
  assert(val_source.dimension(1)==Objet_U::dimension);
  ToDo_Kokkos("critical");
  for (int i=0; i<nb_v; i++)
    {
      Vecteur3 v1,v2,normal,val;
      int s0=elems(i,0);
      int s1=elems(i,1);
      int s2=elems(i,2);
      if (Objet_U::dimension!=3)
        {
          Cerr<<"not implemented "<<finl;
          assert(0);
          Process::exit();
        }
      v1.set(coord(s1,0)-coord(s0,0),coord(s1,1)-coord(s0,1),coord(s1,2)-coord(s0,2));
      v2.set(coord(s2,0)-coord(s0,0),coord(s2,1)-coord(s0,1),coord(s2,2)-coord(s0,2));
      val.set(val_source(i,0),val_source(i,1),val_source(i,2));

      Vecteur3::produit_vectoriel(v1,v2,normal);
      double norm=Vecteur3::produit_scalaire(normal,normal);
      double inv=1./sqrt(norm);
      normal*=inv;


      double prod_scal= Vecteur3::produit_scalaire(normal,val);
      valeurs_espace(i)=prod_scal;
    }

}
//Construction of the storage space returned by this post-processing field
//-Characterization of the storage space returned by this post-processing field (typing ...)
//-Construction of the positions array containing the coordinates of the computation points (support)
//-Interpolation of source values on the retained support
//-Evaluation of the storage space values as a function of (les_fct)

const Champ_base& Champ_Generique_Transformation::get_champ_without_evaluation(OWN_PTR(Champ_base)&) const
{
  if (!espace_stockage_)
    creer_espace_stockage(nature_ch,nb_comp_,espace_stockage_);
  else
    espace_stockage_->changer_temps(get_source(0).get_time());
  return espace_stockage_;
}
const Champ_base& Champ_Generique_Transformation::get_champ(OWN_PTR(Champ_base)&) const
{
  const Domaine_dis_base& domaine_dis = get_ref_domaine_dis_base();
  if (!espace_stockage_)
    creer_espace_stockage(nature_ch,nb_comp_,espace_stockage_);
  else
    espace_stockage_->changer_temps(get_source(0).get_time());
  DoubleTab& valeurs_espace = espace_stockage_->valeurs();
  const Domaine_VF& zvf = ref_cast(Domaine_VF,domaine_dis);
  int nb_elem_tot = zvf.nb_elem_tot();
  int nb_som_tot = get_ref_domain().nb_som_tot();
  const Motcle directive = get_directive_pour_discr();
  bool champ_normal_faces = 0;
  if (get_source(0).get_discretisation().is_vdf() || get_source(0).get_discretisation().is_poly_family())
    champ_normal_faces = 1;

  //Construction of the positions array containing the coordinates
  //of the points (support) where the values of the storage
  //space will be evaluated
  DoubleTrav positions;
  if (localisation_ == "elem")
    {
      if (zvf.xp().nb_dim() != 2) /* xp() not initialized */
        zvf.domaine().calculer_centres_gravite(positions);
      else
        zvf.get_position(positions);
    }
  else if (localisation_ == "som")
    positions = get_ref_domain().coord_sommets();
  else if (localisation_ == "faces")
    positions = zvf.xv();
  else if (localisation_ == "elem_som")
    {
      const DoubleTab& som = get_ref_domain().coord_sommets();
      int dim = som.dimension(1);
      positions = zvf.xp();
      positions.resize(nb_elem_tot+nb_som_tot,dim);
      ToDo_Kokkos("critical loop");
      for (int i1=nb_elem_tot; i1<nb_elem_tot+nb_som_tot; i1++)
        for (int i2=0; i2<dim; i2++)
          positions(i1,i2) = som(i1-nb_elem_tot,i2);
    }
  else
    {
      Cerr<<"Problem in " <<que_suis_je()<<finl;
      Cerr<<"Item "<<localisation_<<" is not recognized."<<finl;
      exit();
    }

  //Interpolation of source values on the retained support and storage in sources_val

  int nb_sources = get_nb_sources();
  int nb_pos = positions.dimension(0);
  assert(nb_pos>0||(methode_=="formule")||(methode_=="composante_normale")||(methode_=="vecteur"));
  using DoubleTravs = TRUST_Vector<DoubleTrav>; // remplace VECT(DoubleTrav)
  const int max_nb_sources = 14;
  if (nb_sources>max_nb_sources)
    {
      Cerr << "Increase max_nb_sources to " << nb_sources << " in Champ_base& Champ_Generique_Transformation::get_champ() !" << finl;
      Process::exit();
    }
  DoubleTravs sources_val(nb_sources);
  IntTrav nb_comps(nb_sources);
  Noms nom_source(nb_sources);
  int dim_compo = 2*dimension;
  Noms compo(dim_compo);
  double temps;
  temps = get_time();
  int nb_compso = 0; // parameter offset so that it is recognized by all methods under test (see the composante method below)
  for (int so=0; so<nb_sources; so++)
    {
      const Champ_Generique_base& source = get_source(so);
      OWN_PTR(Champ_base) stockage_so;
      const Champ_base& source_so = source.get_champ(stockage_so);
      const DoubleTab& source_so_val = source_so.valeurs();
      const Motcle directive_so = source.get_directive_pour_discr();

      //int nb_compso = source_so.nb_comp();
      nb_compso = source_so.nb_comp();
      nb_comps(so) = nb_compso;
      const Noms nom_champ = source.get_property("nom");
      nom_source[so] = nom_champ[0];
      if ((Motcle(methode_)=="produit_scalaire") || (Motcle(methode_)=="norme"))
        {
          const Noms compo_ch = source.get_property("composantes");
          for (int comp=0; comp<dimension; comp++)
            compo[so*dimension+comp] = compo_ch[comp];
        }

      //If VDF and processing a dot product or a norm, then interpolate
      //to retrieve an array with as many components as the problem dimension
      if (directive!=directive_so)
        {
          sources_val[so].resize(nb_pos,nb_compso);
          if (directive == "CHAMP_FACE") source_so.valeur_aux_faces(sources_val[so]);
          else if (directive == "CHAMP_ELEM" && nb_pos==domaine_dis.domaine().nb_elem()) source_so.valeur_aux_centres_de_gravite(domaine_dis.domaine(), sources_val[so]);
          else if (directive == "CHAMP_FONC_QUAD_DG" and directive_so == "CHAMP_ELEM_DG")
            {
              int nelem = valeurs_espace.dimension(0);
              int npoints = valeurs_espace.dimension(1);
              if (methode_=="vecteur") npoints /= dimension; //TODO DG will be cleaner when quadrature is in Kernel
              sources_val[so].resize(nelem,npoints);
              if (!fictive_source_) source_so.eval_elem(sources_val[so]);
            }
          else source_so.valeur_aux(positions,sources_val[so]);
        }
      else
        {
          if (source_so_val.line_size()==1)
            {
              int nn=source_so_val.dimension_tot(0);
              sources_val[so].resize(nn,1);
              sources_val[so] = source_so_val;
            }
          else
            sources_val[so] = source_so_val;
        }

      if ((nb_pos==0)&&(valeurs_espace.dimension_tot(0)!=0))
        {
          int nb_pos_tmp=valeurs_espace.dimension(0);
          if (nb_pos_tmp > sources_val[so].size())
            {
              Cerr << "\nError in Champ_Generique_Transformation::get_champ" << finl;
              Cerr << nom_source[so] << " : Wrong number of elements." << finl;
              exit();
            }
        }
    }
  //End of filling sources_val

  //Computation of storage space values according to the selected methode_
  if ((Motcle(methode_)=="produit_scalaire") || (Motcle(methode_)=="norme"))
    {
      ParserView parser(fxyz[0]);
      parser.parseString();
      int dim = dimension;
      Kokkos::Array<CDoubleTabView, max_nb_sources> sources;
      for (int so=0; so<nb_sources; so++)
        sources[so] = sources_val[so].view_ro();
      DoubleArrView valeurs = static_cast<ArrOfDouble&>(valeurs_espace).view_wo();
      if (directive == "champ_fonc_quad_dg") //This is for DG
        {
          ToDo_Kokkos("Code but check test!");
          int nb_elem = valeurs_espace.dimension(0);
          IntTab nb_points, ind_integ_points;
          zvf.get_ind_integ_points(ind_integ_points);
          zvf.get_nb_integ_points(nb_points);
          int nb_pts_integ_max = zvf.get_max_nb_integ_points();
          CIntArrView ind_integ_points_w = static_cast<const ArrOfInt&>(ind_integ_points).view_ro();
          CIntArrView nb_points_w = static_cast<const ArrOfInt&>(nb_points).view_ro();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_elem, KOKKOS_LAMBDA(const int i)
          {
            for (int pt=0; pt<nb_points_w(i); pt++)
              {
                int threadId = parser.acquire();
                int k = ind_integ_points_w(i)+pt;
                for (int so=0; so<nb_sources; so++)
                  for (int j=0; j<dim; j++)
                    parser.setVar(so*dim+j,sources[so](i,j*nb_pts_integ_max+pt), threadId);
                valeurs(k) = parser.eval(threadId);
                parser.release(threadId);
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else
        {
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_pos, KOKKOS_LAMBDA(const int i)
          {
            int threadId = parser.acquire();
            for (int so=0; so<nb_sources; so++)
              for (int j=0; j<dim; j++)
                parser.setVar(so*dim+j,sources[so](i,j), threadId);
            valeurs(i) = parser.eval(threadId);
            parser.release(threadId);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
  else if (Motcle(methode_)=="vecteur")
    {
      if ((champ_normal_faces) && (localisation_=="faces"))
        {
          assert(nb_comp_==dimension);
          int dim = dimension;
          ParserView fxyz0(fxyz[0]);
          ParserView fxyz1(fxyz[1]);
          ParserView fxyz2(fxyz[dim-1]); // Trick to support 2D/3D
          fxyz0.parseString();
          fxyz1.parseString();
          if (dim==3) fxyz2.parseString();
          Kokkos::Array<CDoubleTabView, max_nb_sources> sources;
          for (int so=0; so<nb_sources; so++)
            sources[so] = sources_val[so].view_ro();
          CDoubleTabView face_normales = zvf.face_normales().view_ro();
          CDoubleArrView surface = zvf.face_surfaces().view_ro();
          CDoubleTabView pos = positions.view_ro();
          //CIntArrView nb_comp_sources = static_cast<const ArrOfInt&>(nb_comps).view_ro();
          DoubleArrView valeurs = static_cast<ArrOfDouble&>(valeurs_espace).view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_pos, KOKKOS_LAMBDA(const int i)
          {
            int threadId = fxyz0.acquire();
            for (int j=0; j<dim; j++)
              {
                double pos_i_j = pos(i, j);
                fxyz0.setVar(j, pos_i_j, threadId);
                fxyz1.setVar(j, pos_i_j, threadId);
                if (dim==3) fxyz2.setVar(j, pos_i_j, threadId);
              }
            fxyz0.setVar(3, temps, threadId);
            fxyz1.setVar(3, temps, threadId);
            if (dim==3) fxyz2.setVar(3, temps, threadId);
            for (int so=0; so<nb_sources; so++)
              {
                double var = sources[so](i,0);
                fxyz0.setVar(so+4,var, threadId);
                fxyz1.setVar(so+4,var, threadId);
                if (dim==3) fxyz2.setVar(so+4,var, threadId);
              }
            valeurs(i)  = fxyz0.eval(threadId) * face_normales(i, 0) / surface(i);
            valeurs(i) += fxyz1.eval(threadId) * face_normales(i, 1) / surface(i);
            if (dim==3) valeurs(i) += fxyz2.eval(threadId) * face_normales(i, 2) / surface(i);
            fxyz0.release(threadId);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else
        {
          if (directive=="champ_fonc_quad_dg") //This is for DG
            {
              ToDo_Kokkos("critical");
              IntTab nb_points, ind_integ_points;
              int nb_elem = valeurs_espace.dimension(0);
              zvf.get_ind_integ_points(ind_integ_points);
              zvf.get_nb_integ_points(nb_points);
              int nb_pts_integ_max = zvf.get_max_nb_integ_points();

              for (int i=0; i<nb_elem; i++)
                {
                  for (int pt=0; pt<nb_points[i]; pt++)
                    {
                      int k = ind_integ_points[i]+pt;
                      double x = positions(k,0);
                      double y = positions(k,1);
                      double z = (dimension>2 ? positions(k,2) : 0);

                      for (int j=0; j<nb_comp_; j++)
                        {
                          fxyz[j].setVar(0,x);
                          fxyz[j].setVar(1,y);
                          fxyz[j].setVar(2,z);
                          fxyz[j].setVar(3,temps);
                          for (int so=0; so<nb_sources; so++)
                            {
                              const DoubleTab& source_so_val = sources_val[so];
                              fxyz[j].setVar(so+4,source_so_val(i,pt));
                            }
                          int l = nb_pts_integ_max*j + pt;
                          valeurs_espace(i,l) = fxyz[j].eval();
                        }
                    }
                }
            }
          else
            {
              Kokkos::Array<CDoubleTabView, max_nb_sources> sources;
              for (int so=0; so<nb_sources; so++)
                sources[so] = sources_val[so].view_ro();
              int dim = dimension;
              CDoubleTabView pos = positions.view_ro();
              DoubleTabView valeurs = valeurs_espace.view_wo();
              for (int j=0; j<nb_comp_; j++)
                {
                  ParserView fxyzj(fxyz[j]);
                  fxyzj.parseString();
                  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_pos, KOKKOS_LAMBDA(const int i)
                  {
                    double x = pos(i,0);
                    double y = pos(i,1);
                    double z = (dim>2 ? pos(i,2) : 0);
                    int threadId = fxyzj.acquire();
                    fxyzj.setVar(0,x,threadId);
                    fxyzj.setVar(1,y,threadId);
                    fxyzj.setVar(2,z,threadId);
                    fxyzj.setVar(3,temps,threadId);
                    for (int so=0; so<nb_sources; so++)
                      fxyzj.setVar(so+4,sources[so](i,0),threadId);
                    valeurs(i,j) = fxyzj.eval(threadId);
                    fxyzj.release(threadId);
                  });
                  end_gpu_timer(__KERNEL_NAME__);
                }
            }
        }
    }
  else if (Motcle(methode_)=="composante_normale")
    {
      const DoubleTab& source_so_val = sources_val[0];
      projette(valeurs_espace,source_so_val,zvf,localisation_,champ_normal_faces);
    }
  else if (Motcle(methode_)=="composante")
    {
      int num_compo = les_fct.size()-1;
      if ((num_compo<0)||(num_compo>=nb_compso))
        {
          Cerr<<"Component number must be >=0 and < "<<nb_compso<<" for composante option."<<finl;
          exit();
        }

      CDoubleTabView source_so_val = sources_val[0].view_ro();
      DoubleArrView v = static_cast<ArrOfDouble&>(valeurs_espace).view_wo();
      if ((champ_normal_faces) && (localisation_=="faces"))
        {
          CDoubleTabView face_normales = zvf.face_normales().view_ro();
          CDoubleArrView surface = zvf.face_surfaces().view_ro();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_pos, KOKKOS_LAMBDA(const int i)
          {
            v(i) = source_so_val(i, 0) * face_normales(i, num_compo) / surface(i);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else
        {
          if (directive=="champ_fonc_quad_dg") //This is for DG
            {
              ToDo_Kokkos("Code but check test!");
              IntTab nb_points, ind_integ_points;
              int nb_elem = valeurs_espace.dimension(0);
              zvf.get_ind_integ_points(ind_integ_points);
              zvf.get_nb_integ_points(nb_points);
              int nb_pts_integ_max = zvf.get_max_nb_integ_points();
              CIntArrView ind_integ_points_w = static_cast<const ArrOfInt&>(ind_integ_points).view_ro();
              CIntArrView nb_points_w = static_cast<const ArrOfInt&>(nb_points).view_ro();

              Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_elem, KOKKOS_LAMBDA(const int i)
              {
                for (int pt=0; pt<nb_points_w[i]; pt++)
                  {
                    int k = ind_integ_points_w[i]+pt;
                    v(k) = source_so_val(i, num_compo*nb_pts_integ_max+pt);
                  }
              });
              end_gpu_timer(__KERNEL_NAME__);
            }
          else
            {
              Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_pos, KOKKOS_LAMBDA(const int i)
              {
                v(i) = source_so_val(i, num_compo);
              });
              end_gpu_timer(__KERNEL_NAME__);
            }
        }
    }
  else if (Motcle(methode_)=="formule")
    {
      //Evaluation of the storage space values as a function of (les_fct)
      ParserView parser(fxyz[0]);
      parser.parseString();
      int special=0;
      if (nb_pos==0 && valeurs_espace.dimension_tot(0)!=0)
        {
          special = 1;
          nb_pos = valeurs_espace.dimension(0);
        }

      if (directive == "champ_fonc_quad_dg") //This is for DG
        {
          ToDo_Kokkos("Code but check test!");
          int dim = dimension;
          int nb_elem = valeurs_espace.dimension(0);
          Kokkos::Array<CDoubleTabView, max_nb_sources> sources;
          for (int so=0; so<nb_sources; so++)
            sources[so] = sources_val[so].view_ro();
          IntTab nb_points, ind_integ_points;
          zvf.get_ind_integ_points(ind_integ_points);
          zvf.get_nb_integ_points(nb_points);
          int nb_pts_integ_max = zvf.get_max_nb_integ_points();

          CIntArrView ind_integ_points_w = static_cast<const ArrOfInt&>(ind_integ_points).view_ro();
          CIntArrView nb_points_w = static_cast<const ArrOfInt&>(nb_points).view_ro();
          CDoubleTabView pos = positions.view_ro();
          DoubleTabView valeurs = valeurs_espace.view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_elem, KOKKOS_LAMBDA(const int i)
          {
            for (int pt=0; pt<nb_points_w(i); pt++)
              {
                int threadId = parser.acquire();
                int k = ind_integ_points_w(i)+pt;
                double x = special ? 1e38 : pos(k,0);
                double y = special ? 1e38 : pos(k,1);
                double z = special ? 1e38 : (dim>2 ? pos(k,2) : 0);
                parser.setVar(0,x,threadId);
                parser.setVar(1,y,threadId);
                parser.setVar(2,z,threadId);
                parser.setVar(3,temps,threadId);

                for (int d = 0; d < nb_comp_; d++)
                  {
                    int j = nb_pts_integ_max*d + pt;
                    for (int so=0; so<nb_sources; so++)
                      parser.setVar(so+4,sources[so](i,j),threadId);
                    valeurs(i, j) = parser.eval(threadId);
                  }
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else
        {
          int line_size = valeurs_espace.line_size();
          int dim = dimension;
          int nb_comp = nb_comp_;
          Kokkos::Array<CDoubleTabView, max_nb_sources> sources;
          for (int so=0; so<nb_sources; so++)
            sources[so] = sources_val[so].view_ro();
          CDoubleTabView pos = positions.view_ro();
          CIntArrView nb_comp_sources = static_cast<const ArrOfInt&>(nb_comps).view_ro();
          DoubleTabView valeurs = valeurs_espace.view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_pos, KOKKOS_LAMBDA(const int i)
          {
            double x = special ? 1e38 : pos(i,0);
            double y = special ? 1e38 : pos(i,1);
            double z = special ? 1e38 : (dim>2 ? pos(i,2) : 0);
            int threadId = parser.acquire();
            parser.setVar(0,x,threadId);
            parser.setVar(1,y,threadId);
            parser.setVar(2,z,threadId);
            parser.setVar(3,temps,threadId);
            if (line_size == 1)
              {
                for (int so=0; so<nb_sources; so++)
                  parser.setVar(so+4,sources[so](i,0),threadId);
                valeurs(i, 0) = parser.eval(threadId);
              }
            else // line_size > 1
              {
                for (int j=0; j<nb_comp; j++)
                  {
                    for (int so=0; so<nb_sources; so++)
                      {
                        int nbcomp_loc = nb_comp_sources(so);
                        if (nbcomp_loc == nb_comp)
                          parser.setVar(so+4,sources[so](i,j),threadId);
                        else if (nbcomp_loc == 1)
                          parser.setVar(so+4,sources[so](i,0),threadId);
                        else
                          Process::Kokkos_exit("The arrays of values don't have compatibles dimensions in Champ_Generique_Transformation::get_champ()");
                      }
                    valeurs(i,j) = parser.eval(threadId);
                  }
              }
            parser.release(threadId);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
  // PL: Removal of a costly synchronization that is very often unnecessary
  // See Champ_Generique_Interpolation (localization = som) for the deferred echange_espace_virtuel
  // valeurs_espace.echange_espace_virtuel();
  return espace_stockage_;
}

//The elem, som and faces localizations can be retained for post-processing
//The elem_som localization is used only for evaluating the values of the storage space
//but cannot be retained as a post-processing localization

Entity Champ_Generique_Transformation::get_localisation(const int index) const
{
  Entity loc;
  if (localisation_=="elem")
    loc = Entity::ELEMENT;
  else if (localisation_=="som")
    loc = Entity::NODE;
  else if (localisation_=="faces")
    loc = Entity::FACE;
  else
    {
      Cerr << "Error of type : localisation should be specified to elem or som or faces for the defined field " << nom_post_ << finl;
      throw Champ_Generique_erreur("INVALID");
    }
  return loc;
}

const Noms Champ_Generique_Transformation::get_property(const Motcle& query) const
{
  Motcles motcles(2);
  motcles[0] = "composantes";
  motcles[1] = "unites";
  int rang = motcles.search(query);
  switch(rang)
    {
    case 0:
      {
        Noms compo(nb_comp_);
        for (int i=0; i<nb_comp_; i++)
          {
            Nom nume(i);
            compo[i] = nom_post_+nume;
          }
        return compo;
      }
    case 1:
      {
        Noms unites(nb_comp_);
        for (int comp=0; comp<nb_comp_; comp++)
          unites[comp] = unite_;
        return unites;
      }
    }
  return Champ_Gen_de_Champs_Gen::get_property(query);
}

const Motcle Champ_Generique_Transformation::get_directive_pour_discr() const
{
  Motcle directive;
  if (localisation_=="elem")
    {
      OWN_PTR(Champ_base) source_espace_stockage;
      const Champ_base& source = get_source(0).get_champ_without_evaluation(source_espace_stockage);
      directive = (source.is_basis_function() || source.is_quadrature()) ? "champ_fonc_quad_dg" : "champ_elem";
    }
  else if (localisation_=="som")
    {
      directive = "champ_sommets";
    }
  else if (localisation_=="faces")
    {
      directive = "champ_face";
    }
  else if (localisation_=="elem_som")
    {
      directive = "pression";
    }
  else
    {
      Cerr<<"Localisation "<<localisation_<<" is not recognized by "<<que_suis_je()<<finl;
      exit();
    }
  return directive;
}

//Names the field as default source
//"Combinaison_"+"les_fct"
void Champ_Generique_Transformation::nommer_source()
{
  if (nom_post_=="??")
    {
      Nom nom_post_source;
      nom_post_source =  "Combinaison_";

      // PDI doesn't understand mathematical symbol inside variable name, so we replace them
      std::string fonction = les_fct[0].getString();
      auto replace = [&](std::string oldS, std::string newS)
      {
        size_t start_pos = 0;
        while((start_pos = fonction.find(oldS, start_pos)) != std::string::npos)
          {
            fonction.replace(start_pos, oldS.length(), newS);
            start_pos += newS.length();
          }
      };
      replace("*", "x");
      replace("+", "PLUS");
      replace("-", "MINUS");
      replace("/", "DIV");
      replace("(", "_BEGINPAR_");
      replace(")", "_ENDPAR_");
      nom_post_source += fonction;
      nommer(nom_post_source);
    }
}

//Returns the number of components to manipulate in the expression
//Verifies the consistency of the data for the chosen method
int Champ_Generique_Transformation::preparer_macro()
{
  int nb_sources = get_nb_sources();
  int nb_var = 0;
  int erreur = 0;
  Nom msg;

  if ((Motcle(methode_)=="produit_scalaire") || (Motcle(methode_)=="norme"))
    {
      for (int i=0; i<nb_sources; i++)
        {
          OWN_PTR(Champ_base) source_espace_stockage;
          const Champ_base& source = get_source(i).get_champ(source_espace_stockage);
          int nb_comp = source.nb_comp();
          if (!source.is_vectorial() && source.nature_du_champ()!=multi_scalaire)
            {
              msg = "At least one of the source fields is not of vector nature.";
              erreur = 1;
            }
          nb_var += nb_comp;
        }

      if ((Motcle(methode_)=="produit_scalaire") && (nb_sources!=2))
        {
          msg = "Exactly two vector source fields must be specified.";
          erreur = 1;
        }
      else if ((Motcle(methode_)=="norme") && (nb_sources!=1))
        {
          msg = "Exactly one vector source field must be specified.";
          erreur = 1;
        }
      // the rest is done if there is no error
      if ( erreur==0 )
        {
          nb_comp_ = 1;
          nature_ch = scalaire;
          // The call to creer_expression_macro is done if we have a vector because it uses the components!
          creer_expression_macro();
          fxyz.dimensionner(1);
          fxyz[0].setNbVar(nb_var);
          fxyz[0].setString(les_fct[0]);
          for (int i=0; i<nb_sources; i++)
            {
              OWN_PTR(Champ_base) source_espace_stockage;
              const Champ_base& source = get_source(i).get_champ(source_espace_stockage);
              int nb_comp = source.nb_vect_comp();
              const Noms compo = get_source(i).get_property("composantes");
              for (int comp=0; comp<nb_comp; comp++)
                {
                  fxyz[0].addVar(compo[comp]);
                  Journal()<<" Source compo "<<compo[comp]<<finl;
                }
            }
          fxyz[0].parseString();
        }
    }
  else if (Motcle(methode_)=="vecteur")
    {
      nb_comp_ = les_fct.size();
      nature_ch = vectoriel;
      fxyz.dimensionner(nb_comp_);

      for (int comp=0; comp<nb_comp_; comp++)
        {
          fxyz[comp].setNbVar(4+nb_sources);
          fxyz[comp].setString(les_fct[comp]);
          fxyz[comp].addVar("x");
          fxyz[comp].addVar("y");
          fxyz[comp].addVar("z");
          fxyz[comp].addVar("t");
          for (int i=0; i<nb_sources; i++)
            {
              const Noms nom = get_source(i).get_property("nom");
              fxyz[comp].addVar(nom[0]);
              Journal()<<" Name of source "<<nom[0]<<finl;
            }
          fxyz[comp].parseString();
        }
    }
  else if ((Motcle(methode_)=="composante")|| ((Motcle(methode_)=="composante_normale") ))
    {
      if (nb_sources!=1)
        {
          msg = "Exactly one vector source field must be specified.";
          erreur = 1;
        }
      // The rest is done if there is no error !!
      if (erreur==0) // 1 source, no need to loop
        {
          OWN_PTR(Champ_base) source_espace_stockage;
          const Champ_base& source = get_source(0).get_champ_without_evaluation(source_espace_stockage);
          if (source.nature_du_champ()!=vectoriel && source.nature_du_champ()!=multi_scalaire)
            {
              msg = "The source field is not of vector nature.";
              erreur = 1;
            }
        }
      nb_comp_ = 1;
      nature_ch = scalaire;
    }
  else if (Motcle(methode_)=="formule")
    nb_var = -1;

  if ((nature_ch==scalaire) && (localisation_!="elem"))
    {
      if (get_source(0).get_discretisation().is_vdf())
        {
          msg = "The nature of the generated storing field will be scalar.\n";
          msg += "In that case the location elem is required since VDF discretization is considered.";
          erreur = 1;
        }
    }

  if (erreur)
    {
      Cerr<<"Please check the syntax of the Transformation generic field "<<nom_post_<<finl;
      Cerr<<"for wich the selected option is : "<<methode_<<finl;
      Cerr<<msg<<finl;
      exit();
    }
  return nb_var;
}

//Creates expressions for specific operations:
//dot product and norm
void Champ_Generique_Transformation::creer_expression_macro()
{
  if ((Motcle(methode_)=="produit_scalaire") || (Motcle(methode_)=="norme"))
    {
      const Noms compo0 = get_source(0).get_property("composantes");

      if (Motcle(methode_)=="norme")
        {
          les_fct[0] = "sqrt(";
          les_fct[0] += compo0[0];
          les_fct[0] += "*";
          les_fct[0] += compo0[0];
          les_fct[0] += "+";
          les_fct[0] += compo0[1];
          les_fct[0] += "*";
          les_fct[0] += compo0[1];
          if (dimension==3)
            {
              les_fct[0] += "+";
              les_fct[0] += compo0[2];
              les_fct[0] += "*";
              les_fct[0] += compo0[2];
            }
          les_fct[0] +=")";
        }

      if (Motcle(methode_)=="produit_scalaire")
        {
          const Noms compo1 = get_source(1).get_property("composantes");
          les_fct[0] = compo0[0];
          les_fct[0] += "*";
          les_fct[0] += compo1[0];
          les_fct[0] += "+";
          les_fct[0] += compo0[1];
          les_fct[0] += "*";
          les_fct[0] += compo1[1];
          if (dimension==3)
            {
              les_fct[0] += "+";
              les_fct[0] += compo0[2];
              les_fct[0] += "*";
              les_fct[0] += compo1[2];
            }
        }
    }
}
