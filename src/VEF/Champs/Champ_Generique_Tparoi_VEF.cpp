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

#include <Champ_Generique_Tparoi_VEF.h>
#include <Probleme_base.h>
#include <Nom.h>
#include <TRUSTTabs.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Equation_base.h>
#include <Modele_turbulence_scal_base.h>
#include <Neumann_paroi.h>
#include <Neumann_homogene.h>
#include <Discretisation_base.h>
#include <Synonyme_info.h>

Implemente_instanciable(Champ_Generique_Tparoi_VEF,"Tparoi_VEF",Champ_Gen_de_Champs_Gen);
// XD tparoi_vef champ_post_de_champs_post tparoi_vef INHERITS_BRACE This keyword is used to post process (only for VEF
// XD_CONT discretization) the temperature field with a slight difference on boundaries with Neumann condition where law
// XD_CONT of the wall is applied on the temperature field. nom_pb is the problem name and field_name is the selected
// XD_CONT field name. A keyword (temperature_physique) is available to post process this field without using
// XD_CONT Definition_champs.

Add_synonym(Champ_Generique_Tparoi_VEF,"Champ_Post_Tparoi_VEF");

Sortie& Champ_Generique_Tparoi_VEF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

//cf Champ_Gen_de_Champs_Gen::readOn
Entree& Champ_Generique_Tparoi_VEF::readOn(Entree& s )
{
  Champ_Gen_de_Champs_Gen::readOn(s);
  return s ;
}

const Champ_base& Champ_Generique_Tparoi_VEF::get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const
{

  OWN_PTR(Champ_base) source_espace_stockage;
  const Champ_base& source = get_source(0).get_champ_without_evaluation(source_espace_stockage);

  const Domaine_dis_base& domaine_dis = get_ref_domaine_dis_base();
  Nature_du_champ nature_source = source.nature_du_champ();
  Noms noms;
  Noms unites;
  noms.add("bidon");
  unites.add("K");
  int nb_comp = 1;
  double temps;
  temps=0.;

  OWN_PTR(Champ_Fonc_base)  espace_stockage_fonc;
  const Discretisation_base&  discr = get_discretisation();
  Motcle directive = get_directive_pour_discr();
  discr.discretiser_champ(directive,domaine_dis,nature_source,noms,unites,nb_comp,temps,espace_stockage_fonc);
  espace_stockage = espace_stockage_fonc;


  return espace_stockage;
}
const Champ_base& Champ_Generique_Tparoi_VEF::get_champ(OWN_PTR(Champ_base)& espace_stockage) const
{

  OWN_PTR(Champ_base) source_espace_stockage;
  const Champ_base& source = get_source(0).get_champ(source_espace_stockage);

  const Domaine_dis_base& domaine_dis = get_ref_domaine_dis_base();
  Nature_du_champ nature_source = source.nature_du_champ();
  Noms noms;
  Noms unites;
  noms.add("bidon");
  unites.add("K");
  int nb_comp = 1;
  double temps;
  temps=0.;

  OWN_PTR(Champ_Fonc_base)  espace_stockage_fonc;
  const Discretisation_base&  discr = get_discretisation();
  Motcle directive = get_directive_pour_discr();
  discr.discretiser_champ(directive,domaine_dis,nature_source,noms,unites,nb_comp,temps,espace_stockage_fonc);
  espace_stockage = espace_stockage_fonc;

  DoubleTab& valeurs_espace = espace_stockage->valeurs();
  const DoubleTab& source_so_val = source.valeurs();
  const DoubleTab& inconnue = source_so_val;
  valeurs_espace = source_so_val;

  // For all faces where flux is imposed, if the turbulent flux is not computed, recompute the temperature
  // retrieve the temperature equation
  const Equation_base& my_eqn=ref_cast(Champ_Inc_base,source).equation();
  const RefObjU& modele_turbulence = my_eqn.get_modele(TURBULENCE);
  if ( modele_turbulence && sub_type(Modele_turbulence_scal_base,modele_turbulence.valeur()))
    {
      const Modele_turbulence_scal_base& mod_turb_scal = ref_cast(Modele_turbulence_scal_base,modele_turbulence.valeur());
      const Turbulence_paroi_scal_base& loiparth = mod_turb_scal.loi_paroi();

      if( loiparth.use_equivalent_distance() )
        {
          // const Paroi_scal_hyd_base_VEF& paroi_scal_vef = ref_cast(Paroi_scal_hyd_base_VEF,loiparth.valeur());

          const Domaine_VEF& domaine_VEF=ref_cast(Domaine_VEF,my_eqn.domaine_dis());
          const IntTab& elem_faces = domaine_VEF.elem_faces();
          const DoubleVect& vol = domaine_VEF.volumes();
          const IntTab& face_voisins = domaine_VEF.face_voisins();
          const DoubleTab& face_normale = domaine_VEF.face_normales();
          const Domaine_Cl_VEF& domaine_Cl_VEF=ref_cast(Domaine_Cl_VEF,my_eqn.domaine_Cl_dis());
          int nb_dim_pb=Objet_U::dimension;
          DoubleVect le_mauvais_gradient(nb_dim_pb);
          int nb_front=domaine_VEF.nb_front_Cl();
          for (int n_bord=0; n_bord<nb_front; n_bord++)
            {
              const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
              const Cond_lim_base& cl_base=la_cl.valeur();
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());

              // Wall laws only apply when the BC is of imposed-temperature type,
              // because in other cases (imposed flux, adiabatic) the wall flux is known and fixed.
              int ldp_appli=0;
              if ((sub_type(Neumann_paroi,cl_base))||(sub_type(Neumann_homogene,cl_base)))
                {
                  ldp_appli=1;
                  // if (paroi_scal_vef.get_flag_calcul_ldp_en_flux_impose() )
                  if (loiparth.get_flag_calcul_ldp_en_flux_impose() )
                    ldp_appli=0;
                }

              if (ldp_appli)
                {
                  // const DoubleVect& d_equiv = paroi_scal_vef.equivalent_distance(n_bord);
                  const DoubleVect& d_equiv = loiparth.equivalent_distance(n_bord);
                  // d_equiv contains the equivalent distance for the boundary.
                  // In d_equiv, for faces that are not fixed-wall (e.g. periodic, symmetry, etc.)
                  // the geometric distance is used, set during initialization in the wall law.

                  int num1 = 0;
                  int num2 = le_bord.nb_faces_tot();
                  ToDo_Kokkos("boundary");
                  for (int ind_face=num1; ind_face<num2; ind_face++)
                    {
                      // Tf is the average fluid temperature in the first element,
                      // without accounting for the wall temperature.
                      double Tf=0.;
                      // double bon_gradient=0.; // norm of the temperature gradient normal to the wall, computed using the wall law.
                      le_mauvais_gradient=0.;
                      int num_face = le_bord.num_face(ind_face);
                      int elem1 = face_voisins(num_face,0);
                      if (elem1==-1) elem1 = face_voisins(num_face,1);
                      double surface_face = domaine_VEF.face_surfaces(num_face);

                      int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
                      for (int i=0; i<nb_faces_elem; i++)
                        {
                          int j = elem_faces(elem1,i);
                          if ( j != num_face )
                            {
                              double surface_pond = 0.;
                              for (int kk=0; kk<nb_dim_pb; kk++)
                                surface_pond -= (face_normale(j,kk)*domaine_VEF.oriente_normale(j,elem1)*face_normale(num_face,kk)*
                                                 domaine_VEF.oriente_normale(num_face,elem1))/(surface_face*surface_face);
                              Tf+=inconnue(j)*surface_pond;
                            }

                          for(int kk=0; kk<nb_dim_pb; kk++)
                            le_mauvais_gradient(kk)+=inconnue(j)*face_normale(j,kk)*domaine_VEF.oriente_normale(j,elem1);
                        }
                      le_mauvais_gradient/=vol(elem1);
                      // mauvais_gradient = le_mauvais_gradient.n
                      double mauvais_gradient=0;
                      for(int kk=0; kk<nb_dim_pb; kk++)
                        mauvais_gradient+=le_mauvais_gradient(kk)*face_normale(num_face,kk)/surface_face;

                      // valeurs_espace(num_face) is the wall temperature: Tw.
                      // The sign of the good gradient does not matter since it is the norm of the
                      // temperature gradient in the element.
                      // It will then be multiplied by the normal vector to the wall face,
                      // which has the correct signs.
                      //bon_gradient=(Tf-inconnue(num_face))/d_equiv(ind_face)*(-domaine_VEF.oriente_normale(num_face,elem1));
                      // we must have bon=mauvais...
                      // compute so that bon_grad = mauvais
                      valeurs_espace(num_face)=Tf-mauvais_gradient*d_equiv(ind_face)*(-domaine_VEF.oriente_normale(num_face,elem1));
                    }
                }
            }
        }
    }
  DoubleTab& espace_valeurs = espace_stockage->valeurs();
  espace_valeurs.echange_espace_virtuel();

  return espace_stockage;
}

const Noms Champ_Generique_Tparoi_VEF::get_property(const Motcle& query) const
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

// Name the field as source by default
// "Combinaison_"+nom_champ_source
void Champ_Generique_Tparoi_VEF::nommer_source()
{
  if (nom_post_=="??")
    {
      Nom nom_post_source, nom_champ_source;
      const Noms nom = get_source(0).get_property("nom");
      nom_champ_source = nom[0];
      nom_post_source =  "Tparoi_";
      nom_post_source +=  nom_champ_source;
      nommer(nom_post_source);
    }
}
