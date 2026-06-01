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

#include <Champ_Input_P0_Composite.h>
#include <Discretisation_tools.h>
#include <Discretisation_base.h>
#include <Champ_Fonc_Tabule.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Champ_Fonc_MED.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <EChaine.h>
#include <Device.h>
#include <Param.h>


Implemente_base(Milieu_base,"Milieu_base",Objet_U);
// XD milieu_base objet_u milieu_base INHERITS_BRACE Basic class for medium (physics properties of medium).

Sortie& Milieu_base::printOn(Sortie& os) const
{
  os << "{" << finl;
  ecrire(os);
  os << "}" << finl;
  return os;
}

/*! @brief Writes a medium object to an output stream.
 *
 * Writes the characteristics of the medium:
 *         - density
 *         - conductivity
 *         - heat capacity
 *         - beta_th
 *         - porosity
 *
 * @param (Sortie& os) the output stream for writing
 */
void  Milieu_base::ecrire(Sortie& os) const
{
  os << "rho " << ch_rho_<< finl;
  os << "lambda " << ch_lambda_ << finl;
  os << "Cp " << ch_Cp_ << finl;
  os << "beta_th " << ch_beta_th_ << finl;
  os << "porosite " << ch_porosites_ << finl;
}

/*! @brief Reading of a medium from an input stream.
 *
 * Format:
 *     {
 *      physical_quantity field_type field_reading_block
 *     }
 *  cf set_param method to know the understood keywords
 *  cf Param class to know possible options reading
 *
 * @param (Entree& is) the input stream for reading the medium parameters
 * @return (Entree&) the modified input stream
 * @throws opening brace expected
 */
Entree& Milieu_base::readOn(Entree& is)
{
  Cerr << "Reading of data for a " << que_suis_je() << " medium" << finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  check_gravity_vector();
  creer_champs_non_lus();
  if (ch_g_) champs_don_.add(ch_g_.valeur());
  if (ch_alpha_) champs_don_.add(ch_alpha_.valeur());
  if (ch_lambda_) champs_don_.add(ch_lambda_.valeur());
  if (ch_Cp_) champs_don_.add(ch_Cp_.valeur());
  if (ch_beta_th_) champs_don_.add(ch_beta_th_.valeur());
  if (ch_porosites_) champs_don_.add(ch_porosites_.valeur());
  if (ch_diametre_hyd_) champs_don_.add(ch_diametre_hyd_.valeur());
  return is;
}

void Milieu_base::set_param(Param& param) const
{
  param.ajouter("rho", &ch_rho_);   // XD_ADD_P field_base
  // XD_CONT Density (kg.m-3).
  param.ajouter("lambda", &ch_lambda_); // XD_ADD_P field_base
  // XD_CONT Conductivity (W.m-1.K-1).
  param.ajouter("Cp", &ch_Cp_);  // XD_ADD_P field_base
  // XD_CONT Specific heat (J.kg-1.K-1).
  param.ajouter("beta_th", &ch_beta_th_);
  set_additional_params(param);
}

// Useful method for F5! F5 does not call Milieu_base::set_param but Milieu_base::set_additional_params ...
void Milieu_base::set_additional_params(Param& param) const
{
  param.ajouter("diametre_hyd_champ", &ch_diametre_hyd_); // XD_ADD_P field_base
  // XD_CONT Hydraulic diameter field (optional).
  param.ajouter("porosites_champ", &ch_porosites_); // XD_ADD_P field_base
  // XD_CONT The porosity is given at each element and the porosity at each face, Psi(face), is calculated by the
  // XD_CONT average of the porosities of the two neighbour elements Psi(elem1), Psi(elem2) :
  // XD_CONT Psi(face)=2/(1/Psi(elem1)+1/Psi(elem2)). This keyword is optional.
  param.ajouter("porosites", &porosites_); // XD_ADD_P porosites
  // XD_CONT Porosities.
  // for F5, gravity is added here ...
  param.ajouter("gravite", &ch_g_); // XD_ADD_P field_base
  // XD_CONT Gravity field (optional).
}

int Milieu_base::lire_motcle_non_standard(const Motcle& mot_lu, Entree& is)
{
  return -1;
}

void Milieu_base::discretiser(const Probleme_base& pb, const  Discretisation_base& dis)
{
  Cerr << "Medium discretization." << finl;
  const Domaine_dis_base& domaine_dis=pb.domaine_dis();

  // PL: no time to make it cleaner, I do as in Fluide_Incompressible::discretiser to handle a conductivity read from a MED file. Test: Reprise_grossier_fin_VEF
  // ToDo: rewrite these two discretiser methods

  // E. Saikali : The thermal conductivity and diffusivity fields are considered as multi_scalaire fields, sure if the number of components read
  // in the data file for lambda > 1, i.e: case of anisotropic diffusion !

  int lambda_nb_comp = 0;

  if(ch_lambda_)
    {
      // Returns number of components of lambda field
      lambda_nb_comp = ch_lambda_->nb_comp( );
      if (sub_type(Champ_Fonc_MED,ch_lambda_.valeur()))
        {
          double temps=ch_lambda_->temps();
          Cerr<<"Convert Champ_fonc_MED lambda in Champ_Don..."<<finl;
          OWN_PTR(Champ_Don_base) ch_lambda_prov;
          dis.discretiser_champ("champ_elem",domaine_dis,"neant","neant",lambda_nb_comp,temps,ch_lambda_prov);
          ch_lambda_prov->affecter(ch_lambda_.valeur());
          ch_lambda_.detach();
          ch_alpha_.detach();
          dis.discretiser_champ("champ_elem",domaine_dis,"neant","neant",lambda_nb_comp,temps,ch_lambda_);
          ch_lambda_->valeurs()=ch_lambda_prov->valeurs();
        }

      if(lambda_nb_comp >1) // For anisotropic case
        ch_lambda_->fixer_nature_du_champ(multi_scalaire);

      dis.nommer_completer_champ_physique(domaine_dis,"conductivite","W/m/K",ch_lambda_.valeur(),pb);

      // the real name will be given later
      if (sub_type(Champ_Fonc_Tabule,ch_lambda_.valeur()))
        {
          double temps=ch_lambda_->temps();
          dis.discretiser_champ("champ_elem",domaine_dis,"neant","neant",lambda_nb_comp,temps,ch_alpha_);
          dis.discretiser_champ("champ_elem",domaine_dis,"neant","neant",lambda_nb_comp,temps,ch_alpha_fois_rho_);
        }
      champs_compris_.ajoute_champ(ch_lambda_.valeur());
    }
  if (!ch_alpha_&&(bool(ch_lambda_)))
    {
      double temps=ch_lambda_->temps();
      // ch_alpha (i.e. diffusivite_thermique) will have same component number as ch_lambda
      dis.discretiser_champ("champ_elem",domaine_dis,"neant","neant",lambda_nb_comp,temps,ch_alpha_);
      dis.discretiser_champ("champ_elem",domaine_dis,"neant","neant",lambda_nb_comp,temps,ch_alpha_fois_rho_);
    }
  if (ch_alpha_)
    {
      dis.nommer_completer_champ_physique(domaine_dis,"diffusivite_thermique","m2/s",ch_alpha_.valeur(),pb);
      dis.nommer_completer_champ_physique(domaine_dis,"alpha_fois_rho","kg/ms",ch_alpha_fois_rho_.valeur(),pb);
      champs_compris_.ajoute_champ(ch_alpha_.valeur());
      champs_compris_.ajoute_champ(ch_alpha_fois_rho_.valeur());
    }
  if (ch_beta_th_)
    {
      dis.nommer_completer_champ_physique(domaine_dis,"dilatabilite","K-1",ch_beta_th_.valeur(),pb);
      champs_compris_.ajoute_champ(ch_beta_th_.valeur());
    }
  if  (ch_Cp_)
    {
      dis.nommer_completer_champ_physique(domaine_dis,"capacite_calorifique","J/kg/K",ch_Cp_.valeur(),pb);
      champs_compris_.ajoute_champ(ch_Cp_.valeur());
    }
  if  (ch_rho_)
    {
      dis.nommer_completer_champ_physique(domaine_dis,"masse_volumique","kg/m^3",ch_rho_,pb);
      champs_compris_.ajoute_champ(ch_rho_);
    }
  if (ch_rho_ && ch_Cp_)
    {
      assert (ch_rho_->nb_comp() == ch_Cp_->nb_comp());
      if(!ch_rho_Cp_comme_T_)
        {
          const double temps = pb.schema_temps().temps_courant();
          // E. Saikali
          // XXX : I reset nb_comp = 1, otherwise it blocks in Solveur_Masse_base => tab_divide_any_shape
          // because we do not have line_size % line_size_vx == 0 (case nb_comp >1 for rho and cp)
          // TODO : FIXME : a generic case must be coded in DoubleVect::tab_divide_any_shape... good luck
          dis.discretiser_champ(dis.is_dg() ? "champ_elem" : "temperature", domaine_dis, "rho_cp_comme_T", "J/m^3/K", 1 /* rho->nb_comp() */, temps, ch_rho_Cp_comme_T_);
          dis.discretiser_champ( "champ_elem", domaine_dis,    "rho_cp_elem", "J/m^3/K", 1 /* rho->nb_comp() */, temps,    ch_rho_Cp_elem_);
        }
      champs_compris_.ajoute_champ(ch_rho_Cp_comme_T_.valeur());
      champs_compris_.ajoute_champ(ch_rho_Cp_elem_.valeur());
    }

  discretiser_porosite(pb,dis);
  discretiser_diametre_hydro(pb,dis);
}

// Useful method for F5! F5 does not call Milieu_base::discretiser but Milieu_base::discretiser_porosite ...
void Milieu_base::discretiser_porosite(const Probleme_base& pb, const Discretisation_base& dis)
{
  if (!zdb_) zdb_ = pb.domaine_dis();
  const double temps = pb.schema_temps().temps_courant();
  Nom fld_name = "porosite_volumique", fld_unit = "rien";

  // We build porosite_face_ with a parallel descriptor
  const MD_Vector& md = ref_cast(Domaine_VF, zdb_.valeur()).md_vector_faces();
  if (!porosite_face_.get_md_vector())
    {
      MD_Vector_tools::creer_tableau_distribue(md, porosite_face_, RESIZE_OPTIONS::NOCOPY_NOINIT);
      assert (ref_cast(Domaine_VF, zdb_.valeur()).nb_faces_tot() == porosite_face_.size_totale());
    }
  porosite_face_ = 1.;

  if (ch_porosites_) // Read via porosites_champ
    {
      assert (!is_user_porosites());
      if (porosites_.is_read())
        {
          Cerr << "WHAT ?? You can not define in your medium both porosites_champ & porosites ! Remove one of them !" << finl;
          Process::exit();
        }

      is_field_porosites_ = true;

      if (sub_type(Champ_input_P0, ch_porosites_.valeur()))
        {
          Cerr << "To control the porosity field from ICoCo, please use Champ_Input_P0_Composite and not Champ_input_P0 !" << finl;
          Process::exit();
        }
      else if (sub_type(Champ_Input_P0_Composite, ch_porosites_.valeur()))
        {
          Champ_Input_P0_Composite& ch_in = ref_cast(Champ_Input_P0_Composite, ch_porosites_.valeur());
          if (!ch_in.is_initialized())
            {
              Cerr << "To control the porosity field from ICoCo, please define the initial field in your Champ_Input_P0_Composite !" << finl;
              Process::exit();
            }

          ch_porosites_->fixer_unite(fld_unit);
          ch_porosites_->valeurs() = ch_in.initial_values(); // Initialize!
        }
      else if (sub_type(Champ_Fonc_MED,ch_porosites_.valeur()))
        {
          Cerr<<"Convert Champ_fonc_MED " << fld_name << " to a OWN_PTR(Champ_Don_base) ..."<<finl;
          OWN_PTR(Champ_Don_base) tmp_fld;
          dis.discretiser_champ("champ_elem",zdb_.valeur(),"neant",fld_unit,1,temps,tmp_fld);
          tmp_fld->affecter(ch_porosites_.valeur()); // interpolate ...
          ch_porosites_.detach();
          dis.discretiser_champ("champ_elem",zdb_.valeur(),fld_name,fld_unit,1,temps,ch_porosites_);
          ch_porosites_->valeurs() = tmp_fld->valeurs();
        }
      else if (sub_type(Champ_Uniforme,ch_porosites_.valeur())) // blabla ...
        {
          const double val = ch_porosites_->valeurs()(0,0);
          ch_porosites_.detach();
          dis.discretiser_champ("champ_elem",zdb_.valeur(),fld_name,fld_unit,1,temps,ch_porosites_);
          ch_porosites_->valeurs() = val;
        }
      else
        dis.nommer_completer_champ_physique(zdb_.valeur(), fld_name, fld_unit, ch_porosites_.valeur(), pb);
    }
  else if (porosites_.is_read()) // via porosites
    {
      assert (!is_field_porosites());
      if (ch_porosites_)
        {
          Cerr << "WHAT ?? You can not define in your medium both porosites_champ & porosites ! Remove one of them !" << finl;
          Process::exit();
        }

      is_user_porosites_ = true;
      // We will now use porosites_champ!
      dis.discretiser_champ("champ_elem", zdb_.valeur(), fld_name, fld_unit, 1, temps, ch_porosites_);
      Domaine_VF& zvf = ref_cast_non_const(Domaine_VF, zdb_.valeur());
      ch_porosites_->valeurs() = 1.; // Initialize to 1 ...
      porosites_.remplir_champ(zvf, ch_porosites_->valeurs(), porosite_face_);
    }
  else // Not defined by the user
    {
      // We will now use porosites_champ!
      dis.discretiser_champ("champ_elem", zdb_.valeur(), fld_name, fld_unit, 1, temps, ch_porosites_);
      ch_porosites_->valeurs() = 1.; // Initialize to 1 ...
    }

  // We add for all cases
  if (sub_type(Champ_Input_P0_Composite, ch_porosites_.valeur()))
    {
      Champ_input_P0& ch_in = ref_cast(Champ_input_P0,ref_cast(Champ_Input_P0_Composite, ch_porosites_.valeur()).input_field());
      ch_in.add_synonymous(fld_name); // So that it can be known also ;-)
      champs_compris_.ajoute_champ(ch_in);
    }
  else champs_compris_.ajoute_champ(ch_porosites_.valeur());

  verifie_champ_porosites();
  if (is_field_porosites()) calculate_face_porosity(); /* otherwise it is already filled ... */
  fill_section_passage_face();
}

void Milieu_base::discretiser_diametre_hydro(const Probleme_base& pb, const Discretisation_base& dis)
{
  if (!zdb_) zdb_ = pb.domaine_dis();
  const double temps = pb.schema_temps().temps_courant();
  Nom fld_name = "diametre_hydraulique", fld_unit = "m";

  // We build porosite_face_ with a parallel descriptor
  const MD_Vector& md = ref_cast(Domaine_VF, zdb_.valeur()).md_vector_faces();
  if (!diametre_hydraulique_face_.get_md_vector())
    {
      MD_Vector_tools::creer_tableau_distribue(md, diametre_hydraulique_face_, RESIZE_OPTIONS::NOCOPY_NOINIT);
      assert (ref_cast(Domaine_VF, zdb_.valeur()).nb_faces_tot() == diametre_hydraulique_face_.size_totale());
    }
  diametre_hydraulique_face_ = 0.; /* hydraulic diameters are initialized to 0 */

  if (ch_diametre_hyd_) // Read via diametre_hyd_champ
    {
      has_hydr_diam_ = true;
      if (sub_type(Champ_Fonc_MED, ch_diametre_hyd_.valeur()))
        {
          Cerr << "Convert Champ_fonc_MED " << fld_name << " to a OWN_PTR(Champ_Don_base) ..." << finl;
          OWN_PTR(Champ_Don_base) tmp_fld;
          dis.discretiser_champ("champ_elem", zdb_.valeur(), "neant", fld_unit, 1, temps, tmp_fld);
          tmp_fld->affecter(ch_diametre_hyd_.valeur()); // interpolate ...
          ch_diametre_hyd_.detach();
          dis.discretiser_champ("champ_elem", zdb_.valeur(), fld_name, fld_unit, 1, temps, ch_diametre_hyd_);
          ch_diametre_hyd_->valeurs() = tmp_fld->valeurs();
        }
      else if (sub_type(Champ_Uniforme, ch_diametre_hyd_.valeur())) // blabla ...
        {
          const double val = ch_diametre_hyd_->valeurs()(0, 0);
          ch_diametre_hyd_.detach();
          dis.discretiser_champ("champ_elem", zdb_.valeur(), fld_name, fld_unit, 1, temps, ch_diametre_hyd_);
          ch_diametre_hyd_->valeurs() = val;
        }
      else
        dis.nommer_completer_champ_physique(zdb_.valeur(), fld_name, fld_unit, ch_diametre_hyd_.valeur(), pb);
    }
  else // Not defined by the user
    {
      dis.discretiser_champ("champ_elem", zdb_.valeur(), fld_name, fld_unit, 1, temps, ch_diametre_hyd_);
      ch_diametre_hyd_->valeurs() = 0.; // Initialize to 0 ...
    }

  champs_compris_.ajoute_champ(ch_diametre_hyd_.valeur());

  if (has_hydr_diam_) calculate_face_hydr_diam(); /* otherwise it is already filled ... */
}

void Milieu_base::calculate_face_hydr_diam()
{
  assert(has_hydr_diam_);
  const Domaine_VF& zvf = ref_cast(Domaine_VF, zdb_.valeur());
  const IntTab& f_e = zvf.face_voisins();
  const int nb_face_tot = zvf.nb_faces_tot();
  assert(diametre_hydraulique_face_.size_totale() == nb_face_tot);
  int e;
  for (int f = 0; f < nb_face_tot; f++)
    {
      double nv = 0.0;
      for (int i = 0; i < 2; i++)
        if ((e = f_e(f, i)) > -1)
          {
            diametre_hydraulique_face_(f) += ch_diametre_hyd_->valeurs()(e);
            nv += 1.0;
          }
      diametre_hydraulique_face_(f) /= nv;
    }
  // diametre_hydraulique_face_.echange_espace_virtuel(); // Elie: to be checked if useful... not used for now
}

void Milieu_base::verifie_champ_porosites()
{
  // XXX : Elie Saikali : Previously we tested this:
  // assert(mp_min_vect(porosites_champ->valeurs()) >= 0. && mp_max_vect(porosites_champ->valeurs()) <= 1.);
  // encountered an F5 case with printf("%.9g\n", 1.0 - mp_max_vect) = -2.88657986e-15
  // tried comparing with std::numeric_limits<double>::epsilon() but the overflow is larger !!
  // so I clean the field like this for now ... if it is reasonable!
  const double min_por = mp_min_vect(ch_porosites_->valeurs()), max_por = mp_max_vect(ch_porosites_->valeurs());

  if (min_por >= 0.0 && max_por <= 1.0) { /* do nothing */ }
  else if (min_por >= -1.e-12 && max_por <= 1. + 1.e-12 ) nettoie_champ_porosites();
  else
    {
      Cerr << "WHAT ?? Your porosity field contains values < 0 or > 1 !!!! You should do something !" << finl;
      Process::exit();
    }
}

void Milieu_base::nettoie_champ_porosites()
{
  Cerr << "************************************************************************" << finl;
  Cerr << " We detected some element-porosity values which are slighly < 0 or > 1 !" << finl;
  Cerr << finl;
  printf("Overflow : 1 - MAX = %.9g\n", 1.0 - mp_max_vect(ch_porosites_->valeurs()));
  printf("Overflow : MIN - 0 = %.9g\n", mp_min_vect(ch_porosites_->valeurs()) - 0.);
  Cerr << finl;
  Cerr << " We will clean the field to prevent overflows and numerical issues ... " << finl;
  Cerr << "************************************************************************" << finl;

  for (int i = 0; i < ch_porosites_->valeurs().dimension_tot(0); i++)
    {
      if (ch_porosites_->valeurs()(i) < 0.) ch_porosites_->valeurs()(i) = 0.;
      if (ch_porosites_->valeurs()(i) > 1.) ch_porosites_->valeurs()(i) = 1.;
    }

  assert(mp_min_vect(ch_porosites_->valeurs()) >= 0.0 && mp_max_vect(ch_porosites_->valeurs()) <= 1.0);
}

void Milieu_base::verifier_coherence_champs(int& err,Nom& msg)
{
  if (err==1)
    {
      Cerr<<"Error while reading the physical properties of a "<<que_suis_je()<<" medium."<<finl;
      Cerr<<msg<<finl;
      Process::exit();
    }
  else if (err==2)
    {
      Cerr<<"Warning while reading the physical properties of a "<<que_suis_je()<<" medium."<<finl;
      Cerr<<msg<<finl;
    }
  else
    Cerr<<"The physical properties of a "<<que_suis_je()<<" medium have been successfully checked."<<finl;
}

void Milieu_base::preparer_calcul()
{
  // does nothing!!!
}

void Milieu_base::check_gravity_vector() const
{
  if (ch_g_)
    if(Objet_U::dimension != ch_g_->nb_comp())
      {
        Cerr << "The dimension is " << Objet_U::dimension << " and you create a gravity vector with " << ch_g_->nb_comp() << " components." << finl;
        Process::exit();
      }
}

void Milieu_base::creer_champs_non_lus()
{
  if (ch_rho_ && ch_lambda_ && ch_Cp_)
    creer_alpha();
}

void Milieu_base::warn_old_syntax()
{
  if (que_suis_je() != "Fluide_Diphasique") // for FT for now ...
    {
      Cerr << "YOU ARE USING AN OLD SYNTAX IN YOUR DATA FILE AND THIS IS NO MORE SUPPORTED !" << finl;
      Cerr << "STARTING FROM TRUST-v1.9.3 : GRAVITY SHOULD BE READ INSIDE THE MEDIUM AND NOT VIA ASSOSCIATION ... " << finl;
      Cerr << "HAVE A LOOK TO ANY TRUST TEST CASE TO SEE HOW IT SHOULD BE DONE ($TRUST_ROOT/tests/) ... " << finl;
      Cerr << "OR RUN -convert_data OPTION OF YOUR APPLICATION SCRIPT, FOR TRUST FOR EXAMPLE:" << finl;
      Cerr << "   trust -convert_data " << Objet_U::nom_du_cas() << ".data" << finl;
      Process::exit();
    }
}

/*! @brief Associates gravity by dynamically checking the type of the object to associate.
 *
 *     If the object is of type OWN_PTR(Champ_Don_base) or Champ_Don_base
 *     the association succeeds, otherwise it fails.
 *
 * @param (Objet_U& ob) a TRUST object representing a gravity field
 * @return (int) 1 if the association succeeded, 0 otherwise.
 */
int Milieu_base::associer_(Objet_U& ob)
{
  if (sub_type(Champ_Don_base,ob))
    {
      warn_old_syntax();
      via_associer_ = true;
      associer_gravite(ref_cast(Champ_Don_base, ob));
      return 1;
    }
  if (dynamic_cast <Champ_Don_base*>(&ob) != nullptr)
    {
      warn_old_syntax();
      via_associer_ = true;
      associer_gravite(dynamic_cast<Champ_Don_base&>(ob));
      return 1;
    }
  return 0;
}

/*! @brief Associates (assigns) a gravity field to the medium.
 *
 * @param (Champ_Don_base& gravite) given field representing gravity
 */
void Milieu_base::associer_gravite(const Champ_Don_base& la_gravite)
{
  // We verify that gravity has the correct dimension
  if (Objet_U::dimension!=la_gravite.nb_comp())
    {
      Cerr << "The dimension is " << Objet_U::dimension << " and you create a gravity vector with " << la_gravite.nb_comp() << " components." << finl;
      exit();
    }
  g_via_associer_ = la_gravite;

  if (via_associer_ && ch_g_)
    {
      assert (g_via_associer_);
      Cerr << "WHAT ?? Remove the associer gravity line from your jdd because it is already in the medium !!!" << finl;
      Process::exit();
    }
}

/*! @brief Computes alpha=lambda/(rho*Cp).
 *
 * Depending on lambda, rho and Cp, alpha can be typed as
 *     a uniform field or a piecewise uniform field.
 *
 * @throws if a precondition is violated
 * @throws if alpha cannot be computed because the types of the fields
 * lambda, rho, Cp are incompatible or not handled.
 */
void Milieu_base::calculer_alpha()
{
  if(ch_lambda_)
    {
      DoubleTab& tabalpha = ch_alpha_->valeurs(), &tab_lambda_sur_cp = ch_alpha_fois_rho_->valeurs();

      tabalpha = ch_lambda_->valeurs();
      tab_lambda_sur_cp = ch_lambda_->valeurs();

      // [ABN]: allows variable rho, Cp at this level (will be used by Solide_Milieu_Variable for instance).
      if (sub_type(Champ_Uniforme,ch_rho_.valeur()))
        tabalpha /= ch_rho_->valeurs()(0,0);
      else
        tab_divide_any_shape(tabalpha,ch_rho_->valeurs());

      if (sub_type(Champ_Uniforme,ch_Cp_.valeur()))
        {
          tabalpha /= ch_Cp_->valeurs()(0,0);
          tab_lambda_sur_cp /= ch_Cp_->valeurs()(0,0);
        }
      else
        {
          tab_divide_any_shape(tabalpha,ch_Cp_->valeurs());
          tab_divide_any_shape(tab_lambda_sur_cp,ch_Cp_->valeurs());
        }
    }
  else
    Cerr << "alpha calculation is not possible since lambda is not known." << finl;
}

void Milieu_base::calculate_face_porosity()
{
  assert(is_field_porosites());
  const Domaine_VF& zvf = ref_cast(Domaine_VF, zdb_.valeur());

  // update porosite_faces
  const int nb_face_tot = zvf.nb_faces_tot();
  assert (nb_face_tot == porosite_face_.size_totale());
  const IntTab& face_voisins = zvf.face_voisins();

  for (int face = 0; face < nb_face_tot; face++)
    {
      const int elem1 = face_voisins(face, 1), elem2 = face_voisins(face, 0);
      if ((elem1 != -1) && (elem2 != -1))
        porosite_face_(face) = 2. / (1. / ch_porosites_->valeurs()(elem1) + 1. / ch_porosites_->valeurs()(elem2));
      else if (elem1 != -1)
        porosite_face_(face) = ch_porosites_->valeurs()(elem1);
      else
        porosite_face_(face) = ch_porosites_->valeurs()(elem2);
    }
  porosite_face_.echange_espace_virtuel();
}

void Milieu_base::mettre_a_jour(double temps)
{
  //Cerr << que_suis_je() << "Milieu_base::mettre_a_jour" << finl;
  if (ch_rho_) ch_rho_->mettre_a_jour(temps);

  if (ch_g_) ch_g_->mettre_a_jour(temps);

  if (g_via_associer_) g_via_associer_->mettre_a_jour(temps);

  if (ch_lambda_) ch_lambda_->mettre_a_jour(temps);

  if (ch_Cp_) ch_Cp_->mettre_a_jour(temps);

  if (ch_beta_th_) ch_beta_th_->mettre_a_jour(temps);

  if ( (bool(ch_lambda_)) && (bool(ch_Cp_)) && (bool(ch_rho_)) )
    {
      calculer_alpha();
      ch_alpha_->changer_temps(temps);
      ch_alpha_fois_rho_->changer_temps(temps);
    }

  if (ch_rho_Cp_comme_T_) update_rho_cp(temps);

  mettre_a_jour_porosite(temps); // for F5!
}

// Utility method for F5! F5 does not call Milieu_base::mettre_a_jour but Milieu_base::mettre_a_jour_porosite ...
void Milieu_base::mettre_a_jour_porosite(double temps)
{
  assert(ch_porosites_ && ch_diametre_hyd_);
  if (is_field_porosites())
    if (sub_type(Champ_Input_P0_Composite, ch_porosites_.valeur()))
      {
        Cerr << "Updating porosity values since Champ_input_P0 !! We update also the face_porosity & section_passage fields ..." << finl;
        verifie_champ_porosites();
        calculate_face_porosity();
        fill_section_passage_face();
        ch_porosites_->changer_temps(temps);
        /* not needed I think but leaving as comment ;-) */
        // porosites_champ->valeurs().echange_espace_virtuel();
      }

  ch_porosites_->mettre_a_jour(temps); /* ne fait rien si Champ_Input_P0_Composite */
  ch_diametre_hyd_->mettre_a_jour(temps);
}

void Milieu_base::update_rho_cp(double temps)
{
  // If the unknown is on the device, copy the data there too:
  if (equation_.size() && (*(equation_.begin()->second)).inconnue().valeurs().isDataOnDevice())
    {
      mapToDevice(ch_rho_Cp_elem_->valeurs());
      mapToDevice(ch_rho_Cp_comme_T_->valeurs());
    }
  ch_rho_Cp_elem_->changer_temps(temps);
  ch_rho_Cp_elem_->changer_temps(temps);
  DoubleTab& rho_cp=ch_rho_Cp_elem_->valeurs();
  if (sub_type(Champ_Uniforme,ch_rho_.valeur()))
    rho_cp=ch_rho_->valeurs()(0,0);
  else
    {
      // AB: rho_cp = rho->valeurs() turns rho_cp into a 2 dimensional array with 1 compo. We want to stay mono-dim:
      rho_cp = 1.;
      tab_multiply_any_shape(rho_cp, ch_rho_->valeurs());
    }
  if (sub_type(Champ_Uniforme,ch_Cp_.valeur()))
    rho_cp*=ch_Cp_->valeurs()(0,0);
  else
    tab_multiply_any_shape(rho_cp,ch_Cp_->valeurs());
  ch_rho_Cp_comme_T_->changer_temps(temps);
  ch_rho_Cp_comme_T_->changer_temps(temps);
  const MD_Vector& md_som = ch_rho_Cp_elem_->domaine_dis_base().domaine().md_vector_sommets(),
                   &md_faces = ref_cast(Domaine_VF,ch_rho_Cp_elem_->domaine_dis_base()).md_vector_faces();
  if (ch_rho_Cp_comme_T_->valeurs().get_md_vector() == ch_rho_Cp_elem_->valeurs().get_md_vector())
    ch_rho_Cp_comme_T_->valeurs() = rho_cp;
  else if (ch_rho_Cp_comme_T_->valeurs().get_md_vector() == md_som)
    Discretisation_tools::cells_to_nodes(ch_rho_Cp_elem_,ch_rho_Cp_comme_T_);
  else if (ch_rho_Cp_comme_T_->valeurs().get_md_vector() == md_faces)
    Discretisation_tools::cells_to_faces(ch_rho_Cp_elem_,ch_rho_Cp_comme_T_);
  else
    {
      Cerr<< que_suis_je()<<(int)__LINE__<<finl;
      Process::exit();
    }
}

void Milieu_base::abortTimeStep()
{
  if (ch_rho_) ch_rho_->abortTimeStep();
}

void Milieu_base::resetTime(double time)
{
  if (ch_rho_)
    {
      if (sub_type(Champ_Don_base, ch_rho_.valeur()))
        ch_rho_->mettre_a_jour(time);
      else if (sub_type(Champ_Inc_base, ch_rho_.valeur()))
        ch_rho_->resetTime(time);
      else
        throw;
    }
}

void Milieu_base::creer_alpha()
{
  Cerr << "Milieu_base::creer_alpha (champ non lu)" << finl;
  assert(ch_lambda_);
  assert(ch_rho_);
  assert(ch_Cp_);
  ch_alpha_ = ch_lambda_;
  ch_alpha_fois_rho_ = ch_lambda_;
  ch_alpha_->nommer("alpha");
  ch_alpha_fois_rho_->nommer("alpha_fois_rho");
}

/*! @brief Returns the gravity of the medium if it has been associated, raises an error otherwise (const version).
 *
 * @return (Champ_Don_base&) the field representing the gravity of the medium
 * @throws if no gravity is associated with the medium
 */
const Champ_Don_base& Milieu_base::gravite() const
{
  if (!ch_g_ && !g_via_associer_)
    {
      Cerr << "The gravity has not been associated with the medium" << finl;
      Process::exit();
    }

  return ch_g_ ? ch_g_.valeur() : g_via_associer_.valeur();
}

/*! @brief Returns the gravity of the medium if it has been associated, raises an error otherwise.
 *
 * @return (Champ_Don_base&) the field representing the gravity of the medium
 */
Champ_Don_base& Milieu_base::gravite()
{
  if (!ch_g_ && !g_via_associer_)
    {
      Cerr << "The gravity has not been associated with the medium" << finl;
      Process::exit();
    }

  return ch_g_ ? ch_g_.valeur() : g_via_associer_.valeur();
}

int Milieu_base::initialiser(const double temps)
{
  Cerr << que_suis_je() << " Milieu_base:::initialiser" << finl;
  if (sub_type(Champ_Don_base, ch_rho_.valeur())) ref_cast(Champ_Don_base, ch_rho_.valeur()).initialiser(temps);

  if (ch_g_) ch_g_->initialiser(temps);

  if (g_via_associer_) g_via_associer_->initialiser(temps);

  if (ch_lambda_) ch_lambda_->initialiser(temps);

  if (ch_Cp_) ch_Cp_->initialiser(temps);

  if (ch_beta_th_) ch_beta_th_->initialiser(temps);

  if ( (bool(ch_lambda_)) && (bool(ch_Cp_)) && (bool(ch_rho_)) )
    {
      calculer_alpha();
      ch_alpha_->changer_temps(temps);
      ch_alpha_fois_rho_->changer_temps(temps);
    }

  if (ch_rho_Cp_comme_T_) update_rho_cp(temps);

  int err=0;
  Nom msg;
  verifier_coherence_champs(err,msg);

  return initialiser_porosite(temps);
}

// method useful for F5! F5 does not call Milieu_base::initialiser but Milieu_base::initialiser_porosite ...
int Milieu_base::initialiser_porosite(const double temps)
{
  // TODO : XXX : to check for ICoCo ? it should be initialized in main ?
  assert(ch_porosites_ && ch_diametre_hyd_);
  ch_porosites_->initialiser(temps);
  ch_diametre_hyd_->initialiser(temps);
  return 1;
}

void Milieu_base::fill_section_passage_face()
{
  const Domaine_VF& zvf = ref_cast(Domaine_VF, zdb_.valeur());
  const DoubleVect& fs = zvf.face_surfaces();
  section_passage_face_.resize(fs.size_array());
  for (int i = 0; i < fs.size_array(); i++) section_passage_face_[i] = fs[i] * porosite_face_[i];
}

/*! @brief Returns the mass density of the medium (const version).
 *
 * @return (Champ_base&) the field representing the mass density
 */
const Champ_base& Milieu_base::masse_volumique() const
{
  return ch_rho_.valeur();
}

/*! @brief Returns the mass density of the medium.
 *
 * @return (Champ_base&) the field representing the mass density
 */
Champ_base& Milieu_base::masse_volumique()
{
  return ch_rho_.valeur();
}

/*! @brief Returns the diffusivity of the medium (const version).
 *
 * @return (Champ_Don_base&) the field representing the diffusivity
 */
const Champ_Don_base& Milieu_base::diffusivite() const
{
  return ch_alpha_.valeur();
}

/*! @brief Returns the diffusivity of the medium.
 *
 * @return (Champ_Don_base&) the field representing the diffusivity
 */
Champ_Don_base& Milieu_base::diffusivite()
{
  return ch_alpha_.valeur();
}

const Champ_Don_base& Milieu_base::diffusivite_fois_rho() const
{
  return ch_alpha_fois_rho_.valeur();
}

Champ_Don_base& Milieu_base::diffusivite_fois_rho()
{
  return ch_alpha_fois_rho_.valeur();
}

/*! @brief Returns the conductivity of the medium (const version).
 *
 * @return (Champ_Don_base&) the field representing the conductivity
 */
const Champ_Don_base& Milieu_base::conductivite() const
{
  return ch_lambda_.valeur();
}

/*! @brief Returns the conductivity of the medium.
 *
 * @return (Champ_Don_base&) the field representing the conductivity
 */
Champ_Don_base& Milieu_base::conductivite()
{
  return ch_lambda_.valeur();
}

/*! @brief Returns the heat capacity of the medium (const version).
 *
 * @return (Champ_Don_base&) the field representing the heat capacity
 */
const Champ_Don_base& Milieu_base::capacite_calorifique() const
{
  return ch_Cp_.valeur();
}

/*! @brief Returns the heat capacity of the medium.
 *
 * @return (Champ_Don_base&) the field representing the heat capacity
 */
Champ_Don_base& Milieu_base::capacite_calorifique()
{
  return ch_Cp_.valeur();
}

/*! @brief Returns beta_t of the medium (const version).
 *
 * @return (Champ_Don_base&) the field representing beta_t
 */
const Champ_Don_base& Milieu_base::beta_t() const
{
  return ch_beta_th_.valeur();
}

/*! @brief Returns beta_t of the medium.
 *
 * @return (Champ_Don_base&) the field representing beta_t
 */
Champ_Don_base& Milieu_base::beta_t()
{
  return ch_beta_th_.valeur();
}

/*! @brief Returns 1 if gravity has been initialized.
 *
 * @return (int) 1 if g is non-null
 */
int Milieu_base::a_gravite() const
{
  return (ch_g_ || g_via_associer_) ? 1 : 0;
}

const Champ_base& Milieu_base::get_champ(const Motcle& nom) const
{
  return champs_compris_.get_champ(nom);
}

bool Milieu_base::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  return champs_compris_.has_champ(nom, ref_champ);
}

bool Milieu_base::has_champ(const Motcle& nom) const
{
  return champs_compris_.has_champ(nom);
}

void Milieu_base::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  if (opt == DESCRIPTION)
    Cerr << que_suis_je() << " : " << champs_compris_.liste_noms_compris() << finl;
  else
    nom.add(champs_compris_.liste_noms_compris());
}

/*! @brief Returns 0 if the medium is already associated with a problem, 1 otherwise.
 *
 * @return (int)
 */
int Milieu_base::est_deja_associe()
{
  if (deja_associe_==1)
    return 0;
  deja_associe_=1;
  return 1;
}

void Milieu_base::associer_equation(const Equation_base *eqn) const
{
  std::string nom_inco(eqn->inconnue().le_nom().getString());
  // E. Saikali
  // At the initialization step, FT problem can have several equations with same unknown name "concentration"
  if (equation_.count(nom_inco) && eqn->probleme().que_suis_je() != "Probleme_FT_Disc_gen")
    {
      Cerr << que_suis_je() << " multiple equations solve the unknown " << eqn->inconnue().le_nom() << " !" << finl;
      Process::exit();
    }
  equation_[nom_inco] = eqn;
}

const Equation_base& Milieu_base::equation(const std::string& inco) const
{
  return *equation_.at(inco);
}

void Milieu_base::set_id_composite(const int i)
{
  id_composite_ = i;
}

void Milieu_base::nommer(const Nom& nom)
{
  nom_ = nom;
}
const Nom& Milieu_base::le_nom() const
{
  return nom_;
}

