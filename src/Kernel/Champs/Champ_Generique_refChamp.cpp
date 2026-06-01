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

#include <Champ_Generique_refChamp.h>
#include <Champ_Inc_P0_base.h>
#include <Champ_Inc_P1_base.h>
#include <Champ_Fonc_P0_base.h>
#include <Champ_Fonc_P1_base.h>
#include <Domaine_VF.h>
#include <Champ_Uniforme.h>
#include <Champ_Inc_Q1_base.h>
#include <Champ_Fonc_Q1_base.h>
#include <Equation_base.h>
#include <Synonyme_info.h>
#include <Param.h>
#include <Postraitement.h>

Implemente_instanciable(Champ_Generique_refChamp,"refChamp",Champ_Generique_base);
// XD refchamp champ_generique_base refchamp INHERITS_BRACE Field of prolem

Add_synonym(Champ_Generique_refChamp,"Champ_Post_refChamp");

Sortie& Champ_Generique_refChamp::printOn(Sortie& os) const
{
  return os;
}

Entree& Champ_Generique_refChamp::readOn(Entree& is)
{
  return Champ_Generique_base::readOn(is);
}

/*! @brief pb_champ :   triggers the reading of the problem name (nom_pb_) to which the discrete field belongs and the name of that discrete field (nom_champ_)
 *
 *   nom_source : option to name the field as a source (otherwise named by default)
 *
 */
void Champ_Generique_refChamp::set_param(Param& param) const
{
  param.ajouter_non_std("nom_source",(this)); // XD attr nom_source chaine nom_source OPT The alias name for the field
  param.ajouter_non_std("Pb_champ",(this),Param::REQUIRED); // XD attr pb_champ deuxmots pb_champ REQ { Pb_champ nom_pb
  // XD_CONT nom_champ } : nom_pb is the problem name and nom_champ is the selected field name.
}

int Champ_Generique_refChamp::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="Pb_champ")
    {
      // Read the problem name and the field name
      is >> nom_pb_ >> nom_champ_;
      // Search for the problem among the objects known to the interpreter
      const Objet_U& ob = interprete().objet(nom_pb_);
      if (!sub_type(Probleme_base, ob))
        {
          Cerr << "Error in Champ_Generique_refChamp ::lire(keyword=\"pb_champ\"\n"
               << " The object " << nom_pb_ << " is not a Probleme_base" << finl;
          exit();
        }
      Probleme_base& pb = ref_cast_non_const(Probleme_base, ob);
      ref_pb_ = pb;
      // Search for the field "nom_champ" in the problem:
      OBS_PTR(Champ_base) ref_champ;
      Noms liste_noms;
      pb.get_noms_champs_postraitables(liste_noms);
      pb.creer_champ(nom_champ_);
      ref_champ = pb.get_champ(nom_champ_);
      ref_champ->corriger_unite_nom_compo();
      set_ref_champ(ref_champ.valeur());
      return 1;
    }
  else if (mot=="nom_source")
    {
      Nom id_source;
      is >> id_source;
      nommer(id_source);
      return 1;
    }
  else
    return Champ_Generique_base::lire_motcle_non_standard(mot,is);
}

/*! @brief Initialize the class with the given field.
 *
 * We take a reference to this field (it must remain valid afterwards).
 *
 */
void Champ_Generique_refChamp::initialize(const Champ_base& champ)
{
  ref_champ_ = champ;
}

/*! @brief Returns the number of coordinates of each vertex of the domain.
 *
 * See GenericField_base::get_dimension()
 *
 */
int Champ_Generique_refChamp::get_dimension() const
{
  // We do not use Objet_U::dimension intentionally as we hope to remove
  // this static variable soon.
  const DoubleTab coords = get_ref_coordinates();
  const int dim = coords.dimension(1);
  return dim;
}

void Champ_Generique_refChamp::get_property_names(Motcles& list) const
{
  list.add("nom");
  list.add("nom_cible");
  list.add("unites");
  list.add("composantes");
  list.add("synonyms");
}

const Noms Champ_Generique_refChamp::get_property(const Motcle& query) const
{

  Motcles motcles(5);
  motcles[0] = "nom";
  motcles[1] = "nom_cible";
  motcles[2] = "unites";
  motcles[3] = "composantes";
  motcles[4] = "synonyms";

  int rang = motcles.search(query);
  switch(rang)
    {
    case 0:
      {
        Noms mots(1);
        mots[0] = nom_post_;
        return mots;
      }
    case 1:
      {
        Noms mots(1);
        mots[0] = nom_champ_;
        return mots;
      }
    case 2 :
      {
        const Noms mots0= get_ref_champ_base().unites();
        ref_cast_non_const(Champ_base,get_ref_champ_base()).corriger_unite_nom_compo();
        const Noms mots = get_ref_champ_base().unites();
        if (mots.size()!=mots0.size())
          {
            Cerr<<"iuuuuuu"<<mots<<" "<<mots0<<finl;
            exit();
          }
        for (int i=0; i<mots.size(); i++)
          {
            if (mots0[i]!=mots[i])
              {
                Cerr <<" iiiiiiiiii"<<mots0[i]<< " "<<mots[i]<<finl;
                exit();
              }
          }
        return mots;
      }
    case 4 :
      {
        if (syno_.size()>0)
          return syno_;

        const Noms mots = get_ref_champ_base().get_synonyms();

        return mots;
      }
    case 3 :
      {
        if (compo_.size()>1)
          return compo_;

        const Noms mots = get_ref_champ_base().noms_compo();
        int nb_comp = mots.size();

        Noms compo(nb_comp);
        for (int i=0; i<nb_comp; i++)
          {
            Nom nume(i);
            compo[i] = nom_post_+nume;
          }
        return compo;
      }
    default :
      {
        Cerr<<"The identifiable properties are : "<<motcles<<finl;
        exit();
      }
    }
  //For compilation
  // We never reach here
  return get_property(query);

}

/*! @brief If the field is not a discrete field: exception Champ_Generique_erreur("INVALID") Otherwise, returns the localisation of the field for support "index".
 *
 *   If index == -1 (default value), raises an exception if the field has multiple supports.
 *   Otherwise, raises an exception if the index exceeds the number of field localisations.
 *
 */
Entity Champ_Generique_refChamp::get_localisation(const int index) const
{
  Entity loc;
  //For initialization
  loc =Entity::NODE;

  const Champ_base& ch = get_ref_champ_base();
  const Domaine_dis_base& z_dis_base = get_ref_domaine_dis_base();

  // Discrete fields with a single localisation:
  if ((sub_type(Champ_Inc_P0_base, ch) || sub_type(Champ_Fonc_P0_base, ch)) && index <= 0)
    {
      loc = Entity::ELEMENT;
    }
  else if ((sub_type(Champ_Inc_P1_base, ch)|| sub_type(Champ_Fonc_P1_base, ch)) && index <= 0)
    {
      loc = Entity::NODE;
    }
  else if ((ch.que_suis_je().debute_par("Champ_Face_PolyMAC")
            || ch.que_suis_je().debute_par("Champ_Fonc_Face_PolyMAC")
            || ch.valeurs().dimension(0) == ref_cast(Domaine_VF,z_dis_base).nb_faces()) && index <= 0)
    {
      loc = Entity::FACE;
    }
  else
    {
      // Discrete fields with multiple localisations
      Nom message="Invalid localization used for postprocessing the field ";
      message+=ch.le_nom()+". Change your data file.";
      throw Champ_Generique_erreur(message);
    }
  return loc;
}

/*! @brief Verifies that the field is indeed a discrete field and returns the value array.
 *
 * Otherwise, raises the exception Champ_Generique_erreur("NO_REF")
 *
 */
const DoubleTab& Champ_Generique_refChamp::get_ref_values() const
{
  // Call to get_localisation to verify that the field is indeed a discrete field
  // (multi-support or not)
  get_localisation(0);
  // Returns the field values
  const DoubleTab& val = get_ref_champ_base().valeurs();
  return val;
}

/*! @brief Creates a copy of the value array. See GenericField_base::get_copy_values()
 *
 */
void Champ_Generique_refChamp::get_copy_values(DoubleTab& values) const
{
  const DoubleTab& val = get_ref_values();
  // Creates a copy of the array
  values = val;
}

/*! @brief call to Champ_base::valeur_aux()
 *
 */
void Champ_Generique_refChamp::get_xyz_values(const DoubleTab& coords, DoubleTab& values, ArrOfBit& validity_flag) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

const Domaine_Cl_dis_base& Champ_Generique_refChamp::get_ref_zcl_dis_base() const
{
  const Champ_base& ch = get_ref_champ_base();
  if (sub_type(Champ_Inc_base,ch))
    return ref_cast(Champ_Inc_base,ch).equation().domaine_Cl_dis();
  else
    {
      Cerr<<"No zcl_dis is available for the field "<<ch.que_suis_je()<<finl;
      exit();
    }

  //For compilation
  return get_ref_zcl_dis_base();
}

const DoubleTab& Champ_Generique_refChamp::get_ref_coordinates() const
{
  const DoubleTab& coord = get_ref_domain().coord_sommets();
  return coord;
}

void Champ_Generique_refChamp::get_copy_coordinates(DoubleTab& coordinates) const
{
  const DoubleTab& coord = get_ref_coordinates();
  coordinates = coord;
}

const IntTab& Champ_Generique_refChamp::get_ref_connectivity(Entity index1, Entity index2) const
{
  const Champ_base& ch = get_ref_champ_base();
  const Domaine_dis_base& domaine_dis_base = ch.domaine_dis_base();
  const Domaine& domaine = domaine_dis_base.domaine();
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis_base);

  switch(index1)
    {
    case Entity::ELEMENT:
      {
        switch(index2)
          {
          case Entity::NODE:
            return domaine.les_elems();
          case Entity::FACE:
            return domaine_vf.elem_faces();
          default :
            {
              exit();
            }
          }
        break;
      }
    case Entity::FACE:
      {
        switch(index2)
          {
          case Entity::NODE:
            return domaine_vf.face_sommets();
          case Entity::ELEMENT:
            return domaine_vf.face_voisins();
          default :
            {
              exit();
            }
          }
        break;
      }
    default :
      {
        exit();
      }
    }
  Nom message="Invalid localization used for postprocessing the field ";
  message+=ch.le_nom()+". Change your data file.";
  throw Champ_Generique_erreur(message);
}

void Champ_Generique_refChamp::get_copy_connectivity(Entity index1, Entity index2, IntTab& tab) const
{
  const IntTab& connectivity = get_ref_connectivity(index1, index2);
  tab = connectivity;
}

// Returns the problem that carries the target field
const Probleme_base& Champ_Generique_refChamp::get_ref_pb_base() const
{
  return ref_pb_.valeur();
}

/*! @brief Returns the underlying champ_base.
 *
 * Tests that the field has been properly associated. Eventually, this method will be removed from GenericField_base but remains
 *   here (to test if the field has been associated).
 *  Exceptions:
 *   Champ_Generique_erreur("NOT_INITIALIZED")
 *
 */
const Champ_base& Champ_Generique_refChamp::get_ref_champ_base() const
{
  if (!ref_champ_)
    throw Champ_Generique_erreur("NOT_INITIALIZED");
  return ref_champ_.valeur();
}

void Champ_Generique_refChamp::reset()
{
  ref_champ_.reset();
  localisation_ = Motcle();
}

void Champ_Generique_refChamp::completer(const Postraitement_base& post)
{
  nommer_source(post);
}

/*! @brief See Champ_Generique_base::mettre_a_jour If the field is a champ_inc in the equation, it must already have
 *
 *   been updated by the equation.
 *  If it is a calculated field, the update is done in get_champ
 *
 */
void Champ_Generique_refChamp::mettre_a_jour(double temps)
{

}

/*! @brief See Champ_Generique_base::get_champ.
 *
 * Here, the storage space is not used, the field already exists
 *
 */
const Champ_base& Champ_Generique_refChamp::get_champ(OWN_PTR(Champ_base)& espace_stockage) const
{
  const Nom& nom_cible = get_ref_champ_base().le_nom();
  ref_pb_->get_champ(nom_cible);
  return get_ref_champ_base();
}

const Champ_base& Champ_Generique_refChamp::get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const
{
  return get_champ(espace_stockage);
}

/*! @brief Associates the field and determines its localisation.
 *
 */
void Champ_Generique_refChamp::set_ref_champ(const Champ_base& champ)
{
  assert(!ref_champ_);
  ref_champ_ = champ;

  // Determination of the field localisation
  // We handle a few fields of known types...
  // For the others, they must be interpolated somewhere.
  // The ideal would be for the localisation to be a property of the discrete field itself.

  Nom type = champ.que_suis_je();
  type.majuscule();

  if (type.debute_par("CHAMP_P0"))
    {
      localisation_ = "ELEMENTS";
    }
  else if (type == "CHAMP_P1")
    {
      localisation_ = "SOMMETS";
    }
  else if (type.debute_par("CHAMP_FACE"))
    {
      localisation_ = "FACES";
    }
  else if (type == "CHAMP_P1NC")
    {
      localisation_ = "FACES";
    }
  else
    {
      localisation_ = "LOCALISATION_INCONNUE";
    }
}

//Returns the time of the target field
double Champ_Generique_refChamp::get_time() const
{
  double temps;
  temps = get_ref_champ_base().temps();
  return temps;
}

//Returns the directive (champ_elem, champ_sommets, champ_face or pression)
//to launch the discretisation of the storage space returned by
//the get_champ() method of the Champ_Generique_base that triggered the call
const Motcle Champ_Generique_refChamp::get_directive_pour_discr() const
{
  Motcle directive;
  const Champ_base& ch = get_ref_champ_base();

  // Discrete fields with a single localisation:
  if (sub_type(Champ_Inc_P0_base,ch) || sub_type(Champ_Fonc_P0_base,ch))
    directive = ch.is_basis_function() ? "champ_elem_DG" : (ch.is_quadrature() ? "champ_fonc_quad_DG" : "champ_elem");
  else if (sub_type(Champ_Inc_P1_base,ch) || sub_type(Champ_Fonc_P1_base,ch)
           || sub_type(Champ_Inc_P1_base,ch) || sub_type(Champ_Inc_Q1_base,ch)
           || sub_type(Champ_Fonc_Q1_base,ch))
    {
      directive = "champ_sommets";
      //   assert(localisation_=="SOMMETS");
    }
  else
    {
      const Nom& type = ch.que_suis_je();
      if ((type.debute_par("Champ_Face")) || (type=="Champ_P1NC") || (type=="Champ_Q1NC") ||
          (type=="Champ_Fonc_Face") || (type=="Champ_Fonc_P1NC") || (type=="Champ_Fonc_Q1NC"))
        {
          directive = "champ_face";
          //      assert(localisation_=="FACES");
        }
      else if ((type=="Champ_P1_isoP1Bulle") || (type=="Champ_Fonc_P1_isoP1Bulle"))
        {
          directive = "pression";
        }
      else if (sub_type(Champ_Uniforme,ch))
        {
          directive = "champ_uniforme";
        }

      else if (sub_type(Champ_Don_base,ch))
        {
          directive = "champ_don";
        }

      else
        {
          Cerr<<"No directive is available to create a storing field for the source field of type "<<ch.que_suis_je()<<finl;
          exit();
        }
    }

  return directive;

}

void Champ_Generique_refChamp::set_nom_champ(const Motcle& nom)
{
  nom_champ_=nom;
}

//Name the field as a source by default
//nom_champ_base + "_natif_" + nom_dom_natif
void Champ_Generique_refChamp::nommer_source(const Postraitement_base& post)
{
  if (nom_post_=="??")
    {
      Nom nom_post_source, nom_champ_base, nom_dom_natif;
      nom_champ_base = get_ref_champ_base().le_nom();
      if (ref_cast_non_const(Postraitement, post).domaine())
        {
          nom_post_source =  nom_champ_base + "_natif_" + ref_cast_non_const(Postraitement, post).domaine()->le_nom();
        }
      else
        {
          nom_dom_natif = get_ref_domain().le_nom();
          nom_post_source = nom_champ_base + "_natif_" + nom_dom_natif;
        }
      nommer(nom_post_source);
    }
}

int Champ_Generique_refChamp::get_info_type_post() const
{
  return 0;
}

const Noms& Champ_Generique_refChamp::fixer_noms_synonyms(const Noms& noms)
{
  return syno_ = noms;
}
const Noms& Champ_Generique_refChamp::fixer_noms_compo(const Noms& noms)
{
  return compo_ = noms;
}
