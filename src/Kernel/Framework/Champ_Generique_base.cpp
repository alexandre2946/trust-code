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

#include <Champ_Generique_base.h>
#include <Champs_compris.h>
#include <Probleme_base.h>
#include <Interprete.h>
#include <TRUST_Ref.h>
#include <Param.h>

Implemente_base(Champ_Generique_base,"Champ_Generique_base",Objet_U);
// XD champ_generique_base objet_u champ_generique_base BRACE not_set
// XD listchamp_generique listobj nul BRACE champ_generique_base COMMA XXX

// XD definition_champ objet_lecture nul NO_BRACE Keyword to create new complex field for advanced postprocessing.
// XD attr name chaine name REQ The name of the new created field.
// XD attr champ_generique champ_generique_base champ_generique REQ not_set

// XD definition_champs listobj nul BRACE definition_champ NO_COMMA List of definition champ

Sortie& Champ_Generique_base::printOn(Sortie& os) const
{
  return os;
}

Entree& Champ_Generique_base::readOn(Entree& is)
{
  Journal()<<"Reading data for a "<<que_suis_je()<<" field." <<finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is;
}

int Champ_Generique_base::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  return -1;
}

/*! @brief Returns the dimension of the space in which the field is defined.
 *
 * This is the number of components of the coordinates that must be provided in get_xyz_values.
 *   (for example, a field defined on a surface can be of dimension 3 if the
 *    vertex coordinates are 3D, or dimension 2 if it is curvilinear coordinates
 *    along the surface)
 *
 */
int Champ_Generique_base::get_dimension() const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
  // We never reach here
}

/*! @brief Returns the time of the Champ_Generique_base.
 *
 */
double Champ_Generique_base::get_time() const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

/*! @brief Returns the problem that carries the target field
 *
 */
const Probleme_base& Champ_Generique_base::get_ref_pb_base() const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

/*! @brief Returns the discretization associated with the problem
 *
 */
const Discretisation_base& Champ_Generique_base::get_discretisation() const
{
  const Objet_U& ob = interprete().objet(nom_pb_);
  const Probleme_base& pb = ref_cast(Probleme_base,ob);
  const Discretisation_base& discr = pb.discretisation();
  return discr;
}

/*! @brief Returns the directive (champ_elem, champ_sommets, champ_face or pression)
 * to launch the discretization of the storage space returned by get_champ()
 *
 */
const Motcle Champ_Generique_base::get_directive_pour_discr() const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

void Champ_Generique_base::nommer(const Nom& nom)
{
  nom_post_ = nom;
}

const Nom& Champ_Generique_base::get_nom_post() const
{
  return nom_post_;
}

/*! @brief Returns the list of possible "queries" for the field.
 *
 */
void Champ_Generique_base::get_property_names(Motcles& list) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

/*! @brief Returns the requested property.
 *
 * Examples: "DISCRETISATION" : type of discrete field (P0, P1, P1NC, etc...)
 *   "ELEMENT_TYPE" : type of element with the largest dimension ({TRIANGLE}, {TETRAHEDRA}, {QUAD}, {HEXA}, etc)
 *   "DYNAMIC_MESH" : is the mesh dynamic or not ({STATIC}, {DYNAMIC})
 *   "NAME" : name of the field
 *   "COMPONENT_NAMES" : name of the field components ({K,EPSILON} or {VITESSE_X,VITESSE_Y,VITESSE_Z})
 *   "BOUNDARY_NAMES" : name of the boundaries
 *   "COORDINATES" : coordinate system of the nodes ({X}, {X,Y}, {X,Y,Z}, {R,THETA}, etc...)
 *  Exceptions:
 *   - GenericFieldError::INVALID : query not understood
 *
 */
const Noms Champ_Generique_base::get_property(const Motcle& query) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}


/*! @brief Returns the type of geometric entities on which the discrete values are attached (NODE for a P1 field, FACE for a P1NC field, ELEMENT for a
 *
 *   P0 field, etc). It is recommended to use the get_localisation() syntax without
 *   parameter, unless you know what to do for multi-support fields.
 *   @sa get_nb_localisations()
 *  Parameter : index
 *  Meaning : index of the requested localization (for multi-support fields).
 *   If index = -1 : if the field is multi-support an exception is raised, otherwise the support is returned.
 *   If index >= 0 : returns the i-th support.
 *  Default value : -1
 *  Exceptions:
 *   - GenericFieldError::INVALID : the field is not discretized on these geometric entities
 *     (it can be an analytical field or a multi-localization field), or the field
 *     is multi-support while index = -1 was requested, or the support "index" does not exist.
 *
 */
Entity Champ_Generique_base::get_localisation(const int index) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

int Champ_Generique_base::get_nb_localisations() const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

/*! @brief Returns a reference to the array of discrete values if it exists in memory.
 *
 * The reference is valid at least until the next call to a non-const method.
 *   The mesh on which these values are defined is accessible through the methods
 *   get_node_coordinates(), get_connectivity(), get_field_type().
 *  Exceptions:
 *   - GenericFieldError::XYZ_ONLY : values are accessible only through get_xyz_values
 *   - GenericFieldError::NO_REF : values are not stored in memory,
 *      must use get_copy_values();
 *   - GenericFieldError::MESH_ONLY : the field carries no values, it serves only to describe a geometry
 *
 */
const DoubleTab& Champ_Generique_base::get_ref_values() const
{
  // Default implementation : exception NO_REF
  throw Champ_Generique_erreur("NO_REF");
}

/*! @brief Fills the values array with the discrete values of the field (creates a copy).
 *
 * Exceptions:
 *   - GenericFieldError::XYZ_ONLY : values are accessible only through get_xyz_values
 *   - GenericFieldError::MESH_ONLY : the field carries no values, it serves only to describe a geometry
 *
 */
void Champ_Generique_base::get_copy_values(DoubleTab& values) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

/*! @brief Computes the point value of the field at the coordinates given in coords and puts them in values.
 *
 * validity_flag is filled with 1 if the value is valid
 *   (coordinates inside the domain), 0 otherwise.
 *  Warning in parallel:
 *   Each processor processes the coords array provided to it: one can either have
 *   the same array computed by everyone: in this case, validity_flag indicates on
 *   each processor which values each processor was able to compute, or one knows in
 *   advance which coordinates are computable by each processor and gives the processor
 *   locally only the coordinates it owns. If the same coordinate is requested
 *   from multiple processors, there is no guarantee that all give the same result. In general,
 *   only one will have the validity_flag set for this coordinate.
 *  Exceptions:
 *   - GenericFieldError::NOT_IMPLEMENTED : the lazy one has not coded the interpolation methods
 *    must work with the mesh and the discrete values.
 *   - GenericFieldError::MESH_ONLY : the field carries no values, it serves only to describe a geometry
 *
 */
void Champ_Generique_base::get_xyz_values(const DoubleTab& coords, DoubleTab& values, ArrOfBit& validity_flag) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

/*! @brief Returns a reference to the array of coordinates of the vertices of the mesh supporting the field, if it exists.
 *
 * The array always has two dimensions:
 *    dimension(0) = number of real vertices
 *    dimension(1) = get_dimension() (dimension of the space in which the field is defined)
 *   In parallel, the array is a distributed array with common items.
 *  Exceptions:
 *   - GenericFieldError::INVALID : the array does not exist (analytical field ...)
 *   - GenericFieldError::NO_REF : values are not stored in memory,
 *      must use get_copy_coordinates();
 *
 */
const DoubleTab& Champ_Generique_base::get_ref_coordinates() const
{
  // Default implementation : exception NO_REF
  throw Champ_Generique_erreur("NO_REF");
}


void Champ_Generique_base::get_copy_coordinates(DoubleTab&) const
{
  throw Champ_Generique_erreur("NO_REF");
}

/*! @brief Returns the connectivity array between the geometric entity index1 and entity index2.
 *
 * For example
 *    get_ref_connectivity(ELEM, NODE) = Domaine::mes_elems
 *    get_ref_connectivity(ELEM, FACE) = DomaineVF::elem_faces_
 *    get_ref_connectivity(FACE, ELEM) = DomaineVF::face_voisins_
 *   The array always has two dimensions:
 *    dimension(0) = number of real entities "index1"
 *    dimension(1) = number of entities "index2" connected to each entity "index1"
 *   In parallel, the array is a distributed array with common items
 *  Exceptions:
 *   - GenericFieldError::INVALID : the array does not exist (analytical field ...)
 *   - GenericFieldError::NO_REF : values are not stored in memory,
 *      must use get_copy_connectivity();
 *
 */
const IntTab& Champ_Generique_base::get_ref_connectivity(Entity index1, Entity index2) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

void Champ_Generique_base::get_copy_connectivity(Entity index1, Entity index2, IntTab&) const
{
  throw Champ_Generique_erreur("NOT_IMPLEMENTED");
}

/*! @brief Returns a ref to the domain on which the storage space will be evaluated.
 *
 * By default the domain associated with the problem.
 *
 */
const Domaine& Champ_Generique_base::get_ref_domain() const
{
  const Objet_U& ob = interprete().objet(nom_pb_);
  const Probleme_base& pb = ref_cast(Probleme_base,ob);
  const Domaine& dom = pb.domaine_dis().domaine();
  return dom;
}


/*! @brief Creates a copy of the domain on which the storage space will be evaluated.
 *
 * By default the domain associated with the problem.
 *
 */
void Champ_Generique_base::get_copy_domain(Domaine& domain) const
{
  const Domaine& dom = get_ref_domain();
  domain = dom;
}


/*! @brief Returns a ref to the discretized domain on which the storage space will be evaluated.
 *
 *  By default the domain associated with the problem.
 *
 */
const Domaine_dis_base& Champ_Generique_base::get_ref_domaine_dis_base() const
{
  const Objet_U& ob = interprete().objet(nom_pb_);
  const Probleme_base& pb = ref_cast(Probleme_base,ob);
  const Domaine_dis_base& domaine_dis = pb.domaine_dis();
  return domaine_dis;
}

/*! @brief Returns a ref to the discretized boundary conditions domain of the equation carrying the target field.
 *
 */
const Domaine_Cl_dis_base& Champ_Generique_base::get_ref_zcl_dis_base() const
{
  throw Champ_Generique_erreur("INVALID");
}

bool Champ_Generique_base::has_champ_post(const Motcle& nom) const
{
  Motcle nom_champ;

  const Noms nom_champ_post = get_property("nom");
  nom_champ = Motcle(nom_champ_post[0]);
  if (nom_champ == nom)
    return true;

  const Noms syno = get_property("synonyms");
  for (int i = 0; i < syno.size(); i++)
    {
      nom_champ = Motcle(syno[i]);
      if (nom_champ == nom)
        return true;
    }

  const Noms composantes = get_property("composantes");
  for (const auto &itr : composantes)
    {
      nom_champ = Motcle(itr);
      if (nom_champ == nom)
        return true;
    }

  return false; /* nothing found */
}

const Champ_Generique_base& Champ_Generique_base::get_champ_post(const Motcle& nom) const
{
  OBS_PTR(Champ_Generique_base) ref_champ;

  Motcle nom_champ;
  const Noms nom_champ_post = get_property("nom");
  nom_champ = Motcle(nom_champ_post[0]);

  if (nom_champ==nom)
    {
      ref_champ = *this;
      ref_champ->fixer_identifiant_appel(nom);
      return ref_champ;
    }
  const Noms syno = get_property("synonyms");
  int nb_syno = syno.size();
  for (int i=0; i<nb_syno; i++)
    {
      nom_champ = Motcle(syno[i]);
      if (nom_champ==nom)
        {
          ref_champ = *this;
          ref_champ->fixer_identifiant_appel(nom);
          return ref_champ;
        }
    }

  {
    const Noms composantes = get_property("composantes");
    for (const auto& itr : composantes)
      {
        nom_champ = Motcle(itr);
        if (nom_champ==nom)
          {
            ref_champ = *this;
            ref_champ->fixer_identifiant_appel(nom);
            return ref_champ;
          }
      }
  }

  throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));
}

int Champ_Generique_base::comprend_champ_post(const Motcle& identifiant) const
{
  Motcle nom_champ;
  const Noms nom_champ_post = get_property("nom");
  nom_champ = Motcle(nom_champ_post[0]);

  if (nom_champ==identifiant)
    return 1;
  else
    {
      const Noms composantes = get_property("composantes");
      for (const auto& itr : composantes)
        {
          nom_champ = Motcle(itr);
          if (nom_champ==identifiant)
            return 1;
        }
    }
  return 0;
}
