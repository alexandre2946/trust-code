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

#ifndef Format_Post_base_included
#define Format_Post_base_included

#include <TRUSTTabs_forward.h>
#include <TRUST_Deriv.h>
#include <Champ_base.h>
#include <TRUST_Ref.h>
#include <Domaine.h>

class Domaine_dis_base;
class Motcle;
class Param;

/*! @brief Base class for post-processing output formats for fields (lata, med, cgns, lml, single_lata).
 *
 *  Using the class through the generic interface:
 *  - type a post-processing format object
 *  - initialize the post-processing format object
 *  - perform completion and verification operations
 *  - write a domain
 *  - ecrire_temps
 *  - ecrire_champ (supported field types: at elements or at nodes,
 *                  or at faces if the domain faces have been written
 *                  and if faces are supported by the post-processing format)
 *
 */

// Parameters involved (possibly) in a post-processing operation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Nom                    nom_fichier                 : name of the file being written
// int           champ                  : 1 if field post-processing, 0 otherwise
// int            stat                          : 1 if statistics post-processing, 0 otherwise
// double            dt_ch                  : period for field post-processing
// double            dt_stat                  : period for statistics post-processing
// int            reprise                : 1 if restart was performed, 0 otherwise
// int            axi                        : 1 if axisymmetric calculation, 0 otherwise
// int            est_le_premier_post         : 1 if first post-processing for a given output file, 0 otherwise
// int            est_le_dernier_post  : 1 if last post-processing for a given output file, 0 otherwise
// Domaine            dom                        : computation domain
// Nom                   id_domaine                : domain name
// IntVect           faces_som                : face-to-vertex connectivity
// IntVect           elem_faces                : element-to-face connectivity
// int           nb_som                : number of domain vertices
// int            nb_faces                : number of domain faces

// double            t_init                : initial time of the computation
// double            temps_courant          : current time of the computation
// double           temps_champ                : target field time
// double            temps_post                : time of the post-processing object (write if temps_post<temps_courant)

// Nom                    id_champ_post        : identifier of the generic field (name or component)
// Nature_du_champ nature_champ                : scalar, multi_scalar, vector
// int            nb_compo                : number of field components
// Noms            nom_compos                : vector containing the component names of the field
// Noms            unites                : vector containing the units of the field to write
// Motcle            loc_post                : location: ELEM, SOM, FACES
// int            ncomp                : component number (-1 if not a component)
// DoubleTab           data                        : array of values to write to file
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Note: data
// If post-processing the full field and not a component of a field:
//        If the field value array has a single entry, data must be dimensioned as data(nb_ddl,1)
// If post-processing a component of a field:
//        If the supplied value array has multiple entries, data must be dimensioned as data(nb_ddl)


// Using the generic interface to perform a post-processing operation.
// The sequence of interface method calls to post-process an array of values of a discrete field at a given instant
// is presented below (not every method is necessarily required depending on the format used)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// format_post.initialize_by_default(nom_fichier)
// format_post.ecrire_entete(temps_courant,reprise,est_le_premier_post)
// format_post.preparer_post(id_domaine,est_le_premier_post,reprise,t_init)
// format_post.completer_post(dom,axi,nature_champ,nb_compo,noms_compo,loc_post,id_champ_post)
// format_post.ecrire_domaine(dom,est_le_premier_post)
// format_post.ecrire_item_int("FACES",id_domaine,id_domaine,"FACES","SOMMETS",faces_som,nb_som);
// format_post.ecrire_item_int("ELEM_FACES",id_domaine,id_domaine,"ELEMENTS","FACES",elem_faces,nb_faces);
// format_post.ecrire_temps(temps_courant)
// format_post.init_ecriture(temps_courant,temps_post,est_le_premier_post,dom)
// format_post.ecrire_champ(dom,unites,noms_compo,ncomp,temps_champ,temps_courant,id_champ_post,id_domaine,loc_post,data)
// format_post.finir_ecriture(temps_courant)
// format_post.finir(est_le_dernier_post)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Format_Post_base : public Objet_U
{
  Declare_base(Format_Post_base);
public:
  // Resets the object to the initial state obtained after the default constructor
  virtual void reset() = 0;
  virtual void resetTime(double t, const std::string dirname);
  virtual void set_param(Param& param) const override=0;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  virtual int initialize_by_default(const Nom& file_basename);
  virtual int initialize(const Nom& file_basename, const int format, const Nom& option_para);
  virtual int modify_file_basename(const Nom file_basename, bool for_restart, const double tinit);
  virtual int ecrire_entete(const double temps_courant, const int reprise, const int est_le_premier_post);
  virtual int completer_post(const Domaine& dom, const int axi, const Nature_du_champ& nature, const int nb_compo, const Noms& noms_compo, const Motcle& loc_post, const Nom& le_nom_champ_post);

  virtual int preparer_post(const Nom& id_du_domaine, const int est_le_premier_post, const int reprise, const double t_init);

  virtual int init_ecriture(double temps_courant, double temps_post, int est_le_premier_postraitement_pour_nom_fich_, const Domaine& domaine);

  virtual int finir_ecriture(double temps_courant);
  virtual int finir(const int est_le_dernier_post);

  virtual int ecrire_domaine(const Domaine& domaine,const int est_le_premier_post);
  virtual void ecrire_domaine_dual(const Domaine& domaine,const int est_le_premier_post);
  virtual int ecrire_domaine_dis(const Domaine& domaine,const OBS_PTR(Domaine_dis_base)& domaine_dis_base,const int est_le_premier_post);
  virtual int ecrire_temps(const double temps);

  virtual int ecrire_champ(const Domaine& domaine, const Noms& unite_, const Noms& noms_compo, int ncomp, double temps_, const Nom& id_du_champ, const Nom& id_du_domaine,
                           const Nom& localisation, const Nom& nature, const DoubleTab& data);

  virtual int ecrire_item_int(const Nom& id_item, const Nom& id_du_domaine, const Nom& id_domaine, const Nom& localisation, const Nom& reference, const IntVect& data, const int reference_size);

  virtual void set_single_lata_option(const bool ) { /* Do nothing */ }
  virtual void set_postraiter_domain() { /* Do nothing */ }
  virtual void set_deformable_domain() { /* Do nothing */ }
  virtual void set_lagrangian_domain() { /* Do nothing */ }
  virtual void set_discr_type(const Nom&) { /* Do nothing */ }
  virtual void set_loc_vector(const std::vector<std::string>&) { /* Do nothing */ }

protected:
  OBS_PTR(Domaine_dis_base) domaine_dis_; ///< Reference to the discretized domain - used for face fields.
};

#endif /* Format_Post_base_included */
