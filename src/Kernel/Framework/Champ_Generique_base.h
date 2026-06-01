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

#ifndef Champ_Generique_base_included
#define Champ_Generique_base_included

#include <TRUST_Deriv.h>
#include <TRUST_Ref.h>
#include <Domaine.h>
#include <ArrOfBit.h>
#include <YAML_data.h>

class Postraitement_base;
class Discretisation_base;
class Domaine_Cl_dis_base;
class Param;
class Champ_base;

enum class Entity { NODE, SEGMENT, FACE, ELEMENT };

/*! @brief class Champ_Generique_base
 *
 * Base class of generic fields for importing a discrete field and elementary actions on this field
 *
 *   (post-processing, etc)
 *  Warning: all methods are PARALLEL, they must be called
 *   simultaneously on all processors (get_domain() can for example
 *   build the parallel domain before returning it).
 *
 */
class Champ_Generique_base : public Objet_U
{
  Declare_base(Champ_Generique_base);
public:

  /* XXX Elie Saikali : put it here and not in Objet_U */
  virtual std::vector<YAML_data> data_a_sauvegarder() const { return std::vector<YAML_data>(); };
  int sauvegarder(Sortie& os) const override { return 0; }
  int reprendre(Entree& is) override { return 1; }

  virtual void set_param(Param& param) const override=0;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  virtual int  get_dimension() const;
  virtual double  get_time() const; //returns the time of the encapsulated field
  virtual const   Probleme_base& get_ref_pb_base() const;
  virtual const   Discretisation_base&  get_discretisation() const;
  virtual const   Motcle                get_directive_pour_discr() const;

  void nommer(const Nom& nom) override;
  virtual const Nom& get_nom_post() const;
  //returns -1 if identifier is the name, component number otherwise
  static inline int composante(const Nom& nom_test,const Nom& nom, const Noms& composantes,const Noms& synonyms);

  virtual void              get_property_names(Motcles& list) const;
  virtual const Noms        get_property(const Motcle& query) const;  //answers requests for name, target_name, units and components

  virtual Entity            get_localisation(const int index = -1) const; //localization ELEMENT or NODE are post-processable
  virtual int               get_nb_localisations() const;

  virtual const DoubleTab&  get_ref_values() const;
  virtual void              get_copy_values(DoubleTab&) const;
  virtual void              get_xyz_values(const DoubleTab& coords, DoubleTab& values, ArrOfBit& validity_flag) const;

  virtual const Domaine&    get_ref_domain() const; //returns a reference to the domain associated with the field
  virtual void              get_copy_domain(Domaine&) const;
  virtual const Domaine_dis_base&  get_ref_domaine_dis_base() const; //returns the discretized domain linked to the domain
  virtual const Domaine_Cl_dis_base&  get_ref_zcl_dis_base() const; //returns the discretized boundary conditions linked to the equation carrying the target field

  virtual const DoubleTab&  get_ref_coordinates() const;
  virtual void              get_copy_coordinates(DoubleTab&) const;
  virtual const IntTab&     get_ref_connectivity(Entity index1, Entity index2) const;
  virtual void              get_copy_connectivity(Entity index1, Entity index2, IntTab&) const;

  // Resets the object to the state obtained by the default constructor
  virtual void  reset() = 0;
  virtual void completer(const Postraitement_base& post) = 0; //Completes the operator possibly carried by the field
  //and names the sources by default
  virtual void  mettre_a_jour(double temps) = 0;              //Updates an operator possibly carried by the field

  //  The get_champ() method is in particular called by the post-processing class when
  //  the field must be post-processed (dt_post elapsed).
  //  The caller must provide an untyped "espace_stockage" field as parameter
  //   "espace_stockage".
  //
  //  Either it returns an existing field (see Champ_Generique_refChamp)
  //   and does not use espace_stockage.
  //  Or it builds a new field stored in espace_stockage,
  //   and the return value is espace_stockage.
  //  The caller retrieves the computation result in the return value,
  //   knowing that it may possibly reference espace_stockage
  //   (so, do not destroy espace_stockage too early).

  // The steps to create the storage space are:
  // espace_stockage.typer(type_champ)
  // espace_stockage.associer_domaine_dis_base(un_domaine_dis)
  // espace_stockage.fixer_nb_comp(nb_comp);
  // espace_stockage.fixer_nb_valeurs_nodales(nb_ddl);
  // Computation of values by an instruction of the form
  // espace_stockage.valeurs() = Operateur.calculer(source.valeurs())
  // espace_stockage.valeurs().echange_espace_virtuel()
  //return espace_stockage

  virtual const Champ_base& get_champ(OWN_PTR(Champ_base) &espace_stockage) const = 0;
  virtual const Champ_base& get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const=0;

  //get_champ_post() returns the field if the identifier passed as parameter designates
  //the field name or one of its components
  virtual const Champ_Generique_base& get_champ_post(const Motcle& nom) const;
  virtual bool has_champ_post(const Motcle& nom) const;

  //returns 1 if the field is identified, 0 otherwise
  virtual int comprend_champ_post(const Motcle& identifiant) const;

  //get_info_type_post() returns 0 if a table must be post-processed, 1 for a tensor
  virtual int get_info_type_post() const = 0;

  //sets the identifiant_appel_ attribute of the field to indicate whether the request
  //was launched by the field name or one of its components (cf Champ_Generique_Interpolation)
  inline void fixer_identifiant_appel(const Nom& identifiant)
  {
    identifiant_appel_ = identifiant;
  }

  virtual const Noms& fixer_noms_compo(const Noms& noms)
  {
    Cerr << "The method " << __func__ << " is not overloaded in " << que_suis_je() << finl;
    throw;
  }
  virtual const Noms& fixer_noms_synonyms(const Noms& noms)
  {
    Cerr << "The method " << __func__ << " is not overloaded in " << que_suis_je() << finl;
    throw;
  }

protected:
  static void               assert_parallel(int);

  Nom nom_post_;
  Nom identifiant_appel_;
  Nom nom_pb_;
  OBS_PTR(Probleme_base) ref_pb_;
};

inline int Champ_Generique_base::composante(const Nom& nom_test,const Nom& nom,const Noms& composantes, const Noms& syno)
{
  Motcle motlu(nom_test);
  if (motlu == Motcle(nom))
    return -1;
  for (auto& itr : syno)
    if (itr==motlu) return -1;
  int n = composantes.size();
  Motcles les_noms_comp(n);
  for (int i=0; i<n; i++)
    les_noms_comp[i] = composantes[i];
  int ncomp = les_noms_comp.search(motlu);
  if (ncomp == -1)
    {
      Cerr << "Error TRUST, the identifier : " << nom_test << finl
           << "did not designate the name of the field tested nor any of its components\n";
      Cerr <<" the components of field named "<<nom <<" are :"<< finl;
      for (int ii=0; ii<n; ii++)
        Cerr << composantes[ii] << " ";
      Cerr << finl;
      exit();
    }
  return ncomp;
}

/*
 * @brief class Champ_Generique_erreur
 *
 * Class Champ_Generique_erreur
 */
class Champ_Generique_erreur
{
public:
  Nom mot1;

  Champ_Generique_erreur(const Nom& mot2)
  {
    mot1 = mot2;
    Cerr<<"Error of type : "<<mot1<<finl;
    Process::exit();
  }
};

#endif /* Champ_Generique_base_included */
