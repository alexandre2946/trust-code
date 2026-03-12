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

#ifndef Param_included
#define Param_included

#include <TRUSTTabs_forward.h>
#include <Objet_a_lire.h>

class Param;
class Entree;
class Motcle;
class ptrParam;     // Defined below in this file

/*! @brief cette classe permet de factoriser les readOn des Objet_U.
 *
 * On peut lui ajouter un int, un double, un Objet_U, un flag
 *   associer a un nom. ajouter_non_std permet d'appeler la methode
 *   lire_motcle_non_standard de l'objet passe en parametre.
 *   To give synonyms or translation for parameters, use the following syntax:
 *   ajouter("keyword|french_keyword|japan_keyword",...)
 *   Example: Schema_Temps_base.pp
 *
 */
class Param
{
public:
  enum Nature { OPTIONAL = 0, REQUIRED = 1 };
  Param(const char *);
  // ajout d'argument
  void ajouter(const char *, const int* ,Param::Nature nat = Param::OPTIONAL);
#if INT_is_64_ == 2
  void ajouter(const char *, const trustIdType* ,Param::Nature nat = Param::OPTIONAL);
#endif
  void ajouter(const char *, const double* ,Param::Nature nat = Param::OPTIONAL);
  void ajouter(const char *, const std::string* ,Param::Nature nat = Param::OPTIONAL);

  void ajouter(const char *, const Objet_U* ,Param::Nature nat = Param::OPTIONAL);
  void ajouter_arr_size_predefinie(const char *, const ArrOfInt* ,Param::Nature nat = Param::OPTIONAL);
  void ajouter_arr_size_predefinie(const char *, const ArrOfDouble* ,Param::Nature nat = Param::OPTIONAL);

  // vectors
  void ajouter(const char *, const std::vector<int>*, Param::Nature nat = Param::OPTIONAL ,int size = -1);
  void ajouter(const char *, const std::vector<double>*,Param::Nature nat = Param::OPTIONAL ,int size = -1);
  void ajouter(const char *, const std::vector<std::string>*, Param::Nature nat = Param::OPTIONAL, int size = -1);

  template<typename T>
  void ajouter(const char *, const std::vector<TRUST_Deriv<T>>*, Param::Nature nat = Param::OPTIONAL, int size = -1);

  // maps
  void ajouter(const char *, const std::map<std::string, int>* ,Param::Nature nat = Param::OPTIONAL);
  void ajouter(const char *, const std::map<std::string, double>* ,Param::Nature nat = Param::OPTIONAL);
  void ajouter(const char *, const std::map<std::string, std::string>* ,Param::Nature nat = Param::OPTIONAL);

  template<typename T>
  void ajouter(const char *, const std::map<std::string, TRUST_Deriv<T>>* ,Param::Nature nat = Param::OPTIONAL);


  void ajouter_flag(const char *, const bool*);
  Param& ajouter_param(const char *, Param::Nature nat = Param::OPTIONAL);
  void ajouter_non_std(const char *,const Objet_U* ,Param::Nature nat = Param::OPTIONAL);
  void ajouter_condition(const char* condition, const char* message,const char*  name=0);
  void supprimer(const char *);
  void supprimer_condition(const char*  name);
  void dictionnaire(const char *, int);
  Param& dictionnaire_param(const char *, int);

  // ajout type (pour eli)
  inline void ajouter_int(const char * c, const int* val,Param::Nature nat = Param::OPTIONAL) { ajouter(c,val,nat); }
  inline void ajouter_double(const char * c, const double* val,Param::Nature nat = Param::OPTIONAL) { ajouter(c,val,nat); }
  inline void ajouter_objet(const char *c , const Objet_U* obj,Param::Nature nat = Param::OPTIONAL) { ajouter(c,obj,nat); }

  int lire_avec_accolades_depuis(Entree& is);
  int lire_sans_accolade(Entree& is);
  inline int lire_avec_accolades(Entree& is) { return lire_avec_accolades_depuis(is); }

  int read(Entree& is,int with_acco=1);
  void print(Sortie& s) const;

  inline const LIST(Nom)& get_list_mots_lus() const { return list_parametre_lu_ ; }

  double get_value(const Nom& mot_lu) const;
  int check();
protected:
  Param();
  Objet_a_lire& create_or_get_objet_a_lire(const char *);

  LIST(Objet_a_lire) list_parametre_a_lire_;
  Nom proprietaire_;
  LIST(Nom) list_parametre_lu_,list_conditions_,list_message_erreur_conditions_, list_nom_conditions_;

};


template<typename T>
void Param::ajouter(const char * mot, const std::vector<TRUST_Deriv<T>>* quoi, Param::Nature nat ,int size)
{

  Objet_a_lire& obj = create_or_get_objet_a_lire(mot);

  auto natc = nat == Param::REQUIRED ? Objet_a_lire::REQUIRED :  Objet_a_lire::OPTIONAL;
  obj.set_nature(natc);
  obj.set_vec_expected_size(size);
  auto ptr = const_cast<std::vector<TRUST_Deriv<T>>*>(quoi);
  // captured by copy for error msg
  std::string attr_name = mot;
  std::string prop = proprietaire_.getString();


  // lambda that will set the values of objects in the map
  auto vec_initializer = [ptr, attr_name, prop](std::vector<DerObjU>& vec)
  {
    for (const auto& ref: vec)
      {

        if (sub_type(T, ref.valeur()))
          {
            const T& cast_obj = ref_cast(T, ref.valeur());
            ptr->push_back(cast_obj);
          }
        else
          {
            Cerr <<"When reading '" << prop << "'" << finl;
            Cerr <<"In keyword '" << attr_name << "', wrong type in vector:" << finl;
            Cerr <<ref.valeur().le_type() << " is not a subtype of " << T::info_obj.name() << finl;
            Process::exit();
          }

      }
  };
  obj.set_vec_obj_initializer(vec_initializer);
}

template<typename T>
void Param::ajouter(const char * mot, const std::map<std::string, TRUST_Deriv<T>>* quoi ,Param::Nature nat)
{

  Objet_a_lire& obj = create_or_get_objet_a_lire(mot);

  auto natc = nat == Param::REQUIRED ? Objet_a_lire::REQUIRED :  Objet_a_lire::OPTIONAL;
  obj.set_nature(natc);
  auto ptr = const_cast<std::map<std::string, TRUST_Deriv<T>>*>(quoi);

  // captured by copy for error msg
  std::string attr_name = mot;
  std::string prop = proprietaire_.getString();

  // lambda that will set the values of objects in the map
  auto map_initializer = [ptr, attr_name, prop](std::map<std::string, DerObjU>& map)
  {
    for (const auto& [key, o]: map)
      {

        if (sub_type(T, o.valeur()))
          {
            const T& cast_obj = ref_cast(T, o.valeur());
            (*ptr)[key] = cast_obj;
            // name the object with the map key
            (*ptr)[key]->nommer(key);
          }
        else
          {
            Cerr <<"When reading '" << prop << "'" << finl;
            Cerr <<"In keyword '" << attr_name << "', wrong type at key " <<  key << finl;
            Cerr <<o.valeur().le_type() << " is not a subtype of " << T::info_obj.name() << finl;
            Process::exit();
          }

      }
  };
  obj.set_map_obj_initializer(map_initializer);
}

#endif
