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


#ifndef Objet_a_lire_included
#define Objet_a_lire_included

#include <TRUST_Deriv.h>
#include <TRUSTArray.h>
#include <TRUST_List.h>
#include <ptrParam.h>
#include <Objet_U.h>
#include <map>
#include <string>
#include <functional>

class Param;
class Objet_a_lire :public Objet_U
{
  Declare_instanciable(Objet_a_lire);
public:
  enum Type { INTEGER = 0, TRUSTID, DOUBLE, STRING, OBJECT, FLAG, NON_STD, DERIV, ArrOfInt_size_imp,
              ArrOfDouble_size_imp, PARAM, MAP_INT, MAP_DOUBLE, MAP_STRING, MAP_OBJET_U, VEC_INT, VEC_DOUBLE, VEC_STRING, VEC_OBJET_U
            };

  enum Nature { OPTIONAL = 0, REQUIRED = 1 };

  void set_entier(int*);
  void set_tid(trustIdType*);
  void set_double(double*);
  void set_string(std::string*);
  void set_objet(Objet_U*);
  void set_arrofint(ArrOfInt*);
  void set_arrofdouble(ArrOfDouble*);

  void set_vec_expected_size(int s) {expected_vec_size_=s;};
  void set_vec_int(std::vector<int>*);
  void set_vec_dbl(std::vector<double>*);
  void set_vec_str(std::vector<std::string>*);

  void set_map_int(std::map<std::string, int>*);
  void set_map_dbl(std::map<std::string, double>*);
  void set_map_str(std::map<std::string, std::string>*);

  using vec_obj_initializer_t = std::function<void(std::vector<DerObjU>&)>;
  using map_obj_initializer_t = std::function<void(std::map<std::string, DerObjU>&)>;
  void set_vec_obj_initializer(vec_obj_initializer_t);
  void set_map_obj_initializer(map_obj_initializer_t);

  template<typename _CLASSE_>
  void set_deriv(TRUST_Deriv<_CLASSE_> *quoi, const char *prefixe)
  {
    obj_a_lire = quoi;
    type = DERIV;
    prefixe_deriv = prefixe;
  }

  Param& create_param(const char*);
  void set_flag(bool*);
  void set_non_std(Objet_U*);
  ptrParam& add_dict(const char*,int, const char* =0);
  void set_nature(Objet_a_lire::Nature n);
  void read(const Motcle& keyword,Entree& is);
  void print(Sortie& s) const;
  const Nom& get_name() const;
  int comprend_name(Motcle& mot) const;
  Nom get_names_message() const;
  void set_name(const LIST(Nom)& n);
  bool is_optional() const;
  bool is_type_simple() const ;
  double get_value() const;

protected:
  Nom name;
  LIST(Nom) names;
  Objet_a_lire::Type type = INTEGER;
  Objet_a_lire::Nature nature = OPTIONAL;
  int *int_a_lire;
  trustIdType *tid_a_lire;
  double *double_a_lire;
  std::string *string_a_lire;
  Objet_U *obj_a_lire, *objet_lu;
  ArrOfInt *arrofint_a_lire;
  ArrOfDouble *arrofdouble_a_lire;
  bool *flag_a_lire;
  Motcles dictionnaire_noms;
  ArrOfInt dictionnaire_valeurs;
  LIST(ptrParam) dictionnaire_params;
  Motcle prefixe_deriv;
  ptrParam param_interne;


  int expected_vec_size_ = -1; // -1 means no size restriction.
  std::vector<int>* vec_int_a_lire = nullptr;
  std::vector<double>* vec_double_a_lire = nullptr;
  std::vector<std::string>* vec_str_a_lire = nullptr;
  vec_obj_initializer_t vec_obj_initializer;

  std::map<std::string, int>* map_int_a_lire = nullptr;
  std::map<std::string, double>* map_double_a_lire = nullptr;
  std::map<std::string, std::string>* map_str_a_lire = nullptr;
  map_obj_initializer_t map_obj_initializer;

};

#endif
