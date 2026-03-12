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

#include <Entree_complete.h>
#include <Objet_a_lire.h>
#include <Param.h>
#include <string>

Implemente_instanciable_sans_constructeur(Objet_a_lire,"Objet_a_lire",Objet_U);

Objet_a_lire::Objet_a_lire() : int_a_lire(nullptr), tid_a_lire(nullptr), double_a_lire(nullptr), obj_a_lire(nullptr), objet_lu(nullptr),
  arrofint_a_lire(nullptr), arrofdouble_a_lire(nullptr), flag_a_lire(nullptr) { }

Entree& Objet_a_lire::readOn(Entree& is)
{
  return Objet_U::readOn(is);
}

Sortie& Objet_a_lire::printOn(Sortie& os) const
{
  return Objet_U::printOn(os);
}

void Objet_a_lire::set_entier(int *quoi)
{
  type = INTEGER;
  int_a_lire = quoi;
}

void Objet_a_lire::set_tid(trustIdType *quoi)
{
  type = TRUSTID;
  tid_a_lire = quoi;
}

void Objet_a_lire::set_double(double *quoi)
{
  type = DOUBLE;
  double_a_lire = quoi;
}

void Objet_a_lire::set_string(std::string *quoi)
{
  type = STRING;
  string_a_lire = quoi;
}

void Objet_a_lire::set_objet(Objet_U *quoi)
{
  obj_a_lire = quoi;
  type = OBJECT;
}

void Objet_a_lire::set_arrofint(ArrOfInt *quoi)
{
  arrofint_a_lire = quoi;
  type = ArrOfInt_size_imp;
}

void Objet_a_lire::set_arrofdouble(ArrOfDouble *quoi)
{
  arrofdouble_a_lire = quoi;
  type = ArrOfDouble_size_imp;
}


void Objet_a_lire::set_vec_int(std::vector<int>*quoi)
{
  type = VEC_INT;
  vec_int_a_lire = quoi;
}
void Objet_a_lire::set_vec_dbl(std::vector<double>*quoi)
{
  type = VEC_DOUBLE;
  vec_double_a_lire = quoi;
}
void Objet_a_lire::set_vec_str(std::vector<std::string>*quoi)
{
  type = VEC_STRING;
  vec_str_a_lire = quoi;
}

void Objet_a_lire::set_vec_obj_initializer(vec_obj_initializer_t initializer)
{
  type = VEC_OBJET_U;
  vec_obj_initializer = initializer;
}

void Objet_a_lire::set_map_int(std::map<std::string, int>*quoi)
{
  type = MAP_INT;
  map_int_a_lire = quoi;
}
void Objet_a_lire::set_map_dbl(std::map<std::string, double>*quoi)
{
  type = MAP_DOUBLE;
  map_double_a_lire = quoi;
}
void Objet_a_lire::set_map_str(std::map<std::string, std::string>*quoi)
{
  type = MAP_STRING;
  map_str_a_lire = quoi;
}

void Objet_a_lire::set_map_obj_initializer(map_obj_initializer_t initializer)
{
  type = MAP_OBJET_U;
  map_obj_initializer = initializer;
}

void Objet_a_lire::set_flag(bool *quoi)
{
  flag_a_lire = quoi;
  // initialisation du flag a false
  *flag_a_lire = false;
  type = FLAG;
}

void Objet_a_lire::set_non_std(Objet_U *quoi)
{
  objet_lu = quoi;
  type = NON_STD;
}
Param& Objet_a_lire::create_param(const char *aname)
{
  type = PARAM;
  param_interne.create(aname);
  return param_interne.valeur();
}

bool Objet_a_lire::is_optional() const
{
  return (nature == OPTIONAL);
}

void Objet_a_lire::set_name(const LIST(Nom) &n)
{
  const auto& list = n.get_stl_list();
  auto itr = list.begin();
  name = *itr;
  ++itr;

  for ( ; itr != list.end(); ++itr) names.add(*itr);
}

int Objet_a_lire::comprend_name(Motcle& mot) const
{
  if (mot == name) return 1;
  for (const auto& itr : names)
    if (mot == itr)
      {
        mot = name;
        return 1;
      }
  return 0;
}

Nom Objet_a_lire::get_names_message() const
{
  Nom titi(name);
  int size = names.size();
  for (int i = size - 1; i >= 0; i--)
    {
      titi += Nom("|");
      titi += names(i);
    }
  return titi;
}

const Nom& Objet_a_lire::get_name() const
{
  return name;
}

ptrParam& Objet_a_lire::add_dict(const char *nom_option, int valeur, const char *aname)
{
  // Le dictionnaire ne fonctionne que pour des parametres de type int:
  assert(int_a_lire != 0);
  // L'option ne doit pas encore exister dans le dictionnaire:
  assert(dictionnaire_noms.search(nom_option) < 0);

  dictionnaire_noms.add(nom_option);

  dictionnaire_valeurs.append_array(valeur);
  dictionnaire_params.add(ptrParam());
  ptrParam& ptr = dictionnaire_params[dictionnaire_params.size() - 1];
  if (aname != 0)
    {
      ptr.create(aname);
    }
  return ptr;
}
void Objet_a_lire::set_nature(Objet_a_lire::Nature n)
{
  nature = n;
}

bool Objet_a_lire::is_type_simple() const
{
  return ((type == INTEGER) || (type == DOUBLE) || (type == FLAG));
}

double Objet_a_lire::get_value() const
{
  switch(type)
    {
    case INTEGER:
      return (*int_a_lire);
    case DOUBLE:
      return (*double_a_lire);
    case FLAG:
      return (*flag_a_lire);
    default:
      Cerr << "get_value not coded for this case" << finl;
      Process::exit();
      break;
    }
  // pour les compilos
  return 0.;
}




namespace
{

template<typename T>
T convert(std::string val);

template <> int convert<int>(std::string val)
{
  return std::stoi(val);
}
template <> double convert<double>(std::string val)
{
  return std::stod(val);
}
template <> std::string convert<std::string>(std::string val)
{
  return val;
}

void expect_word(Entree& is, std::string word)
{
  std::string read;
  is>> read;
  if (read != word)
    {
      Cerr << Nom("From expect_word in Objet_a_lire.cpp: Expected to read: ") + word << endl;
      Cerr << Nom("Found: ") + read << endl;
      Process::exit();
    }

}

template<typename T>
std::vector<T> read_vec_impl(
  Motcle const& motcle,
  Entree& is,
  std::function<T(Entree&, std::string /*first_token*/)> read_value,
  int expected_size = -1)
{
  std::vector<T> vec;

  std::string first;
  is >> first;

  if (first == "[")
    {
      // Bracket syntax: [ a, b, c ]
      std::string token;
      is >> token;

      while (token != "]")
        {
          try
            {
              vec.push_back(read_value(is, token));
            }
          catch (std::exception& e)
            {
              Cerr << "In keyword " << motcle
                   << ", invalid value led to exception: " << e.what() << endl;
              Process::exit();

            }
          is >> token;

          if (token == ",")
            {
              is >> token; // either next element or "]"
            }
          else if (token != "]")
            {
              Cerr << "In keyword " << motcle
                   << ", expected ',' or ']', found: " << token << endl;
              Process::exit();
            }
        }
    }
  else
    {
      // Size-prefixed syntax: N a b c ...
      int size=-1; // for non init warning. this default value CANNOT be used
      try
        {
          size = std::stoi(first);
        }
      catch (...)
        {
          Cerr << "In keyword " << motcle
               << ", expected '[' or an integer size, found: " << first << endl;
          Process::exit();
        }

      if (expected_size >= 0 && size != expected_size)
        {
          Cerr << "In keyword " << motcle << ", declared size " << size
               << " does not match expected size " << expected_size << endl;
          Process::exit();
        }

      vec.reserve(static_cast<size_t>(size));
      for (int i = 0; i < size; i++)
        {
          std::string token;
          is >> token;
          try
            {
              vec.push_back(read_value(is, token));
            }
          catch (std::exception& e)
            {
              Cerr << "In keyword " << motcle
                   << ", invalid value led to exception: " << e.what() << endl;
              Process::exit();

            }
        }
    }

  // Check final size against expected_size (applies to bracket syntax too)
  if (expected_size >= 0 && (int)vec.size() != expected_size)
    {
      Cerr << "In keyword " << motcle << ", read " << vec.size()
           << " elements but expected " << expected_size << endl;
      Process::exit();
    }

  return vec;
}

template<typename T>
std::vector<T> read_vec_any(Motcle const& motcle, Entree& is, int expected_size = -1)
{
  return read_vec_impl<T>(motcle, is,
                          [](Entree& /*is*/, std::string first_token) -> T
  {
    return convert<T>(first_token);
  }, expected_size);
}

std::vector<DerObjU> read_vec_derObjU(Motcle const& motcle, Entree& is, int expected_size = -1)
{
  return read_vec_impl<DerObjU>(motcle, is,
                                [](Entree& is_in, std::string type_str) -> DerObjU
  {
    DerObjU obj;
    obj.typer(type_str.c_str());
    is_in >> obj.valeur();
    return obj;
  }, expected_size);
}

template<typename T>
std::map<std::string, T> read_map_impl(
  Motcle const& motcle,
  Entree& is,
  std::function<T(Entree&, std::string /*first_token*/)> read_value)
{
  std::map<std::string, T> map;

  expect_word(is, "{");

  std::string key;
  is >> key;

  // case of an empty map
  if (key == "}") {return map;}

  bool with_colon_and_comma = false;

  std::string first_token;
  is >> first_token;

  bool has_colon = (first_token == ":");
  if (has_colon)
    {
      with_colon_and_comma = true;
      std::string actual_first;
      is >> actual_first;
      first_token = actual_first;
    }

  T val = read_value(is, first_token);

  while (true)
    {
      if (map.contains(key))
        {
          Cerr << "Duplicate key '" << key << "' found for param: " << motcle << endl;
          Process::exit();
        }
      map[key] = val;

      if (with_colon_and_comma)
        {
          // At this point, the next word should be either a comma or a closing brace
          std::string sep;
          is >> sep;

          // because comma is not needed at the end, we do that:
          // if after a 'key : value' pair, there is a closing brace, stop
          // if there is a comma, read next word into key
          // if it was none of them, syntax is not valid
          // then, if after this the key is a closing brace, stop reading
          // this allows trailing comma at the end

          if (sep == "}") break;
          else if (sep == ",") is >> key;
          else
            {
              Cerr << "In keyword " << motcle
                   << ", expected to read either a comma (,) or closing brace (}). Found "
                   << sep << endl;
              Process::exit();
            }

          if (key == "}") break;

          // after a key, there must be a colon
          expect_word(is, ":");

          // then we can read a value and fill the map
          std::string next_first;
          is >> next_first;
          try
            {
              val = read_value(is, next_first);
            }
          catch (std::exception& e)
            {
              Cerr << "In keyword " << motcle
                   << ", invalid value led to exception: " << e.what() << endl;
              Process::exit();

            }
        }
      else
        {
          // this case is much simpler: less chars, no handling of trailing commas
          is >> key;
          if (key == "}") break;

          std::string next_first;
          is >> next_first;

          try
            {
              val = read_value(is, next_first);
            }
          catch (std::exception& e)
            {
              Cerr << "In keyword " << motcle
                   << ", invalid value led to exception: " << e.what() << endl;
              Process::exit();

            }
        }
    }

  return map;
}
template<typename T>
std::map<std::string, T> read_map_any(Motcle const& motcle, Entree& is)
{
  return read_map_impl<T>(motcle, is,
                          [](Entree& /*is*/, std::string first_token) -> T
  {
    return convert<T>(first_token);
  });
}

std::map<std::string, DerObjU> read_map_derObjU(Motcle const& motcle, Entree& is)
{
  return read_map_impl<DerObjU>(motcle, is,
                                [](Entree& is_in, std::string type_str) -> DerObjU
  {
    DerObjU obj;
    obj.typer(type_str.c_str());
    is_in >> obj.valeur();
    return obj;
  });
}

}




void Objet_a_lire::read(Motcle const& motcle, Entree& is)
{
  int ret = -1;
  switch(type)
    {
    case INTEGER:
      if (dictionnaire_noms.size() == 0)
        {
          is >> *int_a_lire;
        }
      else
        {
          Motcle motlu;
          is >> motlu;
          const int rang = dictionnaire_noms.search(motlu);
          if (rang < 0)
            {
              Cerr << "Error while reading parameter " << name << "\n Found: " << motlu << "\n Expected one on these keywords: " << dictionnaire_noms << finl;
              barrier();
              exit();
            }
          *int_a_lire = dictionnaire_valeurs[rang];
          ptrParam& ptr = dictionnaire_params[rang];
          if (ptr.non_nul())
            {
              ptr->lire_avec_accolades_depuis(is);
            }
        }
      break;
    case TRUSTID:
      is >> (*tid_a_lire);
      break;
    case DOUBLE:
      is >> (*double_a_lire);
      break;
    case STRING:
      is >> (*string_a_lire);
      break;
    case OBJECT:
      is >> (*obj_a_lire);
      break;
    case FLAG:
      (*flag_a_lire) = true;
      break;
    case NON_STD:
      ret = (*objet_lu).lire_motcle_non_standard(motcle, is);
      if (ret < 0)
        {
          Cerr << "Error while reading keyword: '" << motcle << "'. Not recognized!" << finl;
          Process::exit(-1);
        }
      break;
    case DERIV:
      {
        Motcle type_complet(prefixe_deriv);
        Motcle the_type;
        is >> the_type;
        type_complet += the_type;
        Entree_complete is2(type_complet, is);
        is2 >> (*obj_a_lire);
        break;
      }
    case ArrOfInt_size_imp:
      {
        ArrOfInt& arr = *arrofint_a_lire;
        int size = arr.size_array();
        // on ne veut pas de taille nulle pour etre sur que l'utilisateur a bien fait dimensionner...
        assert(size > 0);
        for (int i = 0; i < size; i++)
          is >> arr[i];
        break;
      }
    case ArrOfDouble_size_imp:
      {
        ArrOfDouble& arr = *arrofdouble_a_lire;
        int size = arr.size_array();
        // on ne veut pas de taille nulle pour etre sur que l'utilisateur a bien fait dimensionner...
        assert(size > 0);
        for (int i = 0; i < size; i++)
          is >> arr[i];
        break;
      }
    case PARAM:
      {
        param_interne->lire_avec_accolades_depuis(is);
        break;
      }
    case VEC_INT:
      {
        *vec_int_a_lire = read_vec_any<int>(motcle, is, expected_vec_size_);
        break;
      }
    case VEC_DOUBLE:
      {
        *vec_double_a_lire = read_vec_any<double>(motcle, is, expected_vec_size_);
        break;
      }
    case VEC_STRING:
      {
        *vec_str_a_lire = read_vec_any<std::string>(motcle, is, expected_vec_size_);
        break;
      }
    case VEC_OBJET_U:
      {
        auto read = read_vec_derObjU(motcle, is, expected_vec_size_);
        vec_obj_initializer(read);
        break;
      }
    case MAP_INT:
      {
        *map_int_a_lire = read_map_any<int>(motcle, is);
        break;
      }
    case MAP_DOUBLE:
      {
        *map_double_a_lire = read_map_any<double>(motcle, is);
        break;
      }
    case MAP_STRING:
      {
        *map_str_a_lire = read_map_any<std::string>(motcle, is);
        break;
      }
    case MAP_OBJET_U:
      {
        auto read = read_map_derObjU(motcle, is);
        map_obj_initializer(read);
        break;
      }
    default:
      Cerr << "Invalid parameter type" << finl;
      Process::exit();
      break;
    }
}

