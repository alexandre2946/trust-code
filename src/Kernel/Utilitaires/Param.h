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

/*!
 * @brief Helper class to factorize the `readOn` method of `Objet_U` classes.
 *
 * `Param` provides a declarative, keyword-based way to parse a block of parameters from an
 * input stream (typically a `.data` file enclosed in `{ ... }`). Instead of manually reading
 * tokens in a `readOn` method, the developer registers each expected keyword along with a
 * pointer to the variable that should receive its value, then calls
 * `lire_avec_accolades_depuis()` to perform the actual parsing.
 *
 * Supported registrations include:
 *   - simple scalar types: `int`, `double`, `std::string`, `bool` (flag)
 *   - `Objet_U` references (standard and non-standard via `lire_motcle_non_standard`)
 *   - fixed-size arrays: `ArrOfInt`, `ArrOfDouble` (size must be set before reading)
 *   - STL containers: `std::vector<T>` and `std::map<std::string, T>` with
 *     `T` being `int`, `double`, `std::string`, or `TRUST_Deriv<U>` (see below
 *     for the exact syntaxes supported)
 *   - nested `Param` blocks (see `ajouter_param`)
 *   - dictionaries mapping string options to integer values (see `dictionnaire`)
 *   - logical post-conditions evaluated after the read (see `ajouter_condition`)
 *
 * Each registered keyword can be marked as `OPTIONAL` (default) or `REQUIRED`. Missing
 * required keywords cause the read to fail in `check()`.
 *
 * Synonyms (translations) for a keyword can be declared by separating alternative spellings
 * with a pipe character. For example:
 * @code
 *   param.ajouter("keyword|french_keyword|japan_keyword", &my_var);
 * @endcode
 * Any of the three spellings will be accepted in the input file. See `Schema_Temps_base.cpp`
 * for a real-world example.
 *
 * @par Typical usage
 * @code
 *   // Inside some Objet_U::readOn(Entree& is)
 *   Param param(que_suis_je());
 *   int a = 0, b = 0;
 *   param.ajouter("a", &a, Param::REQUIRED);
 *   param.ajouter("b", &b);                // optional
 *   param.lire_avec_accolades_depuis(is);  // expects { a 1 b 2 }
 * @endcode
 *
 * @par Vector and map syntaxes
 * Every `ajouter` overload that registers a `std::vector<T>` or a
 * `std::map<std::string, T>` (with `T` being `int`, `double`, `std::string` or
 * `TRUST_Deriv<U>`) accepts two interchangeable syntaxes. Spaces are always
 * required between all tokens.
 *
 *   - **Historical TRUST syntax**
 *     - Vectors declare the element count first, followed by the elements:
 *       @code
 *         <keyword> 2 v1 v2
 *       @endcode
 *     - Maps are enclosed in `{}` with alternating key/value pairs:
 *       @code
 *         <keyword> { key1 val1 key2 val2 }
 *       @endcode
 *
 *   - **JSON/Python-like syntax**
 *     - Vectors are enclosed in brackets with comma-separated elements:
 *       @code
 *         <keyword> [ v1 , v2 , ... ]
 *       @endcode
 *     - Maps use `key : value` pairs separated by commas (a trailing comma is
 *       allowed):
 *       @code
 *         <keyword> { key1 : val1 , key2 : val2 }
 *       @endcode
 *
 * For `TRUST_Deriv<U>` (OWN_PTR-like) element types, each individual object must
 * specify its concrete type first, followed by whatever content its `readOn`
 * expects. In maps of such objects, each value is automatically named using its
 * map key. For example, given a `std::map<std::string, TRUST_Deriv<Nom>>`:
 * @code
 *   my_map { nom1 : Nom foo , nom2 : Nom bar , nom3 : Motcle foobar }
 * @endcode
 * will produce a map where `my_map.at("nom1")` is named `"nom1"` and
 * `my_map.at("nom3")` can be `ref_cast` to a `Motcle`.
 *
 * Duplicate map keys are rejected with a fatal error, and for vectors whose
 * expected size has been set via the `size` argument, a mismatch between the
 * declared and the read length is also fatal.
 *
 * For detailed examples, see the unit tests in `unit_params.cpp`.
 */
class Param
{
public:
  /*! @brief Whether a registered parameter must be present in the input or may be omitted. */
  enum Nature { OPTIONAL = 0, REQUIRED = 1 };

  /*!
   * @brief Build a `Param` that will parse the parameters of an object named @p owner_name.
   *
   * @param owner_name Human-readable name of the owning object, used to produce informative
   *                   error messages when parsing fails.
   */
  Param(const char * owner_name);

  // ---------------------------------------------------------------------------
  // Scalar overloads of ajouter(...)
  // ---------------------------------------------------------------------------

  /*!
   * @brief Register an integer parameter.
   *
   * Expected syntax in the input file:
   * @code
   *   <keyword> <integer_value>
   * @endcode
   * For example, if the keyword is `nb_pas_dt_max`, the input would contain:
   * @code
   *   nb_pas_dt_max 100
   * @endcode
   *
   * @param keyword Name of the keyword (optionally with synonyms separated by `|`).
   * @param value   Pointer to the `int` that will be written when the keyword is read.
   *                The pointer must remain valid until the read operation completes.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const int* value, Param::Nature nat = Param::OPTIONAL);

#if INT_is_64_ == 2
  /*!
   * @brief Register a `trustIdType` (64-bit id) parameter.
   *
   * Only available when TRUST is built in mixed 32/64-bit mode (`INT_is_64_ == 2`).
   *
   * Expected syntax:
   * @code
   *   <keyword> <integer_value>
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `trustIdType` that will be written when the keyword is read.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const trustIdType* value, Param::Nature nat = Param::OPTIONAL);
#endif

  /*!
   * @brief Register a floating-point parameter.
   *
   * Expected syntax:
   * @code
   *   <keyword> <double_value>
   * @endcode
   * For example:
   * @code
   *   tmax 1.5e-3
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `double` that will be written when the keyword is read.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const double* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register a string parameter.
   *
   * Expected syntax:
   * @code
   *   <keyword> <single_token>
   * @endcode
   * The value is read as a single whitespace-delimited token.
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::string` that will be written when the keyword is read.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const std::string* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register an `Objet_U` parameter read with its standard `readOn`.
   *
   * The content that follows the keyword is forwarded to the object's `operator>>`
   * (i.e. its `readOn` method). The exact expected syntax therefore depends on the
   * concrete type of the object.
   *
   * Expected syntax:
   * @code
   *   <keyword> <whatever readOn expects for this object>
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `Objet_U` subclass that will be filled when the keyword
   *                is read. Must remain valid until the read completes.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const Objet_U* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register an `ArrOfInt` whose size has already been fixed.
   *
   * The array must have been resized to its expected number of elements @em before
   * the read takes place. Exactly that many integers will be consumed from the stream.
   *
   * Expected syntax (for an array of size N):
   * @code
   *   <keyword> <v_1> <v_2> ... <v_N>
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `ArrOfInt`; its size must be set before calling
   *                `lire_avec_accolades_depuis`.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter_arr_size_predefinie(const char * keyword, const ArrOfInt* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register an `ArrOfDouble` whose size has already been fixed.
   *
   * Same contract as `ajouter_arr_size_predefinie(const char*, const ArrOfInt*, ...)`
   * but for doubles.
   *
   * Expected syntax (for an array of size N):
   * @code
   *   <keyword> <v_1> <v_2> ... <v_N>
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `ArrOfDouble`; its size must be set before calling
   *                `lire_avec_accolades_depuis`.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter_arr_size_predefinie(const char * keyword, const ArrOfDouble* value, Param::Nature nat = Param::OPTIONAL);

  // ---------------------------------------------------------------------------
  // Vector overloads of ajouter(...)
  //
  // See the class-level documentation for the two supported syntaxes.
  // ---------------------------------------------------------------------------

  /*!
   * @brief Register a vector of integers.
   *
   * Expected syntax (either of the two below):
   * @code
   *   <keyword> <N> <v_1> <v_2> ... <v_N>              // historical TRUST syntax
   *   <keyword> [ <v_1> , <v_2> , ... , <v_N> ]        // JSON/Python-like syntax
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::vector<int>` that will be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   * @param size    If non-negative, the expected number of elements. A read whose
   *                length differs from this value is treated as an error. Use `-1`
   *                (the default) to accept any length.
   */
  void ajouter(const char * keyword, const std::vector<int>* value, Param::Nature nat = Param::OPTIONAL, int size = -1);

  /*!
   * @brief Register a vector of doubles.
   *
   * Expected syntax (either of the two below):
   * @code
   *   <keyword> <N> <v_1> <v_2> ... <v_N>              // historical TRUST syntax
   *   <keyword> [ <v_1> , <v_2> , ... , <v_N> ]        // JSON/Python-like syntax
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::vector<double>` that will be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   * @param size    If non-negative, the expected number of elements. Use `-1` (the
   *                default) to accept any length.
   */
  void ajouter(const char * keyword, const std::vector<double>* value, Param::Nature nat = Param::OPTIONAL, int size = -1);

  /*!
   * @brief Register a vector of strings.
   *
   * Each string element is read as a single whitespace-delimited token.
   *
   * Expected syntax (either of the two below):
   * @code
   *   <keyword> <N> <s_1> <s_2> ... <s_N>              // historical TRUST syntax
   *   <keyword> [ <s_1> , <s_2> , ... , <s_N> ]        // JSON/Python-like syntax
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::vector<std::string>` that will be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   * @param size    If non-negative, the expected number of elements. Use `-1` (the
   *                default) to accept any length.
   */
  void ajouter(const char * keyword, const std::vector<std::string>* value, Param::Nature nat = Param::OPTIONAL, int size = -1);

  /*!
   * @brief Register a vector of `OWN_PTR`-like objects (`TRUST_Deriv<T>`).
   *
   * Each element must specify its concrete type first, followed by whatever content
   * its `readOn` expects. At read time, each entry is type-checked: a value whose
   * concrete type is not a subtype of `T` triggers a fatal error.
   *
   * Expected syntax (either of the two below):
   * @code
   *   // Historical TRUST syntax
   *   <keyword> <N> <Type_1> <readOn content_1> <Type_2> <readOn content_2> ...
   *
   *   // JSON/Python-like syntax
   *   <keyword> [ <Type_1> <readOn content_1> , <Type_2> <readOn content_2> , ... ]
   * @endcode
   *
   * @tparam T      Base type of the polymorphic objects stored in the vector; each
   *                read element must be a subtype of `T`.
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::vector<TRUST_Deriv<T>>` that will be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   * @param size    If non-negative, the expected number of elements. Use `-1` (the
   *                default) to accept any length.
   */
  template<typename T>
  void ajouter(const char * keyword, const std::vector<TRUST_Deriv<T>>* value, Param::Nature nat = Param::OPTIONAL, int size = -1);

  // ---------------------------------------------------------------------------
  // Map overloads of ajouter(...)
  //
  // See the class-level documentation for the two supported syntaxes.
  // ---------------------------------------------------------------------------

  /*!
   * @brief Register a map from string keys to integer values.
   *
   * Expected syntax (either of the two below):
   * @code
   *   <keyword> { <key_1> <v_1> <key_2> <v_2> ... }              // historical TRUST syntax
   *   <keyword> { <key_1> : <v_1> , <key_2> : <v_2> , ... }      // JSON/Python-like syntax
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::map<std::string, int>` that will be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const std::map<std::string, int>* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register a map from string keys to double values.
   *
   * Expected syntax (either of the two below):
   * @code
   *   <keyword> { <key_1> <v_1> <key_2> <v_2> ... }              // historical TRUST syntax
   *   <keyword> { <key_1> : <v_1> , <key_2> : <v_2> , ... }      // JSON/Python-like syntax
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::map<std::string, double>` that will be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const std::map<std::string, double>* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register a map from string keys to string values.
   *
   * Each string value is read as a single whitespace-delimited token.
   *
   * Expected syntax (either of the two below):
   * @code
   *   <keyword> { <key_1> <s_1> <key_2> <s_2> ... }              // historical TRUST syntax
   *   <keyword> { <key_1> : <s_1> , <key_2> : <s_2> , ... }      // JSON/Python-like syntax
   * @endcode
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::map<std::string, std::string>` that will be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter(const char * keyword, const std::map<std::string, std::string>* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register a map from string keys to `OWN_PTR`-like objects (`TRUST_Deriv<T>`).
   *
   * Each value must specify its concrete type first, followed by whatever content
   * its `readOn` expects. After reading, each object is named using the key it is
   * registered under (via `nommer`), and each value is type-checked against `T`:
   * a value whose concrete type is not a subtype of `T` triggers a fatal error.
   *
   * Expected syntax (either of the two below):
   * @code
   *   // Historical TRUST syntax
   *   <keyword> { <key_1> <Type_1> <readOn content_1> <key_2> <Type_2> <readOn content_2> ... }
   *
   *   // JSON/Python-like syntax
   *   <keyword> { <key_1> : <Type_1> <readOn content_1> , <key_2> : <Type_2> <readOn content_2> , ... }
   * @endcode
   *
   * See the class-level documentation for a concrete example.
   *
   * @tparam T      Base type of the polymorphic objects stored as map values; each
   *                read value must be a subtype of `T`.
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `std::map<std::string, TRUST_Deriv<T>>` that will
   *                be filled.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  template<typename T>
  void ajouter(const char * keyword, const std::map<std::string, TRUST_Deriv<T>>* value, Param::Nature nat = Param::OPTIONAL);

  // ---------------------------------------------------------------------------
  // Other registrations
  // ---------------------------------------------------------------------------

  /*!
   * @brief Register a boolean flag whose mere presence switches it to `true`.
   *
   * The flag must hold `false` when it is registered (otherwise a fatal error is
   * raised): flags are toggled on only when the keyword is encountered in the
   * input. Flags are always optional.
   *
   * Expected syntax (the keyword appears alone, no associated value):
   * @code
   *   <keyword>
   * @endcode
   *
   * @param keyword Name of the flag.
   * @param value   Pointer to the `bool` that will be set to `true` if the keyword
   *                is present. Must be `false` at the moment of registration.
   */
  void ajouter_flag(const char * keyword, const bool* value);

  /*!
   * @brief Register a nested `Param` block and return a reference to it so it can
   *        be populated in turn.
   *
   * The nested block has its own `{ ... }` scope in the input file.
   *
   * @param keyword Name of the nested block.
   * @param nat     Whether the block is @ref REQUIRED or @ref OPTIONAL (default).
   * @return A reference to the newly-created nested `Param`, on which
   *         `ajouter(...)` can be called to declare its own sub-parameters.
   */
  Param& ajouter_param(const char * keyword, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Register a keyword handled by `Objet_U::lire_motcle_non_standard`.
   *
   * Use this when the content that follows the keyword cannot be described by any
   * of the standard overloads of `ajouter` and requires ad hoc parsing in the
   * object itself.
   *
   * @param keyword Name of the keyword.
   * @param value   Pointer to the `Objet_U` whose `lire_motcle_non_standard` will
   *                be called with the keyword and the input stream.
   * @param nat     Whether the keyword is @ref REQUIRED or @ref OPTIONAL (default).
   */
  void ajouter_non_std(const char * keyword, const Objet_U* value, Param::Nature nat = Param::OPTIONAL);

  /*!
   * @brief Declare a post-read logical condition that must hold on the parameter values.
   *
   * The condition is a formula expression, evaluated after the whole block has been
   * read. Registered parameters contribute two kinds of variables to the expression:
   *   - `is_read_<keyword>`: `1` if the keyword was read, `0` otherwise
   *   - `value_of_<keyword>`: the numerical value, for simple-type parameters
   *     (`int`, `double`, `bool` flag) only
   *
   * If the condition evaluates to false, the read fails and @p message is reported.
   *
   * @param condition Formula string to evaluate.
   * @param message   Human-readable explanation displayed when the condition fails.
   * @param name      Optional identifier for the condition, used by
   *                  `supprimer_condition`. If `nullptr`, a name is generated.
   */
  void ajouter_condition(const char* condition, const char* message, const char* name = 0);

  /*!
   * @brief Remove a previously-registered keyword from this `Param`.
   *
   * @param keyword Name of the keyword to remove. Must have been registered earlier.
   */
  void supprimer(const char * keyword);

  /*!
   * @brief Remove a previously-registered condition from this `Param`.
   *
   * @param name Name of the condition, as passed to `ajouter_condition`.
   */
  void supprimer_condition(const char* name);

  /*!
   * @brief Add an (option name, integer value) entry to the dictionary attached
   *        to a previously registered integer parameter.
   *
   * Once a dictionary has been attached to an integer parameter, the input file
   * must provide one of the option names (instead of a raw integer); the
   * registered `int` receives the associated value.
   *
   * @param option_name Textual option name as it will appear in the input file.
   * @param value       Integer value assigned to the parameter when @p option_name
   *                    is encountered.
   */
  void dictionnaire(const char * option_name, int value);

  /*!
   * @brief Same as `dictionnaire`, but also attaches a nested `Param` block that
   *        will be read when this option is selected.
   *
   * @param option_name Textual option name.
   * @param value       Integer value assigned to the parameter when @p option_name
   *                    is encountered.
   * @return A reference to the nested `Param` associated with this option.
   */
  Param& dictionnaire_param(const char * option_name, int value);

  // --- Aliases kept for ELI (script generation) -----------------------------

  /*! @brief Alias of `ajouter(const char*, const int*, Nature)`. */
  inline void ajouter_int(const char * c, const int* val, Param::Nature nat = Param::OPTIONAL) { ajouter(c,val,nat); }
  /*! @brief Alias of `ajouter(const char*, const double*, Nature)`. */
  inline void ajouter_double(const char * c, const double* val, Param::Nature nat = Param::OPTIONAL) { ajouter(c,val,nat); }
  /*! @brief Alias of `ajouter(const char*, const Objet_U*, Nature)`. */
  inline void ajouter_objet(const char *c, const Objet_U* obj, Param::Nature nat = Param::OPTIONAL) { ajouter(c,obj,nat); }

  // ---------------------------------------------------------------------------
  // Read / print / query
  // ---------------------------------------------------------------------------

  /*!
   * @brief Parse the parameter block `{ ... }` from @p is.
   *
   * Expects an opening brace, then alternating keywords and values (or blocks),
   * then a closing brace. Required keywords that are missing, unknown keywords,
   * and condition failures all trigger a fatal error.
   *
   * @par Example
   * To read `{ a 1 b 2 }`:
   * @code
   *   Param param(que_suis_je());
   *   int a, b;
   *   param.ajouter("a", &a);
   *   param.ajouter("b", &b);
   *   param.lire_avec_accolades_depuis(is);
   * @endcode
   *
   * @param is Input stream positioned just before the opening brace.
   * @return `1` on success; fatal on error.
   */
  int lire_avec_accolades_depuis(Entree& is);

  /*!
   * @brief Read all required keywords in the order they were registered, without
   *        enclosing braces.
   *
   * Only works when every registered parameter is @ref REQUIRED. Optional
   * parameters are not supported in this mode.
   *
   * @param is Input stream.
   * @return `1` on success; fatal on error.
   */
  int lire_sans_accolade(Entree& is);

  /*! @brief Alias of @ref lire_avec_accolades_depuis. */
  inline int lire_avec_accolades(Entree& is) { return lire_avec_accolades_depuis(is); }

  /*!
   * @brief Low-level entry point used by `lire_avec_accolades_depuis`.
   *
   * @param is         Input stream.
   * @param with_acco  Must be `1` (no other mode currently supported).
   * @return `1` on success; fatal on error.
   */
  int read(Entree& is, int with_acco = 1);

  /*!
   * @brief Print the current state of the registered parameters as a `{ ... }` block.
   *
   * Not implemented for every parameter type; falls back to a fatal error for
   * unsupported types.
   *
   * @param s Output stream to print into.
   */
  void print(Sortie& s) const;

  /*! @brief List of keywords that have actually been read from the input. */
  inline const LIST(Nom)& get_list_mots_lus() const { return list_parametre_lu_; }

  /*!
   * @brief Retrieve the value of a simple-type parameter by its `value_of_<keyword>` name.
   *
   * Only works for parameters of simple types (`int`, `double`, `bool` flag).
   *
   * @param mot_lu Name in the form `value_of_<keyword>`.
   * @return The parameter value as a `double`; fatal if the parameter is unknown.
   */
  double get_value(const Nom& mot_lu) const;

  /*!
   * @brief Validate that every required keyword was read and every condition holds.
   *
   * Called automatically by `lire_avec_accolades_depuis` and `lire_sans_accolade`.
   *
   * @return `1` if every check succeeded, `0` otherwise.
   */
  int check();

protected:
  /*! @brief Default-constructed `Param` with no owner; only for derived classes. */
  Param();

  /*!
   * @brief Look up (or create) the `Objet_a_lire` associated with a keyword string.
   *
   * The keyword string may contain synonyms separated by `|`. If any of the spellings
   * is already used, a fatal error is raised.
   *
   * @param keyword Keyword string, possibly containing synonyms.
   * @return A reference to the underlying `Objet_a_lire` to be configured by the caller.
   */
  Objet_a_lire& create_or_get_objet_a_lire(const char * keyword);

private:
  LIST(Objet_a_lire) list_parametre_a_lire_;                 ///< All registered parameters.
  Nom proprietaire_;                                         ///< Name of the owning object, used for error messages.
  LIST(Nom) list_parametre_lu_;                              ///< Keywords actually read from the input.
  LIST(Nom) list_conditions_;                                ///< Registered post-conditions (formulas).
  LIST(Nom) list_message_erreur_conditions_;                 ///< Error messages associated with conditions.
  LIST(Nom) list_nom_conditions_;                            ///< Identifiers of the conditions.
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

  Objet_a_lire& obj_a_lire = create_or_get_objet_a_lire(mot);

  auto natc = nat == Param::REQUIRED ? Objet_a_lire::REQUIRED :  Objet_a_lire::OPTIONAL;
  obj_a_lire.set_nature(natc);
  auto ptr = const_cast<std::map<std::string, TRUST_Deriv<T>>*>(quoi);

  // captured by copy for error msg
  std::string attr_name = mot;
  std::string prop = proprietaire_.getString();

  // lambda that will set the values of objects in the map
  auto map_initializer = [ptr, attr_name, prop](std::map<std::string, DerObjU>& map)
  {
    for (const auto& key_to_objU: map)
      {
        const auto& key = key_to_objU.first;
        const auto& objU = key_to_objU.second;
        if (sub_type(T, objU.valeur()))
          {
            const T& cast_obj = ref_cast(T, objU.valeur());
            (*ptr)[key] = cast_obj;
            // name the object with the map key
            (*ptr)[key]->nommer(key);
          }
        else
          {
            Cerr <<"When reading '" << prop << "'" << finl;
            Cerr <<"In keyword '" << attr_name << "', wrong type at key " <<  key << finl;
            Cerr <<objU.valeur().le_type() << " is not a subtype of " << T::info_obj.name() << finl;
            Process::exit();
          }

      }
  };
  obj_a_lire.set_map_obj_initializer(map_initializer);
}

#endif
