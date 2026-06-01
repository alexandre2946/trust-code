/****************************************************************************
* Copyright (c) 2023, CEA
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

#ifndef Type_info_included
#define Type_info_included
class Objet_U;
class Nom;
class Noms;
class Sortie;
class Synonyme_info;


/*! @brief Models type information for Objet_U objects.
 *
 * @sa Objet_U Nom
 */
class Type_info
{
public:
  ~Type_info();
  Type_info(const char* name, Objet_U* (*f)(), int nb_bases=0, const Type_info** bases=0);
  Type_info(const char* name,                  int nb_bases=0, const Type_info** bases=0);

  inline const Nom& name() const
  {
    return *name_;
  };
  inline int      can_cast(const Type_info* p) const;

  int     same(const Type_info* p) const;
  int     same(const Nom&) const;
  int     has_base(const Type_info* p, int direct=0) const;
  int     has_base(const Nom& , int direct=0) const;
  Sortie&    bases(Sortie&) const;
  Objet_U*   instance() const;
  int     instanciable() const;

  // Static methods:
  static Sortie&           hierarchie(Sortie&) ;
  static int            est_un_type(const char*) ;
  static int            les_sous_types(const Nom&, Noms& sous_types);
  static int            les_sous_types(const Type_info&, Noms& sous_types);
  static const Type_info * type_info_from_name(const char * type_name);
  static Objet_U*          instance(const char* typ);

protected:

private:
  Type_info(Type_info&) {}; // Copy constructor disabled
  Type_info&    operator=(Type_info&); // Assignment operator disabled
  void   ajouter_type(const Type_info& type_info);
  static int search_type_info_name(const char *nom, int& index);

  // Possible names (eg: A|B)
  const char* names_ = "rien";
  // Name and its synonym
  mutable Nom * name_ = nullptr; // (eg: A)
  mutable Nom * synonym_name_ = nullptr; // (eg: B)
  // Object synonym:
  Synonyme_info* synonym_= nullptr; // Synonym

  // Number of base classes of this class
  int nb_bases_ = -1;
  // List of Type_info for the base classes of this class
  const Type_info** b= nullptr;
  // Pointer to the static "cree_instance" method of the class
  // (null if the class is not instantiable)
  Objet_U* (*cree_instance)()= nullptr;

  // List of Type_info for classes declared by declare_base/declare_instanciable
  // The list is sorted alphabetically (case-insensitive)
  static const Type_info** les_types;
  // For each type registered in "les_types", if multiple classes share the same name,
  // then types_homonymes != 0.
  static int * types_homonymes;
  // Number of classes registered in "les_types" and "types_homonymes"
  static int nb_classes;
  // Memory size of the "les_types" and "types_homonymes" arrays
  // (arrays resized by blocks)
  static int les_types_memsize;
};

/*! @brief Traverses the type hierarchy. Returns 1 if p points to a subtype of the current type.
 *
 * @param p Pointer to the type to test.
 * @return 1 if p points to a subtype of the current type, 0 otherwise.
 */
inline int Type_info::can_cast(const Type_info* p) const
{
  return ( (same(p)) || (p->has_base(this)) );
}

#endif
