/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Type_info.h>
#include <Noms.h>
#include <Synonyme_info.h>

// B.Mathieu, 08/2004
//  The initialization of these static members is very important: they MUST be initialized
//  BEFORE the first call to the constructor Type_info::Type_info(...).
//  That constructor is called when initializing the static member info_obj of all Objet_U objects.
//  Risk of "static initialization order fiasco"
//    (see http://www.parashift.com/c++-faq-lite/ctors.html   [10.11])
//  Currently OK because initialization is done with a constant value.

// Array of pointers to types registered during construction of Type_info objects.
// If multiple types share the same name (Type_info::n), only one is registered in Type_info::les_types.
const Type_info** Type_info::les_types=0;
// For each element of the "les_types" array, this value is 1 if the type name
// is shared by multiple types, 0 otherwise. See "ajouter_type".
int * Type_info::types_homonymes=0;

int Type_info::nb_classes=0;
int Type_info::les_types_memsize=0;


// [ABN] I don't dare replacing this with STL for efficiency purpose - here we convert only what's needed for comparison
static inline int strcmp_uppercase(const char *n1, const char *n2)
{
  int i = 0;
  unsigned char c1, c2;
  int delta;
  do
    {
      c1 = (unsigned char) ::toupper(n1[i]);
      c2 = (unsigned char) ::toupper(n2[i]);
      delta = c1 - c2;
      i++;
    }
  while ((delta == 0) && (c1 != 0) && (c2 != 0));
  return delta;
}

// GF: to correctly free memory, at minimum the Nom must be destroyed.
Type_info::~Type_info()
{
  // Find where to remove the type in the array:
  int index;
  int existe_deja = search_type_info_name(name(), index);
  if (existe_deja)
    {
      for (int i = index; i<nb_classes-1; i++)
        {
          les_types[i] = les_types[i+1];
          types_homonymes[i] = types_homonymes[i+1];
        }
      nb_classes--;
    }
  if (name_)
    {
      delete name_;
      name_=0;
      delete synonym_;
      synonym_=0;
      delete synonym_name_;
      synonym_name_=0;
    }
  if ((nb_classes==0)&& (les_types_memsize!=0))
    {
      delete [] les_types;
      delete [] types_homonymes;
      les_types_memsize=0;
    }
}

/*! @brief Searches for the type named "nom" in the list of registered types using binary search.
 *
 * Strings are compared after conversion to uppercase.
 * On return, "index" holds the index of the found type, or the index of the first type after it
 * if not found (i.e. les_types[index-1]->n < nom < les_types[index]->n).
 * Returns 1 if found, 0 otherwise.
 *
 * @param nom The type name to search for.
 * @param index On return, the index of the type or the insertion point.
 * @return 1 if the type was found, 0 otherwise.
 */
int Type_info::search_type_info_name(const char *nom, int& index)
{
  assert(nom != 0);
  // [imin..imax] is the interval where the searched index lies
  int imin = 0;
  int imax = nb_classes;
  while (imax > imin)
    {
      // milieu is always < imax
      int milieu = (imin + imax) / 2;
      int comparaison = strcmp_uppercase(nom, les_types[milieu]->name());
      if (comparaison == 0)
        {
          index = milieu;
          return 1;
        }
      if (comparaison < 0)
        {
          // nom < les_types[milieu]: searched index is <= milieu
          imax = milieu;
        }
      else
        {
          // nom > les_types[milieu]: searched index is strictly > milieu
          imin = milieu + 1;
        }
    }
  index = imax;
  return 0;
}
/*! @brief Constructor from a name and an array of base types.
 *
 * @param un_nom The name of the type to create.
 * @param nb_base Number of base types in the bases array.
 * @param the_bases Array specifying the base (parent) types of the type to create.
 * @throws Exits with an error if the name is null.
 */
Type_info::Type_info(const char* un_nom, int nb_base, const Type_info** the_bases) :
  names_(un_nom),
  name_((Nom*)0),
  synonym_name_((Nom*) 0),
  nb_bases_(nb_base),
  b(the_bases),
  cree_instance(0)
{
  if(un_nom == 0)
    {
      Cerr << "Type_info::Type_info(const char* nom,Objet_U* (*f)()...)\n";
      Cerr << " Error : name == 0" << finl;
      assert(0);
      Process::exit();
    }
  ajouter_type(*this);
}

/*! @brief Constructor from a name, a factory function, and an array of base types.
 *
 * The function is used to create an instance of the appropriate type.
 *
 * @param un_nom The name of the type to create.
 * @param f Factory function that creates an instance of this type.
 * @param nb_base Number of base types in the bases array.
 * @param the_bases Array specifying the base (parent) types of the type to create.
 * @throws Exits with an error if the name is null.
 */
Type_info::Type_info(const char* un_nom,
                     Objet_U* (*f)(),
                     int nb_base,
                     const Type_info** the_bases) :
  names_(un_nom),
  name_((Nom*) 0),
  synonym_name_((Nom*) 0),
  nb_bases_(nb_base),
  b(the_bases),
  cree_instance(f)
{
  if(un_nom == 0)
    {
      Cerr << "Type_info::Type_info(const char* nom,Objet_U* (*f)()...)\n";
      Cerr << " Error : name == 0" << finl;
      assert(0);
      Process::exit();
    }
  ajouter_type(*this);
}

/*! @brief Static method called by Type_info constructors to add a new type to the list of registered types.
 *
 * Verifies that the type name does not already exist.
 *
 * @param type_info The Type_info to register.
 */
void Type_info::ajouter_type(const Type_info& type_info)
{
  // Split type_info.names_ into A if | is found (eg: n=A|B)
  // and add a synonym B
  Nom A("");
  Nom B("");
  int synonym_found=0;
  int i = 0;
  unsigned char c;
  do // Start loop
    {
      c = type_info.names_[i];
      // Find a |
      if (c==124)
        {
          if (synonym_found==1)
            {
              Cerr << "More than 1 synonym found in " << type_info.names_ << finl;
              Cerr << "Not supported yet !" << finl;
              Process::exit();
            }
          else
            {
              synonym_found=1;
              i++;
              c = type_info.names_[i];
              if (c==0)
                {
                  Cerr << "Error in a classname which can't finished by | : " << type_info.names_ << finl;
                  Process::exit();
                }
            }
        }
      // Build the synonym name:
      if (synonym_found)
        B+=c;
      else
        A+=c;
      i++;
    }
  while (c!=0); // End loop

  name_ = new Nom(A);
  if (synonym_found)
    {
      //Commented cause too verbose:
      //Cerr << "Keyword " << A << " has a synonym: " << B << finl;
      synonym_name_ = new Nom(B);
      synonym_ = new Synonyme_info(synonym_name_->getChar(),name_->getChar());
    }
  // Check that there is enough space in the array:
  if (les_types_memsize <= nb_classes + 1)
    {
      static const int INCREMENT = 512;
      // Not enough space in the array: resize it.
      les_types_memsize += INCREMENT;
      const Type_info** nouveau = new const Type_info*[les_types_memsize];
      for (int j = 0; j < nb_classes; j++)
        nouveau[j] = les_types[j];
      delete[] les_types;
      les_types = nouveau;

      int * temp = new int[les_types_memsize];
      for (int j = 0; j < nb_classes; j++)
        temp[j] = types_homonymes[j];
      delete[] types_homonymes;
      types_homonymes = temp;
    }

  // Find where to insert the type in the array:
  int existe_deja=Synonyme_info::est_un_synonyme(type_info.name());
  if (existe_deja)
    {
      Cerr<<" class "<<type_info.name()<<" already exists as a synonym which is forbidden!!!!"<<finl;
      Process::exit();
    }
  int index;
  existe_deja = search_type_info_name(type_info.name(), index);
  if (existe_deja)
    {
      types_homonymes[index] = 1;
      // GF: if we have a homonym and the string_macro_trio macro works,
      // there is a problem except for iterators
      if (strcmp(string_macro_trio("VECT",titi),"VECT"))
        {
          if (strncmp(type_info.name(),"Iterateur_",10))
            {
              Cerr<<" type "<<type_info.name()<<" is in double and it is not allowed!!!!"<<finl;
              Process::exit();
            }
        }
    }
  else
    {
      // Insert the type in the array at position "index":
      for (int j = nb_classes; j > index; j--)
        {
          les_types[j] = les_types[j-1];
          types_homonymes[j] = types_homonymes[j-1];
        }
      les_types[index] = &type_info;
      types_homonymes[index] = 0;
      nb_classes++;
    }
}

/*! @brief Writes the base types of the current type to an output stream.
 *
 * @param os Output stream.
 * @return Reference to the modified output stream.
 */
Sortie& Type_info::bases(Sortie& os) const
{
  int i= nb_bases_;
  while(i--)
    os << b[i]->name() << " ";
  return os << finl;
}

/*! @brief Writes the full hierarchy of the considered type to an output stream.
 *
 * @param (Sortie& os) output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Type_info::hierarchie(Sortie& os)
{
  os << "There is " << nb_classes << " classes:" << finl;
  int i= nb_classes;
  while(i--)
    {
      os << les_types[i]->name() << " inherits from ";
      les_types[i]->bases(os);
    }
  return os << flush;
}

/*! @brief Instantiates an Objet_U of the given type. If a class T whose Type_info has the name typ exists,
 *
 *      instance returns a pointer to a new instance of T.
 *      Returns the null pointer otherwise.
 *
 * @param (const char* typ) string associated with a type
 * @return (Objet_U*) pointer to a new Objet_U of type typ
 */
Objet_U* Type_info::instance(const char* typ)
{
  const Type_info * le_type = type_info_from_name(typ);
  Objet_U * instance;
  if (le_type)
    instance = le_type->instance();
  else
    instance = 0;
  return instance;
}

/*! @brief Creates an instance of the class associated with the type_info.
 *
 */
Objet_U* Type_info::instance() const
{
  if (cree_instance == 0)
    {
      Cerr << "Error in Type_info::instance()\n";
      Cerr << " The type " << name() << " is not instantiable" << finl;
      assert(0);
      Process::exit();
    }
  Objet_U * ainstance = (*cree_instance)();
  return ainstance;
}

/*! @brief Tests whether a class of the given type exists. If a class T whose Type_info has the name nom exists,
 *
 *      est_un_type returns 1, null pointer otherwise.
 *
 * @param (const char* nom) string associated with a type
 * @return (int) return code (0 or 1)
 */
int Type_info::est_un_type(const char* nom)
{
  const Type_info * type = type_info_from_name(nom);
  return (type != 0);
}

/*! @brief Tests whether a type belongs to the base types of the considered type. If direct == 0,
 *
 *      returns 1 if (*p) is among the bases of (*this), 0 otherwise.
 *      If direct != 0, returns 1 if (*p) is among the bases of (*this)
 *      or any direct or indirect parent of (*this), 0 otherwise.
 *
 * @param (const Type_info* p) pointer to the type to search for
 * @param (int direct) 0 to search the entire base hierarchy, non-zero for a direct search
 * @return (int) return code (0 or 1)
 */
int Type_info::has_base(const Type_info* p, int direct) const
{
  // Search for p->name() in b
  // if found return 1
  // else if not direct return 0;
  //      else search in the bases of b
  // B. Mathieu modification: test only on the address of type_info,
  // not on the type name...
  if (p == 0)
    {
      return 0;
    }
  else
    {
      if (p == this)
        {
          return 1;
        }
      else
        {
          for (int i = 0; i < nb_bases_; i++)
            if (b[i] == p)
              return 1;
          if (!direct)
            {
              // Verify ancestors at higher levels
              for (int i = 0; i < nb_bases_; i++)
                if (b[i]->has_base(p, direct)) return 1;
            }
        }
    }
  return 0;
}

/*! @brief Tests whether a type belongs to the base types of the considered type. The type to search for is identified by its name.
 *
 *      If direct == 0, returns 1 if the type named name is among the bases of (*this), 0 otherwise.
 *      If direct != 0, returns 1 if the type named name is among the bases of (*this)
 *      or any direct or indirect parent of (*this), 0 otherwise.
 *
 * @param (const Nom& name) the name of the type to search for
 * @param (int direct) 0 to search the entire base hierarchy, non-zero for a direct search
 * @return (int) return code (0 or 1)
 */
int Type_info::has_base(const Nom& aname, int direct) const
{
  // Search for aname in b
  // if found return 1
  // else if not direct return 0;
  //      else search in the bases of b
  // B. Mathieu modification: compare Type_info addresses only,
  //                          not the type name.

  const Type_info * type = type_info_from_name(aname);
  int resultat = has_base(type, direct);
  return resultat;
}

/*! @brief Comparison on the name of a type. Returns 1 if the name string of the considered type and the given name are identical.
 *
 *      Returns 0 otherwise.
 *
 */
int Type_info::same(const Nom& other_name) const
{
  return strcmp(name(),other_name)==0;
}

/*! @brief Returns 1 if this==p, 0 otherwise.
 *
 */
int Type_info::same(const Type_info* p) const
{
  return (this == p);
}

/*! @brief Returns the names of the subtypes for a given parent type.
 *
 * @param (const Type_info& mere) the parent type to search subtypes of
 * @param (Noms& les_sous_types) the names of the subtypes
 * @return (int) number of subtypes returned
 */
int Type_info::les_sous_types(const Type_info& mere, Noms& les_sous_types)
{
  int compteur=0;
  int i= nb_classes;
  // Modif B. Mathieu: name() no longer returns a static.
  const Nom& nom_mere = mere.name();
  // Cerr << "---------" << (const char*) nom_mere << finl;
  while(i--)
    {
      if( les_types[i]->has_base(nom_mere, 0) )
        if (!les_types[i]->same(nom_mere))
          {
            compteur++;
          }
    }
  // Cerr << compteur << finl ;
  if(compteur==0) return 0;
  les_sous_types.dimensionner(compteur);
  compteur=0;
  i= nb_classes;
  while(i--)
    {
      if( les_types[i]->has_base(nom_mere, 0) )
        if (!les_types[i]->same(nom_mere))
          {
            les_sous_types[compteur++]=les_types[i]->name();
          }
    }
  return compteur;
}

/*! @brief Returns the names of the subtypes for a parent type identified by name.
 *
 * @param (const Nom& type) the name of the parent type to search subtypes of
 * @param (Noms& les_sous_types) the names of the subtypes
 * @return (int) number of subtypes returned
 * @throws Exits with an error if the given name does not correspond to a type known to TRUST
 */
int Type_info::les_sous_types(const Nom& type, Noms& sous_types)
{
  if (!est_un_type(type))
    {
      Cerr << type << "is not a type recognized by TRUST Unitaire " << finl;
      Process::exit();
    }
  int i= nb_classes;
  while(i--)
    {
      const Type_info& mere = *les_types[i];
      if( les_types[i]->same(type) )
        {
          return les_sous_types(mere, sous_types);
        }
    }
  return 0;
}

/*! @brief Static method that returns a pointer to the Type_info whose name is "type_name".
 *
 * If type_name is not a known type, returns a null pointer.
 *
 */

const Type_info * Type_info::type_info_from_name(const char * type_name)
{
  const Type_info * type_info = nullptr;
  if (type_name != 0)
    {
      int index;
      if (search_type_info_name(type_name, index))
        {
          if (types_homonymes[index] == 0)
            {
              type_info = les_types[index];
            }
          else
            {
              // The type is registered but the name corresponds
              // to multiple types...
              Cerr << "const Type_info * Type_info::type_info_from_name(const char * type_name)\n";
              Cerr << " The type " << type_name << " has several homonymous\n";
              Cerr << " We doing as if the type is unknown..." << finl;
            }
        }
    }
  return type_info;
}

/*! @brief Returns 1 if the associated type is instantiable (cree_instance is non-null), 0 otherwise.
 *
 */
int Type_info::instanciable() const
{
  return (cree_instance != 0);
}

