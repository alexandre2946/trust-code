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

#include <Synonyme_info.h>
#include <Noms.h>

// B.Mathieu, 08/2004
//  The initialization process of these static members is very
//  important: they absolutely must be initialized BEFORE
//  the first call to the constructor Synonyme_info::Synonyme_info(...).
//  This constructor is called to initialize the static member info_obj
//  of all Objet_U objects.
//  Risk of "static initialization order fiasco"
//    (see http://www.parashift.com/c++-faq-lite/ctors.html   [10.11])
//  For now it is ok because we initialize with a constant...

// Array of pointers to synonyms registered during the construction
// of Synonyme_info objects. If multiple synonyms share the same name (Synonyme_info::n), then
// only one is registered in Synonyme_info::les_synonymes.
const Synonyme_info** Synonyme_info::les_synonymes=0;
// For each element of the array "les_synonymes", this number is 1 if the
// synonym name is shared with multiple synonyms, 0 otherwise.
// See "ajouter_synonyme"

int Synonyme_info::nb_classes=0;
int Synonyme_info::les_synonymes_memsize=0;


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

// GF to correctly free the memory we need to at least destroy the name
Synonyme_info::~Synonyme_info()
{
  retirer_synonyme(nom());
  if ((nb_classes==0) && (les_synonymes_memsize!=0))
    {
      delete [] les_synonymes;
      les_synonymes_memsize=0;
    }
}

/*! @brief Searches for the synonym named "nom" in the list of registered synonyms using binary search.
 *
 *   Strings are compared converted to uppercase.
 *   Stores in "index" the index of the synonym if found,
 *   otherwise the index of the next synonym after it (in that case,
 *    les_synonymes[index-1]->n < nom < les_synonymes[index]->n )
 *   Returns 1 if the synonym was found, 0 otherwise.
 *
 */
int Synonyme_info::search_synonyme_info_name(const char *nom, int& index)
{
  assert(nom != 0);
  // [imin..imax] is the interval where the searched index lies
  int imin = 0;
  int imax = nb_classes;
  while (imax > imin)
    {
      // We always have milieu < imax
      int milieu = (imin + imax) / 2;
      int comparaison = strcmp_uppercase(nom, les_synonymes[milieu]->n);
      if (comparaison == 0)
        {
          index = milieu;
          return 1;
        }
      if (comparaison < 0)
        {
          // nom < les_synonymes[milieu]
          // the searched index is therefore less than or equal to "milieu"
          imax = milieu;
        }
      else
        {
          // nom > les_synonymes[milieu]
          // the searched index is therefore strictly greater than "milieu"
          imin = milieu + 1;
        }
    }
  index = imax;
  return 0;
}
/*! @brief Constructor from a name and an array of parent classes.
 *
 * @param (const char* nom) the name of the synonym to create
 * @param (int nb_base) the number of parent classes in the bases[] array
 * @param (const Synonyme_info* bases[]) the array specifying the base synonyms (parents) of the synonym to create
 * @throws Exits on error if the name is not defined (null)
 */
Synonyme_info::Synonyme_info(const char* un_nom, const char* org_name) :
  n(un_nom),org(org_name)
{
  if((un_nom == 0)||(org_name==0))
    {
      Cerr << "Synonyme_info::Synonyme_info(const char* nom,Objet_U* (*f)()...)\n";
      Cerr << " Error : nom==0" << finl;
      assert(0);
      Process::exit();
    }
  ajouter_synonyme(*this);
}

/*! @brief Static method called by Synonyme_info constructors to remove a synonym from the list of registered synonyms.
 *
 *   Verifies that the synonym name exists.
 *
 */
void Synonyme_info::retirer_synonyme(const char* nom)
{
  // Search where to find the synonym in the array:
  int index;
  int existe_deja = search_synonyme_info_name(nom, index);
  if (!existe_deja)
    {
      Cerr<<"A synonym is suppressed whereas it doesn't exist !!!!!"<<finl;
      Process::exit();
    }
  else
    {
      // Remove the synonym from the array at index "index":
      for (int i = index; i < nb_classes-1; i++)
        les_synonymes[i] = les_synonymes[i+1];
      nb_classes--;
    }
}

/*! @brief Static method called by Synonyme_info constructors to add a new synonym to the list of registered synonyms.
 *
 *   Verifies that the synonym name does not already exist.
 *
 */
void Synonyme_info::ajouter_synonyme(const Synonyme_info& synonyme_info)
{
  // Check that there is enough room in the array:
  if (les_synonymes_memsize <= nb_classes + 1)
    {
      static const int INCREMENT = 512;
      // Not enough room in the array: resize it.
      les_synonymes_memsize += INCREMENT;
      const Synonyme_info** nouveau = new const Synonyme_info*[les_synonymes_memsize];
      for (int i = 0; i < nb_classes; i++)
        nouveau[i] = les_synonymes[i];
      delete[] les_synonymes;
      les_synonymes = nouveau;
    }

  // Search where to insert the synonym in the array:
  int index;
  int existe_deja=Type_info::est_un_type(synonyme_info.n);
  if (existe_deja)
    {
      Cerr<<" The synonym "<<synonyme_info.n<<" exists as a class which is forbidden !!!!"<<finl;
      Process::exit();
    }
  existe_deja = search_synonyme_info_name(synonyme_info.n, index);
  if (existe_deja)
    {
      Cerr<<" Synonym "<<synonyme_info.n<<" already exists, which is forbidden!!!!"<<finl;
      Process::exit();
    }
  else
    {
      // Add the synonym to the array at index "index":
      for (int i = nb_classes; i > index; i--)
        les_synonymes[i] = les_synonymes[i-1];
      les_synonymes[index] = &synonyme_info;

      nb_classes++;
    }
}

/*! @brief Writes the entire hierarchy of the considered synonym to an output stream.
 *
 * @param (Sortie& os) output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Synonyme_info::hierarchie(Sortie& os)
{
  os << "There is " << nb_classes << " synomyms:" << finl;
  int i= nb_classes;
  while(i--)
    os << les_synonymes[i]->nom() << " <=> "<< les_synonymes[i]->org_name_()<<finl;
  return os << finl<<flush;
}

/*! @brief Tests whether a class synonym exists: if there is a class T whose Synonyme_info has
 *
 *      the name "nom", then est_un_synonyme returns 1,
 *      returns null pointer otherwise.
 *
 * @param (const char* nom) character string associated with a synonym
 * @return (int) return code (0 or 1)
 */
int Synonyme_info::est_un_synonyme(const char* nom)
{
  const Synonyme_info * synonyme = synonyme_info_from_name(nom);
  return (synonyme != 0);
}

/*! @brief Static method that returns a pointer to the Synonyme_info whose name is "synonyme_name".
 *
 * If synonyme_name is not a synonym,
 *   returns a null pointer.
 *
 */
const Synonyme_info * Synonyme_info::synonyme_info_from_name(const char * synonyme_name)
{
  const Synonyme_info * synonyme_info = 0;
  if (synonyme_name != 0)
    {
      int index;
      if (search_synonyme_info_name(synonyme_name, index))
        synonyme_info = les_synonymes[index];
    }
  return synonyme_info;
}
