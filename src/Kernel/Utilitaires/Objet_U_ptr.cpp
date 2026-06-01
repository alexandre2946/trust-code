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

#include <Synonyme_info.h>
#include <Objet_U_ptr.h>
#include <Memoire.h>
#include <Nom.h>

Implemente_base_sans_constructeur_ni_destructeur(Objet_U_ptr,"Objet_U_ptr",Objet_U);

Sortie& Objet_U_ptr::printOn(Sortie& os) const
{
  const Objet_U * objet = get_Objet_U_ptr();
  if (objet)
    {
      os << objet->le_type() << finl;
      os << (*objet);
    }
  else os << "vide" << finl;

  return os;
}

Entree& Objet_U_ptr::readOn(Entree& is)
{
  detach(); // Clear the existing object
  static Nom nom_type; // static to avoid creating an Objet_U on every call
  is >> nom_type;
  Objet_U * objet = nullptr;
  if (nom_type != "vide")
    {
      objet = typer(nom_type);
      if (! objet) Process::exit();
    }

  set_Objet_U_ptr(objet);
  if (objet) is >> (*objet); // Read the object

  return is;
}

/*! @brief Destructor.
 *
 * It does not destroy the referenced object.
 *
 */
Objet_U_ptr::~Objet_U_ptr()
{
  cle_ = -2; // Paranoia: invalidate the pointer
}

/*! @brief Builds a null pointer (key set to -1).
 *
 */
Objet_U_ptr::Objet_U_ptr() : cle_(-1), ptr_object_id_(-1) { }

void Objet_U_ptr::detach()
{
  Objet_U * ptr = get_Objet_U_ptr();
  if (ptr) delete ptr;
  set_Objet_U_ptr(nullptr);
}

int Objet_U_ptr::associer_(Objet_U& objet)
{
  Objet_U * ptr = get_Objet_U_ptr_check();
  assert(ptr != nullptr);
  int resu = ptr->associer_(objet);
  return resu;
}


/*! @brief Checks if the pointer is valid.
 *
 * The pointer is valid if cle_==-1
 *      or if la_memoire().objet_u(cle_) has the same object_id_ as
 *         the one registered in ptr_object_id_.
 *    If the pointer is not valid, the program stops.
 *  Returns the address of the pointed object (of type Objet_U).
 *
 */
Objet_U * Objet_U_ptr::get_Objet_U_ptr_check() const
{
  if (cle_ != -1 || ptr_object_id_ != -1)
    {
      Objet_U& objet = Memoire::Instance().objet_u(cle_);
      Objet_U * addr = &objet;
      const int id = objet.get_object_id();
      if (id == ptr_object_id_)
        return addr;
      else
        {
          Cerr << "(PE" << me() << ") ";
          Cerr << "Error in Objet_U_ptr::get_Objet_U_ptr_check() : id != ptr_object_id_\n";
          Cerr << " Pointer type " << le_type();
          Cerr << "\n cle_           = " << cle_;
          Cerr << "\n ptr_object_id_ = " << ptr_object_id_;
          Cerr << "\n id             = " << id;
          std::cerr << "\n &la_memoire().objet_u(cle_) = " << addr << std::endl;
          // If it crashes at this point, it means the referenced object
          // has been destroyed but the reference is still in use.
          exit();
        }
    }
  return 0;
}

/*! @brief Verifies that the object pointed to by ptr is of an acceptable type for the pointer (via get_info_ptr).
 *
 */
int Objet_U_ptr::check_Objet_U_ptr_type(const Objet_U * ptr) const
{
  if (ptr == nullptr) return 1; // The null pointer is valid

  const Objet_U& objet = *ptr;
  // Check that the object is of the correct type:
  // type accepted by the pointer:
  const Type_info& type_info_ptr = get_info_ptr();
  // type of the object:
  const Type_info& type_info_obj = *(objet.get_info());
  // Can we cast type_info_obj to type_info_ptr?
  if (! type_info_ptr.can_cast(&type_info_obj))
    {
      Cerr << "(PE" << me() << ") ";
      Cerr << "Error in Objet_U_ptr::get_Objet_U_ptr_check() : Type error\n";
      Cerr << " Pointer type " << le_type();
      Cerr << "\n cle_           = " << cle_;
      Cerr << "\n ptr_object_id_ = " << ptr_object_id_;
      std::cerr << "\n &memoire.objet_u(cle_) = " << ptr;
      Cerr << "\n Type accepted by the pointer : " << type_info_ptr.name();
      Cerr << "\n Object type in reference : " << type_info_obj.name();
      Process::exit();
    }
  return 1;
}



/*! @brief Updates keys when Objet_U objects have been renumbered.
 *
 * @param (const int* const new_ones) array of new numbering
 * @return (int) the new key of the pointer
 */
int Objet_U_ptr::change_num(const int* const new_ones)
{
  Objet_U::change_num(new_ones);

  // Do not call a function that does "verifie"
  // because the memory is currently being modified.
  if (cle_ > -1)
    cle_ = new_ones[cle_];
  return cle_;
}

/*! @brief Duplicates the Objet_U obj and then changes the pointer to this copy.
 *
 * @param (const Objet_U& obj) reference to the Objet_U to copy
 */
void Objet_U_ptr::recopie(const Objet_U& obj)
{
  int cle = obj.duplique();
  Memoire& memoire = Memoire::Instance();
  Objet_U& objet = memoire.objet_u(cle);
  set_Objet_U_ptr(& objet);
}

/*! @brief Returns a pointer to the associated Objet_U. WARNING: the address may be null (if the pointer is null).
 *
 */
Objet_U * Objet_U_ptr::get_Objet_U_ptr() const
{
  const Objet_U * objet;
  if (cle_ < 0)
    {
      objet = nullptr;
    }
  else
    {
      Memoire& memoire = Memoire::Instance();
      objet = & memoire.objet_u(cle_);
    }
  assert(objet == get_Objet_U_ptr_check());
  return (Objet_U*) objet;
}

/*! @brief Makes *this point to the object *ptr. The address may be null (null pointer).
 *
 */
void Objet_U_ptr::set_Objet_U_ptr(Objet_U * ptr)
{
  if (ptr != nullptr)
    {
      cle_ = ptr->numero();
      ptr_object_id_ = ptr->get_object_id();
    }
  else
    {
      cle_ = -1;
      ptr_object_id_ = -1;
    }
  // It is sufficient to check the type here: if the type is correct here and object_id_ does not change afterwards, then the type remains correct.
  assert(check_Objet_U_ptr_type(ptr));
  assert(get_Objet_U_ptr_check() == ptr);
}

/*! @brief Tries to create an instance of type "type".
 *
 * If type is not a type or type is not instantiable => stop.
 *    If type is not a subtype of the pointer type => return 0.
 *    If ok, returns the address of the created object.
 */
Objet_U * Objet_U_ptr::typer(const char * type)
{
  const Type_info * type_info = Type_info::type_info_from_name(type); // Type_info of the requested type
  const Type_info& type_base = get_info_ptr(); // Base type of the OWN_PTR

  if ( get_Objet_U_ptr()) detach();

  Objet_U * instance = nullptr;

  if (type_info == 0)
    {
      const Synonyme_info* syn_info= Synonyme_info::synonyme_info_from_name(type);

      if (syn_info!=0) return typer(syn_info->org_name_());
      else
        {
          Cerr << "Error in Deriv_::typer_(const char* const type)" << finl << type << " is not a type." << finl;
          Cerr << "Type required : derived from " << type_base.name() << finl << finl;
          Cerr << type << " is not a recognized keyword." << finl << "Check your data set." << finl;
          Nom nompb = type;
          if (nompb.find("TURBULENT") != -1 || nompb.find("TURBULENCE") != -1)
            Cerr << finl << "*** NOTE :: Since TRUST V1.8.0, turbulence models are in TrioCFD and not anymore in TRUST.\nTry using TrioCFD executable or contact TRUST support." << finl;
          Process::exit();
        }
    }
  if (! type_info->instanciable())
    {
      Cerr << "Error in Deriv_::typer_(const char* const type).\n" << type << " is not instanciable." << finl;
      Process::exit();
    }

  if (! type_info->has_base(&type_base))
    Cerr << "Error in Deriv_::typer_(const char* const type).\n " << type << " is not a subtype of " << type_base.name() << finl;
  else
    instance = type_info->instance(); // Creates an instance of the type described in type_info

  set_Objet_U_ptr(instance);
  return instance;
}

