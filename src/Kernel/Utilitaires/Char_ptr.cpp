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

#include <Char_ptr.h>

#include <Motcle.h>
#include <string>



/*! @brief Default constructor.
 *
 * Creates the string "??"
 *
 */
Char_ptr::Char_ptr()
{
  nom_ = new char[3];
  nom_[0] = '?';
  nom_[1] = '?';
  nom_[2] = 0;
}



/*! @brief Constructs a name from a character string. The string is copied.
 *
 * @param (const char* nom) the character string to use
 */
Char_ptr::Char_ptr(const char* nom)
{
  nom_ = 0;
  operator=(nom);
}

/*! @brief Copy constructor of a name.
 *
 * @param (const Char_ptr& nom) the name to use
 */
Char_ptr::Char_ptr(const Char_ptr& nom)
{
  nom_ = 0;
  operator=(nom);
}

/*! @brief Destructor.
 *
 */
Char_ptr::~Char_ptr()
{
  if(nom_)
    delete[] nom_;
}


/*! @brief Returns the number of characters in the Char_ptr string including the terminating null character.
 *
 *     Example: Char_ptr("hello").longueur() == 6
 *
 */
int Char_ptr::longueur() const
{
  return ((int)strlen(nom_)+1);
}

/*! @brief Copies the string nom.
 *
 * Modified by BM so that nom can point to a sub-part of nom_
 *
 */
Char_ptr& Char_ptr::operator=(const char* const nom)
{
  if (nom_ == nom)
    return *this;
  char *old = nom_;
  const char *n = nom;
  if (!n)
    n = "??";
  nom_ = new char[strlen(n)+1];
  strcpy(nom_, n);
  // Delete the old after copying the new one, in case nom is a part of nom_
  delete [] old;
  return *this;
}

/*! @brief Copies the Char_ptr nom.
 *
 * @param (const Char_ptr& nom) the name to copy
 * @return (Char_ptr&) reference to this, representing the string of Char_ptr nom
 */
Char_ptr& Char_ptr::operator=(const Char_ptr& nom)
{
  operator=(nom.nom_);
  return *this;
}

/*! @brief Returns a pointer to the character string of the name.
 *
 * @return (char*) pointer to the character string of the name
 */
Char_ptr::operator char*() const
{
  return nom_;
}


void Char_ptr::allocate(int n)
{
  if (nom_)
    delete [] nom_;
  nom_=new char[n+1];
  for (int i=0; i<n; i++)
    nom_[i]=' ';
  nom_[n]='\0';
}
