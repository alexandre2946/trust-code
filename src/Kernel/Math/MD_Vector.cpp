/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <MD_Vector.h>
#include <Memoire.h>

//Implemente_base_sans_constructeur_ni_destructeur(MD_Vector_base,"MD_Vector_base",Objet_U);

/*! @brief Constructs an MD_Vector object by copying an existing object.
 *
 * This is the recommended method for creating an MD_Vector object (other than
 * by copying another MD_Vector).
 */
void MD_Vector::copy(const MD_Vector_base& src)
{
  int num_obj = src.duplique();
  MD_Vector_base * p = dynamic_cast<MD_Vector_base *>(Memoire::Instance().objet_u_ptr(num_obj));
  assert(p!= nullptr);
  ptr_.reset(p);
}

/*! @brief Returns 1 if the structures are identical, 0 otherwise.
 *
 */
int MD_Vector::operator==(const MD_Vector& md) const
{
  // For now, very strict test: the two structures are
  //  identical if and only if the pointers are identical,
  //  (i.e. the second was created by copying the first).
  // If we want to remove the multi-reference system and duplicate
  //  the structures on copy, a full equality test
  //  over the entire structure must be done!
  return ptr_ == md.ptr_;
}

/*! @brief Inverse of ==.
 */
int MD_Vector::operator!=(const MD_Vector& md) const
{
  return !operator==(md);
}
