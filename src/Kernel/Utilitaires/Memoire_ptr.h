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

#ifndef Memoire_ptr_included
#define Memoire_ptr_included

#include <Objet_U.h>
#include <assert.h>


/*! @brief Pointer within the TRUST memory for an Objet_U.
 *
 * @sa Objet_U Memoire
 */
class Memoire_ptr
{
public :

  int next;

  Memoire_ptr(Objet_U* ptr=0) ;
  inline int libre() const;
  inline void set(Objet_U* ptr);
  inline Objet_U& obj();
  inline Memoire_ptr& operator=(const Memoire_ptr&);
private :
  Objet_U* o_ptr;
};

/*! @brief Indicates whether the memory pointer is free, i.e. whether it points to a non-null Objet_U.
 *
 * @return (int) 1 if the pointer is free
 */
inline int Memoire_ptr::libre() const
{
  return o_ptr==0;
}

/*! @brief Assigns an Objet_U to a memory pointer.
 *
 * @param (Objet_U* ptr) pointer to an Objet_U
 */
inline void Memoire_ptr::set(Objet_U* ptr)
{
  o_ptr=ptr;
}

/*! @brief Returns a reference to the Objet_U pointed to by the memory pointer.
 *
 * @return (Objet_U&) reference to the pointed-to Objet_U
 */
inline Objet_U& Memoire_ptr::obj()
{
  assert(o_ptr!=0);
  return *o_ptr;
}


/*! @brief Assignment operator between memory pointers. In the case A=B, the Objet_U pointed to by A becomes the Objet_U pointed to by B.
 *
 * @param (const Memoire_ptr& mptr) the memory pointer B
 * @return (Memoire_ptr&) the memory pointer A
 */
inline Memoire_ptr& Memoire_ptr::operator=(const Memoire_ptr& mptr)
{
  o_ptr=mptr.o_ptr;
  next=mptr.next;
  return *this;
}

#endif
