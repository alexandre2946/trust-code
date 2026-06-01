/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Check_espace_virtuel_H
#define Check_espace_virtuel_H

#include <Comm_Group.h>
#include <TRUSTVect.h>

// Returns 1 if the virtual space of v is up to date, 0 otherwise
int check_espace_virtuel_vect(const DoubleVect& v);
int check_espace_virtuel_vect(const IntVect& v);
void assert_invalide_items_non_calcules(DoubleVect& v, double valeur = 0.);

template <typename _TYPE_>
extern void remplir_items_non_calcules_(TRUSTVect<_TYPE_>& v, _TYPE_ valeur);

/*! @brief Fills the "non-computed items" of the array with an invalid value.
 *
 * These are all items not listed in v.get_md_vector().valeurs().get_blocs_items_to_compute(). (items not computed by default vector operations, generally the virtual items)
 *   It is recommended to apply this method at the end of functions that do not return an up-to-date virtual space (using declare_espace_virtuel_invalide(...))
 *   so as to trigger an error if the virtual space is used.
 *
 */
template<typename _TYPE_>
inline void remplir_items_non_calcules(TRUSTVect<_TYPE_>& v, _TYPE_ valeur = 0)
{
  remplir_items_non_calcules_(v, valeur);
}

/*! @brief In comm_check_enabled() mode, checks if the virtual space of the vector is up to date; if not, calls exit().
 *
 * This test is only performed in comm_check_enabled() mode because it requires communications.
 *
 */
template<typename _TYPE_>
inline void assert_espace_virtuel_vect(const TRUSTVect<_TYPE_>& v)
{
  if (Comm_Group::check_enabled())
    {
      if (! check_espace_virtuel_vect(v))
        {
          Cerr << "Fatal error in assert_espace_virtuel_vect: virtual space of this vector is not up to date." << finl;
          Process::barrier();
          Process::exit();
        }
    }
}

/*! @brief When compiled without NDEBUG or when running with check_enabled flag, fills all "not computed" items with an invalid value (ie, items that are not VECT_REAL_ITEMS, usually, these are the virtual items).
 *
 *   You should call this method whenever you compute some field values but you don't compute or update the virtual space.
 *
 */
template<typename _TYPE_>
inline void declare_espace_virtuel_invalide(TRUSTVect<_TYPE_>& v)
{
  const _TYPE_ valeur = (std::is_same<_TYPE_,double>::value) ? -98765.4321 : -1999999999;
#ifdef NDEBUG
  if (Comm_Group::check_enabled())
    remplir_items_non_calcules(v, valeur);
#else
  remplir_items_non_calcules(v, valeur);
#endif
}

#endif /* Check_espace_virtuel_H */
