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

#ifndef Array_tools_included
#define Array_tools_included

#include <TRUSTTab.h>

/*! @brief Utility method to remove duplicates from an array.
 *
 */
template <typename _TYPE_, typename _SIZE_>
void array_trier_retirer_doublons(TRUSTArray<_TYPE_,_SIZE_>& array);

/*! @brief Utility method to compute the intersection of two integer lists.
 *
 */
template <typename _TYPE_, typename _SIZE_>
void array_calculer_intersection(TRUSTArray<_TYPE_,_SIZE_>& liste1, const TRUSTArray<_TYPE_,_SIZE_>& liste2);

/*! @brief Utility method to compute the difference between two sorted integer lists.
 *
 */
void array_retirer_elements(ArrOfInt& sorted_array, const ArrOfInt& sorted_elements_list);

/*! @brief Utility method to search for a value in a sorted array.
 *
 */
int array_bsearch(const ArrOfInt& tab, int valeur);

/*! @brief Lexicographic sort of an array.
 *
 */
template <typename _TYPE_, typename _SIZE_>
int tri_lexicographique_tableau(TRUSTTab<_TYPE_,_SIZE_>& tab);

/*! @brief Indirect sort (sorts the index array which contains row numbers in tab).
 *
 */
template <typename _TYPE_, typename _SIZE_>
int tri_lexicographique_tableau_indirect(const TRUSTTab<_TYPE_,_SIZE_>& tab, ArrOfInt_T<_SIZE_>& index);

/*! @brief Utility method to remove duplicates from an array.
 *
 */
template <typename _SIZE_>
void tableau_trier_retirer_doublons(IntTab_T<_SIZE_>& tab);

#if INT_is_64_ == 2
// BigIntTab = TRUSTTab<int, trustIdType>: doesn't match IntTab_T<_SIZE_> since value/size types differ
void tableau_trier_retirer_doublons(BigIntTab& tab);
#endif

/*! @brief Utility method to find duplicates (allows removing duplicates without changing the order of elements).
 *
 */
void calculer_renum_sans_doublons(const IntTab& tab, ArrOfInt& renum, ArrOfInt& items_a_garder);


template<typename _TYPE_, typename _SIZE_>
inline void append_array_to_array(TRUSTArray<_TYPE_,_SIZE_>& dest, const TRUSTArray<_TYPE_,_SIZE_>& src)
{
  const int n1 = dest.size_array(), n2 = src.size_array();
  dest.resize_array(n1+n2);
  dest.inject_array(src, n2 /* nb elem */, n1 /* dest index */, 0 /* src index */);
}

#endif /* Array_tools_included */
