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

#include <Array_tools.h>

/*! @brief Sorts the array in ascending order and removes duplicates.
 *
 */
template <typename _TYPE_, typename _SIZE_>
void array_trier_retirer_doublons(TRUSTArray<_TYPE_,_SIZE_>& array)
{
  // IntVect is not handled correctly because we use resize_array() instead of resize().
  assert(typeid(array) != typeid(TRUSTVect<_TYPE_, _SIZE_>));
  const _SIZE_ size = array.size_array();
  if (size == 0)
    return;
  // Sort in ascending order
  array.ordonne_array();

  // Remove duplicates
  auto last = std::unique(array.addr(), array.addr()+size);
  _SIZE_ new_size =  static_cast<_SIZE_>(std::distance(array.addr(), last));
  array.resize_array(new_size);
}


/*! @brief Computes the intersection of the two integer lists liste1 and liste2.
 *
 * The result is stored in liste1.
 *   Both lists must be sorted and without duplicates. liste1 is
 *   sorted on output.
 *
 */
template <typename _TYPE_, typename _SIZE_>
void array_calculer_intersection(TRUSTArray<_TYPE_,_SIZE_>& liste1, const TRUSTArray<_TYPE_,_SIZE_>& liste2)
{
  const _SIZE_ sz1 = liste1.size_array();
  const _SIZE_ sz2 = liste2.size_array();
  _SIZE_ j = 0; // Read pointer in liste2
  _SIZE_ k = 0; // Write pointer
  for (_SIZE_ i = 0; i < sz1; i++)
    {
      // Check that the lists are sorted in ascending order
      assert((i >= sz1-1) || (liste1[i] < liste1[i+1]));
      assert((j >= sz2-1) || (liste2[j] < liste2[j+1]));
      const _TYPE_ valeur_i = liste1[i];
      // Advance in liste2 until we find or exceed liste1[i]
      while (j < sz2 && liste2[j] < valeur_i)
        j++;
      if (j == sz2)
        break; // End of list, we are done
      if (liste2[j] == valeur_i)
        {
          liste1[k] = valeur_i;
          k++; // Keep this value
        }
    }
  liste1.resize_array(k);
}

/*! @brief Removes from "sorted_array" the elements that appear in "sorted_elements".
 *
 * Both arrays must initially be sorted in ascending order.
 *   Example:
 *    Input:  sorted_array=[1,4,9,10,12,18], sorted_elements=[3,5,9,10,18,25]
 *    Output: sorted_array=[1,4,12]
 *
 */
void array_retirer_elements(ArrOfInt& sorted_array, const ArrOfInt& sorted_elements_list)
{
  int i_read;      // Index in sorted_array (reading)
  int i_write = 0; // Index in sorted_array (writing position)
  int j = 0;       // Index in sorted_elements
  const int n = sorted_array.size_array();
  const int m = sorted_elements_list.size_array();
  if (m == 0)
    return;

  int j_value = sorted_elements_list[j];
  for (i_read = 0; i_read < n; i_read++)
    {
      // Array sorted?
      assert(i_read == 0 || sorted_array[i_read] > sorted_array[i_read-1]);
      const int i_value = sorted_array[i_read];

      // Advance in the sorted_elements list until we find or exceed
      // element i_value
      while ((j_value < i_value) && (j < m))
        {
          j++;
          if (j == m)
            break;
          assert(sorted_elements_list[j] > j_value); // Array sorted?
          j_value = sorted_elements_list[j];
        }

      if (j == m || j_value != i_value)
        {
          // i_value does not appear in the sorted_elements array, keep it
          sorted_array[i_write] = i_value;
          i_write++;
        }
    }
  sorted_array.resize_array(i_write);
}

template <typename _SIZE_>
static inline int same_line(const IntTab_T<_SIZE_>& v, _SIZE_ i, _SIZE_ j)
{
  const int ls = v.line_size();
  for (int k = 0; k < ls; k++)
    if (v(i,k) != v(j,k))
      return 0;
  return 1;
}

/*! @brief Lexicographic sort of the array tab (ascending order of the first column; if the first column is identical, ascending order
 *
 *   of the second column, etc.).
 *   The array must not be a distributed array.
 *   Return value: number of columns of the array (product of tab.dimension(i) for i>0).
 *
 */
template <typename _TYPE_, typename _SIZE_>
int tri_lexicographique_tableau(TRUSTTab<_TYPE_,_SIZE_>& tab)
{
  // Check that the array is not a distributed array:
  assert(!tab.get_md_vector());

  const _SIZE_ nb_lignes = tab.dimension(0);
  const int nb_colonnes = tab.line_size();
  if (nb_lignes != 0)
    {
      tab.ensureDataOnHost();
      if (nb_colonnes == 1)
        tab.ordonne_array();
      else if (nb_colonnes == 2)
        {
          using pairs = std::array<_TYPE_, 2>;
          _TYPE_ *ptr = tab.addr();
          pairs* tmp = reinterpret_cast<pairs*>(ptr);
          std::sort(tmp, tmp+nb_lignes);
        }
      else if (nb_colonnes == 3)
        {
          using triplets = std::array<_TYPE_, 3>;
          _TYPE_ *ptr = tab.addr();
          triplets* tmp = reinterpret_cast<triplets*>(ptr);
          std::sort(tmp, tmp+nb_lignes);
        }
      else if (nb_colonnes == 4)
        {
          using quadruplets = std::array<_TYPE_, 4>;
          _TYPE_ *ptr = tab.addr();
          quadruplets* tmp = reinterpret_cast<quadruplets*>(ptr);
          std::sort(tmp, tmp+nb_lignes);
        }
      else
        {
          Cerr << "tri_lexicographique_tableau not supported for TRUST tabs with more than 4 columns" << finl;
          Process::exit();
        }
    }
  return nb_colonnes;
}

/*! @brief Same as tri_lexicographique_tableau but sorts the index array which contains row indices of tab such that tab(index[i], *) is
 *
 *   increasing as i increases. Sorts all indices in index...
 *   If the index array has zero size, one of size tab.dimension_tot(0) is created.
 *   Otherwise it is assumed to already contain row indices into tab.
 *   Return value: number of columns of the array (product of tab.dimension(i) for i>0).
 *
 */
template <typename _TYPE_, typename _SIZE_>
int tri_lexicographique_tableau_indirect(const TRUSTTab<_TYPE_,_SIZE_>& tab, ArrOfInt_T<_SIZE_>& index)
{
  using int_t = _SIZE_;
  // Check that the array is not a distributed array:
  assert(!tab.get_md_vector());

  const int_t dimtab = tab.dimension_tot(0);
  if (index.size_array() == 0 && dimtab > 0)
    {
      index.resize_array(dimtab, RESIZE_OPTIONS::NOCOPY_NOINIT);
      for (int_t i = 0; i < dimtab; i++)
        index[i] = i;
    }

  const int_t nb_lignes = index.size_array();
  const int nb_colonnes = tab.line_size();
  const int nb_dim = tab.nb_dim();
  const double epsilon = Objet_U::precision_geom;

  if (nb_lignes != 0)
    {
      std::sort(index.begin(), index.end(), [&](int_t a, int_t b)
      {
        if(nb_dim == 1)
          return ( tab(a)<tab(b) );
        for (int i = 0; i < nb_colonnes; i++)
          {
            if ( std::fabs(tab(a,i)-tab(b,i)) > epsilon )
              return ( tab(a,i)<tab(b,i) );
          }
        return false;
      });
    }
  return nb_colonnes;
}

/*! @brief Sorts the array tab in lexicographic order and removes duplicates (note: [1,2] is not equal to [2,1]).
 *
 */
template <typename _SIZE_>
void tableau_trier_retirer_doublons(IntTab_T<_SIZE_>& tab)
{
  const _SIZE_ nb_lignes = tab.dimension(0);
  if (nb_lignes == 0) return;
  int nb_colonnes = tab.line_size();

  if (nb_colonnes == 1)
    array_trier_retirer_doublons(tab);
  else
    {
      nb_colonnes = tri_lexicographique_tableau(tab);
      if (nb_colonnes == 2)
        {
          _SIZE_ j = 1; // Array size after removal
          _SIZE_ last_x = tab(0, 0);
          _SIZE_ last_y = tab(0, 1);
          for (_SIZE_ i = 1; i < nb_lignes; i++)
            {
              const _SIZE_ x = tab(i, 0);
              const _SIZE_ y = tab(i, 1);
              if (x != last_x || y != last_y)
                {
                  tab(j, 0) = last_x = x;
                  tab(j, 1) = last_y = y;
                  j++;
                }
            }
          tab.resize_dim0(j);
        }
      else
        {
          _SIZE_ j = 0; // Last retained row
          for (_SIZE_ i = 1; i < nb_lignes; i++)
            {
              // If row i differs from row j, keep it:
              if (!same_line(tab, i, j))
                {
                  j++;
                  for (int k = 0; k < nb_colonnes; k++)
                    tab(j, k) = tab(i, k);
                }
            }
          tab.resize_dim0(j+1);
        }
    }
}

// Explicit instanciations
//template const IntVect_T<int> *fct_qsort_tab_ptr<int>;
template int tri_lexicographique_tableau_indirect(const TRUSTTab<int,int>& tab, ArrOfInt_T<int>& index);
template int tri_lexicographique_tableau_indirect(const TRUSTTab<double,int>& tab, ArrOfInt_T<int>& index);
template int tri_lexicographique_tableau(TRUSTTab<int,int>& tab);
template int tri_lexicographique_tableau(TRUSTTab<double,int>& tab);
template void tableau_trier_retirer_doublons(IntTab_T<int>& tab);
template void array_calculer_intersection(TRUSTArray<int,int>& liste1, const TRUSTArray<int,int>& liste2);
template void array_trier_retirer_doublons(TRUSTArray<int,int>& array);
template void array_trier_retirer_doublons(TRUSTArray<double,int>& array);

#if INT_is_64_ == 2
//template const IntVect_T<trustIdType> *fct_qsort_tab_ptr<trustIdType>;
template int tri_lexicographique_tableau_indirect(const TRUSTTab<int,trustIdType>& tab, ArrOfInt_T<trustIdType>& index);
template int tri_lexicographique_tableau_indirect(const TRUSTTab<trustIdType,trustIdType>& tab, ArrOfInt_T<trustIdType>& index);
template int tri_lexicographique_tableau_indirect(const TRUSTTab<double,trustIdType>& tab, ArrOfInt_T<trustIdType>& index);
template int tri_lexicographique_tableau(TRUSTTab<int,trustIdType>& tab);
template int tri_lexicographique_tableau(TRUSTTab<trustIdType,int>& tab);
template int tri_lexicographique_tableau(TRUSTTab<trustIdType,trustIdType>& tab);
template int tri_lexicographique_tableau(TRUSTTab<double,trustIdType>& tab);
template void tableau_trier_retirer_doublons(IntTab_T<trustIdType>& tab);
template void array_calculer_intersection(TRUSTArray<int,trustIdType>& liste1, const TRUSTArray<int,trustIdType>& liste2);
template void array_trier_retirer_doublons(TRUSTArray<int,trustIdType>& array);
template void array_trier_retirer_doublons(TRUSTArray<trustIdType,trustIdType>& array);
template void array_trier_retirer_doublons(TRUSTArray<double,trustIdType>& array);

// BigIntTab = TRUSTTab<int, trustIdType>: value type int, size type trustIdType - doesn't fit IntTab_T<_SIZE_>
static inline int same_line_big(const BigIntTab& v, trustIdType i, trustIdType j)
{
  const int ls = v.line_size();
  for (int k = 0; k < ls; k++)
    if (v(i,k) != v(j,k))
      return 0;
  return 1;
}

void tableau_trier_retirer_doublons(BigIntTab& tab)
{
  const trustIdType nb_lignes = tab.dimension(0);
  if (nb_lignes == 0) return;
  int nb_colonnes = tab.line_size();

  if (nb_colonnes == 1)
    array_trier_retirer_doublons(tab);
  else
    {
      nb_colonnes = tri_lexicographique_tableau(tab);
      if (nb_colonnes == 2)
        {
          trustIdType j = 1;
          int last_x = tab(0, 0);
          int last_y = tab(0, 1);
          for (trustIdType i = 1; i < nb_lignes; i++)
            {
              const int x = tab(i, 0);
              const int y = tab(i, 1);
              if (x != last_x || y != last_y)
                {
                  tab(j, 0) = last_x = x;
                  tab(j, 1) = last_y = y;
                  j++;
                }
            }
          tab.resize_dim0(j);
        }
      else
        {
          trustIdType j = 0;
          for (trustIdType i = 1; i < nb_lignes; i++)
            {
              if (!same_line_big(tab, i, j))
                {
                  j++;
                  for (int k = 0; k < nb_colonnes; k++)
                    tab(j, k) = tab(i, k);
                }
            }
          tab.resize_dim0(j+1);
        }
    }
}
#endif


/*! @brief Finds identical rows in "tab" by lexicographic sort and initializes the sizes and contents of renum and renum_inverse.
 *
 *    renum has size tab.dimension_tot(0).
 *     renum[i] will contain the index of row i in the reduced sorted array (containing unique rows).
 *    renum_inverse contains, for each row of the reduced sorted array, the smallest index of the
 *     corresponding row in tab.
 *     The reduced sorted array can be constructed by extracting rows tab( renum_inverse[i], ...).
 *
 */
void calculer_renum_sans_doublons(const IntTab& tab, ArrOfInt& renum, ArrOfInt& renum_inverse)
{
  // MODIF ELI LAUCOIN 31/01/2012:
  // Completely rewriting this function

  // index allows traversal of tab in order
  ArrOfInt index;
  tri_lexicographique_tableau_indirect(tab, index);
  const int n = index.size_array();

  // resize renum and renum_index
  renum.resize_array(n, RESIZE_OPTIONS::NOCOPY_NOINIT);
  renum_inverse.resize_array(n, RESIZE_OPTIONS::NOCOPY_NOINIT);

  int count  = -1; // row counter in the reduced array
  int latest = -1; // index in the original array of the last row added to the reduced array

  for (int i=0; i<n; ++i)
    {
      // traverse tab in ascending order given by index.
      // if the current row differs from the last row added to the reduced array,
      // increment count and update latest.
      if ( ( latest < 0 ) || ( !(same_line(tab,index[i],index[latest]))) )
        {
          ++count;
          latest=i;
        }

      // record in renum where the current row is in the reduced array
      renum[index[i]]      = count;

      // add the last processed row to the reduced array
      renum_inverse[count] = index[latest];
    }

  // resize renum_inverse to the size of the reduced array

  renum_inverse.resize_array(count+1);
  // END MODIF ELI LAUCOIN 31/01/2012
}

/*! @brief Searches for "valeur" in the array tab by binary search. The array tab must be sorted in ascending order.
 *
*   If not found, returns -1 (including if tab is empty),
 *   otherwise returns an index i such that tab[i] == valeur
 *   (if the value appears multiple times in the array, the first occurrence
 *   is not necessarily returned).
 *
 *  Used in Trio!
 */
int array_bsearch(const ArrOfInt& tab, int valeur)
{
  // attention: all details matter!
  int i = 0;
  int j = tab.size_array(); // j = end of array + 1 (important)
  while (j > i)
    {
      // The array must be sorted
      assert(j == tab.size_array() || tab[i] <= tab[j]);
      const int milieu = (i + j) / 2;
      const int val = tab[milieu];
      if (val > valeur)
        j = milieu; // take the middle value, not milieu - 1
      else if (val < valeur)
        i = milieu + 1; // take milieu + 1, not milieu
      else
        return milieu;
    }
  // If we reach here, i==j, so either:
  // - j == end of array + 1 and the value was not found
  // - tab[j] was tested and is not equal to valeur
  // In both cases, valeur is not in the array
  return -1;
}
