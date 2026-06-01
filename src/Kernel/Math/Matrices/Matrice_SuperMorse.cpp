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
#include <Matrice_SuperMorse.h>

/*! @brief Computes "resu += MATRICE * x" and a dot product (this is a building block for the conjugate gradient, see class Solv_GCP)
 *
 *  Return value:
 *    local contribution of this processor to "(MATRICE * x) dot x"
 *    (note: the dot product counts all rows of the matrix;
 *      shared items are not removed!)
 *    (note: this is different from resu dot x!)
 *
 */
double Matrice_SuperMorse::ajouter_mult_vect_et_prodscal(const DoubleVect& x, DoubleVect& resu) const
{
  assert(resu.size() == x.size());

  const int nb_lignes = lignes_non_vides_.size_array();
  const int *tab_lignes = lignes_non_vides_.addr();
  // The first index of tab1_ we care about is the second element of the array
  // (the first element is always 1)
  assert(tab1_[0] == 1);
  const auto *tab1_ptr = tab1_.addr() + 1;
  assert(tab1_.size_array() == nb_lignes + 1);
  const int *tab2_ptr = tab2_.addr();
  const double *coeff_ptr = coeff_.addr();
  const double *x_ptr = x.addr() - 1; // offset by 1 because we index with Fortran indices
  double *resu_ptr = resu.addr() - 1; // same
  int n = 1; // index of the current coefficient in tab2 and coeff

  double prod_scal_local = 0;

  for (int i = 0; i < nb_lignes; i++)
    {
      // Row index in the matrix:
      const int i_ligne = *(tab_lignes++);
      assert(i_ligne >= 1 && i_ligne <= resu.size_array());
      double r = 0.;
      // End index of the coefficients for this row in tab2 and coeff
      const auto n_fin = *(tab1_ptr++);
      for (; n < n_fin; n++)
        {
          const int colonne = *(tab2_ptr++);
          assert(colonne >= 1 && colonne <= x.size_array());
          const double coef = *(coeff_ptr++);
          const double xb = x_ptr[colonne];
          r += coef * xb;
        }
      resu_ptr[i_ligne] += r;
      // note: we only want to compute the dot product between x and the part added to resu
      //  (not with the existing values of resu before the call)
      prod_scal_local += r * x_ptr[i_ligne];
    }
  return prod_scal_local;
}
