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
#ifndef Matrice_SuperMorse_included
#define Matrice_SuperMorse_included

#include <Matrice_Base.h>

/*! @brief : Matrix with an even sparser storage than Matrice_Morse: only non-empty rows are stored (saving on the size of tab1_,
 *
 *   tab2_ and coeff_ remain identical)
 *
 */
class Matrice_SuperMorse
{
public :
  const auto& tab1(int i) const { return tab1_[i]; }  // i from 0 to n
  auto& tab1(int i) { return tab1_[i]; }
  const int& tab2(int i) const { return tab2_[i]; } // i from 0 to nnz-1
  int& tab2(int i) { return tab2_[i]; }
  const double& coeff(int i) const { return coeff_[i]; } // i from 0 to nnz-1
  double& coeff(int i) { return coeff_[i]; }

  auto& get_set_tab1() { return tab1_ ; }
  auto& get_set_tab2() { return tab2_ ; }
  auto& get_set_coeff() { return coeff_ ; }

  const auto& get_tab1() const { return tab1_ ; }
  const auto& get_tab2() const { return tab2_ ; }
  const auto& get_coeff() const { return coeff_ ; }

  double ajouter_mult_vect_et_prodscal(const DoubleVect& x, DoubleVect& resu) const;
  // Array containing the indices of non-empty rows (Fortran indices)
  ArrOfInt lignes_non_vides_;
  // tab1_ has size lignes_non_vides_.size_array()+1


protected :
#ifdef TRUST_USE_GPU
  ArrOfTID tab1_;
  BigArrOfInt tab2_;
  BigDoubleVect coeff_;
#else
  IntVect tab1_;
  IntVect tab2_;
  DoubleVect coeff_;
#endif
};

#endif
