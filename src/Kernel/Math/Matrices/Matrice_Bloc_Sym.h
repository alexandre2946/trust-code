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

#ifndef Matrice_Bloc_Sym_included
#define Matrice_Bloc_Sym_included



/*! @brief Symmetric block matrix class. Derives from both Matrice_Bloc and Matrice_Sym.
 *
 *
 *
 * @sa Matrice_Bloc Matrice_Sym
 */

#include <Matrice_Bloc.h>
#include <Matrice_Sym.h>

////////////////////////////////////////////////////////////
//
// CLASS : Matrice_Bloc_Sym
//
////////////////////////////////////////////////////////////

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*     class Matrice_Bloc_Sym: storage of blocks of a block matrix A(N,M)    */
/*                                                                            */
/*                         storage is done row by row                        */
/*                         in a vector of type VECT(Matrice).                */
/*     Uses:                                                                  */
/*              blocs_    = vector of matrices Aij                            */
/*              N_        = 1st dim of blocs_                                 */
/*              M_        = 2nd dim of blocs_                                 */
/*              nb_blocs_ = total number of blocks (= N_ * M_)               */
/*                                                                            */
/*     Matrix form:                                                           */
/*                                                                            */
/*              [A11 A12 ... A1M]                                             */
/*              [... A22 ... A2M]                                             */
/*          A = [... ... ... ...]  blocs_=[A11,A12,...,A1M,A21,A22,...,A2M,   */
/*              [... ... ... ANM]                     ,...,AN1,AN2,...,ANM]   */
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/

class Matrice_Morse_Sym;

class Matrice_Bloc_Sym : public Matrice_Bloc, public Matrice_Sym
{
  Declare_instanciable_sans_constructeur(Matrice_Bloc_Sym);

public:

  // Constructeurs :
  Matrice_Bloc_Sym(int N=0, int M=0);

  // Printing
  Sortie& imprimer(Sortie& s) const override;
  Sortie& imprimer_formatte(Sortie& s) const override;

  //Multiplications
  DoubleVect& ajouter_multvect_(const DoubleVect& x, DoubleVect& y) const override;
  DoubleVect& ajouter_multvectT_(const DoubleVect& x, DoubleVect& y) const override;
  DoubleTab& ajouter_multTab_(const DoubleTab& x, DoubleTab& y) const override;

  //Conversions:
  void BlocSymToMatMorseSym(Matrice_Morse_Sym& mat) const;

  // Sizing
  void dimensionner(int N, int M) override; // sizing of blocs_

  // Block access: returns the block Aij with A(N,M)
  const Matrice& get_bloc(int i, int j) const override; // (0<=i<N , i<=j<M)
  Matrice& get_bloc(int i, int j) override;

  void get_stencil(Stencil& stencil) const override;

  void get_symmetric_stencil(Stencil& stencil) const override;

  void get_stencil_and_coefficients(Stencil& stencil, StencilCoeffs& coefficients) const override;

  void get_symmetric_stencil_and_coefficients(Stencil& stencil, StencilCoeffs& coefficients) const override;

  bool check_symmetric_block_matrix_structure() const;
  void assert_check_symmetric_block_matrix_structure() const;
};

#endif
