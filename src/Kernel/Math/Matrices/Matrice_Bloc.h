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

#ifndef Matrice_Bloc_included
#define Matrice_Bloc_included

#include <TRUSTTabs_forward.h>
#include <Matrice_Base.h>
#include <Matrice.h>
#include <TRUSTLists.h>
#include <vector>
#include <TRUST_Vector.h>

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*     class Matrice_Bloc: storage of the blocks of a block matrix A(N,M)    */
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
/*              [A21 A22 ... A2M]                                             */
/*          A = [... ... ... ...]  blocs_=[A11,A12,...,A1M,A21,A22,...,A2M,   */
/*              [AN1 AN2 ... ANM]                     ,...,AN1,AN2,...,ANM]   */
/*                                                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/

class Matrice_Morse;

class Matrice_Bloc : public Matrice_Base
{
  Declare_instanciable_sans_constructeur(Matrice_Bloc);

public :
  int ordre() const override;
  int nb_lignes() const override;
  int nb_colonnes() const override;

  // Methods for computing r+=Ax, implemented in derived classes
  DoubleVect& ajouter_multvect_( const DoubleVect& x, DoubleVect& r ) const override;
  DoubleVect& ajouter_multvectT_( const DoubleVect& x, DoubleVect& r ) const override;
  DoubleTab& ajouter_multTab_( const DoubleTab& x, DoubleTab& r ) const override;

  // multiplication by a scalar
  void scale( const double x ) override;
  // zero out the matrix values
  void clean() override;

  void get_stencil( Stencil& stencil ) const override;
  void get_stencil_and_coefficients(Stencil& stencil, StencilCoeffs& coefficients) const override;
  void get_stencil_and_coeff_ptrs(Stencil& stencil, std::vector<const double *>& coeff_ptr) const override;


  // Printing
  Sortie& imprimer( Sortie& s ) const override;
  Sortie& imprimer_formatte( Sortie& s ) const override;

  // Sizing
  virtual void dimensionner( int N, int M );

  // Block access
  virtual const Matrice& get_bloc( int i, int j ) const;
  virtual Matrice& get_bloc( int i, int j );

  void build_stencil() override;

public :
  // Constructeurs :
  Matrice_Bloc( int N=0, int M=0 );

  // Access to the blocs_ vector characteristics
  int dim( int d ) const;                // if d=0 => N_   if d=1 => M_
  int nb_bloc_lignes() const;            // returns N_
  int nb_bloc_colonnes(void ) const;           // returns M_

  // Fill from a symmetric Morse matrix
  void remplir(const IntLists& voisins, const DoubleLists& valeurs, const DoubleVect& terme_diag, const int i, const int n);

  // Fill from a Morse matrix
  void remplir(const IntLists& voisins, const DoubleLists& valeurs, const int i, const int n, const int j, const int m);

  // Fill from a symmetric or non-symmetric Morse matrix
  void remplir(const IntLists& voisins, const DoubleLists& valeurs, const DoubleVect& terme_diag, const int i, const int n, const int j, const int m);

  // Conversion to a Matrice_Morse
  void block_to_morse( Matrice_Morse& matrix ) const;
  void block_to_morse_with_ptr( Matrice_Morse& result, std::vector<const double *>& coeffs) const;

  void BlocToMatMorse( Matrice_Morse& matrix ) const;


  Matrice_Bloc& operator *=( double x);

  bool check_block_matrix_structure() const;

  void assert_check_block_matrix_structure() const;

protected :
  VECT(Matrice) blocs_;           // les blocs de la matrices source A
  std::vector<Matrice_Base*> blocs_non_nuls_;     // les blocs non nuls
  int N_;                       // 1ere dim de A
  int M_;                       // 2eme dim de A
  int nb_blocs_;                   // nb total des blocs de A (= N_ * M_)

  ArrOfInt offsets_;
  std::vector<int> line_offsets_;
  std::vector<int> column_offsets_;

  template<typename _TAB_T_, typename _VAL_T_>
  void get_stencil_coeff_templ( Stencil& stencil, _TAB_T_& coeff_sp) const;
};

#endif
