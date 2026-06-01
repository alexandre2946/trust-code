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

#ifndef Matrice_Sym_included
#define Matrice_Sym_included

#define _DEGRE_POLY_ 5
#define _SEUIL_GCP_ 1e-12

#include <TRUSTTabs_forward.h>
#include <Matrice_Base.h>

/*! @brief Matrice_Sym class - Base class for the representation of symmetric matrices.
 *
 *     This class is in fact an "interface" class that provides access to
 *     methods for solving linear systems with a symmetric matrix:
 *     conjugate gradient (preconditioned), SSOR solver, polynomial preconditioner.
 *     This class contains no data members (other than those inherited from Process)
 *     because it is used via multiple inheritance.
 *     It "accesses" the matrix through the matrix-vector multiplication method
 *     DoubleVect& multvect(const DoubleVect&, DoubleVect& resu) const
 *     which is an abstract method.
 *
 * @sa Matrice_Morse_Sym, This class does not inherit from Objet_U because it is used via multiple inheritance with other classes already inheriting from Objet_U., Abstract class
 */

class Matrice_Sym
{
public :
  virtual ~Matrice_Sym() {};
  Matrice_Sym():est_definie_(0) {}

  int get_est_definie() const;
  void set_est_definie(int);
  void unsymmetrize_stencil(const int nb_lines, const Stencil& symmetric_stencil, Stencil& stencil) const;
  void unsymmetrize_stencil_and_coefficients(const int nb_lines, const Stencil& symmetric_stencil, const StencilCoeffs& symmetric_coefficients, Stencil& stencil, StencilCoeffs& coefficients) const;

protected :
  virtual DoubleTab& ajouter_multTab_(const DoubleTab&, DoubleTab& ) const=0 ;
  virtual DoubleVect& ajouter_multvect_(const DoubleVect&, DoubleVect& ) const=0 ;
  virtual DoubleVect& ajouter_multvectT_(const DoubleVect&, DoubleVect& ) const=0 ;

private :
  int est_definie_;
};

#endif
