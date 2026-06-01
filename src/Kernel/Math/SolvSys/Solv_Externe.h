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

#ifndef Solv_Externe_included
#define Solv_Externe_included

#include <SolveurSys_base.h>
#include <Solv_tools.h>

class Matrice_Morse;
class Matrice_Morse_Sym;

/*! @brief Common stuff for several external solvers.
 *
 * Note: here we use trustIdType for potentially big identifiers, this maps to PetscInt type in Solv_Petsc class
 * (type equality between the both is checked when creating the solver).
 */
class Solv_Externe : public SolveurSys_base, public Solv_tools
{
  Declare_base_sans_constructeur_ni_destructeur(Solv_Externe);
public:
  Solv_Externe() : SolveurSys_base::SolveurSys_base(),
    matrice_symetrique_(-1)
  {}
  ~Solv_Externe() {}

protected:
  void construit_matrice_morse_intermediaire(const Matrice_Base&, Matrice_Morse& );
  void MorseSymToMorse(const Matrice_Morse_Sym& MS, Matrice_Morse& M);
  void Create_lhs_rhs_onDevice();
  public_for_cuda
  template<typename ExecSpace>
  void Update_lhs_rhs(const DoubleVect& b, DoubleVect& x);
  template<typename ExecSpace>
  void Update_solution(DoubleVect& x);
protected:
  const ArrOfInt& indice_coeff_to_keep(const Matrice_Morse&);

  int matrice_symetrique_;      // Flag for matrix symmetry
  ArrOfDouble lhs_;             // Left-hand side without shared items
  ArrOfDouble rhs_;             // Right-hand side without shared items
private:
  ArrOfInt indice_coeff_to_keep_; // CSR matrix coefficients to keep in the TRUST matrix
};


#endif //TRUST_SOLV_EXTERNE_H
