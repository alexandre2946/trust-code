/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef SSOR_included
#define SSOR_included

#include <Precond_base.h>
#include <TRUSTTab.h>

class Matrice_Morse_Sym;
class Matrice_Bloc_Sym;

class SSOR : public Precond_base
{
  Declare_instanciable_sans_constructeur(SSOR);
public:
  SSOR();

  double get_omega() { return omega_; }
  void   set_omega(double x)
  {
    omega_ = x;
    reinit(REINIT_COEFF);
  }

  // SSOR does not use the virtual space of the input vector:
  int get_flag_updated_input() const override { return 0; }

protected:
  int preconditionner_(const Matrice_Base&, const DoubleVect& secmem, DoubleVect& solution) override;
  void prepare_(const Matrice_Base&, const DoubleVect& secmem) override;

  void ssor(const Matrice_Morse_Sym&, DoubleVect&);
  void ssor(const Matrice_Bloc_Sym&, DoubleVect&);

  double omega_;
  int algo_fortran_, avec_assert_;
  // Members initialised by the prepare() method
  // Flags for items to process (same structure as the RHS) (shared items handled if algo_items_communs_)
  IntTab items_a_traiter_;
  // For each part of the vector: are there any shared items?
  int algo_items_communs_;
  // Descriptor of the RHS (used to verify that the algorithm has been initialised)
  MD_Vector md_secmem_;
  int line_size_;
  // Precomputed omega divided by the diagonal coefficient of the matrix
  ArrOfDouble omega_diag_;
};

#endif /* SSOR_included */
