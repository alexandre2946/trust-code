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

#ifndef Solv_GCP_included
#define Solv_GCP_included

#include <Matrice_SuperMorse.h>
#include <Matrice_Morse_Sym.h>
#include <solv_iteratif.h>
#include <Precond_base.h>
#include <TRUST_Deriv.h>

class Solv_GCP: public solv_iteratif
{
  Declare_instanciable_sans_constructeur(Solv_GCP);
public:
  Solv_GCP();
  int resoudre_systeme(const Matrice_Base&, const DoubleVect&, DoubleVect&) override;
  int resoudre_systeme(const Matrice_Base&, const DoubleVect&, DoubleVect&, int) override;
  inline const OWN_PTR(Precond_base)& get_precond() const { return le_precond_; }
  inline OWN_PTR(Precond_base)& get_precond() { return le_precond_; }
  inline void set_precond(const OWN_PTR(Precond_base)& pre) { le_precond_ = pre; }
  void reinit() override;
  int supporte_matrice_morse_sym() override
  {
    return !le_precond_ || le_precond_->supporte_matrice_morse_sym();
  }
  // GCP does not need that b has an updated virtual space...
  int get_flag_updated_input() const override { return 0; }

protected:
  void prepare_data(const Matrice_Base& matrice, const DoubleVect& secmem, DoubleVect& solution);
  int resoudre_(const Matrice_Base&, const DoubleVect&, DoubleVect&, int);

  bool optimized_ = false;
  OWN_PTR(Precond_base) le_precond_;
  // Data file parameter: should we apply global diagonal preconditioning?
  // In that case, the matrix is copied and multiplied on the left and right by 1/sqrt(diagonal)
  // and the right-hand side and solution are similarly scaled.
  // (Warning: this changes the inner product metric and thus the interpretation of the threshold):
  //     A * X = B
  // <=> D * A * X = D * B    (with D = 1 / sqrt(diagonal))
  // <=> (D * A * D) * Y = D * B, and X = D * Y
  // Property: the diagonal entries of D * A * D are equal to 1
  bool precond_diag_ = false;
  // A vector with virtual items
  DoubleVect tmp_p_avec_items_virt_;
  // Four vectors without virtual items (one could add virtual spaces to all vectors,
  // but this saves space in the cache).
  DoubleVect resu_;
  DoubleVect residu_;
  // tmp_p_ points to the same memory region as tmp_p_avec_items_virt_
  // (this alias is created because vector operations check that
  //  sizes and parallel structures are identical,
  //  e.g. for the preconditioner...)
  DoubleVect tmp_p_;
  DoubleVect tmp_solution_;
  // Matrix of real-real coefficients
  Matrice_Morse_Sym tmp_mat_;
  // Matrix of real-virtual coefficients (storage without empty rows)
  Matrice_SuperMorse tmp_mat_virt_;
  // Memory block holding temporary vectors and possibly matrices
  // (do not read from it directly: it also contains integers!)
  // Using double to ensure correct memory alignment
  ArrOfDouble tmp_data_block_;
  // Renumbering table between the RHS and the temporary vectors
  // (removal of shared and unused virtual items)
  IntVect renum_;
  // If diagonal preconditioning is to be applied,
  // inverse of the square root of the diagonal coefficients of the matrix
  // Warning: the virtual space of this vector is updated because it is needed
  // for D*A*D
  DoubleVect inv_sqrt_diag_;
  int reinit_ = 0; // 0=> nothing ready, 1=> memory allocated, matrix coeffs to copy, 2=> ok
  int nb_it_max_ = -1;
};

#endif /* Solv_GCP_included */
