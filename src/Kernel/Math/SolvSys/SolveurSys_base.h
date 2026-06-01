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

#ifndef SolveurSys_base_included
#define SolveurSys_base_included

#include <TRUSTTabs_forward.h>
#include <Objet_U.h>
#include <Nom.h>
class Matrice_Base;

class SolveurSys_base : public Objet_U
{

  Declare_base_sans_constructeur(SolveurSys_base);

public :

  SolveurSys_base();
  // Each solver should return either 0 for an iterative method, or 1 for a direct method:
  virtual int solveur_direct() const = 0;
  // Search for "x" such that "a*x = b"
  virtual int resoudre_systeme(const Matrice_Base& a, const DoubleVect& b, DoubleVect& x) = 0;
  virtual int resoudre_systeme(const Matrice_Base& a, const DoubleVect& b, DoubleVect& x, int niter_max);
  virtual void fixer_schema_temps_limpr(int l) { schema_temps_limpr_=l; }
  void fixer_limpr(int l) { limpr_=l; }
  void set_return_on_error(int ret) { return_on_error_ = ret; }
  // Print solver convergence if impr option set in the solver AND if time scheme gives authorization
  int limpr() const
  {
    if (schema_temps_limpr_ == 0)
      return 0;
    return  limpr_;
  }

  virtual void reinit() { nouvelle_matrice_=true; }
  bool nouvelle_matrice() const { return nouvelle_matrice_; }
  void fixer_nouvelle_matrice(bool i) { nouvelle_matrice_ = i; }

  // By default all solvers accept Matrice_Morse_Sym (override otherwise)
  virtual int supporte_matrice_morse_sym() { return 1; };

  // Call this to know if "b" (right hand side) must have an updated virtual space before calling resoudre(m, b, x)
  // (this flag cannot be set, it is a property of the solver)
  virtual int get_flag_updated_input() const { return 1; }
  // Call this to tell the solver if you want "x" to have an updated virtual space (default is yes)
  void set_flag_updated_result(bool flag) { echange_ev_resu_ = flag; }
  int get_flag_updated_result() const { return echange_ev_resu_; }

  void save_matrice_secmem_conditionnel(const Matrice_Base& la_matrice, const DoubleVect& secmem,  const DoubleVect& solution, int binaire=1);
  // Methods to set/get reuse_preconditioner_
  void set_reuse_preconditioner(bool flag) { reuse_preconditioner_=flag; }
  bool reuse_preconditioner() { return reuse_preconditioner_; }
  const Nom& get_chaine_lue() const { return chaine_lue_ ; };
  inline bool read_matrix() const { return read_matrix_; };
  void set_read_matrix(bool flag) { read_matrix_ = flag; };
  inline int save_matrix() const { return save_matrice_; };
  void set_save_matrix(int flag) { save_matrice_ = flag; };
  inline const Nom& le_nom() const override { return nom_; };
  inline void nommer(const Nom& nom) override { nom_ = nom; };

protected :
  bool nouvelle_matrice_ = true; // Flag indicating whether storage or factorisation must be recomputed

  // Whether to save the matrix (integer because multiple formats are possible: 1: TRUST 2: PETSc format 3: Matrix Market)
  // WARNING: not a flag in the Param::ajouter_flag sense !!!
  int save_matrice_ = 0;

  bool return_on_error_ = false; // Flag indicating whether to call exit() or return -1 if resoudre_ fails

  // For reading/storing solver parameters:
  Nom chaine_lue_;
  void lecture(Entree&);
private:
  Nom nom_;
  bool read_matrix_ = false; // Flag indicating whether a matrix should be read from file
  int limpr_ = 0;            // Flag to enable or disable printing of solver convergence
  int schema_temps_limpr_ = 1; // Authorization printing flag set by the time scheme
  bool echange_ev_resu_ = true;  // User set flag to tell if the solver must do echange_espace_virtuel() on the result.
  bool reuse_preconditioner_ = false; // Flag to reuse previous preconditioner (default false)
};

#endif

