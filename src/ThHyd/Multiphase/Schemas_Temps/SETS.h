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

#ifndef SETS_included
#define SETS_included

#include <Interface_blocs.h>
#include <Solv_Petsc.h>
#include <Simpler.h>
#include <utility>
#include <vector>
#include <set>

/*! @brief SETS scheme (semi-implicit + stabilisation steps, TRACE-style).
 *
 * SETS ("Stability-Enhancing Two-Step")
 *
 *  Ref : J. H. MAHAFFY, "A stability-enhancing two-step method for fluid flow calculations," Journal of Computational Physics, 46, 3, 329 (1982).
 *
 * @sa Simpler Piso
 */
class SETS: public Simpler
{
  Declare_instanciable_sans_constructeur(SETS);
public:
  SETS();
  Entree& lire(const Motcle&, Entree&) override; /* keyword "criteres_convergence" */
  int nb_valeurs_temporelles_pression() const override /* number of temporal values for the pressure field */
  {
    return 3; /* same as other variables */
  }

  bool iterer_eqn(Equation_base& equation, const DoubleTab& inconnue, DoubleTab& result, double dt, int numero_iteration, int& ok) override;
  void iterer_NS(Equation_base&, DoubleTab& current, DoubleTab& pression, double, Matrice_Morse&, double, DoubleTrav&, int nb_iter, int& converge, int& ok) override;

  /* block elimination of a linear system */
  // input : order -> groups of (variables, block index) to be successively eliminated: e.g. { { {"vitesse", 0 } }, { {"alpha", 0 }, {"temperature", 0 } } }
  //         inco_p -> name of the principal unknown (e.g.: "pression")
  //         mats, sec : matrices / right-hand sides of the linear system mats.d{incos} = {secs}
  //
  //
  // output : A_p / b_p s.t. d{inco} = A_p.d{inco_p} + b_p
  //          return value -> 1 if elimination succeeded, 0 if singularity encountered
  // constraints: - unknowns in the same block must share a common MD_Vector
  //              - for each block { i_1, i_2 }, the matrix { mats[i_j][i_k] } must be block-diagonal w.r.t. this MD_Vector
  //              - outside this diagonal, unknowns in a block may only depend on preceding blocks and on inco_p
  static int eliminer(const std::vector<std::set<std::pair<std::string, int>>> ordre, const std::string inco_p, const std::map<std::string, matrices_t>& mats, const ptabs_t& sec,
                      std::map<std::string, Matrice_Morse>& A_p, tabs_t& b_p);

  /* assembly of a system in inco_p from the expressions d.{inco} : A_p[inco].d{inco_p} + b_p[inco] */
  // input : - inco_p -> the principal unknown
  //         - A_p, b_p -> expressions for the other unknowns computed by eliminer()
  //         - mats, sec -> linear system containing an equation on inco_p (right-hand side in sec[inco_p], Jacobian in mats[inco_p])
  //
  // output : system matrice.d{inco_p} = secmem
  //
  // constraints : all other unknowns must be expressed in A_p / b_p
  static void assembler(const std::string inco_p, const std::map<std::string, Matrice_Morse>& A_p, const tabs_t& b_p, const std::map<std::string, matrices_t>& mats, const ptabs_t& sec,
                        Matrice_Morse& matrice, DoubleTab& secmem, int p_degen);

  double get_default_growth_factor() const override /* time step growth factor */
  {
    return 1.2; /* in case of a failed time step, we recover slowly */
  }

  int iteration_ = 0;  // current iteration number (for evanescence operators)
  int p_degen_ = -1;    // 1 if pressure is degenerate (incompressible medium + no imposed pressure BCs)
  int sets_;      // 1 if the velocity prediction step is performed

  double unknown_positivation(const DoubleTab& uk, DoubleTab& incr); // brings to 0 unknowns that should stay positive

#ifdef PETSCKSP_H
  /* context for the convergence test */
  struct cv_test_t
  {
    void *defctx; // context of the standard convergence test
    SETS *obj;    // the object
    double eps_alpha; // convergence criterion in alpha
    Vec t, v; // Petsc vectors
  };
  DoubleVect norm, residu; // each line equals norm * sum alpha, storage for the residual
  ArrOfTID ix; // indices to retrieve the residual
  cv_test_t *cv_ctx = nullptr;
  void init_cv_ctx(const DoubleTab& secmem, const DoubleVect& norm);
#if PETSC_VERSION_GE(3,24,0)
  static PetscErrorCode destroy_cvctx(void **mctx);
#else
  static PetscErrorCode destroy_cvctx(void *mctx);
#endif
  static PetscErrorCode convergence_test(KSP ksp, PetscInt it, PetscReal rnorm, KSPConvergedReason *reason,void *mctx);
#endif

  double facsec_diffusion_for_sets() const
  {
    return facsec_diffusion_for_sets_;
  }

protected:

  int iter_min_ = 1, iter_max_ = 10; // min/max number of iterations for the nonlinear step
  int first_call_ = 1; // at the very first call, P can be very poor -> velocity is not predicted in SETS
  int pressure_reduction_ = 1; // do we perform pressure reduction?

  /* convergence criteria per unknown (in Linf norm), modifiable via the "criteres_convergence" keyword */
  std::map<std::string, double> crit_conv_;

  /* matrices of the semi-implicit solve */
  std::map<std::string, matrices_t> mats_; // matrices: mats[equation unknown name][other unknown name] = matrix
  Matrice_Bloc mat_semi_impl_; // storage for the mats matrices
  MD_Vector mdv_semi_impl_;    // associated MD_Vector
  std::map<std::string, Matrice_Morse> mat_pred_; // prediction matrices

  /* pressure elimination: dv = A_p[variable name]. dp + b_p[variable name] */
  std::map<std::string, Matrice_Morse> A_p_;

  /* pressure matrix */
  Matrice_Morse matrice_pression_;

  /* Newton out file */
  bool header_written_ = false;
  std::vector<double> incr_var_convergence_;
};

/*! @brief ICE scheme (semi-implicit ICE, CATHARE 3D style).
 *
 * @sa SETS
 */
class ICE: public SETS
{
  Declare_instanciable(ICE);
public:
  bool est_compatible_avec_th_mono() const override /* is this solver compatible with a monolithic thermal solve? */
  {
    return 0; /* no: ICE is explicit in the thermal part */
  }
  double get_default_facsec_max() const override /* recommended facsec_max */
  {
    return 1; /* SETS is semi-implicit */
  }
};

#endif /* SETS_included */
