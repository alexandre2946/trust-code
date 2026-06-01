/****************************************************************************
* Copyright (c) 2024, CEA
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


#ifndef Simple_included
#define Simple_included

#include <TRUSTTabs_forward.h>
#include <Simpler_Base.h>

class Operateur_Grad ;

//Description

//Ref. International Journal For Numerical Methods In Fluids Vol. 12 P. 81-92
// C. T. Shaw

// A = (M/dt + C(Uk) + D)
// Bt and -B denote the gradient and divergence operators respectively

// A solution (U,P) is sought in the form:
// U = U* + u'
// P = P* + p'
// where (U*,P*) satisfies the momentum equation and (u',p') is a correction
// applied to satisfy the continuity equation.

// (U*,P*) satisfies the momentum equation:
//        A[Uk-1]U*k = -BtP* + Sv + Ss + (M/dt)Uk-1                -> U*k
//
// Subtracting this from the full momentum equation for (U,P) gives an equation on u' and p':
//        A[Uk-1]u' = -Btp'

// Keeping only the diagonal part (Da) of A, this becomes:
//        u' = -Da-1Btp'

// The full velocity field U satisfies the continuity equation, giving:
//        -Bu' = BU*
// hence (BDa-1Bt)p' = BU*                                         -> p' then u'

// At the end of the iteration, the solution is:
//        U = U* + beta_u u'
//        P = P* + beta_p p'
// beta_u and beta_p are relaxation coefficients between 0 and 1.
// In practice, relaxation is applied only to the pressure.

// The algorithm can be repeated until convergence ||Uk-Uk-1|| < seuil_convergence_implicite_.
// In practice, only one iteration is needed (seuil_convergence_implicite_ = 1e6).


class Simple : public Simpler_Base
{
  Declare_instanciable_sans_constructeur(Simple);

public :

  Simple();
  bool iterer_eqn(Equation_base& equation, const DoubleTab& inconnue, DoubleTab& result, double dt, int numero_iteration, int& ok) override;
  void iterer_NS(Equation_base&, DoubleTab& current, DoubleTab& pression, double, Matrice_Morse&, double, DoubleTrav&,int nb_iter,int& converge, int& ok) override;
  bool iterer_eqs(LIST(OBS_PTR(Equation_base)) eqs, int compteur, int& ok) override;

  /* memoization of iterer_eqs: public so that the iterated power of TRUST-NK can share it */
  using list_of_eq_ptr_t = std::vector<intptr_t>;  // a list of pointers to equations.
  std::map<list_of_eq_ptr_t, Matrice_Bloc> mbloc;

protected :

  DoubleTab Ustar_old;        // U* = alpha_ U*_new + (1-alpha_)*U*_old   in practice alpha = 1
  double alpha_,beta_;  // beta_ : pressure relaxation coefficient  P = P* + beta_*P'  0<beta_<=1
  int with_d_rho_dt_;

  Entree& lire(const Motcle&, Entree&) override;
  void calculer_correction_en_vitesse(const DoubleTrav& correction_en_pression, DoubleTrav& gradP, DoubleTrav& correction_en_vitesse,const Matrice_Morse& matrice, const Operateur_Grad& gradient );

};

void diviser_par_rho_np1_face(Equation_base& eqn, DoubleTab& tab);

#endif
