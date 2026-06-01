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

#ifndef Piso_included
#define Piso_included

#include <Navier_Stokes_std.h>
#include <Simpler.h>

//Description

// Ref. Journal of Computational Physics 62 P. 40-65
// R. I. Issa

// A = (M/dt + C(Uk) + D)
// Bt and -B denote the gradient and divergence operators respectively
// delta_x denotes the Laplacian operator
// H denotes the convection+diffusion operator

// The PISO (Pressure Implicit Splitting Operator) algorithm is non-iterative
// and consists of three steps * ** and ***.
// The steps of the algorithm are:

// -PREDICTION STEP
//
// Find the velocity field U* satisfying the momentum equation (with pressure Pn)
//        (rho/dt)*(U*-Un) = H(U*) - BtP + S
// U* does not satisfy the continuity equation

// The prediction step is performed by solving the following system:
//        AU* = -BtPn + Sv + Ss + (M/dt)Un                -> U*


// -FIRST CORRECTION STEP
//
// Find U** and P* satisfying the momentum equation and continuity equation (from U* and Un)
//        (rho/dt)*(U**-Un) = H(U*) - BtP* + S
//        delta_x U** = 0

// To stabilise the system, the diagonal part of convection-diffusion is treated implicitly
//        (rho/dt-Ao)*U** - (rho/dt)*Un = H'(U*) -BtP* + S with H(U) = H'(U) + AoU
// Subtracting the prediction equation:
//        (rho/dt-Ao)*(U**-U*) = -Bt(P*-Pn)        and on the other hand delta_x U** = 0
//
// The first correction step is performed by solving (Da = rho/dt-Ao):
//        (BDa-1Bt)P' = BU*                                -> P' -> P* = Pn + P'
// then   Da[Un]U' = -BtP'                                -> U' -> U** = U* + U'


//-SECOND CORRECTION STEP
//
// Find U*** and P** satisfying the momentum equation and continuity equation (from U** and Un)
//        (rho/dt)*(U***-Un) = H(U**) - BtP** + S
//        delta_x U*** = 0

// To stabilise the system, the diagonal part of convection-diffusion is treated implicitly
//        (rho/dt-Ao)*U*** - (rho/dt)*Un = H'(U**) -BtP** + S
// Subtracting the reformulated first correction equation from the reformulated second correction equation:
//        (rho/dt-Ao)*(U***-U**) = H'(U**-U*) -Bt(P**-P*)        and on the other hand delta_x U*** = 0 and delta_x U** = 0
//
// The second correction step is performed by solving (E = H'):
//        (BDa-1Bt)P'' =  (BDa-1E)U'                        -> P'' -> P** = P* + P''
//         DaU'' = EU' -BtP''                                -> U'' -> U*** = U** + U''

// Note: Additional correction steps can be performed.
// The implementation here continues corrections (compt_max-1 maximum)
// unless the residual increases compared to the previous correction.

class Piso : public Simpler
{
  Declare_instanciable(Piso);
public :
  void iterer_NS(Equation_base&, DoubleTab& current, DoubleTab& pression, double, Matrice_Morse&, double, DoubleTrav&,int nb_iter,int& converge, int& ok) override;
protected :

  int nb_corrections_max_ = 21; // maximum number of corrections to refine the projection
  int avancement_crank_ = 0;   // not pure PISO but rather Crank-Nicolson
  int with_sources_ = 0;   // include source terms in the pressure matrix -> more implicit terms, fewer PISO terms

  Entree& lire(const Motcle&, Entree&) override;

private:
  virtual void iterer_NS_PolyMAC_CDO(Navier_Stokes_std& eqn,DoubleTab& current,DoubleTab& pression, double dt, Matrice_Morse& matrice, int& ok);
  // IBM stuff
  void add_penality_term(Navier_Stokes_std& , DoubleTrav& resu , DoubleTrav& gradP);
  void correct_incr_pressure(Navier_Stokes_std& , DoubleTrav& secmem);
  void correct_gradP(Navier_Stokes_std& , DoubleTrav& gradP);
  void correct_pressure(Navier_Stokes_std& , DoubleTab& , DoubleTab& );
};

//Description
// Ref. G. Fauchet

// The IMPLICITE algorithm is non-iterative and consists of two steps * and **.
// The steps of the algorithm are:

// -PREDICTION STEP
//
// Find the velocity field U* satisfying the momentum equation (with pressure Pn)
// As in the PISO algorithm, the diffusion term is expressed implicitly
// and the convection term semi-implicitly.

//         (U*-Un)/dt = H(U*) - BtPn + Sv + Ss
// U* does not satisfy the continuity equation

// The prediction step is performed by solving the following system:
//        AU* = -BtPn + Sv + Ss + (M/dt)Un                -> U*

// -CORRECTION STEP
//
// The momentum equation at step n+1 can be written:
//         (Un+1-Un)/dt = H(Un+1) - BtPn+1 + Sv + Ss
// Subtracting the momentum equation written for U* and neglecting
// the convection and diffusion terms on U' = Un+1 - U*, we get:
//         (Un+1-U*)/dt = -BtP'        with P' = Pn+1 - Pn
//
// The correction step is performed by solving (M = mass matrix):
//        (BM-1Bt)dt*P' = BU*                                -> dt*P' -> P'
//         Un+1 = U* -dt*BtP'                                -> Un+1
//
// Note: The mass matrix is constant, therefore the system (BM-1Bt) is assembled only once.
class Implicite : public Piso
{
  Declare_instanciable(Implicite);
};

#endif /* Piso_included */
