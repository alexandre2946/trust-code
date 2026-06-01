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


#ifndef Simpler_included
#define Simpler_included
#include <Simple.h>

//Description

//Ref. Numerical heat Transfer Vol. 10 P. 209-228
// D. S. Jang - R. Jetli - S. Acharya

// Pressure convergence in the Simple algorithm (see class Simple) is slow
// because the off-diagonal terms are neglected when solving the Poisson equation
// for the pressure correction.
// The goal of the Simpler algorithm is to compute a more accurate pressure estimate
// and then apply the Simple algorithm with it.

// A = (M/dt + C(Uk) + D)
// Bt and -B denote the gradient and divergence operators respectively

// -Pressure estimation
// First evaluate the velocity field UPk by solving (Da = diagonal part of A, E = off-diagonal part):
//        Da[Uk-1]UPk = E[Uk-1]Uk-1 + Sv + Ss                        -> UPk
//
// The momentum equation can then be written as:
//        Da[Uk-1]Uk = Da[Uk-1]UPk - BtPk
// Combining this with the continuity equation (-BUk = 0):
//        (BDa-1Bt)Pk = BUPk                                -> Pk
//
//
// -Applying the Simple algorithm to find the solution (Uk,Pk)
//
// (U*k, Pk) satisfies the following momentum equation:
//        A[Uk-1]U*k = -BtPk + Sv + Ss + (M/dt)Uk-1        -> U*k
//
// p'k is evaluated by solving:
//     (BDa-1Bt)p'k = BU*k                                -> p'k
//
// The velocity correction u'k is deduced by solving:
//        Da[Uk-1] (Uk-U*k) = -Btp'k                        -> U'k
//
// The pressure field Pk is not modified by the correction p'k,
// while the velocity is updated by: Uk = U*k + U'k
//
// The algorithm can be repeated until convergence ||Uk-Uk-1|| < seuil_convergence_implicite_.
// In practice, only one iteration is needed (seuil_convergence_implicite_ = 1e6).

// The algorithm coded in this class (iterer_NS) differs slightly from the above description:
// In the UPk evaluation step:
//   - the pressure gradient is taken into account and the relation becomes:
//         Da[Uk-1]UPk = E[Uk-1]Uk-1 + Sv + Ss - BtPk
//
//   - solving (BDa-1Bt)Pcor = BUPk gives a pressure correction
//     and the pressure field is updated: Pk = Pk-1 + Pcor
//     (Pcor replaces the notation Pk used above)
//
//   - resu is readjusted (BtPk = BtPk-1 + BtPcor) to contain -BtPk

class Simpler : public Simple
{

  Declare_instanciable(Simpler);

public :
  void iterer_NS(Equation_base&, DoubleTab& current, DoubleTab& pression, double, Matrice_Morse&, double, DoubleTrav&,int nb_iter,int& converge, int& ok) override;



protected :


};

int inverser_par_diagonale(const Matrice_Morse& matrice,const DoubleTrav& resu,const DoubleTab& present,DoubleTrav& correction_en_vitesse);

#endif

