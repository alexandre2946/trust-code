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

#ifndef Simpler_Base_included
#define Simpler_Base_included

#include <Solveur_non_lineaire.h>
#include <Parametre_implicite.h>
#include <TRUSTTabs_forward.h>
#include <Nom.h>

class Solveur_Implicite_base;
class Schema_Temps_base;
class Equation_base;
class Matrice_Morse;
class Matrice;
class Sources;
class Motcle;

//Description
//
// Base class for implicit equation solvers.
// The equation is discretised as:
// (M/dt + C(U) + D)Xn+1 = Sv + Ss + (M/dt)Xn
// where
// C(U) convection matrix
// D    diffusion matrix
// M    mass matrix
// Sv   volumetric source term
// Ss   surface source term (Neumann condition)
// Xn and Xn+1 denote the unknown at times tn and tn+1 respectively
//
// At each time step (tn) we solve the following system:
//         A[Xk] = b[Xk]
// with
// A[Xk] = (M/dt + C(Uk) + D)
// b[Xk] = Sv + Ss + (M/dt)Xk
// k denotes an iteration index for the iterative resolution of the matrix system (tn fixed)

// The resolution method used to determine X is a fixed-point method:
// Set f(Xk) = A[Xk]Xk - b[Xk] so that the solution Xsol satisfies f(Xsol) = 0
// A first-order Taylor expansion gives:
// Xk - Xk+1 = f(Xk)/f'(Xk)
// Making the hypothesis f' = A (not strictly true for the convection part)
// the system becomes: A[Xk](Xk-Xk+1) = f(Xk) or equivalently
//                     A[Xk]Xk+1 = A[Xk]Xk - (A[Xk]Xk-Ss) + Sv +(M/dt)Xk
//
//
// Method iterer_eqn() for equations other than Navier-Stokes:
// - builds the matrix (matrice) and right-hand side (resu) from the equation (assembler_avec_inertie(...))
// - triggers the matrix system resolution (le_solveur_.resoudre_systeme(...))

// Method iterer_NS() called by iterer_eqn() for the specific case of the Navier-Stokes equation.
// Available algorithms: Simple - Simpler - Piso - Implicite
// See subclasses of Simpler_base for descriptions of the respective algorithms.


class Simpler_Base : public  Solveur_non_lineaire
{
  Declare_base_sans_constructeur(Simpler_Base);

public :

  inline int max_iter_implicite();
  inline int max_iter_implicite() const;

  bool iterer_eqn(Equation_base& equation, const DoubleTab& inconnue, DoubleTab& result, double dt, int numero_iteration, int& ok) override =0;

  virtual void iterer_NS(Equation_base&, DoubleTab& current, DoubleTab& pression, double, Matrice_Morse&, double, DoubleTrav&,int nb_iter,int& converge, int& ok)=0;

  void assembler_matrice_pression_implicite(Equation_base& eqn_NS,const Matrice_Morse& matrice,Matrice& matrice_en_pression_2);

protected :

  int no_qdm_;
  double facsec_diffusion_for_sets_ = -1.;

  Entree& lire(const Motcle&, Entree&) override;
  int get_controle_residu() const override { return controle_residu_; }
};

#endif
