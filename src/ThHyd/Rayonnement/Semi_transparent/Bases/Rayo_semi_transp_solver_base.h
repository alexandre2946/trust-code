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

#ifndef Rayo_semi_transp_solver_base_included
#define Rayo_semi_transp_solver_base_included

#include <TRUST_Ref.h>
#include <Objet_U.h>

class Eq_rayo_semi_transp;
class Operateur_Diff;
class Matrice_Morse;

class Rayo_semi_transp_solver_base: public Objet_U
{
  Declare_base(Rayo_semi_transp_solver_base);
public:
  void associer_equation_rayo(const Eq_rayo_semi_transp& );

  virtual int nb_colonnes_tot()=0;
  virtual int nb_colonnes()=0;

  virtual void modifier_matrice()=0;
  virtual void assembler_matrice()=0;

  virtual void resoudre(double temps)=0;
  virtual void evaluer_cl_rayonnement(double temps)=0;

protected:
  OBS_PTR(Eq_rayo_semi_transp) eq_rayo_semi_transp_;
};

#endif /* Rayo_semi_transp_solver_base_included */

