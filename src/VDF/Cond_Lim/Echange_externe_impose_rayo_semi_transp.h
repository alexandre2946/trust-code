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

#ifndef Echange_externe_impose_rayo_semi_transp_included
#define Echange_externe_impose_rayo_semi_transp_included

#include <Cond_lim_rayo_semi_transp.h>
#include <Echange_externe_impose.h>

/*! @brief classe Echange_externe_impose_rayo_semi_transp cette classe est utilisee pour imposer une temperature de paroi imposee
 *
 *    uniquement pour une discretisation VDF.
 *
 *
 */
class Echange_externe_impose_rayo_semi_transp: public Cond_lim_rayo_semi_transp, public Echange_externe_impose
{
  Declare_instanciable(Echange_externe_impose_rayo_semi_transp);

public :
  const Cond_lim_base& la_cl() const override;
  Champ_front_base& temperature_bord();

  // La temperature de paroi etant directement donnee par le champ_front T_ext, il n'y a rien a calculer ici
  void calculer_temperature_bord(double temps) { }
  int compatible_avec_eqn(const Equation_base&) const override { return 1; }
  void completer() override;

  inline bool is_bc_rayo_semi_transp(Cond_lim_rayo_semi_transp*& la_cl_rayo) override
  {
    la_cl_rayo = static_cast<Cond_lim_rayo_semi_transp*>(this);
    return true;
  }
};

#endif /* Echange_externe_impose_rayo_semi_transp_included */
