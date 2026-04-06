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

#ifndef Temperature_imposee_paroi_rayo_transp_included
#define Temperature_imposee_paroi_rayo_transp_included

#include <Temperature_imposee_paroi.h>
#include <Cond_lim_rayo_milieu_transp.h>

class Temperature_imposee_paroi_rayo_transp: public Cond_lim_rayo_milieu_transp, public Temperature_imposee_paroi
{
  Declare_instanciable(Temperature_imposee_paroi_rayo_transp);
public:

  int initialiser(double temps) override;
  void completer() override;
  void mettre_a_jour(double temps) override { calculer_Teta_i(temps); }
  void calculer_Teta_i(double temps);

  inline bool is_bc_rayo_milieu_transp(Cond_lim_rayo_milieu_transp *& la_cl_rayo) override
  {
    la_cl_rayo = static_cast<Cond_lim_rayo_milieu_transp*>(this);
    return true;
  }
};

#endif /* Temperature_imposee_paroi_rayo_transp_included */
