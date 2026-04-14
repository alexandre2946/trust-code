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

#ifndef Interface_Baer_Nunziato_included
#define Interface_Baer_Nunziato_included

#include <Interface_base.h>

class Interface_Baer_Nunziato : public Interface_base
{
  Declare_instanciable(Interface_Baer_Nunziato);
public:
  void sigma_(const SpanD T, const SpanD P, SpanD res, int ncomp = 1, int ind = 0) const override
  {
    Process::exit("Dont call Interface_Baer_Nunziato::sigma_ !! \n");
  }
  void sigma_h_(const SpanD H, const SpanD P, SpanD res, int ncomp = 1, int ind = 0) const override
  {
    Process::exit("Dont call Interface_Baer_Nunziato::sigma_h_ !! \n");
  }
  void mettre_a_jour(double ) override { }
  void set_param(Param& param) override;
  int id_phase_vitesse_inter() const { return id_vitesse_interface_; }
  int id_phase_pression_inter() const { return id_pression_interface_; }

private:
  int id_vitesse_interface_ = -123, id_pression_interface_ = -123;
};

#endif /* Interface_Baer_nunziato_included */
