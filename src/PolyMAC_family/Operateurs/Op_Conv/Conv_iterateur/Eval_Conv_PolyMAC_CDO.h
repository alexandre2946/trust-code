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

#ifndef Eval_Conv_PolyMAC_CDO_included
#define Eval_Conv_PolyMAC_CDO_included

class Champ_Face_PolyMAC_CDO;
#include <TRUST_Ref.h>
#include <Evaluateur_PolyMAC_CDO.h>

class Champ_Inc_base;

class Eval_Conv_PolyMAC_CDO : public Evaluateur_PolyMAC_CDO
{
public:
  inline Eval_Conv_PolyMAC_CDO();
  inline Eval_Conv_PolyMAC_CDO(const Eval_Conv_PolyMAC_CDO&);
  void associer(const Champ_Face_PolyMAC_CDO& );
  void mettre_a_jour( );
  const Champ_Inc_base& vitesse() const;
  Champ_Inc_base& vitesse();

protected:
  OBS_PTR(Champ_Face_PolyMAC_CDO) vitesse_;
  DoubleTab dt_vitesse;
};

/*! @brief Default constructor.
 *
 */
inline Eval_Conv_PolyMAC_CDO::Eval_Conv_PolyMAC_CDO()
{}

inline Eval_Conv_PolyMAC_CDO::Eval_Conv_PolyMAC_CDO(const Eval_Conv_PolyMAC_CDO& eval)
  :Evaluateur_PolyMAC_CDO(eval), vitesse_(eval.vitesse_)
{
  dt_vitesse.ref(eval.dt_vitesse);
}

#endif
