/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Champ_Elem_PolyMAC_HFV_included
#define Champ_Elem_PolyMAC_HFV_included

#include <Op_Diff_PolyMAC_HFV_base.h>
#include <Domaine_PolyMAC_HFV.h>
#include <Champ_Elem_PolyMAC_CDO.h>
#include <Operateur.h>

/*! @brief : class Champ_Elem_PolyMAC_HFV
 *
 *  @brief Field corresponding to a scalar unknown (e.g. temperature or pressure).
 *  Degrees of freedom: value at elements + flux at faces.
 */
class Champ_Elem_PolyMAC_HFV : public Champ_Elem_PolyMAC_CDO
{
  Declare_instanciable(Champ_Elem_PolyMAC_HFV);
public :
  Champ_base& affecter_(const Champ_base& ch) override;

  int fixer_nb_valeurs_nodales(int n) override; //values at elements

  const Domaine_PolyMAC_HFV& domaine_PolyMAC_HFV() const;
  int nb_valeurs_nodales() const override;
  virtual void init_auxiliary_variables(); //request additional auxiliary unknowns (face values)
  int reprendre(Entree& fich) override;
};

#endif /* Champ_Elem_PolyMAC_HFV_included */
