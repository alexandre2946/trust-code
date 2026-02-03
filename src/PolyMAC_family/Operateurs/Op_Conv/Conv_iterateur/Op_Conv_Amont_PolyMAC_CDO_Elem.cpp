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

#include <Op_Conv_Amont_PolyMAC_CDO_Elem.h>
#include <Champ_Elem_PolyMAC_CDO.h>
#include <Champ_Face_PolyMAC_CDO.h>

Implemente_instanciable_sans_constructeur(Op_Conv_Amont_PolyMAC_CDO_Elem, "Op_Conv_Amont_PolyMAC_CDO_Elem", Op_Conv_PolyMAC_CDO_iterateur_base);

Sortie& Op_Conv_Amont_PolyMAC_CDO_Elem::printOn(Sortie& s) const
{
  return s << que_suis_je();
}

Entree& Op_Conv_Amont_PolyMAC_CDO_Elem::readOn(Entree& s)
{
  return s;
}

void Op_Conv_Amont_PolyMAC_CDO_Elem::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_cl_dis, const Champ_Inc_base& ch_transporte)
{
  const Domaine_PolyMAC_CDO& zvdf = ref_cast(Domaine_PolyMAC_CDO, domaine_dis);
  const Domaine_Cl_PolyMAC_family& zclvdf = ref_cast(Domaine_Cl_PolyMAC_family, domaine_cl_dis);
  const Champ_Elem_PolyMAC_CDO& inco = ref_cast(Champ_Elem_PolyMAC_CDO, ch_transporte);

  iter_->associer(zvdf, zclvdf, *this);

  Eval_Amont_PolyMAC_CDO_Elem& eval_conv = (Eval_Amont_PolyMAC_CDO_Elem&) iter_->evaluateur();
  eval_conv.associer_domaines(zvdf, zclvdf);          // Evaluateur_PolyMAC_CDO::associer
  eval_conv.associer_inconnue(inco);        // Eval_PolyMAC_CDO_Elem::associer_inconnue
}

/*! @brief associe le champ de vitesse a l'evaluateur
 *
 */
void Op_Conv_Amont_PolyMAC_CDO_Elem::associer_vitesse(const Champ_base& ch_vit)
{
  const Champ_Face_PolyMAC_CDO& vit = ref_cast(Champ_Face_PolyMAC_CDO, ch_vit);

  Eval_Amont_PolyMAC_CDO_Elem& eval_conv = (Eval_Amont_PolyMAC_CDO_Elem&) iter_->evaluateur();
  eval_conv.associer(vit);                // Eval_Conv_PolyMAC_CDO::associer
}

const Champ_base& Op_Conv_Amont_PolyMAC_CDO_Elem::vitesse() const
{
  Eval_Amont_PolyMAC_CDO_Elem& eval_conv = (Eval_Amont_PolyMAC_CDO_Elem&) iter_->evaluateur();
  return eval_conv.vitesse();
}

Champ_base& Op_Conv_Amont_PolyMAC_CDO_Elem::vitesse()
{
  Eval_Amont_PolyMAC_CDO_Elem& eval_conv = (Eval_Amont_PolyMAC_CDO_Elem&) iter_->evaluateur();
  return eval_conv.vitesse();
}

Op_Conv_Amont_PolyMAC_CDO_Elem::Op_Conv_Amont_PolyMAC_CDO_Elem() :
  Op_Conv_PolyMAC_CDO_iterateur_base(Iterateur_PolyMAC_CDO_Elem<Eval_Amont_PolyMAC_CDO_Elem>())
{
}
