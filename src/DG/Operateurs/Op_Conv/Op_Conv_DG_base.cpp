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

#include <Op_Conv_DG_base.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <EcrFicPartage.h>
#include <Probleme_base.h>


Implemente_base(Op_Conv_DG_base, "Op_Conv_DG_base", Operateur_Conv_base);

Sortie& Op_Conv_DG_base::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Op_Conv_DG_base::readOn(Entree& s) { return s; }

double Op_Conv_DG_base::calculer_dt_stab() const { return 1e8; }

void Op_Conv_DG_base::completer() { Operateur_base::completer(); }

void Op_Conv_DG_base::associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl)
{
  la_zcl_dg_ = ref_cast(Domaine_Cl_DG, zcl);
}

void Op_Conv_DG_base::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base&)
{
  le_dom_dg_ = ref_cast(Domaine_DG, domaine_dis);
  la_zcl_dg_ = ref_cast(Domaine_Cl_DG, zcl);
}

int Op_Conv_DG_base::impr(Sortie& os) const
{
  return 1;
}

void Op_Conv_DG_base::associer_vitesse(const Champ_base& ch)
{
  vitesse_ = ch;
}
