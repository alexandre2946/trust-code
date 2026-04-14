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

#include <Op_NConserv_Coloc_base.h>
#include <Domaine_Coloc.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_P0_base.h>
#include <Conservation_Euler.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
#include <Milieu_composite_Euler.h>
#include <Fluide_reel_base.h>


Implemente_base(Op_NConserv_Coloc_base,"Op_NConserv_Coloc_base",Operateur_NConserv_base);
Implemente_instanciable(Op_NConserv_Coloc_base_Elem,"Op_NConserv_Coloc_base_Elem",Op_NConserv_Coloc_base);
Implemente_instanciable(Op_NConserv_Coloc_base_Vect,"Op_NConserv_Coloc_base_Vect",Op_NConserv_Coloc_base);

Sortie& Op_NConserv_Coloc_base::printOn(Sortie& os) const { return Operateur_NConserv_base::printOn(os); }
Entree& Op_NConserv_Coloc_base::readOn(Entree& is) { Operateur_NConserv_base::readOn(is);  return is; }
Sortie& Op_NConserv_Coloc_base_Elem::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_Coloc_base_Elem::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is);  return is; }
Sortie& Op_NConserv_Coloc_base_Vect::printOn(Sortie& os) const { return Op_NConserv_Coloc_base::printOn(os); }
Entree& Op_NConserv_Coloc_base_Vect::readOn(Entree& is) { Op_NConserv_Coloc_base::readOn(is);  return is; }
void Op_NConserv_Coloc_base::completer()
{
  Operateur_base::completer();
  assert(le_dom_poly_.non_nul());
}

void Op_NConserv_Coloc_base::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& inc)
{
  le_dom_poly_ = ref_cast(Domaine_Coloc, domaine_dis);
  la_zcl_poly_ = ref_cast(Domaine_Cl_Coloc, zcl);
  le_champ_inco = ref_cast(Champ_Inc_base,inc);
}

void Op_NConserv_Coloc_base::associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl)
{
  la_zcl_poly_ = ref_cast(Domaine_Cl_Coloc, zcl);
}




