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

#include <Op_NConserv_negligeable.h>
#include <Op_NConserv_Coloc_base.h>
#include <Milieu_composite_Euler.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_base.h>
#include <Domaine_Coloc.h>
#include <Pb_Euler.h>

Implemente_base(Op_NConserv_Coloc_base, "Op_NConserv_Coloc_base", Operateur_NConserv_base);

Sortie& Op_NConserv_Coloc_base::printOn(Sortie& os) const { return Operateur_NConserv_base::printOn(os); }

Entree& Op_NConserv_Coloc_base::readOn(Entree& is) { return Operateur_NConserv_base::readOn(is); }

void Op_NConserv_Coloc_base::completer()
{
  if (!sub_type(Pb_Euler, equation().probleme()))
    {
      Cerr << "WHAT !! Operator " << que_suis_je() << " is only available for Pb_Euler not " << equation().probleme().que_suis_je() << " !! " << finl;
      Process::exit();
    }

  // Seulement operateur negligeable si mono-phasique
  if (!sub_type(Op_NConserv_negligeable, *this))
    {
      const Milieu_composite_Euler& mil = ref_cast(Milieu_composite_Euler, equation().probleme().milieu());
      if (mil.noms_phases().size() == 1)
        {
          Cerr << "You are simulating a Single-Phase Euler problem. You need to use a negligible non-conservative operator !!!" << finl;
          Cerr << "Please remove the non-conservative operator " << que_suis_je() << " from your your equation " << equation().que_suis_je() << " !!" << finl;
          Process::exit();
        }
    }

  Operateur_base::completer();
  assert(le_dom_coloc_);
}

void Op_NConserv_Coloc_base::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& inc)
{
  le_dom_coloc_ = ref_cast(Domaine_Coloc, domaine_dis);
  le_dcl_coloc_ = ref_cast(Domaine_Cl_Coloc, zcl);
  le_champ_inco = ref_cast(Champ_Inc_base, inc);
}

void Op_NConserv_Coloc_base::associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl)
{
  le_dcl_coloc_ = ref_cast(Domaine_Cl_Coloc, zcl);
}
