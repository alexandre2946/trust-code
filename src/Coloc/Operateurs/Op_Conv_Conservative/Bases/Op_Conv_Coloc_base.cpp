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

#include <Op_Conv_Coloc_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_base.h>
#include <Domaine_Coloc.h>
#include <Pb_Euler.h>

Implemente_base(Op_Conv_Coloc_base,"Op_Conv_Coloc_base",Operateur_Conv_base);

Sortie& Op_Conv_Coloc_base::printOn(Sortie& os) const { return Operateur_Conv_base::printOn(os); }
Entree& Op_Conv_Coloc_base::readOn(Entree& is) { Operateur_Conv_base::readOn(is);  return is; }

void Op_Conv_Coloc_base::completer()
{
  if (!sub_type(Pb_Euler, equation().probleme()))
    {
      Cerr << "WHAT !! Operator " << que_suis_je() << " is only available for Pb_Euler not " << equation().probleme().que_suis_je() << " !! " << finl;
      Process::exit();
    }

  Operateur_base::completer();
  assert(le_dom_coloc_);
}

void Op_Conv_Coloc_base::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& inc)
{
  le_dom_coloc_ = ref_cast(Domaine_Coloc, domaine_dis);
  le_dcl_coloc_ = ref_cast(Domaine_Cl_Coloc, zcl);
  le_champ_inco = ref_cast(Champ_Inc_base,inc);
}

void Op_Conv_Coloc_base::associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl)
{
  le_dcl_coloc_ = ref_cast(Domaine_Cl_Coloc, zcl);
}

void Op_Conv_Coloc_base::ajouter_blocs(matrices_t mats, DoubleTab& secmem, const tabs_t& semi_impl) const
{

  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();
  const int n = secmem.line_size();
  const int nb_faces = domaine.nb_faces();
  DoubleTrav num_flux(nb_faces, n);
  Riemann_solver(num_flux);

  for (int f = 0; f < nb_faces; f++)
    for (int i = 0; i < 2; i++)
      {
        const int e = f_e(f, i);
        if (e >= 0 && e < domaine.nb_elem())
          {
            for (int k = 0; k < n; k++)
              secmem(e, k) -= (i ? -1 : 1) * num_flux(f, k) * fs(f);
          }
      }
}
