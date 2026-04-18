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

#include <Op_Conv_Rusanov_Coloc_Vect.h>
#include <Milieu_composite_Euler.h>
#include <Conservation_Euler.h>
#include <Champ_Inc_P0_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Fluide_reel_base.h>
#include <Momentum_Euler.h>
#include <Domaine_Coloc.h>
#include <Pb_Euler.h>

Implemente_instanciable(Op_Conv_Rusanov_Coloc_Vect,"Op_Conv_Rusanov_Coloc_Vect",Op_Conv_Coloc_Vect_base);

Sortie& Op_Conv_Rusanov_Coloc_Vect::printOn(Sortie& os) const { return Op_Conv_Coloc_base::printOn(os); }
Entree& Op_Conv_Rusanov_Coloc_Vect::readOn(Entree& is) { Op_Conv_Coloc_base::readOn(is); return is;}

inline void Op_Conv_Rusanov_Coloc_Vect::scheme(DoubleTab& num_flux, const int f, DoubleTab& flux_l, DoubleTab& flux_r) const
{
  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const DoubleTab& w = le_champ_inco->valeurs();
  const IntTab& f_e = domaine.face_voisins();
  const Momentum_Euler& eq = ref_cast(Momentum_Euler, equation());
  const DoubleTab& vit_n = eq.vitesse_normale();
  const DoubleTab& c = eq.vitesse_son();
  const int nb_phase = ref_cast(Pb_Euler,eq.probleme()).nb_phases();
  const int el = f_e(f, 0), er = f_e(f, 1);
  eq.flux(f, 0, flux_l);
  eq.flux(f, 1, flux_r);

  for (int n = 0; n < nb_phase; n++)
    {
      double un_l = vit_n(f, n), un_r = vit_n(f, n + nb_phase), c_l = c(el, n), c_r = c(er, n);
      double s = std::max(fabs(un_l - c_l), fabs(un_l + c_l));
      s = std::max(s, std::max(fabs(un_r - c_r), fabs(un_r + c_r)));
      for (int d = 0; d < dimension; d++)
        num_flux(f, n + nb_phase * d) = 0.5 * (flux_l(n + nb_phase * d) + flux_r(n + nb_phase * d)) - s * 0.5 * (w(er, n + nb_phase * d) - w(el, n + nb_phase * d));
    }
}
