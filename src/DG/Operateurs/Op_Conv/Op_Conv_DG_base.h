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

#ifndef Op_Conv_DG_base_included
#define Op_Conv_DG_base_included

#include <Operateur_Conv.h>
#include <TRUST_Ref.h>
#include <SFichier.h>
#include <Domaine_DG.h>
#include <Domaine_Cl_DG.h>
#include <Champ_base.h>

class Op_Conv_DG_base: public Operateur_Conv_base
{
  Declare_base(Op_Conv_DG_base);
public:

  void completer() override;
  double calculer_dt_stab() const override;
  inline DoubleTab& calculer(const DoubleTab& inco, DoubleTab& resu) const override;

  int impr(Sortie& os) const override;

  void associer_domaine_cl_dis(const Domaine_Cl_dis_base&) override;
  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;
  void associer_vitesse(const Champ_base&) override;

protected:
  OBS_PTR(Domaine_DG) le_dom_dg_;
  OBS_PTR(Domaine_Cl_DG) la_zcl_dg_;
  OBS_PTR(Champ_base) vitesse_;

  mutable SFichier Flux, Flux_moment, Flux_sum;
};

/*! @brief Computes the convection contribution and stores it in resu, then returns resu.
 *
 */
inline DoubleTab& Op_Conv_DG_base::calculer(const DoubleTab& inco, DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(inco, resu);
}

#endif
