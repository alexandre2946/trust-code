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

#ifndef Op_NConserv_negligeable_included
#define Op_NConserv_negligeable_included

#include <Operateur_NConserv_base.h>
#include <Operateur_negligeable.h>

class Champ_base;

class Op_NConserv_negligeable: public Operateur_negligeable, public Operateur_NConserv_base
{
  Declare_instanciable(Op_NConserv_negligeable);
public :
  inline void contribuer_au_second_membre(DoubleTab& ) const override { }
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override { }
  inline void associer_domaine_cl_dis(const Domaine_Cl_dis_base&) override { }
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override { }
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override { }

  /* interface {dimensionner,ajouter}_blocs -> ne font rien */
  int  has_interface_blocs() const override { return 1; }
  const Champ_base& vitesse() const;

  void ajouter_flux(const DoubleTab& inconnue, DoubleTab& contribution) const override { }
  void calculer_flux(const DoubleTab& inconnue, DoubleTab& flux) const override;

  void check_multiphase_compatibility() const override { }

  inline void mettre_a_jour(double temps) override
  {
    Operateur_negligeable::mettre_a_jour(temps);
  }

  inline void associer(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch) override
  {
    Operateur_negligeable::associer(z, zcl, ch);
  }

protected :
  OBS_PTR(Champ_base) la_vitesse;
};

#endif /* Op_NConserv_negligeable_included */
