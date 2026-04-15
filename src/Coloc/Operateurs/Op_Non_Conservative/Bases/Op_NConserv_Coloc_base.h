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

#ifndef Op_NConserv_Coloc_base_included
#define Op_NConserv_Coloc_base_included

#include <Operateur_NConserv_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Domaine_Coloc.h>
#include <TRUST_Ref.h>
#include <SFichier.h>

class Op_NConserv_Coloc_base : public Operateur_NConserv_base
{
  Declare_base(Op_NConserv_Coloc_base) ;
public:
  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;
  void associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl) override;
  void mettre_a_jour(double temps) override {};
  int has_interface_blocs() const override {  return 1; }
  void completer() override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override {Process::exit();};
  virtual void Abgral_scheme(DoubleTab& num_flux_left, DoubleTab& num_flux_right) const {Process::exit();};
protected:
  OBS_PTR(Domaine_Coloc) le_dom_poly_;
  OBS_PTR(Domaine_Cl_Coloc) la_zcl_poly_;
};

#endif /*Op_NConserv_Coloc_base_included*/

