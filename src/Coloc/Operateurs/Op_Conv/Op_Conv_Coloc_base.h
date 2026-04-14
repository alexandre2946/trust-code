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

#ifndef Op_Conv_Coloc_base_included
#define Op_Conv_Coloc_base_included

#include <Op_Conv_PolyMAC_base.h>
#include <Domaine_Coloc.h>


class Op_Conv_Coloc_base : public Op_Conv_PolyMAC_base
{

  Declare_instanciable(Op_Conv_Coloc_base) ;
public:
  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;
  void associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl) override;
  void mettre_a_jour(double temps) override {};
  int has_interface_blocs() const override {  return 1; }
  void completer() override;
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;
  virtual void Riemann_solver(DoubleTab& num_flux) const { Process::exit(); };
  virtual void scheme(DoubleTab&, const int&) const {Process::exit();};

  //void preparer_calcul() override {Op_Conv_PolyMAC_base::preparer_calcul();};
  //double calculer_dt_stab() const override;
  //void set_incompressible(const int) override {};
//  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override
//  {
//    Operateur_base::get_noms_champs_postraitables(nom, opt);
//  };
//  void creer_champ(const Motcle& motlu) override {};
  /* interface ajouter_blocs */
  // elle sont utiliser dans le mathode ajouter puis calculer ( qui probablement sert a construire la matrice ?)
  //void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override {}
  //double waves_speed(const Domaine_Coloc& domaine, const int& f) const;
  //const DoubleTab& waves_speed() const {return waves_speed;};
  //std::vector<OWN_PTR(Champ_Fonc_base)> waves_speed; // a utiliser au lieu d appeler la fonction a chaque fois

};


class Op_Conv_Coloc_base_Elem : public Op_Conv_Coloc_base
{
  Declare_instanciable(Op_Conv_Coloc_base_Elem) ;
public:
  void Riemann_solver(DoubleTab& num_flux) const override;
};


class Op_Conv_Coloc_base_Vect : public Op_Conv_Coloc_base
{
  Declare_instanciable(Op_Conv_Coloc_base_Vect) ;
public:
  void Riemann_solver(DoubleTab& num_flux) const override;
};



#endif /*Op_Conv_Coloc_base_included*/

