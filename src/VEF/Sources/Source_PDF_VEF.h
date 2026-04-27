/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Source_PDF_VEF_included
#define Source_PDF_VEF_included


#include <Navier_Stokes_std.h>
#include <Source_PDF_base.h>

#include <TRUST_Ref.h>
#include <Matrice.h>

class Probleme_base;
class Domaine_Cl_VEF;
class Domaine_VEF;

class Source_PDF_VEF : public Source_PDF_base
{

  Declare_instanciable(Source_PDF_VEF);

public:
  void associer_pb(const Probleme_base& ) override;
  DoubleTab& ajouter_(const DoubleTab&, DoubleTab&) const override;
  DoubleTab& ajouter_(const DoubleTab&, DoubleTab&, const int) const override ;
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  void verif_ajouter_contrib(const DoubleTab& variable, Matrice_Morse& matrice) const ;
  void multiply_coeff_volume(DoubleTab&) const override;
  void test(Matrice&) const;
  void updateChampRho();

  void correct_pressure(const DoubleTab&,DoubleTab&,const DoubleTab&) const override ;
  void correct_incr_pressure(const DoubleTab&,DoubleTab&) const override ;

  int impr(Sortie&) const override;
  // void ouvrir_fichier_partage(EcrFicPartage&, const Nom&, const Nom&) const;
  // void imprimer_ustar_yplus__mean_only(Sortie&, const Nom& ) const override;

  // Methodes de l interface des champs postraitables
  void creer_champ(const Motcle& motlu) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  void filtre_CLD(DoubleTab&) const override;

protected:
  void calculer_variable_imposee_hybrid()  override;
  void calculer_variable_imposee_elem_fluid() override;
  void calculer_variable_imposee_mean_grad()  override;
  void calculer_vitesse_imposee_power_law_tbl() override;
  void calculer_vitesse_imposee_power_law_tbl_u_star() override;
  void calculer_temperature_imposee_wall_law() override;
//  void calculer_temperature_imposee_mean_wall_law() override;

  OBS_PTR(Domaine_VEF) le_dom_VEF;
  OBS_PTR(Domaine_Cl_VEF) le_dom_Cl_VEF;
  void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base&) override;
  void compute_indicateur_nodal_champ_aire() override;

  DoubleVect tab_u_star_ibm_;                //!< valeurs des u* IBM calculees localement
  DoubleVect tab_y_plus_ibm_;                //!< valeurs des d+ IBM calculees localement

  mutable OWN_PTR(Champ_Fonc_base)  champ_u_star_ibm_;          //!< Champ pour postraitement
  mutable OWN_PTR(Champ_Fonc_base)  champ_y_plus_ibm_;          //!< Champ pour postraitement
};

#endif
