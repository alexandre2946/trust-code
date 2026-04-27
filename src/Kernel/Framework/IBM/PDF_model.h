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

#ifndef PDF_model_included
#define PDF_model_included
#include <Ecrire_MED.h>

#include <TRUSTTabs_forward.h>

#include <Parser_U.h>
#include <Domaine_VF.h>
#include <Probleme_base.h>
#include <Motcle.h>

class Schema_Temps_base;

/*! @brief : class PDF_model
 *
 *  <Description of class PDF_model>
 */
class PDF_model : public Objet_U
{
  Declare_instanciable(PDF_model) ;
public :
  double get_variable_imposee(ArrOfDouble&,int);
  void affecter_variable_imposee(Domaine_VF&, const DoubleTab&);
  void affecter_vitesse_shape_IBM(Domaine_VF&, const DoubleTab&, double);
  double eta() const { return eta_; }
  int pdf_bilan() const { return pdf_bilan_; }
  bool get_PDF_mobile() const { return PDF_mobile_; }
  bool get_use_pseudo_level_set_moving_PDF() const { return use_pseudo_level_set_moving_PDF_; }
  inline void set_PDF_mobile(bool the_flag) { PDF_mobile_ = the_flag; }
  void discretiser_vitesse_shape_IBM(const Probleme_base&);
  inline const DoubleTab& get_vitesse_shape_IBM() const { return vitesse_shape_IBM_->valeurs(); }
  inline void set_vitesse_shape_IBM(DoubleTab& it) { vitesse_shape_IBM_->valeurs() = it; }
  inline int is_vitesse_PDF_donnee( ) const { return vitesse_PDF_donnee_; }
  double raid() const { return raid_; }

protected :
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  OWN_PTR(Champ_Don_base) variable_imposee_lu_;
  OWN_PTR(Champ_Don_base) variable_imposee_;
  int dim_variable_=-1;
  double eta_ = -100.;
  double coefku_= -100.;
  double temps_relax_=1.0e+12;
  double echelle_relax_=5.0e-2;
  double regul_coeff_PDF_=1.;
  int type_variable_imposee_= -1;
  bool local_ = false;
  int pdf_bilan_ = 0; // 0: terme pdf; 1: termes pdf + temps; 2: termes pdf + temps + conv
  bool PDF_mobile_ = false;
  bool use_pseudo_level_set_moving_PDF_ = false;
  int vitesse_PDF_donnee_ = false;
  OWN_PTR(Champ_Don_base) vitesse_shape_IBM_;
  double raid_ = 0.;

  friend class Source_PDF_base;
  friend class Source_PDF_EF;
  friend class Source_PDF_VEF;
  friend class Source_PDF_VDF;

private:
  VECT(Parser_U) parsers_;
  VECT(Parser_U) parsers_vitesse_shape_;

  Ecrire_MED ecr_med_vitesse_shape_IBM_;
};

#endif /* PDF_model */
