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

#ifndef Pb_Couple_Optimisation_IBM_included
#define Pb_Couple_Optimisation_IBM_included

#include <Probleme_Couple.h>
#include <Source_PDF_base.h>
#include <Interpolation_IBM_base.h>
#include <TRUST_Ref.h>
#include <Champ_Don_base.h>
#include <Ecrire_MED.h>

class Cond_lim_base;

class Schema_Temps_base;
class Discretisation_base;

class Pb_Couple_Optimisation_IBM: public Probleme_Couple
{
  Declare_instanciable(Pb_Couple_Optimisation_IBM);

public:
  void initialize() override;
  inline  const Interpolation_IBM_base& le_modele_interpolation_IBM()  const { return my_interpolation_opt_;}
  inline  const Source_PDF_base& la_source_PDF_IBM()  const { return my_source_PDF_opt_;}
  inline  const Source_PDF_base& la_source_PDF_IBM_adjt()  const { return my_source_PDF_opt_adjt_;}
  inline  const Champ_Don_base& le_champ_aire_IBM()  const { return my_source_PDF_opt_->get_champ_aire();}
  inline  const Champ_Don_base& le_champ_rotation_IBM()  const { return my_source_PDF_opt_->get_champ_rotation();}
  double alpha() const { return alpha_; }
  bool solveTimeStep() override;
  void calcul_derivee_forme_IBM();
  bool initTimeStep(double) override;

protected:
  void Save_Med_File_fonction_cout(DoubleTab&);
  Ecrire_MED ecr_med_fonction_cout_;
  int visu_cout_ = 0;

  OBS_PTR(Interpolation_IBM_base) my_interpolation_opt_;
  OBS_PTR(Interpolation_IBM_base) my_interpolation_opt_adjt_;
  OBS_PTR(Prepro_IBM_base) my_prepro_opt_;
  OBS_PTR(Source_PDF_base) my_source_PDF_opt_;
  OBS_PTR(Source_PDF_base) my_source_PDF_opt_adjt_;
  OBS_PTR(Probleme_base) pb_etat_opt_;
  OBS_PTR(Probleme_base) pb_adjt_opt_;
  OBS_PTR(Probleme_base) pb_projection_opt_;
  OWN_PTR(Champ_Don_base) source_derivee_forme_; // Second membre scalaire equ projection (scalaire par element)
  OWN_PTR(Champ_Don_base) normal_derivee_forme_; // Second membre scalaire equ projection (vecteur norme par element)
  OWN_PTR(Champ_Don_base) fonction_cout_lu_,fonction_cout_; // Fonction cout scalaire (par element)
  double alpha_ = 1.;
  double pond_shap_deriv_for_proj_ = 1.;
  int regul_PDF_shape_deriv_ = 0 ;
  int Numero_eq_optimis_;
  int Numero_src_deriv_form_ = -1;
  double modif_aire_pc_low_ = -1.0e+6;
  double modif_aire_pc_high_ = 1.0e+6;
  double area_ref_ = 0.;
  int area_ref_set_ =0;
  int nb_save_call_ = 0;
};
#endif
