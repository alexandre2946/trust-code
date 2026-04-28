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

#ifndef Interpolation_IBM_base_included
#define Interpolation_IBM_base_included

#include <Discretisation_base.h>
#include <Domaine_dis_base.h>
#include <Probleme_base.h>
#include <Prepro_IBM_base.h>
#include <TRUSTList.h>
#include <Param.h>

class Interpolation_IBM_base : public Objet_U
{

  Declare_base( Interpolation_IBM_base ) ;

public:
  virtual void discretise(const Discretisation_base&, Domaine_dis_base&);
  virtual void discretise_PDF_mobile(const Discretisation_base&, Domaine_dis_base&);
  virtual void discretise_pseudo_level_set(const Discretisation_base&, Domaine_dis_base&);

  inline bool get_impr() { return impr_; }
  inline int get_N_histo() { return N_histo_; }
  inline const Champ_Don_base& get_solid_points() { return solid_points_; }
  void set_param(Param&) const override;
  void calculer_normal_et_distance_proj_solid();
  void definir_pseudo_level_set();
  void calcul_cluster_pseudo_level_set(IntLists&, int, IntTab&, DoubleTab&, IntList&, int);
  void define_pseudo_level_set_for_one_cut_cell(IntLists&, int, IntTab&, DoubleTab&, IntList&, int);
  inline void set_source(Source_PDF_base& the_source) {my_source_ = the_source ;}
  inline const DoubleTab& get_normal_proj_solid() const { return champ_normal_proj_solid_->valeurs();}
  inline const DoubleTab& get_champ_solid_points() {return solid_points_->valeurs();}
  inline const DoubleTab& get_champ_corresp_elems() {return corresp_elems_->valeurs();}
  inline const DoubleTab& get_champ_dis_proj_solid() {return champ_dis_proj_solid_->valeurs();}
  inline const DoubleTab& get_champ_champ_normal_proj_solid() {return champ_normal_proj_solid_->valeurs();}
  inline const DoubleTab& get_pseudo_level_set() {return pseudo_level_set_->valeurs();}
  inline void set_pseudo_level_set(DoubleTab& it) { pseudo_level_set_->valeurs() = it;}
  virtual void set_fields_from_prepro_to_interp(Prepro_IBM_base&);

protected:

  OBS_PTR(Source_PDF_base) my_source_;
  OWN_PTR(Champ_Don_base) is_dirichlet_lu_;
  OWN_PTR(Champ_Don_base) is_dirichlet_;
  OWN_PTR(Champ_Don_base) solid_points_lu_;
  OWN_PTR(Champ_Don_base) solid_points_;
  OWN_PTR(Champ_Don_base) corresp_elems_lu_;
  OWN_PTR(Champ_Don_base) corresp_elems_;
  bool has_corresp_ = false;
  bool impr_= false;  // Default value
  int N_histo_=10;  // Default value for number of histogram boxes for printed data

  OWN_PTR(Champ_Don_base) champ_dis_proj_solid_;
  OWN_PTR(Champ_Don_base) champ_normal_proj_solid_;
  OWN_PTR(Champ_Don_base) pseudo_level_set_;
  bool solid_points_from_prepro_ = false;
  bool is_dirichlet_from_prepro_ = false;
  bool corresp_elems_from_prepro_ = false;

  friend class Source_PDF_base;
  friend class Source_PDF_EF;
  friend class Source_PDF_VDF;
  friend class Source_PDF_VEF;
};

#endif /* Interpolation_IBM_base_included */
