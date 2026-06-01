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

#ifndef Op_Conv_Muscl_New_VEF_Face_included
#define Op_Conv_Muscl_New_VEF_Face_included

#include <Op_Conv_VEF_Face.h>
#include <TRUSTTabs.h>
#include <TRUSTList.h>
#include <TRUST_Ref.h>
class Matrice_Morse;
class Sous_domaine_VF;

/*! @brief class Op_Conv_Muscl_New_VEF_Face
 *
 *   This class represents the convection operator associated with a scalar transport equation.
 *   The discretization is VEF.
 *   The convected field is a scalar or vector of type Champ_P1NC.
 *   The convection scheme is derived from the paper
 *   "High-resolution FEM-TVD schemes based on a fully multidimensional flux limiter"
 *    D.Kuzmin and S.Turek.
 *   Inherits from Op to recover upwind implicitation.
 *
 *
 * @sa Operateur_Conv_base
 */


class Op_Conv_Muscl_New_VEF_Face : public Op_Conv_VEF_Face
{

  Declare_instanciable(Op_Conv_Muscl_New_VEF_Face);

public:

  //Auxiliary methods
  void remplir_fluent() const override;
  //int is_compressible() const { return 1;/* is_compressible_; */ }

  void completer() override;

  //Methods for explicit assembly
  DoubleTab& ajouter(const DoubleTab& , DoubleTab& ) const override;
  double calculer_dt_stab() const override;

  //Methods for implicit assembly
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  void ajouter_contribution(const DoubleTab&, Matrice_Morse&) const override;

  //test
  void         modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;

  public_for_cuda
  void calculer_flux_bords(const DoubleTab&, const DoubleTab&, const DoubleTab&) const;
  void calculer_coefficients_operateur_centre(DoubleTab&,DoubleTab&,DoubleTab&,DoubleTab&,const int, const DoubleTab& vitesse) const;
  void calculer_flux_operateur_centre(DoubleTab&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const int,const DoubleTab&,const DoubleTab&) const;
  DoubleTab& ajouter_operateur_centre(const DoubleTab&,const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  DoubleTab& ajouter_diffusion(const DoubleTab&,const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  DoubleTab& ajouter_antidiffusion_v1(const DoubleTab&, const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  DoubleTab& ajouter_antidiffusion_v2(const DoubleTab&, const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  void mettre_a_jour_pour_periodicite(const DoubleTab&,const DoubleTab&,DoubleTab&) const;

private :

  //Auxiliary methods

  //Methods for explicit assembly

  void modifier_flux_operateur_centre(DoubleTab&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const int,const DoubleTab&,const DoubleTab&) const;

  DoubleTab& ajouter_antidiffusion(const DoubleTab&, const DoubleTab&, const DoubleTab&, DoubleTab&) const;


  KOKKOS_INLINE_FUNCTION void calculer_senseur(CDoubleTabView3, CDoubleTabView4, CDoubleArrView, const int, const int, const int, CIntTabView, CIntTabView, CIntTabView, double&, double&, double&, double&) const;
  void calculer_data_pour_dirichlet();

  //Class attributes
  IntTab is_element_for_upwinding_;
  IntVect is_dirichlet_faces_;

  //ArrOfDouble alpha_tab;
  //ArrOfDouble beta_; // zero on faces where upwind (Amont) degeneration is desired.
  //  mutable DoubleTab limiteurs_;//array storing for each face the algebraic mean of the limiter

  double max_limiteur_ = 1.;
  int centered_ = 1;
  int upwind_ = 1;
  int stabilized_ = 1;
  int old_centered_ = 0;
  int version_ = 2;
  int facsec_auto_ = 0;

  // bool sous_domaine;  // Sub-domain case to define so that EF_Stab degenerates to upwind (Amont)
  Nom nom_sous_domaine;
  //OBS_PTR(Sous_domaine_VF) le_sous_domaine_dis;

  //double (*limiteur_)(double)=nullptr;

};

inline void Op_Conv_Muscl_New_VEF_Face::contribuer_a_avec(const DoubleTab& inco,
                                                          Matrice_Morse& matrice) const
{
  ajouter_contribution(inco, matrice);
}

#endif
