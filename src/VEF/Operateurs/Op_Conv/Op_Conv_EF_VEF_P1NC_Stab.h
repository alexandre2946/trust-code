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

#ifndef Op_Conv_EF_VEF_P1NC_Stab_included
#define Op_Conv_EF_VEF_P1NC_Stab_included

#include <Op_Conv_VEF_Face.h>
#include <Equation_base.h>
#include <TRUSTTabs.h>
#include <TRUSTList.h>
#include <TRUST_Ref.h>

class Sous_domaine_VF;
class Matrice_Morse;
/*! @brief class Op_Conv_EF_VEF_P1NC_Stab
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


class Op_Conv_EF_VEF_P1NC_Stab : public Op_Conv_VEF_Face
{

  Declare_instanciable(Op_Conv_EF_VEF_P1NC_Stab);

public:

  //Auxiliary methods
  void remplir_fluent() const override;
  int is_compressible() const;
  void completer() override;

  //Methods for explicit scheme
  DoubleTab& ajouter(const DoubleTab& , DoubleTab& ) const override;

  //Methods for implicit scheme
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  void ajouter_contribution(const DoubleTab&, Matrice_Morse&) const override;

  //test
  void         modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;

  public_for_cuda
  void calculer_flux_bords(const DoubleTab&, const DoubleTab&, const DoubleTab&) const;
  void calculer_coefficients_operateur_centre(DoubleTab&,const int, const DoubleTab& vitesse) const;
  DoubleTab& ajouter_partie_compressible(const DoubleTab&, DoubleTab&, const DoubleTab& vitesse) const;
  DoubleTab& ajouter_operateur_centre(const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  DoubleTab& ajouter_diffusion(const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  DoubleTab& ajouter_antidiffusion(const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  void ajouter_contribution_operateur_centre(const DoubleTab&, const DoubleTab&, Matrice_Morse&) const;
  void ajouter_contribution_diffusion(const DoubleTab&, const DoubleTab&, Matrice_Morse&) const;
  void mettre_a_jour_pour_periodicite(DoubleTab&) const;

private :

  //Auxiliary methods

  //Methods for explicit scheme
  void reinit_conv_pour_Cl(const DoubleTab&,const IntList&, const DoubleTabs&, const DoubleTab&, DoubleTab&) const;

  KOKKOS_INLINE_FUNCTION void calculer_senseur(CDoubleTabView3, CDoubleArrView, const int, const int, CIntTabView, CIntTabView, CIntTabView, double*, double*, double*, double*) const;
  inline void calculer_senseur(const DoubleTab&, const DoubleVect&, const int, const int, const IntTab&, const IntTab&, const IntTab&, ArrOfDouble&, ArrOfDouble&, ArrOfDouble&, ArrOfDouble&) const;
  void ajouter_old(const DoubleTab& , DoubleTab&, const DoubleTab& vitesse) const;
  void calculer_data_pour_dirichlet();

  //Methods for implicit scheme
  void ajouter_contribution_antidiffusion(const DoubleTab&,const DoubleTab&,Matrice_Morse&) const;
  void ajouter_contribution_partie_compressible(const DoubleTab&,const DoubleTab&,Matrice_Morse&) const;

  //Test methods
  void test(const DoubleTab&,const DoubleTab&, const DoubleTab& vitesse) const;
  void test_difference_Kij(const DoubleTab&,DoubleTab&,DoubleTab&, const DoubleTab& vitesse) const;
  void test_difference_resu(const DoubleTab&,const DoubleTab&,const DoubleTab&,const DoubleTab&,const DoubleTab& vitesse) const;
  void test_implicite() const;

  //Class attributes
  ArrOfInt elem_nb_faces_dirichlet_;
  IntTab elem_faces_dirichlet_;
  ArrsOfInt elem_faces_frontiere;

  ArrOfDouble alpha_tab_;
  ArrOfDouble beta_; // equals zero for faces where upwind (Amont) degeneration is desired.
  //  mutable DoubleTab limiteurs_;//array storing for each face the algebraic mean of the limiter

  double alpha_ = 1.;

  int is_compressible_ = 0;
  int test_ = 0;
  int old_ = 0;
  int volumes_etendus_ = 1;

  bool sous_domaine = false;  // Sub-domain case to define so that EF_Stab degenerates to upwind (Amont)
  int new_jacobienne_ = 0;
  Nom nom_sous_domaine;
  OBS_PTR(Sous_domaine_VF) le_sous_domaine_dis;

  int nb_ssz_alpha = -1;
  DoubleVect alpha_ssz;
  Noms noms_ssz_alpha;
  bool ssz_alpha = false;
};

inline void Op_Conv_EF_VEF_P1NC_Stab::contribuer_a_avec(const DoubleTab& inco,
                                                        Matrice_Morse& matrice) const
{
  ajouter_contribution(inco, matrice);
}

inline int Op_Conv_EF_VEF_P1NC_Stab::is_compressible() const
{
  return is_compressible_;
}

#endif
