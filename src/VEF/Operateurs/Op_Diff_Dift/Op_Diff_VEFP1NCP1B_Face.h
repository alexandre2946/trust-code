/****************************************************************************
* Copyright (c) 2024, CEA
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


#ifndef Op_Diff_VEFP1NCP1B_Face_included
#define Op_Diff_VEFP1NCP1B_Face_included

#include <Operateur_Diff_base.h>
#include <Matrice_Morse.h>
#include <Op_Diff_VEF_Face.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <TRUSTLists.h>
#include <SolveurSys.h>
#include <ArrOfBit.h>


/*! @brief class Op_Diff_VEF_Face
 *
 *  @brief This class represents the diffusion operator with VEF discretization.
 *         The diffused field is scalar.
 *
 */
class Op_Diff_VEFP1NCP1B_Face : public Op_Diff_VEF_Face
{
  Declare_instanciable_sans_constructeur(Op_Diff_VEFP1NCP1B_Face);

public:

  Op_Diff_VEFP1NCP1B_Face();
  void associer(const Domaine_dis_base& , const Domaine_Cl_dis_base& ,
                const Champ_Inc_base& ) override;
  void completer() override;

  //Methods for explicit scheme
  DoubleTab& ajouter(const DoubleTab& ,  DoubleTab& ) const override;
  DoubleTab& calculer(const DoubleTab& , DoubleTab& ) const override;
  double calculer_dt_stab() const override;

  //Access to domains
  const Domaine_VEF& domaine_vef() const;
  const Domaine_Cl_VEF& domaine_Cl_VEFPreP1B() const;

  //Methods for implicit scheme.
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  void dimensionner(Matrice_Morse&) const override;
  void ajouter_contribution(const DoubleTab&,Matrice_Morse&) const;
  void contribuer_au_second_membre(DoubleTab& ) const override;


protected :

  // Sizing and filling function
  void initialiser();
  void remplir_nu_p1(const DoubleTab&,DoubleTab&) const;
  void remplir_nu_pA(const DoubleTab&,DoubleTab&) const;

  // Functions for explicit scheme
  DoubleVect& calculer_gradient_elem(const DoubleVect&) const;
  DoubleVect& calculer_gradient_som(const DoubleVect&) const;
  DoubleVect& calculer_gradient_aretes(const DoubleVect&) const;
  DoubleVect& calculer_divergence_elem(DoubleVect&) const;
  DoubleVect& calculer_divergence_som(DoubleVect&) const;
  DoubleVect& calculer_divergence_aretes(DoubleVect&) const;
  DoubleVect& corriger_div_pour_Cl(const DoubleVect&,const DoubleTab&,DoubleVect&) const;
  DoubleTab& corriger_pour_diffusivite(const DoubleTab&,DoubleTab&) const;

  void calculer_laplacien_som(const DoubleTab&) const;

  void calculer_flux_bords_elem(const DoubleVect&) const;
  void calculer_flux_bords_som(const DoubleVect&) const;
  void calculer_flux_bords_aretes(const DoubleVect&) const;

  void calculer_dt_stab_elem(const DoubleTab&,DoubleTab&) const;
  void calculer_dt_stab_som(const DoubleTab&,DoubleTab&) const;
  void calculer_dt_stab_aretes(const DoubleTab&,DoubleTab&) const;

  // Functions for implicit scheme
  void ajouter_contribution_elem(const DoubleTab&,const DoubleVect&, const DoubleTab&,Matrice_Morse&) const;
  void ajouter_contribution_som(const DoubleTab&,const DoubleVect&,const DoubleTab&,Matrice_Morse&) const;
  void ajouter_contribution_aretes(const DoubleTab&,const DoubleVect&,const DoubleTab&,Matrice_Morse&) const;

  void coeff_matrice_som(const int,IntVect&,DoubleTab&,DoubleTab&,const DoubleVect&,const DoubleTab&,
                         const DoubleTab&,Matrice_Morse&) const;
  void coeff_matrice_som_CL(const int,IntVect&,DoubleTab&,DoubleTab&,const DoubleVect&,const DoubleTab&,
                            const DoubleTab&,Matrice_Morse&) const;
  void coeff_matrice_som_symetrie(const int,IntVect&,DoubleTab&,DoubleTab&,const DoubleVect&,const DoubleTab&,
                                  const DoubleTab&,Matrice_Morse&) const;
  void coeff_matrice_som_perio(const int,const int,IntVect&,DoubleTab&,DoubleTab&,
                               const DoubleVect&,const DoubleTab&,const DoubleTab&,Matrice_Morse&) const;


  void gradient_som(const int face,const int,const int,const int,const int,const int,DoubleTab&) const;
  void gradient_som(const int,int&,IntVect&,DoubleTab&) const;
  void gradient_som_CL(const int,int&,IntVect&,DoubleTab&) const;

  void isInStencil(int,int,int&,int&,int&,int&) const;
  void isFaceOfSymetry(ArrOfBit&,int&) const;

  void liste_face(IntLists&,int&) const;

  // Test function
  void test() const;
  void corriger_Cl_test(DoubleTab&) const;
  int test_ = 0;

  // Class attributes
  mutable DoubleTab gradient_p0_;
  mutable DoubleTab gradient_p1_;
  mutable DoubleTab gradient_pa_;
  mutable Matrice_Morse laplacien_p1_;
  mutable int is_laplacian_filled_ = 0;
  mutable int is_laplacian_built_ = 0;

  double convexite_ = 1.e-3; //attribute defining the convexity
  double coeff_ = 1.; //coefficient for the lumped matrix
  int alphaE = 1; //attribute for computing the P0 part of the operator
  int alphaS = 1; //attribute for computing the P1 part of the operator
  int alphaA = 0; //attribute for computing the Pa part of the operator
  int dim_ch_ = -1; //attribute giving the number of dimensions of the unknown
  int decentrage_ = 1; //attribute enabling upwinding of the operator
};

/* //Description: */
/* //Assemble the matrix. */

inline void Op_Diff_VEFP1NCP1B_Face::contribuer_a_avec(const DoubleTab& inco,
                                                       Matrice_Morse& matrice) const
{
  ajouter_contribution(inco, matrice);
}

/* //Description: */
/* //add the right-hand side contribution. */

inline void Op_Diff_VEFP1NCP1B_Face::contribuer_au_second_membre(DoubleTab& resu) const
{
  contribue_au_second_membre(resu);
}

inline const Domaine_VEF& Op_Diff_VEFP1NCP1B_Face::domaine_vef() const
{
  return ref_cast(Domaine_VEF,le_dom_vef.valeur());
}

inline const Domaine_Cl_VEF& Op_Diff_VEFP1NCP1B_Face::domaine_Cl_VEFPreP1B() const
{
  return ref_cast(Domaine_Cl_VEF,la_zcl_vef.valeur());
}
#endif
