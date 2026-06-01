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


#ifndef Op_Conv_Amont_old_VEF_Face_included
#define Op_Conv_Amont_old_VEF_Face_included

#include <Op_Conv_VEF_base.h>


/*! @brief class Op_Conv_Amont_old_VEF_Face
 *
 *   This class represents the convection operator associated with a scalar transport equation.
 *   The discretization is VEF.
 *   The convected field is a scalar or vector of type Champ_P1NC.
 *   The convection scheme is of upwind (Amont) type.
 *   Methods for the implicit scheme are implemented.
 *
 *
 * @sa Operateur_Conv_base
 */
class Op_Conv_Amont_old_VEF_Face : public Op_Conv_VEF_base
{

  Declare_instanciable(Op_Conv_Amont_old_VEF_Face);

public:

  DoubleTab& ajouter(const DoubleTab& , DoubleTab& ) const override;

  // Methods for implicit assembly.
  inline void dimensionner(Matrice_Morse& ) const override;
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;
  inline void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  inline void contribuer_au_second_membre(DoubleTab& ) const override;
  void contribue_au_second_membre(DoubleTab& ) const;
  void ajouter_contribution(const DoubleTab&, Matrice_Morse& ) const;

protected:
// DoubleVect porosite_face;
};


/*! @brief Size the matrix using the dimensionner method of class Op_VEF_Face.
 *
 */
inline  void Op_Conv_Amont_old_VEF_Face::dimensionner(Matrice_Morse& matrice) const
{
  Op_VEF_Face::dimensionner(le_dom_vef.valeur(),la_zcl_vef.valeur(), matrice);
}


/*! @brief Modify the right-hand side and the matrix for Dirichlet conditions.
 *
 */
inline void Op_Conv_Amont_old_VEF_Face::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  Op_VEF_Face::modifier_pour_Cl(le_dom_vef.valeur(),la_zcl_vef.valeur(), matrice, secmem);
}


/*! @brief Assemble the implicit unknown matrix.
 *
 */
inline void Op_Conv_Amont_old_VEF_Face::contribuer_a_avec(const DoubleTab& inco,
                                                          Matrice_Morse& matrice) const
{
  ajouter_contribution(inco, matrice);
}

/*! @brief Add the contribution to the right-hand side.
 *
 */
inline void Op_Conv_Amont_old_VEF_Face::contribuer_au_second_membre(DoubleTab& resu) const
{
  contribue_au_second_membre(resu);
}

#endif
