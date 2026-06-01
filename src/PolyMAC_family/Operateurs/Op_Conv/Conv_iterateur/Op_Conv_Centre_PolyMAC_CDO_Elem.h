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

#ifndef Op_Conv_Centre_PolyMAC_CDO_Elem_included
#define Op_Conv_Centre_PolyMAC_CDO_Elem_included

#include <Op_Conv_PolyMAC_CDO_iterateur_base.h>
#include <Eval_centre_PolyMAC_CDO_Elem.h>
#include <Iterateur_PolyMAC_CDO_Elem.h>
#include <Op_PolyMAC_CDO_Elem.h>

/*! @brief class Op_Conv_Centre_PolyMAC_CDO_Elem
 *
 *   Convection operator for a scalar transport equation using the PolyMAC_CDO
 *   discretization. The convected field is scalar, the convection scheme is
 *   centered, the associated iterator is Iterateur_PolyMAC_CDO_Elem, and the
 *   associated evaluator is Eval_centre_PolyMAC_CDO_Elem.
 *
 */

class Op_Conv_Centre_PolyMAC_CDO_Elem: public Op_Conv_PolyMAC_CDO_iterateur_base, public Op_PolyMAC_CDO_Elem
{

  Declare_instanciable_sans_constructeur(Op_Conv_Centre_PolyMAC_CDO_Elem);

public:

  Op_Conv_Centre_PolyMAC_CDO_Elem();
  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;
  void associer_vitesse(const Champ_base&) override;
  const Champ_base& vitesse() const override;
  Champ_base& vitesse() override;
  inline void dimensionner(Matrice_Morse&) const override;
  inline void dimensionner_bloc_vitesse(Matrice_Morse&) const override;
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;

protected:
  inline Op_Conv_Centre_PolyMAC_CDO_Elem(const Iterateur_PolyMAC_CDO_base&);
};

// This constructor allows creating derived classes (e.g., front_tracking)
inline Op_Conv_Centre_PolyMAC_CDO_Elem::Op_Conv_Centre_PolyMAC_CDO_Elem(const Iterateur_PolyMAC_CDO_base& it) :
  Op_Conv_PolyMAC_CDO_iterateur_base(it)
{
}

/*! @brief Sizes the matrix.
 *
 */
inline void Op_Conv_Centre_PolyMAC_CDO_Elem::dimensionner(Matrice_Morse& matrice) const
{
  Op_PolyMAC_CDO_Elem::dimensionner(iter_->domaine(), iter_->domaine_Cl(), matrice);
}

inline void Op_Conv_Centre_PolyMAC_CDO_Elem::dimensionner_bloc_vitesse(Matrice_Morse& matrice) const
{
  Op_PolyMAC_CDO_Elem::dimensionner_bloc_vitesse(iter_->domaine(), iter_->domaine_Cl(), matrice);
}

inline void Op_Conv_Centre_PolyMAC_CDO_Elem::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  Op_PolyMAC_CDO_Elem::modifier_pour_Cl(iter_->domaine(), iter_->domaine_Cl(), matrice, secmem);
}

#endif
