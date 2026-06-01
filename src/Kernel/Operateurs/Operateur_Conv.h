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

#ifndef Operateur_Conv_included
#define Operateur_Conv_included

#include <Operateur_Conv_base.h>
#include <TRUST_Deriv.h>
#include <Operateur.h>
#include <TRUST_Ref.h>

class Champ_base;

/*! @brief Operateur_Conv Generic class of the hierarchy of operators representing a convection term.
 *
 *     An Operateur_Conv object can reference any object
 *     derived from Operateur_Conv_base.
 *
 * @sa Operateur_Conv_base Operateur
 */
class Operateur_Conv  : public Operateur, public OWN_PTR(Operateur_Conv_base)
{
  Declare_instanciable(Operateur_Conv);
public :

  inline Operateur_base& l_op_base() override;
  inline const Operateur_base& l_op_base() const override;
  DoubleTab& ajouter(const DoubleTab&, DoubleTab& ) const override;
  DoubleTab& calculer(const DoubleTab&, DoubleTab& ) const override;
  inline void associer_vitesse(const Champ_base&) ;
  inline const Champ_base&  vitesse() const ;
  void typer() override;
  inline void typer(const Nom&);
  inline int op_non_nul() const override;
  void associer_norme_vitesse(const Champ_base&);
protected :

  OBS_PTR(Champ_base) la_vitesse;
};

/*! @brief Returns the underlying object upcast to Operateur_base.
 *
 * @return (Operateur_base&) the underlying object upcast to Operateur_base
 */
inline Operateur_base& Operateur_Conv::l_op_base()
{
  return valeur();
}
/*! @brief Returns the underlying object upcast to Operateur_base (const version).
 *
 * @return (Operateur_base&) the underlying object upcast to Operateur_base
 */
inline const Operateur_base& Operateur_Conv::l_op_base() const
{
  return valeur();
}

/*! @brief Associates the velocity (as the transporting velocity) to the convection operator.
 *
 * @param (Champ_Inc_base& vit) the unknown field representing the velocity
 * @return the unknown field representing the transporting velocity
 */
inline void Operateur_Conv::associer_vitesse(const Champ_base& vit)
{
  la_vitesse = vit;
}

/*! @brief Returns the transporting velocity of the operator.
 *
 * @return (Champ_Inc_base&) the unknown field representing the transporting velocity
 */
inline const Champ_base& Operateur_Conv::vitesse() const
{
  return la_vitesse.valeur();
}



/*! @brief Types the operator.
 *
 * @param (Nom& typ) the name representing the type of the operator
 */
inline void Operateur_Conv::typer(const Nom& a_type)
{
  OWN_PTR(Operateur_Conv_base)::typer(a_type);
}

inline int Operateur_Conv::op_non_nul() const
{
  return this->operator bool();
}
#endif
