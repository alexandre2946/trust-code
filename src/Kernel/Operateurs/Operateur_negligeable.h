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

#ifndef Operateur_negligeable_included
#define Operateur_negligeable_included

#include <TRUSTTab.h>

class Domaine_Cl_dis_base;
class Domaine_dis_base;
class Champ_Inc_base;

/*! @brief Classe Opnegligeable This class defines the interface of a negligible operator
 *
 *     by defining ajouter, calculer and mettre_a_jour methods
 *     that perform no computation.
 *     This class is used when defining negligible operators
 *     through multiple inheritance.
 *
 * @sa Op_Diff_negligeable Op_Diff_K_Eps_negligeable Op_Conv_negligeable, Interface class, outside the TRUST operator hierarchy., Must be part of a multiple inheritance to be, useful for defining an instantiable negligible operator.
 */
class Operateur_negligeable
{
public :

  inline DoubleTab& ajouter(const DoubleTab&, DoubleTab& ) const;
  inline DoubleTab& calculer(const DoubleTab&, DoubleTab& ) const;
  inline void mettre_a_jour(double );

protected :
  inline void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&);
};



/*! @brief Adds the contribution of a negligible operator to an array.
 *
 *     DOES NOTHING, returns the array parameter unmodified.
 *
 * @param (DoubleTab&) the array on which the operator is applied
 * @param (DoubleTab& x)
 * @return (DoubleTab&) the input parameter x unmodified
 */
inline DoubleTab& Operateur_negligeable::ajouter(const DoubleTab&, DoubleTab& x) const
{
  return x;
}


/*! @brief Initializes the array parameter with the contribution of the negligible operator: initializes the array to ZERO.
 *
 * @param (DoubleTab&) the array on which the operator is applied
 * @param (DoubleTab& x)
 * @return (DoubleTab&) the input array set to zero
 */
inline DoubleTab& Operateur_negligeable::calculer(const DoubleTab&,
                                                  DoubleTab& x) const
{
  return x=0.0;
}


/*! @brief Time update of a negligible operator: DOES NOTHING
 *
 * @param (double)
 */
inline void Operateur_negligeable::mettre_a_jour(double )
{
}


/*! @brief Associates various objects to a negligible operator: DOES NOTHING
 *
 * @param (Domaine_dis_base&)
 * @param (Domaine_Cl_dis_base&)
 * @param (Champ_Inc_base&)
 */
inline void Operateur_negligeable::associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&)
{
}
#endif
