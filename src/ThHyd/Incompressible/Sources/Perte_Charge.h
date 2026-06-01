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

#ifndef Perte_Charge_included
#define Perte_Charge_included


/*! @brief Classe Perte_Charge This class represents a pressure drop term that is introduced
 *
 *     in the Navier-Stokes equations to model the presence of a
 *     certain type of obstacle in the flow.
 *     The notion of pressure drop is used for obstacles on which
 *     friction is not computed (obstacles internal to cells, or
 *     grids for example).
 *     An object of type Perte_Charge applies in
 *     a single spatial direction (direction_perte_charge() >= 0)
 *     or in all directions (direction_perte_charge() == -1).
 *
 * @sa Does not derive from Objet_U
 */
class Perte_Charge
{

public :

  inline int direction_perte_charge() const;

protected :

  int direction_perte_charge_;
};


/*! @brief Returns the pressure drop direction.
 *
 * @brief Returns the pressure drop direction (X, Y or Z), or -1 for all directions.
 * @return (int) the pressure drop direction
 */
inline int Perte_Charge::direction_perte_charge() const
{
  return direction_perte_charge_;
}

#endif
