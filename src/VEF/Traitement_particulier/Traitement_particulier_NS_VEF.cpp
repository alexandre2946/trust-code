/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Traitement_particulier_NS_VEF.h>

Implemente_instanciable(Traitement_particulier_NS_VEF,"Traitement_particulier_NS_VEF",Traitement_particulier_NS_base);


/*! @brief Prints the equation to an output stream.
 *
 * @brief Simple call to Equation_base::printOn(Sortie&).
 *
 * @param is output stream
 * @return modified output stream
 */
Sortie& Traitement_particulier_NS_VEF::printOn(Sortie& is) const
{
  return is;
}


/*! @brief Read the Navier-Stokes equation specifications from an input stream.
 *
 *     Simple call to Navier_Stokes_std::readOn(Entree&)
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 * @throws no turbulence model specified
 */
Entree& Traitement_particulier_NS_VEF::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_VEF::lire(Entree& is)
{
  return is;
}
void Traitement_particulier_NS_VEF::associer_eqn(const Equation_base& eqn)
{
  Traitement_particulier_NS_base::associer_eqn(eqn);
}
