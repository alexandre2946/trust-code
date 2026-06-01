/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Symetrie.h>

Implemente_instanciable(Symetrie,"Symetrie",Navier);
// XD symetrie condlim_base symetrie INHERITS_BRACE 1). For Navier-Stokes equations, this keyword is used to designate a
// XD_CONT symmetry condition concerning the velocity at the boundary called bord (edge) (normal velocity at the edge
// XD_CONT equal to zero and tangential velocity gradient at the edge equal to zero); 2). For scalar transport equation,
// XD_CONT this keyword is used to set a symmetry condition on scalar on the boundary named bord (edge).



/*! @brief Writes the type of the object to an output stream.
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Symetrie::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

/*! @brief Types the boundary field as "Champ_front_uniforme". Does not read anything from the input stream passed as parameter.
 *
 * @param (Entree& s) an input stream
 * @return (Entree& s) the input stream
 */
Entree& Symetrie::readOn(Entree& s )
{
  le_champ_front.typer("Champ_front_uniforme");
  le_champ_front->fixer_nb_comp(0);
  return s ;
}

