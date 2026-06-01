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

#include <Bord.h>

Implemente_instanciable_32_64(Bord_32_64,"Bord",Frontiere_32_64<_T_>);

// XD defbord objet_lecture nul NO_BRACE Class to define an edge.

// XD defbord_2 defbord nul INHERITS_BRACE 1-D edge (straight line) in the 2-D space.
// XD attr dir chaine(into=["X","Y"]) dir REQ Edge is perpendicular to this direction.
// XD attr eq chaine(into=["="]) eq REQ Equality sign.
// XD attr pos floattant pos REQ Position value.
// XD attr pos2_min floattant pos2_min REQ Minimal value.
// XD attr inf1 chaine(into=["<="]) inf1 REQ Less than or equal to sign.
// XD attr dir2 chaine(into=["X","Y"]) dir2 REQ Edge is parallel to this direction.
// XD attr inf2 chaine(into=["<="]) inf2 REQ Less than or equal to sign.
// XD attr pos2_max floattant pos2_max REQ Maximal value.

// XD defbord_3 defbord nul INHERITS_BRACE 2-D edge (plane) in the 3-D space.
// XD attr dir chaine(into=["X","Y","Z"]) dir REQ Edge is perpendicular to this direction.
// XD attr eq chaine(into=["="]) eq REQ Equality sign.
// XD attr pos floattant pos REQ Position value.
// XD attr pos2_min floattant pos2_min REQ Minimal value.
// XD attr inf1 chaine(into=["<="]) inf1 REQ Less than or equal to sign.
// XD attr dir2 chaine(into=["X","Y"]) dir2 REQ Edge is parallel to this direction.
// XD attr inf2 chaine(into=["<="]) inf2 REQ Less than or equal to sign.
// XD attr pos2_max floattant pos2_max REQ Maximal value.
// XD attr pos3_min floattant pos3_min REQ Minimal value.
// XD attr inf3 chaine(into=["<="]) inf3 REQ Less than or equal to sign.
// XD attr dir3 chaine(into=["Y","Z"]) dir3 REQ Edge is parallel to this direction.
// XD attr inf4 chaine(into=["<="]) inf4 REQ Less than or equal to sign.
// XD attr pos3_max floattant pos3_max REQ Maximal value.

// XD bord bord_base bord NO_BRACE The block side is not in contact with another block and boundary conditions are
// XD_CONT applied to it.
// XD attr nom chaine nom REQ Name of block side.
// XD attr defbord defbord defbord REQ Definition of block side.

/*! @brief Simple call to: Frontiere::printOn(Sortie&)
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the modified output stream
 */
template <typename _SIZE_>
Sortie& Bord_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  return Frontiere_32_64<_SIZE_>::printOn(s) ;
}

/*! @brief Simple call to: Frontiere::readOn(Entree&)
 *
 * @param (Entree& s) an input stream
 * @return (Entree&) the modified input stream
 */
template <typename _SIZE_>
Entree& Bord_32_64<_SIZE_>::readOn(Entree& s)
{
  return Frontiere_32_64<_SIZE_>::readOn(s) ;
}


template class Bord_32_64<int>;
#if INT_is_64_ == 2
template class Bord_32_64<trustIdType>;
#endif
