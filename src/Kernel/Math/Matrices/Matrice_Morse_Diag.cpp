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

#include <Matrice_Morse_Diag.h>

Implemente_instanciable_sans_constructeur(Matrice_Morse_Diag,"Matrice_Morse_Diag",Matrice_Morse_Sym);

/*! @brief Writes the three arrays of the Morse storage structure to an output stream.
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie& s) the modified output stream
 */
Sortie& Matrice_Morse_Diag::printOn(Sortie& s) const
{
  return  Matrice_Morse::printOn(s);
}


/*! @brief NOT IMPLEMENTED
 *
 * @param (Entree& s) an input stream
 * @return (Entree& s) the input stream
 * @throws NOT IMPLEMENTED
 */
Entree& Matrice_Morse_Diag::readOn(Entree& s)
{
  return  Matrice_Morse::readOn(s) ;
}

/*! @brief The purpose of this matrix is only to compute the diagonal.
 *
 */
DoubleVect& Matrice_Morse_Diag::ajouter_multvect_(const DoubleVect& x, DoubleVect& y) const
{
  return Matrice_Morse::ajouter_multvect_(x,y);
}

/*! @brief The purpose of this matrix is only to compute the diagonal.
 *
 */
DoubleVect& Matrice_Morse_Diag::ajouter_multvectT_(const DoubleVect& x, DoubleVect& y) const
{
  Cerr<<" not implemented "<<finl;
  exit();
  return y;
}
/*! @brief The purpose of this matrix is only to compute the diagonal.
 *
 */
DoubleTab& Matrice_Morse_Diag::ajouter_multTab_(const DoubleTab& x, DoubleTab& y) const
{
  Cerr<<" not implemented "<<finl;
  exit();
  return y;
}

void Matrice_Morse_Diag::dimensionne_diag(int nb_case_tot)
{
  IntTab indice(nb_case_tot,2);
  for( int c=0; c<nb_case_tot; c++) indice(c,0)=indice(c,1)=c;
  dimensionner( nb_case_tot, 0);
  dimensionner(indice);
}
