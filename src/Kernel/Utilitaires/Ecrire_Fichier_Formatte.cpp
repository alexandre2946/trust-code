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
#include <Ecrire_Fichier_Formatte.h>

Implemente_instanciable(Ecrire_Fichier_Formatte,"Ecrire_Fichier_Formatte",Ecrire_Fichier);
// XD ecrire_fichier_formatte write_file ecrire_fichier_formatte INHERITS_BRACE Keyword to write the object of name
// XD_CONT name_obj to a file filename in ASCII format.


/*! @brief Call to the printOn method of the Interprete class.
 *
 */
Sortie& Ecrire_Fichier_Formatte::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}

/*! @brief Call to the readOn method of the Interprete class.
 *
 */
Entree& Ecrire_Fichier_Formatte::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

/*! @brief Sequential interpreter.
 *
 * Syntax: Ecrire_Fichier_Formatte object_name file_name
 *
 */
Entree& Ecrire_Fichier_Formatte::interpreter(Entree& is)
{
  Ecrire_Fichier::interpreter_(is, 0 /* formatte */);
  return is;
}
