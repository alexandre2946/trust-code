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

#include <Lire_Fichier_Bin.h>

Implemente_instanciable(Lire_Fichier_Bin,"Lire_Fichier_Bin|Read_File_Binary",Lire_Fichier);
// XD read_file_bin read_file lire_fichier_bin INHERITS_BRACE Keyword to read an object name_obj in the unformatted type
// XD_CONT file filename.

/*! @brief Calls the printOn method of the Interprete class.
 *
 */
Sortie& Lire_Fichier_Bin::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}

/*! @brief Calls the readOn method of the Interprete class.
 *
 */
Entree& Lire_Fichier_Bin::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

/*! @brief Reads a binary-format file. With 2 arguments nom1 and nom2, reads the object from file nom2 into object nom1.
 *
 *     With a single argument nom1, interprets the file named nom1.
 *
 * @param (Entree& is)
 * @return (Entree&)
 */
Entree& Lire_Fichier_Bin::interpreter(Entree& is)
{
  return Lire_Fichier::interpreter(is);
}
