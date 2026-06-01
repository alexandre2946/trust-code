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

#include <Bidim_Axi.h>

Implemente_instanciable(Bidim_Axi,"Bidim_Axi",Interprete);
// XD bidim_axi interprete bidim_axi INHERITS_BRACE Keyword allowing a 2D calculation to be executed using axisymetric
// XD_CONT coordinates (R, Z). If this instruction is not included, calculations are carried out using Cartesian
// XD_CONT coordinates.

/*! @brief Simple call to: Interprete::printOn(Sortie&)
 *
 * @param (Sortie& os) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Bidim_Axi::printOn(Sortie& os) const
{
  Cerr << finl;
  Cerr << "we choose the calculation 2D_Axi, reading"<< finl;
  return Interprete::printOn(os);
}


/*! @brief Simple call to: Interprete::readOn(Entree&)
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 */
Entree& Bidim_Axi::readOn(Entree& is)
{
  Cerr << finl;
  Cerr << "we choose the calculation 2D_Axi, reading"<< finl;


  return Interprete::readOn(is);
}

/*! @brief Main function of the Axi interpreter Sets the 2D_axi variable to 1.
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the input stream
 */
Entree& Bidim_Axi::interpreter(Entree& is)
{
  Cerr << finl;
  Cerr << "we choose the calculation 2D_Axi bidim_axi =1  "<< finl;

  bidim_axi = 1;
  return is;
}
