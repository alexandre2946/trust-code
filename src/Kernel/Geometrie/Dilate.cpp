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

#include <Dilate.h>
#include <Domaine.h>

Implemente_instanciable(Dilate,"Dilate",Interprete_geometrique_base);
// XD dilate interprete dilate INHERITS_BRACE Keyword to multiply the whole coordinates of the geometry.
// XD attr domain_name ref_domaine domain_name REQ Name of domain.
// XD attr alpha floattant alpha REQ Value of dilatation coefficient.


/*! @brief Simple call to: Interprete::printOn(Sortie&)
 *
 * @param (Sortie& os) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Dilate::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}


/*! @brief Simple call to: Interprete::readOn(Entree&)
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 */
Entree& Dilate::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

/*! @brief Main function of the Dilate interpreter Data set structure (in dimension 2):
 *
 *     Dilate dom alpha
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the input stream
 * @throws the object to be meshed is not of Domaine type
 */
Entree& Dilate::interpreter_(Entree& is)
{
  double alpha;
  associer_domaine(is);
  is >> alpha;
  if (Process::is_parallel())
    {
      Cerr << "Dilate can not be use in parallel." << finl;
      Cerr << "Put rather this interpreter in the " << finl;
      Cerr << "data file of the mesh splitter." << finl;
      exit();
    }

  DoubleTab& coord=domaine().les_sommets();
  coord*=alpha;
  Cerr << "Expands the domain " << domaine().le_nom() << " of a factor " << alpha << finl;

  return is;
}

