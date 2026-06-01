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

#include <Lire.h>
#include <Synonyme_info.h>
#include <Interprete_bloc.h>

Implemente_instanciable(Lire,"Lire|Read",Interprete);
// XD read interprete lire NO_BRACE Interpretor to read the a_object objet defined between the braces.
// XD attr a_object chaine a_object REQ Object to be read.
// XD attr bloc chaine bloc REQ Definition of the object.


/*! @brief Calls the printOn method of the Interprete class.
 *
 */
Sortie& Lire::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}

/*! @brief Calls the readOn method of the Interprete class.
 *
 */
Entree& Lire::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

/*! @brief Read an object.
 *
 */
Entree& Lire::interpreter(Entree& is)
{
  // Example in a data file:
  // Lire name { ... }
  Nom name;
  is >> name; // The object name is read from the input stream is

  if (objet_existant(name)) //name of an existing object -> read it
    {
      Objet_U& object=objet(name);
      return is >> object; // Then "{ ... }" is read by the object readOn
    }
  else //not the name of an existing object -> read the type, create the object, then read it
    {
      DerObjU ref;
      Nom type;
      is >> type;
      ref.typer(type);
      //add the object to the current interpreter...
      Objet_U& obj = Interprete_bloc::interprete_courant().ajouter(name, ref);
      return is >> obj;           // then read it
    }
}

