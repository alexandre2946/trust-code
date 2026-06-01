/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Interprete_included
#define Interprete_included

#include <Objet_U.h>

class Interprete;

/*! @brief Base class for "interpreter" objects.
 *
 * These objects define actions to be performed when they are
 *    encountered in the data set. The action is triggered by
 *    a call to the interpreter() method.
 *    The interpreter can then read additional parameters
 *    from the input to perform its task. In general, the interpreter
 *    acts on other objects declared in the data set. It can
 *    access them through the objet() and objet_existant() methods.
 *    See for example the Lire or Associer class.
 *
 * @sa Interprete_bloc which reads a series of instructions
 *    to execute in the data set.
 */
class Interprete : public Objet_U
{
  Declare_base(Interprete);
public:
  virtual Entree& interpreter(Entree&) = 0 ;

  static Objet_U& objet(const Nom&);
  static int objet_existant(const Nom&);
};
#endif
