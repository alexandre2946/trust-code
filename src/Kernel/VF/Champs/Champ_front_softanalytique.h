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

#ifndef Champ_front_softanalytique_included
#define Champ_front_softanalytique_included

#include <Ch_front_var_stationnaire.h>
#include <TRUST_Vector.h>
#include <Parser_U.h>

/*! @brief Champ_front_softanalytique Class derived from Champ_front_var representing boundary
 *
 *      fields that are analytical in space and constant in
 *      time.
 *    New boundary field class.
 *    It allows setting a boundary field parameter that depends on a
 *    function. This function is entered in the data file as
 *    a character string. No recompilation needed,
 *    unlike the Champ_front_analytique class which uses a
 *    hard-coded function.
 *
 * @sa Champ_front_base Champ_front_var
 */
class Champ_front_softanalytique : public Ch_front_var_stationnaire
{
  Declare_instanciable(Champ_front_softanalytique);

public:
  int initialiser(double temps, const Champ_Inc_base& inco) override;
  void valeur_a(DoubleVect& position, DoubleVect& valeur);

private :
  VECT(Parser_U) fxyz;

};

#endif

