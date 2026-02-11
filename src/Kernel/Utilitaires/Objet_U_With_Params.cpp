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

#include <Objet_U_With_Params.h>
#include <Param.h>


Implemente_base(Objet_U_With_Params,"Objet_U_With_Params",Objet_U);

/*! @brief Called in the readOn of ObjetUWithParams. Must be overriden by user.
 *
 * In most cases, it should call the method of the base class to inherit the params.
 * This is not always the case. If not done, reasons must be explained clearly in a comment
 *
 * Should fill the param argument with parameters to be read, using method ajouter
 *
 * Is const, because the modifications are done later in the readOn, with param.lire_avec_accolades_depuis
 *
 */
void Objet_U_With_Params::set_param(Param& param) const
{
}

/*! @brief Reading of the input stream (dataset) with the helper class Param
 *
 * This method is 'final'. The goal is to enforce a specific flow (see details below)
 * of the reading/initialization of TRUST Objects.
 *
 *
 * Unfolding of the method:
 *  1. param is filled by the const method set_param.
 *  2. params are read with param.lire_avec_accolades_depuis
 *  3. params are validated by the user in the const method validate_params
 *
 *
 * @param (Entree& is) an input stream from which param will be read
 * @return (Entree&) the input stream, after reading a block { ... }
 * @throws wrong input format, missing brackets
 */
Entree& Objet_U_With_Params::readOn(Entree& is) // final, see header
{
  return __readOn_Impl_DO_NOT_OVERRIDE_UNLESS_THE_WORLD_IS_ENDING(is);
}

Sortie& Objet_U_With_Params::printOn(Sortie& os) const
{
  return os;
}

/*! @brief FORBIDDEN KNOWLEDGE
 *
 * Provided if absolutely necessary to change the impl of a readOn.
 *
 */
Entree& Objet_U_With_Params::__readOn_Impl_DO_NOT_OVERRIDE_UNLESS_THE_WORLD_IS_ENDING(Entree& is)
{

  Cerr<<"Reading of data for a "<<que_suis_je()<<" equation"<<finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  validate_params();
  return is;

}
