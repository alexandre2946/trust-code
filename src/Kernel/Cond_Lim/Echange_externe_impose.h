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

#ifndef Echange_externe_impose_included
#define Echange_externe_impose_included

#include <Echange_impose_base.h>

/*! @brief Classe Echange_externe_impose: This class represents the special case of the class
 *
 *     Echange_impose_base where the total heat exchange coefficient is computed using
 *     the wall heat exchange coefficient provided by the user.
 *
 *      h_total  : total exchange coefficient
 *      h_imp    : wall exchange coefficient (user input)
 *      e/lambda : internal exchange coefficient
 *                 where: lambda, conductivity in the fluid
 *                        e = d,  in laminar regime
 *                        e = d',  in turbulent regime
 *      (where d' is the equivalent distance computed during application
 *        of the wall laws)
 *           1/h_total = (1/h_imp) + (e/lambda)
 *     It is not the Echange_externe_impose class that manages the computation
 *     of h_total. The class only provides h_imp and it is the diffusive
 *     flux evaluator that computes the e/lambda and h_total terms.
 *     The classes Echange_global_impose and Echange_externe_impose have
 *     exactly the same interface; the Echange_externe_impose class
 *     therefore only serves to signal to the diffusive flux evaluator
 *     that it must compute a total exchange coefficient.
 *
 * @sa Echange_impose_base Echange_global_impose
 */
class Echange_externe_impose: public Echange_impose_base
{
  Declare_instanciable(Echange_externe_impose);
  void verifie_ch_init_nb_comp() const override;
};

#endif
