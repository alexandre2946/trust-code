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

#ifndef Mailler_included
#define Mailler_included

#include <Interprete_geometrique_base.h>

/*! @brief Class Mailler A mesher by agglomeration of domains (paves in 2D and blocks in 3D).
 *
 *     Data set structure (in dimension 2):
 *     Mailler dom
 *     {
 *     [Epsilon eps]
 *     Pave nompave1
 *     {
 *        Origine OX OY
 *           Longueurs LX LY
 *           Nombre_de_Noeuds NX NY
 *     }
 *     {
 *        Bord nom_bord1 X = X0 Y0 <= Y <= Y1
 *        Bord nom_bord2 X = X0 Y1 <= Y <= Y2
 *        Bord nom_bord2 Y = Y0 X1 <= X <= X2
 *        ...
 *     } ,
 *     Pave nompave2 ...
 *     }
 *     Two points will be merged as soon as the distance between them is
 *     less than Epsilon.
 *
 * @sa Interprete Pave, Currently the only object type recognized by TRUST for meshing a domain is Pave
 */
template <typename _SIZE_>
class Mailler_32_64 : public Interprete_geometrique_base_32_64<_SIZE_>
{
  Declare_instanciable_32_64(Mailler_32_64);
public :
  using Domaine_t = Domaine_32_64<_SIZE_>;

  Entree& interpreter_(Entree&) override;
};
#endif
