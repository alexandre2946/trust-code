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

#ifndef Extraire_plan_included
#define Extraire_plan_included

#include <Interprete_geometrique_base.h>

class Nom;
#include <Domaine_forward.h>

/*! @brief Class Extraire_plan Reading a file
 *
 * @sa Interprete
 */
class Extraire_plan : public Interprete_geometrique_base
{
  Declare_instanciable(Extraire_plan);
public :
  Entree& interpreter_(Entree&) override;
};

KOKKOS_INLINE_FUNCTION
void calcul_normal(const double* origine,const double* point1, const double* point2, double* normal)
{
  normal[0]=(point1[1]-origine[1])*(point2[2]-origine[2])-((point2[1]-origine[1])*(point1[2]-origine[2]));
  normal[1]=(point1[2]-origine[2])*(point2[0]-origine[0])-((point2[2]-origine[2])*(point1[0]-origine[0]));
  normal[2]=(point1[0]-origine[0])*(point2[1]-origine[1])-((point2[0]-origine[0])*(point1[1]-origine[1]));
}

#endif /* Extraire_plan_included */
