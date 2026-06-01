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

#ifndef Tetraedriser_par_prisme_included
#define Tetraedriser_par_prisme_included

#include <Triangulation_base.h>
#include <Domaine_forward.h>

/*! @brief Class Tetra_par_prisme This class is an interpreter that serves to read and execute
 *
 *     the Tetra_par_prisme directive:
 *         Tetra_par_prisme domain_name
 *     This directive is to be used in VEF discretization to obtain
 *     a tetrahedral mesh (via prisms) from a mesh made up of blocks.
 *
 * @sa Interprete Pave Tetraedre, This class is usable in 3D
 */
class Tetraedriser_par_prisme : public Triangulation_base
{
  Declare_instanciable(Tetraedriser_par_prisme);

public :

  void trianguler(Domaine&) const override;
  inline int dimension_application() const override;
};

inline int Tetraedriser_par_prisme::dimension_application() const
{
  return 3;
}

#endif
