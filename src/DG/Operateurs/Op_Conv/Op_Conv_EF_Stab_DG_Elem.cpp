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

#include <Op_Conv_EF_Stab_DG_Elem.h>
#include <Dirichlet_homogene.h>
#include <Masse_ajoutee_base.h>
#include <Schema_Temps_base.h>
#include <Pb_Multiphase.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTLists.h>
#include <Dirichlet.h>
#include <Param.h>
#include <cmath>

Implemente_instanciable( Op_Conv_EF_Stab_DG_Elem, "Op_Conv_EF_Stab_DG_Elem", Op_Conv_DG_base );
Implemente_instanciable( Op_Conv_Amont_DG_Elem, "Op_Conv_Amont_DG_Elem", Op_Conv_EF_Stab_DG_Elem );
Implemente_instanciable( Op_Conv_Centre_DG_Elem, "Op_Conv_Centre_DG_Elem", Op_Conv_EF_Stab_DG_Elem );

// XD Op_Conv_EF_Stab_DG_Elem interprete Op_Conv_EF_Stab_DG_Elem BRACE Class Op_Conv_EF_Stab_DG_Elem

Sortie& Op_Conv_EF_Stab_DG_Elem::printOn(Sortie& os) const { return Op_Conv_DG_base::printOn(os); }
Sortie& Op_Conv_Amont_DG_Elem::printOn(Sortie& os) const { return Op_Conv_DG_base::printOn(os); }
Sortie& Op_Conv_Centre_DG_Elem::printOn(Sortie& os) const { return Op_Conv_DG_base::printOn(os); }

Entree& Op_Conv_EF_Stab_DG_Elem::readOn(Entree& is) { return Op_Conv_DG_base::readOn(is); }

Entree& Op_Conv_Amont_DG_Elem::readOn(Entree& is)
{
  return Op_Conv_DG_base::readOn(is);
}

Entree& Op_Conv_Centre_DG_Elem::readOn(Entree& is)
{
  return Op_Conv_DG_base::readOn(is);
}

void Op_Conv_EF_Stab_DG_Elem::completer()
{
  Op_Conv_DG_base::completer();
}

double Op_Conv_EF_Stab_DG_Elem::calculer_dt_stab() const
{
  double dt = 1e10;
  return Process::mp_min(dt);
}

void Op_Conv_EF_Stab_DG_Elem::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  throw; //TODO DG
}

void Op_Conv_EF_Stab_DG_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  throw; //TODO DG
}

void Op_Conv_EF_Stab_DG_Elem::set_incompressible(const int flag)
{
  if (flag == 0)
    {
      Cerr << "Compressible form of operator \"" << que_suis_je() << "\" :" << finl;
      Cerr << "Discretization of \u2207(inco \u2297 v) - v \u2207.(inco)" << finl;
    }
  incompressible_ = flag;
}
