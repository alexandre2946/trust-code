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

#include <Operateur_NConserv.h>

Implemente_instanciable(Operateur_NConserv, "Operateur_NConserv", OWN_PTR(Operateur_NConserv_base));
// XD op_non_conservativtifs_deriv objet_lecture op_non_conservativtifs_deriv NO_BRACE not_set


// XD termes_non_conservatifs_negligeable op_non_conservativtifs_deriv negligeable NO_BRACE For Coloc discretization. Suppresses the non_conservative operator.

// XD termes_non_conservatifs_hll op_non_conservativtifs_deriv hll NO_BRACE Keyword for Coloc discretization. Activates the HLL non-conservative scheme.

// XD convection_hll convection_deriv hll NO_BRACE Keyword for Coloc discretization. Activates the HLL conservative scheme.

// XD convection_rusanov convection_deriv rusanov NO_BRACE Keyword for Coloc discretization. Activates the Rusanov conservative scheme.

Sortie& Operateur_NConserv::printOn(Sortie& os) const { return Operateur::ecrire(os); }

Entree& Operateur_NConserv::readOn(Entree& is)
{
  Operateur::lire(is);
  return is;
}

void Operateur_NConserv::typer()
{
  if (Motcle(typ) == Motcle("negligeable"))
    {
      OWN_PTR(Operateur_NConserv_base)::typer("Op_NConserv_negligeable");
    }
  else
    {
      Equation_base& eqn = mon_equation.valeur();
      Nom nom_type = eqn.discretisation().get_name_of_type_for(que_suis_je(), typ, eqn);
      OWN_PTR(Operateur_NConserv_base)::typer(nom_type);
    }
  Cerr << valeur().que_suis_je() << finl;
}

DoubleTab& Operateur_NConserv::ajouter(const DoubleTab& donnee, DoubleTab& resu) const
{
  DoubleTab& tmp = valeur().ajouter(donnee, resu);
  return tmp;
}

DoubleTab& Operateur_NConserv::calculer(const DoubleTab& donnee, DoubleTab& resu) const
{
  DoubleTab& tmp = valeur().calculer(donnee, resu);
  return tmp;
}

