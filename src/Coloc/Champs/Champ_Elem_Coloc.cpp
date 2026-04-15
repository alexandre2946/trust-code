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

#include <Champ_Elem_Coloc.h>
#include <Domaine_Coloc.h>
#include <Pb_Euler.h>

Implemente_instanciable(Champ_Elem_Coloc, "Champ_Elem_Coloc", Champ_Inc_P0_base);

Sortie& Champ_Elem_Coloc::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Champ_Elem_Coloc::readOn(Entree& s)
{
  lire_donnees(s);
  return s;
}

int Champ_Elem_Coloc::nb_valeurs_nodales() const
{
  return domaine_Coloc().nb_elem();
}

int Champ_Elem_Coloc::fixer_nb_valeurs_nodales(int n)
{
  assert (n == domaine_dis_base().domaine().nb_elem());
  creer_tableau_distribue(domaine_dis_base().domaine().md_vector_elements());
  return n;
}

const Domaine_Coloc& Champ_Elem_Coloc::domaine_Coloc() const
{
  return ref_cast(Domaine_Coloc, le_dom_VF.valeur());
}

int Champ_Elem_Coloc::reprendre(Entree& fich)
{
  Process::exit("Champ_Elem_Coloc::reprendre not coded ... \n");
  return 1;
}
