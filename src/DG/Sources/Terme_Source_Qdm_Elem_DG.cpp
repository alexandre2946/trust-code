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

#include <Terme_Source_Qdm_Elem_DG.h>
#include <Neumann_sortie_libre.h>
#include <Domaine_Cl_DG.h>
#include <Neumann_homogene.h>
#include <Domaine_DG.h>
#include <Champ_Uniforme.h>

#include <Equation_base.h>
#include <Milieu_base.h>
#include <Neumann.h>

Implemente_instanciable(Terme_Source_Qdm_Elem_DG, "Source_Qdm_Elem_DG", Source_base);

Sortie& Terme_Source_Qdm_Elem_DG::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Terme_Source_Qdm_Elem_DG::readOn(Entree& s)
{
  s >> la_source;
  if (la_source->nb_comp() != dimension)
    {
      Cerr << "Erreur a la lecture du terme source de type " << que_suis_je() << finl;
      Cerr << "le champ source doit avoir " << dimension << " composantes" << finl;
      Process::exit();
    }
  return s ;
}

void Terme_Source_Qdm_Elem_DG::associer_domaines(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_DG = ref_cast(Domaine_DG, domaine_dis);
  le_dom_Cl_DG = ref_cast(Domaine_Cl_DG, domaine_Cl_dis);
}

DoubleTab& Terme_Source_Qdm_Elem_DG::ajouter(DoubleTab& resu) const
{
  if (has_interface_blocs()) return Source_base::ajouter(resu);

  return resu;
}

DoubleTab& Terme_Source_Qdm_Elem_DG::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Terme_Source_Qdm_Elem_DG::mettre_a_jour(double temps)
{
  la_source->mettre_a_jour(temps);
}
