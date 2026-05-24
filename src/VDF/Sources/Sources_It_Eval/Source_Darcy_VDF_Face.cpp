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

#include <Source_Darcy_VDF_Face.h>
#include <Fluide_Incompressible.h>
#include <Probleme_base.h>
#include <Domaine_Cl_VDF.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Source_Darcy_VDF_Face,"Darcy_VDF_Face",Terme_Source_VDF_base);

// Hand-maintained XD: the dataset keyword 'darcy' diverges from the C++
// 'Darcy_VDF_Face' name above, and its declared parent 'source_base' diverges
// from the actual 'Terme_Source_VDF_base' — preserved verbatim from the
// legacy TRAD_2.org to keep current schema behavior. Body is parsed opaquely
// via bloc_lecture so the divergence is invisible to the parser.
// XD darcy source_base darcy NO_BRACE Class for calculation in a porous media with source term of Darcy -nu/K*V. This
// XD_CONT keyword must be used with a permeability model. For the moment there are two models : permeability constant
// XD_CONT or Ergun's law. Darcy source term is available for quasi compressible calculation. A new keyword is aded for
// XD_CONT porosity (porosite).
// XD attr bloc bloc_lecture bloc REQ Description.

Sortie& Source_Darcy_VDF_Face::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Source_Darcy_VDF_Face::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  set_fichier("Source_Darcy");
  set_description("Darcy term = Integral(-nu/K*vitesse*dv) [m/s2]");
  return is;
}

void Source_Darcy_VDF_Face::set_param(Param& param) const
{
  param.ajouter_non_std("modele_K",(this),Param::REQUIRED);
  param.ajouter_non_std("porosite",(this),Param::REQUIRED);
}

int Source_Darcy_VDF_Face::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="modele_K")
    {
      Motcle motlu;
      is >> motlu;
      eval().modK_.typer(motlu);
      is >> eval().modK_.valeur();
      return 1;
    }
  else if (mot=="porosite")
    {
      is >> eval().getPorosite();
      return 1;
    }
  else
    {
      Cerr << "Unknown keyword in Source_Darcy_VDF_Face" << finl;
      Process::exit();
    }
  return -1;
}

void Source_Darcy_VDF_Face::associer_domaines(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_cl_dis)
{
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_dis);
  const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF,domaine_cl_dis);
  iter_->associer_domaines(zvdf, zclvdf);
  eval().associer_domaines(zvdf, zclvdf);
}

void Source_Darcy_VDF_Face::associer_pb(const Probleme_base& pb)
{
  const Champ_Don_base& diffu = ref_cast(Fluide_base,pb.milieu()).viscosite_cinematique();
  const Champ_Inc_base& vit = pb.equation(0).inconnue();
  eval().associer(diffu);
  eval().associer(vit);
}
