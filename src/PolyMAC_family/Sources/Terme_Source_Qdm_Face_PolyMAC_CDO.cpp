/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Terme_Source_Qdm_Face_PolyMAC_CDO.h>
#include <Neumann_sortie_libre.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Neumann_homogene.h>
#include <Domaine_PolyMAC_CDO.h>
#include <Champ_Uniforme.h>

#include <Equation_base.h>
#include <Milieu_base.h>
#include <Neumann.h>

Implemente_instanciable(Terme_Source_Qdm_Face_PolyMAC_CDO, "Source_Qdm_face_PolyMAC_CDO", Source_base);

Sortie& Terme_Source_Qdm_Face_PolyMAC_CDO::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Terme_Source_Qdm_Face_PolyMAC_CDO::readOn(Entree& s)
{
  s >> la_source;
  if (la_source->nb_comp() != equation().inconnue().nb_comp())
    {
      Cerr << "Error reading the source term of type " << que_suis_je() << finl;
      Cerr << "the source field must have " << dimension << " components" << finl;
      Process::exit();
    }
  return s ;
}

void Terme_Source_Qdm_Face_PolyMAC_CDO::associer_domaines(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_PolyMAC_CDO = ref_cast(Domaine_PolyMAC_CDO, domaine_dis);
  le_dom_Cl_PolyMAC_CDO = ref_cast(Domaine_Cl_PolyMAC_family, domaine_Cl_dis);
}

DoubleTab& Terme_Source_Qdm_Face_PolyMAC_CDO::ajouter(DoubleTab& resu) const
{
  if (has_interface_blocs()) return Source_base::ajouter(resu);

  const Domaine_PolyMAC_CDO& domaine_PolyMAC_CDO = le_dom_PolyMAC_CDO.valeur();
  const Domaine_Cl_PolyMAC_family& Domaine_Cl_PolyMAC_family = le_dom_Cl_PolyMAC_CDO.valeur();

  /* 1. boundary faces: contribute only to Neumann faces */
  for (int n_bord = 0; n_bord < domaine_PolyMAC_CDO.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = Domaine_Cl_PolyMAC_family.les_conditions_limites(n_bord);
      if (!sub_type(Neumann, la_cl.valeur()) && !sub_type(Neumann_homogene, la_cl.valeur()) && !sub_type(Neumann_val_ext, la_cl.valeur()))
        continue;
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      for (int f = le_bord.num_premiere_face(); f < le_bord.num_premiere_face() + le_bord.nb_faces(); f++)
        {
          int e = domaine_PolyMAC_CDO.face_voisins(f, 0);
          double fac = equation().milieu().porosite_face(f) * domaine_PolyMAC_CDO.face_surfaces(f);
          for (int r = 0; r < dimension; r++)
            resu(f) += fac * la_source->valeurs()(sub_type(Champ_Uniforme,la_source.valeur()) ? 0 : e, r) * (domaine_PolyMAC_CDO.xv(f, r) - domaine_PolyMAC_CDO.xp(e, r));
        }
    }
  /* 2. internal faces -> upstream/downstream contributions */
  for (int f = domaine_PolyMAC_CDO.premiere_face_int(); f < domaine_PolyMAC_CDO.nb_faces(); f++)
    {
      double fac = equation().milieu().porosite_face(f) * domaine_PolyMAC_CDO.face_surfaces(f);
      for (int i = 0; i < 2; i++)
        for (int r = 0, e = domaine_PolyMAC_CDO.face_voisins(f, i); r < dimension; r++)
          resu(f) += fac * la_source->valeurs()(sub_type(Champ_Uniforme,la_source.valeur()) ? 0 : e, r) * (domaine_PolyMAC_CDO.xv(f, r) - domaine_PolyMAC_CDO.xp(e, r)) * (i ? -1 : 1);
    }
  return resu;
}

DoubleTab& Terme_Source_Qdm_Face_PolyMAC_CDO::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Terme_Source_Qdm_Face_PolyMAC_CDO::mettre_a_jour(double temps)
{
  la_source->mettre_a_jour(temps);
}
