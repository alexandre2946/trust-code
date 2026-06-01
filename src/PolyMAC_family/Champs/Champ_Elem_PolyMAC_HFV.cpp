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

#include <Champ_Elem_PolyMAC_HFV.h>
#include <Domaine_PolyMAC_HFV.h>
#include <Connectivite_som_elem.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <MD_Vector_base.h>

#include <Pb_Multiphase.h>

Implemente_instanciable(Champ_Elem_PolyMAC_HFV,"Champ_Elem_PolyMAC_HFV",Champ_Elem_PolyMAC_CDO);

Sortie& Champ_Elem_PolyMAC_HFV::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Champ_Elem_PolyMAC_HFV::readOn(Entree& s) { return Champ_Elem_PolyMAC_CDO::readOn(s) ; }

const Domaine_PolyMAC_HFV& Champ_Elem_PolyMAC_HFV::domaine_PolyMAC_HFV() const
{
  return ref_cast(Domaine_PolyMAC_HFV, le_dom_VF.valeur());
}

int Champ_Elem_PolyMAC_HFV::fixer_nb_valeurs_nodales(int n)
{
  assert (n == domaine_dis_base().domaine().nb_elem() || n < 0); //accept both VEF and VDF conventions
  creer_tableau_distribue(domaine_dis_base().domaine().md_vector_elements());
  return n;
}

int Champ_Elem_PolyMAC_HFV::nb_valeurs_nodales() const
{
  return le_dom_VF->nb_elem(); //ignore auxiliary variables
}

void Champ_Elem_PolyMAC_HFV::init_auxiliary_variables()
{
  const Domaine_PolyMAC_HFV& domaine = domaine_PolyMAC_HFV();
  const IntTab& f_e = domaine.face_voisins();
  for (int n = 0; n < nb_valeurs_temporelles(); n++)
    if (futur(n).size_reelle_ok())
      {
        DoubleTab& vals = futur(n);
        vals.set_md_vector(MD_Vector()); //remove the MD_Vector...
        vals.resize_dim0(domaine.mdv_elems_faces->get_nb_items_tot()); //...resize to the correct size...
        vals.set_md_vector(domaine.mdv_elems_faces); //...and restore the correct MD_Vector
        /* initialize face variables: using the value of the upstream element */
        for (int f = 0, ne_tot = domaine.nb_elem_tot(); f < domaine.nb_faces(); f++)
          for (int m = 0, e = f_e(f, 0); m < vals.dimension(1); m++)
            vals(ne_tot + f, m) = vals(e, m);
        vals.echange_espace_virtuel();
      }
}

int Champ_Elem_PolyMAC_HFV::reprendre(Entree& fich)
{
  if (! via_ch_fonc_reprise()) return Champ_Inc_base::reprendre(fich); /* ie: resume last time ! */

  const Pb_Multiphase * pbm = mon_equation_non_nul() ? (sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : nullptr) : nullptr;
  if (pbm) return Champ_Inc_base::reprendre(fich);

  // otherwise do this ...
  init_auxiliary_variables();
  return Champ_Inc_base::reprendre(fich);
}

Champ_base& Champ_Elem_PolyMAC_HFV::affecter_(const Champ_base& ch)
{
  const Domaine_PolyMAC_HFV& domaine = domaine_PolyMAC_HFV();
  if (domaine.nb_elem_tot() > 0 && ch.valeurs().dimension_tot(0) > domaine.nb_elem_tot())
    init_auxiliary_variables(), valeurs() = ch.valeurs();
  else
    {
      DoubleTab_parts part(valeurs());
      for (int i = 0; i < part.size(); i++) ch.valeur_aux(i ? domaine.xv() : domaine.xp(), part[i]);
      valeurs().echange_espace_virtuel();
    }
  return *this;
}
