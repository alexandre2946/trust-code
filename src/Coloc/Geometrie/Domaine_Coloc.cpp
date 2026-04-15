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

#include <Poly_geom_base.h>
#include <Domaine_Coloc.h>

Implemente_instanciable(Domaine_Coloc, "Domaine_Coloc", Domaine_Poly_base);

Sortie& Domaine_Coloc::printOn(Sortie& os) const { return Domaine_Poly_base::printOn(os); }

Entree& Domaine_Coloc::readOn(Entree& is) { return Domaine_Poly_base::readOn(is); }

void Domaine_Coloc::discretiser()
{
  Domaine_Poly_base::discretiser();
  calculer_h_carre();
// discretiser_aretes();
}

void Domaine_Coloc::calculer_h_carre()
{
  // Calcul de h_carre
  h_carre = 1.e30;
  h_carre_.resize(nb_faces());
  // Calcul des surfaces
  Elem_geom_base& elem_geom = domaine().type_elem().valeur();
  int is_polyedre = sub_type(Poly_geom_base, elem_geom) ? 1 : 0;
  const ArrOfInt PolyIndex = is_polyedre ? ref_cast(Poly_geom_base, domaine().type_elem().valeur()).getElemIndex() : ArrOfInt(0);
  const DoubleVect& surfaces=face_surfaces();
  const int nbe=nb_elem();
  for (int num_elem=0; num_elem<nbe; num_elem++)
    {
      double surf_max = 0;
      const int nb_faces_elem = is_polyedre ? PolyIndex[num_elem+1] - PolyIndex[num_elem] : domaine().nb_faces_elem();
      for (int i=0; i<nb_faces_elem; i++)
        {
          double surf = surfaces(elem_faces(num_elem,i));
          surf_max = (surf > surf_max)? surf : surf_max;
        }
      double vol = volumes(num_elem)/surf_max;
      vol *= vol;
      h_carre_(num_elem) = vol;
      h_carre = ( vol < h_carre )? vol : h_carre;
    }
}

void Domaine_Coloc::modifier_pour_Cl(const Conds_lim& conds_lim)
{
  Cerr << "Le Domaine_Coloc a ete rempli avec succes" << finl;
  for (int num_cond_lim=0; num_cond_lim<conds_lim.size(); num_cond_lim++)
    {
      const Cond_lim_base& cl = conds_lim[num_cond_lim].valeur();
      if (sub_type(Periodique, cl))
        Process::exit("Domaine_Coloc::modifier_pour_Cl not coded for Periodic BC \n");
    }

  // PQ : 10/10/05 : les faces periodiques etant a double contribution
  //          l'appel a marquer_faces_double_contrib s'effectue dans cette methode
  //          afin de pouvoir beneficier de conds_lim.
  Domaine_VF::marquer_faces_double_contrib(conds_lim);
}
