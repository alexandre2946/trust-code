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
#include <BasisFunction.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Neumann.h>
#include <Field_base.h>

Implemente_instanciable(Terme_Source_Qdm_Elem_DG, "Source_Qdm_Elem_DG", Source_base);

Sortie& Terme_Source_Qdm_Elem_DG::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Terme_Source_Qdm_Elem_DG::readOn(Entree& s)
{
  s >> la_source;
  const int nb_comp = la_source->nb_comp();

  Noms noms;
  Noms unites;
  noms.add("Terme_Source_DG");
  unites.add("kg.m.s^-1");

  equation().probleme().discretisation().discretiser_champ("champ_fonc_quad_dg", equation().domaine_dis(), Nature_du_champ::vectoriel, noms, unites,nb_comp,0., la_source_DG);
  la_source_DG->valeurs() = 0.;
  la_source_DG->affecter(la_source);

  return s ;
}

void Terme_Source_Qdm_Elem_DG::associer_domaines(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_DG = ref_cast(Domaine_DG, domaine_dis);
  le_dom_Cl_DG = ref_cast(Domaine_Cl_DG, domaine_Cl_dis);
}

void Terme_Source_Qdm_Elem_DG::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_DG& dom = ref_cast(Domaine_DG, le_dom_DG.valeur());

  int order = Option_DG::Get_order_for("vitesse");

  const BasisFunction& bfunc = dom.get_basisFunction(order);
  const int nb_bfunc = bfunc.nb_bfunc();

  const Quadrature_base& quad = dom.get_quadrature(5);
  int nb_pts_integ_max = quad.nb_pts_integ_max();

  DoubleTab product(nb_pts_integ_max);

  DoubleTab fbase(nb_bfunc, nb_pts_integ_max);

  const int dim = Objet_U::dimension;

  for (int elem = 0; elem < dom.nb_elem(); elem++)
    {
      bfunc.eval_bfunc(quad, elem, fbase);
      for (int d = 0; d < dim; d++)
        for (int fb = 0; fb < nb_bfunc; fb++)
          {
            for (int k = 0; k < quad.nb_pts_integ(elem) ; k++)
              product(k) =  la_source_DG->valeurs()(sub_type(Champ_Uniforme,la_source.valeur()) ? 0 : elem, k + d*nb_pts_integ_max) * fbase(fb, k);

            secmem(elem, fb + d*nb_bfunc) += quad.compute_integral_on_elem(elem, product);
          }
    }
}

void Terme_Source_Qdm_Elem_DG::mettre_a_jour(double temps)
{
  la_source->mettre_a_jour(temps);
}
