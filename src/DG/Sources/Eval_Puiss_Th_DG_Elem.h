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

#ifndef Eval_Puiss_Th_DG_Elem_included
#define Eval_Puiss_Th_DG_Elem_included

#include <Evaluateur_Source_Elem.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Champ_Don_base.h>
#include <TRUST_Ref.h>
#include <TRUSTTab.h>
#include <Champ_Elem_DG.h>
#include <Domaine_DG.h>
#include <BasisFunction.h>

/**
 * @brief Evaluator computing the volumetric heat source contribution for DG element unknowns.
 *
 * For each element e, this evaluator computes the projection of the volumetric heat
 * power density Q onto each basis function phi_i:
 *   S(i) = integral of Q(x) * phi_i(x) dV
 * using a fixed order-5 quadrature rule. The power field Q is expected to be
 * pre-evaluated at the quadrature points and stored in the puissance array
 * (shape: nb_elem x nb_pts_integ_max), as set up by associer_champs().
 *
 * calculer_terme_source() is templated on Type_Double to allow both scalar and
 * small-vector output types
 *
 * @sa Terme_Puissance_Thermique_DG_Elem, Evaluateur_Source_Elem
 */
class Eval_Puiss_Th_DG_Elem: public Evaluateur_Source_Elem
{
public:
  void mettre_a_jour() override { }
  void associer_champs(const Champ_Don_base&);

  template <typename Type_Double>
  inline void calculer_terme_source(const int, Type_Double&) const;

protected:
  OBS_PTR(Champ_Don_base) la_puissance;
  DoubleTab puissance;
};

/**
 * @brief Computes the heat source projection onto the local DG basis for element e.
 *
 * @details Evaluates the integral S(fb) = integral of Q * phi_fb over element e
 * using a fixed order-5 quadrature. The integrand at each quadrature point k is:
 *   product(k) = puissance(e, k) * fbase(fb, k)
 * where puissance contains Q pre-sampled at quadrature points and fbase contains
 * the basis function values from eval_bfunc().
 *
 * @tparam Type_Double Output type (typically a fixed-size array of nb_bfunc values).
 * @param e  Element index.
 * @param S  Output array of size nb_bfunc, filled with the projected source values.
 */
template <typename Type_Double>
inline void Eval_Puiss_Th_DG_Elem::calculer_terme_source(const int e, Type_Double& S) const
{
  const Champ_Elem_DG& ch = ref_cast(Champ_Elem_DG, la_zcl->inconnue());
  const Nom& nom_inc = ch.le_nom();

  const Domaine_DG& dom = ref_cast(Domaine_DG, le_dom.valeur());

  int order = Option_DG::Get_order_for(nom_inc);

  const BasisFunction& bfunc = dom.get_basisFunction(order);
  const int nb_bfunc = bfunc.nb_bfunc();

  const Quadrature_base& quad = dom.get_quadrature(5);
  int nb_pts_integ_max = quad.nb_pts_integ_max();

  DoubleTab product(nb_pts_integ_max);

  DoubleTab fbase(nb_bfunc, nb_pts_integ_max);
  bfunc.eval_bfunc(quad, e, fbase);

  for (int fb = 0; fb < nb_bfunc; fb++)
    {
      for (int k = 0; k < quad.nb_pts_integ(e) ; k++)
        product(k) = puissance(e,k) * fbase(fb, k);

      S(fb) = quad.compute_integral_on_elem(e, product);
    }
}

#endif /* Eval_Puiss_Th_DG_Elem_included */
