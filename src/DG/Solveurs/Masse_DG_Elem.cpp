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

#include <Masse_DG_Elem.h>
#include <Champ_Elem_DG.h>
#include <TRUSTTab_parts.h>
#include <Equation_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <BasisFunction.h>

Implemente_instanciable(Masse_DG_Elem, "Masse_DG_Elem", Masse_DG_base);

Sortie& Masse_DG_Elem::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Masse_DG_Elem::readOn(Entree& s) { return s; }

/**
 * @brief Applies the inverse mass matrix M^{-1} to the right-hand side vector sm in-place.
 *
 * @details Two strategies depending on gram_schmidt:
 *
 *  - **Orthonormal basis** (gram_schmidt == true):
 *    Since M = diag(volume[e]), M^{-1} * sm reduces to dividing each element's
 *    DOF values by the element volume, handled by tab_divide_any_shape().
 *
 *  - **Non-orthonormal basis** (gram_schmidt == false):
 *    For each element, the local inverse mass matrix M^{-1} is computed on-the-fly
 *    by BasisFunction::eval_invMassMatrix(). It is then applied independently to
 *    each spatial component d as a dense matrix-vector product:
 *      sm[e, d*nb_bfunc : (d+1)*nb_bfunc] = M^{-1} * sm[e, d*nb_bfunc : (d+1)*nb_bfunc]
 *    using ref_array slices to avoid data copies.
 *
 *  After inversion, the ghost cell values are synchronized via echange_espace_virtuel().
 *
 * @param sm The right-hand side vector to transform in-place into M^{-1} * sm.
 * @return A reference to sm.
 */
DoubleTab& Masse_DG_Elem::appliquer_impl(DoubleTab& sm) const
{
  if (le_dom_dg_->gram_schmidt())
    {
      const DoubleVect& volume = le_dom_dg_->volumes();

      tab_divide_any_shape(sm, volume);
    }
  else
    {
      const Nom& nom_inco = equation().inconnue().le_nom();
      int order = Option_DG::Get_order_for(nom_inco);
      int dim = nom_inco.debute_par("vitesse") ? Objet_U::dimension : 1;

      const BasisFunction& bfunc = le_dom_dg_->get_basisFunction(order);
      const int nb_bfunc = bfunc.nb_bfunc();

      const int quad_order = bfunc.get_default_quadrature_order();
      const Quadrature_base& quad = le_dom_dg_->get_quadrature(quad_order);

      DoubleTab res, loc;
      DoubleTab invMsm, temp_sm;
      invMsm.copy(sm, RESIZE_OPTIONS::COPY_INIT);
      temp_sm.copy(sm, RESIZE_OPTIONS::COPY_INIT);

      for (int num_elem = 0; num_elem < le_dom_dg_->nb_elem(); num_elem++)
        {
          Matrice_Dense invM = bfunc.eval_invMassMatrix(quad, num_elem);
          for (int d = 0 ; d<dim; d++)
            {
              res.ref_array(temp_sm, num_elem*dim*nb_bfunc + d*nb_bfunc, nb_bfunc);
              loc.ref_array(invMsm, num_elem*dim*nb_bfunc + d*nb_bfunc, nb_bfunc);

              invM.multvect_(loc, res);
            }
        }
      sm.copy(temp_sm, RESIZE_OPTIONS::COPY_INIT);
    }

  sm.echange_espace_virtuel();

  return sm;
}
