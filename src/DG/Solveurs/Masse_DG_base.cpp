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

#include <Masse_DG_base.h>
#include <Domaine_Cl_DG.h>
#include <Domaine_DG.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <BasisFunction.h>

Implemente_base(Masse_DG_base, "Masse_DG_base", Solveur_Masse_base);
Sortie& Masse_DG_base::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }
Entree& Masse_DG_base::readOn(Entree& s) { return s; }

void Masse_DG_base::associer_domaine_dis_base(const Domaine_dis_base& le_dom_dis_base)
{
  le_dom_dg_ = ref_cast(Domaine_DG, le_dom_dis_base);
}

void Masse_DG_base::associer_domaine_cl_dis_base(const Domaine_Cl_dis_base& le_dom_Cl_dis_base)
{
  le_dom_Cl_dg_ = ref_cast(Domaine_Cl_DG, le_dom_Cl_dis_base);
}

/**
 * @brief Multiplies the mass coefficient vector element-wise by an optional temporal field.
 *
 * @details If a temporal coefficient has been registered (has_coefficient_temporel_ == true),
 * the field named name_of_coefficient_temporel_ is retrieved from the equation and its
 * values are applied to coef via tab_multiply_any_shape(). Three field types are handled:
 *  - Champ_Inc_base: uses the first part of the field's values (ConstDoubleTab_parts[0]).
 *  - Champ_Fonc_base: uses the field's values directly.
 *  - Champ_Don_base: evaluates the field at the unknown's node coordinates.
 * If no temporal coefficient is registered, coef is left unchanged.
 *
 * @param coef The coefficient vector to scale in-place (typically initialized to 1).
 */
void Masse_DG_base::appliquer_coef(DoubleVect& coef) const
{
  if (has_coefficient_temporel_)
    {
      OBS_PTR(Champ_base) ref_coeff;
      ref_coeff = equation().get_champ(name_of_coefficient_temporel_);

      DoubleTab values;
      if (sub_type(Champ_Inc_base,ref_coeff.valeur()))
        {
          const Champ_Inc_base& coeff = ref_cast(Champ_Inc_base,ref_coeff.valeur());
          ConstDoubleTab_parts val_parts(coeff.valeurs());
          values.ref(val_parts[0]);

        }
      else if (sub_type(Champ_Fonc_base,ref_coeff.valeur()))
        {
          const Champ_Fonc_base& coeff = ref_cast(Champ_Fonc_base,ref_coeff.valeur());
          values.ref(coeff.valeurs());

        }
      else if (sub_type(Champ_Don_base,ref_coeff.valeur()))
        {
          DoubleTab nodes;
          equation().inconnue().remplir_coord_noeuds(nodes);
          ref_coeff->valeur_aux(nodes,values);
        }
      tab_multiply_any_shape(coef, values, VECT_REAL_ITEMS);
    }
}

/**
 * @brief Builds the sparsity pattern of the mass matrix block.
 *
 * @details The pattern depends on whether the basis is orthonormalized:
 *  - **Orthonormal basis**: purely diagonal — one non-zero per DOF
 *    (indice(k,0) == indice(k,1) for every k).
 *  - **Non-orthonormal basis**: full local block — every pair (i,j) within
 *    the same element and same spatial component is coupled.
 *
 * The number of spatial components (dim) is set to Objet_U::dimension for
 * velocity unknowns and 1 for all other fields (scalar).
 *
 * @param matrices  Map of matrix name → Matrice_Morse pointer to be sized.
 * @param semi_impl Unused here; kept for interface compatibility.
 */
void Masse_DG_base::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Nom& nom_inco = equation().inconnue().le_nom();
  if (!matrices.count(nom_inco.getString())) return; //nothing to do

  int order = Option_DG::Get_order_for(nom_inco);
  int dim = nom_inco.debute_par("vitesse") ? Objet_U::dimension : 1;

  const BasisFunction& bfunc = le_dom_dg_->get_basisFunction(order);
  const int nb_bfunc = bfunc.nb_bfunc();

  Matrice_Morse& mat = *matrices.at(nom_inco.getString());
  const int nb_elem_tot = le_dom_dg_->nb_elem_tot();

  int size_indices = le_dom_dg_->gram_schmidt() ? nb_elem_tot*dim*nb_bfunc : nb_elem_tot*dim*nb_bfunc*nb_bfunc;

  IntTab indice(size_indices, 2);
  int current_indice = 0;
  if (le_dom_dg_->gram_schmidt())
    {
      for (int e = 0; e < nb_elem_tot; e++)
        {
          for (int d = 0; d<dim; d++)
            {
              for (int i = 0; i < nb_bfunc; i++ )
                {
                  indice(current_indice+i, 0) = current_indice+i;
                  indice(current_indice+i, 1) = current_indice+i;
                }
              current_indice+=nb_bfunc;
            }
        }
    }
  else
    {
      int index = 0;
      for (int e = 0; e < nb_elem_tot; e++)
        {
          for (int d = 0; d<dim; d++)
            {
              for (int i = 0; i < nb_bfunc; i++ )
                for (int j = 0; j < nb_bfunc; j++ )
                  {
                    indice(index, 0) = current_indice+i;
                    indice(index, 1) = current_indice+j;
                    index++;
                  }
              current_indice+=nb_bfunc;
            }
        }
    }
  mat.dimensionner(indice);
}

/**
 * @brief Assembles the mass matrix contribution (M/dt) into the matrix and right-hand side.
 *
 * @details Computes and accumulates the term (M/dt) * u into the global system,
 * where M is the L2 mass matrix and dt is the time step. The assembly strategy
 * depends on whether the basis is orthonormalized:
 *
 *  - **Orthonormal basis** (gram_schmidt == true):
 *    The mass matrix is diagonal with entries coef[e] * volume[e]. For each DOF:
 *      mat(dof, dof) += coef[e] * volume[e] / dt
 *      secmem(e, dof) += coef[e] * volume[e] * (u^n - delta * u^{n+1}) / dt
 *    where delta = resoudre_en_increments (1 if solving for the increment, 0 otherwise).
 *
 *  - **Non-orthonormal basis, order 0**: reduces to a single scalar per element
 *    (one DOF), equivalent to the cell-average finite volume mass term.
 *
 *  - **Non-orthonormal basis, order > 0**: the full local mass matrix M_ij is
 *    assembled by Gaussian quadrature (Ern, Finite Elements II, 2021, p.71):
 *      M_ij = integral of phi_i * phi_j
 *    accumulated as mat(i,j) += coef[e] * M_ij / dt and the corresponding RHS term.
 *
 * The temporal coefficient (e.g., density) is applied via appliquer_coef() before
 * the loop, and the medium porosity is used as the base coefficient array.
 *
 * @param matrices              Map of matrix name → Matrice_Morse pointer to accumulate into.
 * @param secmem                Right-hand side to accumulate into.
 * @param dt                    Current time step size.
 * @param semi_impl             Map of semi-implicit field values; if the unknown is present,
 *                              its values are used as u^n instead of the stored past values.
 * @param resoudre_en_increments If 1, the system is solved for the increment (u^{n+1} - u^n);
 *                              the current solution is subtracted from the RHS accordingly.
 */
void Masse_DG_base::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, double dt, const tabs_t& semi_impl, int resoudre_en_increments) const
{
  const Nom& nom_inco = equation().inconnue().le_nom();
  const std::string& nom_inco_str = equation().inconnue().le_nom().getString();
  const DoubleTab& passe = semi_impl.count(nom_inco_str) ? semi_impl.at(nom_inco_str) : equation().inconnue().passe();
  const DoubleTab& inco = equation().inconnue().valeurs();
  Matrice_Morse *mat = matrices.count(nom_inco_str) ? matrices.at(nom_inco_str) : nullptr;

  int order = Option_DG::Get_order_for(nom_inco);
  int dim = nom_inco.debute_par("vitesse") ? Objet_U::dimension : 1;

  const BasisFunction& bfunc = le_dom_dg_->get_basisFunction(order);
  const int nb_bfunc = bfunc.nb_bfunc();

  assert(nb_bfunc*dim == inco.line_size());

  const int quad_order = bfunc.get_default_quadrature_order();
  const Quadrature_base& quad = le_dom_dg_->get_quadrature(quad_order);

  const DoubleVect& volume = le_dom_dg_->volumes();
  const int nb_elem_tot = inco.dimension_tot(0);

  int nb_pts_integ_max = quad.nb_pts_integ_max();
  const IntTab& tab_pts_integ= quad.get_tab_nb_pts_integ();

  DoubleTab fbase(nb_bfunc, nb_pts_integ_max);
  DoubleTab product(nb_pts_integ_max);

  DoubleTrav coef(equation().milieu().porosite_elem());
  coef = 1.;
  appliquer_coef(coef);

  int current_indice = 0;
  if (le_dom_dg_->gram_schmidt())
    {
      for (int e = 0; e < nb_elem_tot; e++)
        {
          for (int i=0; i<nb_bfunc; i++)
            {
              for (int d = 0; d<dim; d++)
                {
                  if (mat)
                    (*mat)(current_indice+i+d*nb_bfunc, current_indice+i+d*nb_bfunc) += coef[e]*volume[e] / dt;
                  secmem(e,i+d*nb_bfunc) += coef[e]*volume[e]*(passe(e,i+d*nb_bfunc) - resoudre_en_increments*inco(e,i+d*nb_bfunc))/ dt;
                }
            }
          current_indice+=nb_bfunc*dim;
        }
    }
  else
    {
      for (int e = 0; e < nb_elem_tot; e++)
        {
          if (order == 0)
            for (int d = 0; d<dim; d++)
              {
                if (mat)
                  (*mat)(current_indice+d, current_indice+d) += coef[e]*volume[e] / dt;
                secmem(e,d) += coef[e]*volume[e]*(passe(e,d)-resoudre_en_increments*inco(e,d)) / dt;
                current_indice+=nb_bfunc;
              }
          else
            {
              bfunc.eval_bfunc(quad, e, fbase);
              /****************************************************************/
              /* Quadrature formula: Ern, Finite Elements II, 2021, p 71      */
              /****************************************************************/
              for (int i=0; i<nb_bfunc; i++)
                {
                  for (int j=0; j<nb_bfunc; j++)
                    {
                      product = 0.;
                      for (int k = 0; k < tab_pts_integ(e) ; k++)
                        product(k) = fbase(i, k) * fbase(j, k);

                      double integral = quad.compute_integral_on_elem(e, product);

                      for (int d = 0; d<dim; d++)
                        {
                          if (mat)
                            (*mat)(current_indice+i+d*nb_bfunc, current_indice+j+d*nb_bfunc) += coef[e]*integral / dt;
                          secmem(e,i+d*nb_bfunc) += coef[e]*integral*(passe(e,j+d*nb_bfunc) - resoudre_en_increments*inco(e,j+d*nb_bfunc)) / dt;
                        }
                    }
                }
              current_indice+=nb_bfunc*dim;
            }
        }
    }
}
