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

#include <Op_Grad_DG.h>
#include <Check_espace_virtuel.h>
#include <Domaine_Cl_DG.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
#include <EcrFicPartage.h>
#include <Matrice_Morse.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <SFichier.h>
#include <Debog.h>
#include <BasisFunction.h>

Implemente_instanciable(Op_Grad_DG, "Op_Grad_DG", Operateur_Grad_base);

Sortie& Op_Grad_DG::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Op_Grad_DG::readOn(Entree& s) { return s; }

/**
 * @brief Associates the operator with a DG domain and its boundary conditions.
 * @param domaine_dis    The discretized domain, expected to be a Domaine_DG.
 * @param domaine_Cl_dis The boundary condition container, expected to be a Domaine_Cl_DG.
 */
void Op_Grad_DG::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base&)
{
  le_dom_DG = ref_cast(Domaine_DG, domaine_dis);
  le_dcl_DG = ref_cast(Domaine_Cl_DG, domaine_Cl_dis);
}

/**
 * @brief Sizes the velocity-pressure matrix block for the gradient operator.
 *
 * @details Builds the sparsity pattern of the rectangular matrix mapping pressure DOFs
 * to velocity DOFs. Each velocity DOF of element T is coupled to all pressure DOFs of
 * T and its face-neighbours, as given by the pre-computed stencil. The resulting matrix
 * has size_v rows (velocity global DOFs) and size_p columns (pressure global DOFs).
 *
 * Returns immediately if "pression" is absent from matrices or is treated semi-implicitly.
 *
 * @param matrices  Map of matrix name → Matrice_Morse pointer to be sized.
 * @param semi_impl Map of semi-implicit field names; if "pression" is present, returns immediately.
 */
void Op_Grad_DG::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  if (!matrices.count("pression")) return; //rien a faire
  if (semi_impl.count("pression"))
    return; // semi-implicite -> rien a dimensionner

  Matrice_Morse *mat = matrices["pression"], mat2;

  const Domaine_DG& domaine = le_dom_DG.valeur();

  int order_v = Option_DG::Get_order_for("vitesse");
  int order_p = Option_DG::Get_order_for("pression");

  const BasisFunction& bfunc_v = domaine.get_basisFunction(order_v);
  const int nb_bfunc_v = bfunc_v.nb_bfunc();

  const BasisFunction& bfunc_p = domaine.get_basisFunction(order_p);
  const int nb_bfunc_p = bfunc_p.nb_bfunc();

  int dim = Objet_U::dimension;
  const IntTab& indices_glob_elem_v = bfunc_v.indices_glob_elem(dim);
  const IntTab& indices_glob_elem_p = bfunc_p.indices_glob_elem();

  int nb_elem_tot = domaine.nb_elem_tot();

  int size_row = indices_glob_elem_v(nb_elem_tot);
  int size_col = indices_glob_elem_p(nb_elem_tot);

  mat2.dimensionner(size_row, size_col, 0);

  auto& tab1 = mat2.get_set_tab1();
  auto& tab2 = mat2.get_set_tab2();
  auto& coeff = mat2.get_set_coeff();
  coeff = 0;

  const Stencil& stencil_sorted = domaine.get_stencil_sorted();
  const int nb_stencil_max = stencil_sorted.dimension(1);

  int nb_indices_line;
  int row, col, indice;

  tab1(0) = 1;
  for (int nelem = 0; nelem < nb_elem_tot; nelem++)
    {
      nb_indices_line = 0;
      for (int k = 0; k < nb_stencil_max; k++)
        {
          if (stencil_sorted(nelem, k) < 0)
            break;
          nb_indices_line += nb_bfunc_p;
        }
      for (int k = 0; k < nb_bfunc_v * dim; k++)
        tab1(indices_glob_elem_v(nelem) + k + 1) = nb_indices_line + tab1(indices_glob_elem_v(nelem) + k);
    }

  mat2.dimensionner(size_row, size_col, tab1(size_row) - 1);

  for (int nelem = 0; nelem < nb_elem_tot; nelem++)
    {
      row = tab1[indices_glob_elem_v(nelem)] - 1;
      nb_indices_line = tab1[indices_glob_elem_v(nelem) + 1] - tab1[indices_glob_elem_v(nelem)];
      indice = 0;

      for (int i = 0; i < nb_bfunc_v*dim; i++)
        {
          for (int k = 0; k < nb_stencil_max; k++)
            {
              if (stencil_sorted(nelem, k) < 0)
                break;
              col = indices_glob_elem_p(stencil_sorted(nelem, k)) + 1;
              for (int j = 0 ;  j < nb_bfunc_p; j++)
                tab2[row + indice + j + k*nb_bfunc_p] = col + j;
            }
          indice += nb_indices_line;
        }
    }
  mat2.is_sorted_stencil();
  assert(mat2.is_sorted_stencil());
  mat->nb_colonnes() ? *mat += mat2 : *mat = mat2;
}

/**
 * @brief Assembles the DG gradient operator into the matrix and right-hand side.
 *
 * @details The assembly has two stages:
 *
 * **1. Volume term** — for each element:
 *      integral of grad(phi_p_j) . phi_v_i  over the element,
 * accumulated as (*mat)(v_dof, p_dof) += coeff and secmem(elem, v_dof) -= coeff * p(elem, p_dof).
 *
 * **2. Internal face term** — for each internal face shared by elem0 and elem1:
 *   The average-jump coupling integral of {{phi_p}} * [phi_v . n] is expanded into
 *   four pointwise products at face quadrature points:
 *    - coeff00: -0.5 * phi_p0 * n_d * phi_v0  (elem0 pressure, elem0 velocity)
 *    - coeff01: +0.5 * phi_p1 * n_d * phi_v0  (elem1 pressure, elem0 velocity)
 *    - coeff10: -0.5 * phi_p0 * n_d * phi_v1  (elem0 pressure, elem1 velocity)
 *    - coeff11: +0.5 * phi_p1 * n_d * phi_v1  (elem1 pressure, elem1 velocity)
 *
 *   The sign convention follows the outward normal of elem0: the jump [v.n] on the
 *   face is (v0 - v1).n/|f|, hence the minus sign on elem0-side contributions.
 *
 * @param matrices  Map of matrix name → Matrice_Morse pointer to accumulate into.
 * @param secmem    Right-hand side (momentum residual) to accumulate into.
 * @param semi_impl Map of semi-implicit field values; if "pression" is present,
 *                  its values are used instead of the current pressure field.
 */
void Op_Grad_DG::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{

  const DoubleTab& inco_p = semi_impl.count("pression") ? semi_impl.at("pression") :  ref_cast(Navier_Stokes_std, equation()).pression().valeurs(); // NB : is this working ?
  Matrice_Morse *mat = matrices.count("pression") ? matrices.at("pression") : nullptr; // pression for the stabilisation term if np==nv

  const Domaine_DG& domaine = le_dom_DG.valeur();

  int order_v = Option_DG::Get_order_for("vitesse");
  int order_p = Option_DG::Get_order_for("pression");

  const BasisFunction& bfunc_v = domaine.get_basisFunction(order_v);
  const int nb_bfunc_v = bfunc_v.nb_bfunc();

  const BasisFunction& bfunc_p = domaine.get_basisFunction(order_p);
  const int nb_bfunc_p = bfunc_p.nb_bfunc();

  const int dim = Objet_U::dimension;
  const IntTab& indices_glob_elem_v = bfunc_v.indices_glob_elem(dim);
  const IntTab& indices_glob_elem_p = bfunc_p.indices_glob_elem();

  const int quad_order = bfunc_p.get_default_quadrature_order(); //TODO DG should we choose the max or the p order quadrature ?
  const Quadrature_base& quad = domaine.get_quadrature(quad_order);  // Same quadrature for all champs


  int nb_pts_integ_max = quad.nb_pts_integ_max();
  double coeff;
  DoubleTab grad_fbase_elem(nb_bfunc_p, nb_pts_integ_max, Objet_U::dimension);
  DoubleTab f_base_v(nb_bfunc_v, nb_pts_integ_max);
  DoubleTab scalar_product_dim(nb_pts_integ_max);

  // Loop over elements to compute \int q_h div(u_h) dV
  for (int elem = 0; elem < domaine.nb_elem(); elem++)
    {
      int ind_elem_v = indices_glob_elem_v(elem);
      int ind_elem_p = indices_glob_elem_p(elem);
      bfunc_p.eval_grad_bfunc(quad, elem, grad_fbase_elem);
      bfunc_v.eval_bfunc(quad, elem, f_base_v);
      for (int d = 0; d < dim; d++)
        for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
          for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
            {
              for (int k = 0; k < quad.nb_pts_integ(elem); k++)
                scalar_product_dim(k) = grad_fbase_elem(pressure_index, k, d) * f_base_v(velocity_index, k);
              coeff = quad.compute_integral_on_elem(elem, scalar_product_dim);
              if (mat)
                (*mat)(ind_elem_v + velocity_index + d * nb_bfunc_v, ind_elem_p + pressure_index) += coeff;
              secmem(elem, velocity_index + d * nb_bfunc_v) -= coeff * inco_p(elem, pressure_index);
            }
    }

  const DoubleTab& face_normales = domaine.face_normales();
  const DoubleVect& face_surfaces = domaine.face_surfaces();
  int nb_pts_int_fac = quad.nb_pts_integ_facets();
  const IntTab& face_voisins = domaine.face_voisins();

  double coeff00, coeff10, coeff01, coeff11;

  DoubleTab eval_jump_on_facet00(nb_pts_int_fac);
  DoubleTab eval_jump_on_facet01(nb_pts_int_fac);
  DoubleTab eval_jump_on_facet10(nb_pts_int_fac);
  DoubleTab eval_jump_on_facet11(nb_pts_int_fac);

  DoubleTab f_base_v0(nb_bfunc_v, nb_pts_int_fac);
  DoubleTab f_base_v1(nb_bfunc_v, nb_pts_int_fac);
  DoubleTab f_base_p0(nb_bfunc_p, nb_pts_int_fac);
  DoubleTab f_base_p1(nb_bfunc_p, nb_pts_int_fac);

  int premiere_face_int = domaine.premiere_face_int();

  // Loop over facets to compute \int_f [u_h.n]_F {{q_h}} dS
  for (int face = premiere_face_int; face < domaine.nb_faces(); face++)
    {
      int elem0 = face_voisins(face,0);
      int elem1 = face_voisins(face,1);
      double sur_f = face_surfaces(face);
      int ind_elem0_v = indices_glob_elem_v(elem0);
      int ind_elem1_v = indices_glob_elem_v(elem1);
      int ind_elem0_p = indices_glob_elem_p(elem0);
      int ind_elem1_p = indices_glob_elem_p(elem1);
      bfunc_v.eval_bfunc_on_facets(quad, elem0, face, f_base_v0);
      bfunc_v.eval_bfunc_on_facets(quad, elem1, face, f_base_v1);
      bfunc_p.eval_bfunc_on_facets(quad, elem0, face, f_base_p0);
      bfunc_p.eval_bfunc_on_facets(quad, elem1, face, f_base_p1);
      for (int d = 0; d < Objet_U::dimension; d++)
        {
          for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
            {
              for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
                {
                  for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                    {
                      eval_jump_on_facet00(k) = -f_base_p0(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v0(velocity_index, k) / sur_f;
                      eval_jump_on_facet01(k) = +f_base_p1(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v0(velocity_index, k) / sur_f;
                      eval_jump_on_facet10(k) = -f_base_p0(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v1(velocity_index, k) / sur_f;
                      eval_jump_on_facet11(k) = +f_base_p1(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v1(velocity_index, k) / sur_f;
                    }
                  coeff00 = quad.compute_integral_on_facet(face, eval_jump_on_facet00);
                  coeff01 = quad.compute_integral_on_facet(face, eval_jump_on_facet01);
                  coeff10 = quad.compute_integral_on_facet(face, eval_jump_on_facet10);
                  coeff11 = quad.compute_integral_on_facet(face, eval_jump_on_facet11);

                  if (mat)
                    {
                      (*mat)(ind_elem0_v + velocity_index + d * nb_bfunc_v, ind_elem0_p + pressure_index) += coeff00;
                      (*mat)(ind_elem0_v + velocity_index + d * nb_bfunc_v, ind_elem1_p + pressure_index) += coeff01;
                      (*mat)(ind_elem1_v + velocity_index + d * nb_bfunc_v, ind_elem0_p + pressure_index) += coeff10;
                      (*mat)(ind_elem1_v + velocity_index + d * nb_bfunc_v, ind_elem1_p + pressure_index) += coeff11;
                    }
                  secmem(elem0, velocity_index + d * nb_bfunc_v) -= coeff00 * inco_p(elem0, pressure_index);
                  secmem(elem0, velocity_index + d * nb_bfunc_v) -= coeff01 * inco_p(elem1, pressure_index);
                  secmem(elem1, velocity_index + d * nb_bfunc_v) -= coeff10 * inco_p(elem0, pressure_index);
                  secmem(elem1, velocity_index + d * nb_bfunc_v) -= coeff11 * inco_p(elem1, pressure_index);
                }
            }
        }
    }
}

DoubleTab& Op_Grad_DG::calculer(const DoubleTab& pressure, DoubleTab& grad) const
{
  grad = 0.;
  return ajouter(pressure, grad);
}


int Op_Grad_DG::impr(Sortie& os) const
{
  return 0;
}
