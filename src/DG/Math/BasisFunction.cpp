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

#include <BasisFunction.h>
#include <TRUSTTab_parts.h>
#include <Matrix_tools.h>
#include <Array_tools.h>

/**
 * @brief Initializes the basis function object for a given DG domain and polynomial order.
 *
 * @details Sets up the global DOF index table indices_glob_elem_, selects the default
 * quadrature order (3 for order 1, 5 for order 2), and, if orthonormalization is
 * requested, allocates and builds the block-diagonal transition matrix via
 * allocate_transition_matrix() and build_transition_matrix(). Finally computes the
 * element and face stabilization parameters via compute_stab_param().
 *
 * @param dom          The DG domain providing mesh geometry and connectivity.
 * @param order        Polynomial order of the basis (0, 1, or 2).
 * @param gram_schmidt If true, the basis is L2-orthonormalized via Gram-Schmidt.
 */
void BasisFunction::initialize(const Domaine_DG& dom, const int& order, const bool& gram_schmidt)
{
  dom_ = dom;
  order_ = order;
  nb_bfunc_ = Option_DG::Nb_col_from_order(order_);
  is_orthonormalized_ = false;
  is_diagonal_ = dom.gram_schmidt();

  int nb_elem_tot = dom_->nb_elem_tot();
  indices_glob_elem_.resize(nb_elem_tot+1);
  indices_glob_elem_(0)=0;
  for (int e = 0; e < nb_elem_tot; e++)
    indices_glob_elem_(e+1) = indices_glob_elem_(e) +  nb_bfunc_;
  default_quad_order_ = 3*(order_==1)+5*(order_==2);

  if (is_diagonal_)
    {
      allocate_transition_matrix();
      build_transition_matrix();
    }



  // computation of the stabilization parameters
  compute_stab_param();
}

/**
 * @brief Allocates the sparsity pattern of the block-diagonal transition matrix.
 *
 * @details Builds a stencil with a full nb_bfunc x nb_bfunc dense block for each
 * element, then calls Matrix_tools::allocate_morse_matrix() to size transition_matrix_.
 * The transition matrix is later filled by build_transition_matrix() with the
 * Gram-Schmidt change-of-basis coefficients.
 */
void BasisFunction::allocate_transition_matrix()
{
  Stencil indice(0, 2);
  int nb_elem_tot = dom_->nb_elem_tot();

  int current_indice = 0;
  for (int e = 0; e < nb_elem_tot; e++)
    {
      for (int i = 0; i < nb_bfunc_; i++ )
        for (int j = 0; j < nb_bfunc_; j++ )
          indice.append_line( current_indice+i, current_indice+j);
      current_indice+=nb_bfunc_;
    }

  int size_inc = indices_glob_elem_(nb_elem_tot);

  tableau_trier_retirer_doublons(indice);
  Matrix_tools::allocate_morse_matrix(size_inc, size_inc, indice, transition_matrix_);
}

/**
 * @brief Computes and stores the Gram-Schmidt orthonormalization coefficients.
 *
 * @details Iterates over all elements, evaluates the raw monomial basis at the element
 * quadrature points, then calls gramSchmidt() to orthonormalize the basis in-place and
 * record the change-of-basis coefficients in transition_matrix_. Sets is_orthonormalized_
 * to true upon completion.
 */
void BasisFunction::build_transition_matrix()
{
  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());

  const DoubleVect& ve = domaine.volumes();
  int current_indice = 0;
  const Quadrature_base& quad = domaine.get_quadrature(default_quad_order_);
  const IntTab& tab_pts_integ= quad.get_tab_nb_pts_integ();
  int nb_pts_integ_max = quad.nb_pts_integ_max();
  DoubleTab fbase(nb_bfunc_, nb_pts_integ_max);

  for (int e = 0; e < domaine.nb_elem_tot(); e++)
    {
      int nb_pts_integ = tab_pts_integ(e);

      eval_bfunc(quad, e, fbase);
      gramSchmidt(fbase, quad, e, current_indice, nb_pts_integ, ve(e), 0);

      current_indice+=nb_bfunc_;
    }
  is_orthonormalized_ = true;

}

/**
 * @brief Recursively orthonormalizes the local basis using the modified Gram-Schmidt process.
 *
 * @details At step `index`, the function:
 *  1. Projects fbase[index] onto all previously orthonormalized vectors fbase[0..index-1]
 *     and subtracts those projections, making fbase[index] orthogonal to all previous ones.
 *  2. Normalizes fbase[index] by its L2 norm (integral divided by element volume).
 *  3. Records each operation as a row in transition_matrix_ so that orthonormalize()
 *     can re-apply the same transform to any future raw basis evaluation without
 *     repeating the quadrature.
 *  4. Recurses to process index+1.
 *
 * @param fbase           Raw basis values at quadrature points (nb_bfunc x nb_pts_integ),
 *                        modified in-place to become orthonormal.
 * @param quad            Quadrature rule used to compute inner products.
 * @param num_elem        Global element index (used to query the quadrature).
 * @param current_indice  First global DOF index of this element in transition_matrix_.
 * @param nb_pts_integ    Number of active quadrature points for this element.
 * @param volume          Volume of the element, used for normalization.
 * @param index           Current basis function index being orthonormalized (0-based).
 */
void BasisFunction::gramSchmidt(DoubleTab& fbase, const Quadrature_base& quad, const int& num_elem, const int& current_indice, const int& nb_pts_integ, const double& volume, int index)
{
  if (index >= nb_bfunc_) return;

  DoubleTab product(nb_pts_integ);
  transition_matrix_(current_indice+index, current_indice+index) = 1.;

  for (int j = 0; j < index; ++j)
    {
      for (int k = 0; k < nb_pts_integ ; k++)
        product(k) = fbase(index, k) * fbase(j, k);

      double numerateur = quad.compute_integral_on_elem(num_elem, product);

      for (int k = 0; k < nb_pts_integ ; k++)
        product(k) = fbase(j, k) * fbase(j, k);

      double denominateur = quad.compute_integral_on_elem(num_elem, product);

      double projection = numerateur / denominateur;

      for (int k = 0; k < nb_pts_integ ; k++)
        fbase(index,k) -= projection * fbase(j,k);
      for (int l = 0; l < j+1; l++)
        transition_matrix_(current_indice+index, current_indice+l) -= projection*transition_matrix_(current_indice+j, current_indice+l);
    }

  // Normalize the basis
  for (int k = 0; k < nb_pts_integ ; k++)
    product(k) = fbase(index, k) * fbase(index, k);
  double norm = quad.compute_integral_on_elem(num_elem, product)/volume;
  for (int j = 0; j < index+1 ; j++)
    transition_matrix_(current_indice+index, current_indice+j) /= sqrt(norm);
  for (int k = 0; k < nb_pts_integ ; k++)
    fbase(index, k) /= sqrt(norm);

  // Recursive for the next index
  gramSchmidt(fbase, quad, num_elem, current_indice, nb_pts_integ, volume, index + 1);
}

/**
 * @brief Computes the local L2 mass matrix M_ij = integral of phi_i * phi_j on element nelem.
 *
 * @param quad   Quadrature rule used for integration.
 * @param nelem  Index of the element.
 * @return A nb_bfunc x nb_bfunc dense matrix containing the local mass matrix.
 */
const Matrice_Dense BasisFunction::build_local_mass_matrix(const Quadrature_base& quad, const int nelem) const
{
  int nb_pts_integ_max = quad.nb_pts_integ_max();
  const IntTab& tab_pts_integ= quad.get_tab_nb_pts_integ();
  DoubleTab fbase(nb_bfunc_, nb_pts_integ_max);
  DoubleTab product(nb_pts_integ_max);
  Matrice_Dense loc_mass_mat(nb_bfunc_, nb_bfunc_);

  eval_bfunc(quad, nelem, fbase);
  for (int i=0; i<nb_bfunc_; i++)
    {
      for (int j=0; j<nb_bfunc_; j++)
        {
          product = 0.;
          for (int k = 0; k < tab_pts_integ(nelem) ; k++)
            product(k) = fbase(i, k) * fbase(j, k);

          loc_mass_mat(i, j) = quad.compute_integral_on_elem(nelem, product);
        }
    }
  return loc_mass_mat;
}

/**
 * @brief Evaluates all basis functions at the element quadrature points.
 *
 * @details Fills fbasis(i, k) with the value of the i-th basis function at the k-th
 * quadrature point of element nelem, using the scaled monomial basis centered at the
 * element barycenter. If the basis is orthonormalized, the transition matrix is applied
 * via orthonormalize().
 *
 * @param quad    Quadrature rule providing integration point coordinates.
 * @param nelem   Element index.
 * @param fbasis  Output array of shape (nb_bfunc, nb_pts_integ_max), filled in-place.
 *
 * @note Only 2D and orders 0-2 are implemented. Throws for order > 2 or 3D.
 */
void BasisFunction::eval_bfunc(const Quadrature_base& quad, const int& nelem, DoubleTab& fbasis) const
{
  const DoubleTab& integ_points = quad.get_integ_points();

  assert(fbasis.dimension(0) == nb_bfunc_ && fbasis.dimension(1) == quad.nb_pts_integ_max());


  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());
  const DoubleTab& xp = domaine.xp(); // barycentre elem

  if (Objet_U::dimension == 2)
    {
      double invh = 1./sqrt(domaine.carre_pas_maille(nelem));

      for (int pt = 0; pt < quad.nb_pts_integ(nelem); pt++)
        {
          fbasis(0, pt) = 1;
          if (order_ == 0) continue;

          fbasis(1, pt) = (integ_points(quad.ind_pts_integ(nelem) + pt, 0) - xp(nelem,0))*invh;
          fbasis(2, pt) = (integ_points(quad.ind_pts_integ(nelem) + pt, 1) - xp(nelem,1))*invh;

          if (order_ == 1) continue;
          fbasis(3, pt) = (integ_points(quad.ind_pts_integ(nelem) + pt, 0) - xp(nelem,0))*invh*(integ_points(quad.ind_pts_integ(nelem) + pt, 1) - xp(nelem,1))*invh; //xy
          fbasis(4, pt) = (integ_points(quad.ind_pts_integ(nelem) + pt, 0) - xp(nelem,0))*invh*(integ_points(quad.ind_pts_integ(nelem) + pt, 0) - xp(nelem,0))*invh; //x2
          fbasis(5, pt) = (integ_points(quad.ind_pts_integ(nelem) + pt, 1) - xp(nelem,1))*invh*(integ_points(quad.ind_pts_integ(nelem) + pt, 1) - xp(nelem,1))*invh; //y2
          if (order_ == 2) continue;
          throw;
        }
    }
  else if (Objet_U::dimension == 3)
    throw; //TODO
  else
    Process::exit();

  if (is_orthonormalized_)
    orthonormalize(nelem, quad.nb_pts_integ(nelem), fbasis);
}

/**
 * @brief Applies the pre-computed Gram-Schmidt transition matrix to a raw basis evaluation.
 *
 * @details Iterates over basis functions in reverse order (exploiting the upper-triangular
 * structure of the transition matrix) and replaces each raw value with the linear
 * combination given by the corresponding row of transition_matrix_. Handles both 2D
 * arrays (scalar basis: fbasis(i, k)) and 3D arrays (gradient basis: fbasis(i, k, d)).
 *
 * @param nelem         Element index, used to locate the block in transition_matrix_.
 * @param nb_pts_integ  Number of quadrature points to transform.
 * @param fbasis        Basis (or gradient) array to transform in-place.
 */
void BasisFunction::orthonormalize(const int& nelem, const int& nb_pts_integ, DoubleTab& fbasis) const
{
  int current_indice = indices_glob_elem_(nelem);

  double multvect;

  for (int i=nb_bfunc_-1; i>=0; i--) // reverse to take advantage of the triangular matrix
    {
      for (int k = 0; k < nb_pts_integ; k++)
        {
          if (fbasis.nb_dim() == 2)
            {
              multvect = 0.;
              for (int j=0; j<i+1; j++)
                multvect += transition_matrix_(current_indice+i,current_indice+j)*fbasis(j,k);
              fbasis(i,k) = multvect;
            }
          else if (fbasis.nb_dim() == 3) // grad cases
            {
              for (int l = 0; l < Objet_U::dimension ; l++)
                {
                  multvect = 0.;
                  for (int j=0; j<i+1; j++)
                    multvect += transition_matrix_(current_indice+i,current_indice+j)*fbasis(j,k,l);
                  fbasis(i,k,l) = multvect;
                }
            }
        }
    }
}

/**
 * @brief Evaluates all basis functions at a set of arbitrary coordinates.
 *
 * @details Same polynomial evaluation as the quadrature-based overload but uses
 * a caller-provided coordinate array instead of the quadrature point table.
 * Useful for post-processing or point-wise evaluations.
 *
 * @param coords  Input coordinates array of shape (nb_points, dimension).
 * @param nelem   Element index (used for barycenter and mesh size).
 * @param fbasis  Output array of shape (nb_bfunc, nb_points), filled in-place.
 *
 * @note Only 2D and orders 0-2 are implemented. Throws for order > 2 or 3D.
 */
void BasisFunction::eval_bfunc(const DoubleTab& coords, const int& nelem, DoubleTab& fbasis) const
{

  assert(fbasis.dimension(0) == nb_bfunc_);

  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());
  const DoubleTab& xp = domaine.xp(); // barycentre elem

  const Quadrature_base& quad = domaine.get_quadrature();
  int nb_points = quad.nb_pts_integ(nelem);

  if (Objet_U::dimension == 2)
    {
      double invh = 1./sqrt(domaine.carre_pas_maille(nelem));

      for (int pt = 0; pt < nb_points; pt++)
        {
          fbasis(0, pt) = 1;
          if (order_ == 0) continue;

          fbasis(1, pt) = (coords(pt, 0) - xp(nelem,0))*invh;
          fbasis(2, pt) = (coords(pt, 1) - xp(nelem,1))*invh;

          if (order_ == 1) continue;

          fbasis(3, pt) = (coords(pt, 0) - xp(nelem,0))*invh*(coords(pt, 1) - xp(nelem,1))*invh;
          fbasis(4, pt) = (coords(pt, 0) - xp(nelem,0))*invh*(coords(pt, 0) - xp(nelem,0))*invh;
          fbasis(5, pt) = (coords(pt, 1) - xp(nelem,1))*invh*(coords(pt, 1) - xp(nelem,1))*invh;

          if (order_ == 2) continue;
          throw;
        }
    }
  else if (Objet_U::dimension == 3)
    throw; //TODO
  else
    Process::exit();

  if (is_orthonormalized_)
    orthonormalize(nelem, nb_points, fbasis);
}

/* @brief Evaluation of the divergence of the basis functions (need to indicate which scalar component we are using) on integration points for elements
 *
 * @param quad Quadrature used to evaluate the basis functions
 * @param nelem Index of the element
 * @param dim Scalar component for which the divergence is computed
 * @param div_fbasis Output tab containing the divergence of the basis functions in dimension dim at integration points
 */
void BasisFunction::eval_div_bfunc(const Quadrature_base& quad, const int& nelem, DoubleTab& div_fbasis) const
{
  int nb_pts_integ_max = div_fbasis.dimension(2);
  DoubleTab grad_fbase_elem(nb_bfunc_, nb_pts_integ_max, Objet_U::dimension);
  eval_grad_bfunc(quad, nelem, grad_fbase_elem);
  for(int dim=0; dim<Objet_U::dimension; dim++)
    for (int i=0; i<nb_bfunc_; i++)
      for (int j=0; j<quad.nb_pts_integ(nelem); j++)
        div_fbasis(dim,i,j) = grad_fbase_elem(i,j,dim);

}

/**
 * @brief Evaluates the gradients of all basis functions at the element quadrature points.
 *
 * @details Fills grad_fbasis(i, k, d) with the d-th component of grad(phi_i) at the
 * k-th quadrature point of element nelem. The constant basis function (index 0) has a
 * zero gradient (left implicitly as zero by the caller's initialization). If the basis
 * is orthonormalized, the transition matrix is applied via orthonormalize().
 *
 * @param quad        Quadrature rule providing integration point coordinates.
 * @param nelem       Element index.
 * @param grad_fbasis Output array of shape (nb_bfunc, nb_pts_integ_max, dimension), filled in-place.
 *
 * @note Only 2D and orders 1-2 are implemented. Throws for order > 2 or 3D.
 */
void BasisFunction::eval_grad_bfunc(const Quadrature_base& quad, const int& nelem, DoubleTab& grad_fbasis) const
{
//  const DoubleTab& integ_points_on_facets = quad.get_integ_points_facets();
  assert(grad_fbasis.dimension(0) == nb_bfunc_ && grad_fbasis.dimension(1) == quad.nb_pts_integ_max() && grad_fbasis.dimension(2) == Objet_U::dimension);

  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());

  int ind_elem=quad.ind_pts_integ(nelem);
  const DoubleTab& integ_points = quad.get_integ_points();
  const DoubleTab& xp = domaine.xp(); // barycentre elem

  if (Objet_U::dimension == 2)
    {
      double invh = 1./sqrt(domaine.carre_pas_maille(nelem));

      for (int pt = 0; pt < quad.nb_pts_integ(nelem); pt++)
        {

          grad_fbasis(1, pt, 0) = invh;
          grad_fbasis(1, pt, 1) = 0.;
          grad_fbasis(2, pt, 0) = 0.;
          grad_fbasis(2, pt, 1) = invh;

          if (order_ == 1) continue;
          grad_fbasis(3, pt, 0) = invh*(integ_points(ind_elem+pt, 1) - xp(nelem,1))*invh;
          grad_fbasis(3, pt, 1) = invh*(integ_points(ind_elem+pt, 0) - xp(nelem,0))*invh;
          grad_fbasis(4, pt, 0) = 2*invh*(integ_points(ind_elem+pt, 0) - xp(nelem,0))*invh;
          grad_fbasis(4, pt, 1) = 0.;
          grad_fbasis(5, pt, 0) = 0.;
          grad_fbasis(5, pt, 1) = 2*invh*(integ_points(ind_elem+pt, 1) - xp(nelem,1))*invh;

          if (order_ == 2) continue;

          throw;
        }
    }
  else if (Objet_U::dimension == 3)
    throw; //TODO
  else
    Process::exit();

  if (is_orthonormalized_)
    orthonormalize(nelem, quad.nb_pts_integ(nelem), grad_fbasis);
}

/**
 * @brief Evaluates all basis functions of element nelem at the quadrature points of face num_face.
 *
 * @details Same polynomial as eval_bfunc() but uses the face quadrature point coordinates
 * instead of the element interior points. The basis is still centered at the element
 * barycenter, so the evaluation is consistent with the interior values across the face.
 * If the basis is orthonormalized, the transition matrix is applied via orthonormalize().
 *
 * @param quad      Quadrature rule providing face integration point coordinates.
 * @param nelem     Element index (provides barycenter and mesh size).
 * @param num_face  Face index (selects the row in the face quadrature point table).
 * @param fbasis    Output array of shape (nb_bfunc, nb_pts_integ_facets), filled in-place.
 *
 * @note Only 2D and orders 0-2 are implemented. Throws for order > 2 or 3D.
 */
void BasisFunction::eval_bfunc_on_facets(const Quadrature_base& quad, const int& nelem, const int& num_face, DoubleTab& fbasis) const
{
  const DoubleTab& integ_points_on_facets = quad.get_integ_points_facets();
  int nb_pts_integ_on_facets = quad.nb_pts_integ_facets();

  assert(fbasis.dimension(0) == nb_bfunc_ && fbasis.dimension(1) == nb_pts_integ_on_facets);

  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());
  const DoubleTab& xp = domaine.xp(); // barycentre elem

  if (Objet_U::dimension == 2)
    {
      double invh = 1./sqrt(domaine.carre_pas_maille(nelem));

      for (int pt = 0; pt < nb_pts_integ_on_facets; pt++)
        {
          fbasis(0, pt) = 1;
          if (order_ == 0) continue;

          fbasis(1, pt) = (integ_points_on_facets(num_face, pt, 0) - xp(nelem,0))*invh;
          fbasis(2, pt) = (integ_points_on_facets(num_face, pt, 1) - xp(nelem,1))*invh;

          if (order_ == 1) continue;
          fbasis(3, pt) = (integ_points_on_facets(num_face, pt, 0) - xp(nelem,0))*invh*(integ_points_on_facets(num_face, pt, 1) - xp(nelem,1))*invh;
          fbasis(4, pt) = (integ_points_on_facets(num_face, pt, 0) - xp(nelem,0))*invh*(integ_points_on_facets(num_face, pt, 0) - xp(nelem,0))*invh;
          fbasis(5, pt) = (integ_points_on_facets(num_face, pt, 1) - xp(nelem,1))*invh*(integ_points_on_facets(num_face, pt, 1) - xp(nelem,1))*invh;
          if (order_ == 2) continue;
          throw;
        }
    }
  else if (Objet_U::dimension == 3)
    throw; //TODO
  else
    Process::exit();

  if (is_orthonormalized_)
    orthonormalize(nelem, nb_pts_integ_on_facets, fbasis);
}

/**
 * @brief Evaluates the gradients of all basis functions of element nelem at the quadrature points of face num_face.
 *
 * @details Same gradient formulas as eval_grad_bfunc() but evaluated at face quadrature
 * points. Used in the SIP consistency and symmetry terms where the normal flux
 * { nu * grad(phi_i) } . n must be integrated over a face. If the basis is
 * orthonormalized, the transition matrix is applied via orthonormalize().
 *
 * @param quad        Quadrature rule providing face integration point coordinates.
 * @param nelem       Element index (provides barycenter and mesh size).
 * @param num_face    Face index (selects the row in the face quadrature point table).
 * @param grad_fbasis Output array of shape (nb_bfunc, nb_pts_integ_facets, dimension), filled in-place.
 *
 * @note Only 2D and orders 1-2 are implemented. Throws for order > 2 or 3D.
 */
void BasisFunction::eval_grad_bfunc_on_facets(const Quadrature_base& quad, const int& nelem, const int& num_face, DoubleTab& grad_fbasis) const
{

  const DoubleTab& integ_points_on_facets = quad.get_integ_points_facets();
  int nb_pts_integ_on_facets = quad.nb_pts_integ_facets();
  assert(grad_fbasis.dimension(0) == nb_bfunc_ && grad_fbasis.dimension(1) == nb_pts_integ_on_facets && grad_fbasis.dimension(2) == Objet_U::dimension );

  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());
  const DoubleTab& xp = domaine.xp(); // barycentre elem


  if (Objet_U::dimension == 2)
    {
      double invh = 1./sqrt(domaine.carre_pas_maille(nelem));

      for (int pt = 0; pt < nb_pts_integ_on_facets; pt++)
        {

          grad_fbasis(1, pt, 0) = invh;
          grad_fbasis(1, pt, 1) = 0.;
          grad_fbasis(2, pt, 0) = 0.;
          grad_fbasis(2, pt, 1) = invh;

          if (order_ == 1) continue;

          grad_fbasis(3, pt, 0) = invh*(integ_points_on_facets(num_face, pt, 1) - xp(nelem,1))*invh;  //xy
          grad_fbasis(3, pt, 1) = invh*(integ_points_on_facets(num_face, pt, 0) - xp(nelem,0))*invh;  //xy
          grad_fbasis(4, pt, 0) = 2*invh*(integ_points_on_facets(num_face, pt, 0) - xp(nelem,0))*invh;  //xx
          grad_fbasis(4, pt, 1) = 0.;  //xy
          grad_fbasis(5, pt, 0) = 0.;  //xy
          grad_fbasis(5, pt, 1) = 2*invh*(integ_points_on_facets(num_face, pt, 1) - xp(nelem,1))*invh;  //xy
          if (order_ == 2) continue;
          throw;
        }
    }
  else if (Objet_U::dimension == 3)
    throw; //TODO
  else
    Process::exit();

  if (is_orthonormalized_)
    orthonormalize(nelem, nb_pts_integ_on_facets, grad_fbasis);
}

/**
 * @brief Evaluates the divergence of the basis functions at the quadrature points of face num_face.
 *
 * @details Computes the gradient via eval_grad_bfunc_on_facets() and then rearranges the
 * result so that div_fbasis(d, i, k) = d(phi_i)/dx_d at the k-th face quadrature point.
 * This permuted layout is used when assembling divergence-based operators.
 *
 * @param quad       Quadrature rule providing face integration point coordinates.
 * @param nelem      Element index.
 * @param num_face   Face index.
 * @param div_fbasis Output array of shape (nb_bfunc, nb_pts_integ_facets), filled in-place.
 */
void BasisFunction::eval_div_bfunc_on_facets(const Quadrature_base& quad, const int& nelem, const int& num_face,  DoubleTab& div_fbasis) const
{
  assert(nb_bfunc_ == div_fbasis.dimension(0));
  int nb_pts_integ_max = div_fbasis.dimension(1);
  DoubleTab grad_fbase_elem(nb_bfunc_, nb_pts_integ_max, Objet_U::dimension);
  eval_grad_bfunc_on_facets(quad, nelem, num_face, grad_fbase_elem);
  for(int dim=0; dim<Objet_U::dimension; dim++)
    for (int i=0; i<nb_bfunc_; i++)
      for (int j=0; j<quad.nb_pts_integ(nelem); j++)
        div_fbasis(dim,i,j) = grad_fbase_elem(i,j,dim);
}

/**
 * @brief Computes the inverse of the local L2 mass matrix for element nelem.
 *
 * @details Two strategies are used depending on the polynomial order:
 *  - **Order 1 (2D)**: The 3x3 mass matrix has an analytic block structure. The (0,0)
 *    entry is 1/volume, and the lower-right 2x2 block (linear DOFs) is inverted
 *    analytically using the determinant of the moment integrals sum_x2, sum_y2, sum_xy.
 *  - **Order 2 (2D)**: The full 6x6 mass matrix is assembled by build_local_mass_matrix()
 *    and then numerically inverted via Matrice_Dense::inverse().
 *
 * @param quad   Quadrature rule used to compute the mass matrix entries.
 * @param nelem  Element index.
 * @return A nb_bfunc x nb_bfunc dense matrix containing M^{-1}.
 *
 * @note 3D is not yet implemented and throws at runtime.
 */
const Matrice_Dense BasisFunction::eval_invMassMatrix(const Quadrature_base& quad, const int& nelem) const
{

  // TODo DG adapt to high order quadrature ? possible ? or inverse mass matrix with specific quadrature ?

  const DoubleTab& weights = quad.get_weights();
  int nb_pts_integ_max = quad.nb_pts_integ_max();

  Matrice_Dense matrice(nb_bfunc_, nb_bfunc_);

  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());
  const DoubleVect& ve = domaine.volumes();

  DoubleTab fbase(nb_bfunc_, nb_pts_integ_max);

  eval_bfunc(quad, nelem, fbase);

  if (Objet_U::dimension == 2 && order_==1) // TODO: When general order, make it local to the cell
    {
      double invV = 1./ve(nelem);
      matrice(0, 0) = invV;

      matrice(0, 1) = 0.;
      matrice(0, 2) = 0.;

      matrice(1, 0) = 0.;
      matrice(2, 0) = 0.;

      double sumx2 = 0.;
      double sumy2 = 0.;
      double sumxy = 0.;

      for (int f = 0; f < quad.nb_pts_integ(nelem); f++)
        {
          sumx2 += weights(quad.ind_pts_integ(nelem)+f)*ve(nelem)*fbase(1,f)*fbase(1,f);
          sumxy += weights(quad.ind_pts_integ(nelem)+f)*ve(nelem)*fbase(1,f)*fbase(2,f);;
          sumy2 += weights(quad.ind_pts_integ(nelem)+f)*ve(nelem)*fbase(2,f)*fbase(2,f);
        }

      double inv_det = 1./(sumy2*sumx2 - sumxy*sumxy);

      matrice(1, 1) = sumy2*inv_det;
      matrice(2, 2) = sumx2*inv_det;

      matrice(1, 2) = -sumxy*inv_det;
      matrice(2, 1) = -sumxy*inv_det;
    }
  else if (Objet_U::dimension == 2 && order_==2)
    {
      matrice=build_local_mass_matrix(quad,  nelem); // Creation of the local mass matrix.
      matrice.inverse(); // Inversion of the local mass matrix.
    }
  else if (Objet_U::dimension == 3)
    throw; //TODO
  else
    Process::exit();

  return matrice;
}

/**
 * @brief Computes the SIP penalty stabilization parameters eta_elem and eta_facet.
 *
 * @details Element parameters eta_elem(e) depend on the polynomial order and the
 * element geometry:
 *  - **Triangles**: eta = (6/pi) * sigma^2 * (p+1)*(p+2), where sigma is the element
 *    shape parameter from Domaine_DG::get_sig(). This is a theoretically grounded
 *    lower bound for coercivity of the SIP bilinear form on triangles.
 *  - **Other polygons/polyhedra**: eta = 10*p^2, a conservative estimate pending a
 *    more geometry-aware formula.
 *
 * Face parameters eta_facet(f) are derived from the adjacent element values:
 *  - **Boundary faces**: eta_facet = eta_elem of the single adjacent element.
 *  - **Internal faces**: eta_facet = harmonic mean of the two adjacent element values,
 *    providing a balanced penalty that accounts for differing element sizes or orders
 *    on each side.
 */
void BasisFunction::compute_stab_param()
{
  const Domaine_DG& domaine = ref_cast(Domaine_DG,dom_.valeur());
  int nb_elem_tot = domaine.nb_elem_tot();
  eta_elem.resize(nb_elem_tot);
  eta_facet.resize(domaine.nb_faces());
  const IntTab& nfaces_elem = domaine.get_nfaces_elem();  // IntTab that indicate the number of facet of each elem

  //          Computation of the stabilisation parameters for triangle
  const DoubleTab& sig=domaine.get_sig();
  for (int e = 0; e < nb_elem_tot; e++)
    {
      double ordre = get_order();
      if(nfaces_elem(e)==3) // triangle
        {
          eta_elem(e) = 6. / M_PI * (sig(e)*sig(e))*(ordre + 1.)*(ordre + 2.);
        }
      else   // Polyhedra
        {
          eta_elem(e) = 10*ordre*ordre; //  TODO: calculate a more optimized pen. coeff.
        }
    }

  for (unsigned int f = 0; f < (unsigned int) domaine.premiere_face_int(); f++) // For the boundary
    {
      int elem0 = domaine.face_voisins(f, 0);
      eta_facet(f) = eta_elem(elem0); // EtaF=EtaT
    }

  for (unsigned int f = domaine.premiere_face_int(); f < (unsigned int) domaine.nb_faces(); f++)
    {
      unsigned int elem0 = domaine.face_voisins(f, 0);
      unsigned int elem1 = domaine.face_voisins(f, 1);
      double eta_t0 = eta_elem(elem0);
      double eta_t1 = eta_elem(elem1);
      eta_facet(f) = 2. * eta_t0 * eta_t1 / (eta_t0 + eta_t1); // Harmonique mean
    }

  /*for (int e = 0; e < domaine.nb_elem(); e++)
    {
      eta_elem(e) = 50.;     // Todo DG : replace for generic order
    }
  for (unsigned int f = 0; f < (unsigned int) domaine.nb_faces(); f++) // For the boundary
    {
      eta_facet(f) = 10.; // EtaF=EtaT
    }*/
}

