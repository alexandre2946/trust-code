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

// Recursive Gram-Schmidt orthogonalization function with transition matrix
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

/*! @brief Compute the mass matrix of cell nelem
 *
 * @param quad quadature used to compute mass matrix
 * @param nelem index of the cell
 *
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
    orthonormalize(nelem, fbasis);
}

void BasisFunction::orthonormalize(const int& nelem, DoubleTab& fbasis) const
{
  int current_indice = indices_glob_elem_(nelem);
  const int nb_pts_integ= fbasis.dimension(1);

  double multvect;

  for (int i=nb_bfunc_-1; i>=0; i--) // reverse to take advantage of the triangular matrix
    {
      for (int k = 0; k < nb_pts_integ ; k++)
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
    orthonormalize(nelem, fbasis);
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
  int nb_pts_integ_max = div_fbasis.dimension(1);
  DoubleTab grad_fbase_elem(nb_bfunc_, nb_pts_integ_max, Objet_U::dimension);
  eval_grad_bfunc(quad, nelem, grad_fbase_elem);
  for(int dim=0; dim<Objet_U::dimension; dim++)
    for (int i=0; i<nb_bfunc_; i++)
      for (int j=0; j<quad.nb_pts_integ(nelem); j++)
        div_fbasis(dim,i,j) = grad_fbase_elem(i,j,dim);

}

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
    orthonormalize(nelem, grad_fbasis);
}


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
    orthonormalize(nelem, fbasis);
}

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
    orthonormalize(nelem, grad_fbasis);
}

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

