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

#ifndef BasisFunction_included
#define BasisFunction_included

#include <Quadrature_base.h>
#include <Option_DG.h>
#include <Domaine_DG.h>
#include <Matrice_Dense.h>
#include <Matrice_Morse.h>

class BasisFunction
{
public:
  BasisFunction() = default;
  virtual ~BasisFunction() {}

  void initialize(const Domaine_DG& dom, const int& order, const bool& gram_schmidt);

  inline const int& get_order() const { return order_; }
  inline const int& get_default_quadrature_order() const { return default_quad_order_; }
  inline const IntTab& indices_glob_elem(const int dim = 1) const
  {
    switch (dim)
      {
      case 1:
        return indices_glob_elem_;
      case 2:
        if (indices_glob_elem_2D_.size() == 0)
          {
            indices_glob_elem_2D_ = indices_glob_elem_;
            indices_glob_elem_2D_ *= 2;
          }
        return indices_glob_elem_2D_;
      case 3:
        if (indices_glob_elem_3D_.size() == 0)
          {
            indices_glob_elem_3D_ = indices_glob_elem_;
            indices_glob_elem_3D_ *= 2;
          }
        return indices_glob_elem_3D_;
      default:
        Cerr << "bad dimension indices_glob_elem" << finl;
        throw;
      }
    return indices_glob_elem_;
  }
//  inline const int& indices_glob_elem(int elem) const { return indices_glob_elem_(elem); }
  inline const int& nb_bfunc() const { return nb_bfunc_; }

  //Evaluation of the basis functions on integration points for elements and facets
  void eval_bfunc(const Quadrature_base& quad, const int& nelem, DoubleTab& fbasis) const;
  void eval_bfunc_on_facets(const Quadrature_base& quad, const int& nelem, const int& num_face, DoubleTab& grad_fbasis) const;
  void eval_bfunc(const DoubleTab& coord, const int& nelem, DoubleTab& fbasis) const;

  //Evaluation of the gradient of the basis functions on integration points for elements and facets
  void eval_grad_bfunc(const Quadrature_base& quad, const int& nelem, DoubleTab& fbasis) const;
  void eval_grad_bfunc_on_facets(const Quadrature_base& quad, const int& nelem, const int& num_face, DoubleTab& grad_fbasis) const;

  //Evaluation of the divergence of the basis functions (need to indicate which scalar component we are using) on integration points for elements and facets
  void eval_div_bfunc_on_facets(const Quadrature_base& quad, const int& nelem, const int& num_face, DoubleTab& div_fbasis) const;
  void eval_div_bfunc(const Quadrature_base& quad, const int& nelem, DoubleTab& div_fbasis) const;

  const Matrice_Dense eval_invMassMatrix(const Quadrature_base& quad, const int& nelem) const;

  inline const DoubleTab& get_eta_elem() const { return eta_elem; }
  inline const DoubleTab& get_eta_facet() const { return eta_facet; }

protected:
  /*! Compute the mass matrix
   */
  void allocate_mass_matrix();
  void allocate_transition_matrix();
  void compute_stab_param();

  const Matrice_Dense build_local_mass_matrix(const Quadrature_base& quad, const int nelem) const;

  void build_mass_matrix();
  void build_transition_matrix();
  void orthonormalize(const Quadrature_base& quad, const int& nelem, const int& nb_pts_integ, DoubleTab& fbasis) const;

  void gramSchmidt(DoubleTab& fbase, const Quadrature_base& quad, const int& num_elem, const int& current_indice, const int& nb_pts_integ, const double& volume, int index);

  OBS_PTR(Domaine_DG) dom_;
  int order_;
  bool is_orthonormalized_;
  bool is_diagonal_;
  int nb_bfunc_;
  int default_quad_order_;

  IntTab indices_glob_elem_;
  mutable IntTab indices_glob_elem_2D_;
  mutable IntTab indices_glob_elem_3D_;


  Matrice_Morse transition_matrix_;

  DoubleTab eta_elem;
  DoubleTab eta_facet;

};


#endif /* BasisFunction_included */
