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

void Op_Grad_DG::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base&)
{
  le_dom_DG = ref_cast(Domaine_DG, domaine_dis);
  le_dcl_DG = ref_cast(Domaine_Cl_DG, domaine_Cl_dis);
}

void Op_Grad_DG::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const std::string nom_inco = equation().inconnue().le_nom().getString();
  if (semi_impl.count("pression"))
    return; // semi-implicite -> rien a dimensionner

  int n_ext = 1; //TODO DG what is op_ext in this case ?
  std::vector<Matrice_Morse *> mat(n_ext);
  for (int i = 0; i < n_ext; i++)
    {
      std::string nom_mat = i ? nom_inco + "/" + (this)->equation().probleme().le_nom().getString() : nom_inco; //TODO DG is that correspond ?
      mat[i] = matrices.count(nom_mat) ? matrices.at(nom_mat) : nullptr;
      if (!mat[i])
        continue;
      Matrice_Morse mat2;
      if (i == 0)
        dimensionner(mat2);
      else
        throw; // TODO DG for dimensionner_terme_croises

      mat[i]->nb_colonnes() ? *mat[i] += mat2 : *mat[i] = mat2;
    }
}

void Op_Grad_DG::dimensionner(Matrice_Morse& matrice) const
{
  if (has_interface_blocs())
    {
      Operateur_base::dimensionner(matrice);
      return;
    }

  const Domaine_DG& domaine_DG = le_dom_DG.valeur();
  int nb_faces = domaine_DG.nb_faces();
  int nb_faces_tot = domaine_DG.nb_faces_tot();
  int nb_elem_tot = domaine_DG.nb_elem_tot();
  IntTab stencil(0, 2);

  const IntTab& face_voisins = domaine_DG.face_voisins();

  int nb_coef = 0;
  for (int face = 0; face < nb_faces; face++)
    {
      for (int dir = 0; dir < 2; dir++)
        {
          const int elem = face_voisins(face, dir);
          if (elem != -1)
            {
              stencil.resize(nb_coef + 1, 2);
              stencil(nb_coef, 0) = elem;
              stencil(nb_coef, 1) = face;
              nb_coef++;
            }
        }
    }
  tableau_trier_retirer_doublons(stencil);
  Matrix_tools::allocate_morse_matrix(nb_elem_tot, nb_faces_tot, stencil, matrice);
}

void Op_Grad_DG::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{

  const DoubleTab& inco_p = semi_impl.count("pression") ? semi_impl.at("pression") : equation().inconnue().valeurs(); // NB : is this working ?
  Matrice_Morse *mat = matrices.count("pression") ? matrices.at("pression") : nullptr; // pression for the stabilisation term if np==nv

  const Domaine_DG& domaine = le_dom_DG.valeur();
  const IntTab& face_voisins = domaine.face_voisins();

  int order_v = Option_DG::Get_order_for("vitesse");
  int order_p = Option_DG::Get_order_for("pression");

  const BasisFunction& bfunc_v = domaine.get_basisFunction(order_v);
  const int nb_bfunc_v = bfunc_v.nb_bfunc();

  const BasisFunction& bfunc_p = domaine.get_basisFunction(order_p);
  const int nb_bfunc_p = bfunc_p.nb_bfunc();

  const int quad_order = bfunc_p.get_default_quadrature_order(); //TODO DG should we choose the max or the p order quadrature ?
  const Quadrature_base& quad = domaine.get_quadrature(quad_order);  // Same quadrature for all champs
  int nb_pts_integ_max = quad.nb_pts_integ_max();
  double coeff, coeff00, coeff10, coeff01, coeff11;

  DoubleTab grad_fbase_elem(nb_bfunc_p, nb_pts_integ_max, Objet_U::dimension);
  DoubleTab f_base_v(nb_bfunc_v, nb_pts_integ_max);
  DoubleTab scalar_product_dim(nb_pts_integ_max);

  const IntTab& indices_glob_elem_v = bfunc_v.indices_glob_elem();
  const IntTab& indices_glob_elem_p = bfunc_p.indices_glob_elem();

  const int dim = Objet_U::dimension;

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
              scalar_product_dim = 0.;
              for (int k = 0; k < quad.nb_pts_integ(elem); k++)
                scalar_product_dim(k) += grad_fbase_elem(pressure_index, k, d) * f_base_v(velocity_index, k);
              coeff = quad.compute_integral_on_elem(elem, scalar_product_dim);
              if (mat)
                (*mat)(ind_elem_v*dim + velocity_index + d * nb_bfunc_v, ind_elem_p + pressure_index) -= coeff;
              secmem(elem, velocity_index + d * nb_bfunc_v) += coeff * inco_p(elem, pressure_index);
            }
    }

  const DoubleTab& face_normales = domaine.face_normales();
  const DoubleVect& face_surfaces = domaine.face_surfaces();
  int nb_pts_int_fac = quad.nb_pts_integ_facets();

  DoubleTab eval_jump_on_facet00(nb_pts_int_fac);
  DoubleTab eval_jump_on_facet01(nb_pts_int_fac);
  DoubleTab eval_jump_on_facet10(nb_pts_int_fac);
  DoubleTab eval_jump_on_facet11(nb_pts_int_fac);
  //DoubleTab mean_v(nb_pts_int_fac);

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
      int ind_elem0_v = bfunc_v.indices_glob_elem(elem0);
      int ind_elem1_v = bfunc_v.indices_glob_elem(elem1);
      int ind_elem0_p = bfunc_p.indices_glob_elem(elem0);
      int ind_elem1_p = bfunc_p.indices_glob_elem(elem1);
      bfunc_v.eval_bfunc_on_facets(quad, elem0, face, f_base_v0);
      bfunc_v.eval_bfunc_on_facets(quad, elem1, face, f_base_v1);
      bfunc_p.eval_bfunc_on_facets(quad, elem0, face, f_base_p0);
      bfunc_p.eval_bfunc_on_facets(quad, elem1, face, f_base_p1);
      for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
        {
          //for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
          //  mean_v(k) = 0.5*(f_base_v0(velocity_index, k) + f_base_v1(velocity_index, k));

          for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
            {
              for (int d = 0; d < Objet_U::dimension; d++)
                {
                  eval_jump_on_facet00 = 0.;
                  eval_jump_on_facet01 = 0.;
                  eval_jump_on_facet10 = 0.;
                  eval_jump_on_facet11 = 0.;
                  for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                    {
                      eval_jump_on_facet00(k) -= f_base_p0(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v0(velocity_index, k) / sur_f;
                      eval_jump_on_facet01(k) += f_base_p1(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v0(velocity_index, k) / sur_f;
                      eval_jump_on_facet10(k) -= f_base_p0(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v1(velocity_index, k) / sur_f;
                      eval_jump_on_facet11(k) += f_base_p1(pressure_index, k) * face_normales(face, d)  * 0.5 * f_base_v1(velocity_index, k) / sur_f;
                    }
                  coeff00 = quad.compute_integral_on_facet(face, eval_jump_on_facet00);
                  coeff01 = quad.compute_integral_on_facet(face, eval_jump_on_facet01);
                  coeff10 = quad.compute_integral_on_facet(face, eval_jump_on_facet10);
                  coeff11 = quad.compute_integral_on_facet(face, eval_jump_on_facet11);

                  if (mat)
                    {
                      (*mat)(ind_elem0_v*dim + velocity_index + d * nb_bfunc_v, ind_elem0_p + pressure_index) += coeff00;
                      (*mat)(ind_elem0_v*dim + velocity_index + d * nb_bfunc_v, ind_elem1_p + pressure_index) += coeff01;
                      (*mat)(ind_elem1_v*dim + velocity_index + d * nb_bfunc_v, ind_elem0_p + pressure_index) += coeff10;
                      (*mat)(ind_elem1_v*dim + velocity_index + d * nb_bfunc_v, ind_elem1_p + pressure_index) += coeff11;
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
