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

#include <Modele_turbulence_scal_base.h>
#include <Echange_externe_impose.h>
#include <Op_Diff_DG_Elem.h>
#include <Dirichlet_homogene.h>
#include <Domaine_Cl_DG.h>
#include <Champ_Elem_DG.h>
#include <Schema_Temps_base.h>
#include <Champ_front_calc.h>
#include <Domaine_DG.h>
#include <communications.h>
#include <Synonyme_info.h>
#include <Probleme_base.h>
#include <Neumann_paroi.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTLists.h>
#include <Dirichlet.h>
#include <cmath>
#include <Quadrature_base.h>
#include <Champ_front_txyz.h>
#include <Champ_front_softanalytique.h>
#include <BasisFunction.h>

Implemente_instanciable(Op_Diff_DG_Elem, "Op_Diff_DG_Elem", Op_Diff_DG_base);

Sortie& Op_Diff_DG_Elem::printOn(Sortie& os) const { return Op_Diff_DG_base::printOn(os); }

Entree& Op_Diff_DG_Elem::readOn(Entree& is) { return Op_Diff_DG_base::readOn(is); }


/**
 * @brief Finalizes operator setup after all associations have been made.
 *
 * @details Calls the parent completer(), then:
 *  - Casts the unknown field to Champ_Elem_DG and checks that the ghost cell layer
 *    is at least 1 element thick (required for face-neighbour communication in DG).
 *  - Initializes the face boundary flux array flux_bords_ with the correct number
 *    of components (2 for Transport_K_Epsilon, otherwise the line size of the unknown).
 *  - If the operator is a turbulent diffusion operator (name starts with "Op_Dift"),
 *    retrieves the turbulent conductivity from the turbulence model and registers it
 *    as the turbulent diffusivity via associer_diffusivite_turbulente().
 */
void Op_Diff_DG_Elem::completer()
{
  Op_Diff_DG_base::completer();
  const Champ_Elem_DG& ch = ref_cast(Champ_Elem_DG, equation().inconnue());
  const Domaine_DG& domaine = le_dom_dg_.valeur();
  if (domaine.domaine().nb_joints() && domaine.domaine().joint(0).epaisseur() < 1)
    Cerr << "Op_Diff_DG_Elem :  ghost cell layer width too small, minimum 1" << finl, Process::exit();
  ch.fcl();
  int nb_comp = (equation().que_suis_je() == "Transport_K_Epsilon") ? 2 : ch.valeurs().line_size();
  flux_bords_.resize(domaine.premiere_face_int(), nb_comp);

  if (!que_suis_je().debute_par("Op_Dift"))
    return;

  const RefObjU& modele_turbulence = equation().get_modele(TURBULENCE);
  const Modele_turbulence_scal_base& mod_turb = ref_cast(Modele_turbulence_scal_base, modele_turbulence.valeur());
  const Champ_Fonc_base& lambda_t = mod_turb.conductivite_turbulente();
  associer_diffusivite_turbulente(lambda_t);
}

/**
 * @brief Builds the sparsity pattern of the DG diffusion matrix in a Matrice_Morse.
 *
 * @details The matrix couples each element to all its face-neighbours as given by the
 * pre-computed sorted stencil (Domaine_DG::get_stencil_sorted()). The global index
 * space is built from the BasisFunction index map indices_glob_elem, so that each
 * element contributes a block of nb_basis_func * dim rows and the column range of a
 * row spans nb_basis_func columns per neighbour element (including itself).
 *
 * The method fills tab1 (row pointers) in a first pass, then tab2 (column indices)
 * in a second pass, and finally marks the stencil as sorted.
 *
 * @param la_matrice The Matrice_Morse whose sparsity pattern is to be set.
 */
void Op_Diff_DG_Elem::dimensionner(Matrice_Morse& la_matrice) const // TODO a remonter dans Op_DG_Elem
{

  const Nom& nom_inco = equation().inconnue().le_nom();
  int nordre = Option_DG::Get_order_for(nom_inco);
  int dim = equation().inconnue().is_vectorial() ? Objet_U::dimension : 1;

  const Domaine_DG& domaine = le_dom_dg_.valeur();

  const BasisFunction& bfunc = domaine.get_basisFunction(nordre);
  const int nb_basis_func = bfunc.nb_bfunc();

  const IntTab& indices_glob_elem =bfunc.indices_glob_elem(dim);

  int nb_elem_tot = le_dom_dg_->nb_elem_tot();
  int size_inc = indices_glob_elem(nb_elem_tot);

  const Stencil& stencil_sorted = domaine.get_stencil_sorted();
  const int nb_stencil_max = stencil_sorted.dimension(1);

  la_matrice.dimensionner(size_inc, size_inc, 0);

  auto& tab1 = la_matrice.get_set_tab1();
  auto& tab2 = la_matrice.get_set_tab2();
  auto& coeff = la_matrice.get_set_coeff();
  coeff = 0;

  int col, indice;

  tab1(0) = 1;
  for (int nelem = 0; nelem < nb_elem_tot; nelem++)
    {
      int nb_indices_line = 0;
      for (int k = 0 ; k < nb_stencil_max; k++)
        {
          if (stencil_sorted(nelem, k) < 0)
            break;
          nb_indices_line += nb_basis_func;
        }
      for (int k = 0; k < nb_basis_func * dim; k++)
        tab1(indices_glob_elem(nelem) + k + 1) = nb_indices_line + tab1(indices_glob_elem(nelem) + k);
    }

  la_matrice.dimensionner(size_inc, tab1(size_inc) - 1);

  for (int nelem = 0; nelem < nb_elem_tot; nelem++)
    {
      auto row = tab1[indices_glob_elem(nelem)]-1 ;
      auto nb_indices_line = tab1[indices_glob_elem(nelem)+1] - tab1[indices_glob_elem(nelem)];
      indice = 0;
      for (int d = 0; d < dim; d++)
        {
          for (int i = 0; i < nb_basis_func; i++)
            {
              for (int k = 0; k < nb_stencil_max; k++)
                {
                  if (stencil_sorted(nelem, k) < 0)
                    break;
                  col = indices_glob_elem(stencil_sorted(nelem, k)) + 1;
                  for (int j = 0 ;  j < nb_basis_func; j++)
                    tab2[row + indice + j + k*nb_basis_func] = col + j + d*nb_basis_func;
                }
              indice += nb_indices_line;
            }
        }
    }
  la_matrice.is_sorted_stencil();
  assert(la_matrice.is_sorted_stencil());
}

/**
 * @brief Sizes the block matrices used by the interface_blocs assembly mechanism.
 *
 * @details If the unknown is treated semi-implicitly (its name appears in semi_impl),
 * no matrix needs to be dimensioned and the method returns immediately.
 * Otherwise, for each external operator registered in op_ext (used for monolithic
 * thermal coupling), the corresponding sub-matrix is sized by calling dimensionner().
 * Cross-problem coupling (i > 0) is not yet implemented and throws at runtime.
 *
 * @param matrices Map of matrix name → Matrice_Morse pointer to be sized.
 * @param semi_impl Map of semi-implicit field names to their current values.
 */
void Op_Diff_DG_Elem::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const std::string nom_inco = equation().inconnue().le_nom().getString();
  if (semi_impl.count(nom_inco))
    return; // semi-implicite -> rien a dimensionner

  init_op_ext();                  // TODO DG a completer
  int n_ext = (int)op_ext.size(); // pour la thermique monolithique

  std::vector<Matrice_Morse *> mat(n_ext);
  for (int i = 0; i < n_ext; i++)
    {
      std::string nom_mat = i ? nom_inco + "/" + op_ext[i]->equation().probleme().le_nom().getString() : nom_inco;
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

/**
 * @brief Assembles the SIP diffusion operator into the matrix and right-hand side.
 *
 * @details This is the core assembly routine. It proceeds in three stages:
 *
 * **1. Volume integrals (stiffness term)**
 * For each element, the term  integral of nu * grad(phi_i) . grad(phi_j)  is integrated
 * using the element quadrature rule and accumulated into the diagonal block of the
 * matrix and into secmem.
 *
 * **2. Internal face integrals (SIP terms)**
 * For each internal face shared by elem0 and elem1:
 *  - *Penalty term*: gamma * (eta_F/h_T) * integral of phi_i * phi_j, added to the
 *    diagonal blocks (same-element test and trial) and subtracted from the off-diagonal
 *    blocks (cross-element), enforcing the jump penalization symmetrically.
 *  - *Consistency + symmetry terms*: 0.5 * integral of { nu * grad(phi_i) } . n * phi_j,
 *    assembled into all four (elem0/elem1) x (elem0/elem1) block combinations with the
 *    appropriate sign to yield a symmetric bilinear form.
 *
 * **3. Boundary face integrals (Dirichlet enforcement)**
 * For boundary faces flagged as Dirichlet (fcl flag > 5), the same penalty and
 * consistency terms are applied to the single adjacent element, imposing the boundary
 * condition weakly in the SIP sense. The Dirichlet value contribution to secmem is
 * handled separately by contribuer_au_second_membre().
 *
 * Both isotropic and anisotropic diffusivities are supported: for anisotropic cases,
 * nu_F is computed as the normal projection of the diffusivity tensor onto the face.
 *
 * The DOF ordering within an element is: phi0.ex, phi1.ex, ..., phi0.ey, phi1.ey, ...
 * for vector fields (nb_bfunc DOFs per spatial direction).
 *
 * @param matrices  Map of matrix name → Matrice_Morse pointer to accumulate into.
 * @param secmem    Right-hand side array to accumulate into.
 * @param semi_impl Map of semi-implicit field values (matrix assembly is skipped if the
 *                  unknown is present here).
 */
void Op_Diff_DG_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{

  const Nom& nom_inco = equation().inconnue().le_nom();
  const std::string& nom_inco_str = equation().inconnue().le_nom().getString();
  const DoubleTab& inco = semi_impl.count(nom_inco_str) ? semi_impl.at(nom_inco_str) : equation().inconnue().valeurs();
  Matrice_Morse *mat = matrices.count(nom_inco_str) ? matrices.at(nom_inco_str) : nullptr;

  update_nu();

  const Domaine_DG& domaine = le_dom_dg_.valeur();
  const IntTab& face_voisins = domaine.face_voisins();

  const Champ_Elem_DG& ch = ref_cast(Champ_Elem_DG, equation().inconnue());
  const int dim = ch.is_vectorial() ? Objet_U::dimension : 1;

  int order = Option_DG::Get_order_for(nom_inco);

  const BasisFunction& bfunc = le_dom_dg_->get_basisFunction(order);
  const int nb_bfunc = bfunc.nb_bfunc();
  const IntTab& indices_glob_elem = bfunc.indices_glob_elem(dim);

  const DoubleTab& eta_F = bfunc.get_eta_facet(); // Compute the penalisation coefficient

  const int quad_order = bfunc.get_default_quadrature_order();
  const Quadrature_base& quad = domaine.get_quadrature(quad_order);
  int nb_pts_integ_max = quad.nb_pts_integ_max();
  double coeff;

  DoubleTab grad_fbase_elem(nb_bfunc, nb_pts_integ_max, Objet_U::dimension);
  DoubleTab diffusion(nb_pts_integ_max);

  for (int e = 0; e < le_dom_dg_->nb_elem(); e++)
    {
      bfunc.eval_grad_bfunc(quad, e, grad_fbase_elem);
      int ind_elem = indices_glob_elem(e);
      for (int i = 0; i < nb_bfunc; i++)
        for (int j = 0; j < nb_bfunc; j++)
          {
            diffusion = 0.;
            for (int k = 0; k < quad.nb_pts_integ(e); k++)
              for (int d = 0; d < Objet_U::dimension; d++)
                {
                  bool ori = is_aniso_ ? d : 0;
                  diffusion(k) += nu(e, ori) * grad_fbase_elem(i, k, d) * grad_fbase_elem(j, k, d);
                }

            coeff = quad.compute_integral_on_elem(e, diffusion);
            for (int d_base = 0; d_base < dim; d_base++)
              {
                if (mat)
                  (*mat)(ind_elem + i + d_base * nb_bfunc, ind_elem + j + d_base * nb_bfunc) += coeff;
                secmem(e, i + d_base * nb_bfunc) -= coeff * inco(e, j + d_base * nb_bfunc);
              }
          }
    }

  int nb_pts_int_fac = quad.nb_pts_integ_facets();
  int premiere_face_int = domaine.premiere_face_int();
  const DoubleTab& face_normales = domaine.face_normales();

  int elem0, elem1;

  DoubleTab product(nb_pts_int_fac);
  DoubleTab scalar_product(nb_pts_int_fac);

  DoubleTab fbase0(nb_bfunc, nb_pts_int_fac);
  DoubleTab fbase1(nb_bfunc, nb_pts_int_fac);

  DoubleTab grad_fbase0(nb_bfunc, nb_pts_int_fac, Objet_U::dimension);
  DoubleTab grad_fbase1(nb_bfunc, nb_pts_int_fac, Objet_U::dimension);

  for (int f = premiere_face_int; f < domaine.nb_faces(); f++)
    {

      elem0 = face_voisins(f, 0);
      elem1 = face_voisins(f, 1);

      int ind_elem0 = indices_glob_elem(elem0);
      int ind_elem1 = indices_glob_elem(elem1);

      double sur_f = domaine.face_surfaces(f);

      double h_T = sqrt(std::min(domaine.carre_pas_maille(elem0), domaine.carre_pas_maille(elem1))); // TODO possibilite de prendre moyenne harmonique (stabilite)
      double invh_T = 1. / h_T;

      double nu_0 = 0., nu_1 = 0.;
      if (is_aniso_)
        {
          for (int d = 0; d < Objet_U::dimension; d++)
            {
              nu_0 += nu(elem0, d) * face_normales(f, d) * face_normales(f, d);
              nu_1 += nu(elem1, d) * face_normales(f, d) * face_normales(f, d);
            }
          nu_0 /= sur_f * sur_f;
          nu_1 /= sur_f * sur_f;
        }
      else
        {
          nu_0 = nu(elem0, 0);
          nu_1 = nu(elem1, 0);
        }
      double gamma = 2 * nu_0 * nu_1 / (nu_0 + nu_1);

      //*****************//
      // penalizing term //
      //*****************//
      for (int i_elem = 0; i_elem < 2; i_elem++)
        {
          int elem = face_voisins(f, i_elem);
          int ind_elem = indices_glob_elem(elem);

          bfunc.eval_bfunc_on_facets(quad, elem, f, fbase0);

          for (int i = 0; i < nb_bfunc; i++)
            for (int j = 0; j < nb_bfunc; j++)
              {
                for (int k = 0; k < nb_pts_int_fac; k++)
                  product(k) = fbase0(i, k) * fbase0(j, k); // TODO DG kronecker ?

                coeff = gamma * eta_F(f) * invh_T * quad.compute_integral_on_facet(f, product);
                for (int d_base = 0; d_base < dim; d_base++)
                  {
                    if (mat)
                      (*mat)(ind_elem + i + d_base * nb_bfunc, ind_elem + j + d_base * nb_bfunc) += coeff;
                    secmem(elem, i + d_base * nb_bfunc) -= coeff * inco(elem, j + d_base * nb_bfunc);
                  }
              }
        }

      // crossed_term
      bfunc.eval_bfunc_on_facets(quad, elem0, f, fbase0);
      bfunc.eval_bfunc_on_facets(quad, elem1, f, fbase1);

      for (int i = 0; i < nb_bfunc; i++)
        for (int j = 0; j < nb_bfunc; j++)
          {
            for (int k = 0; k < nb_pts_int_fac; k++)
              product(k) = fbase0(i, k) * fbase1(j, k);

            double integral = quad.compute_integral_on_facet(f, product);
            coeff = gamma * eta_F(f) * invh_T * integral;
            for (int d_base = 0; d_base < dim; d_base++)
              {
                if (mat)
                  {
                    (*mat)(ind_elem0 + i + d_base * nb_bfunc, ind_elem1 + j + d_base * nb_bfunc) -= coeff;
                    (*mat)(ind_elem1 + j + d_base * nb_bfunc, ind_elem0 + i + d_base * nb_bfunc) -= coeff; // symmetry
                  }
                secmem(elem0, i + d_base * nb_bfunc) += coeff * inco(elem1, j + d_base * nb_bfunc);
                secmem(elem1, j + d_base * nb_bfunc) += coeff * inco(elem0, i + d_base * nb_bfunc); // symmetry
              }
          }

      //****************//
      // symmetric term //
      //****************//
      bfunc.eval_grad_bfunc_on_facets(quad, elem0, f, grad_fbase0);
      bfunc.eval_grad_bfunc_on_facets(quad, elem1, f, grad_fbase1);

      for (int i = 0; i < nb_bfunc; i++)
        {
          scalar_product = 0.;
          for (int k = 0; k < nb_pts_int_fac; k++)
            for (int d = 0; d < Objet_U::dimension; d++)
              {
                bool ori = is_aniso_ ? d : 0;
                scalar_product(k) += nu(elem0, ori) * face_normales(f, d) / sur_f * grad_fbase0(i, k, d);
              }

          for (int j = 0; j < nb_bfunc; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase0(j, k);
              double integral = quad.compute_integral_on_facet(f, product);
              for (int d_base = 0; d_base < dim; d_base++)
                {
                  if (mat)
                    {
                      (*mat)(ind_elem0 + i + d_base * nb_bfunc, ind_elem0 + j + d_base * nb_bfunc) -= 0.5 * integral;
                      (*mat)(ind_elem0 + j + d_base * nb_bfunc, ind_elem0 + i + d_base * nb_bfunc) -= 0.5 * integral; // symmetry
                    }
                  secmem(elem0, i + d_base * nb_bfunc) += 0.5 * integral * inco(elem0, j + d_base * nb_bfunc);
                  secmem(elem0, j + d_base * nb_bfunc) += 0.5 * integral * inco(elem0, i + d_base * nb_bfunc); // symmetry
                }
            }

          for (int j = 0; j < nb_bfunc; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase1(j, k);
              double integral = quad.compute_integral_on_facet(f, product);
              for (int d_base = 0; d_base < dim; d_base++)
                {
                  if (mat)
                    {
                      (*mat)(ind_elem0 + i + d_base * nb_bfunc, ind_elem1 + j + d_base * nb_bfunc) += 0.5 * integral;
                      (*mat)(ind_elem1 + j + d_base * nb_bfunc, ind_elem0 + i + d_base * nb_bfunc) += 0.5 * integral; // symmetry
                    }
                  secmem(elem0, i + d_base * nb_bfunc) -= 0.5 * integral * inco(elem1, j + d_base * nb_bfunc);
                  secmem(elem1, j + d_base * nb_bfunc) -= 0.5 * integral * inco(elem0, i + d_base * nb_bfunc); // symmetry
                }
            }

          scalar_product = 0.;
          for (int k = 0; k < nb_pts_int_fac; k++)
            for (int d = 0; d < Objet_U::dimension; d++)
              {
                bool ori = is_aniso_ ? d : 0;
                scalar_product(k) += nu(elem1, ori) * face_normales(f, d) / sur_f * grad_fbase1(i, k, d);
              }

          for (int j = 0; j < nb_bfunc; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase1(j, k);
              double integral = quad.compute_integral_on_facet(f, product);
              for (int d_base = 0; d_base < dim; d_base++)
                {
                  if (mat)
                    {
                      (*mat)(ind_elem1 + i + d_base * nb_bfunc, ind_elem1 + j + d_base * nb_bfunc) += 0.5 * integral;
                      (*mat)(ind_elem1 + j + d_base * nb_bfunc, ind_elem1 + i + d_base * nb_bfunc) += 0.5 * integral; // symmetry
                    }
                  secmem(elem1, i + d_base * nb_bfunc) -= 0.5 * integral * inco(elem1, j + d_base * nb_bfunc);
                  secmem(elem1, j + d_base * nb_bfunc) -= 0.5 * integral * inco(elem1, i + d_base * nb_bfunc); // symmetry
                }
            }

          for (int j = 0; j < nb_bfunc; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase0(j, k);
              double integral = quad.compute_integral_on_facet(f, product);
              for (int d_base = 0; d_base < dim; d_base++)
                {
                  if (mat)
                    {
                      (*mat)(ind_elem1 + i + d_base * nb_bfunc, ind_elem0 + j + d_base * nb_bfunc) -= 0.5 * integral;
                      (*mat)(ind_elem0 + j + d_base * nb_bfunc, ind_elem1 + i + d_base * nb_bfunc) -= 0.5 * integral; // symmetry
                    }
                  secmem(elem1, i + d_base * nb_bfunc) += 0.5 * integral * inco(elem0, j + d_base * nb_bfunc);
                  secmem(elem0, j + d_base * nb_bfunc) += 0.5 * integral * inco(elem1, i + d_base * nb_bfunc); // symmetry
                }
            }
        }
    }

  /* Treatment of the boundary conditions */

  for (int f = 0; f < premiere_face_int; f++) // For the boundary
    {

      if (ch.fcl()(f, 0) > 5)
        {
          int elem = face_voisins(f, 0); // The cell that have one facet on the boundary
          int ind_elem = indices_glob_elem(elem);

          bfunc.eval_bfunc_on_facets(quad, elem, f, fbase0);
          bfunc.eval_grad_bfunc_on_facets(quad, elem, f, grad_fbase0);

          double h_T = sqrt(domaine.carre_pas_maille(elem));
          double invh_T = 1. / h_T; // TODO regarder penalisation remplacer h_T par h_F
          double sur_f = domaine.face_surfaces(f);
          double nu_F = 0.;
          if (is_aniso_)
            {
              for (int d = 0; d < Objet_U::dimension; d++)
                nu_F += nu(elem, d) * face_normales(f, d) * face_normales(f, d);
              nu_F /= sur_f * sur_f;
            }
          else
            nu_F = nu(elem, 0);

          for (int i = 0; i < nb_bfunc; i++)
            {
              scalar_product = 0.;
              for (int k = 0; k < nb_pts_int_fac; k++)
                for (int d = 0; d < Objet_U::dimension; d++)
                  {
                    bool ori = is_aniso_ ? d : 0;
                    scalar_product(k) += nu(elem, ori) * face_normales(f, d) / sur_f * grad_fbase0(i, k, d);
                  }

              for (int j = 0; j < nb_bfunc; j++)
                {
                  for (int k = 0; k < nb_pts_int_fac; k++)
                    product(k) = fbase0(i, k) * fbase0(j, k); // TODO DG kronecker ?

                  coeff = nu_F * eta_F(f) * invh_T * quad.compute_integral_on_facet(f, product);
                  for (int d_base = 0; d_base < dim; d_base++)
                    {
                      if (mat)
                        (*mat)(ind_elem + i + d_base * nb_bfunc, ind_elem + j + d_base * nb_bfunc) += coeff;
                      secmem(elem, i + d_base * nb_bfunc) -= coeff * inco(elem, j + d_base * nb_bfunc);
                    }
                  for (int k = 0; k < nb_pts_int_fac; k++)
                    product(k) = scalar_product(k) * fbase0(j, k);

                  double integral = quad.compute_integral_on_facet(f, product);
                  for (int d_base = 0; d_base < dim; d_base++)
                    {
                      if (mat)
                        {
                          (*mat)(ind_elem + i + d_base * nb_bfunc, ind_elem + j + d_base * nb_bfunc) -= integral;
                          (*mat)(ind_elem + j + d_base * nb_bfunc, ind_elem + i + d_base * nb_bfunc) -= integral;
                        }
                      secmem(elem, i + d_base * nb_bfunc) += integral * inco(elem, j + d_base * nb_bfunc);
                      secmem(elem, j + d_base * nb_bfunc) += integral * inco(elem, i + d_base * nb_bfunc);
                    }
                }
            }
        }
    }
  contribuer_au_second_membre(secmem); // TODO DG a integrer proprement dans la boucle
}

void Op_Diff_DG_Elem::dimensionner_termes_croises(Matrice_Morse& matrice, const Probleme_base& autre_pb, int nl, int nc) const
{
  // TODO pour problemes croises
  throw;
}

void Op_Diff_DG_Elem::ajouter_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, DoubleTab& resu) const
{
  throw;
  // TODO idem above
}

void Op_Diff_DG_Elem::contribuer_termes_croises(const DoubleTab& inco, const Probleme_base& autre_pb, const DoubleTab& autre_inco, Matrice_Morse& matrice) const
{
  // TODO idem above
  throw;
}

/**
 * @brief Adds boundary condition contributions to the right-hand side.
 *
 * @details Loops over all boundary faces and, depending on the boundary condition type,
 * adds the appropriate weak enforcement term to resu:
 *
 *  - **Neumann / Neumann_paroi**: the imposed flux g_N is integrated against each basis
 *    function on the boundary face:
 *      resu(elem, i) += integral of phi_i * g_N
 *    Point-wise flux values are retrieved either from a Champ_front_var_instationnaire
 *    (using valeur_au_temps_et_au_point) or directly from Neumann::flux_impose().
 *
 *  - **Dirichlet** (weak SIP enforcement): the boundary value g_D enters two terms:
 *      - Symmetry:  -integral of { nu * grad(phi_i) } . n * g_D
 *      - Penalty:   (nu_F * eta_F / h_T) * integral of phi_i * g_D
 *    The Dirichlet value is retrieved either point-wise (Champ_front_var_instationnaire)
 *    or face-wise (Dirichlet::val_imp_au_temps()).
 *
 *  - **Dirichlet_homogene**: no contribution (homogeneous condition, nothing to add).
 *
 * Both isotropic and anisotropic diffusivities are handled for the Dirichlet terms.
 * The use of Champ_front_softanalytique is explicitly forbidden and triggers an error.
 *
 * @param resu  The right-hand side array to accumulate boundary contributions into.
 */
void Op_Diff_DG_Elem::contribuer_au_second_membre(DoubleTab& resu) const
{

  const Domaine_DG& domaine = le_dom_dg_.valeur();

  int nb_bords = domaine.nb_front_Cl();

  const IntTab& face_voisins = domaine.face_voisins();
  const DoubleTab& face_normales = domaine.face_normales();

  const Nom& nom_inco = equation().inconnue().le_nom();
  int order = Option_DG::Get_order_for(nom_inco);
  int dim = equation().inconnue().is_vectorial() ? Objet_U::dimension : 1;

  const BasisFunction& bfunc = le_dom_dg_->get_basisFunction(order);
  const int nb_bfunc = bfunc.nb_bfunc();

  assert(dim*nb_bfunc == equation().inconnue().valeurs().line_size());

  const int quad_order = bfunc.get_default_quadrature_order();
  const Quadrature_base& quad = domaine.get_quadrature(quad_order);
  const DoubleTab& integ_points_facets = quad.get_integ_points_facets();
  int nb_pts_int_fac = integ_points_facets.dimension(1);

  DoubleTab fbase(nb_bfunc, nb_pts_int_fac);
  DoubleTab grad_fbase(nb_bfunc, nb_pts_int_fac, Objet_U::dimension);
  DoubleTab scalar_product_dim(dim, nb_pts_int_fac); // DoubleTab used for storing scalar products of the RHS and the basis functions in x, y, (z)
  DoubleTab scalar_product(nb_pts_int_fac);          // DoubleTab used for reftab scalar_product_dim for a given dimension
  // Les conditions aux limites pour le second membre
  const DoubleTab& eta_F = bfunc.get_eta_facet(); // Compute the penalisation coefficient
  int ind_face;
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = la_zcl_dg_->les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      int num1f = 0;
      int num2f = le_bord.nb_faces();
      double xk = 0., yk = 0., zk = 0.;
      double temps = equation().schema_temps().temps_courant();
      double sur_f = 0.;

      if (sub_type(Champ_front_softanalytique, la_cl.valeur().champ_front()))
        {
          Cerr << "You have to use a Champ_front_fonc_txyz and not " << la_cl.valeur().champ_front().que_suis_je() << finl;
          exit();
        }

      bool avec_valeur_aux_points = false;
      if (sub_type(Champ_front_var_instationnaire, la_cl.valeur().champ_front()))
        {
          const Champ_front_var_instationnaire& ch_txyz = ref_cast(Champ_front_var_instationnaire, la_cl.valeur().champ_front());
          avec_valeur_aux_points = ch_txyz.valeur_au_temps_et_au_point_disponible();
        }

      if (sub_type(Neumann, la_cl.valeur()))
        {
          if (avec_valeur_aux_points)
            {
              if (sub_type(Champ_front_var_instationnaire, la_cl.valeur().champ_front()))
                {
                  const Champ_front_var_instationnaire& champ_front =
                    ref_cast(Champ_front_var_instationnaire, la_cl.valeur().champ_front());

                  for (int ind_faceb = num1f; ind_faceb < num2f; ind_faceb++)
                    {
                      ind_face = le_bord.num_face(ind_faceb);
                      int elem = face_voisins(ind_face, 0); // The cell that have one facet on the boundary

                      bfunc.eval_bfunc_on_facets(quad, elem, ind_face, fbase);

                      for (int i = 0; i < nb_bfunc; i++)
                        {
                          scalar_product_dim = 0.;
                          for (int k = 0; k < nb_pts_int_fac; k++)
                            {

                              // Coordonnees des points d'integration
                              xk = integ_points_facets(ind_face, k, 0);
                              yk = integ_points_facets(ind_face, k, 1);
                              if (dimension == 3)
                                zk = integ_points_facets(ind_face, k, 2);
                              for (int d_base = 0; d_base < dim; d_base++)
                                {
                                  double flux_impose_k = champ_front.valeur_au_temps_et_au_point(temps, 0, xk, yk, zk, d_base);
                                  scalar_product_dim(d_base, k) = fbase(i, k) * flux_impose_k;
                                }
                            }
                          for (int d_base = 0; d_base < dim; d_base++)
                            {
                              scalar_product.ref_array(scalar_product_dim, d_base*nb_pts_int_fac, nb_pts_int_fac);
                              resu(elem, i + d_base * nb_bfunc) += quad.compute_integral_on_facet(ind_face, scalar_product);
                            }
                        }
                    }
                }
            }
          else
            {
              const Neumann& neumann = ref_cast(Neumann, la_cl.valeur());

              for (int ind_faceb = num1f; ind_faceb < num2f; ind_faceb++)
                {
                  ind_face = le_bord.num_face(ind_faceb);
                  int elem = face_voisins(ind_face, 0); // The cell that have one facet on the boundary

                  bfunc.eval_bfunc_on_facets(quad, elem, ind_face, fbase);

                  for (int i = 0; i < nb_bfunc; i++)
                    {
                      scalar_product_dim = 0.;
                      for (int k = 0; k < nb_pts_int_fac; k++)
                        {
                          for (int d_base = 0; d_base < dim; d_base++)
                            {
                              double flux_impose_k = neumann.flux_impose(ind_faceb, d_base);
                              scalar_product_dim(d_base, k) = fbase(i, k) * flux_impose_k;
                            }
                        }
                      for (int d_base = 0; d_base < dim; d_base++)
                        {
                          scalar_product.ref_array(scalar_product_dim, d_base*nb_pts_int_fac, nb_pts_int_fac);
                          resu(elem, i + d_base * nb_bfunc) += quad.compute_integral_on_facet(ind_face, scalar_product);
                        }
                    }
                }
            }
        }
      else if (sub_type(Dirichlet_homogene, la_cl.valeur()))
        {
          // On ne fait rien et c'est normal
        }
      else if (sub_type(Dirichlet, la_cl.valeur()))
        {

          if (avec_valeur_aux_points)
            {
              if (sub_type(Champ_front_var_instationnaire, la_cl.valeur().champ_front()))
                {
                  const Champ_front_var_instationnaire& champ_front =
                    ref_cast(Champ_front_var_instationnaire, la_cl.valeur().champ_front());

                  for (int ind_faceb = num1f; ind_faceb < num2f; ind_faceb++)
                    {

                      ind_face = le_bord.num_face(ind_faceb);

                      int elem = face_voisins(ind_face, 0); // The cell that have one facet on the boundary

                      bfunc.eval_bfunc_on_facets(quad, elem, ind_face, fbase);
                      bfunc.eval_grad_bfunc_on_facets(quad, elem, ind_face, grad_fbase);

                      sur_f = domaine.face_surfaces(ind_face);

                      double h_T = sqrt(domaine.carre_pas_maille(elem));
                      double invh_T = 1. / h_T; // TODO regarder penalisation remplacer h_T par h_F
                      double nu_F = 0.;
                      if (is_aniso_)
                        {
                          for (int d = 0; d < Objet_U::dimension; d++)
                            nu_F += nu(elem, d) * face_normales(ind_face, d) * face_normales(ind_face, d);
                          nu_F /= sur_f * sur_f;
                        }
                      else
                        nu_F = nu(elem, 0);

                      for (int i = 0; i < nb_bfunc; i++)
                        {
                          double u_bord_k = 0.;
                          scalar_product_dim = 0.;
                          for (int k = 0; k < nb_pts_int_fac; k++)
                            {
                              // Coordonnees des points d'integration
                              xk = integ_points_facets(ind_face, k, 0);
                              yk = integ_points_facets(ind_face, k, 1);
                              if (dimension == 3)
                                zk = integ_points_facets(ind_face, k, 2);

                              for (int d_base = 0; d_base < dim; d_base++)
                                {
                                  u_bord_k = champ_front.valeur_au_temps_et_au_point(temps, 0, xk, yk, zk, d_base);
                                  for (int d = 0; d < Objet_U::dimension; d++)
                                    {
                                      bool ori = is_aniso_ ? d : 0;
                                      scalar_product_dim(d_base, k) -= nu(elem, ori) * face_normales(ind_face, d) / sur_f * grad_fbase(i, k, d) * u_bord_k;
                                    }
                                  scalar_product_dim(d_base, k) += nu_F * eta_F(ind_face) * invh_T * u_bord_k * fbase(i, k); // \eta/H_F \int g \vvec_h
                                }
                            }
                          for (int d_base = 0; d_base < dim; d_base++)
                            {
                              scalar_product.ref_array(scalar_product_dim, d_base*nb_pts_int_fac, nb_pts_int_fac);
                              resu(elem, i + d_base * nb_bfunc) += quad.compute_integral_on_facet(ind_face, scalar_product);
                            }
                        }
                    }
                }
            }
          else
            {
              const Dirichlet& dirichlet = ref_cast(Dirichlet, la_cl.valeur());

              for (int ind_faceb = num1f; ind_faceb < num2f; ind_faceb++)
                {

                  ind_face = le_bord.num_face(ind_faceb);

                  int elem = face_voisins(ind_face, 0); // The cell that have one facet on the boundary

                  bfunc.eval_bfunc_on_facets(quad, elem, ind_face, fbase);
                  bfunc.eval_grad_bfunc_on_facets(quad, elem, ind_face, grad_fbase);

                  sur_f = domaine.face_surfaces(ind_face);

                  double h_T = sqrt(domaine.carre_pas_maille(elem));
                  double invh_T = 1. / h_T; // TODO regarder penalisation remplacer h_T par h_F
                  double nu_F = 0.;
                  if (is_aniso_)
                    {
                      for (int d = 0; d < Objet_U::dimension; d++)
                        nu_F += nu(elem, d) * face_normales(ind_face, d) * face_normales(ind_face, d);
                      nu_F /= sur_f * sur_f;
                    }
                  else
                    nu_F = nu(elem, 0);

                  for (int i = 0; i < nb_bfunc; i++)
                    {
                      scalar_product_dim = 0.;
                      for (int d_base = 0; d_base < dim; d_base++)
                        {
                          double u_bord = dirichlet.val_imp_au_temps(temps, ind_faceb, d_base);
                          for (int k = 0; k < nb_pts_int_fac; k++)
                            {
                              for (int d = 0; d < Objet_U::dimension; d++)
                                {
                                  bool ori = is_aniso_ ? d : 0;
                                  scalar_product_dim(d_base, k) -= nu(elem, ori) * face_normales(ind_face, d) / sur_f * grad_fbase(i, k, d) * u_bord;
                                }
                              scalar_product_dim(d_base, k) += nu_F * eta_F(ind_face) * invh_T * u_bord * fbase(i, k); // \eta/H_F \int g \vvec_h
                            }
                          scalar_product.ref_array(scalar_product_dim, d_base*nb_pts_int_fac, nb_pts_int_fac);
                          resu(elem, i + d_base * nb_bfunc) += quad.compute_integral_on_facet(ind_face, scalar_product);
                        }
                    }
                }
            }
        }
      else
        {
        }
    }
}
