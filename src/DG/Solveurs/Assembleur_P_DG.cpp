/****************************************************************************
* Copyright (c) 2025, CEA
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
/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Neumann_sortie_libre.h>
#include <Assembleur_P_DG.h>
#include <Domaine_Cl_DG.h>
#include <Champ_Elem_DG.h>
#include <Matrice_Diagonale.h>
#include <Matrice_Morse_Sym.h>
#include <Matrice_Bloc_Sym.h>
#include <Static_Int_Lists.h>
#include <Domaine_DG.h>
#include <TRUSTTab_parts.h>
#include <Operateur_Grad.h>
#include <Matrix_tools.h>
#include <Milieu_base.h>
#include <Array_tools.h>
#include <Dirichlet.h>
#include <Debog.h>
#include <Perf_counters.h>
#include <BasisFunction.h>
#include <Navier_Stokes_std.h>

Implemente_instanciable(Assembleur_P_DG,"Assembleur_P_DG",Assembleur_base);

Sortie& Assembleur_P_DG::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Assembleur_P_DG::readOn(Entree& s) { return Assembleur_base::readOn(s); }

int Assembleur_P_DG::assembler(Matrice& la_matrice)
{
  DoubleVect rien;
  return assembler_mat(la_matrice, rien, 1, 1);
}

int Assembleur_P_DG::assembler_rho_variable(Matrice& la_matrice, const Champ_Don_base& rho)
{
  abort();
  return 0;
}

/**
 * @brief Core routine that builds the SIP pressure Laplacian matrix.
 *
 * @details The assembly follows the same three-term SIP structure as Op_Diff_DG_Elem
 * but with unit diffusivity (no nu weighting) and acting on the pressure unknown only.
 * It proceeds in three stages:
 *
 * **1. Sparsity pattern construction**
 * tab1 (row pointers) and tab2 (column indices) are filled in two passes using the
 * pre-computed sorted face-neighbour stencil, yielding a pressure x pressure matrix
 * of size size_p x size_p where each element block couples to all its face-neighbours.
 *
 * **2. Volume stiffness term**
 * For each element:
 *   mat(i, j) += integral of grad(phi_i).grad(phi_j)
 *
 * **3. Internal face SIP terms**
 * For each internal face shared by elem0 and elem1:
 *  - *Penalty term*: (eta_F/h_T) * integral of phi_i * phi_j, added to the diagonal
 *    blocks and subtracted from the off-diagonal (cross-element) blocks.
 *  - *Consistency + symmetry terms*: 0.5 * integral of grad(phi_i).n * phi_j,
 *    assembled into all four block combinations with signs ensuring global symmetry.
 *
 * **4. Boundary face SIP terms** (Dirichlet faces: fcl flag 6 or 7)
 * The same penalty and consistency terms are applied to the single adjacent element,
 * weakly enforcing the Dirichlet pressure condition in the SIP sense.
 *
 * @param la_matrice       The output matrix (typed as Matrice_Morse).
 * @param diag             Unused diagonal coefficient vector (kept for interface compatibility).
 * @param incr_pression    If 1, the solver works on pressure increments.
 * @param resoudre_en_u    If 1, the solver works in velocity units.
 * @return Always 1.
 */
int Assembleur_P_DG::assembler_mat(Matrice& la_matrice, const DoubleVect& diag, int incr_pression, int resoudre_en_u)
{
  set_resoudre_increment_pression(incr_pression);
  set_resoudre_en_u(resoudre_en_u);
  Cerr << "[DG] Starting the pressure matrix assembly ... ";
  statistics().begin_count(STD_COUNTERS::matrix_assembly,statistics().get_last_opened_counter_level()+1);
  la_matrice.typer("Matrice_Morse");
  Matrice_Morse& mat = ref_cast(Matrice_Morse, la_matrice.valeur());

  const Domaine_DG& domaine = ref_cast(Domaine_DG, le_dom_dg_.valeur());
  const Champ_Elem_DG& ch = ref_cast(Champ_Elem_DG, ref_cast(Navier_Stokes_std, equation()).pression());

  int nordre = Option_DG::Get_order_for("pression");

  const BasisFunction& bfunc = domaine.get_basisFunction(nordre);
  const int nb_basis_func = bfunc.nb_bfunc();

  const IntTab& indices_glob_elem = bfunc.indices_glob_elem();

  int nb_elem_tot = le_dom_dg_->nb_elem_tot();
  int size_inc = indices_glob_elem(nb_elem_tot);

  const Stencil& stencil_sorted = domaine.get_stencil_sorted();
  const int nb_stencil_max = stencil_sorted.dimension(1);

  mat.dimensionner(size_inc, size_inc, 0);
  auto& tab1 = mat.get_set_tab1();
  auto& tab2 = mat.get_set_tab2();

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
          nb_indices_line += nb_basis_func;
        }
      for (int k = 0; k < nb_basis_func; k++)
        tab1(indices_glob_elem(nelem) + k + 1) = nb_indices_line + tab1(indices_glob_elem(nelem) + k);
    }

  mat.dimensionner(size_inc, tab1(size_inc) - 1);

  for (int nelem = 0; nelem < nb_elem_tot; nelem++)
    {
      row = tab1[indices_glob_elem(nelem)] - 1;
      nb_indices_line = tab1[indices_glob_elem(nelem) + 1] - tab1[indices_glob_elem(nelem)];
      indice = 0;
      for (int k = 0; k < nb_stencil_max; k++)
        {
          if (stencil_sorted(nelem,k) < 0) break;
          col = indices_glob_elem(stencil_sorted(nelem,k))+1;
          for (int j=0; j<nb_basis_func; j++)
            for (int i=0; i<nb_basis_func; i++)
              tab2[row+indice+j+nb_indices_line*i] = col+j;
          indice += nb_basis_func;
        }
    }
  mat.sort_stencil();

  const DoubleTab& eta_F = bfunc.get_eta_facet(); // Compute the penalisation coefficient

  const int quad_order = bfunc.get_default_quadrature_order();
  const Quadrature_base& quad = domaine.get_quadrature(quad_order);
  int nb_pts_integ_max = quad.nb_pts_integ_max();
  double coeff;

  DoubleTab grad_fbase_elem(nb_basis_func, nb_pts_integ_max, Objet_U::dimension);
  DoubleTab divergence(nb_pts_integ_max);

  for (int e = 0; e < le_dom_dg_->nb_elem(); e++)
    {
      bfunc.eval_grad_bfunc(quad, e, grad_fbase_elem);
      int ind_elem = indices_glob_elem(e);
      for (int i = 0; i < nb_basis_func; i++)
        for (int j = 0; j < nb_basis_func; j++)
          {
            divergence = 0.;
            for (int k = 0; k < quad.nb_pts_integ(e) ; k++)
              for (int d=0; d<Objet_U::dimension; d++)
                divergence(k) += grad_fbase_elem(i,k,d) * grad_fbase_elem(j,k,d);

            coeff = quad.compute_integral_on_elem(e, divergence);
            mat(ind_elem+i, ind_elem+j) += coeff;
          }
    }

  int nb_pts_int_fac = quad.nb_pts_integ_facets();
  const IntTab& face_voisins = domaine.face_voisins();

  int premiere_face_int = domaine.premiere_face_int();
  const DoubleTab& face_normales = domaine.face_normales();

  int elem0, elem1;

  DoubleTab product(nb_pts_int_fac);
  DoubleTab scalar_product(nb_pts_int_fac);

  DoubleTab fbase0(nb_basis_func, nb_pts_int_fac);
  DoubleTab fbase1(nb_basis_func, nb_pts_int_fac);

  DoubleTab grad_fbase0(nb_basis_func, nb_pts_int_fac, Objet_U::dimension);
  DoubleTab grad_fbase1(nb_basis_func, nb_pts_int_fac, Objet_U::dimension);

  for (int f = premiere_face_int; f < domaine.nb_faces(); f++)
    {

      elem0 = face_voisins(f, 0);
      elem1 = face_voisins(f, 1);

      int ind_elem0 = indices_glob_elem(elem0);
      int ind_elem1 = indices_glob_elem(elem1);

      double sur_f = domaine.face_surfaces(f);

      double h_T = sqrt(std::min(domaine.carre_pas_maille(elem0), domaine.carre_pas_maille(elem1))); // TODO possibility to use harmonic mean (stability)
      double invh_T = 1. / h_T;

      //*****************//
      // penalizing term //
      //*****************//
      for (int i_elem = 0; i_elem < 2; i_elem++)
        {
          int elem = face_voisins(f, i_elem);
          int ind_elem = indices_glob_elem(elem);

          bfunc.eval_bfunc_on_facets(quad, elem, f, fbase0);

          for (int i = 0; i < nb_basis_func; i++)
            for (int j = 0; j < nb_basis_func; j++)
              {
                for (int k = 0; k < nb_pts_int_fac; k++)
                  product(k) = fbase0(i, k) * fbase0(j, k); // TODO DG kronecker ?

                coeff = eta_F(f) * invh_T * quad.compute_integral_on_facet(f, product);
                mat(ind_elem + i, ind_elem + j) += coeff;
              }
        }

      // crossed_term
      bfunc.eval_bfunc_on_facets(quad, elem0, f, fbase0);
      bfunc.eval_bfunc_on_facets(quad, elem1, f, fbase1);

      for (int i = 0; i < nb_basis_func; i++)
        for (int j = 0; j < nb_basis_func; j++)
          {
            for (int k = 0; k < nb_pts_int_fac; k++)
              product(k) = fbase0(i, k) * fbase1(j, k);

            double integral = quad.compute_integral_on_facet(f, product);
            coeff =  eta_F(f) * invh_T *integral;
            mat(ind_elem0 + i, ind_elem1 + j) -= coeff;
            mat(ind_elem1 + j, ind_elem0 + i) -= coeff; //symmetry
          }

      //****************//
      // symmetric term //
      //****************//
      bfunc.eval_grad_bfunc_on_facets(quad, elem0, f, grad_fbase0);
      bfunc.eval_grad_bfunc_on_facets(quad, elem1, f, grad_fbase1);

      for (int i = 0; i < nb_basis_func; i++)
        {
          scalar_product = 0.;
          for (int k = 0; k < nb_pts_int_fac; k++)
            for (int d = 0; d < Objet_U::dimension; d++)
              scalar_product(k) += face_normales(f, d) / sur_f * grad_fbase0(i, k, d);

          for (int j = 0; j < nb_basis_func; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase0(j, k);
              double integral = quad.compute_integral_on_facet(f, product);

              mat(ind_elem0 + i, ind_elem0 + j) -= 0.5 *integral;
              mat(ind_elem0 + j, ind_elem0 + i) -= 0.5 *integral; //symmetry
            }

          for (int j = 0; j < nb_basis_func; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase1(j, k);
              double integral = quad.compute_integral_on_facet(f, product);

              mat(ind_elem0 + i, ind_elem1 + j) += 0.5 *integral;
              mat(ind_elem1 + j, ind_elem0 + i) += 0.5 *integral; //symmetry

            }

          scalar_product = 0.;
          for (int k = 0; k < nb_pts_int_fac; k++)
            for (int d = 0; d < Objet_U::dimension; d++)
              scalar_product(k) += face_normales(f, d) / sur_f * grad_fbase1(i, k, d);

          for (int j = 0; j < nb_basis_func; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase1(j, k);
              double integral = quad.compute_integral_on_facet(f, product);
              mat(ind_elem1 + i, ind_elem1 + j) += 0.5 *integral;
              mat(ind_elem1 + j, ind_elem1 + i) += 0.5 *integral; //symmetry
            }

          for (int j = 0; j < nb_basis_func; j++)
            {
              for (int k = 0; k < nb_pts_int_fac; k++)
                product(k) = scalar_product(k) * fbase0(j, k);
              double integral = quad.compute_integral_on_facet(f, product);

              mat(ind_elem1 + i, ind_elem0 + j) -= 0.5 *integral;
              mat(ind_elem0 + j, ind_elem1 + i) -= 0.5 *integral; //symmetry

            }
        }
    }

  for (int f = 0; f < premiere_face_int; f++) // For the boundary
    {

      if ((ch.fcl()(f, 0)==6)||(ch.fcl()(f, 0)==7))
        {
          int elem = face_voisins(f, 0); // The cell that have one facet on the boundary
          int ind_elem = indices_glob_elem(elem);

          bfunc.eval_bfunc_on_facets(quad, elem, f, fbase0);
          bfunc.eval_grad_bfunc_on_facets(quad, elem, f, grad_fbase0);

          double h_T = sqrt(domaine.carre_pas_maille(elem));
          double invh_T = 1. / h_T; //TODO review penalization: replace h_T by h_F
          double sur_f = domaine.face_surfaces(f);

          for (int i = 0; i < nb_basis_func; i++)
            {
              scalar_product = 0.;
              for (int k = 0; k < nb_pts_int_fac; k++)
                for (int d = 0; d < Objet_U::dimension; d++)
                  scalar_product(k) += face_normales(f, d) / sur_f * grad_fbase0(i, k, d);

              for (int j = 0; j < nb_basis_func; j++)
                {
                  for (int k = 0; k < nb_pts_int_fac; k++)
                    product(k) = fbase0(i, k) * fbase0(j, k); // TODO DG kronecker ?

                  coeff = eta_F(f) * invh_T * quad.compute_integral_on_facet(f, product);
                  mat(ind_elem + i, ind_elem + j) += coeff;

                  for (int k = 0; k < nb_pts_int_fac; k++)
                    product(k) = scalar_product(k) * fbase0(j, k);

                  double integral = quad.compute_integral_on_facet(f, product);
                  mat(ind_elem + i, ind_elem + j) -= integral;
                  mat(ind_elem + j, ind_elem + i) -= integral;
                }
            }
        }
    }


  Cerr << statistics().get_time_since_last_open(STD_COUNTERS::matrix_assembly) << " s" << finl;
  statistics().end_count(STD_COUNTERS::matrix_assembly);
  return 1;
}

/*! @brief Assembles the pressure matrix for a quasi-compressible fluid: laplacian(P) is replaced by div(grad(P)/rho).
 *
 * @param tab_rho The density field.
 * @return Always 1.
 */
int Assembleur_P_DG::assembler_QC(const DoubleTab& tab_rho, Matrice& matrice)
{
  Cerr << "[DG] Starting the pressure matrix assembly for Quasi Compressible" << finl;
  assembler(matrice);
  set_resoudre_increment_pression(1);
  set_resoudre_en_u(0);
  abort();
  Matrice_Bloc& matrice_bloc = ref_cast(Matrice_Bloc, matrice.valeur());
  Matrice_Morse_Sym& la_matrice = ref_cast(Matrice_Morse_Sym, matrice_bloc.get_bloc(0, 0).valeur());
  if ((la_matrice.get_est_definie() != 1) && (1))
    {
      Cerr << "[DG] No imposed pressure  --> P(0)=0" << finl;
//      if (je_suis_maitre())
//        la_matrice(0, 0) *= 2; //TODO dg a adapter
      la_matrice.set_est_definie(1);
    }

  Cerr << "[DG] End of pressure matrix assembly" << finl;
  return 1;
}

int Assembleur_P_DG::modifier_secmem(DoubleTab& secmem)
{
  return 1;
}

/**
 * @brief Removes the arbitrary pressure constant by pinning the minimum pressure to zero.
 *
 * @details When no Dirichlet pressure reference is imposed (has_P_ref == 0), the
 * pressure field is defined only up to a constant. This method makes the solution
 * unique by subtracting the global minimum pressure value (computed consistently
 * across MPI processes via mp_min) from all elements, then exchanges ghost values.
 * If a pressure reference is already set (has_P_ref == 1), the method is a no-op.
 *
 * @param pression The pressure field to modify in-place.
 * @return Always 1.
 */
int Assembleur_P_DG::modifier_solution(DoubleTab& pression)
{
  // Projection :
  double press_0;
  if(!has_P_ref)
    {
      // Use the minimum pressure as the reference pressure
      // to ensure the same pressure reference in sequential and parallel runs
      press_0=DMAXFLOAT;
      int nb_elem=le_dom_dg_->domaine().nb_elem();
      for(int n=0; n<nb_elem; n++)
        if (pression(n,0) < press_0)
          press_0 = pression(n,0);
      press_0 = mp_min(press_0);
      for(int n=0; n<nb_elem; n++)
        pression(n,0) -=press_0;
      pression.echange_espace_virtuel();
    }
  return 1;
}

const Domaine_dis_base& Assembleur_P_DG::domaine_dis_base() const
{
  return le_dom_dg_.valeur();
}

const Domaine_Cl_dis_base& Assembleur_P_DG::domaine_Cl_dis_base() const
{
  return le_dom_Cl_dg_.valeur();
}

void Assembleur_P_DG::associer_domaine_dis_base(const Domaine_dis_base& le_dom_dis)
{
  le_dom_dg_ = ref_cast(Domaine_DG, le_dom_dis);
}

void Assembleur_P_DG::associer_domaine_cl_dis_base(const Domaine_Cl_dis_base& le_dom_Cl_dis)
{
  le_dom_Cl_dg_ = ref_cast(Domaine_Cl_DG, le_dom_Cl_dis);
}

void Assembleur_P_DG::completer(const Equation_base& Eqn)
{
  mon_equation = Eqn;
  stencil_done = 0;
}
