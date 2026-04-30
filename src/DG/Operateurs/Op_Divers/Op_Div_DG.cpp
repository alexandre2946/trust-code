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

#include <Op_Div_DG.h>
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
#include <Champ_front_txyz.h>
#include <Champ_front_softanalytique.h>
#include <Neumann_paroi.h>
#include <Dirichlet.h>

Implemente_instanciable(Op_Div_DG, "Op_Div_DG", Operateur_Div_base);

Sortie& Op_Div_DG::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Op_Div_DG::readOn(Entree& s) { return s; }

void Op_Div_DG::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base&)
{
  le_dom_DG = ref_cast(Domaine_DG, domaine_dis);
  le_dcl_DG = ref_cast(Domaine_Cl_DG, domaine_Cl_dis);
}

void Op_Div_DG::completer()
{
  Operateur_base::completer();
  op_diff_ = ref_cast(Op_Diff_DG_base, equation().operateur(0).l_op_base());
}

void Op_Div_DG::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{

  if (!matrices.count("vitesse"))
    return; // rien a faire
  if (semi_impl.count("vitesse"))
    return; // semi-implicite -> rien a dimensionner

  Matrice_Morse *matv = matrices.count("vitesse") ? matrices["vitesse"] : nullptr,
                 *matp = matrices.count("pression") ? matrices["pression"] : nullptr,
                  matv2, matp2;

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

  int size_row = indices_glob_elem_p(nb_elem_tot);
  int size_col = indices_glob_elem_v(nb_elem_tot);

  const IntTab& stencil_sorted = domaine.get_stencil_sorted();
  const int nb_stencil_max = stencil_sorted.dimension(1);

  int nb_indices_line;
  int row, col, indice;

  if (matv)
    {
      matv2.dimensionner(size_row, size_col, 0);

      IntVect& tabv1 = matv2.get_set_tab1();
      IntVect& tabv2 = matv2.get_set_tab2();
      DoubleVect& coeff = matv2.get_set_coeff();
      coeff = 0;
      tabv1(0) = 1;
      for (int nelem = 0; nelem < nb_elem_tot; nelem++)
        {
          nb_indices_line = 0;
          for (int k = 0; k < nb_stencil_max; k++)
            {
              if (stencil_sorted(nelem, k) < 0)
                break;
              nb_indices_line += nb_bfunc_v * dim;
            }
          for (int k = 0; k < nb_bfunc_p; k++)
            tabv1(indices_glob_elem_p(nelem) + k + 1) = nb_indices_line + tabv1(indices_glob_elem_p(nelem) + k);
        }

      matv2.dimensionner(size_row, size_col, tabv1(size_row) - 1);

      for (int nelem = 0; nelem < nb_elem_tot; nelem++)
        {
          row = tabv1[indices_glob_elem_p(nelem)] - 1;
          nb_indices_line = tabv1[indices_glob_elem_p(nelem) + 1] - tabv1[indices_glob_elem_p(nelem)];
          indice = 0;

          for (int i = 0; i < nb_bfunc_p; i++)
            {
              for (int k = 0; k < nb_stencil_max; k++)
                {
                  if (stencil_sorted(nelem, k) < 0)
                    break;
                  col = indices_glob_elem_v(stencil_sorted(nelem, k)) + 1;
                  for (int j = 0; j < nb_bfunc_v * dim; j++)
                    tabv2[row + indice + j + k * nb_bfunc_v * dim] = col + j;
                }
              indice += nb_indices_line;
            }
        }
      matv2.is_sorted_stencil();
      assert(matv2.is_sorted_stencil());
      matv->nb_colonnes() ? *matv += matv2 : *matv = matv2;
    }
  if (matp && order_v == order_p) // Stabilization term for equal order interpolation of velocity and pressure
    {
      const int size_p = indices_glob_elem_p(nb_elem_tot);

      matp2.dimensionner(size_p, size_p, 0);

      IntVect& tabp1 = matp2.get_set_tab1();
      IntVect& tabp2 = matp2.get_set_tab2();
      DoubleVect& coeffp = matp2.get_set_coeff();
      coeffp = 0.;

      tabp1(0) = 1;

      // Nombre de colonnes non nulles par ligne :
      // pour chaque ddl pression d'un élément, on couple avec tous les ddl pression
      // de l'élément courant + de son stencil.
      for (int nelem = 0; nelem < nb_elem_tot; nelem++)
        {
          nb_indices_line = 0;
          for (int k = 0; k < nb_stencil_max; k++)
            {
              if (stencil_sorted(nelem, k) < 0)
                break;
              nb_indices_line += nb_bfunc_p;
            }

          for (int i = 0; i < nb_bfunc_p; i++)
            tabp1(indices_glob_elem_p(nelem) + i + 1) = nb_indices_line + tabp1(indices_glob_elem_p(nelem) + i);
        }

      matp2.dimensionner(size_p, size_p, tabp1(size_p) - 1);

      for (int nelem = 0; nelem < nb_elem_tot; nelem++)
        {
          nb_indices_line = tabp1(indices_glob_elem_p(nelem) + 1) - tabp1(indices_glob_elem_p(nelem));

          for (int i = 0; i < nb_bfunc_p; i++)
            {
              row = tabp1(indices_glob_elem_p(nelem) + i) - 1;

              for (int k = 0; k < nb_stencil_max; k++)
                {
                  if (stencil_sorted(nelem, k) < 0)
                    break;

                  col = indices_glob_elem_p(stencil_sorted(nelem, k)) + 1;
                  for (int j = 0; j < nb_bfunc_p; j++)
                    tabp2[row + j + k * nb_bfunc_p] = col + j;
                }
            }
        }
    }
  else // no stabilization term, but we still need to dimension the matrix
    {
      matp2.dimensionner(size_row, 1);
      IntVect& tabp1 = matp2.get_set_tab1();
      tabp1 = 2;
      tabp1(0) = 1;
      IntVect& tabp2 = matp2.get_set_tab2();
      tabp2(0) = 1;
      DoubleVect& coeffp = matp2.get_set_coeff();
      coeffp = 0;
    }
  matp2.is_sorted_stencil();
  assert(matp2.is_sorted_stencil());

  matp->nb_colonnes() ? *matp += matp2 : *matp = matp2;
}

void Op_Div_DG::ajouter_blocs_ext(const DoubleTab& vit, matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  int order_v = Option_DG::Get_order_for("vitesse");
  int order_p = Option_DG::Get_order_for("pression");

  Matrice_Morse *matv = matrices.count("vitesse") ? matrices["vitesse"] : nullptr;

  bool stabilisation = ((order_v == order_p) && matrices.count("pression")); // Calculate the stabilization term if the same order is used for velocity and pressure, and if the matrix for pressure is allocated
  Matrice_Morse *matp = matrices.count("pression") ? (stabilisation ? matrices["pression"] : nullptr) : nullptr;
  const DoubleTab& inco_p = semi_impl.count("pression") ? semi_impl.at("pression") : ref_cast(Navier_Stokes_std, equation()).pression().valeurs();

  const Domaine_DG& domaine = le_dom_DG.valeur();
  const IntTab& face_voisins = domaine.face_voisins();

  const BasisFunction& bfunc_v = domaine.get_basisFunction(order_v);
  const int nb_bfunc_v = bfunc_v.nb_bfunc();

  const BasisFunction& bfunc_p = domaine.get_basisFunction(order_p);
  const int nb_bfunc_p = bfunc_p.nb_bfunc();

  const int dim = Objet_U::dimension;
  const IntTab& indices_glob_elem_v = bfunc_v.indices_glob_elem(dim);
  const IntTab& indices_glob_elem_p = bfunc_p.indices_glob_elem();

  {
    // div part
    const int quad_order = bfunc_v.get_default_quadrature_order();
    const Quadrature_base& quad = domaine.get_quadrature(quad_order); // Same quadrature for all champs

    int nb_pts_integ_max = quad.nb_pts_integ_max();
    double coeff;
    DoubleTab grad_fbase_elem(nb_bfunc_v, nb_pts_integ_max, Objet_U::dimension);
    DoubleTab f_base_p(nb_bfunc_p, nb_pts_integ_max);
    DoubleTab scalar_product_dim(nb_pts_integ_max);

    // Loop over elements to compute \int q_h div(u_h) dV
    for (int elem = 0; elem < domaine.nb_elem(); elem++)
      {
        int ind_elem_v = indices_glob_elem_v(elem);
        int ind_elem_p = indices_glob_elem_p(elem);
        bfunc_v.eval_grad_bfunc(quad, elem, grad_fbase_elem);
        bfunc_p.eval_bfunc(quad, elem, f_base_p);
        for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
          for (int d = 0; d < dim; d++)
            for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
              {
                for (int k = 0; k < quad.nb_pts_integ(elem); k++)
                  scalar_product_dim(k) = grad_fbase_elem(velocity_index, k, d) * f_base_p(pressure_index, k);
                coeff = quad.compute_integral_on_elem(elem, scalar_product_dim);
                if (matv)
                  (*matv)(ind_elem_p + pressure_index, ind_elem_v + velocity_index + d * nb_bfunc_v) += coeff;
                secmem(elem, pressure_index) -= coeff * vit(elem, velocity_index + d * nb_bfunc_v);
              }
      }

    const DoubleTab& face_normales = domaine.face_normales();
    const DoubleVect& face_surfaces = domaine.face_surfaces();
    int nb_pts_int_fac = quad.nb_pts_integ_facets();

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

    // Loop over facets to compute \int_f [u_h.n]_F q_h dS
    for (int face = premiere_face_int; face < domaine.nb_faces(); face++)
      {
        int elem0 = face_voisins(face, 0);
        int elem1 = face_voisins(face, 1);
        double sur_f = face_surfaces(face);
        int ind_elem0_v = indices_glob_elem_v(elem0);
        int ind_elem1_v = indices_glob_elem_v(elem1);
        int ind_elem0_p = indices_glob_elem_p(elem0);
        int ind_elem1_p = indices_glob_elem_p(elem1);
        bfunc_v.eval_bfunc_on_facets(quad, elem0, face, f_base_v0);
        bfunc_v.eval_bfunc_on_facets(quad, elem1, face, f_base_v1);
        bfunc_p.eval_bfunc_on_facets(quad, elem0, face, f_base_p0);
        bfunc_p.eval_bfunc_on_facets(quad, elem1, face, f_base_p1);
        for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
          {
            for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
              {
                for (int d = 0; d < Objet_U::dimension; d++)
                  {
                    for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                      {
                        eval_jump_on_facet00(k) = -f_base_v0(velocity_index, k) * face_normales(face, d) * 0.5 * f_base_p0(pressure_index, k) / sur_f;
                        eval_jump_on_facet01(k) = +f_base_v1(velocity_index, k) * face_normales(face, d) * 0.5 * f_base_p0(pressure_index, k) / sur_f;
                        eval_jump_on_facet10(k) = -f_base_v0(velocity_index, k) * face_normales(face, d) * 0.5 * f_base_p1(pressure_index, k) / sur_f;
                        eval_jump_on_facet11(k) = +f_base_v1(velocity_index, k) * face_normales(face, d) * 0.5 * f_base_p1(pressure_index, k) / sur_f;
                      }
                    coeff00 = quad.compute_integral_on_facet(face, eval_jump_on_facet00);
                    coeff01 = quad.compute_integral_on_facet(face, eval_jump_on_facet01);
                    coeff10 = quad.compute_integral_on_facet(face, eval_jump_on_facet10);
                    coeff11 = quad.compute_integral_on_facet(face, eval_jump_on_facet11);

                    if (matv)
                      {
                        (*matv)(ind_elem0_p + pressure_index, ind_elem0_v + velocity_index + d * nb_bfunc_v) += coeff00;
                        (*matv)(ind_elem0_p + pressure_index, ind_elem1_v + velocity_index + d * nb_bfunc_v) += coeff01;
                        (*matv)(ind_elem1_p + pressure_index, ind_elem0_v + velocity_index + d * nb_bfunc_v) += coeff10;
                        (*matv)(ind_elem1_p + pressure_index, ind_elem1_v + velocity_index + d * nb_bfunc_v) += coeff11;
                      }
                    secmem(elem0, pressure_index) -= coeff00 * vit(elem0, velocity_index + d * nb_bfunc_v);
                    secmem(elem0, pressure_index) -= coeff01 * vit(elem1, velocity_index + d * nb_bfunc_v);
                    secmem(elem1, pressure_index) -= coeff10 * vit(elem0, velocity_index + d * nb_bfunc_v);
                    secmem(elem1, pressure_index) -= coeff11 * vit(elem1, velocity_index + d * nb_bfunc_v);
                  }
              }
          }
      }

    // Treatment of the boundary conditions
    for (int face = 0; face < premiere_face_int; face++)
      {
        int elem = face_voisins(face, 0); // The cell that have one facet on the boundary
        double sur_f = face_surfaces(face);

        int ind_elem_v = indices_glob_elem_v(elem);
        int ind_elem_p = indices_glob_elem_p(elem);
        bfunc_v.eval_bfunc_on_facets(quad, elem, face, f_base_v0);
        bfunc_p.eval_bfunc_on_facets(quad, elem, face, f_base_p0);
        for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
          {
            for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
              {
                for (int d = 0; d < Objet_U::dimension; d++)
                  {
                    //eval_jump_on_facet00 = 0.;
                    for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                      eval_jump_on_facet00(k) = -f_base_v0(velocity_index, k) * face_normales(face, d) * f_base_p0(pressure_index, k) / sur_f;
                    coeff00 = quad.compute_integral_on_facet(face, eval_jump_on_facet00);
                    if (matv)
                      (*matv)(ind_elem_p + pressure_index, ind_elem_v + velocity_index + d * nb_bfunc_v) += coeff00;
                    secmem(elem, pressure_index) -= coeff00 * vit(elem, velocity_index + d * nb_bfunc_v);
                  }
              }
          }
      }

    DoubleTab u_bord_k(nb_pts_int_fac, dim); // Dirichlet projection
    const DoubleTab& integ_points_facets = quad.get_integ_points_facets();
    for (int num_cl = 0; num_cl < le_dom_DG->nb_front_Cl(); num_cl++)
      {
        const Cond_lim& la_cl = le_dcl_DG->les_conditions_limites(num_cl);
        const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
        int num1f = 0;
        int num2f = le_bord.nb_faces();
        double xk = 0., yk = 0., zk = 0.;
        double temps = equation().schema_temps().temps_courant();

        if (sub_type(Champ_front_softanalytique, la_cl.valeur().champ_front()))
          {
            Cerr << " Il faut utiliser Champ_front_fonc_txyz et non " << la_cl.valeur().champ_front().que_suis_je() << finl;
            exit();
          }

        bool avec_valeur_aux_points = false;
        if (sub_type(Champ_front_var_instationnaire, la_cl.valeur().champ_front()))
          {
            const Champ_front_var_instationnaire& ch_txyz = ref_cast(Champ_front_var_instationnaire, la_cl.valeur().champ_front());
            avec_valeur_aux_points = ch_txyz.valeur_au_temps_et_au_point_disponible();
          }

        if (sub_type(Dirichlet, la_cl.valeur()))
          {
            if (avec_valeur_aux_points)
              {
                if (sub_type(Champ_front_var_instationnaire, la_cl.valeur().champ_front()))
                  {
                    const Champ_front_var_instationnaire& champ_front =
                      ref_cast(Champ_front_var_instationnaire, la_cl.valeur().champ_front());

                    for (int ind_faceb = num1f; ind_faceb < num2f; ind_faceb++)
                      {
                        u_bord_k = 0.;
                        int face = le_bord.num_face(ind_faceb);
                        int elem = face_voisins(face, 0); // The cell that have one facet on the boundary
                        double sur_f = face_surfaces(face);

                        bfunc_v.eval_bfunc_on_facets(quad, elem, face, f_base_v0);
                        bfunc_p.eval_bfunc_on_facets(quad, elem, face, f_base_p0);

                        for (int k = 0; k < nb_pts_int_fac; k++)
                          {
                            // Coordonnees des points d'integration
                            xk = integ_points_facets(face, k, 0);
                            yk = integ_points_facets(face, k, 1);
                            if (dimension == 3)
                              zk = integ_points_facets(face, k, 2);

                            for (int d = 0; d < Objet_U::dimension; d++)
                              u_bord_k(k, d) = champ_front.valeur_au_temps_et_au_point(temps, 0, xk, yk, zk, d);
                          }

                        for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
                          {
                            for (int d = 0; d < Objet_U::dimension; d++)
                              {
                                //eval_jump_on_facet01 = 0.;
                                for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                                  eval_jump_on_facet01(k) = +u_bord_k(k, d) * face_normales(face, d) * f_base_p0(pressure_index, k) / sur_f;
                                coeff01 = quad.compute_integral_on_facet(face, eval_jump_on_facet01);
                                secmem(elem, pressure_index) -= coeff01;
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
                    int face = le_bord.num_face(ind_faceb);
                    int elem = face_voisins(face, 0); // The cell that have one facet on the boundary
                    double sur_f = face_surfaces(face);

                    bfunc_v.eval_bfunc_on_facets(quad, elem, face, f_base_v0);
                    bfunc_p.eval_bfunc_on_facets(quad, elem, face, f_base_p0);

                    for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
                      {
                        for (int d = 0; d < Objet_U::dimension; d++)
                          {
                            //eval_jump_on_facet01 = 0.;
                            double u_bord = dirichlet.val_imp_au_temps(temps, ind_faceb, d);
                            for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                              eval_jump_on_facet01(k) = +u_bord * face_normales(face, d) * f_base_p0(pressure_index, k) / sur_f;
                            coeff01 = quad.compute_integral_on_facet(face, eval_jump_on_facet01);
                            secmem(elem, pressure_index) -= coeff01;
                          }
                      }
                  }
              }
          }
      }
  }


  if (!matp) return; // No matrix allocated for the stabilization term, we skip the calculation
  {
    //stabilization part
    op_diff_->update_nu();

    int premiere_face_int = domaine.premiere_face_int();

    const DoubleVect& face_surfaces = domaine.face_surfaces();
    const int quad_order = bfunc_p.get_default_quadrature_order();
    const Quadrature_base& quad = domaine.get_quadrature(quad_order); // pressure quad

    int nb_pts_int_fac = quad.nb_pts_integ_facets();
    DoubleTab f_base_p0(nb_bfunc_p, nb_pts_int_fac);
    DoubleTab f_base_p1(nb_bfunc_p, nb_pts_int_fac);
    DoubleTab eval_jump_on_facet00(nb_pts_int_fac);
    DoubleTab eval_jump_on_facet01(nb_pts_int_fac);
    DoubleTab eval_jump_on_facet10(nb_pts_int_fac);
    DoubleTab eval_jump_on_facet11(nb_pts_int_fac);
    double coeff00, coeff01, coeff11;

    // Loop over facets to compute \int_f [u_h.n]_F q_h dS
    for (int face = premiere_face_int; face < domaine.nb_faces(); face++)
      {
        int elem0 = face_voisins(face, 0);
        int elem1 = face_voisins(face, 1);
        double sur_f = face_surfaces(face);
        int ind_elem0_p = indices_glob_elem_p(elem0);
        int ind_elem1_p = indices_glob_elem_p(elem1);
        bfunc_p.eval_bfunc_on_facets(quad, elem0, face, f_base_p0);
        bfunc_p.eval_bfunc_on_facets(quad, elem1, face, f_base_p1);

        double nu0 = op_diff_->nu(elem0, 0);
        double nu1 = op_diff_->nu(elem1, 0);
        double nu_f = 2. * nu0 * nu1 / (nu0 + nu1); // harmonic mean

        for (int pressure_index_l = 0; pressure_index_l < nb_bfunc_p; pressure_index_l++) // Loop for basis function
          {
            for (int pressure_index_r = 0; pressure_index_r < nb_bfunc_p; pressure_index_r++) // Loop for test function
              {
                // elem 0 with elem 0
                for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                  eval_jump_on_facet00(k) = nu_f * f_base_p0(pressure_index_l, k) * f_base_p0(pressure_index_r, k) * sur_f;
                coeff00 = quad.compute_integral_on_facet(face, eval_jump_on_facet00);
                (*matp)(ind_elem0_p + pressure_index_l, ind_elem0_p + pressure_index_r) += coeff00;
                secmem(elem0, pressure_index_l) -= coeff00 * inco_p(elem0, pressure_index_r);

                // elem 0 with elem 1
                for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                  eval_jump_on_facet01(k) = -nu_f * f_base_p0(pressure_index_l, k) * f_base_p1(pressure_index_r, k) * sur_f;
                coeff01 = quad.compute_integral_on_facet(face, eval_jump_on_facet01);
                (*matp)(ind_elem0_p + pressure_index_l, ind_elem1_p + pressure_index_r) += coeff01;
                secmem(elem0, pressure_index_l) -= coeff01 * inco_p(elem1, pressure_index_r);
                // elem 1 with elem 0
                (*matp)(ind_elem1_p + pressure_index_r, ind_elem0_p + pressure_index_l) += coeff01;
                secmem(elem1, pressure_index_r) -= coeff01 * inco_p(elem0, pressure_index_l);

                // elem 1 with elem 1
                for (int k = 0; k < quad.nb_pts_integ_facets(); k++)
                  eval_jump_on_facet11(k) = nu_f * f_base_p1(pressure_index_l, k) * f_base_p1(pressure_index_r, k) * sur_f;
                coeff11 = quad.compute_integral_on_facet(face, eval_jump_on_facet11);
                (*matp)(ind_elem1_p + pressure_index_l, ind_elem1_p + pressure_index_r) += coeff11;
                secmem(elem1, pressure_index_l) -= coeff11 * inco_p(elem1, pressure_index_r);
              }
          }
      }
  }
}

DoubleTab& Op_Div_DG::calculer(const DoubleTab& vit, DoubleTab& div) const
{
  div = 0.;
  return ajouter(vit, div);
}

int Op_Div_DG::impr(Sortie& os) const
{
  return 1;
}

void Op_Div_DG::volumique(DoubleTab& div) const
{
  const Domaine_DG& domaine_DG = le_dom_DG.valeur();
  const DoubleVect& vol = domaine_DG.volumes();
  const int nb_elem = domaine_DG.domaine().nb_elem_tot();

  for (int num_elem = 0; num_elem < nb_elem; num_elem++)
    div(num_elem, 0) /= vol(num_elem); // TODO DFG ici c'est n'importe quoi, trouve comment avoir une valeur coherente !!!
}
