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

#include <Check_espace_virtuel.h>
#include <Domaine_Cl_DG.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>
#include <Op_Div_DG.h>
#include <Probleme_base.h>
#include <EcrFicPartage.h>
#include <Matrice_Morse.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <SFichier.h>
#include <Debog.h>
#include <BasisFunction.h>

Implemente_instanciable(Op_Div_DG, "Op_Div_DG", Operateur_Div_base);

Sortie& Op_Div_DG::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Op_Div_DG::readOn(Entree& s) { return s; }

void Op_Div_DG::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base&)
{
  le_dom_DG = ref_cast(Domaine_DG, domaine_dis);
  le_dcl_DG = ref_cast(Domaine_Cl_DG, domaine_Cl_dis);
}



void Op_Div_DG::ajouter_blocs([[maybe_unused]] matrices_t matrices, [[maybe_unused]] DoubleTab& secmem, const tabs_t& semi_impl) const
{

  /* const DoubleTab& inco_v = semi_impl.count("vitesse") ? semi_impl.at("vitesse") : equation().inconnue().valeurs(); // NB : is this working ?
  Matrice_Morse *mat = matrices.count("vitesse") ? matrices.at("vitesse") : nullptr; // pression for the stabilisation term if np==nv

  const Domaine_DG& domaine = le_dom_DG.valeur();
  const IntTab& face_voisins = domaine.face_voisins();

  int order_v = Option_DG::Get_order_for("vitesse");
  int order_p = Option_DG::Get_order_for("pression");

  const BasisFunction& bfunc_v = le_dom_dg_->get_basisFunction(order_v);
  const int nb_bfunc_v = bfunc_v.nb_bfunc();

  const BasisFunction& bfunc_p = le_dom_dg_->get_basisFunction(order_p);
  const int nb_bfunc_p = bfunc_p.nb_bfunc();

  const int quad_order = bfunc_v.get_default_quadrature_order();
  const Quadrature_base& quad = domaine.get_quadrature(quad_order);  // Same quadrature for all champs
  int nb_pts_integ_max = quad.nb_pts_integ_max();
  double coeff, coeff0, coeff1;

  DoubleTab Div_fbase(Objet_U::dimension, nb_bfunc_v, nb_pts_integ_max);
  DoubleTab f_base_p(nb_bfunc_p, nb_pts_integ_max);
  DoubleTab scalar_product_dim(nb_pts_integ_max);

  // Loop over elements to compute \int q_h div(u_h) dV
  for (int elem = 0; elem < domaine.nb_elem(); elem++)
    {
      int ind_elem_v = bfunc_v.indices_glob_elem_v(elem);
      int ind_elem_p = bfunc_p.indices_glob_elem_p(elem);
      bfunc_v.eval_div_bfunc(quad, elem, Div_fbase);
      bfunc_p.eval_bfunc(quad, elem, f_base_p);
      for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
        for (int d = 0; d < Objet_U::dimension; d++)
          for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
            {
              scalar_product_dim = 0.;
              for (int k = 0; k < quad.nb_pts_integ(elem); k++)
                scalar_product_dim(k) += Div_fbase(d, velocity_index, k) * f_base_p(pressure_index, k);
              coeff = quad.compute_integral_on_elem(elem, scalar_product_dim);
              if (mat)
                (*mat)(ind_elem_p + pressure_index, ind_elem_v*dim + velocity_index + d * nb_bfunc_v) += coeff;
              secmem(elem, pressure_index) -= coeff * inco_v(elem, velocity_index + d * nb_bfunc_v);
            }
    }

  const DoubleTab& face_normales = domaine.face_normales();
  int nb_pts_int_fac = quad.nb_pts_integ_facets();

  DoubleTab eval_jump_on_facet0(nb_pts_int_fac);
  DoubleTab eval_jump_on_facet1(nb_pts_int_fac);

  DoubleTab f_base_v0(nb_bfunc_v, nb_pts_int_fac);
  DoubleTab f_base_v1(nb_bfunc_v, nb_pts_int_fac);
  DoubleTab f_base_p0(nb_bfunc_p, nb_pts_int_fac);
  DoubleTab f_base_p1(nb_bfunc_p, nb_pts_int_fac);

  // Loop over facets to compute \int_f [u_h.n]_F q_h dS
  for (int face = 0; face < domaine.nb_faces(); face++)
    {
      int elem0 = face_voisins(face,0);
      int elem1 = face_voisins(face,1);
      if (elem1 != -1) // internal face
        {
          int ind_elem0_v = bfunc_v.indices_glob_elem_v(elem0);
          int ind_elem1_v = bfunc_v.indices_glob_elem_v(elem1);
          int ind_elem0_p = bfunc_p.indices_glob_elem_p(elem0);
          int ind_elem1_p = bfunc_p.indices_glob_elem_p(elem1);
          bfunc_v.eval_bfunc_on_facets(quad, elem0, face, f_base_v0);
          bfunc_v.eval_bfunc_on_facets(quad, elem1, face, f_base_v1);
          bfunc_p.eval_bfunc_on_facets(quad, elem0, face, f_base_p0);
          bfunc_p.eval_bfunc_on_facets(quad, elem1, face, f_base_p1);
          for (int pressure_index = 0; pressure_index < nb_bfunc_p; pressure_index++)
            for (int d = 0; d < Objet_U::dimension; d++)
              for (int velocity_index = 0; velocity_index < nb_bfunc_v; velocity_index++)
                {
                  eval_jump_on_facet0 = 0.;
                  eval_jump_on_facet1 = 0.;
                  for (int k = 0; k < quad.nb_pts_integ_on_facet(elem0, face); k++)
                    {
                      double mean_P = 0.5*(f_base_p0(pressure_index, k) + f_base_p1(pressure_index, k));
                      for (int dim = 0 ; dim < Objet_U::dimension ; dim++)
                        {
                          eval_jump_on_facet0(k) += f_base_v0(velocity_index, k) * face_normales(face, dim)  * mean_P;
                          eval_jump_on_facet1(k) -= f_base_v1(velocity_index, k) * face_normales(face, dim)  * mean_P;
                        }
                    }
                  coeff0 = quad.compute_integral_on_facet(elem0, face, eval_jump_on_facet0);
                  coeff1 = quad.compute_integral_on_facet(elem1, face, eval_jump_on_facet1);

                  if (mat)
                    {
                      (*mat)(ind_elem0_p + pressure_index, ind_elem0_v*dim + velocity_index + d * nb_bfunc_v) -= coeff0;
                      (*mat)(ind_elem1_p + pressure_index, ind_elem1_v*dim + velocity_index + d * nb_bfunc_v) += coeff1;
                    }
                  secmem(elem0, pressure_index) -= coeff0 * inco_v(elem0, velocity_index + d * nb_bfunc_v);
                  secmem(elem1, pressure_index) += coeff1 * inco_v(elem1, velocity_index + d * nb_bfunc_v);
                }
        }
    }

  */

}



void Op_Div_DG::dimensionner(Matrice_Morse& matrice) const
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

DoubleTab& Op_Div_DG::calculer(const DoubleTab& vit, DoubleTab& div) const
{
  div = 0.;
  return ajouter(vit, div);
}

int Op_Div_DG::impr(Sortie& os) const
{
  const int impr_bord = (le_dom_DG->domaine().bords_a_imprimer().est_vide() ? 0 : 1);
  ouvrir_fichier(Flux_div, "", je_suis_maitre());
  EcrFicPartage Flux_face;
  ouvrir_fichier_partage(Flux_face, "", impr_bord);
  const Schema_Temps_base& sch = equation().probleme().schema_temps();
  const double temps = sch.temps_courant();
  if (je_suis_maitre())
    Flux_div.add_col(temps);

  const int nb_compo = flux_bords_.dimension(1);

  // On parcours les frontieres pour sommer les flux par frontiere dans le tableau flux_bord
  DoubleVect flux_bord(nb_compo);
  DoubleVect bilan(nb_compo);
  bilan = 0.;

  for (int num_cl = 0; num_cl < le_dom_DG->nb_front_Cl(); num_cl++)
    {
      flux_bord = 0;
      const Cond_lim& la_cl = le_dcl_DG->les_conditions_limites(num_cl);
      const Front_VF& frontiere_dis = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = frontiere_dis.num_premiere_face();
      const int nfin = ndeb + frontiere_dis.nb_faces();
      for (int face = ndeb; face < nfin; face++)
        for (int k = 0; k < nb_compo; k++)
          flux_bord(k) += flux_bords_(face, k);

      for (int k = 0; k < nb_compo; k++)
        flux_bord(k) = Process::mp_sum(flux_bord(k));

      if (je_suis_maitre())
        {
          for (int k = 0; k < nb_compo; k++)
            {
              //Ajout pour impression sur fichiers separes
              Flux_div.add_col(flux_bord(k));
              bilan(k) += flux_bord(k);
            }
        }
    }

  if (je_suis_maitre())
    {
      for (int k = 0; k < nb_compo; k++)
        Flux_div.add_col(bilan(k));

      Flux_div << finl;
    }

  for (int num_cl = 0; num_cl < le_dom_DG->nb_front_Cl(); num_cl++)
    {
      const Frontiere_dis_base& la_fr = le_dcl_DG->les_conditions_limites(num_cl)->frontiere_dis();
      const Cond_lim& la_cl = le_dcl_DG->les_conditions_limites(num_cl);
      const Front_VF& frontiere_dis = ref_cast(Front_VF, la_cl->frontiere_dis());
      int ndeb = frontiere_dis.num_premiere_face();
      int nfin = ndeb + frontiere_dis.nb_faces();
      if (le_dom_DG->domaine().bords_a_imprimer().contient(la_fr.le_nom()))
        {
          Flux_face << "# Flux par face sur " << la_fr.le_nom() << " au temps " << temps << " : " << finl;
          for (int face = ndeb; face < nfin; face++)
            {
              if (dimension == 2)
                Flux_face << "# Face a x= " << le_dom_DG->xv(face, 0) << " y= " << le_dom_DG->xv(face, 1) << " flux=";
              else if (dimension == 3)
                Flux_face << "# Face a x= " << le_dom_DG->xv(face, 0) << " y= " << le_dom_DG->xv(face, 1) << " z= " << le_dom_DG->xv(face, 2) << " flux=";
              for (int k = 0; k < nb_compo; k++)
                Flux_face << " " << flux_bords_(face, k);
              Flux_face << finl;
            }
          Flux_face.syncfile();
        }
    }

  return 1;
}

void Op_Div_DG::volumique(DoubleTab& div) const
{
  const Domaine_DG& domaine_DG = le_dom_DG.valeur();
  const DoubleVect& vol = domaine_DG.volumes();
  const int nb_elem = domaine_DG.domaine().nb_elem_tot();

  for (int num_elem = 0; num_elem < nb_elem; num_elem++)
    div(num_elem) /= vol(num_elem);
}
