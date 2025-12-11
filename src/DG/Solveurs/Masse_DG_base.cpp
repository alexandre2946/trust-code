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

void Masse_DG_base::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Nom& nom_inco = equation().inconnue().le_nom();
  if (!matrices.count(nom_inco.getString())) return; //rien a faire

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
      for (int e = 0; e < nb_elem_tot; e++)
        {
          for (int d = 0; d<dim; d++)
            {
              for (int i = 0; i < nb_bfunc; i++ )
                for (int j = 0; j < nb_bfunc; j++ )
                  {
                    indice(((e*dim+d)*nb_bfunc+i)*nb_bfunc+j, 0) = current_indice+i;
                    indice(((e*dim+d)*nb_bfunc+i)*nb_bfunc+j, 1) = current_indice+j;
                  }
              current_indice+=nb_bfunc;
            }
        }
    }
  mat.dimensionner(indice);
}

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

  assert(nb_bfunc == inco.line_size());

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
          for (int d = 0; d<dim; d++)
            {
              for (int i=0; i<nb_bfunc; i++)
                {
                  if (mat)
                    (*mat)(current_indice+i, current_indice+i) += coef[e]*volume[e] / dt;
                  secmem(e,i) += coef[e]*volume[e]*passe(e,i) / dt;
                }
              current_indice+=nb_bfunc;
            }
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
                secmem(e,d) += coef[e]*volume[e]*passe(e,d) / dt;
                current_indice+=nb_bfunc;
              }
          else
            {
              bfunc.eval_bfunc(quad, e, fbase);
              /****************************************************************/
              /* Formule de quadrature : Ern, Finite Elements II, 2021, p 71  */
              /****************************************************************/
              for (int d = 0; d<dim; d++)
                {
                  for (int i=0; i<nb_bfunc; i++)
                    {
                      for (int j=0; j<nb_bfunc; j++)
                        {
                          product = 0.;
                          for (int k = 0; k < tab_pts_integ(e) ; k++)
                            product(k) = fbase(i, k) * fbase(j, k);

                          double integral = quad.compute_integral_on_elem(e, product);

                          if (mat)
                            (*mat)(current_indice+i+d*nb_bfunc, current_indice+j+d*nb_bfunc) += coef[e]*integral / dt;
                          secmem(e,i) += coef[e]*integral*passe(e,j) / dt;
                        }
                    }
                  current_indice+=nb_bfunc;
                }
            }
        }
    }
}
