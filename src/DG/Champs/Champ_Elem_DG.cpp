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

#include <Champ_Elem_DG.h>
#include <TRUSTTab_parts.h>
#include <Domaine_Cl_dis_base.h>
#include <Domaine_DG.h>
#include <Option_DG.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <BasisFunction.h>


Implemente_instanciable(Champ_Elem_DG, "Champ_Elem_DG", Champ_Inc_P0_base);

Sortie& Champ_Elem_DG::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Champ_Elem_DG::readOn(Entree& s)
{
  lire_donnees(s);
  return s;
}

int Champ_Elem_DG::imprime(Sortie& os, int ncomp) const //TODO DG adapt
{
  const Domaine_dis_base& domaine_dis = domaine_dis_base();
  const Domaine& domaine = domaine_dis.domaine();
  const DoubleTab& coord=domaine.coord_sommets();
  const int nb_som = domaine.nb_som();
  const DoubleTab& val = valeurs();
  int som;
  os << nb_som << finl;
  for (som=0; som<nb_som; som++)
    {
      if (dimension==3)
        os << coord(som,0) << " " << coord(som,1) << " " << coord(som,2) << " " ;
      if (dimension==2)
        os << coord(som,0) << " " << coord(som,1) << " " ;
      if (nb_compo_ == 1)
        os << val(som) << finl;
      else
        os << val(som,ncomp) << finl;
    }
  os << finl;
  Cout << "Champ_Elem_DG::imprime FIN >>>>>>>>>> " << finl;
  return 1;
}

int Champ_Elem_DG::fixer_nb_valeurs_nodales(int n)
{
  creer_tableau_distribue(domaine_dis_base().domaine().md_vector_elements()); // TODO 26/08/2024 nb_ddl a
  return n;
}


void Champ_Elem_DG::associer_domaine_dis_base(const Domaine_dis_base& z_dis)
{
  le_dom_VF = ref_cast(Domaine_VF, z_dis);

  order_ = Option_DG::Get_order_for(nom_);// Todo regler la relation ordre inconnu/quadrature avec dictionnaire ?
  nb_bfunc_ = Option_DG::Nb_col_from_order(order_);
  if (nom_.debute_par("vitesse"))
    is_scalar_ = false;

  const int dim = is_scalar_ ? 1 : Objet_U::dimension;
  int nb_elem_tot = le_dom_VF->nb_elem_tot();

  indices_glob_elem_.resize(dim*nb_elem_tot+1);
  indices_glob_elem_(0)=0;
  for (int e = 0; e < nb_elem_tot; e++)
    indices_glob_elem_(e+1) = indices_glob_elem_(e) +  nb_bfunc_*dim;
}


/*@brief Used to project a Champ_base to a Champ_Elem_DG. Mostly used for initial condition
 *
 */
Champ_base& Champ_Elem_DG::affecter_(const Champ_base& ch)
{
  const Domaine_DG& domaine = ref_cast(Domaine_DG,le_dom_VF.valeur());
  const DoubleVect& volume = domaine.volumes();

  const BasisFunction& bfunc = domaine.get_basisFunction(order_);
  assert(nb_bfunc_ == bfunc.nb_bfunc());

  const int quad_order = bfunc.get_default_quadrature_order();

  const Quadrature_base& quad = domaine.get_quadrature(quad_order); // if Uniforme remplissage automatique
  //sinon interdit d'avoir 1

  //creation d'un DoubleTab intermediaire pour recuperer les valeurs du champ ch sur les points de quadrature ?
  const DoubleTab& integ_points = quad.get_integ_points();
  int nb_pts_integ_max = quad.nb_pts_integ_max();

  if (nom_.debute_par("vitesse"))
    is_scalar_ = false;
  const int dim = is_scalar_ ? 1: Objet_U::dimension;

  int nb_elem = domaine.nb_elem();

  DoubleTab product(nb_pts_integ_max);
  int nb_pts_integ = integ_points.dimension(0);
  DoubleTab values(nb_pts_integ,dim);
  DoubleTab phi_rhs(nb_bfunc_);
  DoubleTab res;

  DoubleTab fbase(nb_bfunc_, nb_pts_integ_max);

  ch.valeur_aux(integ_points, values);

  for (int num_elem = 0; num_elem < nb_elem; num_elem++)
    {
      Matrice_Dense invM;
      bfunc.eval_bfunc(quad, num_elem, fbase);
      if (!domaine.gram_schmidt())
        invM = bfunc.eval_invMassMatrix(quad, num_elem); //to remove and ref on global matrix
      for (int d =0; d<dim; d++)
        {
          for (int fb = 0; fb < nb_bfunc_; fb++)
            {
              for (int k = 0; k < quad.nb_pts_integ(num_elem) ; k++)
                product(k) = values(quad.ind_pts_integ(num_elem) + k,d) * fbase(fb, k);

              phi_rhs(fb) = quad.compute_integral_on_elem(num_elem, product);
            }

          if (domaine.gram_schmidt())
            {
              for (int fb = 0; fb < nb_bfunc_; fb++)
                valeurs()(num_elem,fb+d*nb_bfunc_) = phi_rhs(fb)/volume(num_elem);
            }
          else
            {
              res.ref_array(valeurs(), (num_elem*dim+d)*nb_bfunc_, nb_bfunc_);
              invM.ajouter_multvect_(phi_rhs, res);
            }
        }
    }

  valeurs().echange_espace_virtuel();
  return *this;
}
/*@brief Compute from a Champ_Elem_DG, its value at a certain position. Mostly used for post-processing
 *
 */
DoubleTab& Champ_Elem_DG::valeur_aux(const DoubleTab& positions, DoubleTab& tab_valeurs) const
{
  const Domaine_DG& domaine_DG = ref_cast(Domaine_DG,le_dom_VF.valeur());
  const Domaine& domaine = domaine_dis_base().domaine();
  throw; // /!\ TODO
  const int dim = Objet_U::dimension;

  const BasisFunction& bfunc = domaine_DG.get_basisFunction(order_);
  assert(nb_bfunc_ == bfunc.nb_bfunc());

  const int quad_order = bfunc.get_default_quadrature_order();

  const Quadrature_base& quad = domaine_DG.get_quadrature(quad_order);
  int nb_pts_integ_max = quad.nb_pts_integ_max();

  IntVect les_polys;
  les_polys.resize(tab_valeurs.dimension(0), RESIZE_OPTIONS::NOCOPY_NOINIT);

  domaine.chercher_elements(positions, les_polys); //TODO DG selectionner uniquement la premiere valeur de tab_valeurs (on refait plein de fois le même truc)

  const Champ_base& ch_base = le_champ();
  const DoubleTab& values = ch_base.valeurs();
  int nb_polys = les_polys.size();


  if (nb_polys == 0)
    return tab_valeurs;

  DoubleTab fbase(nb_bfunc_,nb_pts_integ_max);
  DoubleTab coords(nb_pts_integ_max,dim);
  for (int i = 0; i < nb_polys; i++)
    {
      int cell = les_polys(i);
      assert(cell < values.dimension_tot(0));

      if (cell != -1)
        {
          for (int j = 0; j < quad.nb_pts_integ(cell); j++)
            for (int k = 0; k<dim; k++)
              coords(j,k) = positions(quad.ind_pts_integ(cell) +j, k);
//          coords.ref_tab(positions, i*nb_points, nb_points); //pas possible car const

          bfunc.eval_bfunc(coords, cell, fbase);

          for (int j = 0; j < quad.nb_pts_integ(cell) ; j++)
            {
              tab_valeurs(i,j) = 0.;
              for (int l =0; l<nb_bfunc_; l++)
                tab_valeurs(i,j) += values(cell,l) * fbase(l,j); // reconstruction valeurs du champ aux points d'integrations
            }
        }
    }

  return tab_valeurs;
}

DoubleTab& Champ_Elem_DG::eval_elem(DoubleTab& tab_valeurs) const
{
  const Domaine_DG& domaine = ref_cast(Domaine_DG,le_dom_VF.valeur());

  const int nb_elem = domaine.nb_elem();

  const BasisFunction& bfunc = domaine.get_basisFunction(order_);
  assert(nb_bfunc_ == bfunc.nb_bfunc());

  const Quadrature_base& quad = domaine.get_quadrature(5);
  int nb_pts_integ_max = quad.nb_pts_integ_max();

  const int dim = tab_valeurs.dimension(1)/nb_pts_integ_max;

  const Champ_base& ch_base = le_champ();
  const DoubleTab& values = ch_base.valeurs();

  assert(tab_valeurs.dimension(0) == nb_elem && tab_valeurs.dimension(1) == dim*nb_pts_integ_max );

  DoubleTab fbase(nb_bfunc_,nb_pts_integ_max);
  for (int i = 0; i < nb_elem; i++)
    {
      bfunc.eval_bfunc(quad, i, fbase);

      for (int d =0; d<dim; d++)
        {

          for (int j = 0; j < quad.nb_pts_integ(i) ; j++)
            {
              tab_valeurs(i,j) = 0.;
              for (int l =0; l<nb_bfunc_; l++)
                tab_valeurs(i,j+d*nb_pts_integ_max) += values(i,l+d*nb_bfunc_) * fbase(l,j); // reconstruction valeurs du champ aux points d'integrations
            }
        }
    }

  return tab_valeurs;
}

DoubleTab& Champ_Elem_DG::valeur_aux_elems(const DoubleTab& positions, const IntVect& polys, DoubleTab& result) const
{

  const Domaine_DG& domaine = ref_cast(Domaine_DG,le_dom_VF.valeur());

  const DoubleVect& volume = domaine.volumes();

  const BasisFunction& bfunc = domaine.get_basisFunction(order_);
  assert(nb_bfunc_ == bfunc.nb_bfunc());

  const Quadrature_base& quad = domaine.get_quadrature(5);
  int nb_pts_integ_max = quad.nb_pts_integ_max();

  const Champ_base& ch_base = le_champ();
  const DoubleTab& values = ch_base.valeurs();
  int nb_polys = polys.size();

  if (nb_polys == 0)
    return result;

  // TODO : FIXME
  // For FT the resize should be done in its good position and not here ...
  if (result.nb_dim() == 1) result.resize(nb_polys, 1);

  assert(result.line_size() == 1);
  ToDo_Kokkos("critical");

  DoubleTab fbase(nb_bfunc_, nb_pts_integ_max);
  DoubleTab product(nb_pts_integ_max);

  for (int i = 0; i < nb_polys; i++)
    {
      int cell = polys(i);
      assert(cell < values.dimension_tot(0));

      if (cell != -1)
        {
          bfunc.eval_bfunc(quad, cell, fbase);

          product = 0.;
          for (int k = 0; k < quad.nb_pts_integ(cell) ; k++)
            for (int l =0; l<nb_bfunc_; l++)
              product(k) += values(cell,l) * fbase(l,k); // reconstruction valeurs du champ aux points coords

          result(i,0) = quad.compute_integral_on_elem(cell, product);
          result(i,0) /= volume(cell);
        }
    }

  return result;

}
