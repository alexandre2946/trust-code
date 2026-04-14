/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Champ_Vect_Elem_Coloc.h>
//#include <Domaine_Cl_dis_base.h>
//#include <Dirichlet_homogene.h>
//#include <Neumann_homogene.h>
//#include <Pb_Euler.h>
//#include <Periodique.h>
//#include <Dirichlet.h>
//#include <Symetrie.h>
//#include <Neumann.h>
//#include <Navier.h>

#include <Domaine_Cl_dis_base.h>
#include <Dirichlet_homogene.h>
#include <Neumann_homogene.h>
#include <Periodique.h>
#include <Dirichlet.h>
#include <Symetrie.h>
#include <Neumann.h>
#include <Navier.h>

Implemente_instanciable(Champ_Vect_Elem_Coloc, "Champ_Vect_Elem_Coloc", Champ_Elem_Coloc);
Sortie& Champ_Vect_Elem_Coloc::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }
Entree& Champ_Vect_Elem_Coloc::readOn(Entree& s) { return s; }

//void Champ_Vect_Elem_Coloc::init_fcl() const
//{
//  Champ_Inc_P0_base::init_fcl();
//
//}

void Champ_Vect_Elem_Coloc::init_fcl() const
{
  const Domaine_VF& domaine = ref_cast(Domaine_VF,le_dom_VF.valeur());
  const Conds_lim& cls = domaine_Cl_dis().les_conditions_limites();
  int i, f, n;

  fcl_.resize(domaine.nb_faces_tot(), 3);
  for (n = 0; n < cls.size(); n++)
    {
      const Front_VF& fvf = ref_cast(Front_VF, cls[n]->frontiere_dis());
      int idx = sub_type(Neumann, cls[n].valeur())
                + 2 * sub_type(Navier, cls[n].valeur())
                + 3 * sub_type(Dirichlet, cls[n].valeur()) + 3 * sub_type(Neumann_homogene, cls[n].valeur())
                + 4 * sub_type(Dirichlet_homogene, cls[n].valeur())
                + 5 * sub_type(Periodique, cls[n].valeur());
      if (!idx)
        {
          Cerr << "Champ_Vect_Elem_Coloc : CL non codee rencontree!" << finl;
          Process::exit();
        }
      for (i = 0; i < fvf.nb_faces_tot(); i++)
        f = fvf.num_face(i), fcl_(f, 0) = idx, fcl_(f, 1) = n, fcl_(f, 2) = i;
    }
  fcl_init_ = 1;
}



//
//const Domaine_Coloc& Champ_Vect_Elem_Coloc::domaine_Coloc() const
//{
//  return ref_cast(Domaine_Coloc, le_dom_VF.valeur());
//}
//
//int Champ_Vect_Elem_Coloc::nb_valeurs_nodales() const
//{
//  return domaine_Coloc().nb_elem()*nb_compo_; //// pas nb d elem
//}
//
//int Champ_Vect_Elem_Coloc::reprendre(Entree& fich)
//{
//  if (! via_ch_fonc_reprise()) return Champ_Inc_base::reprendre(fich); /* ie: resume last time ! */
//  ////////////////
//  const Pb_Euler * pbm = mon_equation_non_nul() ? (sub_type(Pb_Euler, equation().probleme()) ? &ref_cast(Pb_Euler, equation().probleme()) : nullptr) : nullptr; // pourquoi cela ?????
//  if (pbm) return Champ_Inc_base::reprendre(fich);
//
//  // sinon on fait ca ...
//  const Domaine_Coloc* domaine = le_dom_VF.non_nul() ? &ref_cast( Domaine_Coloc,le_dom_VF.valeur()) : nullptr;
//  valeurs().set_md_vector(MD_Vector()); //on enleve le MD_Vector...
//  valeurs().resize(0);
//  int ret = Champ_Inc_base::reprendre(fich);
//  //et on remet le bon si on peut
//  //if (domaine) valeurs().set_md_vector(valeurs().dimension_tot(0) > domaine->nb_elem_tot() ? domaine->mdv_elems_faces : domaine->domaine().md_vector_elements()); // pas besoin de ce test
//  if (domaine) valeurs().set_md_vector(domaine->domaine().md_vector_elements());
//  return ret;
//}
//
////int Champ_Vect_Elem_Coloc::fixer_nb_valeurs_nodales(int n)
////{
////  //// Non : pas correcte ; tab dim=1
////  assert(n == domaine_dis_base().domaine().nb_elem());
////
////  // Probleme: nb_comp vaut dimension mais on ne veut qu'une dimension !!!
////  // HACK :
////  int old_nb_compo = nb_compo_;
////  if(nb_compo_ != 1) nb_compo_ /= dimension;
////
////  /* variables : valeurs normales aux faces, puis valeurs aux elements par blocs -> pour que line_size() marche */
////  creer_tableau_distribue(domaine_dis_base().domaine().md_vector_elements());
////  nb_compo_ = old_nb_compo;
////  return n;
////}
//
//
//int Champ_Vect_Elem_Coloc::fixer_nb_valeurs_nodales(int n)
//{
//  assert (n == domaine_dis_base().domaine().nb_elem()*nb_compo_); /// a modifier comp$*n
//  creer_tableau_distribue(domaine_dis_base().domaine().md_vector_elements());
//  return n;
//}



