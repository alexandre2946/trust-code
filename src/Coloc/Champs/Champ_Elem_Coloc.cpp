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

#include <Champ_Elem_Coloc.h>
#include <Pb_Euler.h>

Implemente_instanciable(Champ_Elem_Coloc, "Champ_Elem_Coloc", Champ_Inc_P0_base);

Sortie& Champ_Elem_Coloc::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Champ_Elem_Coloc::readOn(Entree& s)
{
  lire_donnees(s);
  return s;
}

//int Champ_Elem_Coloc::imprime(Sortie& os, int ncomp) const
//{
//  const Domaine_dis_base& domaine_dis = domaine_dis_base();
//  const Domaine& domaine = domaine_dis.domaine();
//  const DoubleTab& coord=domaine.coord_sommets();
//  const int nb_som = domaine.nb_som();
//  const DoubleTab& val = valeurs();
//  int som;
//  os << nb_som << finl;
//  for (som=0; som<nb_som; som++)
//    {
//      if (dimension==3)
//        os << coord(som,0) << " " << coord(som,1) << " " << coord(som,2) << " " ;
//      if (dimension==2)
//        os << coord(som,0) << " " << coord(som,1) << " " ;
//      if (nb_compo_ == 1)
//        os << val(som) << finl;
//      else
//        os << val(som,ncomp) << finl;
//    }
//  os << finl;
//  Cout << "Champ_Elem_Coloc::imprime FIN >>>>>>>>>> " << finl;
//  return 1;
//}

int Champ_Elem_Coloc::nb_valeurs_nodales() const
{
  return domaine_Coloc().nb_elem();
}

int Champ_Elem_Coloc::fixer_nb_valeurs_nodales(int n)
{
  assert (n == domaine_dis_base().domaine().nb_elem());
  creer_tableau_distribue(domaine_dis_base().domaine().md_vector_elements());
  return n;
}

const Domaine_Coloc& Champ_Elem_Coloc::domaine_Coloc() const
{
  return ref_cast(Domaine_Coloc, le_dom_VF.valeur());
}

//DoubleTab& Champ_Elem_Coloc::valeur_aux_faces(DoubleTab& dst) const
//{
//  const Domaine_VF& domaine = ref_cast(Domaine_VF, le_dom_VF.valeur());
//  const IntTab& f_e = domaine.face_voisins();
//  const DoubleTab& src = valeurs();
//
//  /* vals doit etre pre-dimensionne */
//  int i, e, f, n, N = (src.nb_dim() == 1 ? 1 : src.dimension(1));
//  assert(dst.dimension(0) == domaine.xv().dimension(0) && N == (dst.nb_dim() == 1 ? 1 : dst.dimension(1)));
//  // e enelver
//  if (src.dimension_tot(0) > domaine.nb_elem_tot()) //on a les valeurs aux faces
//    for (f = 0; f < dst.dimension(0); f++)
//      for (n = 0; n < N; n++) dst(f, n) = src(domaine.nb_elem_tot() + f, n);
//  ///
//  else for (f = 0; f < dst.dimension(0); f++) //on prend (amont + aval) / 2
//      for (i = 0; i < 2 && (e = f_e(f, i)) >= 0; i++)
//        for (n = 0; n < N; n++)
//          dst(f, n) += src(e, n) * (f < domaine.premiere_face_int() ? 1 : 0.5);
//
//  return dst;
//}

int Champ_Elem_Coloc::reprendre(Entree& fich)
{
  Process::exit();
  return 1;
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
}


