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

#include <Op_Grad_DG.h>
#include <Neumann_sortie_libre.h>
#include <Check_espace_virtuel.h>
#include <Domaine_Cl_DG.h>
#include <Champ_Elem_DG.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTTrav.h>
#include <Dirichlet.h>

Implemente_instanciable(Op_Grad_DG, "Op_Grad_DG", Operateur_Grad_base);

Sortie& Op_Grad_DG::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Op_Grad_DG::readOn(Entree& s) { return s; }

void Op_Grad_DG::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base&)
{
  const Domaine_DG& domaine_DG = ref_cast(Domaine_DG, domaine_dis);
  const Domaine_Cl_DG& domaine_Cl_DG = ref_cast(Domaine_Cl_DG, domaine_Cl_dis);
  ref_domaine = domaine_DG;
  ref_dcl = domaine_Cl_DG;
  porosite_surf.ref(equation().milieu().porosite_face());
  face_voisins.ref(domaine_DG.face_voisins());
}

void Op_Grad_DG::dimensionner(Matrice_Morse& mat) const
{
  if (has_interface_blocs())
    {
      Operateur_Grad_base::dimensionner(mat);
      return;
    }

//  const Domaine_DG& domaine_DG = ref_domaine.valeur();
//  IntTab stencil(0, 2);
//
//  for (int f = 0; f < domaine_DG.nb_faces(); f++)
//    for (int i = 0; i < 2; i++)
//      {
//        const int e = domaine_DG.face_voisins(f, i);
//        if (e < 0) continue;
//
//        stencil.append_line(f, e);
//      }
//  tableau_trier_retirer_doublons(stencil);
//  Matrix_tools::allocate_morse_matrix(domaine_DG.nb_faces_tot(), domaine_DG.nb_elem_tot(), stencil, mat);
}

DoubleTab& Op_Grad_DG::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  if (has_interface_blocs()) return Operateur_Grad_base::ajouter(inco, resu);

//  assert_espace_virtuel_vect(inco);
//  const Domaine_DG& domaine_DG = ref_domaine.valeur();
//  const Domaine_Cl_DG& domaine_Cl_DG = ref_dcl.valeur();
//  const DoubleVect& face_surfaces = domaine_DG.face_surfaces();

  resu.echange_espace_virtuel();
  return resu;
}

DoubleTab& Op_Grad_DG::calculer(const DoubleTab& inco, DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(inco, resu);
}

void Op_Grad_DG::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& la_matrice) const
{
  if (has_interface_blocs())
    {
      Operateur_Grad_base::contribuer_a_avec(inco, la_matrice);
      return;
    }

//  assert_espace_virtuel_vect(inco);
//  const Domaine_DG& domaine_DG = ref_domaine.valeur();
//  const Domaine_Cl_DG& domaine_Cl_DG = ref_dcl.valeur();
//  const DoubleVect& face_surfaces = domaine_DG.face_surfaces();

}

int Op_Grad_DG::impr(Sortie& os) const
{
  return 0;
}
