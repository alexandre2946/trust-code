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

#include <Op_Grad_PolyMAC_CDO_Face.h>
#include <Neumann_sortie_libre.h>
#include <Check_espace_virtuel.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Champ_Elem_PolyMAC_CDO.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTTrav.h>
#include <Dirichlet.h>

Implemente_instanciable(Op_Grad_PolyMAC_CDO_Face, "Op_Grad_PolyMAC_CDO_Face", Operateur_Grad_base);

Sortie& Op_Grad_PolyMAC_CDO_Face::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Op_Grad_PolyMAC_CDO_Face::readOn(Entree& s) { return s; }

void Op_Grad_PolyMAC_CDO_Face::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base&)
{
  const Domaine_PolyMAC_CDO& zPolyMAC_HFV = ref_cast(Domaine_PolyMAC_CDO, domaine_dis);
  const Domaine_Cl_PolyMAC_family& zclPolyMAC_HFV = ref_cast(Domaine_Cl_PolyMAC_family, domaine_Cl_dis);
  ref_domaine = zPolyMAC_HFV;
  ref_dcl = zclPolyMAC_HFV;
  porosite_surf.ref(equation().milieu().porosite_face());
  face_voisins.ref(zPolyMAC_HFV.face_voisins());
}

void Op_Grad_PolyMAC_CDO_Face::dimensionner(Matrice_Morse& mat) const
{
  if (has_interface_blocs())
    {
      Operateur_Grad_base::dimensionner(mat);
      return;
    }

  const Domaine_PolyMAC_CDO& zPolyMAC_HFV = ref_domaine.valeur();
  Stencil stencil(0, 2);

  for (int f = 0; f < zPolyMAC_HFV.nb_faces(); f++)
    for (int i = 0; i < 2; i++)
      {
        const int e = zPolyMAC_HFV.face_voisins(f, i);
        if (e < 0) continue;

        stencil.append_line(f, e);
      }
  tableau_trier_retirer_doublons(stencil);
  Matrix_tools::allocate_morse_matrix(zPolyMAC_HFV.nb_faces_tot(), zPolyMAC_HFV.nb_elem_tot(), stencil, mat);
}

DoubleTab& Op_Grad_PolyMAC_CDO_Face::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  if (has_interface_blocs()) return Operateur_Grad_base::ajouter(inco, resu);

  assert_espace_virtuel_vect(inco);
  const Domaine_PolyMAC_CDO& zPolyMAC_HFV = ref_domaine.valeur();
  const Domaine_Cl_PolyMAC_family& zclPolyMAC_HFV = ref_dcl.valeur();
  const DoubleVect& face_surfaces = zPolyMAC_HFV.face_surfaces();

  double coef;
  int n0, n1;

  // Boucle sur les bords pour traiter les conditions aux limites
  int ndeb, nfin, num_face;
  for (int n_bord = 0; n_bord < zPolyMAC_HFV.nb_front_Cl(); n_bord++)
    {

      // pour chaque Condition Limite on regarde son type
      // Si face de Dirichlet ou de Symetrie on ne fait rien
      // Si face de Neumann on calcule la contribution au terme source

      const Cond_lim& la_cl = zclPolyMAC_HFV.les_conditions_limites(n_bord);
      if (sub_type(Neumann_sortie_libre, la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_cl_typee = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face = ndeb; num_face < nfin; num_face++)
            {
              double P_imp = la_cl_typee.flux_impose(num_face - ndeb);
              n0 = face_voisins(num_face, 0);
              if (n0 != -1)
                {
                  coef = face_surfaces(num_face) * porosite_surf(num_face);
                  resu(num_face) += (coef * (P_imp - inco(n0)));
                }
              else
                {
                  n1 = face_voisins(num_face, 1);
                  coef = face_surfaces(num_face) * porosite_surf(num_face);
                  resu(num_face) += (coef * (inco(n1) - P_imp));
                }
            }
        }
      // Fin de la boucle for
    }

  // Boucle sur les faces internes
  for (num_face = zPolyMAC_HFV.premiere_face_int(); num_face < zPolyMAC_HFV.nb_faces(); num_face++)
    {
      n0 = face_voisins(num_face, 0);
      n1 = face_voisins(num_face, 1);
      coef = face_surfaces(num_face) * porosite_surf(num_face);
      resu(num_face) += coef * (inco(n1) - inco(n0));
    }
  resu.echange_espace_virtuel();
  return resu;
}

DoubleTab& Op_Grad_PolyMAC_CDO_Face::calculer(const DoubleTab& inco, DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(inco, resu);
}

void Op_Grad_PolyMAC_CDO_Face::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& la_matrice) const
{
  if (has_interface_blocs())
    {
      Operateur_Grad_base::contribuer_a_avec(inco, la_matrice);
      return;
    }

  assert_espace_virtuel_vect(inco);
  const Domaine_PolyMAC_CDO& zPolyMAC_HFV = ref_domaine.valeur();
  const Domaine_Cl_PolyMAC_family& zclPolyMAC_HFV = ref_dcl.valeur();
  const DoubleVect& face_surfaces = zPolyMAC_HFV.face_surfaces();

  double coef;
  int n0, n1;

  // Boucle sur les bords pour traiter les conditions aux limites
  int ndeb, nfin, num_face;
  for (int n_bord = 0; n_bord < zPolyMAC_HFV.nb_front_Cl(); n_bord++)
    {

      // pour chaque Condition Limite on regarde son type
      // Si face de Dirichlet ou de Symetrie on ne fait rien
      // Si face de Neumann on calcule la contribution au terme source

      const Cond_lim& la_cl = zclPolyMAC_HFV.les_conditions_limites(n_bord);
      if (sub_type(Neumann_sortie_libre, la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (num_face = ndeb; num_face < nfin; num_face++)
            {
              n0 = face_voisins(num_face, 0);
              if (n0 != -1)
                {
                  coef = face_surfaces(num_face) * porosite_surf(num_face);
                  la_matrice(num_face, n0) += coef;
                }
              else
                {
                  n1 = face_voisins(num_face, 1);
                  coef = face_surfaces(num_face) * porosite_surf(num_face);
                  la_matrice(num_face, n1) -= coef;
                }
            }
        }
      // Fin de la boucle for
    }

  // Boucle sur les faces internes
  for (num_face = zPolyMAC_HFV.premiere_face_int(); num_face < zPolyMAC_HFV.nb_faces(); num_face++)
    {
      n0 = face_voisins(num_face, 0);
      n1 = face_voisins(num_face, 1);
      coef = face_surfaces(num_face) * porosite_surf(num_face);
      la_matrice(num_face, n0) += coef;
      la_matrice(num_face, n1) -= coef;
    }
}

int Op_Grad_PolyMAC_CDO_Face::impr(Sortie& os) const
{
  return 0;
}
