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

#include <Assembleur_P_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Domaine_VDF.h>
#include <Periodique.h>
#include <Symetrie.h>
#include <Neumann_sortie_libre.h>
#include <Dirichlet_entree_fluide_leaves.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Matrice_Bloc.h>
#include <Option_VDF.h>
#include <Champ_Fonc_Face_VDF.h>
#include <Matrice_Morse_Sym.h>
#include <Milieu_base.h>
#include <Matrix_tools.h>
#include <Pb_Multiphase.h>

Implemente_instanciable_sans_constructeur(Assembleur_P_VDF,"Assembleur_P_VDF",Assembleur_base);

Assembleur_P_VDF::Assembleur_P_VDF() : has_P_ref(0) { }

Sortie& Assembleur_P_VDF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom() ;
}

Entree& Assembleur_P_VDF::readOn(Entree& s )
{
  return Assembleur_base::readOn(s);
}

/*! @brief Fills the array faces with the list of indices of the periodic faces in the face_voisins array.
 *
 * @brief Each periodic face appears twice in face_voisins (each face corresponds to the opposite face).
 *  Only the one of the two with the smaller index in the list of faces of each periodic boundary
 *  is stored in the array faces.
 *  Return value:
 *  number of periodic faces (equal to the size of the faces array).
 *
 */
int Assembleur_P_VDF::liste_faces_periodiques(ArrOfInt& faces)
{
  // First, largely overestimate the array size:
  // number of boundary faces
  const int nb_faces_bord = le_dom_VDF->nb_faces_bord();
  faces.resize_array(nb_faces_bord);

  // Search for periodic faces in the boundary conditions:
  const Conds_lim& les_cl = le_dom_Cl_VDF->les_conditions_limites();
  const int nb_cl = les_cl.size();
  int nb_faces_periodiques = 0;
  for (int num_cl = 0; num_cl < nb_cl; num_cl++)
    {
      const Cond_lim_base& la_cl = les_cl[num_cl].valeur();
      // Select only the Periodique conditions
      if ( ! sub_type(Periodique,la_cl))
        continue;
      const Periodique& la_cl_perio = ref_cast(Periodique, la_cl);
      const Front_VF&    frontiere = ref_cast(Front_VF, la_cl.frontiere_dis());
      const int nb_faces_cl = frontiere.nb_faces();
      const int num_premiere_face = frontiere.num_premiere_face();
      for (int i = 0; i < nb_faces_cl; i++)
        {
          // Index of the opposite face in the boundary face array:
          const int face_associee = la_cl_perio.face_associee(i);
          if (face_associee > i)
            {
              const int num_face_global = num_premiere_face + i;
              faces[nb_faces_periodiques] = num_face_global;
              nb_faces_periodiques++;
            }
        }
    }

  // Final size of the faces array
  faces.resize_array(nb_faces_periodiques);
  return nb_faces_periodiques;
}

/*! @brief Determines the nonzero entries of the matrix and prepares the storage.
 *
 * @brief Sparse matrix of size nb_elements (rows) * nb_elem_tot (columns)
 *   Stored as a block matrix composed of two Morse matrices:
 *    * Square symmetric matrix nb_elements * nb_elements
 *      (contains the terms M(i,j) where i and j are indices of real elements)
 *    * Rectangular matrix nb_elements * (nb_elem_tot - nb_elem)
 *      (contains the terms M(i,j) where i is real and j is virtual)
 *
 */
int Assembleur_P_VDF::construire(Matrice& la_matrice)
{
  int i;
  const Domaine_VDF& domaine_vdf   = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_vdf.face_voisins();

  // Count the total number of non-zero entries:
  // square matrix: number of internal faces / 2 + nb_elem + nb periodic faces
  //                (each internal face gives one coefficient, there is one diagonal
  //                 element and each periodic face also gives one coefficient)
  // rectangular matrix: number of joint faces


  // First step: count the number of non-zero entries per row
  // For each row of the square matrix, number of non-zero entries
  const int nb_elem     = domaine_vdf.nb_elem();
  const int nb_elem_tot = domaine_vdf.nb_elem_tot();
  ArrOfInt carre_nb_non_zero(nb_elem);
  // Same for the rectangular matrix
  ArrOfInt rect_nb_non_zero(nb_elem);
  // There is one element on the diagonal:
  carre_nb_non_zero = 1;
  rect_nb_non_zero = 0;
  int carre_nb_non_zero_tot = nb_elem;
  int rect_nb_non_zero_tot = 0;

  // Plus one non-zero entry for each internal and periodic face
  // (symmetric matrix, only the entry m(line,col) with col>line is stored)

  ArrOfInt liste_faces_perio;
  const int nb_faces_periodiques = liste_faces_periodiques(liste_faces_perio);
  const int nb_faces_internes = domaine_vdf.nb_faces_internes();
  const int premiere_face_interne = domaine_vdf.premiere_face_int();
  for (i = 0; i < nb_faces_internes + nb_faces_periodiques; i++)
    {
      int face;
      if (i < nb_faces_internes) // Trick to loop over internal and periodic faces
        face = premiere_face_interne + i;
      else
        face = liste_faces_perio[i - nb_faces_internes];

      int elem0 = face_voisins(face,0);
      int elem1 = face_voisins(face,1);
      if (elem0 > elem1)
        {
          int tmp = elem1;
          elem1 = elem0;
          elem0 = tmp;
        }
      if (elem0 < nb_elem)   // elem0 is real
        {
          if (elem1 < nb_elem)      // elem1 real
            {
              carre_nb_non_zero[elem0] ++;
              carre_nb_non_zero_tot ++;
            }
          else                      // elem1 virtual
            {
              rect_nb_non_zero[elem0] ++;
              rect_nb_non_zero_tot ++;
            }
        }
    }

  // Type and size the pressure matrix
  la_matrice.typer("Matrice_Bloc");
  Matrice_Bloc& matrice =ref_cast(Matrice_Bloc , la_matrice.valeur());
  matrice.dimensionner(1,2);
  matrice.get_bloc(0,0).typer("Matrice_Morse_Sym");
  matrice.get_bloc(0,1).typer("Matrice_Morse");
  Matrice_Morse_Sym& carre = ref_cast(Matrice_Morse_Sym ,matrice.get_bloc(0,0).valeur());
  Matrice_Morse&      rect  = ref_cast(Matrice_Morse ,     matrice.get_bloc(0,1).valeur());

  carre.dimensionner(nb_elem, carre_nb_non_zero_tot);
  rect.dimensionner(nb_elem, nb_elem_tot - nb_elem, rect_nb_non_zero_tot);

  {
    const int nb_faces_bord = domaine_vdf.nb_faces_bord();
    les_coeff_pression.resize_array(nb_faces_bord);
  }
  auto& carre_tab1 = carre.get_set_tab1();
  auto& rect_tab1 = rect.get_set_tab1();

  // Sparse matrix, Morse storage with Fortran indices:
  // rows numbered 1..n, columns 1..m
  // The k-th nonzero coefficient on row i (1<=i<=n) is (with 1<=k)
  //   M(i,j) = coeff_[tab1_[k]]     in Fortran
  //   M(i,j) = coeff_[tab1_[k-1]-1] in C
  // The column index j of this coefficient (1<=j<=m) is
  //   j = tab2_[tab1_[k]]     in Fortran
  //   j = tab2_[tab1_[k-1]-1] in C
  //
  // Compute the index of the first coefficient on row i
  // in the Morse index array of the two matrices (tab1_)
  {
    int indice = 1; // tab1_ contains a Fortran index (first element at 1)
    for (i = 0; i < nb_elem; i++)
      {
        carre_tab1[i] = indice;
        indice += carre_nb_non_zero[i];
      }
    carre_tab1[i] = indice;

    indice = 1;
    for (i = 0; i < nb_elem; i++)
      {
        rect_tab1[i] = indice;
        indice += rect_nb_non_zero[i];
      }
    rect_tab1[i] = indice;
  }

  // Second step: fill tab2_ = column index of each non-zero term of the matrix
  auto& carre_tab2 = carre.get_set_tab2();
  auto& rect_tab2 = rect.get_set_tab2();

  carre_tab2 = -1;
  rect_tab2 = -1;

  // Diagonal term:
  for (i = 1; i <= nb_elem; i++)
    carre_tab2[carre_tab1[i-1]-1] = i; // Fortran index 1<=i<=nb_elem

  carre_nb_non_zero = 1; // Number of nonzero coefficients on each row
  rect_nb_non_zero = 0;

  // Off-diagonal terms:
  for (int i_face = 0; i_face < nb_faces_internes + nb_faces_periodiques; i_face++)
    {

      // Compute the index of the face to process
      const int face = (i_face < nb_faces_internes)
                       ? premiere_face_interne + i_face
                       : liste_faces_perio[i_face - nb_faces_internes];

      int elem0 = face_voisins(face,0);
      int elem1 = face_voisins(face,1);
      if (elem0 > elem1)
        {
          int tmp = elem1;
          elem1 = elem0;
          elem0 = tmp;
        }
      assert(elem0 >= 0);            // Verify that we have two neighboring elements
      if (elem0 < nb_elem)                              // elem0 is real
        {
          const int ligne = elem0 + 1;                 // Fortran index
          if (elem1 < nb_elem)                            // elem1 is real too
            {
              const int colonne = elem1 + 1;             // Fortran index
              const int n = carre_nb_non_zero[ligne-1]++;
              const auto index = carre_tab1[ligne-1] + n; // Fortran index in tab2
              carre_tab2[index - 1] = colonne;
            }
          else                                           // elem1 is virtual
            {
              const int colonne = elem1 - nb_elem + 1;  // Fortran index
              const int n = rect_nb_non_zero[ligne-1]++;
              const auto index = rect_tab1[ligne-1] + n; // Fortran index in tab2
              rect_tab2[index - 1] = colonne;
            }
        }
    }

  return 1;
}

/*! @brief Computes the coefficients of the pressure matrix with a rho field.
 *
 * @brief If rho_ptr == 0, compute the matrix -div( porosity * grad P ),
 *   otherwise compute -div( porosity/rho grad P ) and *rho_ptr must be a Champ_Fonc_Face_VDF.
 *
 */

int Assembleur_P_VDF::remplir(Matrice& la_matrice, const DoubleVect& volumes_entrelaces,const Champ_Don_base * rho_ptr)
{
  const Domaine_VDF& domaine_vdf   = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_vdf.face_voisins();
  const DoubleVect& face_surfaces = domaine_vdf.face_surfaces();
  //const DoubleVect & volumes_entrelaces = domaine_vdf.volumes_entrelaces();
  const DoubleVect& porosite_face = le_dom_Cl_VDF->equation().milieu().porosite_face();


  const DoubleVect * valeurs_rho = 0;
  if (rho_ptr)
    {
      assert(sub_type(Champ_Fonc_Face_VDF, *rho_ptr));
      valeurs_rho = & (rho_ptr->valeurs());
    }

  // Shortcuts to the square part (real/real element coefficients)
  // and the rectangular part (real/virtual elements) of the matrix
  Matrice_Bloc& matrice = ref_cast(Matrice_Bloc, la_matrice.valeur());
  Matrice_Morse_Sym& carre = ref_cast(Matrice_Morse_Sym, matrice.get_bloc(0,0).valeur());
  Matrice_Morse&      rect  = ref_cast(Matrice_Morse,     matrice.get_bloc(0,1).valeur());

  const int nb_elem = domaine_vdf.nb_elem();
  ArrOfInt carre_nb_non_zero(nb_elem);
  ArrOfInt rect_nb_non_zero(nb_elem);
  carre_nb_non_zero = 1;
  rect_nb_non_zero = 0;

  auto& carre_tab1 = carre.get_set_tab1();
  auto& rect_tab1 = rect.get_set_tab1();
  auto& carre_coeff = carre.get_set_coeff();
  auto& rect_coeff = rect.get_set_coeff();

  carre_coeff = 0.;
  rect_coeff = 0.;

  // Processing internal and periodic faces:
  // For each face between two elements elem0 and elem1, there are four terms to add:
  //   M(elem0,elem0)
  //   M(elem0,elem1)
  //   M(elem1,elem1)
  //   M(elem1,elem0)  (omitted because the matrix is stored as symmetric)

  // Build the list of periodic faces
  ArrOfInt liste_faces_perio;
  const int nb_faces_periodiques = liste_faces_periodiques(liste_faces_perio);
  const int nb_faces_internes = domaine_vdf.nb_faces_internes();
  const int premiere_face_interne = domaine_vdf.premiere_face_int();
  for (int i_face = 0; i_face < nb_faces_internes + nb_faces_periodiques; i_face++)
    {

      // Compute the index of the face to process
      const int num_face = (i_face < nb_faces_internes)
                           ? premiere_face_interne + i_face
                           : liste_faces_perio[i_face - nb_faces_internes];
      // Compute rho on this face
      const double rho_face = (valeurs_rho) ? (*valeurs_rho)[num_face] : 1.;
      // Compute the coefficient
      const double surface  = face_surfaces[num_face];
      const double volume   = volumes_entrelaces[num_face];
      const double porosite = porosite_face[num_face];
      const double coefficient = surface * surface * porosite / (volume * rho_face);
      // Indices of the two neighboring elements (the smaller one in elem0)
      int elem0 = face_voisins(num_face,0);
      int elem1 = face_voisins(num_face,1);
      if (elem0 > elem1)
        {
          int tmp = elem1;
          elem1 = elem0;
          elem0 = tmp;
        }
      if (elem0 < nb_elem)
        {
          // elem0 is real
          const int ligne = elem0 + 1;   // Fortran index
          // Fortran index of the diagonal element (elem0, elem0)
          const auto index_diag = carre_tab1[ligne-1];
          carre_coeff[index_diag - 1] += coefficient;
          if (elem1 < nb_elem)
            {
              // elem1 is real too
              // Fortran index of the diagonal element (elem1, elem1)
              const auto index_diag1 = carre_tab1[elem1]; // at row elem1+1
              // Fortran index of the off-diagonal element (elem0, elem1)
              const int n = carre_nb_non_zero[ligne-1]++;
              const auto index = index_diag + n;
              // Diagonal coefficient
              carre_coeff[index_diag1 - 1] += coefficient;
              // Off-diagonal coefficient
              carre_coeff[index - 1] = - coefficient;
              assert(carre.get_tab2()(index - 1) == elem1 + 1);
            }
          else
            {
              // elem1 is virtual
              const int n = rect_nb_non_zero[ligne-1]++;
              const auto index = rect_tab1[ligne-1] + n; // Fortran index in tab2
              // Off-diagonal coefficient
              rect_coeff[index - 1] = - coefficient;
              assert(rect.get_tab2()(index - 1) == elem1 - nb_elem + 1);
            }
        }
    }

  // Processing the boundary conditions
  const Conds_lim& les_cl = le_dom_Cl_VDF->les_conditions_limites();
  const int nb_cl = les_cl.size();
  for (int num_cl = 0; num_cl < nb_cl; num_cl++)
    {
      const Cond_lim_base& la_cl = les_cl[num_cl].valeur();
      const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());

      // Test on boundary conditions in 2D RZ (symmetry about the axis of revolution is required)
      if (bidim_axi && !sub_type(Symetrie,la_cl))
        {
          const int ndeb = la_front_dis.num_premiere_face();
          const int nfin = ndeb + la_front_dis.nb_faces();
          if (nfin>ndeb && est_egal(face_surfaces[ndeb],0))
            {
              Cerr << "\nFirst face surface is smaller than PrecisionGeom = " << precision_geom << finl;
              Cerr << "May be you have an error in the definition of the boundary conditions." << finl;
              Cerr << "The axis of revolution for this 2D calculation is along Y." << finl;
              Cerr << "So you must specify symmetry boundary condition (symetrie keyword) for the boundary " << la_front_dis.le_nom() << finl;
              exit();
            }
        }

      // For each boundary face between elem0 and a fictitious exterior element
      // with imposed pressure P0, we have:
      //    grad P = (P(elem0) - P0) * surface / volume_entrelace
      // elem0 is an unknown; P0 is added to the right-hand side in "modifier_secmem".
      if (sub_type(Neumann_sortie_libre,la_cl))
        {
          has_P_ref = 1;
          carre.set_est_definie(1);
          const int ndeb = la_front_dis.num_premiere_face();
          const int nfin = ndeb + la_front_dis.nb_faces();
          for (int num_face = ndeb; num_face < nfin; num_face++)
            {
              // Compute rho on this face
              const double rho_face = (valeurs_rho) ? (*valeurs_rho)[num_face] : 1.;
              // Compute the coefficient to add to the matrix
              const double surface  = face_surfaces[num_face];
              // Note: the staggered volume has a special value at the boundary
              // (see Domaine_VDF::calculer_volumes_entrelaces())
              const double volume   = volumes_entrelaces[num_face];
              const double porosite = porosite_face[num_face];
              const double coefficient = Option_VDF::coeff_P_neumann * surface * surface * porosite / (volume * rho_face);
              assert(coefficient > 0.);
              // Index of the neighboring element (one is -1, the other is a real element)
              const int elem0 = face_voisins(num_face, 0);
              const int elem1 = face_voisins(num_face, 1);
              assert(elem0 == -1 || elem1 == -1);
              const int elem = elem0 + elem1 + 1;
              // Add the coefficient to the matrix
              assert(elem < nb_elem);
              const auto index = carre_tab1[elem]; // Fortran index
              carre_coeff[index - 1] += coefficient;
              les_coeff_pression[num_face] = coefficient;
            }
        }
      else
        {
          // For other boundary conditions, no additional term in the matrix
          // (grad P dot n = 0 on the boundary,
          // or time derivative of grad P dot n = 0 on the boundary)
        }
    }
  has_P_ref = (int)mp_max(has_P_ref);

  // Sanity check: no zero element on the diagonal
  for (int i = 0; i < nb_elem; i++)
    {
      const auto index = carre_tab1[i];
      const double coeff_diagonal = carre_coeff[index - 1];
      if (coeff_diagonal == 0.)
        {
          // Cell i has no neighbor: pressure is arbitrary
          carre_coeff[index - 1] = 1.;
        }
    }

  carre.compacte();
  rect.compacte();
  return 1;
}

/*! @brief Modifies the right-hand side to apply boundary conditions.
 *
 * @brief The supported conditions are:
 *   Neumann_sortie_libre,
 *   Entree_fluide_vitesse_imposee,
 *   Dirichlet_paroi_defilante (nothing to do),
 *   Dirichlet_paroi_fixe (nothing to do),
 *   Symetrie (nothing to do)
 *
 */
int Assembleur_P_VDF::modifier_secmem(DoubleTab& secmem)
{
  const Domaine_Cl_VDF& le_dom_cl = le_dom_Cl_VDF.valeur();
  int nb_cond_lim = le_dom_cl.nb_cond_lim();

  for (int indice_cl = 0; indice_cl < nb_cond_lim; indice_cl++)
    {
      const Cond_lim_base& la_cl_base =
        le_dom_cl.les_conditions_limites(indice_cl).valeur();

      const Front_VF& frontiere_vf = ref_cast(Front_VF, la_cl_base.frontiere_dis());

      if (sub_type(Neumann_sortie_libre, la_cl_base))
        {
          modifier_secmem_pression_imposee(ref_cast( Neumann_sortie_libre, la_cl_base),
                                           frontiere_vf,
                                           secmem);
        }
      else if (sub_type(Entree_fluide_vitesse_imposee, la_cl_base))
        {
          modifier_secmem_vitesse_imposee(ref_cast(Entree_fluide_vitesse_imposee ,la_cl_base),
                                          frontiere_vf,
                                          secmem);
        }
      else if (sub_type(Dirichlet_paroi_defilante, la_cl_base))
        {
          // For a sliding wall, nothing to do.
        }
      else if (sub_type(Dirichlet_paroi_fixe, la_cl_base))
        {
          // Nothing to do either.
        }
      else if (sub_type(Symetrie, la_cl_base))
        {
          // Still nothing to do
        }
      else if (sub_type(Periodique, la_cl_base))
        {
          // Nothing to do
        }
      else
        {
          Cerr << "Error in Assembleur_P_VDF::modifier_secmem\n the boundary condition ";
          Cerr << la_cl_base.que_suis_je() << " is not supported." << finl;
          assert(0);
          exit();
        }
    }
  secmem.echange_espace_virtuel();
  return 1;
}

/*! @brief Modifies the right-hand side of the pressure solver for a "Neumann_sortie_libre" condition.
 *
 *  @brief Computation in "pressure increment" mode:
 *   add the pressure increment, i.e. zero (unsteady boundary condition not supported)
 *  Computation in "pressure" mode:
 *   Add the term Pimpose * surface / volume_entrelace to the right-hand side in the discretization of the
 *   pressure at the boundary (between an element elem0 and a fictitious exterior element with imposed pressure):
 *     grad P = (P(elem0) - Pimpose) * surface / volume_entrelace
 *
 */

void Assembleur_P_VDF::modifier_secmem_pression_imposee(const Neumann_sortie_libre& cond_lim,
                                                        const Front_VF& frontiere_vf,
                                                        DoubleTab& secmem)
{
  const Domaine_VDF& le_dom = le_dom_VDF.valeur();
  const IntTab& face_voisins = le_dom.face_voisins();
  if (get_resoudre_increment_pression())
    {
      /*
        const Champ_front_base & champ_front = cond_lim.champ_front();
        if (sub_type(Champ_front_instationnaire_base, champ_front)
        || sub_type(Champ_front_var_instationnaire, champ_front)) {
        Cerr << "Erreur dans Assembleur_P_VDF::modifier_secmem_pression_imposee\n ";
        Cerr << champ_front.que_suis_je();
        Cerr << " + resoudre_increment_pression non code" << finl;
        assert(0);
        exit();
        } else {
        // Champ stationnaire, on ajoute un increment de pression nul.
        // So nothing to do.
        }
      */
    }
  else
    {
      const int nb_faces = frontiere_vf.nb_faces();
      const int num_premiere_face = frontiere_vf.num_premiere_face();
      for (int i = 0; i < nb_faces; i++)
        {
          const int num_face = num_premiere_face + i;
          const double Pimp = cond_lim.flux_impose(i);
          const double coef = les_coeff_pression[num_face] * Pimp;
          const int elem = face_voisins(num_face, 0) + face_voisins(num_face, 1) + 1;
          secmem[elem] += coef;
        }
    }
}

/*! @brief Modifies the right-hand side of the pressure system for an imposed velocity boundary condition.
 *
 *  @brief If solving in pressure increment mode, ...
 *  otherwise nothing to do.
 *
 */
void Assembleur_P_VDF::modifier_secmem_vitesse_imposee(const Entree_fluide_vitesse_imposee& cond_lim,
                                                       const Front_VF& frontiere_vf,
                                                       DoubleTab& secmem)
{
  const Champ_front_base& champ_front = cond_lim.champ_front();
  const Domaine_VDF& le_dom = le_dom_VDF.valeur();
  const DoubleVect& face_surfaces = le_dom.face_surfaces();
  const IntTab& face_voisins = le_dom.face_voisins();

  if (get_resoudre_en_u())
    {
      if (champ_front.instationnaire())
        {
          const DoubleTab& tab_gpoint = champ_front.derivee_en_temps();
          int nb_dim = tab_gpoint.nb_dim();
          bool ch_unif = (tab_gpoint.nb_dim()==1 || tab_gpoint.dimension(0)==1);
          const int nb_faces = frontiere_vf.nb_faces();
          const int num_premiere_face = frontiere_vf.num_premiere_face();
          for (int i = 0; i < nb_faces; i++)
            {
              const int num_face = num_premiere_face + i;
              const double surface = face_surfaces(num_face);
              const int elem0 = face_voisins(num_face, 0);
              const int elem1 = face_voisins(num_face, 1);
              // gpoint is relative to the face normal (pointing towards elem1)
              // Is the normal inward or outward?
              const double signe = (elem0 < 0) ? 1. : -1.;
              // Index of the element adjacent to the boundary face
              const int elem = elem0 + elem1 + 1;
              const int ori = le_dom.orientation(num_face);
              const double gpoint = nb_dim==1 ? tab_gpoint(ori) : tab_gpoint(ch_unif ? 0 : i, ori);

              secmem[elem] += signe * surface * gpoint;
            }
        }
      else
        {
          // The boundary field is steady, nothing to do.
        }
    }
  else
    {
      // Pressure resolution: the boundary condition is imposed elsewhere
    }
}

int Assembleur_P_VDF::modifier_solution(DoubleTab& pression)
{
  // Projection :
  double press_0;
  if(!has_P_ref)
    {
      // Take the minimum pressure as the reference pressure
      // to have the same reference pressure in sequential and parallel runs
      press_0=DMAXFLOAT;
      int nb_elem=le_dom_VDF->domaine().nb_elem();
      for(int n=0; n<nb_elem; n++)
        if (pression[n] < press_0)
          press_0 = pression[n];
      press_0 = mp_min(press_0);
      pression -=press_0;
      pression.echange_espace_virtuel();
    }
  return 1;
}
int Assembleur_P_VDF::assembler_mat(Matrice& matrice,const DoubleVect& volumes_entrelaces,int incr_pression,int resoudre_en_u)
{
  if (!matrice)
    {
      if (je_suis_maitre())
        Cerr << "Assembling the pressure matrix: Assembleur_P_VDF::assembler" << finl;
      // By default, solve in pressure increment
      construire(matrice);
    }
  set_resoudre_increment_pression(incr_pression);
  set_resoudre_en_u(resoudre_en_u);

  remplir(matrice,volumes_entrelaces, 0);
  return 1;
}

/*! @brief Assembles the pressure matrix M such that M*P = div(porosity * grad(P))
 *
 *  @brief and computes the coefficients for modifier_secmem.
 *
 */
int Assembleur_P_VDF::assembler(Matrice& matrice)
{
  if (je_suis_maitre())
    Cerr << "Assembling the pressure matrix: Assembleur_P_VDF::assembler" << finl;
  // By default, solve in pressure increment
  set_resoudre_increment_pression(1);
  set_resoudre_en_u(1);
  construire(matrice);
  const Domaine_VDF& domaine_vdf   = le_dom_VDF.valeur();

  const DoubleVect& volumes_entrelaces = domaine_vdf.volumes_entrelaces();
  remplir(matrice,volumes_entrelaces, 0);
  return 1;
}

/*! @brief Assembles the pressure matrix M such that M*P = div(porosity/rho * grad(P))
 *
 *  @brief and computes the coefficients for modifier_secmem.
 *
 * @param matrice The matrix to assemble. Constraint: either the matrix has not yet been typed (in which case it is "constructed"), or it is the same as from the previous call.
 * @param rho Density field.
 */
int Assembleur_P_VDF::assembler_rho_variable(Matrice& matrice,
                                             const Champ_Don_base& rho)
{
  // assembler_rho_variable was introduced for front-tracking:
  // must explicitly specify whether we solve in pressure increment
  assert(get_resoudre_increment_pression() >= 0);
  // same for solving in u
  assert(get_resoudre_en_u() >= 0);
  // If the matrix has not yet been typed, it must be constructed:
  if (!matrice)
    {
      if (je_suis_maitre())
        {
          Cerr << "Assembling the pressure matrix: ";
          Cerr << "Assembleur_P_VDF::assembler_rho_variable" << finl;
        }
      construire(matrice);
    }
  const Domaine_VDF& domaine_vdf   = le_dom_VDF.valeur();

  const DoubleVect& volumes_entrelaces = domaine_vdf.volumes_entrelaces();
  remplir(matrice,volumes_entrelaces, & rho);
  return 1;
}

/*! @brief Assembles the pressure matrix for a quasi-compressible fluid.
 *
 * @brief The matrix M is such that M*P = div( porosity * grad(P) ).
 *     The resoudre_increment_pression flag is set to zero if not yet assigned.
 *
 * @param tab_rho Density array.
 * @return Always returns 1.
 */
int Assembleur_P_VDF::assembler_QC(const DoubleTab& tab_rho, Matrice& matrice)
{
  // Default for QC: solve in pressure, not in pressure increment.
  if (get_resoudre_increment_pression() < 0)
    {
      set_resoudre_increment_pression(1);
      set_resoudre_en_u(0);
    }
  if (!matrice)
    {
      if (je_suis_maitre())
        {
          Cerr << "Assembling the pressure matrix: ";
          Cerr << "Assembleur_P_VDF::assembler_QC" << finl;
        }
      construire(matrice);
      const Domaine_VDF& domaine_vdf   = le_dom_VDF.valeur();

      const DoubleVect& volumes_entrelaces = domaine_vdf.volumes_entrelaces();
      remplir(matrice,volumes_entrelaces, 0);

      Matrice_Bloc& matrice_bloc=ref_cast(Matrice_Bloc,matrice.valeur());
      Matrice_Morse_Sym& la_matrice =ref_cast(Matrice_Morse_Sym,matrice_bloc.get_bloc(0,0).valeur());
      if (la_matrice.get_est_definie()!=1)
        {
          if ((je_suis_maitre()) && (la_matrice.nb_lignes()==0) && (la_matrice.nb_colonnes()==0))
            {
              Cerr<<"Pressure matrix will not be defined."<<finl;
              exit();
            }

          if ((la_matrice.nb_lignes()>0) && (la_matrice.nb_colonnes()>0))
            {
              Cerr<<"la_matrice(0,0)"<<la_matrice(0,0)<<finl;
              Cerr<<"No imposed pressure --> P(0)=0"<<finl;
              if (je_suis_maitre())    la_matrice(0,0) *= 2;
            }
          la_matrice.set_est_definie(1);
        }
    }
  return 1;
}

/* equation sum_k alpha_k = 1 en Pb_Multiphase */
void Assembleur_P_VDF::dimensionner_continuite(matrices_t matrices, int aux_only) const
{
  if (aux_only) return; //nothing to do
  int e, n, N = ref_cast(Pb_Multiphase, le_dom_Cl_VDF->equation().probleme()).nb_phases(), ne_tot = le_dom_VDF->nb_elem_tot();
  Stencil stencil(0, 2);

  for (e = 0; e < le_dom_VDF->nb_elem(); e++)
    for (n = 0; n < N; n++) stencil.append_line(e, N * e + n);
  Matrix_tools::allocate_morse_matrix(ne_tot, N * ne_tot, stencil, *matrices.at("alpha"));
}

void Assembleur_P_VDF::assembler_continuite(matrices_t matrices, DoubleTab& secmem, int aux_only) const
{
  if (aux_only) return;
  const DoubleTab& alpha = ref_cast(Pb_Multiphase, le_dom_Cl_VDF->equation().probleme()).equation_masse().inconnue().valeurs();
  Matrice_Morse& mat = *matrices.at("alpha");
  const DoubleVect& ve = le_dom_VDF->volumes(), &pe = le_dom_Cl_VDF->equation().milieu().porosite_elem();
  int e, n, N = alpha.line_size();
  /* right-hand side: multiply by porosity * volume so that the pressure system is symmetric in Cartesian coordinates */
  for (e = 0; e < le_dom_VDF->nb_elem(); e++)
    for (secmem(e) = -pe(e) * ve(e), n = 0; n < N; n++) secmem(e) += pe(e) * ve(e) * alpha(e, n);
  /* matrice */
  for (e = 0; e < le_dom_VDF->nb_elem(); e++)
    for (n = 0; n < N; n++) mat(e, N * e + n) = -pe(e) * ve(e);
}

/* norme pour assembler_continuite */
DoubleTab Assembleur_P_VDF::norme_continuite() const
{
  const DoubleVect& pe = le_dom_Cl_VDF->equation().milieu().porosite_elem(), &ve = le_dom_VDF->volumes();
  DoubleTab norm(le_dom_VDF->nb_elem());
  for (int e = 0; e < le_dom_VDF->nb_elem(); e++) norm(e) = pe(e) * ve(e);
  return norm;
}

const Domaine_dis_base& Assembleur_P_VDF::domaine_dis_base() const
{
  return le_dom_VDF.valeur();
}

const Domaine_Cl_dis_base& Assembleur_P_VDF::domaine_Cl_dis_base() const
{
  return le_dom_Cl_VDF.valeur();
}

void Assembleur_P_VDF::associer_domaine_dis_base(const Domaine_dis_base& le_dom_dis)
{
  le_dom_VDF = ref_cast(Domaine_VDF, le_dom_dis);
}

void Assembleur_P_VDF::associer_domaine_cl_dis_base(const Domaine_Cl_dis_base& le_dom_Cl_dis)
{
  le_dom_Cl_VDF = ref_cast(Domaine_Cl_VDF, le_dom_Cl_dis);
}

void Assembleur_P_VDF::completer(const Equation_base& Eqn)
{
  // CCa 30/04/99: not sure if anything needs to be done here
  ;
}
