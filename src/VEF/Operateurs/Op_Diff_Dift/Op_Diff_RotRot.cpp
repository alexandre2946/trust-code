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

#include <Champ_P1_isoP1Bulle.h>
#include <LecFicDistribueBin.h>
#include <Matrice_Morse_Sym.h>
#include <Champ_Uniforme.h>
#include <Op_Diff_RotRot.h>

#include <Domaine_Cl_VEF.h>

#include <SFichier.h>
#include <Solv_GCP.h>
#include <SSOR.h>

Implemente_instanciable(Op_Diff_RotRot, "Op_Diff_VEF_ROTROT_P1NC", Operateur_Diff_base);

Sortie& Op_Diff_RotRot::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Op_Diff_RotRot::readOn(Entree& s) { return s; }

const Domaine_VEF& Op_Diff_RotRot::domaine_vef() const { return ref_cast(Domaine_VEF, le_dom_vef.valeur()); }

void Op_Diff_RotRot::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base& inco)
{
  const Domaine_VEF& zvef = ref_cast(Domaine_VEF, domaine_dis);
  const Domaine_Cl_VEF& zclvef = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis);
  le_dom_vef = zvef;
  la_zcl_vef = zclvef;

  curl_.associer(domaine_dis, domaine_Cl_dis, inco);
  rot_.associer(domaine_dis, domaine_Cl_dis, inco);

  //////////////////////////////////////////////
  // Define the vorticity field
  vorticite_.typer("Champ_P1_isoP1Bulle");
  Champ_P1_isoP1Bulle& vorticite = ref_cast(Champ_P1_isoP1Bulle, vorticite_.valeur());

  vorticite.associer_domaine_dis_base(zvef);
  vorticite.nommer("vorticite");

  if (dimension == 2)
    vorticite.fixer_nb_comp(1);
  else
    vorticite.fixer_nb_comp(dimension);

  int nb_tot = zvef.nb_elem() + zvef.nb_som() - 1;
  vorticite.fixer_nb_valeurs_nodales(nb_tot);

  vorticite.fixer_unite("s-1");

  //////////////////////////////////////////////////

  solveur_.typer("Solv_GCP");
  OWN_PTR(Precond_base) p;
  p.typer("SSOR");
  ref_cast(Solv_GCP,solveur.valeur()).set_precond(p);
  assembler_matrice(matrice_vorticite_);
  tester();

}

DoubleTab& Op_Diff_RotRot::calculer(const DoubleTab& vitesse, DoubleTab& diffusion) const
{
  diffusion = 0;
  return ajouter(vitesse, diffusion);
}

//Method of the diffusion operator in the form
//curl of the vorticity without const — the only one that
//should be called.
DoubleTab& Op_Diff_RotRot::ajouter(const DoubleTab& vitesse, DoubleTab& diffusion) const
{
  Cerr << "entering OpDiffRotRot" << finl;
  DoubleTab curl(matrice_vorticite_->ordre());

  curl_.calculer(vitesse, curl);
  //curl=-1*curl;
  calculer_vorticite(vorticite_->valeurs(), curl);
  rot_.calculer(vorticite_->valeurs(), diffusion);

  Cerr << "je sors de OpDiffRotRot" << finl;

  return diffusion;

}

int Op_Diff_RotRot::calculer_vorticite(DoubleTab& solution, const DoubleTab& curl) const
{
  const Domaine& domaine = domaine_vef().domaine();
  //static int nb_appel2=0;

  // Vorticity resolution: only sequential case is considered for now.
  // For the parallel case, the approach used in class N_S.cpp should be followed.

  if (Process::is_sequential())
    {
      // The vorticity matrix is treated as a Morse matrix
      // by default, and the inverse() method is applied
      // to solve the linear system.
      const Matrice_Morse_Sym& la_matrice = ref_cast(Matrice_Morse_Sym, matrice_vorticite_.valeur());
      //      Matrice_Morse& la_matrice = (Matrice_Morse&) matrice_vorticite_.valeur();
      DoubleTab solution_temporaire(la_matrice.ordre());

      assert(solution_temporaire.size() == solution.size() - 1);
      assert(curl.size() == solution_temporaire.size());

      // Solve the linear system to compute the vorticity.
      // The vorticity is then stored in the variable solution.
      //      la_matrice.inverse(curl,solution_temporaire,1e-15);
      Solv_GCP& solv = ref_cast_non_const(Solv_GCP, solveur_.valeur());
      solv.set_seuil(1e-17);
      solv.resoudre_systeme(la_matrice, curl, solution_temporaire);

      // Copy values from solution_temporaire into solution.
      // One slot is not filled yet: it is done later.
      for (int i = 0; i < solution_temporaire.size(); i++)
        solution[i] = solution_temporaire[i];

      // Fill in the last entry of "solution"
      int sommet = domaine.nb_som() - 1;

      for (int i = 0; i < curl_.elem_som_size(sommet); i++)
        solution[sommet] += solution_temporaire[curl_.elements_pour_sommet(sommet, i)];

      //nb_appel2++;

    }

  return 1;
}

//////////////////////////////////////////////////////
/* Functions for assembling the vorticity matrix */
/////////////////////////////////////////////////////
/* Compute the vorticity matrix for the triple problem:
 / we consider the case without boundary conditions for the vorticity
 / for the moment, and the parallel part of the algorithm is not coded.
 / REM: currently only works in 2D
 */
int Op_Diff_RotRot::assembler_matrice(Matrice& matrice)
{
  const Domaine& domaine = domaine_vef().domaine();

  int colonne_a_remplir_tab2, colonne_a_remplir_coeff;

  Cerr << "Assembling the vorticity matrix..." << finl;
  matrice.typer("Matrice_Morse_Sym");
  //Matrice_Morse& la_matrice=(Matrice_Morse&) matrice.valeur();

  // Size the matrix appropriately.
  // The matrix is square of size (nb_elem + nb_som).
  // The number of non-zero coefficients does not exceed
  // (2*(dimension+1) + 1)*nb_elem + (nb_som-1)*(nb_som-1).
  // NOTE: a BASIS of our space consists of element indicator functions
  // + hat functions - 1 of those functions.
  // Otherwise the sum of all hat functions - the sum
  // of all indicator functions = 0.
  // By default we remove the last hat function to
  // to form our basis.
  // Cf. Paper in Latex/Vorticity
  int nombre_coeff_non_nuls = (2 * dimension + 3) * domaine.nb_elem() + (domaine.nb_som() - 1) * (domaine.nb_som() - 1);
  //la_matrice.dimensionner(domaine.nb_elem()+domaine.nb_som()-1,nombre_coeff_non_nuls);
  Matrice_Morse la_matrice(domaine.nb_elem() + domaine.nb_som() - 1, nombre_coeff_non_nuls);

  // In case the matrix already exists in a file
  Nom nomfic("Vorticite.sv");
  LecFicDistribueBin vorticite;
  int fic_vorticite_existe;

  if (vorticite.ouvrir(nomfic))
    fic_vorticite_existe = 1;
  else
    fic_vorticite_existe = 0;

  if (fic_vorticite_existe)
    {
      Cerr << "Reading the vorticity matrix from file: " << nomfic << finl;
      vorticite >> la_matrice;
      vorticite.close();
      Cerr << "Done reading the vorticity matrix." << finl;
      return 1;
    }

  Cerr << "Assembling the vorticity matrix " << nomfic << finl;

  // Now fill the matrix row by row.
  // Elements and vertices are numbered starting from 0.
  // However, column and row indices of the arrays
  // in a FORTRAN Morse matrix start at 1.
  // Construction example: matrice[C++ index] = FORTRAN index

  // Key parameters for filling the arrays
  nombre_coeff_non_nuls = 1;
  colonne_a_remplir_tab2 = 0; // for C++
  colonne_a_remplir_coeff = 0; // for C++

  // Start by filling the rows of the sub-matrix
  // of size nb_elem * (nb_elem + nb_som): cf. matrix structure
  for (int numero_elem = 0; numero_elem < domaine.nb_elem(); numero_elem++)
    {
      // For a given element "numero_elem", compute the list
      // of vertices belonging to that element.
      IntList sommets_pour_elem = sommets_pour_element(numero_elem);
      Tri(sommets_pour_elem); //sorted array

      // Fill tab1
      la_matrice.get_set_tab1()(numero_elem) = nombre_coeff_non_nuls;

      // Fill tab2 and coeff

      // The nb_elem first entries of row "numero_elem" (FORTRAN indexing)
      la_matrice.get_set_tab2()(colonne_a_remplir_tab2) = numero_elem + 1;
      la_matrice.get_set_coeff()(colonne_a_remplir_coeff) = remplir_elem_elem_EF(numero_elem);

      // Increment nb_coeff_non_nuls because one non-zero entry was just filled
      nombre_coeff_non_nuls++;

      // Increment column counters since array slots were just filled
      colonne_a_remplir_tab2++;
      colonne_a_remplir_coeff++;

      // The nb_som entries of row "numero_elem"
      for (int i = 0; i < domaine.nb_som_elem(); i++)
        {
          // If "sommet_pour_elem[i]" differs from the index of the last vertex,
          // store the correct coefficient
          if (sommets_pour_elem[i] != domaine.nb_som() - 1)
            {
              la_matrice.get_set_tab2()(colonne_a_remplir_tab2) = domaine.nb_elem() + sommets_pour_elem[i] + 1; //for FORTRAN indexing
              la_matrice.get_set_coeff()(colonne_a_remplir_coeff) = remplir_elem_som_EF(numero_elem, sommets_pour_elem[i]);

              // Increment nb_coeff_non_nuls since one non-zero entry was just filled
              nombre_coeff_non_nuls++;

              // Increment colonne_* to avoid overwriting already stored coefficients
              colonne_a_remplir_tab2++;
              colonne_a_remplir_coeff++;
            }

        }

      // No need to modify the integers nombre_coeff_non_nuls and
      // colonne_*: they are already at the correct value.

    }

  // Fill the rows of the sub-matrix of size
  // nb_som * (nb_elem + nb_som): cf. matrix structure
  for (int numero_som = 0; numero_som < domaine.nb_som() - 1; numero_som++)
    {
      // Store the arrays we need along with their respective sizes
      IntList Elem_pour_sommet = elements_pour_sommet(numero_som);
      IntList Sommets_voisins = sommets_voisins(numero_som, Elem_pour_sommet);
      Tri(Elem_pour_sommet); //liste triee
      Tri(Sommets_voisins); //liste triee

      // Fill tab1 (FORTRAN indexing)
      la_matrice.get_set_tab1()(domaine.nb_elem() + numero_som) = nombre_coeff_non_nuls;

      // Fill tab2 and coeff

      // The nb_elem first entries of row "numero_som"
      for (int i = 0; i < Elem_pour_sommet.size(); i++)
        {
          // FORTRAN indexing
          la_matrice.get_set_tab2()(colonne_a_remplir_tab2) = Elem_pour_sommet[i] + 1;
          la_matrice.get_set_coeff()(colonne_a_remplir_coeff) = remplir_som_elem_EF(Elem_pour_sommet[i], numero_som);

          // Increment nb_coeff_non_nuls since one non-zero entry was just filled
          nombre_coeff_non_nuls++;

          // Increment column indices accordingly
          colonne_a_remplir_tab2++;
          colonne_a_remplir_coeff++;
        }

      // No need to modify nb_coeff_non_nuls and colonne_*: they are at the correct value.

      // The nb_som entries of row "numero_som"
      for (int i = 0; i < Sommets_voisins.size(); i++)
        {
          // Retrieve the global index of the neighbouring vertex
          int numero_som_global = Sommets_voisins[i];

          // If "numero_som_global" differs from the index of the last vertex,
          // store the correct coefficient

          if (numero_som_global != domaine.nb_som() - 1)
            {
              // FORTRAN indexing
              la_matrice.get_set_tab2()(colonne_a_remplir_tab2) = domaine.nb_elem() + Sommets_voisins[i] + 1;
              la_matrice.get_set_coeff()(colonne_a_remplir_coeff) = remplir_som_som_EF(numero_som, Sommets_voisins[i], Elem_pour_sommet);

              // Increment nb_coeff_non_nuls since one non-zero entry was just filled
              nombre_coeff_non_nuls++;

              // Increment indices
              colonne_a_remplir_tab2++;
              colonne_a_remplir_coeff++;
            }

        }

      // No need to increment nombre_coeff_non_nuls and colonne_*:
      // they are already at the correct value.
    }

  // By convention for Morse matrices, the last entry
  // of tab1 is tab1[nb_elem + nb_som + 1] and equals
  // the total number of non-zero coefficients + 1,
  // i.e. with our algorithm: nombre_coeff_non_nuls
  la_matrice.get_set_tab1()(domaine.nb_elem() + domaine.nb_som() - 1) = nombre_coeff_non_nuls;

  //   if(Debog::mode_db==2) Debog::save_matrix_seq(la_matrice);
  //   else if(Debog::mode_db==3)
  //     Debog::save_and_distribute_matrix_seq(la_matrice);

  Matrice_Morse_Sym& la_matrice_sym = ref_cast(Matrice_Morse_Sym, matrice.valeur());
  la_matrice_sym = la_matrice;

  Cerr << "Fin de l'assemblage de la matrice de vorticite. " << finl;

  return 1;
}

/* For a given element "numero_elem", returns */
/* the list of vertices belonging to that element */
IntList Op_Diff_RotRot::sommets_pour_element(int numero_elem) const
{
  IntList resultat;
  int numero_global_sommet = 0;
  const Domaine& domaine = domaine_vef().domaine();

  // As a precaution but normally not needed
  if (!resultat.est_vide())
    resultat.vide();

  for (int i = 0; i < domaine.nb_som_elem(); i++)
    {
      numero_global_sommet = domaine.sommet_elem(numero_elem, i);
      resultat.add_if_not(numero_global_sommet);
    }

  return resultat;
}

/* For a given vertex "numero_sommet", returns */
/* the list of elements containing that vertex */
IntList Op_Diff_RotRot::elements_pour_sommet(int numero_sommet) const
{
  IntList resultat;
  int numero_global_som;
  const Domaine& domaine = domaine_vef().domaine();

  // As a precaution but normally not needed
  if (!resultat.est_vide())
    resultat.vide();

  // Not very efficient but works regardless of the dimension:
  // loop over elements, check for each element
  // its vertices, then compare those vertices to the input parameter;
  // if one of them coincides with the input parameter, store the element.
  for (int numero_elem = 0; numero_elem < domaine.nb_elem(); numero_elem++)
    for (int numero_som = 0; numero_som < domaine.nb_som_elem(); numero_som++)
      {
        numero_global_som = domaine.sommet_elem(numero_elem, numero_som);
        if (numero_sommet == numero_global_som)
          resultat.add_if_not(numero_elem);
      }

  return resultat;

}

/* For a given vertex "numero_sommet", returns */
/* the list of neighbouring vertices of "numero_sommet" */
/* Parameter: the list of elements containing "numero_sommet" */
/* Search within these elements to obtain the result */
/* NOTE: the result list contains the vertex "numero_sommet" itself */
IntList Op_Diff_RotRot::sommets_voisins(int numero_sommet, const IntList& liste) const
{
  IntList resultat;
  int numero_global_som;
  const Domaine& domaine = domaine_vef().domaine();

  //Retrieve the vertices of the elements in "liste"
  //then compare them to "numero_som" and keep them
  //without duplicates.
  for (int numero_elem_loc = 0; numero_elem_loc < liste.size(); numero_elem_loc++)
    for (int numero_som = 0; numero_som < domaine.nb_som_elem(); numero_som++)
      {
        numero_global_som = domaine.sommet_elem(liste[numero_elem_loc], numero_som);
        resultat.add_if_not(numero_global_som); //contains "numero_sommet"
      }

  //   Cerr << "Affichage de sommets_voisins pour numero_som " << numero_sommet
  //        << finl;
  //   for (int i=0;i<resultat.size();i++)
  //     Cerr << resultat[i] << finl;

  return resultat;
}

/* Sorting function for an IntList */
/* The sort is performed in ascending order */
/* REM: there are no duplicates in the sorted list */
void Op_Diff_RotRot::Tri(IntList& liste_a_trier) const
{
  if (liste_a_trier.est_vide())
    {
      Cerr << "Error in Op_Diff_RotRot::Tri()." << finl;
      Cerr << "The list to sort is empty: exiting." << finl;
      Process::exit();
    }

  IntList temporaire;
  int minimum;

  while (!liste_a_trier.est_vide())
    {
      minimum = liste_a_trier[0];

      for (int i = 0; i < liste_a_trier.size(); i++)
        minimum = (minimum <= liste_a_trier[i] ? minimum : liste_a_trier[i]);

      temporaire.add_if_not(minimum);
      liste_a_trier.suppr(minimum);
    }

  for (int i = 0; i < temporaire.size(); i++)
    liste_a_trier.add_if_not(temporaire[i]);

}

/* For element "numero_elem", returns the coefficient */
/* to place in the nb_elem * nb_elem sub-matrix */
/* at row "numero_elem", column "numero_elem" */
/* EF matrix */
double Op_Diff_RotRot::remplir_elem_elem_EF(const int numero_elem) const
{
  return 1. * domaine_vef().volumes(numero_elem);
}

/* For element "numero_elem", returns the coefficient */
/* to place in the nb_elem * nb_som sub-matrix */
/* at row "numero_elem", column "numero_som" */
/* EF matrix */
double Op_Diff_RotRot::remplir_elem_som_EF(const int numero_elem, const int numero_som) const
{
  return (1. * domaine_vef().volumes(numero_elem) / (dimension + 1));
}

/* For vertex "numero_som", returns the coefficient */
/* to place in the nb_som * nb_elem sub-matrix */
/* at row "numero_som", column "numero_elem" */
/* EF matrix */
double Op_Diff_RotRot::remplir_som_elem_EF(const int numero_elem, const int numero_som) const
{
  return (1. * domaine_vef().volumes(numero_elem) / (dimension + 1));
}

/* For vertex "numero_som", returns the coefficient */
/* to place in the nb_som * nb_som sub-matrix */
/* at row "numero_som", column "sommet_voisin" */
/* "elem_voisins" is the array of elements containing "numero_som" */
/* EF matrix */
double Op_Diff_RotRot::remplir_som_som_EF(const int numero_som, const int sommet_voisin, const IntList& elem_voisins) const
{
  double resultat = 0.;
  int test = 0;

  //First test: if numero_som == sommet_voisin
  //then we can immediately return the result
  if (numero_som == sommet_voisin)
    {
      for (int i = 0; i < elem_voisins.size(); i++)
        {
          resultat += domaine_vef().volumes(elem_voisins[i]);
        }

      //Account for the spatial dimension
      resultat *= 2. / ((dimension + 1) * (dimension + 2));

      return resultat;
    }

  //Sinon:
  // Compute the contribution to the result

  for (int i = 0; i < elem_voisins.size(); i++)
    {
      if (sommets_pour_element(elem_voisins[i]).contient(sommet_voisin))
        {
          resultat += domaine_vef().volumes(elem_voisins[i]);
          test++;
        }
    }

  //Account for the spatial dimension
  resultat *= 1. / ((dimension + 1) * (dimension + 2));

  //Verify that sommet_voisin was indeed in the neighbourhood of numero_voisin
  if (test == 0)
    {
      Cerr << "Error in Op_Diff_RotRot::remplir_som_som_EF." << finl;
      Cerr << "sommet_voisin is not in the neighbourhood of numero_voisin." << finl;
      Cerr << "Exiting." << finl;
      Process::exit();
    }

  return resultat;

}

DoubleTab Op_Diff_RotRot::vecteur_normal(const int face, const int elem) const
{
  assert(dimension == 2);

  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  DoubleTab le_vecteur_normal(dimension);

  for (int composante = 0; composante < dimension; composante++)

    le_vecteur_normal(composante) = domaine_VEF.face_normales(face, composante) * domaine_VEF.oriente_normale(face, elem);

  return le_vecteur_normal;
}

int Op_Diff_RotRot::tester() const
{
  //  const Domaine& domaine = domaine_Vef().domaine();

  //   if (vorticite_->nb_valeurs_nodales() !=
  //       domaine.nb_elem()+domaine.nb_som()-1 )
  //     {
  //       Cerr << "Probleme dans la definition de la vorticite." << finl;
  //       Cerr << "The number of registered components is incorrect." << finl;
  //       Process::exit();
  //     }

  //Test the vorticity matrix
  const Matrice_Morse_Sym& la_matrice = ref_cast(Matrice_Morse_Sym, matrice_vorticite_.valeur());
  Solv_GCP& solv = ref_cast_non_const(Solv_GCP, solveur_.valeur());

  DoubleTab resultat(la_matrice.ordre());
  //DoubleTab resultat1(la_matrice.ordre());
  DoubleTab secmem(la_matrice.ordre());
  DoubleTab secmem1(la_matrice.ordre());
  secmem = 0.;
  secmem[0] = -0.25;
  secmem[4] = -0.0833333;
  secmem[5] = -0.0833333;
  secmem1 = 1e-10 + 1e-15;

  SFichier fic("Matrice.test");
  la_matrice.imprimer_formatte(fic);

  resultat = 0.;
  solv.set_seuil(1e-17);
  solv.resoudre_systeme(la_matrice, secmem, resultat);
  Cerr << "Resultat de l'inversion" << finl;
  for (int i = 0; i < 8; i++)
    Cerr << resultat[i] << " ; ";
  Cerr << finl;

  solv.resoudre_systeme(la_matrice, secmem1, resultat);
  Cerr << "Resultat de l'inversion" << finl;
  for (int i = 0; i < 8; i++)
    Cerr << resultat[i] << " ; ";
  Cerr << finl;

  return 1;
}

void Op_Diff_RotRot::associer_diffusivite(const Champ_base& diffu)
{
  diffusivite_ = ref_cast(Champ_Uniforme, diffu);
}

const Champ_base& Op_Diff_RotRot::diffusivite() const
{
  return diffusivite_;
}
