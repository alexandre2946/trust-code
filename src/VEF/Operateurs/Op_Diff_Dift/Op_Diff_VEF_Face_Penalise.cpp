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
#include <Op_Diff_VEF_Face_Penalise.h>
#include <Dirichlet_paroi_fixe.h>


Implemente_instanciable(Op_Diff_VEF_Face_Penalise,"Op_Diff_VEFpenalise_P1NC",Op_Diff_VEF_Face);

/* Mandatory implementation of the printOn() function */
Sortie& Op_Diff_VEF_Face_Penalise::
printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

/* Mandatory implementation of the readOn() function */
Entree& Op_Diff_VEF_Face_Penalise::readOn(Entree& is )
{
  /* IMPORTANT NOTE: this class was created by derivation but
   * it should only be applied to non-turbulent Navier-Stokes equations,
   * hence the following checks.
   */

  Cerr<<"In Op_Diff_VEF_Face_Penalise::readOn()"<<finl;

  /* Check problem dimension since the theory has only been validated
   * in 2D.
   */
  if (dimension != 2)
    {
      Cerr << "Error in Op_Diff_VEF_Face_Penalise::readOn()" << finl;
      Cerr << "Problem dimension must be 2" << finl;
      Process::exit();
    }

  return is ;
}

/*
 *
 *
 *
 *
 *
 */

/*! @brief Method that computes the contribution of the operator.
 *
 */
DoubleTab& Op_Diff_VEF_Face_Penalise::
calculer(const DoubleTab& inconnue, DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(inconnue,resu);
}

/*
 *
 *
 *
 *
 *
 */

/*! @brief Method that computes the velocity at time n+1 when the explicit scheme is used.
 *
 */
DoubleTab& Op_Diff_VEF_Face_Penalise::
ajouter(const DoubleTab& inconnue, DoubleTab& resu) const
{
  //  Cerr << "Entering ajouter() of penalisation" << finl;
  int nb_composante,numero_global_face,local;
  int face_penalisation,face;
  double coeff;
  IntList voisinage_ordre1,voisinage_ordre1_strict;
  IntList voisinage_ordre2,voisinage_ordre2_strict;
  IntList ensemble_faces;

  /* Compute the result due to the penalisation matrix
   * i.e. sum_j U_j P_i,j where P_i,j is the term of the
   * penalisation matrix.
   */
  for (nb_composante = 0; nb_composante < dimension; nb_composante++)
    {
      for (numero_global_face = 0; numero_global_face < domaine_vef().nb_faces();
           numero_global_face++)
        {
          /* Compute the neighbourhoods of numero_global_face */
          voisinage(numero_global_face,voisinage_ordre1);
          voisinage(voisinage_ordre1,voisinage_ordre2);
          reduction(voisinage_ordre1,voisinage_ordre2,voisinage_ordre2_strict);
          //           Cerr << "Nombre iterations " << numero_global_face << finl;
          //           Cerr << "Nombre de faces " << domaine_vef().nb_faces() << finl;

          /* Compute resu(numero_face_global,nb_comp)
           * for voisinage_ordre2_strict
           */
          for (local = 0 ; local < voisinage_ordre2_strict.size(); local++)
            {
              face = voisinage_ordre2_strict[local];
              faces_communes(numero_global_face,face,ensemble_faces);
              //               Cerr << "Taille voisinage_ordre_1 " << voisinage_ordre1.size()
              //                    << finl;
              //               Cerr << "Taille voisinage_ordre2 " << voisinage_ordre2.size()
              //                    << finl;
              //               Cerr << "Taille voisinage_ordre2_strict " << voisinage_ordre2_strict.size() << finl;
              //               Cerr << "Taille ensemble_faces " << ensemble_faces.size() << finl;

              // if (ensemble_faces.size() != 1)
              //                 {
              //                   Cerr << "Erreur dans ajouter() de la penalisation." << finl;
              //                   Process::exit();
              //                 }

              for (int mm = 0; mm < ensemble_faces.size(); mm++)
                {
                  face_penalisation = ensemble_faces[mm];

                  resu(numero_global_face,nb_composante) +=
                    inconnue(face,nb_composante)*
                    signe(numero_global_face,face)*
                    coefficient_penalisation(face_penalisation)*
                    longueur(face_penalisation)*
                    1./3.;
                }

            }

          /* Compute resu(numero_global_face,nb_comp)
           * for voisinage_ordre1_strict
           */
          for (local = 0 ; local < voisinage_ordre1_strict.size(); local++)
            {
              assert(voisinage_ordre1_strict.size() == 2
                     || voisinage_ordre1_strict.size() == 4);

              face = voisinage_ordre1_strict[local];
              face_penalisation = autre_face(numero_global_face,face);

              resu(numero_global_face,nb_composante) +=
                inconnue(face,nb_composante)*
                coefficient_penalisation(face_penalisation)*
                longueur(face_penalisation)*
                (-1./3.);
            }

          /* Case 1: the face under consideration is on the domain boundary
           * and numero_global_face = face.
           */
          if (numero_global_face < domaine_vef().nb_faces_bord())
            for (local = 0; local < voisinage_ordre1.size(); local++)
              {
                assert(voisinage_ordre1.size() == 3);

                face_penalisation = voisinage_ordre1[local];

                if (numero_global_face == face_penalisation)
                  coeff = 1.;
                else
                  coeff = 1./3.;

                resu(numero_global_face,nb_composante) +=
                  inconnue(numero_global_face,nb_composante) *
                  coeff *
                  coefficient_penalisation(face_penalisation)*
                  longueur(face_penalisation);
              }

          /* Case 2: the face is internal and numero_global_face = face */
          else
            for (local = 0; local < voisinage_ordre1_strict.size(); local++)
              {
                assert( voisinage_ordre1_strict.size() == 4);

                face_penalisation = voisinage_ordre1_strict[local];

                resu(numero_global_face,nb_composante) +=
                  inconnue(numero_global_face,nb_composante) *
                  coefficient_penalisation(face_penalisation)*
                  longueur(face_penalisation)*
                  1./3.;
              }
        }
    }


  /* Add the classical diffusion matrix */
  Op_Diff_VEF_Face::ajouter(inconnue,resu);

  //  Cerr << "Exiting calcul_matrice_de_penalisation_" << finl;
  return resu;
}

/*
 *
 *
 *
 *
 *
 */

/*! @brief Method computing the neighbourhood of a face.
 *
 */
void Op_Diff_VEF_Face_Penalise::
voisinage(const int Numero_face, IntList& Voisinage) const
{
  //  Cerr << "Entering voisinage for a face index" << finl;

  /* Clear the list Voisinage to avoid surprises */
  if (! Voisinage.est_vide() ) Voisinage.vide();

  /* Declaration of the main local parameters. */
  int numero_local;

  /* The number of faces per element in the discretisation domain.
   * NOTE: prisms are excluded.
   */
  const int nb_faces_element = domaine().nb_faces_elem();

  /* The neighbouring elements of Numero_face */
  const int voisin1 = domaine_vef().face_voisins(Numero_face,1);
  const int voisin2 = domaine_vef().face_voisins(Numero_face,0);

  /* Retrieve the faces of voisin* and inject them into the list Voisinage.
   * The element must exist, hence the first test.
   */
  if (voisin1 != -1)
    {
      for (numero_local = 0; numero_local < nb_faces_element; numero_local++)
        {
          /* Retrieve the global index of each face of voisin* */
          const int numero_global_face =
            domaine_vef().elem_faces(voisin1,numero_local);

          /* Then store this index in the list Voisinage. */
          Voisinage.add_if_not(numero_global_face);

        }// end for

    }// end if

  if (voisin2 != -1)
    {
      for (numero_local = 0; numero_local < nb_faces_element; numero_local++)
        {
          /* Retrieve the global index of each
           * face of numero_element_*
           */
          const int numero_global_face =
            domaine_vef().elem_faces(voisin2,numero_local);

          /* Then store this index in the list Voisinage. */
          Voisinage.add_if_not(numero_global_face);

        }// end for

    }//end if

  //  Cerr << "Exiting voisinage for a face index" << finl;
}

/*
 *
 *
 *
 *
 *
 */

/*! @brief Method computing the neighbourhood of a list of faces.
 *
 */
void  Op_Diff_VEF_Face_Penalise::
voisinage(const IntList& Ensemble_faces, IntList& Voisinage) const
{
  //  Cerr << "Entering voisinage for a set of faces" << finl;

  /* Clear the list Voisinage to avoid surprises */
  if (! Voisinage.est_vide() ) Voisinage.vide();

  /* Declaration of the main local parameters. */
  int nb_elements_Ensemble_faces,nb_elements_liste_temporaire;
  IntList liste_temporaire;

  for (nb_elements_Ensemble_faces = 0;
       nb_elements_Ensemble_faces < Ensemble_faces.size();
       nb_elements_Ensemble_faces++)
    {
      /* Loop-internal parameter */
      const int numero_face_dans_Ensemble_faces =
        Ensemble_faces[nb_elements_Ensemble_faces];

      /* Clear the temporary list each time we change face */
      if (! liste_temporaire.est_vide() ) liste_temporaire.vide();

      /* Place the neighbourhood of numero_face_dans_Ensemble_faces
       * into the temporary list.
       */
      voisinage(numero_face_dans_Ensemble_faces,liste_temporaire);
      //      Cerr << "Taille liste temporaire " << liste_temporaire.size() << finl;


      for (nb_elements_liste_temporaire = 0;
           nb_elements_liste_temporaire < liste_temporaire.size();
           nb_elements_liste_temporaire++)
        {
          //          Cerr << "les elements " << liste_temporaire[nb_elements_liste_temporaire] << finl;
          /* Second-loop-internal parameter */
          const int numero_face_dans_liste_temporaire =
            liste_temporaire[nb_elements_liste_temporaire];

          /* Finally, store in Voisinage */
          Voisinage.add_if_not(numero_face_dans_liste_temporaire);

        }// end second for

    }// end first for

  // Cerr << "Exiting voisinage for a set of faces" << finl;
}

/*
 *
 *
 *
 *
 *
 */


double  Op_Diff_VEF_Face_Penalise::
signe(const int Face1, const int Face2) const
{
  //  Cerr << "Entering signe" << finl;

  /* Local parameters of the procedure */
  int numero_local;

  /* Retrieve the number of vertices per face in the discretised domain.
   * NOTE: prismes are excluded by convention.
   */
  const int nb_sommets_par_face = domaine_vef().nb_som_face();

  /* Create a list containing the vertices of Face2 */
  IntList sommets_Face2;

  for (numero_local = 0 ; numero_local < nb_sommets_par_face ; numero_local++)
    sommets_Face2.add(domaine_vef().face_sommets(Face2,numero_local));

  /* Check whether any vertices of Face2 belong to Face1.
   * If yes, return 1; otherwise return -1.
   */
  for (numero_local = 0; numero_local < nb_sommets_par_face ; numero_local++)
    if ( sommets_Face2.contient(domaine_vef().face_sommets(Face1,numero_local)))
      return 1.;

  return -1.;

  //  Cerr << "Exiting signe" << finl;
}

/*
 *
 *
 *
 *
 *
 */

/*! @brief Member function returning the penalisation coefficient associated with each face of the primary mesh.
 *
 */
double  Op_Diff_VEF_Face_Penalise::
coefficient_penalisation(const int Numero_face) const
{
  //  Cerr << "Entering coefficient_penalisation" << finl;
  /* Initialisation of local parameters */
  double eta=0.;
  double coefficientpenalisation = 0.;

  const int voisin1 =
    domaine_vef().face_voisins(Numero_face,1);

  const int voisin2 =
    domaine_vef().face_voisins(Numero_face,0);

  if (voisin1 == -1 && voisin2 == -1)
    {
      Cerr << "Error in Op_Dift_standard_Face_VEF_penalise::"
           << "coefficient_penalisation()" << finl;
      Process::exit();
    }

  if (voisin1 != -1 && voisin2 != -1)
    {
      /* Compute the penalisation coefficients */
      coefficientpenalisation = 1./diametre(voisin1);
      eta = 1./diametre(voisin2);

      coefficientpenalisation = std::min(coefficientpenalisation,eta);
    }

  if (voisin1 == -1)
    coefficientpenalisation = 1./diametre(voisin2);

  if (voisin2 == -1)
    coefficientpenalisation = 1./diametre(voisin1);

  return coefficientpenalisation;
  //  Cerr << "Exiting coefficient_penalisation" << finl;
}

/*
 *
 *
 *
 *
 *
 */

/*! @brief Member function returning the list of faces belonging to the neighbourhood of both Face1 AND Face2.
 *
 */
void  Op_Diff_VEF_Face_Penalise::
faces_communes(const int Face1,const int Face2,
               IntList& Faces_communes) const
{

  //  Cerr <<"Entering faces_communes" << finl;
  /* First clear Faces_communes to avoid errors. */
  if (! Faces_communes.est_vide() ) Faces_communes.vide();

  /* Declaration of local parameters of the procedure */
  IntList voisinage_Face1,voisinage_Face2;
  int nb_element_voisinage_Face2;

  /* Compute the neighbourhoods of both Face1 and Face2 */
  voisinage(Face1,voisinage_Face1);
  voisinage(Face2,voisinage_Face2);

  /* Then find the faces common to these 2 neighbourhoods. */
  for (nb_element_voisinage_Face2 = 0;
       nb_element_voisinage_Face2 < voisinage_Face2.size();
       nb_element_voisinage_Face2++)
    {
      const int numero_face_voisinage_Face2 =
        voisinage_Face2[nb_element_voisinage_Face2];

      if (voisinage_Face1.contient(numero_face_voisinage_Face2))
        Faces_communes.add_if_not(numero_face_voisinage_Face2);

    }// end for

  //  Cerr << "Exiting faces_communes" << finl;
}

/*
 *
 *
 *
 *
 *
 */

void  Op_Diff_VEF_Face_Penalise::
reduction(const IntList& Liste1,const IntList& Liste2,
          IntList& Liste_reduite) const
{
  //  Cerr << "Entering reduction" << finl;
  /* Clear Liste_reduite to avoid errors */
  if (! Liste_reduite.est_vide() ) Liste_reduite.vide();

  /* Declaration of local parameters of the procedure */
  const IntList *liste_de_plus_petite_taille,*liste_de_plus_grande_taille;
  int nb_element_dans_liste;

  /* Check the sizes of the lists passed as parameters, then
   * allocate accordingly.
   */
  if (Liste1.size() >= Liste2.size())
    {
      liste_de_plus_petite_taille = &Liste2;
      liste_de_plus_grande_taille = &Liste1;
    }
  else
    {
      liste_de_plus_petite_taille = &Liste1;
      liste_de_plus_grande_taille = &Liste2;
    }

  /* Remove from liste_de_plus_grande_taille the elements
   * of liste_de_plus_petite_taille that are present in it.
   */
  for (int ll = 0 ; ll < (*liste_de_plus_grande_taille).size() ; ll++)
    Liste_reduite.add( (*liste_de_plus_grande_taille)[ll] );

  for (nb_element_dans_liste = 0;
       nb_element_dans_liste < (*liste_de_plus_petite_taille).size();
       nb_element_dans_liste++)
    if
    (Liste_reduite.contient( (*liste_de_plus_petite_taille)
                             [nb_element_dans_liste] ) )
      Liste_reduite.suppr( (*liste_de_plus_petite_taille)
                           [nb_element_dans_liste] );

  //  Cerr << "Exiting reduction" << finl;
}


/*
 *
 *
 *
 *
 *
 */

/*! @brief Member function returning the index of the element containing both Face1 and Face2 if it exists, or -1 otherwise.
 *
 */
int Op_Diff_VEF_Face_Penalise::
element_commun(const int Face1,const int Face2) const
{
  //  Cerr << "Entering element_commun" << finl;

  /* Neighbouring elements of Face1 */
  const int voisin1_Face1 = domaine_vef().face_voisins(Face1,1);
  const int voisin2_Face1 = domaine_vef().face_voisins(Face1,0);

  /* Neighbouring elements of Face2 */
  const int voisin1_Face2 = domaine_vef().face_voisins(Face2,1);
  const int voisin2_Face2 = domaine_vef().face_voisins(Face2,0);

  /* Search for the common element */
  if (voisin1_Face1 != -1)
    if (voisin1_Face1 == voisin1_Face2 || voisin1_Face1 == voisin2_Face2)
      return voisin1_Face1;

  if (voisin2_Face1 != -1)
    if (voisin2_Face1 == voisin1_Face2 || voisin2_Face1 == voisin2_Face2)
      return voisin2_Face1;

  Cerr << " Op_Diff_VEF_Face_Penalise::element_commun()" << finl;
  Cerr << "Warning: face " << Face1 << " and face " << Face2
       << " have no common element" << finl;
  Cerr << "Exiting element_commun" << finl;
  return -1;

}

/*
 *
 *
 *
 *
 *
 */

/*! @brief Member function returning the 3rd face of element Element if Face1 and Face2 belong to the same element.
 *
 *  Returns -1 otherwise.
 *
 */
int Op_Diff_VEF_Face_Penalise::
autre_face(const int Face1, const int Face2)
const
{
  //  Cerr << "Entering autre_face" << finl;
  /* Declaration of local variables */
  int numero_local,lautre_face=-1;
  int elem_commun = element_commun(Face1,Face2);

  /* Number of faces per element in the discretisation domain.
   * NOTE: prismes are excluded.
   */
  if (elem_commun == -1)
    {
      Cerr << "Function element_commun" << finl;
      Cerr << "The 2 faces do not belong to the same element." << finl;
      Process::exit();
      return lautre_face;
    }

  const int nb_faces_element = domaine().nb_faces_elem();

  for (numero_local = 0; numero_local < nb_faces_element; numero_local++)

    {
      /* Retrieve the global index of each face of Element */
      const int numero_global_face =
        domaine_vef().elem_faces(elem_commun,numero_local);

      if ( numero_global_face != Face1 && numero_global_face != Face2)
        {
          lautre_face = numero_global_face;
          break;
        }

    }// end for

  return lautre_face;

}
/*
 *
 *
 *
 *
 *
 */


