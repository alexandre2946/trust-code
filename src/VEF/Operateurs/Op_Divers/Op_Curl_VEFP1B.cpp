/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Op_Curl_VEFP1B.h>
#include <Domaine_Cl_VEF.h>
#include <Domaine_VEF.h>

Implemente_instanciable(Op_Curl_VEFP1B, "Op_Curl_VEFPreP1B_P1NC", Operateur_base);

Sortie& Op_Curl_VEFP1B::printOn(Sortie& s) const { return s << que_suis_je(); }
Entree& Op_Curl_VEFP1B::readOn(Entree& is) { return is; }

inline void add_curl_som(int nps, int sommet, int face, double flux, DoubleTab& curl, const Domaine& domaine)
{
  curl(nps + domaine.get_renum_som_perio(sommet)) += flux;
}

inline void traiter_flux(DoubleTab& curl, double flux, int element1, int element2, int npe)
{
  curl(npe + element1) += flux;
  curl(npe + element2) -= flux;
}

void Op_Curl_VEFP1B::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const Champ_Inc_base& inco)
{
  le_dom_vef = ref_cast(Domaine_VEF, domaine_dis);
  la_zcl_vef = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis);
  elements_pour_sommet();
}

DoubleTab& Op_Curl_VEFP1B::calculer(const DoubleTab& vitesse, DoubleTab& curl) const
{
  curl = 0;
  return ajouter(vitesse, curl);
}

DoubleTab& Op_Curl_VEFP1B::ajouter(const DoubleTab& vitesse, DoubleTab& curl) const
{
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  const Domaine& domaine = domaine_VEF.domaine();
  //int prems=domaine_VEF.premiere_face_int();
  if (dimension != 2)
    {
      Cerr << "Only 2D is supported at the moment. " << finl;
      Process::exit();
    }

  int face0 = 0, face1 = 0, face2 = 0;
  int numero_triangle = 0;
  int face_globale = 0, face_opp = 0;

  DoubleTab vecteur_normal0(dimension);
  DoubleTab vecteur_normal1(dimension);
  DoubleTab vecteur_normal2(dimension);

  // Process internal faces i.e. without boundary conditions
  // NOTE: the vorticity basis consists of: the set of indicator
  // functions of elements + the set of hat functions of P1
  // minus the last of these shape functions

  //P0 part of the vorticity
  for (int numero_elem = 0; numero_elem < domaine.nb_elem(); numero_elem++)
    {

      //REM: this part can be generalized to 3D
      //by using a loop with domaine.nb_faces_element()

      //First we need the indices of the 3 faces
      //belonging to element K
      face0 = domaine_VEF.elem_faces(numero_elem, 0);
      face1 = domaine_VEF.elem_faces(numero_elem, 1);
      face2 = domaine_VEF.elem_faces(numero_elem, 2);

      //Then we need the tangent vectors of
      //these three faces.
      vecteur_normal0 = vecteur_normal(face0, numero_elem);
      vecteur_normal1 = vecteur_normal(face1, numero_elem);
      vecteur_normal2 = vecteur_normal(face2, numero_elem);

      int modulo;
      for (int composante = 0; composante < dimension; composante++)
        {
          //The P0 part has been tested with functions (1,0);
          //(0,1);(x,0) and (0,x).
          //All results are correct

          modulo = (composante + 1) % 2;
          curl(numero_elem) += pow(-1., modulo)
                               * (vitesse(face0, composante) * vecteur_normal0(modulo) + vitesse(face1, composante) * vecteur_normal1(modulo) + vitesse(face2, composante) * vecteur_normal2(modulo));

        }

//      Cerr << "Element curl(" << numero_elem << ") " << curl(numero_elem) << finl;
    }

  //P1 part of the vorticity

  for (int numero_som = 0; numero_som < domaine.nb_som() - 1; numero_som++)
    {
      for (int num_loc_elem = 0; num_loc_elem < elem_som_size(numero_som); num_loc_elem++)

        {
          // for num_loc_elem

          //Retrieve the global index of the triangle
          numero_triangle = elements_pour_sommet(numero_som, num_loc_elem);

          //Retrieve the global index of the face opposite to "numero_som"
          //in triangle "numero_triangle"
          face_opp = domaine_VEF.numero_sommet_local(numero_som, numero_triangle);
          face_opp = domaine_VEF.elem_faces(numero_triangle, face_opp);

          //Retrieve the normal vector of this opposite face.
          vecteur_normal1 = vecteur_normal(face_opp, numero_triangle);

          for (int num_loc_face = 0; num_loc_face < domaine.nb_faces_elem(); num_loc_face++)

            {
              // for num_loc_face

              //Retrieve the global index of face "num_loc_face"
              face_globale = domaine_VEF.elem_faces(numero_triangle, num_loc_face);

              //               //If "face_globale" is an internal edge, perform the correct
              //               //traitement
              //               if (face_globale >= domaine_VEF.premiere_face_int() )
              {
                //Compute the normal vectors associated with these faces.
                vecteur_normal0 = vecteur_normal(face_globale, numero_triangle);

                //Finally compute the curl contribution of each of
                //these faces for each of the 2 triangles.
                int modulo;
                for (int composante = 0; composante < dimension; composante++)
                  {
                    //P1 part tested with functions (1,0);(0,1)
                    //(x,0) and (0,x).
                    //All tests are correct

                    modulo = (composante + 1) % 2;

                    //Partie (lambda_s,curl u)
                    curl(domaine.nb_elem() + numero_som) += -pow(-1., modulo) * 1. / (dimension + 1) * vitesse(face_globale, composante) * vecteur_normal0(modulo);

                    //Partie (rot lambda_s, u)
                    curl(domaine.nb_elem() + numero_som) += pow(-1., modulo) * 1. / (dimension * (dimension + 1)) * vitesse(face_globale, composante) * vecteur_normal1(modulo);
                  }

              } // end if

            } // end for num_loc_face

          /* For now, we only work with H10 velocities */
          /* This is the subject of the thesis */
          /* Consequently, no need to treat boundary faces */

        } // end for num_loc_elem

      Cerr << "Vertex curl(" << numero_som << ") " << curl(domaine.nb_elem() + numero_som) << finl;

    } // end for over vertices

  Cerr << "exiting OpCurl" << finl;

  return curl;
}

DoubleTab Op_Curl_VEFP1B::vecteur_normal(const int face, const int elem) const
{
  assert(dimension == 2);

  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  DoubleTab le_vecteur_normal(dimension);

  for (int composante = 0; composante < dimension; composante++)

    le_vecteur_normal(composante) = domaine_VEF.face_normales(face, composante) * domaine_VEF.oriente_normale(face, elem);

  return le_vecteur_normal;
}

// Array that stores at position "i" all mesh elements containing the vertex with global index "i"
int Op_Curl_VEFP1B::elements_pour_sommet()
{
  const Domaine& domaine = le_dom_vef->domaine();
  int numero_global_som;
  elements_pour_sommet_.dimensionner(domaine.nb_som());

  for (int numero_elem = 0; numero_elem < domaine.nb_elem(); numero_elem++)
    for (int numero_som_loc = 0; numero_som_loc < domaine.nb_som_elem(); numero_som_loc++)
      {
        numero_global_som = domaine.sommet_elem(numero_elem, numero_som_loc);
        elements_pour_sommet_[numero_global_som].add_if_not(numero_elem);
      }

  return 1;
}

// Function returning the global index of the element containing "sommet" located at position "indice" in the "elements_pour_sommet_" list
int Op_Curl_VEFP1B::elements_pour_sommet(const int sommet, const int indice) const
{
  return elements_pour_sommet_[sommet][indice];
}

// Function returning the size of the list at position "sommet" in the "elements_pour_sommet_" array
int Op_Curl_VEFP1B::elem_som_size(const int sommet) const
{
  return elements_pour_sommet_[sommet].size();
}
