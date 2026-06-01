/****************************************************************************
* Copyright (c) 2023, CEA
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


#ifndef Op_Diff_VEF_Face_Penalise_included
#define Op_Diff_VEF_Face_Penalise_included
#include <Op_Diff_VEF_Face.h>
#include <Domaine_Cl_VEF.h>
#include <TRUSTList.h>
#include <Domaine.h>

class Op_Diff_VEF_Face_Penalise : public Op_Diff_VEF_Face
{

  Declare_instanciable(Op_Diff_VEF_Face_Penalise);

public :
  /* Member function computing the penalised diffusion.
   * Overrides the base class version.
   */
  DoubleTab& ajouter(const DoubleTab& inconnue, DoubleTab& resu) const override;
  DoubleTab& calculer(const DoubleTab& inconnue, DoubleTab& resu) const override;

protected:

  /* Member function returning in the list Voisinage the set
   * of face indices neighbouring face Numero_face.
   * Neighbouring faces are those belonging to elements
   * that contain face Numero_face.
   * Note: Voisinage also contains face Numero_face itself.
   */
  void voisinage(const int Numero_face, IntList& Voisinage) const;

  /* Member function returning in the list Voisinage
   * the set of faces constituting the neighbourhood of the list
   * Ensemble_faces.
   * Note: Voisinage also contains the elements of Ensemble_faces.
   */
  void voisinage(const IntList& Ensemble_faces, IntList& Voisinage) const;

  /* Member function returning the appropriate sign for
   * computing the interactions of the basis functions of faces
   * Face1 and Face2 along another face.
   */
  double signe(const int Face1, const int Face2 ) const;

  /* Member function returning the penalisation coefficient
   * to be applied on each edge of the primary mesh.
   */
  double coefficient_penalisation(const int Numero_face) const;

  /* Member function returning the list Faces_communes of faces
   * that belong to the neighbourhood of both Face1 AND Face2.
   */
  void faces_communes(const int Face1,const int Face2,
                      IntList& Face_commune) const;

  /* Member function returning the list Liste_reduite which is
   * the set difference of Liste1 and Liste2.
   * Note: the function dynamically checks the lengths of
   * Liste1 and Liste2 before executing.
   */
  void reduction(const IntList& Liste1,const IntList& Liste2,
                 IntList& Liste_reduite) const;

  /* Member function returning the index of the element shared by
   * Face1 and Face2 if it exists, or -1 otherwise.
   */
  int element_commun(const int Face1,const int Face2) const;

  /* Member function returning the index of the 3rd face
   * of an element, if Face1 and Face2 belong to the same element.
   * Returns -1 otherwise.
   */
  int autre_face(const int Face1, const int Face2) const;

  /* Member function returning the diameter of element Element.
   * Note: assumed to operate in 2D.
   */
  inline double diametre(const int Element) const;

  /* Member function returning the length of face Face.
   * Note: assumed to operate in 2D.
   */
  inline double longueur(const int Face) const;


private:
  /* Member function returning the Domaine_VEF of the domain. */
  inline const Domaine_VEF& domaine_vef() const;

  /* Member function returning the Domaine of the problem. */
  inline const Domaine&  domaine() const;

  /* Member function returning the boundary conditions domain. */
  inline const Domaine_Cl_VEF& domaine_cl() const;

};




inline double Op_Diff_VEF_Face_Penalise::longueur(const int Face) const
{
  double x_sommet1,x_sommet2;
  double y_sommet1,y_sommet2;

  int sommet1 = domaine_vef().face_sommets(Face,0);
  int sommet2 = domaine_vef().face_sommets(Face,1);

  x_sommet1 = domaine().coord(sommet1,0);
  y_sommet1 = domaine().coord(sommet1,1);

  x_sommet2 = domaine().coord(sommet2,0);
  y_sommet2 = domaine().coord(sommet2,1);

  return sqrt( pow(x_sommet1 - x_sommet2,2) + pow(y_sommet1 - y_sommet2,2) );
}

inline double Op_Diff_VEF_Face_Penalise::diametre(const int Element) const
{
  int face1,face2,face3;
  double longueur_face1,longueur_face2,longueur_face3;
  double longueur_max,diametre_element;

  face1 = domaine_vef().elem_faces(Element,0);
  face2 = domaine_vef().elem_faces(Element,1);
  face3 = domaine_vef().elem_faces(Element,2);

  longueur_face1 = longueur(face1);
  longueur_face2 = longueur(face2);
  longueur_face3 = longueur(face3);

  longueur_max = longueur_face1 >= longueur_face2 ? longueur_face1 :
                 longueur_face2;
  diametre_element = longueur_face3 >= longueur_max ? longueur_face3 :
                     longueur_max;

  return diametre_element;
}

inline const Domaine_VEF& Op_Diff_VEF_Face_Penalise::domaine_vef() const
{
  return le_dom_vef.valeur();
}

inline const Domaine_Cl_VEF& Op_Diff_VEF_Face_Penalise::domaine_cl() const
{
  return la_zcl_vef.valeur();
}

inline const Domaine&  Op_Diff_VEF_Face_Penalise::domaine() const
{
  return domaine_vef().domaine();
}


#endif
