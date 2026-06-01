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

#ifndef Op_Diff_RotRot_included
#define Op_Diff_RotRot_included

#include <Operateur_Diff_base.h>
#include <Op_Curl_VEFP1B.h>
#include <Equation_base.h>
#include <Op_Rot_VEFP1B.h>
#include <Champ_base.h>
#include <SolveurSys.h>
#include <TRUST_Ref.h>

class Champ_Uniforme;
class Domaine_VEF;

class Op_Diff_RotRot: public Operateur_Diff_base
{
  Declare_instanciable(Op_Diff_RotRot);
public:
  void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base&) override;
  void associer_diffusivite(const Champ_base&) override;
  const Champ_base& diffusivite() const override;
  DoubleTab& calculer(const DoubleTab&, DoubleTab&) const override;
  DoubleTab& ajouter(const DoubleTab&, DoubleTab&) const override;

  int tester() const;

  /* Function to compute the vorticity from the velocity field. */
  int calculer_vorticite(DoubleTab&, const DoubleTab&) const;

  /* Assemble the vorticity matrix. */
  int assembler_matrice(Matrice&);

  /* For a given element "numero_elem", returns the list of vertices belonging to that element. */
  IntList sommets_pour_element(int numero_elem) const;

  /* For a given vertex "numero_som", returns the list of elements containing that vertex. */
  IntList elements_pour_sommet(int numero_som) const;

  /* For a given vertex "numero_som", returns */
  /* the list of neighboring vertices of "numero_som". */
  /* Parameter: the list of elements containing "numero_som". */
  /* It is sufficient to search in those elements to get the result. */
  /* NOTE: the result list includes the vertex "numero_som" itself. */
  IntList sommets_voisins(int numero_som, const IntList& liste) const;

  /* For element "numero_elem", returns the coefficient */
  /* to place in the nb_elem * nb_elem sub-matrix */
  /* at row "numero_elem". */
  /* EF matrix */
  double remplir_elem_elem_EF(const int numero_elem) const;

  /* For element "numero_elem", returns the coefficient */
  /* to place in the nb_elem * nb_som sub-matrix */
  /* at row "numero_elem". */
  /* EF matrix */
  double remplir_elem_som_EF(const int numero_elem, const int numero_som) const;

  /* For vertex "numero_som", returns the coefficient */
  /* to place in the nb_som * nb_elem sub-matrix */
  /* at row "numero_som". */
  /* EF matrix */
  double remplir_som_elem_EF(const int numero_elem, const int numero_som) const;

  /* For element "numero_elem", returns the coefficient */
  /* to place in the nb_elem * nb_som sub-matrix */
  /* at row "numero_som". */
  /* "IntList" is the array of elements containing "numero_som". */
  /* EF matrix */
  double remplir_som_som_EF(const int numero_som, const int sommet_voisin, const IntList&) const;

  /* Function to sort an IntList */
  /* Sorting is done in ascending order. */
  void Tri(IntList& liste_a_trier) const;

  //Method to return the outward normal vector at "face" of element "elem"
  DoubleTab vecteur_normal(const int face, const int elem) const;

  const Domaine_VEF& domaine_vef() const;

protected:

  mutable OWN_PTR(Champ_Inc_base) vorticite_;
  Matrice matrice_vorticite_;

  Op_Rot_VEFP1B rot_;
  Op_Curl_VEFP1B curl_;

  SolveurSys solveur_;
  //mutable Solv_GCP solveur_;

  OBS_PTR(Domaine_VEF) le_dom_vef;
  OBS_PTR(Domaine_Cl_VEF) la_zcl_vef;
  OBS_PTR(Champ_Uniforme) diffusivite_;
};

#endif
