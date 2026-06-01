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

#ifndef Domaine_VEF_included
#define Domaine_VEF_included

#include <TRUSTArray_kokkos.tpp>
#include <Elem_VEF_base.h>
#include <TRUST_Deriv.h>
#include <Domaine_VF.h>
#include <kokkos++.h>

class VEF_discretisation;
class Geometrie;

/*! @brief class Domaine_VEF
 *
 * @brief Instantiable class derived from Domaine_VF.
 *          This class contains the geometric information required by the Finite Element Volume method (Crouzeix-Raviart element).
 *          The class holds a number of pieces of information about faces.
 *          Among these faces, boundary and joint faces are also included. Two categories of faces are distinguished:
 *            - non-standard faces: on a joint, a boundary, or internal faces
 *              belonging to a boundary element
 *            - standard faces: internal faces not belonging to any boundary element
 *       This distinction corresponds to boundary condition treatment: standard faces do not "see" boundary conditions.
 *       The full set of faces is numbered as follows:
 *            - faces on a Domaine_joint appear first
 *                    (in the order of the les_joints vector)
 *                 - faces on a Domaine_bord appear next
 *                (in the order of the les_bords vector)
 *                - internal non-standard faces appear next
 *            - internal standard faces appear last
 *       All non-standard faces requiring special treatment are therefore grouped at the beginning.
 *       Two element types are distinguished:
 *            - non-standard elements: they have at least one boundary face
 *            - standard elements: they have no boundary face
 *       Standard (resp. non-standard) elements are not stored consecutively in the Domaine object.
 *       The array rang_elem_non_std is used to selectively access either type.
 */
class Domaine_VEF: public Domaine_VF
{
  Declare_instanciable(Domaine_VEF);
public:
  void discretiser() override;
  virtual void discretiser_suite(const VEF_discretisation&);
  void discretiser_arete();

  void construire_ok_arete();
  int lecture_ok_arete();
  void verifie_ok_arete(int) const;
  void construire_renum_arete_perio(const Conds_lim&);

  virtual void creer_tableau_p1bulle(Array_base&, RESIZE_OPTIONS opt = RESIZE_OPTIONS::COPY_INIT) const;

  void swap(int, int, int);
  void modifier_pour_Cl(const Conds_lim&) override;
  void typer_elem(Domaine&) override;
  void calculer_volumes_entrelaces();
  void calculer_h_carre();
  DoubleTab& vecteur_face_facette();

  inline const Elem_VEF_base& type_elem() const { return type_elem_.valeur(); }
  inline int nb_elem_Cl() const { return nb_elem() - nb_elem_std_; }
  inline int nb_faces_joint() const { return 0; }
  inline int nb_faces_std() const { return nb_faces_std_; }
  inline int nb_elem_std() const { return nb_elem_std_; }
  inline int premiere_face_std() const { return nb_faces() - nb_faces_std_; }
  inline int nb_faces_non_std() const { return nb_faces() - nb_faces_std_; }
  inline double carre_pas_du_maillage() const { return h_carre; }
  inline const DoubleVect& carre_pas_maille() const { return h_carre_; }
  inline auto& facette_normales() { return facette_normales_; }
  inline const auto& facette_normales() const { return facette_normales_; }
  inline IntVect& rang_elem_non_std() { return rang_elem_non_std_; }
  inline const IntVect& rang_elem_non_std() const { return rang_elem_non_std_; }

  inline double volume_au_sommet(int som) const { return volumes_som_[som]; }
  inline const DoubleVect& volume_aux_sommets() const { return volumes_som_; }
  inline int get_P1Bulle() const { assert(P1Bulle != -1); return P1Bulle; }
  inline int get_alphaE() const { assert(alphaE != -1); return alphaE; }
  inline int get_alphaS() const { assert(alphaS != -1); return alphaS; }
  inline int get_alphaA() const { assert(alphaA != -1); return alphaA; }
  inline int get_alphaRT() const { assert(alphaRT != -1); return alphaRT; }
  inline int get_modif_div_face_dirichlet() const { assert(modif_div_face_dirichlet != -1); return modif_div_face_dirichlet; }
  inline int get_cl_pression_sommet_faible() const { assert(cl_pression_sommet_faible != -1); return cl_pression_sommet_faible; }
  inline const ArrOfInt& get_renum_arete_perio() const { return renum_arete_perio; }
  inline const IntVect& get_ok_arete() const { return ok_arete; }
  inline const DoubleVect& get_volumes_aretes() const { return volumes_aretes; }

  inline virtual const MD_Vector& md_vector_p1b() const { assert(md_vector_p1b_); return md_vector_p1b_; }

  inline int numero_premier_element() const;
  inline int numero_premier_sommet() const;
  inline int numero_premiere_arete() const;

  inline double dist_face_elem0(int num_face,int n0) const override;
  inline double dist_face_elem1(int num_face,int n1) const override;

private:
  double h_carre = 1.e30;                         // squared mesh step size
  DoubleVect h_carre_;                        // squared cell step size
  OWN_PTR(Elem_VEF_base) type_elem_;                  // type of the discretisation element
  // normals to interlaced volume faces:
#ifdef TRUST_USE_GPU
  BigDoubleTab facette_normales_; // Cause size=nb_elem*6*dim may be > 2^31
#else
  DoubleTab facette_normales_;
#endif
  DoubleTab vecteur_face_facette_;                // vector from face center to facette center
  IntVect orientation_;


  DoubleVect volumes_som_, volumes_aretes;
  ArrOfInt renum_arete_perio;
  IntVect ok_arete;

  int P1Bulle = -1, alphaE = -1, alphaS = -1, alphaA = -1;
  int alphaRT = -1; // for trio stationary
  int modif_div_face_dirichlet = -1;
  int cl_pression_sommet_faible = -1; // determines whether pressure BCs are imposed weakly or strongly -> see divergence and assembler
  // Descriptor for p1b arrays (depending on alphaE, alphaS and alphaA) (built in Domaine_VEF::discretiser())
  MD_Vector md_vector_p1b_;

  Sortie& ecrit(Sortie& os) const;
};

// Out-of-class Kokkos function: otherwise dom_VEF.oriente_normale(...) would copy a Domaine_VEF instance from host to device!
KOKKOS_INLINE_FUNCTION int oriente_normale(int face_opp, int elem2, CIntTabView face_voisins)
{
  return (face_voisins(face_opp, 0) == elem2) ? 1 : -1;
}

// Test method:
void exemple_champ_non_homogene(const Domaine_VEF&, DoubleTab&);

inline int Domaine_VEF::numero_premier_element() const
{
  if (!alphaE)
    return -1;
  else
    return 0;
}
inline int Domaine_VEF::numero_premier_sommet() const
{
  if (!alphaS)
    return -1;
  else if (!alphaE)
    return 0;
  else
    return nb_elem_tot();
}
inline int Domaine_VEF::numero_premiere_arete() const
{
  if (!alphaA)
    return -1;
  else if (!alphaE && !alphaS)
    return 0;
  else if (!alphaE && alphaS)
    return nb_som_tot();
  else if (!alphaS && alphaE)
    return nb_elem_tot();
  else
    return nb_elem_tot() + nb_som_tot();
}

inline double Domaine_VEF::dist_face_elem0(int f,int e) const
{
  return std::fabs(dot(&xp_(e, 0), &face_normales_(f, 0), &xv_(f, 0))) / face_surfaces(f);
}

inline double Domaine_VEF::dist_face_elem1(int f,int e) const
{
  return std::fabs(dot(&xp_(e, 0), &face_normales_(f, 0), &xv_(f, 0))) / face_surfaces(f);
}
#endif
