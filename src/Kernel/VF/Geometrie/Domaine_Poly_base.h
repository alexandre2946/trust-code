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

#ifndef Domaine_Poly_base_included
#define Domaine_Poly_base_included

#include <Domaine_Poly_tools.h>
#include <Static_Int_Lists.h>
#include <Elem_poly_base.h>
#include <Elem_poly_base.h>
#include <TRUST_Deriv.h>
#include <TRUSTLists.h>
#include <Periodique.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Conds_lim.h>
#include <Domaine.h>
#include <Lapack.h>
#include <math.h>
#include <vector>
#include <string>
#include <array>
#include <map>

class Geometrie;
extern bool polymac_flica5;
/*! @brief class Domaine_Poly_base
 *
 *  	Instantiable class derived from Domaine_VF.
 *  	This class contains the geometric information required by the
 *  	Finite Volume Element method (Crouzeix-Raviart element).
 *  	The class holds a number of pieces of information concerning the faces.
 *  	In this set of faces, the boundary and joint faces are also included.
 *       To handle the faces, two categories are distinguished:
 *            - non-standard faces: located on a joint, a boundary, or internal
 *              but belonging to a boundary element
 *            - standard faces: internal faces not belonging to a boundary element
 *       This distinction corresponds to the treatment of boundary conditions:
 *       standard faces do not "see" boundary conditions.
 *       The set of faces is numbered as follows:
 *            - faces on a Domaine_joint appear first
 *     	       (in the order of the les_joints vector)
 *    	     - faces on a Domaine_bord appear next
 * 	       (in the order of the les_bords vector)
 *   	     - non-standard internal faces appear next
 *            - standard internal faces appear last
 *       Consequently, all non-standard faces requiring special treatment are
 *       grouped together at the beginning.
 *       Two types of elements are distinguished:
 *            - non-standard elements: they have at least one boundary face
 *            - standard elements: they have no boundary face
 *       Standard elements (resp. non-standard elements) are not stored
 *       consecutively in the Domaine object. The rang_elem_non_std array
 *       is used to selectively access one or the other type of element.
 *
 */
class Domaine_Poly_base : public Domaine_VF
{
  Declare_base(Domaine_Poly_base);
public :
  void typer_elem(Domaine& domaine_geom) override;
  void discretiser() override;
  virtual void calculer_volumes_entrelaces() { }
  void discretiser_aretes();

  void orthocentrer();

  inline const DoubleVect& longueur_aretes() const { return longueur_aretes_; }
  inline const DoubleTab& ta() const { return ta_; }

  void modifier_pour_Cl(const Conds_lim& ) override;

  inline const Elem_poly_base& type_elem() const { return type_elem_.valeur(); }
  inline int nb_elem_Cl() const { return nb_elem() - nb_elem_std_; }
  inline int nb_faces_joint() const { return 0; /* return nb_faces_joint_;    A FAIRE */ }
  inline int nb_faces_std() const { return nb_faces_std_; }
  inline int nb_elem_std() const { return nb_elem_std_; }
  inline double carre_pas_du_maillage() const { return h_carre; }
  inline double carre_pas_maille(int i) const { return h_carre_(i); }
  inline IntVect& rang_elem_non_std() { return rang_elem_non_std_; }
  inline const IntVect& rang_elem_non_std() const { return rang_elem_non_std_; }

  virtual void calculer_h_carre();

  inline DoubleTab& volumes_entrelaces_dir() { return volumes_entrelaces_dir_; } // returns the interlaced volumes array per side.
  inline const DoubleTab& volumes_entrelaces_dir() const { return volumes_entrelaces_dir_; }

  //equivalent to dot(), but for the product (a - ma).nu.(b - mb)
  inline double nu_dot(const DoubleTab* nu, int e, int n, const double *a, const double *b, const double *ma = nullptr, const double *mb = nullptr) const;

  inline double dist_norm(int num_face) const override;
  inline double dist_norm_bord(int num_face) const override;
  DoubleVect& dist_norm_bord(DoubleVect& , const Nom& nom_bord) const;
  inline double dist_face_elem0(int num_face,int n0) const override;
  inline double dist_face_elem1(int num_face,int n1) const override;
  inline double dist_face_elem0_period(int num_face,int n0,double l) const override;
  inline double dist_face_elem1_period(int num_face,int n1,double l) const override;

  void detecter_faces_non_planes() const;

  //equivalent faces: equiv(f, 0/1, i) = face equivalent to e_f(f_e(f, 0/1), i) on the other side, -1 if there is none
  const IntTab& equiv() const;
  virtual void init_equiv() const = 0;

  //connectivite sommet-elements
  const Static_Int_Lists& som_elem() const;

  //indexing in arrays of type (element, vertex) and (element, edge)
  const IntTab& elem_som_d() const; //entry of vertex les_elems(e, i) of element e: elem_som_d()(e) + i
  const IntTab& elem_arete_d() const; //entry of edge elem_arete(e, i) of element e: elem_arete_d()(e) + i

  //for each element, distribution of its volume among each of its vertices
  const DoubleTab& vol_elem_som() const;
  //for each vertex, product porosity * volume
  const DoubleTab& pvol_som(const DoubleVect& poro) const;

  //som_arete[som1][som2 > som1] -> arete correspondant a (som1, som2)
  std::vector<std::map<int, int> > som_arete;

  //MD_Vectors for Champ_Elem_PolyMAC_CDO (elems + faces) and for Champ_Face_PolyMAC_CDO (faces + edges)
  mutable MD_Vector mdv_elems_faces, mdv_faces_aretes;

  void calculer_infos_aretes();
  void fill_normales();
  void recalculer_xv();

protected:
  void verifier_type_elem() const;
  void corriger_face_voisins_sur_les_faces_virtuelles();

  double h_carre = DMAXFLOAT;			 // square of the mesh step size
  DoubleVect h_carre_;			// square of the mesh cell step size
  OWN_PTR(Elem_poly_base) type_elem_;                  // type of the discretization element

  Sortie& ecrit(Sortie& os) const;

  mutable IntTab equiv_;
  mutable Static_Int_Lists som_elem_;
  mutable IntTab elem_som_d_, elem_arete_d_;
  mutable DoubleTab vol_elem_som_, pvol_som_;

  DoubleVect longueur_aretes_; //edge lengths
  mutable DoubleTab ta_;       //tangent vectors to the edges
};

/* equivalent of dist_norm_bord from VDF */
inline double Domaine_Poly_base::dist_norm_bord(int f) const
{
  assert(face_voisins(f, 1) == -1);
  return std::fabs(dot(&xp_(face_voisins(f, 0), 0), &face_normales_(f, 0), &xv_(f, 0))) / face_surfaces(f);
}

inline double Domaine_Poly_base::dist_norm(int f) const
{
  return std::fabs(dot(&xp_(face_voisins(f, 0), 0), &face_normales_(f, 0), &xp_(face_voisins(f, 1), 0))) / face_surfaces(f);
}

inline double Domaine_Poly_base::dist_face_elem0(int f,int e) const
{
  return std::fabs(dot(&xp_(e, 0), &face_normales_(f, 0), &xv_(f, 0))) / face_surfaces(f);
}

inline double Domaine_Poly_base::dist_face_elem1(int f,int e) const
{
  return std::fabs(dot(&xp_(e, 0), &face_normales_(f, 0), &xv_(f, 0))) / face_surfaces(f);
}

inline double Domaine_Poly_base::dist_face_elem0_period(int num_face,int n0,double l) const
{
  abort();
  return 0;
}

inline double Domaine_Poly_base::dist_face_elem1_period(int num_face,int n1,double l) const
{
  abort();
  return 0;
}

//returns the dot product a.nu.b regardless of the number of components and the type of nu tensor
inline double Domaine_Poly_base::nu_dot(const DoubleTab* nu, int e, int n, const double *a, const double *b, const double *ma, const double *mb) const
{
  if (!nu) return dot(a, b, ma, mb);
  int d, db, D = dimension;
  double resu = 0;
  if (nu->nb_dim() == 2) resu += (*nu)(e, n) * dot(a, b, ma, mb); //isotrope
  else if (nu->nb_dim() == 3)
    for (d = 0; d < D; d++) //anisotrope diagonal
      resu += (*nu)(e, n, d) * (a[d] - (ma ? ma[d] : 0)) * (b[d] - (mb ? mb[d] : 0));
  else for (d = 0; d < D; d++)
      for (db = 0; db < D; db++)
        resu += (*nu)(e, n, d, db) * (a[d] - (ma ? ma[d] : 0)) * (b[db] - (mb ? mb[db] : 0));
  return resu;
}

#endif /* Domaine_Poly_base_included */
