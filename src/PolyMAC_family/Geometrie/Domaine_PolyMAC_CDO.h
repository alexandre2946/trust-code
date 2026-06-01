
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

#ifndef Domaine_PolyMAC_CDO_included
#define Domaine_PolyMAC_CDO_included

#include <Echange_global_impose.h>
#include <Neumann_sortie_libre.h>
#include <Domaine_Poly_base.h>
#include <Matrice_Morse_Sym.h>
#include <Neumann_homogene.h>
#include <SolveurSys.h>
#include <Periodique.h>
#include <Dirichlet.h>
#include <Symetrie.h>

class Domaine_PolyMAC_CDO : public Domaine_Poly_base
{
  Declare_instanciable(Domaine_PolyMAC_CDO);
public :
  void discretiser() override;
  void swap(int, int, int) { }
  void modifier_pour_Cl(const Conds_lim& ) override;
  void init_equiv() const override;

  inline const IntTab& arete_faces() const { return arete_faces_; }
  void calculer_volumes_entrelaces() override;
  void calculer_h_carre() override;

  inline double dot (const double *a, const double *b, const double *ma = nullptr, const double *mb = nullptr) const { return dot(dimension, a, b, ma, mb); }
  KOKKOS_INLINE_FUNCTION double dot (const int dim, const double *a, const double *b, const double *ma = nullptr, const double *mb = nullptr) const;

  IntVect cyclic; // cyclic(i) = 1 if poly i is cyclic

  //which optional structures have been initialized
  mutable std::map<std::string, int> is_init;
  //first-order interpolations of the velocity vector at elements
  void init_ve() const;
  mutable IntTab vedeb, veji; //reconstruction of ve via (veji, veci)[vedeb(e), vedeb(e + 1)[ (faces)
  mutable DoubleTab veci;

  //curl at faces of a field tangent to edges
  void init_rf() const;
  mutable IntTab rfdeb, rfji; //reconstruction of the curl via (rfji, rfci)[rfdeb(f), rfdeb(f + 1)[ (field at edges)
  mutable DoubleTab rfci;

  //stabilization of a mimetic mass matrix in an element: in PolyMAC_CDO -> m1 or m2
  inline void ajouter_stabilisation(DoubleTab& M, DoubleTab& N) const;
  inline int W_stabiliser(DoubleTab& W, DoubleTab& R, DoubleTab& N, int *ctr, double *spectre) const;

  //mimetic matrix of a face field: (normal value at faces) -> (linear integral on broken lines)
  void init_m2() const;
  // ToDo switch these arrays to IntVect and DoubleVect:
  mutable IntTab m2d, m2i, m2j, w2i, w2j; //storage: rows of M_2^e in m2i([m2d(e), m2d(e + 1)[), indices/coeffs of these rows in (m2j/m2c)[m2i(i), m2i(i+1)[
  mutable DoubleTab m2c, w2c;             //          with the diagonal coefficient first (facilitates Echange_contact_PolyMAC_CDO)
  void init_m2solv() const; //to solve m2.v = s
  mutable Matrice_Morse_Sym m2mat;
  mutable SolveurSys m2solv;

  //first-order interpolation at elements of a field defined by its tangential components at edges (e.g. vorticity)
  inline void init_we() const;
  void init_we_2d() const;
  void init_we_3d() const;
  mutable IntTab wedeb, weji; //reconstruction of we via (weji, weci)[wedeb(e), wedeb(e + 1)[ (vertices in 2D, edges in 3D)
  mutable DoubleTab weci;

  //mimetic matrix of an edge field: (tangential value at edges) -> (flux through the union of facets touching the edge)
  void init_m1() const;
  void init_m1_2d() const;
  void init_m1_3d() const;
  mutable IntTab m1deb, m1ji; //reconstruction of m1 via (m1ji(.,0), m1ci)[m1deb(a), m1deb(a + 1)[ (vertices in 2D, edges in 3D); m1ji(.,1) contains the element index
  mutable DoubleTab m1ci;

  //std::map to retrieve the (proc, local item) pair associated with a virtual item for mdv_elem_faces
  void init_virt_ef_map() const;
  mutable std::map<std::array<int, 2>, int> virt_ef_map;

  //local matrices per element (Hodge operators) for performing interpolations:
  void M2(const DoubleTab *nu, int e, DoubleTab& m2) const; //normals at faces -> tangentials at dual faces:  (nu x_ef.v)    = m2 (|f|n_ef.v)
  void W2(const DoubleTab *nu, int e, DoubleTab& w2) const; //tangentials at dual faces -> normals at faces: (nu |f|n_ef.v) = w2 (x_ef.v)

private:
  void init_m2_new() const;
  void init_m2_osqp() const;

  mutable IntTab arete_faces_; //connectivite face -> aretes
};

/* dot product of two vectors */
KOKKOS_INLINE_FUNCTION double Domaine_PolyMAC_CDO::dot(const int dim, const double *a, const double *b, const double *ma, const double *mb) const
{
  double res = 0;
  for (int i = 0; i < dim; i++) res += (a[i] - (ma ? ma[i] : 0)) * (b[i] - (mb ? mb[i] : 0));
  return res;
}

#endif /* Domaine_PolyMAC_CDO_included */
