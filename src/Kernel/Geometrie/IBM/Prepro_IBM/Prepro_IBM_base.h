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

#ifndef Prepro_IBM_base_included
#define Prepro_IBM_base_included

#include <Discretisation_base.h>
#include <Champ_Don_base.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Octree_Double.h>
#include <Champ_base.h>
#include <TRUST_Ref.h>
#include <Param.h>
#include <Debog.h>

using namespace std;

class Domaine_dis_base;
class Source_PDF_base;

class Prepro_IBM_base :public Objet_U
{
  Declare_base(Prepro_IBM_base);

public:
  virtual void associer_pb(const Probleme_base&);
  void discretiser() ;
  inline int verify_results_prepro() { return verify_results_prepro_; };
  void compute_h_max_elem();
  void compute_NeighNode(int);
  inline const DoubleTab& get_champ_aire() {return champ_aire_->valeurs();}
  inline const DoubleTab& get_champ_rotation() {return champ_rotation_->valeurs();}
  inline const DoubleTab& get_champ_barycentre() {return champ_bary_->valeurs();}
  inline const DoubleTab& get_champ_solid_points() {return solid_points_->valeurs();}
  inline const DoubleTab& get_champ_fluid_points() {return fluid_points_->valeurs();}
  inline const DoubleTab& get_champ_solid_elems() {return solid_elems_->valeurs();}
  inline const DoubleTab& get_champ_fluid_elems() {return fluid_elems_->valeurs();}
  inline const DoubleTab& get_champ_corresp_elems() {return corresp_elems_->valeurs();}
  inline const DoubleTab& get_isNodeDirichlet() {return isNodeDirichlet_->valeurs();}
  inline const DoubleTab& get_h_max_elem() {return h_max_elem_->valeurs();}
  inline const DoubleTab& get_h_max_node() {return h_max_node_->valeurs();}
  inline void set_champ_aire(DoubleTab& champ_aire) {champ_aire_->valeurs() = champ_aire;}
  inline void set_champ_rotation(DoubleTab& champ_rotation) {champ_rotation_->valeurs() = champ_rotation;}
  inline void set_champ_barycentre(DoubleTab& champ_bary) {champ_bary_->valeurs() = champ_bary;}
  void calculer_normal_proj_solid(DoubleTab&, DoubleTab&);
  virtual void compute_solid_fluid(int);

protected:
  void set_param(Param&) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void computeLocalFrame(const DoubleTab&, DoubleTab&, DoubleTab&, int);
  void computeMatRot(const DoubleTab&, DoubleTab&, DoubleTab&, int);
  void computeAire2();
  void Save_Med_File();
  void compute_effective_error();
  void intersectPolyPoly2D(MEDCouplingUMesh *, MEDCouplingUMesh *, DoubleTab&, IntTab&, double, int&);
  void intersectSegPoly2D(MEDCouplingUMesh *, DoubleTab&, DoubleTab&, double, MCAuto<MEDCoupling::DataArrayDouble>, int&);
  void intersectSegSeg2D(DoubleTab&, DoubleTab&, DoubleTab&, DoubleTab&, DoubleTab&, double, int&);
  OBS_PTR(Probleme_base) mon_pb_;
  Nom nom_fichier_med_IB_, nom_maillage_IB_ = "??";
  Domaine dom_med_IB_;
  MCAuto<MEDCoupling::MEDCouplingUMesh> aSkinUMesh_ = nullptr; // Mesh MedCoupling IBM
  DoubleTab barySurf_; // Barycenters of the Lagrangian mesh
  DoubleTab normalArr_; // Normals of the Lagrangian mesh
  DoubleTab coordsSur3D_; // Coordinates of the Lagrangian mesh

  double eps_ = 1.0e-12; // geometric precision
  double eps_effec_ = 1.0e-12; // effective geometric precision
  double c_prepro_ = 0.; // multiplicative factor for the search of the fluid point
  IntTab dimTab_ ; // choice of search directions for the fluid point
  bool save_prepro_ = false; // Save results to a MED file
  bool verify_results_prepro_ = false; // Verification of fields computed by prepro versus fields read from MED file (source_PDF)

  // Fields produced by the prepro_IBM
  OWN_PTR(Champ_Don_base) champ_rotation_, champ_aire_; // IBM rotation and area on the Eulerian mesh
  OWN_PTR(Champ_Don_base) champ_bary_; // IBM barycenters on the Eulerian mesh
  OWN_PTR(Champ_Don_base) champ_normal_; // IBM normals on the Eulerian mesh

  OWN_PTR(Champ_Don_base) isNodeDirichlet_ ;
  OWN_PTR(Champ_Don_base) fluid_points_;
  OWN_PTR(Champ_Don_base) fluid_elems_;
  OWN_PTR(Champ_Don_base) solid_points_;
  OWN_PTR(Champ_Don_base) solid_elems_;
  OWN_PTR(Champ_Don_base) corresp_elems_;

  OWN_PTR(Champ_Don_base) h_max_elem_;
  OWN_PTR(Champ_Don_base) h_max_node_;

  IntLists sommets_voisins_;

  Nom nom_fichier_med_Out_;

  int verbose_=0;

  friend class Source_PDF_base;
};

#endif /* Prepro_IBM_base_included */
