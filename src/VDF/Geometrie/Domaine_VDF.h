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

#ifndef Domaine_VDF_included
#define Domaine_VDF_included


#include <math.h>
#include <map>
#include <array>
#include <Domaine_VF.h>
#include <Domaine.h>
class Geometrie;

/*! @brief class Domaine_VDF
 *
 * @brief Instantiable class derived from Domaine_VF.
 *          This class contains the geometric information required by the
 *          Finite Difference Volume (VDF) method.
 *          The class holds a number of face-related data.
 *       All faces are numbered as follows:
 *            - faces on a Domaine_joint appear first
 *                    (in the order of the les_joints vector)
 *                 - faces on a Domaine_bord appear next
 *                (in the order of the les_bords vector)
 *                - internal faces appear last
 *       Each face is associated with an integer indicating its orientation.
 *       Within each face family (boundary, joint, internal) it is assumed that:
 *            - the block of faces with equation x = const (orientation 0) comes first
 *            - the block of faces with equation y = const (orientation 1) comes next
 *            - the block of faces with equation z = const (orientation 2) comes last
 *       For the boundary face block, sub-blocks corresponding to each boundary are preserved.
 *       No particular element numbering is required.
 *       The notion of edge is introduced for computing diffusive and convective fluxes
 *       in the momentum conservation equation.
 *       The Qdm array contains the edge/face connectivity. In this array,
 *       edges appear in the following order:
 *            - block of joint edges
 *            - block of boundary edges
 *            - block of mixed edges
 *            - block of internal edges
 *       Within each block, edges appear in the order: XY edges, XZ edges, YZ edges.
 *
 *
 *
 */

class Domaine_VDF : public Domaine_VF
{

  Declare_instanciable(Domaine_VDF);

public :

  void discretiser() override;
  Faces* creer_faces() override;
  inline int nb_faces_X() const;
  inline int nb_faces_Y() const;
  inline int nb_faces_Z() const;
  inline int nb_aretes() const;
  inline int nb_aretes_joint() const;
  inline int nb_aretes_coin() const;
  inline int premiere_arete_coin() const;
  inline int nb_aretes_bord() const;
  inline int premiere_arete_bord() const;
  inline int nb_aretes_mixtes() const;
  inline int premiere_arete_mixte() const;
  inline int nb_aretes_internes() const;
  inline int premiere_arete_interne() const;
  inline double h_x() const;
  inline double h_y() const;
  inline double h_z() const;

  inline int Qdm(int num_arete,int ) const;
  //inline double porosite_face(int ) const;
  //inline double porosite_elem(int i) const;
  using Domaine_VF::face_normales;
  inline double face_normales(int , int ) const override;
  inline int orientation(int ) const override;
  inline double dist_face(int , int , int k) const;
  inline double dist_norm(int num_face) const override;
  inline double dist_norm_bord(int num_face) const override;
  inline double dist_face_elem0(int ,int ) const override;
  inline double dist_face_elem1(int ,int ) const override;
  inline double dist_face_axi(int , int , int k) const;
  inline double dist_face_period(int , int , int ) const;
  inline double dist_norm_period(int ,double ) const;
  inline double dist_face_elem0_period(int ,int ,double ) const override;
  inline double dist_face_elem1_period(int ,int ,double ) const override;
  inline double dist_norm_axi(int num_face) const ;
  inline double dist_norm_bord_axi(int num_face) const;
  inline double dist_face_elem0_axi(int ,int ) const;
  inline double dist_face_elem1_axi(int ,int ) const;
  inline double distance_face(int , int , int k) const;
  inline double distance_normale(int num_face) const;
  inline double dist_elem(int ,int ,int ) const;
  inline double dist_elem_period(int ,int ,int ) const;
  inline double dim_elem(int ,int ) const;
  inline double dim_face(int ,int ) const;
  inline double delta_C(int ) const;
  inline int amont_amont(int, int ) const;
  inline int face_amont_princ(int ,int ) const;
  inline int face_amont_conj(int ,int ,int ) const;
  inline int face_bord_amont(int ,int ,int ) const;
  inline int elem_voisin(int , int , int ) const;

  inline IntVect& orientation();
  inline const IntVect& orientation() const override;

  // inline DoubleVect& porosite_face();
  // inline const DoubleVect& porosite_face() const;
  // inline DoubleVect& porosite_elem();
  // inline const DoubleVect& porosite_elem() const;

  inline IntTab& Qdm();
  inline const IntTab& Qdm() const;
  void calculer_volumes_entrelaces();
  void modifier_pour_Cl(const Conds_lim& cl) override;
  void creer_elements_fictifs(const Domaine_Cl_dis_base& ) override;
  DoubleVect& dist_norm_bord(DoubleVect& , const Nom& nom_bord) const;

  //std::map to retrieve the (proc, local item) pair associated with a virtual item for mdv_elem
  void init_virt_e_map() const;
  mutable std::map<std::array<int, 2>, int> virt_e_map;

protected:
  void prepare_elem_non_std(Faces&) override;
  void compute_sort_key(Faces&, IntTab& sort_key) override;
  void renumber_faces(Faces& les_faces, IntTab& sort_key) override;

private:

  IntVect orientation_;                    // face orientation
  // 0 if face perpendicular to the X axis
  // 1 if face perpendicular to the Y axis
  // 2 if face perpendicular to the Z axis
  int nb_faces_X_ = -1;                         // number of faces perpendicular to the X axis
  int nb_faces_Y_ = -1;                         // number of faces perpendicular to the Y axis
  int nb_faces_Z_ = -1;                         // number of faces perpendicular to the Z axis
  int nb_aretes_ = -1;                          // total number of edges of all types
  int nb_aretes_joint_ = -1;                    // number of joint edges
  int nb_aretes_coin_ = -1;                     // number of corner edges
  int nb_aretes_bord_ = -1;                     // number of boundary edges
  int nb_aretes_mixtes_ = -1;                   // number of mixed edges
  int nb_aretes_internes_ = -1;                 // number of internal edges
  IntTab Qdm_;                            // edge/face connectivity
  // DoubleVect porosite_elem_;               // volumetric porosities for mass control volumes
  // DoubleVect porosite_face_;               // surface porosities for mass and volumetric momentum control
  double h_x_ = 1.e30 , h_y_ = 1.e30 ,h_z_ = 1.e30;                   // mesh step sizes in the three spatial directions
  // h_x_ (resp. h_y_) is the smallest distance between two face centers
  // with equation X = const (resp. Y = const)

  // void calculer_porosites();
  void genere_aretes();
  void calcul_h();
  void remplir_face_normales();
};

// Inline functions

/*! @brief
 *
 */
inline IntTab& Domaine_VDF::Qdm()
{
  return Qdm_;
}

inline double Domaine_VDF::face_normales(int num_face,int k) const
{
  int ori = orientation(num_face);
  double surf=0;
  if (ori == k)
    surf=face_surfaces(num_face);
  return surf;
}

/*! @brief
 *
 */
inline const IntTab& Domaine_VDF::Qdm() const
{
  return Qdm_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_faces_X() const
{
  return nb_faces_X_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_faces_Y() const
{
  return nb_faces_Y_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_faces_Z() const
{
  return nb_faces_Z_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::Qdm(int num_arete,int i) const
{
  return Qdm_(num_arete,i);
}

/*! @brief inline double Domaine_VDF::porosite_face(int i) const {
 *
 *   return porosite_face_[i];
 *  }
 *
 */

/*! @brief
 *
 */
inline IntVect& Domaine_VDF::orientation()
{
  return orientation_;
}

/*! @brief
 *
 */
inline const IntVect& Domaine_VDF::orientation() const
{
  return orientation_;
}

/*! @brief inline DoubleVect& Domaine_VDF::porosite_face() {
 *
 *   return porosite_face_;
 *  }
 *
 */

/*! @brief inline const DoubleVect& Domaine_VDF::porosite_face() const {
 *
 *    return porosite_face_;
 *  }
 *
 */

/*! @brief inline DoubleVect& Domaine_VDF::porosite_elem() {
 *
 *   return porosite_elem_;
 *  }
 *
 */

/*! @brief inline const DoubleVect& Domaine_VDF::porosite_elem() const {
 *
 *   return porosite_elem_;
 *  }
 *
 */

/*! @brief inline double Domaine_VDF::porosite_elem(int i) const {
 *
 *   return porosite_elem_[i];
 *  }
 *
 */

/*! @brief
 *
 */
inline int Domaine_VDF::orientation(int i) const
{
  return orientation_[i];
}

// Compute function usable only in Cartesian coordinates:
// distance between the centers of 2 faces in direction k.

/*! @brief Returns the distance between two face centers in direction k (Cartesian only).
 *
 * @param fac1 Index of the first face.
 * @param fac2 Index of the second face.
 * @param k Direction index.
 * @return Distance xv(fac2,k) - xv(fac1,k).
 */
inline double Domaine_VDF::dist_face(int fac1, int fac2, int k) const
{
  // Note: this method is no longer called by dist_face in
  // Eval_Diff_VDF_Multi_inco_const.cpp and Eval_Diff_VDF_const.cpp for evaluator optimization
  return xv_(fac2,k) - xv_(fac1,k);
}

// Compute function usable only in cylindrical coordinates:
// distance between the centers of 2 faces in direction k.

/*! @brief Returns the distance between two face centers in direction k (cylindrical coordinates).
 *
 * @param fac1 Index of the first face.
 * @param fac2 Index of the second face.
 * @param k Direction index.
 * @return Arc-length distance (in direction k) between the two face centers.
 */
inline double Domaine_VDF::dist_face_axi(int fac1, int fac2, int k) const
{
  if (k != 1)
    return xv_(fac2,k) - xv_(fac1,k);
  else
    {
      double d_teta = xv_(fac2,1) - xv_(fac1,1);
      if (d_teta < 0)
        d_teta += 2.0*M_PI;
      return d_teta*xv_(fac1,0);
    }
}

// Compute function usable only in Cartesian coordinates:
// normal distance for an internal face.
// For an internal face the normal distance equals the distance between
// the centers of the two neighboring cells.

/*! @brief Returns the normal distance for an internal face (Cartesian coordinates).
 *
 * @param num_face Face index.
 * @return Distance between the centers of the two neighboring cells.
 */
inline double Domaine_VDF::dist_norm(int num_face) const
{
  int n1 = face_voisins_(num_face,0);
  int n2 = face_voisins_(num_face,1);
  int k = orientation_[num_face];
  return (xp_(n2,k) - xp_(n1,k));
}

// Compute function usable only in cylindrical coordinates:
// normal distance for an internal face.

/*! @brief Returns the normal distance for an internal face (cylindrical coordinates).
 *
 * @param num_face Face index.
 * @return Arc-length normal distance between the two neighboring cell centers.
 */
inline double Domaine_VDF::dist_norm_axi(int num_face) const
{
  int n1 = face_voisins_(num_face,0);
  int n2 = face_voisins_(num_face,1);
  int k = orientation_[num_face];
  double dist;
  if (k != 1)
    dist = xp_(n2,k) - xp_(n1,k);
  else
    {
      double d_teta = xp_(n2,1) - xp_(n1,1);
      if (d_teta < 0)
        d_teta += 2.0*M_PI;
      dist = d_teta*xp_(n1,0);
    }
  return dist;
}

// Compute function usable only in Cartesian coordinates:
// normal distance for a boundary face.
// For a boundary face the normal distance equals the distance between
// the center of the neighboring cell and the boundary.

/*! @brief Returns the normal distance for a boundary face (Cartesian coordinates).
 *
 * @param num_face Face index.
 * @return Distance between the neighboring cell center and the boundary face.
 */
inline double Domaine_VDF::dist_norm_bord(int num_face) const
{
  int n1 = face_voisins_(num_face,0);
  int n2 = face_voisins_(num_face,1);
  assert(num_face<nb_faces_bord() || n1==-1 || n2==-1); // Check that num_face is a real or virtual boundary face
  int k = orientation_[num_face];
  if (n1!=-1)
    return (xv_(num_face,k) - xp_(n1,k));
  else
    return (xp_(n2,k) - xv_(num_face,k));
}

// Compute function usable only in cylindrical coordinates:
// normal distance for a boundary face.

/*! @brief Returns the normal distance for a boundary face (cylindrical coordinates).
 *
 * @param num_face Face index.
 * @return Arc-length normal distance between the neighboring cell center and the boundary face.
 */
inline double Domaine_VDF::dist_norm_bord_axi(int num_face) const
{
  int n1 = face_voisins_(num_face,0);
  int n2 = face_voisins_(num_face,1);
  assert(num_face<nb_faces_bord() || n1==-1 || n2==-1); // Check that num_face is a real or virtual boundary face
  int k = orientation_[num_face];
  double dist;
  if (n1!=-1)
    if (k != 1)
      dist = xv_(num_face,k) - xp_(n1,k);
    else
      {
        double d_teta = xv_(num_face,1) - xp_(n1,1);
        if (d_teta < 0)
          d_teta += 2.0*M_PI;
        dist = d_teta*xp_(n1,0);
      }
  else if (k != 1)
    dist = xp_(n2,k) - xv_(num_face,k);
  else
    {
      double d_teta = xp_(n2,1) - xv_(num_face,1);
      if (d_teta < 0)
        d_teta += 2.0*M_PI;
      dist = d_teta*xp_(n2,0);
    }
  return dist;
}

// Compute function for the distance between the centers of 2 faces
// of the same orientation, usable in cylindrical and Cartesian coordinates.

/*! @brief Returns the distance between the centers of two same-orientation faces.
 *
 * @param n1 Index of the first face.
 * @param n2 Index of the second face.
 * @param k Direction index.
 * @return Distance between the two face centers (arc-length in cylindrical, Euclidean in Cartesian).
 */
inline double Domaine_VDF::distance_face(int n1, int n2, int k) const
{
  double dist,d_teta;
  assert ( (orientation_[n1]==orientation_[n2]) );
  if ( (axi!=1) || (k!=1) )
    dist = xv_(n2,k) - xv_(n1,k);
  else
    {
      d_teta = xv_(n2,1) - xv_(n1,1);
      if (d_teta < 0)
        d_teta += 2.0*M_PI;
      dist = d_teta*xv_(n1,0);
    }
  return dist;
}

// Compute function for the normal distance for any face,
// usable in cylindrical and Cartesian coordinates.

/*! @brief Returns the normal distance for any face (cylindrical or Cartesian coordinates).
 *
 * @param num_face Face index.
 * @return Normal distance between the two neighboring cell centers, or between a cell center and the boundary.
 */
inline double Domaine_VDF::distance_normale(int num_face) const
{

  double dist,d_teta;
  int n1 = face_voisins_(num_face,0);
  int n2 = face_voisins_(num_face,1);
  int k = orientation_[num_face];
  if ((n1!=-1) && (n2!=-1))
    {
      if ( (k!=1) || (axi!=1) )
        dist = xp_(n2,k) - xp_(n1,k);
      else
        {
          d_teta = xp_(n2,1) - xp_(n1,1);
          if (d_teta < 0)
            d_teta += 2.0*M_PI;
          dist = d_teta*xp_(n1,0);
        }
    }
  else if (n1!=-1)
    {
      if ( (k!=1) || (axi!=1) )
        dist = (xv_(num_face,k) - xp_(n1,k));
      else
        {
          d_teta = xv_(num_face,1) - xp_(n1,1);
          if (d_teta < 0)
            d_teta += 2.0*M_PI;
          dist = d_teta*xp_(n1,0);
        }
    }
  else
    {
      if ( (k!=1) || (axi!=1) )
        dist = (xp_(n2,k) - xv_(num_face,k));
      else
        {
          d_teta = xp_(n2,1) - xv_(num_face,1);
          if (d_teta < 0)
            d_teta += 2.0*M_PI;
          dist = d_teta*xp_(n2,0);
        }
    }
  return dist;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_aretes_joint() const
{
  return nb_aretes_joint_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_aretes_coin() const
{
  return nb_aretes_coin_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::premiere_arete_coin() const
{
  return nb_aretes_joint_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_aretes_bord() const
{
  return nb_aretes_bord_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::premiere_arete_bord() const
{
  return nb_aretes_joint_+ nb_aretes_coin_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_aretes_mixtes() const
{
  return nb_aretes_mixtes_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::premiere_arete_mixte() const
{
  return nb_aretes_ - nb_aretes_mixtes_ - nb_aretes_internes_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_aretes_internes() const
{
  return nb_aretes_internes_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::nb_aretes() const
{
  return nb_aretes_;
}

/*! @brief
 *
 */
inline int Domaine_VDF::premiere_arete_interne() const
{
  return nb_aretes_ - nb_aretes_internes_;
}

/*! @brief
 *
 */
inline double Domaine_VDF::h_x() const
{
  return h_x_;
}

/*! @brief
 *
 */
inline double Domaine_VDF::h_y() const
{
  return h_y_;
}

/*! @brief
 *
 */
inline double Domaine_VDF::h_z() const
{
  return h_z_;
}

/*! @brief
 *
 */
inline double Domaine_VDF::dist_elem(int n1, int n2, int k) const
{
  return xp_(n2,k)-xp_(n1,k);
}

/*! @brief
 *
 */
inline double Domaine_VDF::dist_elem_period(int n1, int n2, int k) const
{
  return xp_(n2,k) - xv_(elem_faces(n2,k),k)
         + xv_(elem_faces(n1,k+dimension),k) - xp_(n1,k);
}

/*! @brief
 *
 */
inline double Domaine_VDF::dim_elem(int n1, int k) const
{
  return xv_(elem_faces_(n1,k+dimension),k)-xv_(elem_faces_(n1,k),k);
}

/*! @brief
 *
 */
inline double Domaine_VDF::dim_face(int n1, int k) const
{
  int elem = std::max(face_voisins_(n1,0), face_voisins_(n1,1));
  return dim_elem(elem, k);
}

/*! @brief
 *
 */
inline double Domaine_VDF::delta_C(int elem) const
{
  double dist= 1;
  for (int i=0; i<dimension; i++)
    dist *= dim_elem(elem,i);
  return pow(dist,1./3.);
}

/*! @brief
 *
 */
inline int Domaine_VDF::amont_amont(int num_face, int i) const
{
  int k=orientation_[num_face];
  int num_elem = face_voisins_(num_face,i);
  int face = elem_faces_(num_elem,k+i*dimension);
  return face_voisins_(face,i);
}

/*! @brief
 *
 */
inline int Domaine_VDF::face_amont_princ(int num_face, int i) const
{
  int ori=orientation(num_face);
  int elem=face_voisins_(num_face,i);
  if(elem !=-1)
    elem=elem_faces_(elem,ori+i*dimension);
  return elem;
}

/*! @brief
 *
 */
inline int Domaine_VDF::face_amont_conj(int num_face, int k, int i) const
{
  int ori = orientation(num_face);
  int elem = face_voisins_(num_face,1);
  int face_conj=-2,face,elem_bis=-2;

  if(elem != -1)
    {
      face = elem_faces_(elem, k+i*dimension);
      elem_bis = face_voisins_(face,i);
      if (elem_bis != -1)
        face_conj = elem_faces_(elem_bis, ori);
      else
        face_conj = -1;
    }
  if ((elem==-1) || (elem_bis==-1))
    {
      elem = face_voisins_(num_face,0);
      if(elem != -1)
        {
          face = elem_faces_(elem, k+i*dimension);
          elem_bis = face_voisins_(face,i);
          if (elem_bis != -1)
            face_conj = elem_faces_(elem_bis, ori+dimension);
          else
            face_conj = -1;
        }
    }
  assert(face_conj!=-2);
  return face_conj;
}

/*! @brief Returns the neighboring face, accounting for the possibility that it may be a boundary face.
 *
 * @param num_face Face index.
 * @param k Direction index of the conjugate face.
 * @param i Side index (0 or 1).
 * @return Index of the neighboring face element.
 */
inline int Domaine_VDF::face_bord_amont(int num_face, int k, int i) const
{
  int ori = orientation(num_face);
  int elem = face_voisins_(num_face,1);
  if(elem != -1)
    {
      int face = elem_faces_(elem, k+i*dimension);
      int elem_bis = face_voisins_(face,i);
      if (elem_bis != -1)
        elem = elem_faces_(elem_bis, ori);
      else
        elem = -1;
    }
  if (elem == -1)
    {
      elem = face_voisins_(num_face,0);
      if(elem != -1)
        {
          int face = elem_faces_(elem, k+i*dimension);
          int elem_bis = face_voisins_(face,i);
          if (elem_bis != -1)
            elem = elem_faces_(elem_bis, ori+dimension);
          else
            elem = -1;
        }
    }
  return elem;
}

/*! @brief
 *
 */
inline int Domaine_VDF::elem_voisin(int elem, int face , int indic) const
{
  int ori = orientation_(face);
  return face_voisins_(elem_faces_(elem,ori+indic*dimension),indic);
}

/*! @brief Returns the distance between the center of a face and the center of face_voisins(face,0) (Cartesian coordinates only).
 *
 * @param num_face Face index.
 * @param n0 Index of the element on side 0.
 * @return Distance from the face center to the element center.
 */
inline double Domaine_VDF::dist_face_elem0(int num_face,int n0) const

{
  int ori = orientation_[num_face];
  return xv_(num_face,ori) - xp_(n0,ori);
}

/*! @brief Returns the distance between the center of a face and the center of face_voisins(face,1) (Cartesian coordinates only).
 *
 * @param num_face Face index.
 * @param n1 Index of the element on side 1.
 * @return Distance from the element center to the face center.
 */
inline double Domaine_VDF::dist_face_elem1(int num_face,int n1) const
{
  int ori = orientation_[num_face];
  return xp_(n1,ori) - xv_(num_face,ori);
}

/*! @brief Returns the distance between the center of a face and the center of face_voisins(face,0) (cylindrical coordinates only).
 *
 * @param num_face Face index.
 * @param n0 Index of the element on side 0.
 * @return Arc-length or Euclidean distance from the face center to the element center.
 */
inline double Domaine_VDF::dist_face_elem0_axi(int num_face,int n0) const
{
  int ori = orientation_[num_face];
  double dist;
  if (ori!=1)
    dist = xv_(num_face,ori) - xp_(n0,ori);
  else
    {
      double d_teta = xv_(num_face,1) - xp_(n0,1);
      if (d_teta < 0)
        d_teta += 2.0*M_PI;
      dist = d_teta*xp_(n0,0);
    }
  return dist;
}

/*! @brief Returns the distance between the center of a face and the center of face_voisins(face,1) (cylindrical coordinates only).
 *
 * @param num_face Face index.
 * @param n1 Index of the element on side 1.
 * @return Arc-length or Euclidean distance from the element center to the face center.
 */
inline double Domaine_VDF::dist_face_elem1_axi(int num_face,int n1) const
{
  int ori = orientation_[num_face];
  double dist;
  if (ori!=1)
    dist = xp_(n1,ori) - xv_(num_face,ori);
  else
    {
      double d_teta = xp_(n1,1) - xv_(num_face,1);
      if (d_teta < 0)
        d_teta += 2.0*M_PI;
      dist = d_teta*xp_(n1,0);
    }
  return dist;
}

inline double Domaine_VDF::dist_norm_period(int num_face,double l) const
{
  int n1 = face_voisins_(num_face,0);
  int n2 = face_voisins_(num_face,1);
  int k = orientation_[num_face];
  return (xp_(n2,k) + l - xp_(n1,k));
}

inline double Domaine_VDF::dist_face_elem0_period(int num_face,int n0,double l) const
{
  int ori = orientation_[num_face];
  double dist = xv_(num_face,ori) - xp_(n0,ori);
  if (dist > 0)
    return dist;
  else
    return dist + l;
}

inline double Domaine_VDF::dist_face_elem1_period(int num_face,int n1,double l) const
{
  int ori = orientation_[num_face];
  double dist = xp_(n1,ori) - xv_(num_face,ori);
  if (dist > 0)
    return dist;
  else
    return dist + l;
}

inline double Domaine_VDF::dist_face_period(int fac1 , int fac2, int k) const
{
  const Domaine& le_domaine = domaine();
  const DoubleTab& coord_sommets = le_domaine.coord_sommets();
  double dist= std::fabs(coord_sommets(face_sommets(fac1,1),k)-xv_(fac1,k));
  dist += std::fabs(xv_(fac2,k) - coord_sommets(face_sommets(fac2,0),k));
  return dist;

}

#endif




