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

#ifndef Sous_domaine_VF_included
#define Sous_domaine_VF_included

#include <Sous_domaine_dis_base.h>
#include <Domaine_VF.h>
#include <TRUST_Ref.h>

//! This abstract class contains the geometrical subdomain information common to finite-volume methods (VDF and VEF for example)
/**
   We assume that each face has at most two neighbouring elements in the
   domain, which is always the case in a conforming mesh.

   The array les_faces contains all the faces of the Domaine_dis that
   belong to this subdomain, sorted as follows:
   * Faces internal to the subdomain
   * Faces internal to the domain, but whose only first neighbouring element
   belongs to the subdomain
   * Faces internal to the domain, but whose only second neighbouring element
   belongs to the subdomain
   * Faces that have only one element in the domain.
   The separations are indicated respectively by the variables
   premiere_face_bord_0, premiere_face_bord_1, and premiere_face_bord.

   les_faces : nb_dim=1
   dimension(0) = number of faces of the subdomain
   integer value = index in the arrays le_dom_VF->face_voisins_,
   le_dom_VF->face_sommets_,...

   volumes_entrelaces(int face) returns the interlaced volume restricted
   to the subdomain. The face number refers to the array
   les_faces. Only the interlaced volumes that differ from those of the domain
   are stored locally.
*/

class Sous_domaine_VF : public Sous_domaine_dis_base
{

  Declare_instanciable(Sous_domaine_VF);

public:

  // Accessor methods
  inline const IntTab& les_faces() const
  {
    return les_faces_;
  }
  inline IntTab& les_faces()
  {
    return les_faces_;
  }
  inline int premiere_face_bord_0() const
  {
    return premiere_face_bord_0_;
  }
  inline int premiere_face_bord_1() const
  {
    return premiere_face_bord_1_;
  }
  inline int premiere_face_bord() const
  {
    return premiere_face_bord_;
  }
  //! Returns the interlaced volume restricted to the sub-domain. face is the index in the les_faces_ array.
  inline double volumes_entrelaces(int) const;

  // Specific methods

  //! Generates les_faces by iterating over the faces of domaine_dis and identifying which neighbours belong to le_sous_domaine.
  void discretiser() override;

protected:
  OBS_PTR(Domaine_VF) le_dom_VF;
  IntTab les_faces_;
  int premiere_face_bord_0_ = -10;
  int premiere_face_bord_1_= -10;
  int premiere_face_bord_= -10;
  DoubleTab volumes_entrelaces_;
};

inline double Sous_domaine_VF::volumes_entrelaces(int face) const
{
  if (face<premiere_face_bord_0_ || face>=premiere_face_bord_)
    return le_dom_VF->volumes_entrelaces(les_faces_(face));
  else
    return volumes_entrelaces_(face-premiere_face_bord_0_);
}

#endif
