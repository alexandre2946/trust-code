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

#ifndef Echange_contact_PolyMAC_CDO_included
#define Echange_contact_PolyMAC_CDO_included

#include <Echange_externe_impose.h>
#include <TRUSTTabs_forward.h>
#include <TRUST_Ref.h>

class Domaine_PolyMAC_CDO;
class Front_VF;
#include <Domaine_forward.h>

//  In addition to the champ_front representing the wall temperature,
//  this class holds another champ_front with as many time values
//  that represents the temperature in the other problem.
class Echange_contact_PolyMAC_CDO  : public Echange_externe_impose
{
  Declare_instanciable(Echange_contact_PolyMAC_CDO);
public :
  void completer() override;
  int initialiser(double temps) override;
  void calculer_correspondance();
  void update_coeffs();
  void update_delta() const;
  void mettre_a_jour(double ) override;
  inline Champ_front_base& T_autre_pb() { return T_autre_pb_; }
  inline const Champ_front_base& T_autre_pb() const { return T_autre_pb_; }
  inline const Nom& nom_autre_pb() const { return nom_autre_pb_; }

  //item(i, j): index of the j-th item needed for face i of the boundary
  mutable IntTab item;

  //coeff(i, j): coefficient of the face, then coefficient of item(i, j - 1) (element, then other faces) in the flux formula at the face
  //delta(i, j, 0/1) -> same for the nonlinear correction of Le Potier
  mutable DoubleTab coeff, delta_int, delta;
  int monolithic = 0; //1 if solving thermal problem monolithically
protected :
  int stab_ = 0; //1 if using the Le Potier stabilization
  mutable int coeffs_a_jour_ = 0, delta_a_jour_ = 0; //flag: whether coefficients have been updated
  double h_paroi = -123.;
  OWN_PTR(Champ_front_base) T_autre_pb_;
  Nom nom_autre_pb_;
};

#endif
