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

#ifndef Periodique_included
#define Periodique_included

#include <Cond_lim_base.h>
#include <TRUSTTab.h>

/*! @brief class Periodique This class represents a periodic boundary condition.
 *
 *      Periodicity can be specified in X, Y or Z.
 *      All faces of the boundary associated with this condition must
 *      have the same orientation.
 *
 * @sa Cond_lim_base, All faces of the boundary associated with this condition must have the same orientation
 */
class Periodique: public Cond_lim_base
{
  Declare_instanciable(Periodique);
public:
  void mettre_a_jour(double temps) override { }
  int face_associee(int i) const { return face_front_associee_[i]; }
  inline const ArrOfInt& face_associee() const { return face_front_associee_; }
  double distance() const { return distance_; }
  inline const ArrOfDouble& direction_perio() const { return direction_perio_; }
  int direction_periodicite() const;
  inline int est_periodique_selon_un_axe() const { return direction_xyz_ >= 0; }

protected:
  // Array of size nb_faces() + nb_faces_virt()
  // face_front_associee_[i] is the index of the opposite face on this boundary (-1 if the opposite face does not exist for a virtual face)
  ArrOfInt face_front_associee_;
  ArrOfDouble direction_perio_;
  double distance_ = -500.;
  // -1, 0, 1, or 2
  int direction_xyz_ = -100;
  int compatible_avec_eqn(const Equation_base&) const override { return 1; }
  void completer() override;
};

#endif
