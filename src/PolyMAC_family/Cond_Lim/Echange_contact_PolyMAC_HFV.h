/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Echange_contact_PolyMAC_HFV_included
#define Echange_contact_PolyMAC_HFV_included

#include <Op_Diff_PolyMAC_HFV_Elem.h>
#include <Echange_externe_impose.h>
#include <TRUSTTabs_forward.h>
#include <MD_Vector_tools.h>
#include <TRUST_Ref.h>

class Domaine_PolyMAC_HFV;
class Front_VF;
#include <Domaine_forward.h>

/*! @brief class : Echange_contact_PolyMAC_HFV
 *
 *  @brief In addition to the champ_front representing the wall temperature,
 *   this class holds another champ_front with as many time values
 *   that represents the temperature in the other problem.
 *
 */
class Echange_contact_PolyMAC_HFV  : public Echange_externe_impose
{
  Declare_instanciable(Echange_contact_PolyMAC_HFV);
public :
  void init_op() const;
  void mettre_a_jour(double temps) override { } //not used
  void verifie_ch_init_nb_comp() const override { } //no constraint on the number of components on each side

  mutable OBS_PTR(Front_VF) fvf, o_fvf; //boundary in the other problem
  mutable int i_fvf = -1 , i_o_fvf = -1;  //boundary indices on each side
  mutable OBS_PTR(Op_Diff_PolyMAC_HFV_Elem) diff, o_diff; //diffusion operators on each side
  mutable int o_idx = -1; //index of the other operator in the op_ext array of Op_Diff_PolyMAC_HFV_Elem (to be filled by that operator)

  /* faces on the other side of the boundary */
  void init_f_dist() const; //initialization of f_dist
  mutable IntTab f_dist;     //face on the other side of each boundary face
  mutable int f_dist_init_ = 0;

  double invh_paroi = 1e30; //thermal resistance (1 / h) of the wall

protected :
  Nom nom_autre_pb_, nom_bord_, nom_champ_; //name of the remote problem, the boundary, and the field
};

#endif /* Echange_contact_PolyMAC_HFV_included */
