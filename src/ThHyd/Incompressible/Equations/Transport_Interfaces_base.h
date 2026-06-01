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

#ifndef Transport_Interfaces_base_included
#define Transport_Interfaces_base_included

#include <Equation_base.h>


/*! @brief Transport_Interfaces_base This class is the base class for interface transport equations.
 *
 *      Currently only one instantiable derived class: Transport_Interfaces_FT_Disc
 *
 *
 * @sa Transport_Interfaces_FT_Disc, Abstract methods:, void modifier_vpoint_pour_imposer_vit(...), void integrer_ensemble_lagrange(...)
 */

class Transport_Interfaces_base : public Equation_base
{
  Declare_base(Transport_Interfaces_base);

public:

  // The presence of a solid interface in a flow is accompanied by a source term
  // in the momentum equation that cannot be considered as a classical source term
  // but rather as a modification of vpoint to impose the interface velocity on the fluid.

  // Applies the modification of vpoint (momentum equation)
  // to impose the interface velocity on the fluid.
  virtual void modifier_vpoint_pour_imposer_vit(const DoubleTab& inco_val,DoubleTab& vpoint0,DoubleTab& vpoint,
                                                const DoubleTab& rho_faces,DoubleTab& source_val,
                                                const double temps, const double dt,
                                                const int is_explicite = 1, const double eta = 1.) = 0;
  // Returns the tag of the interface mesh.
  virtual int get_mesh_tag() const = 0;
  // Returns the indicator field.
  // TODO : Why virtual here? Why not getting the daughter version from Transport_Interfaces_FT_Disc directly here?
  virtual const Champ_base& get_indicatrice() =0;
  virtual void update_indicatrice() =0;
  virtual void check_indicatrice_is_up_to_date() =0;


  // Performs trajectory integration of point particles
  // intended to mark the fluid.
  // xn+1 = xn + v_interpolated*dt
  virtual void integrer_ensemble_lagrange(const double temps) = 0;

protected:

  // In the future, other elements of Transport_Interfaces_FT_Disc may be factored out here.
  // e.g.: OBS_PTR(Probleme_base) probleme_base_;
  // This refactoring can be done once the recent developments
  // purely related to Front-Tracking are integrated.

private:

};

#endif
