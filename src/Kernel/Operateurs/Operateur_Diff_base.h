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

#ifndef Operateur_Diff_base_included
#define Operateur_Diff_base_included

#include <Support_Champ_Masse_Volumique.h>
#include <Correlation_base.h>
#include <Operateur_base.h>
#include <TRUST_Ref.h>

class Champ_base;

/*! @brief Operateur_Diff_base This class is the base of the hierarchy of operators representing
 *
 *     a diffusion term in an equation. The choice of term depends
 *     on the laminar or turbulent modelling of the flow, on the
 *     discretisation and on the type of the diffusivity field. These variants
 *     give rise to derived classes of Operateur_Diff_base.
 *
 * @sa Operateur_base Operateur_Diff, Abstract class, Abstract method, void associer_diffusivite(const Champ_Don_base& ), const Champ_Don_base& diffusivite() const
 */
class Operateur_Diff_base  : public Operateur_base,
  public Support_Champ_Masse_Volumique
{
  Declare_base(Operateur_Diff_base);
public:
  virtual void associer_diffusivite(const Champ_base&) = 0;
  virtual void associer_diffusivite_pour_pas_de_temps(const Champ_base&);
  virtual void associer_diffusivite_volumique(const Champ_base&);
  virtual const Champ_base& diffusivite() const=0;
  inline virtual void calculer_borne_locale(DoubleVect& ,double,double ) const {};

  //list of Op_Diff from problems solved simultaneously (monolithic thermal)
  mutable std::vector<const Operateur_Diff_base *> op_ext;
  virtual void init_op_ext() const { op_ext = { this }; }    //filling of op_ext (cannot be done in completer(), too early)

  virtual bool is_turb() const { return false; }
  virtual const Correlation_base* correlation_viscosite_turbulente() const { return nullptr; }
  virtual void calculer_von_mises(const DoubleTab& deplacement, DoubleTab& deformation, DoubleTab& contraintes, DoubleTab& von_mises) const { throw; }

protected:
  virtual const Champ_base& diffusivite_pour_pas_de_temps() const;
  OBS_PTR(Champ_base) diffusivite_pour_pas_de_temps_;
};

#endif
