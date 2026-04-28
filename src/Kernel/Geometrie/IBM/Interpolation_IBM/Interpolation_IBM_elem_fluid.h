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

#ifndef Interpolation_IBM_elem_fluid_included
#define Interpolation_IBM_elem_fluid_included

#include <Interpolation_IBM_base.h>

class Interpolation_IBM_elem_fluid : public Interpolation_IBM_base
{
  Declare_instanciable( Interpolation_IBM_elem_fluid ) ;

public :
  void discretise(const Discretisation_base&, Domaine_dis_base&) override;
  void set_fields_from_prepro_to_interp(Prepro_IBM_base&) override;

protected :
  virtual void computeFluidElems(Domaine_dis_base&);
  void set_param(Param&) const override;

  OWN_PTR(Champ_Don_base) fluid_points_lu_;
  OWN_PTR(Champ_Don_base) fluid_points_;

  OWN_PTR(Champ_Don_base) fluid_elems_lu_;
  OWN_PTR(Champ_Don_base) fluid_elems_;

  bool fluid_points_from_prepro_ = false;
  bool fluid_elems_from_prepro_ = false;

  friend class Source_PDF_base;
  friend class Source_PDF_EF;
  friend class Source_PDF_VEF;
  friend class Source_PDF_VDF;
};

#endif /* Interpolation_IBM_elem_fluid_included */
