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

#ifndef Cond_lim_rayo_milieu_transp_included
#define Cond_lim_rayo_milieu_transp_included

#include <Cond_lim_base.h>
#include <TRUST_Ref.h>

class Modele_rayo_transp;

class Cond_lim_rayo_milieu_transp
{
public:
  virtual ~Cond_lim_rayo_milieu_transp() { }

  virtual void completer();
  void preparer_surface(const Frontiere_dis_base&, const Domaine_Cl_dis_base&);

  inline virtual double surface(int numfa) const { return surf_i_[numfa]; }
  inline virtual double teta_i(int numfa) const { return teta_i_[numfa]; }

protected:
  DoubleVect surf_i_, teta_i_;
  OBS_PTR(Modele_rayo_transp) le_modele_rayo_;

  void error_pb_name(const Nom& nom_class, const Nom& nom_pb, const Nom& nom_pb_ray);
  void error_non_rad_bc(const Nom& nom_class, const Nom& nom_pb, const Nom& nom_bord, const Nom& other_type);
};

#endif /* Cond_lim_rayo_milieu_transp_included */
