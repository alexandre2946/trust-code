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

#ifndef Cond_lim_rayo_semi_transp_included
#define Cond_lim_rayo_semi_transp_included

#include <Neumann_sortie_libre.h>
#include <Cond_lim_base.h>
#include <TRUST_Ref.h>

class Pb_rayo_semi_transp;

/*! @brief classe Cond_lim_rayo_semi_transp
 *
 * @sa Ce n'est pas une classe de l'arbre TRUST a elle seule., Cette classe est faite etre une classe mere d'une classe, qui heritera par ailleurs d'Objet_U
 */
class Cond_lim_rayo_semi_transp
{
public:
  virtual ~Cond_lim_rayo_semi_transp() { }

  virtual void associer_pb_rayo_semi_transp(const Pb_rayo_semi_transp&);
  inline const Pb_rayo_semi_transp& pb_rayo_semi_transp() const
  {
    assert(pb_rayo_semi_transp_.non_nul());
    return pb_rayo_semi_transp_.valeur();
  }
  inline Pb_rayo_semi_transp& pb_rayo_semi_transp()
  {
    assert(pb_rayo_semi_transp_.non_nul());
    return pb_rayo_semi_transp_.valeur();
  }

  inline Champ_front_base& emissivite() { return emissivite_; }
  inline const Champ_front_base& emissivite() const { return emissivite_; }
  inline double& A() { return A_; }
  inline const double& A() const { return A_; }

  virtual void recherche_emissivite_et_A();
  virtual const Cond_lim_base& la_cl() const =0;
  virtual void completer_Cl_opposee_si_contact() { }

protected:
  OBS_PTR(Pb_rayo_semi_transp) pb_rayo_semi_transp_;
  OWN_PTR(Champ_front_base) emissivite_;
  double A_ = -123.;
};

#endif /* Cond_lim_rayo_semi_transp_included */
