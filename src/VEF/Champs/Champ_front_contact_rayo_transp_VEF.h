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

#ifndef Champ_front_contact_rayo_transp_VEF_included
#define Champ_front_contact_rayo_transp_VEF_included

#include <Modele_rayo_transp.h>
#include <Champ_front_contact_VEF.h>

class Champ_front_contact_rayo_transp_VEF: public Champ_front_contact_VEF
{
  Declare_instanciable(Champ_front_contact_rayo_transp_VEF);
public:

  Champ_front_base& affecter_(const Champ_front_base& ch) override;
  int initialiser(double temps, const Champ_Inc_base& inco) override;
  void calculer_temperature_bord(double temps);
  void mettre_a_jour(double temps) override;
  void mettre_a_jour_flux_radiatif(); // Le fait calculer par le modele et le stocke
  void calcul_grads_locaux(double temps) override;
  void modifie_gradients_pour_rayonnement(DoubleVect& gradient_num_transf, DoubleVect& gradient_num_transf_autre_pb);
  void calculer_coeffs_echange(double temps) override;

  inline Champ_Inc_base& inconnue1() { return l_inconnue1.valeur(); }
  inline const Champ_Inc_base& inconnue1() const { return l_inconnue1.valeur(); }

  inline Champ_Inc_base& inconnue2() { return l_inconnue2.valeur(); }
  inline const Champ_Inc_base& inconnue2() const { return l_inconnue2.valeur(); }

  inline Nom& nom_prob1() { return nom_pb1; }
  inline const Nom& nom_prob1() const { return nom_pb1; }
  inline Nom& nom_prob2() { return nom_pb2; }
  inline const Nom& nom_prob2() const { return nom_pb2; }

  inline void associer_modele_rayo(Modele_rayo_transp& mod) { le_modele_rayo_ = mod; }
  inline Modele_rayo_transp& modele_rayo() { return le_modele_rayo_; }

protected:
  DoubleVect flux_radiatif_;
  OBS_PTR(Modele_rayo_transp) le_modele_rayo_;
};

#endif /* Champ_front_contact_rayo_transp_VEF_included */
