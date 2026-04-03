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

#ifndef Face_rayo_transp_included
#define Face_rayo_transp_included

#include <Ensemble_faces_rayo_transp.h>
#include <TRUST_Vector.h>

class Face_rayo_transp: public Objet_U
{
  Declare_instanciable(Face_rayo_transp);
public:

  int chercher_ensemble_faces(const Nom&) const;
  void ecrire_temperature_bord() const;

  double calculer_temperature();
  double imprimer_flux_radiatif(Sortie&, Sortie&, Sortie&) const;

  inline double T_face_rayo() const { return T_face_rayo_; }
  inline double emissivite() const { return emissivite_; }
  inline double flux_radiatif() const { return flux_radiatif_; }
  inline double surface_rayo() const { return surf_; }
  inline void mettre_a_jour_flux_radiatif(double J) { flux_radiatif_ = J; }
  inline const Nom& nom_bord_rayo() const { return nom_bord_rayo_; }
  inline const Nom& nom_bord_rayo_lu() const { return nom_bord_rayo_lu_; }
  inline const IntVect& num_faces(int j) const { return num_faces_; }
  inline int nb_ensembles_faces() const { return nb_ensembles_faces_; }

  inline const Ensemble_faces_rayo_transp& ensembles_faces_bord(int j) const { return les_ensembles_faces_bord_[j]; }
  inline Ensemble_faces_rayo_transp& ensembles_faces_bord(int j) { return les_ensembles_faces_bord_[j]; }

private:
  VECT(Ensemble_faces_rayo_transp) les_ensembles_faces_bord_;
  IntVect num_faces_;
  Nom nom_bord_rayo_, nom_bord_rayo_lu_;
  double T_face_rayo_ = -123., emissivite_ = -123., surf_ = -123., flux_radiatif_ = 0.;
  int nb_ensembles_faces_ = 1;
};

#endif /* Face_rayo_transp_included */
