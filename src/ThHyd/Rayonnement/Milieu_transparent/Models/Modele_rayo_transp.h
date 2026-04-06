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

#ifndef Modele_rayo_transp_included
#define Modele_rayo_transp_included

#include <Face_rayo_transp.h>
#include <TRUST_Vector.h>

class Pb_Fluide_base;

class Modele_rayo_transp: public Objet_U
{
  Declare_instanciable(Modele_rayo_transp);
public:
  void associer_pb_fluide_rayo(const Pb_Fluide_base& );
  void mettre_a_jour(double temps);
  void completer();
  void calculer_temperatures();
  void calculer_flux_radiatifs();
  void imprimer_flux_radiatifs(Sortie&) const;
  int postraiter();
  double flux_radiatif(int num_face_global) const; // 0 < face < nb_faces_de_bord

  inline void associer_processeur_rayonnant(int proc) { processeur_rayonnant_ = proc; }
  inline Face_rayo_transp& face_rayonnante(int j) { return les_faces_rayonnantes_[j]; }
  inline const Face_rayo_transp& face_rayonnante(int j) const { return les_faces_rayonnantes_[j]; }
  inline const Nom& nom_pb_rayonnant() const { return nom_pb_rayonnant_; }
  inline double relaxation() const { return relaxation_; }
  inline int processeur_rayonnant() { return processeur_rayonnant_; }
  inline int processeur_rayonnant() const { return processeur_rayonnant_; }
  inline int nb_faces_rayonnantes() const { return nb_faces_rayonnantes_; }
  inline int nb_faces_totales() const { return nb_faces_totales_; }
  inline int ordre_matrice_fac_forme() const { return ordre_mat_forme_; }

private:
  void lire_fichiers(Nom& nom1, Nom& nom2);
  void lire_fichiers(Nom& nom1, Nom& nom2, Nom& nom3);
  void init_matrice_rayo();

  OBS_PTR(Pb_Fluide_base) pb_fluide_rayo_;
  VECT(Face_rayo_transp) les_faces_rayonnantes_;
  int nb_faces_rayonnantes_ = -123, nb_faces_totales_ = -123, ordre_mat_forme_ = -123;
  double temps_ = -123.; // on garde le temps pour les impressions
  mutable int deja_imprime_ = 0;

  DoubleTab les_facteurs_de_forme_, matrice_rayo_, les_flux_radiatifs_;
  mutable IntVect corres_;

  static constexpr double sigma_ = 5.67e-8;
  double relaxation_ = 1.;

  // Par defaut, on suppose qu'il faut inverser la matrice de rayonnement
  bool lire_matrice_inv_ = false, fic_mat_ray_inv_bin_ = false;
  Nom nom_fic_mat_ray_inv_, nom_pb_rayonnant_;
  int processeur_rayonnant_ = -123;
};

#endif /* Modele_rayo_transp_included */
