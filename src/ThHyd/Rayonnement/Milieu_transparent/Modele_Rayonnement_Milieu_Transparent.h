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

#ifndef Modele_Rayonnement_Milieu_Transparent_included
#define Modele_Rayonnement_Milieu_Transparent_included

#include <Face_Rayonnante.h>
#include <Domaine_forward.h>
#include <Matrice_Morse.h>
#include <TRUST_Vector.h>
#include <Cond_lim_base.h>

class Schema_Temps_base;
class Discretisation_base;

constexpr double SIGMA_DEFAULT = 5.67e-8;

class Modele_Rayonnement_Milieu_Transparent: public Objet_U
{
  Declare_instanciable_sans_constructeur(Modele_Rayonnement_Milieu_Transparent);
public:

  Modele_Rayonnement_Milieu_Transparent(double Sigma = SIGMA_DEFAULT) { SIGMA_ = Sigma; }

  void lire_fichiers(Nom& nom1, Nom& nom2);
  void lire_fichiers(Nom& nom1, Nom& nom2, Nom& nom3);
  void discretiser(const Discretisation_base&, const Domaine&);
  void mettre_a_jour(double temps);
  void preparer_calcul();
  void calculer_temperatures();
  void calculer_flux_radiatifs();
  void imprimer_flux_radiatifs(Sortie&) const;
  double flux_radiatif(int num_face_global) const; // 0 < face < nb_faces_de_bord

  inline void associer_processeur_rayonnant(int proc) { processeur_rayonnant_ = proc; }
  inline Face_Rayonnante& face_rayonnante(int j) { return les_faces_rayonnantes_[j]; }
  inline const Face_Rayonnante& face_rayonnante(int j) const { return les_faces_rayonnantes_[j]; }
  inline const Nom& nom_pb_rayonnant() const { return nom_pb_rayonnant_; }
  inline Nom& nom_pb_rayonnant() { return nom_pb_rayonnant_; }
  inline double relaxation() const { return relaxation_; }
  inline int processeur_rayonnant() { return processeur_rayonnant_; }
  inline int processeur_rayonnant() const { return processeur_rayonnant_; }
  inline int nb_faces_rayonnantes() const { return nb_faces_rayonnantes_; }
  inline int nb_faces_totales() const { return nb_faces_totales_; }
  inline int ordre_matrice_fac_forme() const { return ordre_mat_forme_; }

private:
  VECT(Face_Rayonnante) les_faces_rayonnantes_;
  int nb_faces_rayonnantes_ = -123, nb_faces_totales_ = -123, ordre_mat_forme_ = -123;
  double temps_ = -123.; // on garde le temps pour les impressions
  mutable int deja_imprime_ = 0;

  DoubleTab les_facteurs_de_forme_, matrice_rayo_, les_flux_radiatifs_;
  mutable IntVect corres_;

  double SIGMA_ = SIGMA_DEFAULT, relaxation_ = 1.;
  int inversion_debut_ = -123;
  // Par defaut, on suppose qu'il faut inverser la matrice de rayonnement
  int lire_matrice_inv_ = 0, fic_mat_ray_inv_bin_ = -123;
  Nom nom_fic_mat_ray_inv_, nom_pb_rayonnant_ = "non_donne";
  int processeur_rayonnant_ = -123;
};

#endif /* Modele_Rayonnement_Milieu_Transparent_included */
