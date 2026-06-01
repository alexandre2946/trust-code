/****************************************************************************
* Copyright (c) 2024, CEA
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

#ifndef Marqueur_Lagrange_base_included
#define Marqueur_Lagrange_base_included

#include <Champs_compris_interface.h>
#include <Champ_Fonc_base.h>
#include <Champs_compris.h>

class Ensemble_Lagrange_base;
class Probleme_base;
class Discretisation_base;

/*! @brief Marqueur_Lagrange_base The Marqueur_Lagrange_base class is the base class for Lagrangian marker classes.
 *
 *      Currently only one instantiable derived class: Marqueur_FT
 *      A marker is intended to track the fluid motion by integrating the trajectory
 *      of a set of particles from the (interpolated) fluid velocity.
 *      - the set of tracked points is an attribute of Marqueur_FT as it is of type Maillage_FT_Disc
 *      - post-processing is performed on the number of particles per cell (densite_particules_)
 *        or on the point cloud
 *      - integration starts from a time t_debut_integr_ set by the user
 *          or equal to t_init by default
 *
 * @sa Abstract class., Abstract methods:, Ensemble_Lagrange_base&  ensemble_points(), void calculer_valeurs_champs()
 */
class Marqueur_Lagrange_base: public Champs_compris_interface, public Objet_U
{
  Declare_base_sans_constructeur(Marqueur_Lagrange_base);

public :

  Marqueur_Lagrange_base();
  virtual Ensemble_Lagrange_base& ensemble_points() = 0;
  virtual const Ensemble_Lagrange_base& ensemble_points() const = 0;
  virtual void calculer_valeurs_champs() = 0;
  const inline double& temps_debut_integration() const;
  virtual void discretiser(const Probleme_base& pb, const  Discretisation_base& dis);
  virtual void mettre_a_jour(double temps);

  // Methods of the post-processable fields interface
  /////////////////////////////////////////////////////
  void creer_champ(const Motcle& motlu) override { }
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  /////////////////////////////////////////////////////

protected :

  OWN_PTR(Champ_Fonc_base)  densite_particules_; // Expresses the number of particles per cell
  double t_debut_integr_;            // Start time for trajectory integration

private :

  Champs_compris champs_compris_;

};

const double& Marqueur_Lagrange_base::temps_debut_integration() const
{
  return t_debut_integr_;
}

#endif
