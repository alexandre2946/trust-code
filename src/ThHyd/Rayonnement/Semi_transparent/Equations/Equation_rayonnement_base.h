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

#ifndef Equation_rayonnement_base_included
#define Equation_rayonnement_base_included

#include <Operateur_Diff.h>
#include <Matrice_Morse.h>
#include <TRUST_Ref.h>

class Pb_rayo_semi_transp;
class Motcle;
class Milieu_base;
class Fluide_base;

class Equation_rayonnement_base: public Equation_base
{
  Declare_base(Equation_rayonnement_base);
public:

  void set_param(Param& titi) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  bool initTimeStep(double dt) override;
  void discretiser() override;
  void completer() override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  const Discretisation_base& discretisation() const;

  bool solve();

  void associer_pb_base(const Probleme_base& pb) override;
  void associer_milieu_base(const Milieu_base&) override;
  Milieu_base& milieu() override;
  const Milieu_base& milieu() const override;

  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;

  inline const Pb_rayo_semi_transp& pb_rayo_semi_transp() const { return pb_rayo_semi_transp_.valeur(); }
  inline Pb_rayo_semi_transp& pb_rayo_semi_transp() { return pb_rayo_semi_transp_.valeur(); }
  inline const Champ_Inc_base& inconnue() const override { return irradiance_.valeur(); }
  inline Champ_Inc_base& inconnue() override { return irradiance_.valeur(); }
  inline void associer_fluide(const Fluide_base& un_fluide) { le_fluide_ = un_fluide; }
  inline Fluide_base& fluide() { return le_fluide_.valeur(); }
  inline const Fluide_base& fluide() const { return le_fluide_.valeur(); }
  inline int nombre_d_operateurs() const override { return 1; }

  // pas de flux calcule correctement par les operateurs...
  inline int impr(Sortie& os) const override { return 1; }

  void Mat_Morse_to_Mat_Bloc(Matrice& matrice_tmp);
  void dimensionner_Mat_Bloc_Morse_Sym(Matrice& matrice_tmp);

  virtual int nb_colonnes_tot()=0;
  virtual int nb_colonnes()=0;
  virtual void resoudre(double temps)=0;
  virtual void assembler_matrice()=0;
  virtual void modifier_matrice()=0;
  virtual void evaluer_cl_rayonnement(double temps)=0;

protected:
  OBS_PTR(Fluide_base) le_fluide_;
  OBS_PTR(Pb_rayo_semi_transp) pb_rayo_semi_transp_;
  OWN_PTR(Champ_Inc_base) irradiance_;

  Operateur_Diff terme_diffusif_;
  Matrice_Morse la_matrice_;
  SolveurSys solveur_;
};

#endif /* Equation_rayonnement_base_included */
