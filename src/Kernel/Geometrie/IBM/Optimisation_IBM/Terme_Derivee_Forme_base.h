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

#ifndef Terme_Derivee_Forme_base_included
#define Terme_Derivee_Forme_base_included

#include <TRUST_Ref.h>
#include <TRUST_Deriv.h>
#include <TRUSTTabs_forward.h>
#include <Champ_Don_base.h>
#include <Champ_Fonc_base.h>
#include <Source_base.h>

/*! @brief Classe Terme_Derivee_Forme_base Cette classe represente un terme source de l'equation de projection en optimisation de forme
 *
 */
class Terme_Derivee_Forme_base : public Source_base
{

  Declare_base(Terme_Derivee_Forme_base);

public :
  DoubleTab& calculer(DoubleTab& ) const override;
  void mettre_a_jour(double ) override;
  void modify_name_file(Nom& ) const;

  void set_source_derivee_forme(DoubleTab&) const;
  const DoubleTab& get_source_derivee_forme() const { return source_derivee_forme->valeurs(); };

  // Methodes de l interface des champs postraitables
  void creer_champ(const Motcle& motlu) override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;
  const Champ_base& get_champ(const Motcle&) const override;

protected:
  OWN_PTR(Champ_Don_base) source_derivee_forme;
  mutable OWN_PTR(Champ_Fonc_base)  champ_derivee_forme_; //!< Champ pour postraitement
};

#endif
