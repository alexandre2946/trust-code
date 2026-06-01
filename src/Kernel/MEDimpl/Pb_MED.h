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


#ifndef Pb_MED_included
#define Pb_MED_included

#include <Discretisation_base.h>
#include <Probleme_Couple.h>
#include <Probleme_base.h>
#include <Champs_Fonc.h>

/*! @brief class Pb_MED Class for re-reading MED files and post-processing them.
 */
class Pb_MED : public Probleme_base
{
  Declare_instanciable(Pb_MED);

public:

  int nombre_d_equations() const override;
  const Equation_base& equation(int) const override ;
  Equation_base& equation(int) override;
  int comprend_champ(const Motcle& ) const;

  inline const ArrOfDouble& temps_sauv() const   { return temps_sauv_ ;      }
  inline Champs_Fonc& get_champs_fonc_post()     { return champs_fonc_post;  }

  /////////////////////////////////////////////////////
  // Methods of the post-processable fields interface
  /////////////////////////////////////////////////////
  void creer_champ(const Motcle& motlu) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom,Option opt=NONE) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;

  void typer_lire_milieu(Entree& is) override { /* Do nothing */ }

protected :
  Champs_Fonc champs_fonc_post;

private:
  Nom nom_fic;
  Noms nomschampmed;
  Discretisation dis_bidon;
  ArrOfDouble temps_sauv_;

  // A class carrying post-processable fields normally has an attribute
  // champs_compris_ containing a reference to those fields.
  // In the specific case of Pb_MED, this attribute is not declared,
  // so that a very specific get_champ() method can be implemented for this
  // problem type, which will not use the Champs_compris methods.
  // The methods get_noms_champs_postraitables() and creer_champ() will also
  // not manipulate champs_compris_.

  //Champs_compris champs_compris_;
};

class Pbc_MED : public Probleme_Couple
{
  Declare_instanciable(Pbc_MED);
private:
  OWN_PTR(Schema_Temps_base) sch_;
};

#endif
