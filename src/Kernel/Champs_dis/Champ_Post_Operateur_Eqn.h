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

#ifndef Champ_Post_Operateur_Eqn_included
#define Champ_Post_Operateur_Eqn_included

#include <Champ_Generique_Operateur_base.h>
#include <TRUST_Ref.h>
#include <Operateur.h>

class Equation_base;

/*! @brief class Champ_Post_Operateur_Eqn OWN_PTR(Champ_base) intended to post-process the gradient of a generic field
 *
 *  The class carries a statistical operator "gradient"
 *
 *  Syntax to follow for the data file
 *
 *  chosen_name operateur_eqn
 *    {
 *      sources { refchamp { pb_champ nom_pb nom_inconnue } }
 *      numero_source/numero op 1
 *    }
 * The name of the unknown must be entered, not that of the equation
 * The number of source terms is that of the data file, the order must be respected
 * For operators: 0 diffusion, 1 convection
 *
 */

class Champ_Post_Operateur_Eqn : public Champ_Generique_Operateur_base
{

  Declare_instanciable(Champ_Post_Operateur_Eqn);

public:

  const Noms get_property(const Motcle& query) const override;
  Entity  get_localisation(const int index = -1) const override;
  const Champ_base&  get_champ(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const override;
  const Champ_base&  get_champ_compo_without_evaluation(OWN_PTR(Champ_base)& espace_stockage) const;

  const Operateur_base& Operateur() const override;
  Operateur_base& Operateur() override;
  void completer(const Postraitement_base& post) override;
  void nommer_source() override;
  //virtual Entree &   lire(const Motcle & motcle, Entree & is);
  void set_param(Param& param) const override;
  void verification_cas_compo() const;

protected:
  int numero_source_ = -1;
  int numero_op_ = -1;
  int numero_masse_ = -1;
  OBS_PTR(Equation_base) ref_eq_;
  bool sans_solveur_masse_ = false;
  Entity localisation_inco_ = Entity::NODE;
  int compo_ = -1;                            //To identify the component to retrieve
};

#endif
