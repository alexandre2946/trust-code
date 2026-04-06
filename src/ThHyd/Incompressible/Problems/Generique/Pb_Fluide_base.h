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

#ifndef Pb_Fluide_base_included
#define Pb_Fluide_base_included

#include <Modele_rayo_transp.h>
#include <Probleme_base.h>

/*! @brief classe  Pb_Fluide_base Cette classe a pour but de disposer d une classe amont pour
 *
 *      la hierarchie des problemes portant une equation de Navier_Stokes
 *
 * @sa Probleme_base
 */
class Pb_Fluide_base : public Probleme_base
{
  Declare_base(Pb_Fluide_base);
public:
  int expression_predefini(const Motcle& motlu, Nom& expression) override;

  /* Transparent radiation model */
  void preparer_calcul() override;
  int postraiter(int force = 1) override;
  void validateTimeStep() override;

  Entree& lire_radiation_models(Entree& is, Motcle& mot) override final;
  void assoscier_rayo_model_CL() override final;
  inline bool has_mod_rayo_transp() const override final { return mod_rayo_transp_.non_nul(); }

  inline Modele_rayo_transp& get_mod_rayo_transp()
  {
    if(mod_rayo_transp_.est_nul())
      Process::exit("Pb_Fluide_base::get_mod_rayo_transp() -- No transparent radiation model is associated for your problem !!! ");

    return mod_rayo_transp_.valeur();
  }

  inline const Modele_rayo_transp& get_mod_rayo_transp() const
  {
    if(mod_rayo_transp_.est_nul())
      Process::exit("Pb_Fluide_base::get_mod_rayo_transp() -- No transparent radiation model is associated for your problem !!! ");

    return mod_rayo_transp_.valeur();
  }

protected:
  OWN_PTR(Modele_rayo_transp) mod_rayo_transp_;
};

#endif /* Pb_Fluide_base_included */
