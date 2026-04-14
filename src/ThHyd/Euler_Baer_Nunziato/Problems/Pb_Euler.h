/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef Pb_Euler_included
#define Pb_Euler_included

#include <Correlation_base.h>
#include <Fraction_Euler.h>
#include <Pb_Fluide_base.h>
#include <Momentum_Euler.h>
#include <Density_Euler.h>
#include <Energy_Euler.h>
#include <TRUST_Deriv.h>
#include <TRUST_List.h>
#include <Interprete.h>
#include <Verif_Cl.h>

class Pb_Euler : public Pb_Fluide_base
{
  Declare_instanciable(Pb_Euler);
public:
  void associer_milieu_base(const Milieu_base& ) override;
  void typer_lire_milieu(Entree& is) override;
  Entree& lire_equations(Entree& is, Motcle& dernier_mot) override;
  const Equation_base& equation(int) const override ;
  Equation_base& equation(int) override;
  void preparer_calcul() override;
  void mettre_a_jour(double temps) override;

  inline int nb_phases() const { return noms_phases_.size(); }
  inline Momentum_Euler& equation_qdm() { return eq_qdm_; }
  inline const Momentum_Euler& equation_qdm() const { return eq_qdm_; }
  inline Density_Euler& equation_masse() { return eq_masse_; }
  inline const Density_Euler& equation_masse() const { return eq_masse_; }
  inline Energy_Euler& equation_energie() { return eq_energie_; }
  inline const Energy_Euler& equation_energie() const { return eq_energie_; }
  inline Fraction_Euler& equation_fraction() { return eq_fraction_; }
  inline const Fraction_Euler& equation_fraction() const { return eq_fraction_; }
  inline const Nom& nom_phase(int i) const { return noms_phases_[i]; }
  inline const Noms& noms_phases() const { return noms_phases_; }
  inline int verifier() override { return 1; }
  inline int nombre_d_equations() const override { return 4; }

protected:
  Noms noms_phases_;
  Density_Euler eq_masse_;
  Momentum_Euler eq_qdm_;
  Energy_Euler eq_energie_;
  Fraction_Euler eq_fraction_;
};

#endif /* Pb_Euler_included */
