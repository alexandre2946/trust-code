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

#include <Pb_Hydraulique.h>
#include <Fluide_Ostwald.h>

Implemente_instanciable(Pb_Hydraulique,"Pb_Hydraulique",Pb_Fluide_base);
// XD pb_hydraulique Pb_base pb_hydraulique INHERITS_BRACE Resolution of the Navier-Stokes equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr correlations bloc_lecture correlations OPT List of correlations used in specific source terms (e.g.
// XD_CONT wall-friction adaptive law).
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard REQ Navier-Stokes equations.
// XD attr Energie_cinetique_turbulente Energie_cinetique_turbulente Energie_cinetique_turbulente OPT Turbulent kinetic
// XD_CONT Energy conservation equation for a turbulent mono-phase problem (available in TrioCFD).
// XD attr Echelle_temporelle_turbulente Echelle_temporelle_turbulente Echelle_temporelle_turbulente OPT Turbulent
// XD_CONT Dissipation time scale equation for a turbulent mono-phase problem (available in TrioCFD).
// XD attr Taux_dissipation_turbulent Taux_dissipation_turbulent Taux_dissipation_turbulent OPT Turbulent Dissipation
// XD_CONT frequency equation for a turbulent mono-phase problem (available in TrioCFD).

Sortie& Pb_Hydraulique::printOn(Sortie& os) const { return Pb_Fluide_base::printOn(os); }
Entree& Pb_Hydraulique::readOn(Entree& is) { return Pb_Fluide_base::readOn(is); }

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_std si i=0, sort (exit) sinon.
 *
 *     (version const)
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation d'hydraulique de type Navier_Stokes_std
 */
const Equation_base& Pb_Hydraulique::equation(int i) const
{
  if (i == 0) return eq_hydraulique;
  else if (i < 1 + eq_opt_.size() && i > 0) return eq_opt_[i - 1].valeur();
  else
    {
      Cerr << "Pb_Hydraulique::equation() : Wrong equation number" << i << "!" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_std si i=0, sort (exit) sinon.
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation d'hydraulique de type Navier_Stokes_std
 */
Equation_base& Pb_Hydraulique::equation(int i)
{
  if (i == 0) return eq_hydraulique;
  else if (i < 1 + eq_opt_.size() && i > 0) return eq_opt_[i - 1].valeur();
  else
    {
      Cerr << "Pb_Hydraulique::equation() : Wrong equation number" << i << "!" << finl;
      Process::exit();
    }
  return eq_hydraulique;
}



/*! @brief Associe le milieu au probleme.
 *
 * Le milieu doit etre de type fluide incompressible.
 *
 * @param (Milieu_base& mil) le milieu physique a associer au probleme
 * @throws le milieu n'est pas du type Fluide_base
 */
void Pb_Hydraulique::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_base,mil) && ref_cast(Fluide_base, mil).is_incompressible())
    eq_hydraulique.associer_milieu_base(mil);
  else if (sub_type(Fluide_Ostwald,mil))
    {
      Cerr << "Je suis non newto" << finl;
      eq_hydraulique.associer_milieu_base(mil);
    }
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je() << " ne peut etre associe a "
           << finl;
      Cerr << "un probleme de type Pb_Hydraulique " << finl;
      exit();
    }
}

