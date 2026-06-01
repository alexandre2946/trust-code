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

#include <Operateur_Diff_base.h>
#include <Milieu_base.h>
#include <Champ_base.h>

Implemente_base(Operateur_Diff_base,"Operateur_Diff_base",Operateur_base);

Sortie& Operateur_Diff_base::printOn(Sortie& os) const
{
  return os;
}

Entree& Operateur_Diff_base::readOn(Entree& is)
{
  return is;
}

/*! @brief Associates the true diffusivity in m^2/s (in QC for example, the operator is applied to rho*u, and the dynamic viscosity is then associated to
 *
 *   the operator instead of the kinematic viscosity. In that case, it
 *   is necessary to associate the true diffusivity here for the computation
 *   of the stability time step.
 *
 */
void Operateur_Diff_base::associer_diffusivite_pour_pas_de_temps(
  const Champ_base& diffu)
{
  if (je_suis_maitre())
    {
      Cerr << "Operateur_Diff_base::associer_diffusivite_pour_pas_de_temps\n";
      Cerr << " field name : " << diffu.le_nom();
      Cerr << " unit : "        << diffu.unite() << finl;
    }
  diffusivite_pour_pas_de_temps_ = diffu;
}

/*! @brief Returns the field corresponding to the true diffusivity of the medium used for the time step computation.
 *
 * If the operator is applied to rho*v
 *   for example (QC case for instance), the associated diffusivity is the
 *   dynamic viscosity and not the kinematic viscosity. The kinematic
 *   viscosity must be used for the time step computation...
 *   The true diffusivity (in m^2/s) is determined as follows:
 *   * if diffusivite_pour_pas_de_temps_ has been initialised => it is used
 *      (mechanism set up for front-tracking, because the medium
 *       is invalid)
 *   * otherwise diffusivite() is used
 *
 */
const Champ_base& Operateur_Diff_base::diffusivite_pour_pas_de_temps() const
{
  assert( diffusivite_pour_pas_de_temps_ );
  return diffusivite_pour_pas_de_temps_.valeur();
}

void Operateur_Diff_base::associer_diffusivite_volumique(const Champ_base& champ)
{
  Cerr << que_suis_je() << " does not support volumic diffusivity (received field " << champ.le_nom() << ")." << finl;
  Process::exit();
}
