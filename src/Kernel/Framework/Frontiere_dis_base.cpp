/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Frontiere_dis_base.h>
#include <Frontiere.h>
Implemente_base(Frontiere_dis_base,"Frontiere_dis_base",Objet_U);


/*! @brief Overrides Objet_U::printOn(Sortie&) DOES NOTHING
 *
 *     To override in derived classes.
 *     Prints the discretized boundary to an output stream
 *
 * @param (Sortie& os) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Frontiere_dis_base::printOn(Sortie& os ) const
{
  return os;
}


/*! @brief DOES NOTHING - to override in derived classes.
 *
 *     Reads a discretized boundary from an input stream
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 */
Entree& Frontiere_dis_base::readOn(Entree& is )
{
  return is;
}


/*! @brief Associates the geometric boundary object with the discretized boundary.
 *
 * @param (Frontiere& fr) the geometric boundary
 */
void Frontiere_dis_base::associer_frontiere(const Frontiere& fr)
{
  la_frontiere=fr;
}

/*! @brief Returns the associated geometric boundary.
 *
 * (const version)
 *
 * @return (Frontiere&) the associated geometric boundary
 */
const Frontiere& Frontiere_dis_base::frontiere() const
{
  return la_frontiere.valeur();
}

/*! @brief Returns the associated geometric boundary.
 *
 * @return (Frontiere&) the associated geometric boundary
 */
Frontiere& Frontiere_dis_base::frontiere()
{
  return la_frontiere.valeur();
}

/*! @brief Returns the name of the geometric boundary.
 *
 * @return (Nom&) the name of the geometric boundary
 */
const Nom& Frontiere_dis_base::le_nom() const
{
  return la_frontiere->le_nom();
}

void Frontiere_dis_base::associer_Domaine_dis(const Domaine_dis_base& z)
{
  le_dom_dis=z;
}

const Domaine_dis_base& Frontiere_dis_base::domaine_dis() const
{
  return le_dom_dis.valeur();
}
