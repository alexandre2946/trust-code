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

#include <Postraitement.h>
#include <Sondes_Int.h>

Implemente_instanciable(Sondes_Int, "Sondes_Int", LIST(Sonde_Int));

Sortie& Sondes_Int::printOn(Sortie& s ) const { return s ; }

/*! @brief Reads a list of probes from an input stream. Format: { [READ A PROBE AS MANY TIMES AS NEEDED] }
 *
 * @param s an input stream
 * @return the modified input stream
 * @throws opening brace expected
 */
Entree& Sondes_Int::readOn(Entree& s )
{
  assert(mon_post);

  Motcle motlu;
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");

  s >> motlu;

  if (motlu != accolade_ouverte)
    {
      Cerr << "Error while reading the probes in the postprocessing" << finl;
      Cerr << "We expected { to start to read the probes" << finl;
      exit();
    }
  s >> motlu;

  if (motlu == accolade_fermee)
    {
      Cerr << "Error while reading the probes in the postprocessing" << finl;
      Cerr << "You have not defined any probe" << finl;
      exit();
    }
  while (motlu != accolade_fermee)
    {
      Sonde_Int une_sonde(motlu);
      une_sonde.associer_post(mon_post.valeur());
      s >> une_sonde;
      add(une_sonde);
      s >> motlu;
    }
  Cerr << "End of reading probes " << finl;
  return s;
}

/*! @brief Performs post-processing on each probe in the list.
 *
 * @param temps the current time
 */
void Sondes_Int::postraiter(double temps)
{
  for (auto& itr : *this) itr.postraiter(temps);
}

/*! @brief Updates each probe in the list in time.
 *
 * @param temps the update time
 * @param tinit the initial time of the probes
 */
void Sondes_Int::mettre_a_jour(double temps, double tinit)
{
  for (auto& itr : *this) itr.mettre_a_jour(temps,tinit);
}

/*! @brief Associates a post-processing object with the probe list.
 *
 */
void Sondes_Int::associer_post(const Postraitement& post)
{
  mon_post = post;
}
