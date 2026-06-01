/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Champ_front_var.h>
#include <Domaine_VF.h>

Implemente_base(Champ_front_var,"Champ_front_var",Champ_front_base);


/*! @brief Prints the field name to an output stream
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Champ_front_var::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}


/*! @brief DOES NOTHING must be overridden in derived classes
 *
 *     Must in all cases set nb_comp
 *
 * @param (Entree& is)
 * @return (Entree&) the modified input stream
 */
Entree& Champ_front_var::readOn(Entree& s )
{
  fixer_nb_comp(1);
  return s ;
}


/*! @brief Initialization at the beginning of the calculation.
 *
 * Sizes the value array and creates its virtual space.
 *     Derived classes must imperatively call this method.
 *
 * @return (0 in case of error, 1 otherwise.)
 */
int Champ_front_var::initialiser(double temps, const Champ_Inc_base& inco)
{
  const Frontiere& frontiere = frontiere_dis().frontiere();
  const int n = les_valeurs->nb_cases();
  const int nbc = nb_comp();
  for(int i = 0; i < n; i++)
    {
      DoubleTab& tab = les_valeurs[i].valeurs();
      // B.M. the completer method sometimes sizes the array.
      // We have a bit of everything in input: sometimes empty array,
      // sometimes nb_dim==1, etc...
      // Be careful, existing values must be preserved if there are any!
      if (tab.nb_dim()!=2 || tab.dimension(1)!=nbc)
        {
          tab.resize(tab.dimension(0), nbc);
          frontiere.creer_tableau_faces(tab);
        }
      tab.echange_espace_virtuel();
    }
  return 1;
}
