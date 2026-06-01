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

#include <Champ_front_uniforme.h>

Implemente_instanciable(Champ_front_uniforme,"Champ_front_uniforme",Champ_front_base);
// XD champ_front_uniforme front_field_base champ_front_uniforme NO_BRACE Boundary field which is constant in space and
// XD_CONT stationary.
// XD attr val list val REQ Values of field components.


/*! @brief Print the field to an output stream.
 *
 * Prints the field size and the (constant) value on
 *     the boundary.
 *
 * @param (Sortie& os) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Champ_front_uniforme::printOn(Sortie& os) const
{
  const DoubleTab& tab=valeurs();
  os << tab.size() << " ";
  for(int i=0; i<tab.size(); i++)
    os << tab(0,i);
  return os;
}


/*! @brief Read the field from an input stream.
 *
 * Format:
 *       Champ_front_uniforme nb_compo vrel_1 ... [vrel_i]
 *
 * @param (Entree& is) an input stream
 * @return (Entree& is) the modified input stream
 */
Entree& Champ_front_uniforme::readOn(Entree& is)
{
  int dim;
  DoubleTab& tab=les_valeurs->valeurs();
  dim=lire_dimension(is,que_suis_je());
  tab.resize(1,dim);
  fixer_nb_comp(dim);
  for(int i=0; i<dim; i++)
    is >> tab(0,i);
  return is;
}



/*! @brief Returns the object upcast to Champ_front_base&
 *
 * @param (Champ_front_base& ch)
 * @return (Champ_front_base&) (*this) upcast to Champ_front_base&
 */
Champ_front_base& Champ_front_uniforme::affecter_(const Champ_front_base& ch)
{
  return *this;
}

/*! @brief Returns the vector of field values for the given face.
 *
 * @return the field values array
 */
void Champ_front_uniforme::valeurs_face(int face,DoubleVect& var) const
{
  var.resize(nb_compo_);
  int i;
  for (i=0 ; i<nb_compo_ ; i++)
    var(i) = valeurs()(0,i);
}

/*! @brief Returns the values without caring about time since the field is stationary.
 *
 */
DoubleTab& Champ_front_uniforme::valeurs_au_temps(double temps)
{
  return les_valeurs->valeurs();
}

/*! @brief Returns the values without caring about time since the field is stationary.
 *
 */
const DoubleTab& Champ_front_uniforme::valeurs_au_temps(double temps) const
{
  return les_valeurs->valeurs();
}

/*! @brief Advance in time: nothing to do for a stationary field!
 *
 */
int Champ_front_uniforme::avancer(double temps)
{
  return 1;
}

/*! @brief Step back in time: nothing to do for a stationary field!
 *
 */
int Champ_front_uniforme::reculer(double temps)
{
  return 1;
}

/*! @brief Nothing to do for a stationary field!
 *
 */
void Champ_front_uniforme::changer_temps_futur(double temps,int i)
{
}
