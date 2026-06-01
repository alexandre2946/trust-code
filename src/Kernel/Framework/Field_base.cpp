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

#include <Field_base.h>

Implemente_base_sans_constructeur(Field_base,"Field_base",Objet_U);
// XD field_base objet_u champ_base INHERITS_BRACE Basic class of fields.


Field_base::Field_base() : nb_compo_(0), nature_(scalaire) { }

Sortie& Field_base::printOn(Sortie& s ) const { return s; }
Entree& Field_base::readOn(Entree& is ) { return is; }

/*! @brief Returns the name of the field.
 *
 */
const Nom& Field_base::le_nom() const
{
  return nom_;
}

/*! @brief Gives a name to the field
 *
 * @param (Nom& name) the name to give to the field
 */
void Field_base::nommer(const Nom& name)
{
  nom_ = name ;
}

/*! @brief Sets the number of components of the field.
 *
 * The field is vectorial if it has the same dimension as the space.
 *
 * @param (int i) the number of components of the field
 */
void Field_base::fixer_nb_comp(int i)
{
  // Prohibition to change the field nature once the nb_valeurs_nodales is set.
  nb_compo_ = i;
  noms_compo_.dimensionner(i);

  if (i == dimension) fixer_nature_du_champ(vectoriel);
  else if (i > dimension) fixer_nature_du_champ(multi_scalaire);
}

/*! @brief Sets the names of the field components
 *
 * @param (Noms& noms) the array of names to give to the field components
 * @return (Noms&) the array of names of the field components
 */
const Noms& Field_base::fixer_noms_compo(const Noms& noms)
{
  // teo boutin: not sure this assert is in line with the use of the method.
  // Maybe a madman sets the first 2 fields here then uses fixer_nom_compo(2, xxx) for the third one
  // Remove the assert for now, but I think it is necessary. Though a lot of tests may fail
  // assert(noms.size() == nb_comp());
  return noms_compo_ = noms;
}

/*! @brief Returns the array of names of the field components
 *
 */
const Noms& Field_base::noms_compo() const
{
  return noms_compo_;
}

/*! @brief Sets the name of the i-th component of the field
 *
 * @param (int i) the index of the field component whose name we want to specify
 * @param (Nom& nom) the name to give to the i-th component of the field
 * @return (Nom&) the name of i-th component of the field
 * @throws invalid field component index
 */
const Nom& Field_base::fixer_nom_compo(int i, const Nom& nom)
{
  assert(i < nb_comp());
  return noms_compo_[i] = nom;
}

/*! @brief Returns the name of the i-th component of the field
 *
 * @param (int i) the index of the field component whose name we want to specify
 * @return (Nom& nom) the name of i-th component of the field
 * @throws invalid field component index
 */
const Nom& Field_base::nom_compo(int i) const
{
  // teo boutin: those are not intensive operation I think. Maybe check those bounds all the time and exit in case of mistakes.
  // assert should be kept for very intensive operations in big loops, otherwise users who are also developers will make mistakes and never know... because they only run big simulations in optim
  assert(i < nb_comp());

  // it is important to check this bound too...
  // Surprinsingly otherwise does not cause a segfault. only valgrind error...
  assert(i < noms_compo_.size());
  return noms_compo_[i];
}

/*! @brief Sets the name of a scalar field
 *
 * @param (Nom& nom) the name to give to the field (scalar)
 * @return (Nom&) the name of the scalar field
 * @throws the field is not scalar
 */
const Nom& Field_base::fixer_nom_compo(const Nom& nom)
{
  assert(nb_comp() == 1);
  nommer(nom);
  noms_compo_.dimensionner_force(1);
  return noms_compo_[0] = nom;
}

/*! @brief Returns the name of a scalar field
 *
 */
const Nom& Field_base::nom_compo() const
{
  assert(nb_comp() == 1);
  return le_nom();
}

/*! @brief Specifies the units of the field components.
 *
 * These units are specified through an array of Nom and can be different for each component of the field.
 *
 * @param (Noms& noms) the names of the units of the field components
 */
const Noms& Field_base::fixer_unites(const Noms& noms)
{
  unite_.dimensionner_force(nb_comp());
  fixer_nature_du_champ(multi_scalaire);
  return unite_ = noms;
}

/*! @brief Returns the units of the field components
 *
 * @return (Noms&) the names of the units of the field components
 */
const Noms& Field_base::unites() const
{
  return unite_;
}

/*! @brief Specifies the unit of the i-th component of the field Meaning: the index of the field component whose unit we want to specify
 *
 * @param (Nom& nom) the type of unit to specify
 * @return (Nom&) the type of unit of the i-th component of the field
 */
const Nom& Field_base::fixer_unite(int i, const Nom& nom)
{
  fixer_nature_du_champ(multi_scalaire);
  return unite_[i] = nom;
}

/*! @brief Returns the unit of the i-th component of the field
 *
 * @param (int i) the index of the field component whose unit we want to know
 * @return (Nom&) the unit of the i-th component of the field
 */
const Nom& Field_base::unite(int i) const
{
  return unite_[i];
}

/*! @brief Specifies the unit of a scalar field or whose all components have the same unit
 *
 * @param (Nom& nom) the unit to specify
 * @return (Nom&) the unit of the field
 */
const Nom& Field_base::fixer_unite(const Nom& nom)
{
  unite_.dimensionner_force(1);
  return unite_[0] = nom;
}

/*! @brief Returns the unit of a scalar field whose all components have the same unit
 *
 * @return (Nom&) the (common) unit of the field components
 */
const Nom& Field_base::unite() const
{
  assert(unite_.size() == 1);
  return unite_[0];
}

/*! @brief Sets the nature of a field: scalar, multiscalar, vectorial.
 *
 * The type (enum) Nature_du_champ is defined in Ch_base.h.
 *
 * @param (Nature_du_champ n) the nature to assign to the field
 */
Nature_du_champ Field_base::fixer_nature_du_champ(Nature_du_champ n)
{
  return nature_ = n;
}

/*! @brief Returns the order of the basis functions
 *
 */
int Field_base::order_field() const
{
  switch (nature_)
    {
    case scalaire:
    case vectoriel:
    case quadrature_scalaire:
    case quadrature_vectoriel:
    case multi_scalaire:
      return 0;
    case basis_function_order_1_scalar:
    case basis_function_order_1_vectorial:
      return 1;
    case basis_function_order_2_scalar:
    case basis_function_order_2_vectorial:
      return 2;
    }
  return 0;
}

int Field_base::nb_vect_comp() const
{
  switch (nature_)
    {
    case scalaire:
    case multi_scalaire:
      return nb_compo_;
    case quadrature_scalaire:
    case basis_function_order_1_scalar:
    case basis_function_order_2_scalar:
      return 1;
    case vectoriel:
    case quadrature_vectoriel:
    case basis_function_order_1_vectorial:
    case basis_function_order_2_vectorial:
      return Objet_U::dimension;
    }
  return 0;
}
