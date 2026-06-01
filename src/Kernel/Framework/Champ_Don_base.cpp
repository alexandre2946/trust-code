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

#include <Champ_Don_base.h>
#include <Check_espace_virtuel.h>

Implemente_base(Champ_Don_base,"Champ_Don_base",Champ_base);
// XD champ_don_base field_base champ_don_base INHERITS_BRACE Basic class for data fields (not calculated), p.e. physics
// XD_CONT properties.


Sortie& Champ_Don_base::printOn(Sortie& os) const { return Champ_base::printOn(os); }
Entree& Champ_Don_base::readOn(Entree& is) { return Champ_base::readOn(is); }

/*! @brief Sets the number of degrees of freedom per component
 *
 * @param (int nb_noeuds) the number of degrees of freedom per component
 */
int Champ_Don_base::fixer_nb_valeurs_nodales(int nb_noeuds)
{
  valeurs_.resize(nb_noeuds, nb_compo_);
  return nb_noeuds;
}

/*! @brief Causes an error! To be overridden by derived classes! Not a pure virtual for development convenience!
 *
 */
Champ_base& Champ_Don_base::affecter_(const Champ_base&)
{
  Cerr << "Champ_Don_base::affecter_ : " << que_suis_je() << "::affecter_ must be overloaded in derived classes" << finl;
  throw;
}

/*! @brief Causes an error! To be overridden by derived classes! Not a pure virtual for development convenience!
 *
 */
Champ_base& Champ_Don_base::affecter_compo(const Champ_base&, int)
{
  Cerr << "Champ_Don_base::affecter_compo must be overloaded in derived classes" << finl;
  throw;
}

/*! @brief Time update.
 *
 * @param (double) the update time
 */
void Champ_Don_base::mettre_a_jour(double t)
{
  changer_temps(t);
}

/*!
 * See comments in Probleme_base_interface_proto::resetTime_impl()
 * Here we update.
 */
void Champ_Don_base::resetTime(double time)
{
  mettre_a_jour(time);
}

/*! @brief DOES NOTHING.
 *
 * To be overridden in derived classes. Causes the field initialization if necessary
 *
 */
int Champ_Don_base::initialiser(const double un_temps)
{
  mettre_a_jour(un_temps);
  return 1;
}

/*! @brief DOES NOTHING.
 *
 * To be overridden in derived classes
 *
 */
int Champ_Don_base::reprendre(Entree& )
{
  return 1;
}

/*! @brief DOES NOTHING.
 *
 * To be overridden in derived classes
 *
 */
int Champ_Don_base::sauvegarder(Sortie& ) const
{
  return 1;
}

/*! @brief DOES NOTHING.
 *
 * EXIT! To be overridden in derived classes
 *
 */
int Champ_Don_base::imprime(Sortie& os, int ncomp) const
{
  Cerr << que_suis_je() << "::imprime not coded." << finl;
  Process::exit();
  return 1;
}

/*! @brief Sets the number of components and the number of nodal values.
 *
 * @param (int) the number of nodes per component of the field (the number of dof per component)
 * @param (int) the number of components of the field
 */
void Champ_Don_base::dimensionner(int nb_noeuds, int nb_compo)
{
  fixer_nb_comp(nb_compo);
  fixer_nb_valeurs_nodales(nb_noeuds);
}
