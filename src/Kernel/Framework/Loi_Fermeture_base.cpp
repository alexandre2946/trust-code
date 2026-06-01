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


#include <Loi_Fermeture_base.h>
#include <Probleme_base.h>
#include <Param.h>
#include <Champs_compris.h>
#include <Schema_Temps_base.h>

Implemente_base_sans_constructeur(Loi_Fermeture_base, "Loi_Fermeture_base", Objet_U);
// XD loi_fermeture_base objet_u loi_fermeture_base BRACE Class for appends fermeture to problem

Loi_Fermeture_base::Loi_Fermeture_base()
{
  status_ = INITIAL;
}

/*! @brief This method is the first called by the problem to construct the object when it is associated with the problem.
 *
 *   We check that it is not already associated.
 *
 */
void Loi_Fermeture_base::associer_pb_base(const Probleme_base& pb)
{
  if (status_ == PB_ASSOCIE)
    {
      Cerr << "Error associating " << que_suis_je() << " " << le_nom();
      Cerr << " to problem " << pb.que_suis_je();
      Cerr  << ".\n This object is already associated to problem " << mon_probleme().le_nom() << finl;
      barrier();
      exit();
    }
  assert(status_ == INITIAL);
  mon_probleme_ = pb;
  status_ = PB_ASSOCIE;
}

/*! @brief This method is called by the problem after the discretization of equations and medium and before
 *
 *   the call to readOn() for reading parameters.
 *   In derived classes it must at least discretize the fields
 *   that will be required in the readOn() of equations or other
 *   closure laws.
 *
 */
void Loi_Fermeture_base::discretiser(const Discretisation_base&)
{
  assert(status_ == PB_ASSOCIE);
  status_ = DISCRETISE;
}

/*! @brief This method calls the set_param() method to initialize the parameters, then reads the parameters.
 *
 *   In the implementation of derived classes we can simply
 *   call the base class method and check the parameters.
 *
 */
Entree& Loi_Fermeture_base::readOn(Entree& is)
{
  assert(status_ == DISCRETISE);
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  status_ = READON_FAIT;
  return is;
}

/*! @brief For now, exit()
 *
 */
Sortie& Loi_Fermeture_base::printOn(Sortie& os) const
{
  Cerr << "Loi_Fermeture_base::printOn not coded" << finl;
  exit();
  return os;
}

/*! @brief This method is called by the readOn of the class.
 *
 * It must be overridden in derived classes to add to "param" the
 *   various parameters to be read from the dataset and call
 *   the ancestor method. In the base class: no parameters.
 *
 */
void Loi_Fermeture_base::set_param(Param& param) const
{
}

/*! @brief This method is called after reading all equations and closure laws (all fields and boundary conditions
 *
 *   of the problem are available)
 *
 */
void Loi_Fermeture_base::completer()
{
  assert(status_ == READON_FAIT);
  status_ = COMPLET;
}

/*! @brief Returns the problem (this method was created to avoid giving write access to the problem via the REF)
 *
 */
const Probleme_base& Loi_Fermeture_base::mon_probleme() const
{
  return mon_probleme_.valeur();
}

/*! @brief This method is called by the problem after preparer_calcul() of equations and medium.
 *
 * It must update all fields it manages based on the other
 *   fields of the problem.
 *
 */
void Loi_Fermeture_base::preparer_calcul()
{
  assert(status_ == COMPLET);
  const double temps = mon_probleme().schema_temps().temps_courant();
  mettre_a_jour(temps);
}

/*! @brief This method is called by the problem after mettre_a_jour() of equations and medium.
 *
 * It must update all fields it manages based on the other
 *   fields of the problem.
 *
 */
void Loi_Fermeture_base::mettre_a_jour(double temps)
{
  assert(status_ == COMPLET);
}

/*! @brief This method returns the field named "nom" if it is understood by the class, otherwise calls the get_champ method of the ancestor.
 *
 *   In the base class, the Champ_compris_erreur exception is raised.
 *
 */
const Champ_base& Loi_Fermeture_base::get_champ(const Motcle& nom) const
{
  return champs_compris_.get_champ(nom);
}
bool Loi_Fermeture_base::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  return champs_compris_.has_champ(nom, ref_champ);
}

bool Loi_Fermeture_base::has_champ(const Motcle& nom) const
{
  return champs_compris_.has_champ(nom);
}

void Loi_Fermeture_base::get_noms_champs_postraitables(Noms& noms, Option opt) const
{
  if (opt == DESCRIPTION)
    Cerr << que_suis_je() << ": " << champs_compris_.liste_noms_compris() << finl;
  else
    noms.add(champs_compris_.liste_noms_compris());
}
