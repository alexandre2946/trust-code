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

#include <Champ_Fonc_Fonction_txyz_Morceaux.h>
#include <Champ_Fonc_Tabule.h>

Implemente_instanciable(Champ_Fonc_Fonction_txyz_Morceaux,"Champ_Fonc_Fonction_txyz_Morceaux",TRUSTChamp_Morceaux_generique<Champ_Morceaux_Type::FONC_TXYZ>);
// XD champ_fonc_fonction_txyz_morceaux champ_don_base champ_fonc_fonction_txyz_morceaux NO_BRACE Field defined by
// XD_CONT analytical functions in each sub-domaine. On each zone, the value is defined as a function of x,y,z,t and of
// XD_CONT scalar value taken from a parameter field. This values is associated to the variable 'val' in the expression.
// XD attr problem_name ref_Pb_base problem_name REQ Name of the problem.
// XD attr inco chaine inco REQ Name of the field (for example: temperature).
// XD attr nb_comp entier nb_comp REQ Number of field components.
// XD attr data bloc_lecture data REQ { Defaut val_def sous_domaine_1 val_1 ... sous_domaine_i val_i } By default, the
// XD_CONT value val_def is assigned to the field. It takes the sous_domaine_i identifier Sous_Domaine (sub_area) type
// XD_CONT object function, val_i. Sous_Domaine (sub_area) type objects must have been previously defined if the
// XD_CONT operator wishes to use a champ_fonc_fonction_txyz_morceaux type object.

Sortie& Champ_Fonc_Fonction_txyz_Morceaux::printOn(Sortie& os) const { return os << valeurs(); }

/*! @brief Reads the values of the piecewise uniform field from an input stream.
 *
 * Reads the domain name (nom_domaine), the number of field components (nb_comp), the default value
 *     of the field and the values on the subdomains.
 *     Format:
 *      Champ_Fonc_Fonction_txyz_Morceaux pb champ nb_comp { Defaut val_def sous_domaine_1 val_1 ... sous_domaine_i val_i }
 *
 */
Entree& Champ_Fonc_Fonction_txyz_Morceaux::readOn(Entree& is)
{
  int dim;
  Nom nom;
  is >> nom;
  interprete_get_domaine(nom);

  is >> nom_champ_parametre_;
  bool isNum = Champ_Fonc_Tabule::Check_if_int(nom_champ_parametre_);
  if (isNum)
    {
      Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << finl;
      Cerr << "Error in call to " << que_suis_je() << ":" << finl;
      Cerr << "The syntax has changed in version 1.8.2." << finl;
      Cerr << "You should now pass the dimension/number of components AFTER the field/parameter name." << finl;
      Cerr << "Please update your dataset or contact TRUST support team." << finl;
      Process::exit();
    }
  is >> dim;

  creer_tabs(dim);
  is >> nom;

  return complete_readOn(dim,que_suis_je(),is,nom);
}
