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

#include <Coloc_discretisation.h>
#include <Domaine_Coloc.h>
#include <Equation_base.h>
#include <Probleme_base.h>

Implemente_instanciable(Coloc_discretisation, "Coloc", Discret_Thyd);

// XD coloc discretisation_base coloc INHERITS_BRACE Co-localised cell-center discretization

Entree& Coloc_discretisation::readOn(Entree& s) { return s;}

Sortie& Coloc_discretisation::printOn(Sortie& s) const { return s;}

void Coloc_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, int nb_pas_dt,
                                             double temps, OWN_PTR(Champ_Inc_base)& champ, const Nom& sous_type) const
{
  Motcles motcles(6);
  motcles[0] = "vitesse";     // Standard choice for velocity
  motcles[1] = "pression";    // Standard choice for pressure
  motcles[2] = "temperature"; // Standard choice for temperature
  motcles[3] = "divergence_vitesse"; // Field type obtained by computing div v
  motcles[4] = "gradient_pression";  // Field type obtained by computing grad P
  motcles[5] = "champ_elem";    // Create a field at elements (type P0)

  Nom type_champ_vec("Champ_Vect_Elem_Coloc");
  Nom type_champ_elem("Champ_Elem_Coloc");
  Nom type;
  int default_nb_comp = 0; // Default value for the number of components
  int rang = motcles.search(directive);

  switch(rang)
    {
    case 0:
    case 4:
      type = type_champ_vec;
      default_nb_comp = 3;
      break;
    case 1:
    case 2:
    case 3:
    case 5:
      type = type_champ_elem;
      default_nb_comp = 1;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "Coloc_discretisation : " << motcles;

  if (sous_type != NOM_VIDE)
    rang = verifie_sous_type(type, sous_type, directive);

  // If the directive was not understood (or if it is a description request),
  // call the parent class:
  if (rang < 0)
    {
      Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, nb_pas_dt, temps, champ);
      return;
    }

  int nb_ddl = 0;
  if (type.debute_par(type_champ_elem) || type.debute_par(type_champ_vec))
    nb_ddl = z.nb_elem() ;
  else
    assert(0);

  if (nb_comp < 0)
    nb_comp = default_nb_comp;
  assert(nb_comp > 0);
  creer_champ(champ, z, type, noms[0], unites[0], nb_comp, nb_ddl, nb_pas_dt, temps, directive, que_suis_je());

  if (nature == multi_scalaire)
    {
      champ->fixer_nature_du_champ(nature);
      champ->fixer_unites(unites);
      champ->fixer_noms_compo(noms);
    }
}

void Coloc_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                             OWN_PTR(Champ_Fonc_base)& champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

void Coloc_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                             OWN_PTR(Champ_Don_base)& champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

Nom Coloc_discretisation::get_name_of_type_for(const Nom& class_operateur, const Nom& type_operateur, const Equation_base& eqn, const OBS_PTR(Champ_base) &champ_sup) const
{
  if (eqn.probleme().que_suis_je() != "Pb_Euler")
    Process::exit("\nError with the chosen discretization !!!! \nThe colocalised discretization currently works only for the Euler problem !! Please select another one ... \n");

  Nom type;
  Nom type_ch = eqn.inconnue().que_suis_je();
  if (type_ch.debute_par("Champ_Elem"))  type_ch = "_Elem";
  else if (type_ch.debute_par("Champ_Vect")) type_ch = "_Vect";

  if (class_operateur == "Source")
    type = type_operateur + type_ch + "_" + que_suis_je();
  else if (class_operateur == "Solveur_Masse")
    type = Nom("Masse_") + que_suis_je() + type_ch;
  else if (class_operateur == "Operateur_Grad")
    type = Nom("Op_Grad_") + type_operateur + "_" + que_suis_je() + type_ch;
  else if (class_operateur == "Operateur_Diff")
    type = Nom("Op_Diff") + (type_operateur != "" ? "_" : "") + type_operateur + "_" + que_suis_je() + type_ch;
  else if (class_operateur == "Operateur_Conv")
    type = Nom("Op_Conv_") + type_operateur + "_" + que_suis_je() + type_ch;
  else if (class_operateur == "Operateur_NConserv")
    type = Nom("Op_NConserv_") + type_operateur + "_" + que_suis_je() + type_ch;
  else
    return Discret_Thyd::get_name_of_type_for(class_operateur, type_operateur, eqn);

  return type;
}

void Coloc_discretisation::discretiser_champ_fonc_don(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp,
                                                      double temps, Objet_U& champ) const
{
// Two pointers for easy access to champ_don or champ_fonc, depending on the type of the champ object.
  OWN_PTR(Champ_Fonc_base) *champ_fonc = dynamic_cast<OWN_PTR(Champ_Fonc_base)*>(&champ);
  OWN_PTR(Champ_Don_base) *champ_don = dynamic_cast<OWN_PTR(Champ_Don_base)*>(&champ);

  Motcles motcles(6);
  motcles[0] = "pression";    // Standard choice for pressure
  motcles[1] = "temperature"; // Standard choice for temperature
  motcles[2] = "divergence_vitesse"; // Field type obtained by computing div v
  motcles[3] = "champ_elem";  // Create a field at elements (type P0)
  motcles[4] = "vitesse";     // Standard choice for velocity
  motcles[5] = "gradient_pression";  // Field type obtained by computing grad P

// The velocity field type depends on the element type:
  int default_nb_comp = 0,
      rang = motcles.search(directive);

  Nom type_champ_elem("Champ_Fonc_Elem_Coloc"),
      type_champ_vect("Champ_Fonc_Vect_Coloc"),
      type;

  switch(rang)
    {
    case 0:
    case 1:
    case 2:
    case 3:
      type = type_champ_elem;
      default_nb_comp = 1;
      break;
    case 4:
    case 5:
      type = type_champ_vect;
      default_nb_comp = 3;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "Coloc_discretisation : " << motcles;

// If the directive was not understood (or if it is a description request),
// call the parent class:
  if (rang < 0)
    {
      if (champ_fonc)
        Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, temps, *champ_fonc);
      else
        Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, temps, *champ_don);
      return;
    }

// Compute the number of degrees of freedom
  int nb_ddl = 0;
  if ((type == "Champ_Fonc_Elem_Coloc") || (type == "Champ_Fonc_Vect_Coloc"))    // to be replaced by debute_par("Champ_Elem")
    nb_ddl = z.nb_elem();
  else
    assert(0);

// If it is a multi-scalar field, handle accordingly
  if (nb_comp < 0)
    nb_comp = default_nb_comp;
  if (champ_fonc)
    creer_champ(*champ_fonc, z, type, noms[0], unites[0], nb_comp, nb_ddl, temps, directive, que_suis_je());
  else
    creer_champ(*champ_don, z, type, noms[0], unites[0], nb_comp, nb_ddl, temps, directive, que_suis_je());
  if ((nature == multi_scalaire) && (champ_fonc))
    {
      champ_fonc->valeur().fixer_nature_du_champ(nature);
      champ_fonc->valeur().fixer_unites(unites);
      champ_fonc->valeur().fixer_noms_compo(noms);
    }
  else if ((nature == multi_scalaire) && (champ_don))
    {
      Cerr << "There is no field of type OWN_PTR(Champ_Don_base) with a multi_scalaire nature." << finl;
      exit();
    }
}
