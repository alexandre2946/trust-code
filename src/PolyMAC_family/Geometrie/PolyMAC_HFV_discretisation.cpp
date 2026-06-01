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

#include <Champ_Fonc_Tabule_Elem_PolyMAC_CDO.h>
#include <PolyMAC_HFV_discretisation.h>
#include <Op_Diff_PolyMAC_HFV_base.h>
#include <Champ_Fonc_Elem_PolyMAC_CDO.h>
#include <Domaine_PolyMAC_HFV.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Schema_Temps_base.h>
#include <Champ_Fonc_Tabule.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <DescStructure.h>
#include <Milieu_base.h>
#include <Operateur.h>
#include <Motcle.h>

Implemente_instanciable(PolyMAC_HFV_discretisation, "PolyMAC_HFV|PolyMAC_P0P1NC", PolyMAC_CDO_discretisation);
// XD PolyMAC_HFV discretisation_base PolyMAC_P0P1NC INHERITS_BRACE PolyMAC_HFV discretization (previously PolyMAC_CDO
// XD_CONT discretization compatible with pb_multi).

Entree& PolyMAC_HFV_discretisation::readOn(Entree& s) { return s; }

Sortie& PolyMAC_HFV_discretisation::printOn(Sortie& s) const { return s; }

/*! @brief Discretizes a field for PolyMAC_HFV based on a discretization directive.
 *
 * @brief The directive is a Motcle such as "vitesse", "pression",
 *  "temperature", "champ_elem" (creates a P0-type field), ...
 *  This method determines the type of field to create based on the element type
 *  and the discretization directive. It then determines the number of degrees of freedom
 *  and sets all field parameters (type, nb_compo, nb_ddl, nb_pas_dt,
 *  name(s), unit(s), and field nature) and associates the Domaine_dis with the field.
 *  See the code for the correspondence between directives and
 *  the type of field created.
 *
 */
void PolyMAC_HFV_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, int nb_pas_dt,
                                                   double temps, OWN_PTR(Champ_Inc_base) &champ, const Nom& sous_type) const
{
  const Domaine_PolyMAC_HFV& domaine_PolyMAC_HFV = ref_cast(Domaine_PolyMAC_HFV, z);

  Motcles motcles(7);
  motcles[0] = "vitesse";     // Standard choice for velocity
  motcles[1] = "pression";    // Standard choice for pressure
  motcles[2] = "temperature"; // Standard choice for temperature
  motcles[3] = "divergence_vitesse"; // Type of field obtained by computing div v
  motcles[4] = "gradient_pression";  // Type of field obtained by computing grad P
  motcles[5] = "champ_elem";    // Create a field at elements (of type P0)
  motcles[6] = "champ_sommets"; // Create a field at vertices (type P1)

  // The velocity field type depends on the element type:
  int zp1 = false, default_nb_comp = 0, rang = motcles.search(directive);
  Nom type_elem = Nom("Champ_Elem_") + que_suis_je(), type_som = "Champ_Som_PolyMAC_HFV", type_champ_scal = zp1 ? type_som : type_elem, type_champ_vitesse =
                                                                                                              zp1 ? "Champ_Arete_PolyMAC_HFV" : Nom("Champ_Face_") + que_suis_je(), type;
  switch(rang)
    {
    case 0:
    case 4:
      type = type_champ_vitesse;
      default_nb_comp = 3;
      break;
    case 1:
    case 2:
    case 3:
      type = type_champ_scal;
      default_nb_comp = 1;
      break;
    case 5:
      type = type_elem;
      default_nb_comp = 1;
      break;
    case 6:
      type = type_som;
      default_nb_comp = 1;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "PolyMAC_HFV_discretisation : " << motcles;

  if (sous_type != NOM_VIDE)
    rang = verifie_sous_type(type, sous_type, directive);

  // If the directive was not understood (or if it is a demande_description)
  // call the parent class:
  if (rang < 0)
    {
      Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, nb_pas_dt, temps, champ);
      return;
    }

  // Compute the number of degrees of freedom
  int nb_ddl = 0;
  if (type.debute_par(type_elem))
    nb_ddl = z.nb_elem();
  else if (type.debute_par(type_som))
    nb_ddl = domaine_PolyMAC_HFV.nb_som();
  else if (type.debute_par(type_champ_vitesse))
    nb_ddl = zp1 ? domaine_PolyMAC_HFV.domaine().nb_aretes() : domaine_PolyMAC_HFV.nb_faces();
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

/*! @brief Same as PolyMAC_HFV_discretisation::discretiser_champ(.
 *
 * @brief .. , Champ_Inc) Common processing for champ_fonc and champ_don.
 *  This method is private (passing an Objet_U is not clean from
 *  the outside ...)
 *
 */
void PolyMAC_HFV_discretisation::discretiser_champ_fonc_don(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp,
                                                            double temps, Objet_U& champ) const
{
  // Two pointers for easy access to champ_don or champ_fonc, depending on the type of the champ object.
  OWN_PTR(Champ_Fonc_base) *champ_fonc = dynamic_cast<OWN_PTR(Champ_Fonc_base)*>(&champ);
  OWN_PTR(Champ_Don_base) *champ_don = dynamic_cast<OWN_PTR(Champ_Don_base)*>(&champ);

  const Domaine_PolyMAC_HFV& domaine_PolyMAC_HFV = ref_cast(Domaine_PolyMAC_HFV, z);

  Motcles motcles(8);
  motcles[0] = "pression";    // Standard choice for pressure
  motcles[1] = "temperature"; // Standard choice for temperature
  motcles[2] = "divergence_vitesse"; // Type of field obtained by computing div v
  motcles[3] = "champ_elem";  // Create a field at elements (of type P0)
  motcles[4] = "vitesse";     // Standard choice for velocity
  motcles[5] = "gradient_pression";  // Type of field obtained by computing grad P
  motcles[6] = "champ_sommets";  // Create a field at vertices (of type P1)
  motcles[7] = "champ_face";     // Standard choice for velocity

  // The velocity field type depends on the element type:
  int zp1 = false, default_nb_comp = 0, rang = motcles.search(directive);
  Nom type_elem("Champ_Fonc_Elem_PolyMAC_CDO"), type_som("Champ_Fonc_Som_PolyMAC_CDO"), type_scal = zp1 ? type_som : type_elem, type_champ_vitesse(
                                                                                                      zp1 ? "Champ_Fonc_Arete_PolyMAC_HFV" : "Champ_Fonc_Face_PolyMAC_CDO"), type;
  switch(rang)
    {
    case 0:
    case 1:
    case 2:
      type = type_scal;
      default_nb_comp = 1;
      break;
    case 4:
    case 5:
      type = type_champ_vitesse;
      default_nb_comp = 3;
      break;
    case 3:
      type = type_elem;
      default_nb_comp = 1;
      break;
    case 6:
      type = type_som;
      default_nb_comp = 1;
      break;
    case 7:
      type = "Champ_Fonc_Face_PolyMAC_CDO";
      default_nb_comp = 3;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "PolyMAC_HFV_discretisation : " << motcles;

  // If the directive was not understood (or if it is a demande_description)
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
  if (type == "Champ_Fonc_Elem_PolyMAC_CDO")
    nb_ddl = z.nb_elem();
  else if (type == "Champ_Fonc_Face_PolyMAC_CDO")
    nb_ddl = domaine_PolyMAC_HFV.nb_faces();
  else if (type == "Champ_Fonc_Som_PolyMAC_CDO")
    nb_ddl = domaine_PolyMAC_HFV.nb_som();
  else if (type == "Champ_Fonc_Arete_PolyMAC_HFV")
    nb_ddl = domaine_PolyMAC_HFV.domaine().nb_aretes();
  else
    assert(0);

  // If it is a multi-scalar field:
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

void PolyMAC_HFV_discretisation::y_plus(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  Cerr << "Discretization of y plus" << finl; // Used as a global wall-distance model
  Noms noms(1), unites(1);
  noms[0] = Nom("y_plus");
  unites[0] = Nom("adimensionnel");
  discretiser_champ(Motcle("champ_elem"), z, scalaire, noms, unites, 1, 0, ch);
  DoubleTab& tab_y_p = ch->valeurs();
  for (int i = 0; i < tab_y_p.dimension_tot(0); i++)
    for (int n = 0; n < tab_y_p.dimension_tot(1); n++)
      tab_y_p(i, n) = -1.;
}

void PolyMAC_HFV_discretisation::grad_u(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  abort();
}

Nom PolyMAC_HFV_discretisation::get_name_of_type_for(const Nom& class_operateur, const Nom& type_operateur, const Equation_base& eqn, const OBS_PTR(Champ_base) &champ_sup) const
{
  Nom type;
  Nom type_ch = eqn.inconnue().que_suis_je();
  if (type_ch.debute_par("Champ_Elem"))
    type_ch = "_Elem";
  else if (type_ch.debute_par("Champ_Som"))
    type_ch = "_Som";
  else if (type_ch.debute_par("Champ_Face"))
    type_ch = "_Face";
  else if (type_ch.debute_par("Champ_Arete"))
    type_ch = "_Arete";

  if (class_operateur == "Source")
    type = type_operateur + type_ch + "_" + que_suis_je();
  else if (class_operateur == "Solveur_Masse")
    type = Nom("Masse_") + que_suis_je() + type_ch;
  else if (class_operateur == "Operateur_Grad")
    type = Nom("Op_Grad_") + que_suis_je() + "_Face";
  else if (class_operateur == "Operateur_Div")
    type = Nom("Op_Div_") + que_suis_je();
  else if (class_operateur == "Operateur_Diff")
    type = Nom("Op_Diff") + (type_operateur != "" ? "_" : "") + type_operateur + "_" + que_suis_je() + type_ch;
  else if (class_operateur == "Operateur_Conv")
    type = Nom("Op_Conv_") + type_operateur + "_" + que_suis_je() + type_ch;
  else if (class_operateur == "Operateur_Evanescence")
    type = Nom("Op_Evanescence") + (type_operateur != "" ? "_" : "") + type_operateur + "_" + que_suis_je() + type_ch;
  else
    return Discret_Thyd::get_name_of_type_for(class_operateur, type_operateur, eqn);

  return type;
}

void PolyMAC_HFV_discretisation::residu(const Domaine_dis_base& z, const Champ_Inc_base& ch_inco, OWN_PTR(Champ_Fonc_base) &champ) const
{
  Nom ch_name(ch_inco.le_nom());
  ch_name += "_residu";
  Cerr << "Discretization of " << ch_name << finl;

  Nom type_ch = ch_inco.que_suis_je();
  if (type_ch.debute_par("Champ_Face"))
    {
      Motcle loc = "champ_face";
      Noms nom(1), unites(1);
      nom[0] = ch_name;
      unites[0] = "units_not_defined";
      int nb_comp = ch_inco.valeurs().line_size() * dimension;
      discretiser_champ(loc, z, vectoriel, nom, unites, nb_comp, ch_inco.temps(), champ);
    }
  else
    {
      Motcle loc = "champ_elem";
      Noms nom(1), unites(1);
      nom[0] = ch_name;
      unites[0] = "units_not_defined";
      int nb_comp = ch_inco.valeurs().line_size();
      discretiser_champ(loc, z, scalaire, nom, unites, nb_comp, ch_inco.temps(), champ);
    }

  Champ_Fonc_base& ch_fonc = ref_cast(Champ_Fonc_base, champ.valeur());
  DoubleTab& tab = ch_fonc.valeurs();
  tab = -10000.0;
  Cerr << "[Information] Discretisation_base::residu : the residue is set to -10000.0 at initial time" << finl;
}
