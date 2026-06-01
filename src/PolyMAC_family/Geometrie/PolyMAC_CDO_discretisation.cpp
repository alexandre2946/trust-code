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

#include <PolyMAC_CDO_discretisation.h>
#include <Domaine_PolyMAC_CDO.h>
#include <Champ_Fonc_Tabule.h>
#include <Champ_Fonc_Elem_PolyMAC_CDO.h>
#include <Champ_Fonc_Tabule_Elem_PolyMAC_CDO.h>
#include <Milieu_base.h>
#include <Equation_base.h>
#include <Champ_Uniforme.h>
#include <DescStructure.h>
#include <Schema_Temps_base.h>
#include <Motcle.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <grad_U_Champ_Face_PolyMAC_CDO.h>

Implemente_instanciable(PolyMAC_CDO_discretisation, "PolyMAC_CDO|PolyMAC", Discret_Thyd);
// XD PolyMAC_CDO discretisation_base PolyMAC INHERITS_BRACE PolyMAC_CDO discretization (PolyMAC_CDO discretization that
// XD_CONT is not compatible with pb_multi).


Entree& PolyMAC_CDO_discretisation::readOn(Entree& s) { return s; }

Sortie& PolyMAC_CDO_discretisation::printOn(Sortie& s) const { return s; }

/*! @brief Discretizes a field for PolyMAC_CDO based on a discretization directive.
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
void PolyMAC_CDO_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, int nb_pas_dt,
                                                   double temps, OWN_PTR(Champ_Inc_base)& champ, const Nom& sous_type) const
{
  const Domaine_PolyMAC_CDO& domaine_PolyMAC_CDO = ref_cast(Domaine_PolyMAC_CDO, z);

  Motcles motcles(7);
  motcles[0] = "vitesse";     // Standard choice for velocity
  motcles[1] = "pression";    // Standard choice for pressure
  motcles[2] = "temperature"; // Standard choice for temperature
  motcles[3] = "divergence_vitesse"; // Type of field obtained by computing div v
  motcles[4] = "gradient_pression";  // Type of field obtained by computing grad P
  motcles[5] = "champ_elem";    // Create a field at elements (of type P0)
  motcles[6] = "champ_sommets"; // Create a field at vertices (type P1)

  // The velocity field type depends on the element type:
  Nom type_champ_vitesse("Champ_Face_PolyMAC_CDO");
  Nom type_elem("Champ_Elem_PolyMAC_CDO");
  Nom type;
  int default_nb_comp = 0; // Default value of the number of components
  int rang = motcles.search(directive);
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
    case 5:
      type = type_elem;
      default_nb_comp = 1;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "PolyMAC_CDO_discretisation : " << motcles;

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
    nb_ddl = z.nb_elem() + domaine_PolyMAC_CDO.nb_faces();
  else if (type.debute_par(type_champ_vitesse))
    nb_ddl = domaine_PolyMAC_CDO.nb_faces() + (dimension < 3 ? domaine_PolyMAC_CDO.nb_som() : domaine_PolyMAC_CDO.domaine().nb_aretes());
  else if (type.debute_par("Champ_Som_PolyMAC_CDO"))
    nb_ddl = domaine_PolyMAC_CDO.nb_som();
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

/*! @brief Same as PolyMAC_CDO_discretisation::discretiser_champ(.
 *
 * @brief .. , Champ_Inc)
 *
 */
void PolyMAC_CDO_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                                   OWN_PTR(Champ_Fonc_base)& champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

/*! @brief Same as PolyMAC_CDO_discretisation::discretiser_champ(.
 *
 * @brief .. , Champ_Inc)
 *
 */
void PolyMAC_CDO_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                                   OWN_PTR(Champ_Don_base)& champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

/*! @brief Same as PolyMAC_CDO_discretisation::discretiser_champ(.
 *
 * @brief .. , Champ_Inc) Common processing for champ_fonc and champ_don.
 *  This method is private (passing an Objet_U is not clean from
 *  the outside ...)
 *
 */
void PolyMAC_CDO_discretisation::discretiser_champ_fonc_don(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                                            Objet_U& champ) const
{
  // Two pointers for easy access to champ_don or champ_fonc, depending on the type of the champ object.
  OWN_PTR(Champ_Fonc_base)  * champ_fonc = dynamic_cast<OWN_PTR(Champ_Fonc_base)*>(&champ);
  OWN_PTR(Champ_Don_base) * champ_don = dynamic_cast<OWN_PTR(Champ_Don_base)*>(&champ);

  const Domaine_PolyMAC_CDO& domaine_PolyMAC_CDO = ref_cast(Domaine_PolyMAC_CDO, z);

  Motcles motcles(8);
  motcles[0] = "pression";    // Standard choice for pressure
  motcles[1] = "temperature"; // Standard choice for temperature
  motcles[2] = "divergence_vitesse"; // Type of field obtained by computing div v
  motcles[3] = "champ_elem";  // Create a field at elements (of type P0)
  motcles[6] = "champ_sommets";  // Create a field at elements (of type P1)
  motcles[4] = "vitesse";     // Standard choice for velocity
  motcles[7] = "champ_face";     // Standard choice for velocity
  motcles[5] = "gradient_pression";  // Type of field obtained by computing grad P

  // The velocity field type depends on the element type:
  Nom type_champ_vitesse("Champ_Fonc_Face_PolyMAC_CDO");
  Nom type_elem("Champ_Fonc_Elem_PolyMAC_CDO");

  Nom type;
  int default_nb_comp = 0; // Default value of the number of components
  int rang = motcles.search(directive);
  switch(rang)
    {
    case 0:
    case 1:
    case 2:
    case 3:
      type = type_elem;
      default_nb_comp = 1;
      break;
    case 4:
    case 5:
    case 7:
      type = type_champ_vitesse;
      default_nb_comp = 3;
      break;
    case 6:
      type = "Champ_Fonc_Som_PolyMAC_CDO";
      default_nb_comp = 1;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "PolyMAC_CDO_discretisation : " << motcles;

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
  else if (type == type_champ_vitesse)
    nb_ddl = domaine_PolyMAC_CDO.nb_faces();
  else if (type == "Champ_Fonc_Som_PolyMAC_CDO")
    nb_ddl = domaine_PolyMAC_CDO.nb_som();
  else
    assert(0);

  // If it is a multi-scalar field, handle accordingly!
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

void PolyMAC_CDO_discretisation::distance_paroi(const Schema_Temps_base& sch, Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const
{
  Cerr << "Discretization of the wall distance" << finl;
  Domaine_PolyMAC_CDO& domaine_PolyMAC_CDO = ref_cast(Domaine_PolyMAC_CDO, z);
  ch.typer("Champ_Fonc_Elem_PolyMAC_CDO");
  Champ_Fonc_Elem_PolyMAC_CDO& ch_dist_paroi = ref_cast(Champ_Fonc_Elem_PolyMAC_CDO, ch.valeur());
  ch_dist_paroi.associer_domaine_dis_base(domaine_PolyMAC_CDO);
  ch_dist_paroi.nommer("distance_paroi");
  ch_dist_paroi.fixer_nb_comp(1);
  ch_dist_paroi.fixer_nb_valeurs_nodales(domaine_PolyMAC_CDO.nb_elem());
  ch_dist_paroi.fixer_unite("m");
  ch_dist_paroi.changer_temps(sch.temps_courant());
}


/*! @brief Discretizes the incompressible fluid in PolyMAC_CDO, i.e., K and N.
 *
 * @brief Postcondition: the method does not modify the object.
 *
 * @param (Domaine_dis_base&) domain to discretize
 * @param (Fluide_Ostwald&) fluid to discretize
 * @param (Champ_Inc_base&) ch_vitesse
 * @param (Champ_Inc_base&) temperature
 */

void PolyMAC_CDO_discretisation::grad_u(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base)& ch) const
{
  const Champ_Face_PolyMAC_CDO& vit = ref_cast(Champ_Face_PolyMAC_CDO, ch_vitesse);
  const Domaine_PolyMAC_CDO& domaine_poly = ref_cast(Domaine_PolyMAC_CDO, z);
  const Domaine_Cl_PolyMAC_family& domaine_cl_poly = ref_cast(Domaine_Cl_PolyMAC_family, zcl);
  ch.typer("grad_U_Champ_Face_PolyMAC_CDO");
  grad_U_Champ_Face_PolyMAC_CDO& ch_grad_u = ref_cast(grad_U_Champ_Face_PolyMAC_CDO, ch.valeur());
  ch_grad_u.associer_domaine_dis_base(domaine_poly);
  ch_grad_u.associer_domaine_Cl_dis_base(domaine_cl_poly);
  ch_grad_u.associer_champ(vit);
  ch_grad_u.nommer("gradient_vitesse");
  ch_grad_u.fixer_nb_comp(dimension * dimension);

  if (dimension == 2)
    {
      ch_grad_u.fixer_nom_compo(0, "dUdX"); // du/dx
      ch_grad_u.fixer_nom_compo(1, "dUdY"); // du/dy
      ch_grad_u.fixer_nom_compo(2, "dVdX"); // dv/dx
      ch_grad_u.fixer_nom_compo(3, "dVdY"); // dv/dy
    }
  else
    {
      ch_grad_u.fixer_nom_compo(0, "dUdX"); // du/dx
      ch_grad_u.fixer_nom_compo(1, "dUdY"); // du/dy
      ch_grad_u.fixer_nom_compo(2, "dUdZ"); // du/dz
      ch_grad_u.fixer_nom_compo(3, "dVdX"); // dv/dx
      ch_grad_u.fixer_nom_compo(4, "dVdY"); // dv/dy
      ch_grad_u.fixer_nom_compo(5, "dVdZ"); // dv/dz
      ch_grad_u.fixer_nom_compo(6, "dWdX"); // dw/dx
      ch_grad_u.fixer_nom_compo(7, "dWdY"); // dw/dy
      ch_grad_u.fixer_nom_compo(8, "dWdZ"); // dw/dz
    }
  ch_grad_u.fixer_nature_du_champ(vectoriel);
  ch_grad_u.fixer_nb_valeurs_nodales(domaine_poly.nb_elem());
  ch_grad_u.fixer_unite("s-1");
  ch_grad_u.changer_temps(ch_vitesse.temps());
}

void PolyMAC_CDO_discretisation::modifier_champ_tabule(const Domaine_dis_base& domaine_poly, Champ_Fonc_Tabule& lambda_tab, const VECT(OBS_PTR(Champ_base)) &champs_param) const
{
  lambda_tab.typer_champ_tabule_discretise("Champ_Fonc_Tabule_Elem_PolyMAC_CDO");
  Champ_Fonc_base& lambda_tab_dis = lambda_tab.le_champ_tabule_discretise();
  Champ_Fonc_Tabule_Elem_PolyMAC_CDO& ch_tab_lambda_dis = ref_cast(Champ_Fonc_Tabule_Elem_PolyMAC_CDO, lambda_tab_dis);
  //ch_tab_lambda_dis.nommer(nom_champ);
  ch_tab_lambda_dis.associer_domaine_dis_base(domaine_poly);
  ch_tab_lambda_dis.associer_param(champs_param, lambda_tab.table());
  ch_tab_lambda_dis.fixer_nb_comp(lambda_tab.nb_comp());
  ch_tab_lambda_dis.fixer_nb_valeurs_nodales(domaine_poly.nb_elem());
  // ch_tab_lambda_dis.fixer_unite(unite);
  ch_tab_lambda_dis.changer_temps(champs_param[0]->temps());
}

Nom PolyMAC_CDO_discretisation::get_name_of_type_for(const Nom& class_operateur, const Nom& type_operateur, const Equation_base& eqn, const OBS_PTR(Champ_base) &champ_sup) const
{
  Nom type;
  if (class_operateur == "Source")
    {
      type = type_operateur;
      Nom champ = (eqn.inconnue().que_suis_je());
      champ.suffix("Champ");
      type += champ;
      //type+="_PolyMAC_CDO";
      return type;

    }
  else if (class_operateur == "Solveur_Masse")
    {
      Nom type_ch = eqn.inconnue().que_suis_je();
      if (type_ch.debute_par("Champ_Elem"))
        type_ch = "_Elem";

      if (type_ch.debute_par("Champ_Face"))
        type_ch = "_Face";

      type = "Masse_PolyMAC_CDO";
      type += type_ch;
    }
  else if (class_operateur == "Operateur_Grad")
    {
      type = "Op_Grad_PolyMAC_CDO_Face";
    }
  else if (class_operateur == "Operateur_Div")
    {
      type = "Op_Div_PolyMAC_CDO";
    }

  else if (class_operateur == "Operateur_Diff")
    {
      Nom type_ch = eqn.inconnue().que_suis_je();
      if (type_ch.debute_par("Champ_Elem"))
        type_ch = "_Elem";

      if (type_ch.debute_par("Champ_Face"))
        type_ch = "_Face";

      type = "Op_Diff";
      if (type_operateur != "")
        {
          type += "_";
          type += type_operateur;
        }
      type += "_PolyMAC_CDO";
      type += type_ch;
    }
  else if (class_operateur == "Operateur_Conv")
    {
      type = "Op_Conv_";
      type += type_operateur;
      Nom tiret = "_";
      type += tiret;
      type += que_suis_je();
      Nom type_ch = eqn.inconnue().que_suis_je();
      if (type_ch.debute_par("Champ_Elem"))
        type += "_Elem";
      if (type_ch.debute_par("Champ_Face"))
        type += "_Face";
    }

  else
    return Discret_Thyd::get_name_of_type_for(class_operateur, type_operateur, eqn);

  return type;
}

void PolyMAC_CDO_discretisation::distance_paroi_globale(const Schema_Temps_base& sch, Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const
{
  Cerr << "Discretization of global wall distance" << finl;
  Noms noms(1), unites(1);
  noms[0] = Nom("distance_paroi_globale");
  unites[0] = Nom("m");
  discretiser_champ(Motcle("champ_elem"), z, scalaire, noms, unites, 1, 0, ch);
}
