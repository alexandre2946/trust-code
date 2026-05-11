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

#include <DG_discretisation.h>
#include <Domaine_DG.h>
#include <Champ_Fonc_Tabule.h>
#include <Milieu_base.h>
#include <Equation_base.h>
#include <Champ_Uniforme.h>
#include <Champ_Inc_base.h>
#include <Schema_Temps_base.h>
#include <Motcle.h>
#include <Domaine_Cl_DG.h>
#include <Option_DG.h>
#include <Quadrature_base.h>
#include <Quadrature_Ord1_Polygone.h>
#include <Quadrature_Ord3_Polygone.h>
#include <Quadrature_Ord5_Polygone.h>

Implemente_instanciable(DG_discretisation, "DG", Discret_Thyd);
// XD DG discretisation_base DG INHERITS_BRACE DG discretization


Entree& DG_discretisation::readOn(Entree& s)
{
  return Discret_Thyd::readOn(s);
}

Sortie& DG_discretisation::printOn(Sortie& s) const { return s; }

/*! @brief Discretisation of a field for DG formulation
 *
 * TODO, refactor this after stokes. And decide if we want to put this on Option_DG or not
 *
 * La directive est un Motcle comme "vitesse", "pression",
 *  "temperature", "champ_elem" (cree un champ de type P0), ...
 *  Cette methode determine le type du champ a creer en fonction du type d'element
 *  et de la directive de discretisation. Elle determine ensuite le nombre de ddl
 *  et fixe l'ensemble des parametres du champ (type, nb_compo, nb_ddl, nb_pas_dt,
 *  nom(s), unite(s) et nature du champ) et associe la Domaine_dis au champ.
 *  Voir le code pour avoir la correspondance entre les directives et
 *  le type de champ cree.
 *
 */
void DG_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& dom_dis, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, int nb_pas_dt,
                                          double temps, OWN_PTR(Champ_Inc_base)& champ, const Nom& sous_type) const
{
  Motcles motcles(7);
  motcles[0] = "vitesse";     // Choix standard pour la vitesse
  motcles[1] = "pression";    // Choix standard pour la pression
  motcles[2] = "temperature"; // Choix standard pour la temperature
  motcles[3] = "divergence_vitesse"; // Le type de champ obtenu en calculant div v
  motcles[4] = "gradient_pression";  // Le type de champ obtenu en calculant grad P
  motcles[5] = "champ_elem";    // Creer un champ aux elements (de type P0)
  motcles[6] = "champ_sommets"; // Creer un champ aux sommets (type P1)
  // Le type de champ de vitesse depend du type d'element :
//  Nom type_champ_vitesse("Champ_Face_DG");
  Nom type_elem("Champ_Elem_DG");
  Nom type;
  int nb_basis_func = 0; // Valeur par defaut du nombre de composantes
  int rang = motcles.search(directive);
  const int order_DG = Option_DG::Get_order_for(noms[0]);
  switch(rang)
    {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      type = type_elem;
      nb_basis_func = Option_DG::Nb_col_from_order(order_DG);
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "DG discretisation : " << motcles;

  if (sous_type != NOM_VIDE)
    rang = verifie_sous_type(type, sous_type, directive);

  // Si on n'a pas compris la directive (ou si c'est une demande_description)
  // alors on appelle l'ancetre :
  if (rang < 0)
    {
      Discret_Thyd::discretiser_champ(directive, dom_dis, nature, noms, unites, nb_comp, nb_pas_dt, temps, champ);
      return;
    }

  // Calcul du nombre de ddl
  int nb_ddl = 0;
  if (type.debute_par(type_elem))
    nb_ddl = dom_dis.nb_elem();
  else
    assert(0);

  creer_champ(champ, dom_dis, type, noms[0], unites[0], nb_comp*nb_basis_func, nb_ddl, nb_pas_dt, temps, directive, que_suis_je());

  if (nature == multi_scalaire)
    {
      throw;
    }

  if (nb_comp == 1)
    switch(order_DG)
      {
      case 0:
        champ->fixer_nature_du_champ(scalaire);
        break;
      case 1:
        champ->fixer_nature_du_champ(basis_function_order_1_scalar);
        break;
      case 2:
        champ->fixer_nature_du_champ(basis_function_order_2_scalar);
        break;
      default:
        assert(0);
      }
  else
    switch(order_DG)
      {
      case 0:
        champ->fixer_nature_du_champ(vectoriel);
        break;
      case 1:
        champ->fixer_nature_du_champ(basis_function_order_1_vectorial);
        break;
      case 2:
        champ->fixer_nature_du_champ(basis_function_order_2_vectorial);
        break;
      default:
        assert(0);
      }


}

/*! @brief Idem que DG_discretisation::discretiser_champ(.
 *
 * .. , Champ_Inc)
 *
 */
void DG_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                          OWN_PTR(Champ_Fonc_base)& champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

/*! @brief Idem que DG_discretisation::discretiser_champ(.
 *
 * .. , Champ_Inc)
 *
 */
void DG_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                          OWN_PTR(Champ_Don_base)& champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

/*! @brief Idem que DG_discretisation::discretiser_champ(.
 *
 * .. , Champ_Inc) Traitement commun aux champ_fonc et champ_don.
 *  Cette methode est privee (passage d'un Objet_U pas propre vu
 *  de l'exterieur ...)
 *
 */
void DG_discretisation::discretiser_champ_fonc_don(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                                   Objet_U& champ) const
{
  // Deux pointeurs pour acceder facilement au champ_don ou au champ_fonc, suivant le type de l'objet champ.
  OWN_PTR(Champ_Fonc_base) * champ_fonc = dynamic_cast<OWN_PTR(Champ_Fonc_base)*>(&champ);
  OWN_PTR(Champ_Don_base) * champ_don = dynamic_cast<OWN_PTR(Champ_Don_base)*>(&champ);

  const Domaine_DG& domaine_DG = ref_cast(Domaine_DG, z);


  Motcles motcles(8);
  motcles[0] = "pression";    // Choix standard pour la pression
  motcles[1] = "temperature"; // Choix standard pour la temperature
  motcles[2] = "champ_fonc_quad_dg"; // With value on quadrature points
  motcles[5] = "champ_elem_dg"; // With value on quadrature points
  motcles[3] = "champ_elem";  // Creer un champ aux elements (de type P0)
  motcles[6] = "champ_sommets";  // Creer un champ aux elements (de type P1)
  motcles[4] = "vitesse";     // Choix standard pour la vitesse
  motcles[7] = "champ_face";     // Choix standard pour la vitesse

  Nom type;
  int nb_points = 0; // Valeur par defaut du nombre de composantes
  int rang = motcles.search(directive);
  const Quadrature_base& quad = domaine_DG.get_quadrature(5); // TODO: Make this depend from the order of discretization ...
  int nb_pts_integ_max = quad.nb_pts_integ_max();

//  const int order_DG = Option_DG::Get_order_for(directive);
  switch(rang)
    {
    case 0:
    case 1:
    case 2:
    case 4:
    case 5:
    case 7:
      type = "Champ_Fonc_Quad_DG";
      nb_points = nb_pts_integ_max; //Option_DG::Nb_col_from_order(order_DG);;
      break;
    case 3:
      type = "Champ_Fonc_Elem_DG";
      nb_points = 1;
      break;
    case 6:
      type = "Champ_Fonc_Som_DG";
      nb_points = 1;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "DG discretisation : " << motcles;

  // Si on n'a pas compris la directive (ou si c'est une demande_description)
  // alors on appelle l'ancetre :
  if (rang < 0)
    {
      if (champ_fonc)
        Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, temps, *champ_fonc);
      else
        Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, temps, *champ_don);
      return;
    }

  // Calcul du nombre de ddl
  int nb_ddl = 0;
  if (type == "Champ_Fonc_Elem_DG" || type == "Champ_Fonc_Quad_DG")
    nb_ddl = z.nb_elem();
  else if (type == "Champ_Fonc_Som_DG")
    nb_ddl = domaine_DG.nb_som();
  else
    assert(0);

  bool vector = (nature==vectoriel or nature==quadrature_vectoriel or nature == basis_function_order_1_vectorial or nature == basis_function_order_2_vectorial );
  if (vector) nb_comp = dimension;
  else nb_comp = 1;
  if (champ_fonc)
    {
      creer_champ(*champ_fonc, z, type, noms[0], unites[0], nb_comp*nb_points, nb_ddl, temps, directive, que_suis_je());
      if (nb_comp == 1)
        {
          (nb_points == 1) ? champ_fonc->valeur().fixer_nature_du_champ(scalaire)
          : champ_fonc->valeur().fixer_nature_du_champ(quadrature_scalaire);
        }
      else if (nb_comp == dimension)
        {
          (nb_points == 1) ? champ_fonc->valeur().fixer_nature_du_champ(vectoriel)
          : champ_fonc->valeur().fixer_nature_du_champ(quadrature_vectoriel);
        }
      else
        {
          Cerr << "multi_scalaire not implemented for now" << finl;
          exit();
        }
    }
  else
    {
      creer_champ(*champ_don, z, type, noms[0], unites[0], nb_comp*nb_points, nb_ddl, temps, directive, que_suis_je());
      if (nb_comp == 1)
        {
          (nb_points == 1) ? champ_don->valeur().fixer_nature_du_champ(scalaire)
          : champ_don->valeur().fixer_nature_du_champ(quadrature_scalaire);
        }
      else if (nb_comp == dimension)
        {
          (nb_points == 1) ? champ_don->valeur().fixer_nature_du_champ(vectoriel)
          : champ_don->valeur().fixer_nature_du_champ(quadrature_vectoriel);
        }
      else
        {
          Cerr << "multi_scalaire not implemented for now" << finl;
          exit();
        }
    }



  if ((nature == multi_scalaire) && (champ_fonc))
    {
      throw;
    }
  else if ((nature == multi_scalaire) && (champ_don))
    {
      Cerr << "There is no field of type Champ_Don with a multi_scalaire nature." << finl;
      exit();
    }
}

void DG_discretisation::distance_paroi(const Schema_Temps_base& sch, Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const
{
  throw;
}


void DG_discretisation::grad_u(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base)& ch) const
{
  throw;
}


void DG_discretisation::modifier_champ_tabule(const Domaine_dis_base& domaine_poly, Champ_Fonc_Tabule& lambda_tab, const VECT(OBS_PTR(Champ_base)) &champs_param) const
{
  throw;
}

/*! @brief Old copy paste, give the name to fields. Name used to allocate size of matrixes in discretiser
 *
 */
Nom DG_discretisation::get_name_of_type_for(const Nom& class_operateur, const Nom& type_operateur, const Equation_base& eqn, const OBS_PTR(Champ_base) &champ_sup) const
{
  Nom type;
  if (class_operateur == "Source")
    {
      type = type_operateur;
      Nom champ = (eqn.inconnue().que_suis_je());
      champ.suffix("Champ");
      type += champ;
      //type+="_DG";
      return type;

    }
  else if (class_operateur == "Solveur_Masse")
    {
      Nom type_ch = eqn.inconnue().que_suis_je();
      if (type_ch.debute_par("Champ_Elem"))
        type_ch = "_Elem";

      if (type_ch.debute_par("Champ_Face"))
        type_ch = "_Face";

      type = "Masse_DG";
      type += type_ch;
    }
  else if (class_operateur == "Operateur_Grad")
    {
      type = "Op_Grad_DG";
    }
  else if (class_operateur == "Operateur_Div")
    {
      type = "Op_Div_DG";
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
      type += "_DG";
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
      type += "_DG";
    }

  else
    return Discret_Thyd::get_name_of_type_for(class_operateur, type_operateur, eqn);

  return type;
}

void DG_discretisation::distance_paroi_globale(const Schema_Temps_base& sch, Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const
{
  Cerr << "Global wall distance discretisation" << finl;
  Noms noms(1), unites(1);
  noms[0] = Nom("distance_paroi_globale");
  unites[0] = Nom("m");
  discretiser_champ(Motcle("champ_elem"), z, scalaire, noms, unites, 1, 0, ch);
}
