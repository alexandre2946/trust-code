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

#include <VDF_discretisation.h>
#include <Rotationnel_Champ_Face.h>
#include <grad_U_Champ_Face.h>
#include <Correlation_Vec_Sca_VDF.h>
#include <Champ_Fonc_Tabule.h>
#include <Champ_Fonc_Tabule_P0_VDF.h>
#include <Critere_Q_Champ_Face.h>
#include <Fluide_Ostwald.h>
#include <Champ_Ostwald_VDF.h>
#include <Champ_Uniforme.h>
#include <Y_plus_Champ_Face.h>
#include <Reynolds_maille_Champ_Face.h>
#include <Courant_maille_Champ_Face.h>
#include <Taux_cisaillement_P0_VDF.h>
#include <Postraitement_base.h>
#include <Schema_Temps_base.h>
#include <Domaine_Cl_VDF.h>
#include <Navier_Stokes_std.h>
#include <Domaine_VDF.h>
#include <Probleme_base.h>
#include <Champ_P0_VDF.h>
#include <T_paroi_Champ_P0_VDF.h>

Implemente_instanciable(VDF_discretisation, "VDF", Discret_Thyd);
// XD vdf discretisation_base vdf INHERITS_BRACE Finite difference volume discretization.

Entree& VDF_discretisation::readOn(Entree& s)
{
  return Discret_Thyd::readOn(s);
}

Sortie& VDF_discretisation::printOn(Sortie& s) const { return s; }

/*! @brief Discretizes a OWN_PTR(Champ_Inc_base) for VDF based on a discretization directive.
 *
 * The directive is a Motcle such as "vitesse", "pression",
 *  "temperature", "champ_elem" (creates a P0-type field ...), ...
 *  This method determines the type of field to create based on the element type
 *  and the discretization directive. It then determines the number of dofs
 *  and sets all field parameters (type, nb_compo, nb_ddl, nb_pas_dt,
 *  name(s), unit(s), nature of the field and assigns a time) and associates the Domaine_dis with the field.
 *  See the code for the correspondence between directives and
 *  the type of field created.
 *
 */
void VDF_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, int nb_pas_dt, double temps,
                                           OWN_PTR(Champ_Inc_base) &champ, const Nom& sous_type) const
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);

  Motcles motcles(7);
  motcles[0] = "vitesse";     // Standard choice for velocity
  motcles[1] = "pression";    // Standard choice for pressure
  motcles[2] = "temperature"; // Standard choice for temperature
  motcles[3] = "divergence_vitesse"; // Field type obtained by computing div v
  motcles[4] = "gradient_pression";  // Field type obtained by computing grad P
  motcles[5] = "champ_face";   // Create a face field
  motcles[6] = "champ_elem";   // Create an element field (P0 type)

  Nom type;
  int default_nb_comp = -1;
  int rang = motcles.search(directive);
  switch(rang)
    {
    case 0:
      type = "Champ_Face_VDF";
      default_nb_comp = dimension;
      break;
    case 1:
      type = "Champ_P0_VDF";
      default_nb_comp = 1;
      break;
    case 2:
      type = "Champ_P0_VDF";
      default_nb_comp = 1;
      break;
    case 3:
      type = "Champ_P0_VDF";
      default_nb_comp = 1;
      break;
    case 4:
      type = "Champ_Face_VDF";
      default_nb_comp = dimension;
      break;
    case 5:
      type = "Champ_Face_VDF";
      default_nb_comp = 1;
      break;
    case 6:
      type = "Champ_P0_VDF";
      default_nb_comp = 1;
      break;
    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "VDF_discretisation : " << motcles;

  if (sous_type != NOM_VIDE)
    rang = verifie_sous_type(type, sous_type, directive);

  // If the directive was not understood (or if it is a demande_description),
  // call the ancestor:
  if (rang < 0)
    {
      Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, nb_pas_dt, temps, champ);
      return;
    }

  // Compute the number of dofs
  int nb_ddl = 0;
  if (type.debute_par("Champ_P0_VDF"))
    nb_ddl = z.nb_elem();
  else if (type.debute_par("Champ_Face"))
    nb_ddl = domaine_vdf.nb_faces();
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

/*! @brief Same as VDF_discretisation::discretiser_champ(.
 *
 * .. , Champ_Inc) for a Champ_Fonc.
 *
 */
void VDF_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                           OWN_PTR(Champ_Fonc_base) &champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

/*! @brief Same as VDF_discretisation::discretiser_champ(.
 *
 * .. , Champ_Inc) for a Champ_Don.
 *
 */
void VDF_discretisation::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                           OWN_PTR(Champ_Don_base)& champ) const
{
  discretiser_champ_fonc_don(directive, z, nature, noms, unites, nb_comp, temps, champ);
}

/*! @brief Same as VEF_discretisation::discretiser_champ(.
 *
 * .. , Champ_Inc). Common treatment for champ_fonc and champ_don.
 *  This method is private (passing an Objet_U is not clean from
 *  an external viewpoint ...)
 *
 */
void VDF_discretisation::discretiser_champ_fonc_don(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                                    Objet_U& champ) const
{
  // Two pointers to access the champ_don or champ_fonc easily, depending on the type of the champ object.
  OWN_PTR(Champ_Fonc_base) *champ_fonc = dynamic_cast<OWN_PTR(Champ_Fonc_base)*>(&champ);
  OWN_PTR(Champ_Don_base) *champ_don = dynamic_cast<OWN_PTR(Champ_Don_base)*>(&champ);

  Motcles motcles(8);
  motcles[0] = "pression";
  motcles[1] = "temperature";
  motcles[2] = "divergence_vitesse"; // Field type obtained by computing div v
  motcles[3] = "champ_elem";  // Create an element field (P0 type)
  motcles[4] = "vitesse";     // Create a field like velocity
  motcles[5] = "gradient_pression";  // Field type obtained by computing grad P
  motcles[6] = "champ_face";   // Create a face field
  motcles[7] = "champ_sommets";   // Create a vertex field

  Nom type;
  int rang = motcles.search(directive);
  int default_nb_comp = 1;
  switch(rang)
    {
    case 0:
      type = "Champ_Fonc_P0_VDF";
      break;
    case 1:
      type = "Champ_Fonc_P0_VDF";
      break;
    case 2:
      type = "Champ_Fonc_P0_VDF";
      break;
    case 3:
      type = "Champ_Fonc_P0_VDF";
      break;
    case 4:
      type = "Champ_Fonc_Face_VDF";
      default_nb_comp = dimension;
      break;
    case 5:
      type = "Champ_Fonc_Face_VDF";
      default_nb_comp = dimension;
      break;
    case 6:
      type = "Champ_Fonc_Face_VDF";
      default_nb_comp = 1;
      break;
    case 7:
      type = "Champ_Fonc_Q1_VDF";
      default_nb_comp = 1;
      break;

    default:
      assert(rang < 0);
      break;
    }

  if (directive == DEMANDE_DESCRIPTION)
    Cerr << "VDF_discretisation : " << motcles;

  // If the directive was not understood (or if it is a demande_description),
  // call the ancestor:
  if (rang < 0)
    {
      if (champ_fonc)
        Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, temps, *champ_fonc);
      else
        Discret_Thyd::discretiser_champ(directive, z, nature, noms, unites, nb_comp, temps, *champ_don);
      return;
    }

  // Compute the number of dofs
  int nb_ddl = 0;
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  if (type == "Champ_Fonc_P0_VDF")
    nb_ddl = z.nb_elem();
  else if (type == "Champ_Fonc_Face_VDF")
    nb_ddl = domaine_vdf.nb_faces();
  else if (type == "Champ_Fonc_Q1_VDF")
    nb_ddl = domaine_vdf.nb_som();
  else
    assert(0);

  if (nb_comp < 0)
    nb_comp = default_nb_comp; // Default number of components
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

void VDF_discretisation::vorticite(Domaine_dis_base& z, const Champ_Inc_base& ch_vitesse,
                                   OWN_PTR(Champ_Fonc_base) &ch) const
{
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  ch.typer("Rotationnel_Champ_Face");
  Rotationnel_Champ_Face& ch_W = ref_cast(Rotationnel_Champ_Face, ch.valeur());
  ch_W.associer_domaine_dis_base(domaine_vdf);
  ch_W.associer_champ(vit);
  ch_W.nommer("vorticite");
  if (dimension == 2)
    ch_W.fixer_nb_comp(1);
  else
    {
      ch_W.fixer_nb_comp(dimension);
      ch_W.fixer_nom_compo(0, "vorticiteX");
      ch_W.fixer_nom_compo(1, "vorticiteY");
      ch_W.fixer_nom_compo(2, "vorticiteZ");
    }
  ch_W.fixer_nb_valeurs_nodales(domaine_vdf.nb_elem());
  ch_W.fixer_unite("s-1");
  ch_W.changer_temps(ch_vitesse.temps());
}

void VDF_discretisation::critere_Q(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  const Domaine_Cl_VDF& domaine_cl_vdf = ref_cast(Domaine_Cl_VDF, zcl);
  ch.typer("Critere_Q_Champ_Face");
  Critere_Q_Champ_Face& ch_Criter_Q = ref_cast(Critere_Q_Champ_Face, ch.valeur());
  ch_Criter_Q.associer_domaine_dis_base(domaine_vdf);
  ch_Criter_Q.associer_domaine_Cl_dis_base(domaine_cl_vdf);
  ch_Criter_Q.associer_champ(vit);
  ch_Criter_Q.nommer("Critere_Q");
  ch_Criter_Q.fixer_nb_comp(1);
  ch_Criter_Q.fixer_nb_valeurs_nodales(domaine_vdf.nb_elem());
  ch_Criter_Q.fixer_unite("s-2");
  ch_Criter_Q.changer_temps(ch_vitesse.temps());
}

void VDF_discretisation::grad_u(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  const Domaine_Cl_VDF& domaine_cl_vdf = ref_cast(Domaine_Cl_VDF, zcl);
  const int N = ch_vitesse.valeurs().line_size();
  ch.typer("grad_U_Champ_Face");
  grad_U_Champ_Face& ch_grad_u = ref_cast(grad_U_Champ_Face, ch.valeur());
  ch_grad_u.associer_domaine_dis_base(domaine_vdf);
  ch_grad_u.associer_domaine_Cl_dis_base(domaine_cl_vdf);
  ch_grad_u.associer_champ(vit);
  ch_grad_u.nommer("gradient_vitesse");
  ch_grad_u.fixer_nb_comp(dimension * dimension * N);

  for (int n = 0; n < N; n++)
    if (dimension == 2)
      {
        ch_grad_u.fixer_nom_compo(N * 0 + n, "dUdX"); // du/dx
        ch_grad_u.fixer_nom_compo(N * 1 + n, "dUdY"); // du/dy
        ch_grad_u.fixer_nom_compo(N * 2 + n, "dVdX"); // dv/dx
        ch_grad_u.fixer_nom_compo(N * 3 + n, "dVdY"); // dv/dy
      }
    else
      {
        ch_grad_u.fixer_nom_compo(N * 0 + n, "dUdX"); // du/dx
        ch_grad_u.fixer_nom_compo(N * 1 + n, "dUdY"); // du/dy
        ch_grad_u.fixer_nom_compo(N * 2 + n, "dUdZ"); // du/dz
        ch_grad_u.fixer_nom_compo(N * 3 + n, "dVdX"); // dv/dx
        ch_grad_u.fixer_nom_compo(N * 4 + n, "dVdY"); // dv/dy
        ch_grad_u.fixer_nom_compo(N * 5 + n, "dVdZ"); // dv/dz
        ch_grad_u.fixer_nom_compo(N * 6 + n, "dWdX"); // dw/dx
        ch_grad_u.fixer_nom_compo(N * 7 + n, "dWdY"); // dw/dy
        ch_grad_u.fixer_nom_compo(N * 8 + n, "dWdZ"); // dw/dz
      }
  ch_grad_u.fixer_nature_du_champ(vectoriel);
  ch_grad_u.fixer_nb_valeurs_nodales(domaine_vdf.nb_elem());
  ch_grad_u.fixer_unite("s-1");
  ch_grad_u.changer_temps(ch_vitesse.temps());
}

void VDF_discretisation::reynolds_maille(const Domaine_dis_base& z, const Fluide_base& le_fluide, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &champ) const
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  champ.typer("Reynolds_maille_Champ_Face");
  Reynolds_maille_Champ_Face& ch = ref_cast(Reynolds_maille_Champ_Face, champ.valeur());
  ch.associer_domaine_dis_base(domaine_vdf);
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
  const Champ_Don_base& nu = le_fluide.viscosite_cinematique();
  ch.associer_champ(vit, nu);
  ch.nommer("Reynolds_maille");
  ch.fixer_nb_comp(dimension);
  ch.fixer_nom_compo(0, "Reynolds_maille_X");
  ch.fixer_nom_compo(1, "Reynolds_maille_Y");
  if (dimension == 3)
    ch.fixer_nom_compo(2, "Reynolds_maille_Z");
  ch.fixer_nb_valeurs_nodales(domaine_vdf.nb_faces());
  ch.fixer_unite("adimensionnel");
  ch.changer_temps(ch_vitesse.temps());
}

void VDF_discretisation::courant_maille(const Domaine_dis_base& z, const Schema_Temps_base& sch, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &champ) const
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  champ.typer("Courant_maille_Champ_Face");
  Courant_maille_Champ_Face& ch = ref_cast(Courant_maille_Champ_Face, champ.valeur());
  ch.associer_domaine_dis_base(domaine_vdf);
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
  ch.associer_champ(vit, sch);
  ch.nommer("Courant_maille");
  ch.fixer_nb_comp(dimension);
  ch.fixer_nom_compo(0, "Courant_maille_X");
  ch.fixer_nom_compo(1, "Courant_maille_Y");
  if (dimension == 3)
    ch.fixer_nom_compo(2, "Courant_maille_Z");
  ch.fixer_nb_valeurs_nodales(domaine_vdf.nb_faces());
  ch.fixer_unite("adimensionnel");
  ch.changer_temps(ch_vitesse.temps());
}

void VDF_discretisation::taux_cisaillement(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &champ) const
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  const Domaine_Cl_VDF& domaine_cl_vdf = ref_cast(Domaine_Cl_VDF, zcl);
  champ.typer("Taux_cisaillement_P0_VDF");
  Taux_cisaillement_P0_VDF& ch = ref_cast(Taux_cisaillement_P0_VDF, champ.valeur());
  ch.associer_domaine_dis_base(domaine_vdf);
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
  ch.associer_champ(vit, domaine_cl_vdf);
  ch.nommer("Taux_cisaillement");
  int N = vit.valeurs().line_size();
  ch.fixer_nb_comp(N);
  for (int n = 0; n < N; n++)
    {
      Nom phase = Nom(n);
      ch.fixer_nom_compo(n, Nom("Taux_cisaillement_") + phase);
    }
  ch.fixer_nb_valeurs_nodales(domaine_vdf.nb_elem());
  ch.fixer_unite("s-1");
  ch.changer_temps(ch_vitesse.temps());
}

void VDF_discretisation::y_plus(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  const Domaine_Cl_VDF& domaine_cl_vdf = ref_cast(Domaine_Cl_VDF, zcl);
  if (domaine_cl_vdf.equation().probleme().que_suis_je().debute_par("Pb_Multiphase"))
    {
      Cerr << "Discretization of y plus" << finl; // Used as a global wall distance model
      Noms noms(1), unites(1);
      noms[0] = Nom("Y_plus");
      unites[0] = Nom("adimensionnel");
      discretiser_champ(Motcle("champ_elem"), z, scalaire, noms, unites, 1, 0, ch);
      DoubleTab& tab_y_p = ch->valeurs();
      for (int i = 0; i < tab_y_p.dimension_tot(0); i++)
        for (int n = 0; n < tab_y_p.dimension_tot(1); n++)
          tab_y_p(i, n) = -1.;
    }
  else
    {
      ch.typer("Y_plus_Champ_Face");
      Y_plus_Champ_Face& ch_y_plus = ref_cast(Y_plus_Champ_Face, ch.valeur());
      ch_y_plus.associer_domaine_dis_base(domaine_vdf);
      ch_y_plus.associer_domaine_Cl_dis_base(domaine_cl_vdf);
      ch_y_plus.associer_champ(vit);
      ch_y_plus.nommer("Y_plus");
      ch_y_plus.fixer_nb_comp(1);
      ch_y_plus.fixer_nb_valeurs_nodales(domaine_vdf.nb_elem());
      ch_y_plus.fixer_unite("adimensionnel");
      ch_y_plus.changer_temps(ch_vitesse.temps());
    }
}

void VDF_discretisation::t_paroi(const Domaine_dis_base& z,const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_temp, OWN_PTR(Champ_Fonc_base)& ch) const
{
  Cerr << "Discretization of temperature_paroi" << finl;
  const Champ_P0_VDF& temp = ref_cast(Champ_P0_VDF, ch_temp);
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  const Domaine_Cl_VDF& domaine_cl_vdf = ref_cast(Domaine_Cl_VDF, zcl);
  ch.typer("T_paroi_Champ_P0_VDF");
  T_paroi_Champ_P0_VDF& ch_tp = ref_cast(T_paroi_Champ_P0_VDF, ch.valeur());
  ch_tp.associer_domaine_dis_base(domaine_vdf);
  ch_tp.associer_domaine_Cl_dis_base(domaine_cl_vdf);
  ch_tp.associer_champ(temp);
  ch_tp.nommer("temperature_paroi");
  ch_tp.add_synonymous("wall_temperature");
  ch_tp.fixer_nb_comp(ch_temp.valeurs().line_size()); // for multiphase ...
  ch_tp.fixer_nb_valeurs_nodales(domaine_vdf.nb_elem());
  ch_tp.fixer_unite("K-C");
  ch_tp.changer_temps(ch_temp.temps());
}

void VDF_discretisation::modifier_champ_tabule(const Domaine_dis_base& domaine_dis, Champ_Fonc_Tabule& le_champ_tabule, const VECT(OBS_PTR(Champ_base)) &ch_inc) const
{
  le_champ_tabule.typer_champ_tabule_discretise("Champ_Fonc_Tabule_P0_VDF");
  Champ_Fonc_Tabule_P0_VDF& le_champ_tabule_dis = ref_cast(Champ_Fonc_Tabule_P0_VDF, le_champ_tabule.le_champ_tabule_discretise());
  le_champ_tabule_dis.associer_domaine_dis_base(domaine_dis);
  le_champ_tabule_dis.associer_param(ch_inc, le_champ_tabule.table());
  le_champ_tabule_dis.nommer(le_champ_tabule.le_nom()); // We give a name to this field, help for debug
  le_champ_tabule_dis.fixer_nb_comp(le_champ_tabule.nb_comp());
  le_champ_tabule_dis.fixer_nb_valeurs_nodales(domaine_dis.nb_elem());
  le_champ_tabule_dis.changer_temps(ch_inc[0]->temps());
}

/*! @brief Discretizes the incompressible fluid in VDF, including K, e, N.
 *
 * @param z domain to discretize
 * @param le_fluide fluid to discretize
 * @param eqn_hydr Navier-Stokes equation (used to access velocity)
 */
void VDF_discretisation::proprietes_physiques_fluide_Ostwald(const Domaine_dis_base& z, Fluide_Ostwald& le_fluide, const Navier_Stokes_std& eqn_hydr, const Champ_Inc_base&) const
{
  Cerr << "Discretization of fluide_Ostwald" << finl;
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);
  const Champ_Inc_base& ch_vitesse = eqn_hydr.inconnue();
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);

  Champ_Don_base& mu = le_fluide.viscosite_dynamique();
  //  mu is always a champ_Ostwald_VDF, so the following must always be done
  Champ_Ostwald_VDF& ch_mu = ref_cast(Champ_Ostwald_VDF, mu);
  Cerr << "associating domainedisbase" << finl;
  ch_mu.associer_domaine_dis_base(domaine_vdf);
  ch_mu.associer_fluide(le_fluide);
  ch_mu.associer_champ(vit);
  Cerr << "associations done" << finl;
  ch_mu.fixer_nb_comp(1);

  Cerr << "calling fixer_nb_valeurs_nodales" << finl;
  Cerr << "nb_valeurs_nodales = " << domaine_vdf.nb_elem() << finl;
  ch_mu.fixer_nb_valeurs_nodales(domaine_vdf.nb_elem());

  Cerr << "calling changer_temps" << finl;
  ch_mu.changer_temps(vit.temps());

  Cerr << "mu is discretized " << finl;
}

void VDF_discretisation::creer_champ_vorticite(const Schema_Temps_base& sch, const Champ_Inc_base& ch_vitesse,
                                               OWN_PTR(Champ_Fonc_base) &ch) const
{
  if (sub_type(Champ_Face_VDF, ch_vitesse))
    {
      const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, ch_vitesse);
      const Domaine_VDF& domaine_VDF = ref_cast(Domaine_VDF, vit.domaine_dis_base());
      int N = ch_vitesse.valeurs().line_size();
      ch.typer("Rotationnel_Champ_Face");
      Rotationnel_Champ_Face& ch_W = ref_cast(Rotationnel_Champ_Face, ch.valeur());
      ch_W.associer_domaine_dis_base(domaine_VDF);
      ch_W.associer_champ(vit);
      ch_W.nommer("vorticite");
      if (dimension == 2)
        ch_W.fixer_nb_comp(N);
      else
        for (int n = 0; n < N; n++)
          {
            ch_W.fixer_nb_comp(dimension * N);
            ch_W.fixer_nom_compo(N * 0 + n, "vorticiteX");
            ch_W.fixer_nom_compo(N * 1 + n, "vorticiteY");
            ch_W.fixer_nom_compo(N * 2 + n, "vorticiteZ");
          }
      ch_W.fixer_nb_valeurs_nodales(domaine_VDF.nb_elem());
      ch_W.fixer_unite("s-1");
      ch_W.changer_temps(sch.temps_courant());
    }
}

void VDF_discretisation::residu(const Domaine_dis_base& z, const Champ_Inc_base& ch_inco, OWN_PTR(Champ_Fonc_base) &champ) const
{
  Nom ch_name(ch_inco.le_nom());
  ch_name += "_residu";
  Cerr << "Discretization of " << ch_name << finl;

  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, z);

  Motcle loc;
  int nb_comp;
  Nom type_ch = ch_inco.que_suis_je();
  if (type_ch.debute_par("Champ_Face"))
    {
      loc = "champ_face";
      nb_comp = dimension * ch_inco.valeurs().line_size();
    }
  else
    {
      loc = "champ_elem";
      nb_comp = ch_inco.valeurs().line_size();
    }

  Discretisation_base::discretiser_champ(loc, domaine_vdf, ch_name, "units_not_defined", nb_comp, ch_inco.temps(), champ);
  Champ_Fonc_base& ch_fonc = ref_cast(Champ_Fonc_base, champ.valeur());
  DoubleTab& tab = ch_fonc.valeurs();
  tab = -10000.0;
  Cerr << "[Information] Discretisation_base::residu : the residue is set to -10000.0 at initial time" << finl;
}

void VDF_discretisation::distance_paroi_globale(const Schema_Temps_base& sch, Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base) &ch) const
{
  Cerr << "Discretization of global wall distance" << finl;
  Noms noms(1), unites(1);
  noms[0] = Nom("distance_paroi_globale");
  unites[0] = Nom("m");
  discretiser_champ(Motcle("champ_elem"), z, scalaire, noms, unites, 1, 0, ch);
}
