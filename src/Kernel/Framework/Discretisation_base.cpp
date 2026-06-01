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

#include <Discretisation_base.h>
#include <Domaine_dis_cache.h>
#include <Schema_Temps_base.h>
#include <Champ_Fonc_Tabule.h>
#include <Champ_Uniforme.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Interprete.h>
#include <Domaine_VF.h>

Implemente_base(Discretisation_base, "Discretisation_base", Objet_U);

const Motcle Discretisation_base::DEMANDE_DESCRIPTION = Motcle("demande_description");
const Nom Discretisation_base::NOM_VIDE = Nom("-");

Sortie& Discretisation_base::printOn(Sortie& os) const { return os; }

Entree& Discretisation_base::readOn(Entree& is) { return is; }

void Discretisation_base::associer_domaine(const Domaine& dom)
{
  le_domaine_ = dom;
}

// Construction of objects depending on the discretization:
// Fields,
// Pressure matrices,
// ...

// A directive is passed that must be recognized by the type
//  of discretization. Directives are keywords:
//  case is ignored and there should be no spaces.
// Examples of directives:
//  "CHAMP_P0",
//  "CHAMP_FACES",
//  "VITESSE",
//  "PRESSION",
//  "TEMPERATURE",
// Sometimes, the number of components depends on the discretization
// (for example for velocity: 1 component in VDF, 3 in VEF)
// If we set nb_comp = -1, the discretization chooses the number
// appropriate, otherwise it uses the provided value.
void Discretisation_base::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, const Nom& nom, const Nom& unite, int nb_comp, int nb_pas_dt, double temps, OWN_PTR(Champ_Inc_base)& champ,
                                            const Nom& sous_type) const
{
  Noms noms;
  Noms unites;
  noms.add(nom);
  unites.add(unite);

  discretiser_champ(directive, z, scalaire, noms, unites, nb_comp, nb_pas_dt, temps, champ, sous_type);
}

void Discretisation_base::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, const Nom& nom, const Nom& unite, int nb_comp, double temps, OWN_PTR(Champ_Fonc_base)& champ) const
{
  Noms noms;
  Noms unites;
  noms.add(nom);
  unites.add(unite);

  discretiser_champ(directive, z, scalaire, noms, unites, nb_comp, temps, champ);
}

void Discretisation_base::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, const Nom& nom, const Nom& unite, int nb_comp, double temps, OWN_PTR(Champ_Don_base)& champ) const
{
  Noms noms;
  Noms unites;
  noms.add(nom);
  unites.add(unite);

  discretiser_champ(directive, z, scalaire, noms, unites, nb_comp, temps, champ);
}

/*! @brief This function is a tool for the three following methods.
 *
 * ..
 *
 */
void Discretisation_base::test_demande_description(const Motcle& directive, const Nom& type_objet) const
{
  if (directive == Discretisation_base::DEMANDE_DESCRIPTION)
    {
      Cerr << "Discr_base : none directive understood" << finl;
      // We have gone through all descendants, they have written all the directives understood, we stop.
      Process::exit();
    }
  // If we get here, it is because no descendant of the method understood the directive, we cause the display of the directives understood:
  Cerr << "\nError in Discr_base::discretiser_(..., ";
  Cerr << type_objet << ")\n";
  Cerr << " The discretization " << que_suis_je();
  Cerr << " does not understood the following directive :\n  " << directive;
  Cerr << "\n The understood directives are :\n";
}

/*! @brief Discretization of a field depending on the directive and other parameters.
 *
 * See derived classes. This method
 *  does not handle any directive, it displays the list of directives understood
 *  by the derived classes.
 *  See for example VDF_discretisation.cpp and VEF...
 *
 */
void Discretisation_base::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, int nb_pas_dt, double temps,
                                            OWN_PTR(Champ_Inc_base)& champ, const Nom& sous_type) const

{
  test_demande_description(directive, champ.que_suis_je());
  // Recursive call to produce the display of directives: XXX: Elie Saikali -> CODE DOES NOT COMPILE WITHOUT THROW
  if (directive == DEMANDE_DESCRIPTION)
    {
      Process::exit();
      throw;
    }
  discretiser_champ(DEMANDE_DESCRIPTION, z, nature, noms, unites, nb_comp, nb_pas_dt, temps, champ, sous_type);
}

/*! @brief idem
 *
 */
void Discretisation_base::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                            OWN_PTR(Champ_Fonc_base)& champ) const

{
  test_demande_description(directive, champ.que_suis_je());
  // Recursive call to produce the display of directives: XXX: Elie Saikali -> CODE DOES NOT COMPILE WITHOUT THROW
  if (directive == DEMANDE_DESCRIPTION)
    {
      Process::exit();
      throw;
    }
  discretiser_champ(DEMANDE_DESCRIPTION, z, nature, noms, unites, nb_comp, temps, champ);
}

/*! @brief idem
 *
 */
void Discretisation_base::discretiser_champ(const Motcle& directive, const Domaine_dis_base& z, Nature_du_champ nature, const Noms& noms, const Noms& unites, int nb_comp, double temps,
                                            OWN_PTR(Champ_Don_base)& champ) const

{
  test_demande_description(directive, champ.que_suis_je());
  // Recursive call to produce the display of directives: XXX: Elie Saikali -> CODE DOES NOT COMPILE WITHOUT THROW
  if (directive == DEMANDE_DESCRIPTION)
    {
      Process::exit();
      throw;
    }
  discretiser_champ(DEMANDE_DESCRIPTION, z, nature, noms, unites, nb_comp, temps, champ);
}

void Discretisation_base::discretiser_variables() const
{
  Process::exit("Discretisation_base::discretiser_variables() does nothing and must be overloaded !");
}

/*! @brief Tool function to set the common members for all types of fields (used in creer_champ)
 *
 */
void Discretisation_base::champ_fixer_membres_communs(Champ_base& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, double temps)
{
  ch.nommer(nom);
  ch.associer_domaine_dis_base(z);
  ch.fixer_nb_comp(nb_comp);
  if ((nb_comp > 1) && (nb_comp == dimension))
    {
      Noms noms(nb_comp);
      if (!axi)
        {
          noms[0] = nom + "X";
          noms[1] = nom + "Y";
          if (nb_comp > 2)
            noms[2] = nom + "Z";
        }
      else
        {
          noms[0] = nom + "R";
          noms[1] = nom + "teta";
          if (nb_comp > 2)
            noms[2] = nom + "Z";
        }
      ch.fixer_noms_compo(noms);
    }
  ch.fixer_unite(unite);
  ch.fixer_nb_valeurs_nodales(nb_ddl);
  ch.changer_temps(temps);
}

/*! @brief Static method that creates an OWN_PTR(Champ_Inc_base) of the specified type.
 *
 * The parameters "directive" and "nom_discretisation" are
 *  used for display only and are optional
 *
 */
void Discretisation_base::creer_champ(OWN_PTR(Champ_Inc_base)& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, int nb_pas_dt, double temps,
                                      const Nom& directive, const Nom& nom_discretisation)
{
  //Nom nomd = nom_discretisation; // To work around the problem of "static" in type_info::nom()
  ch.typer(type);
  Champ_Inc_base& chb = ch.valeur();
  chb.fixer_nb_valeurs_temporelles(nb_pas_dt);
  champ_fixer_membres_communs(chb, z, type, nom, unite, nb_comp, nb_ddl, temps);
}

/*! @brief Static method that creates an OWN_PTR(Champ_Fonc_base) of the specified type.
 *
 * The parameters "directive" and "nom_discretisation" are
 *  used for display only and are optional
 *
 */
void Discretisation_base::creer_champ(OWN_PTR(Champ_Fonc_base)& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, double temps, const Nom& directive,
                                      const Nom& nom_discretisation)
{
  //Nom nomd = nom_discretisation; // To work around the problem of "static" in type_info::nom()
  ch.typer(type);
  Champ_Fonc_base& chb = ch.valeur();
  champ_fixer_membres_communs(chb, z, type, nom, unite, nb_comp, nb_ddl, temps);
}

/*! @brief Static method that creates an OWN_PTR(Champ_Don_base) of the specified type.
 *
 * The parameters "directive" and "nom_discretisation" are
 *  used for display only and are optional
 *
 */
void Discretisation_base::creer_champ(OWN_PTR(Champ_Don_base)& ch, const Domaine_dis_base& z, const Nom& type, const Nom& nom, const Nom& unite, int nb_comp, int nb_ddl, double temps, const Nom& directive,
                                      const Nom& nom_discretisation)
{
  //Nom nomd = nom_discretisation; // To work around the problem of "static" in type_info::nom()
  ch.typer(type);
  Champ_Don_base& chb = ch.valeur();
  champ_fixer_membres_communs(chb, z, type, nom, unite, nb_comp, nb_ddl, temps);
}

Domaine_dis_base& Discretisation_base::discretiser() const
{
  Nom type = "Domaine_", dis = que_suis_je();
  if (dis == "VEFPreP1B") dis = "VEF";
  type += dis;
  const Domaine& dom = le_domaine_.valeur();
  return Domaine_dis_cache::Build_or_get(type, dom, this);
}

void Discretisation_base::volume_maille(const Schema_Temps_base& sch, const Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const
{
  Cerr << "Discretization of the field 'volume of meshes'" << finl;
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, z);
  discretiser_champ("champ_elem", domaine_VF, "volume_maille", "m3", 1, sch.temps_courant(), ch);
  Champ_Fonc_base& ch_fonc = ref_cast(Champ_Fonc_base, ch.valeur());
  ch_fonc.valeurs().ref(domaine_VF.volumes());
  ch_fonc.valeurs().promote_scalar_to_dim2(); // keep a (nb_elem,1) view even when referencing a DoubleVect
}

void Discretisation_base::mesh_numbering(const Schema_Temps_base& sch, const Domaine_dis_base& z, OWN_PTR(Champ_Fonc_base)& ch) const
{
  Cerr << "Discretization of the field mesh numbering" << finl;
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, z);
  Noms noms(3);
  Noms unit(3);
  noms[0]="nodes"; // Component 0: node numbering
  noms[1]="cells"; // Component 1: element numbering
  noms[2]="faces"; // Component 2: face numbering
  // ToDo: make VisIt recognize component names:
  // Managing fields with components via TRUST + lata plugin is a pain...
  // Looking forward to VTK to avoid this bridge between TRUST and VisIt/Paraview/Tecplot
  discretiser_champ("champ_elem", domaine_VF, multi_scalaire, noms, unit, 3, sch.temps_courant(), ch);
  ch->nommer("mesh_numbering");
  DoubleTab& tab = ch->valeurs();
  tab=0;
  const IntTab& les_elems = z.domaine().les_elems();
  const IntTab& elem_faces = domaine_VF.elem_faces();
  int nb_elem = les_elems.dimension(0);
  int nb_soms_elem = les_elems.dimension(1);
  int nb_faces_elem = elem_faces.dimension(1);
  for (int elem=0; elem<nb_elem; elem++)
    {
      for (int som = 0; som < nb_soms_elem; som++)
        tab(elem, 0) += les_elems(elem, som);
      tab(elem, 0) /= nb_soms_elem;
      tab(elem, 1) = elem;
      for (int face = 0; face < nb_faces_elem; face++)
        tab(elem, 2) += elem_faces(elem, face);
      tab(elem, 2) /= nb_faces_elem;
    }
}

void Discretisation_base::residu(const Domaine_dis_base&, const Champ_Inc_base&, OWN_PTR(Champ_Fonc_base)&) const
{
  Cerr << "Discret_Thyd::residu() does nothing ! " << que_suis_je() << " needs to overload it !" << finl;
  Process::exit();
}

void Discretisation_base::modifier_champ_tabule(const Domaine_dis_base& domaine_dis, Champ_Fonc_Tabule& ch_tab, const VECT(OBS_PTR(Champ_base)) &ch_inc) const
{
  Cerr << que_suis_je() << " must overload Discretisation_base::modifier_champ_tabule !" << finl;
  Process::exit();
}

void Discretisation_base::nommer_completer_champ_physique(const Domaine_dis_base& domaine_dis, const Nom& nom_champ, const Nom& unite, Champ_base& le_champ, const Probleme_base& pb) const
{
  // we name the field and the unit
  le_champ.nommer(nom_champ);
  le_champ.fixer_unite(unite);
  if (sub_type(Champ_Fonc_Tabule, le_champ))
    {
      Noms& noms_variables = ref_cast(Champ_Fonc_Tabule, le_champ).noms_champs_parametre();
      Noms& noms_pbs = ref_cast(Champ_Fonc_Tabule, le_champ).noms_problemes();
      VECT(OBS_PTR(Champ_base)) les_ch_eq;
      for (int i = 0; i < noms_variables.size(); i++)
        {
          OBS_PTR(Champ_base) champ;
          const Probleme_base& pb_ch = noms_pbs.size() == 0 ? pb : ref_cast(Probleme_base, Interprete::objet(noms_pbs[i]));
          champ = pb_ch.get_champ(Motcle(noms_variables[i]));
          les_ch_eq.add(champ);
        }
      modifier_champ_tabule(domaine_dis, ref_cast(Champ_Fonc_Tabule, le_champ), les_ch_eq);
    }
}

/*! @brief Fills the Nom type depending on the class of operator, the type of operator and the equation.
 *
 * @param (class_operateur) This name corresponds to the type of the base class of the object to be constructed. Examples: "source","Op_conv","Op_diff","Op_div","Op_grad","solveur_masse" Example: get_name_of_type_for("Op_conv","amont",eqn,type); get_name_of_type_for("Op_diff"," ",eqn,type,champ_diffusisivite);
 */

// this allows to specify the behavior for each discretization
Nom Discretisation_base::get_name_of_type_for(const Nom& class_operateur, const Nom& type_operateur, const Equation_base& eqn, const OBS_PTR(Champ_base) &champ_sup) const
{
  Nom type;
  if (class_operateur == "Source")
    {
      type = type_operateur;
      Nom disc = eqn.discretisation().que_suis_je();

      int isQC = eqn.probleme().is_dilatable();

      if (isQC && ((eqn.que_suis_je() != "Transport_K_Epsilon") && (eqn.que_suis_je() != "Transport_K_Epsilon_Bas_Reynolds") && (eqn.que_suis_je() != "Transport_K_Epsilon_Realisable") && (eqn.que_suis_je() != "Transport_K_Epsilon_V2")))
        type += "_QC";

      // except for the boussinesq term
      if (disc == "VEFPreP1B")
        {
          if ((Motcle(type_operateur) == "boussinesq_temperature") || (Motcle(type_operateur) == "boussinesq_concentration") || (Motcle(type_operateur) == "boussinesq"))
            disc = "VEFPreP1B";
          else
            disc = "VEF";
        }

      Nom type_ch = eqn.inconnue().que_suis_je();
      if (type_ch == "Champ_Q1NC") type_ch = "Champ_P1NC";

      type_ch.suffix("Champ_");
      type += "_";
      type += disc;
      type += "_";
      type += type_ch;
      return type;
    }
  else if (class_operateur == "Solveur_Masse")
    {
      type = "Masse_";

      Nom discr = eqn.discretisation().que_suis_je();
      if (discr == "VEFPreP1B") discr = "VEF";

      type += discr;

      Nom type_ch = eqn.inconnue().que_suis_je();

      if (type_ch == "Champ_Q1NC") type_ch = "Champ_P1NC";
      if (type_ch.debute_par("Champ_P0_VDF")) type_ch = "Champ_P0_VDF";
      if (type_ch.debute_par("Champ_Face")) type_ch = "Champ_Face";

      type_ch.suffix("Champ");
      type += type_ch;
      return type;
    }
  else if (class_operateur == "Operateur_Grad")
    {
      type = "Op_Grad_";
      Nom type_pb = eqn.probleme().que_suis_je();

      if (type_pb == "Probleme_SG")
        {
          type += (type_pb.suffix("Probleme_"));
          type += "_";
        }

      Nom discr = eqn.discretisation().que_suis_je();

      type += discr;
      type += "_";
      Nom type_inco = eqn.inconnue().que_suis_je();

      if (type_inco == "Champ_Q1NC") type_inco = "Champ_P1NC";

      type += (type_inco.suffix("Champ_"));

      //Test to apply a gradient to a Champ_P1NC with one component in VEF: Typing to review (revision of operators)
      if ((eqn.inconnue().le_nom() != "vitesse") && (eqn.inconnue().que_suis_je() == "Champ_P1NC")) type = "Op_Grad_P1NC_to_P0";

      //Test to apply a gradient to a Champ_P0 with one component in VDF: Typing to review (revision of operators)
      if ((eqn.inconnue().le_nom() != "vitesse") && (eqn.inconnue().que_suis_je() == "Champ_P0_VDF")) type = "Op_Grad_P0_to_Face";

      return type;
    }
  else if (class_operateur == "Operateur_Div")
    {
      type = "Op_Div_";

      Nom discr = eqn.discretisation().que_suis_je();
      Nom type_inco = eqn.inconnue().que_suis_je();
      type += discr;

      type += "_";
      if (type_inco == "Champ_Q1NC") type_inco = "Champ_P1NC";

      type += (type_inco.suffix("Champ_"));
      return type;
    }
  else if (class_operateur == "Operateur_Diff")
    {
      Nom typ(type_operateur);
      if (typ == "standard") typ = "";

      Cerr << "We treat the diffusive operator of : " << eqn.que_suis_je() << finl;
      type = "Op_Diff_";

      Nom nom_discr = que_suis_je();
      Cerr << "The discretization used is : " << nom_discr << finl;
      assert(champ_sup);

      const Champ_base& diffusivite = champ_sup.valeur();

      if (nom_discr == "VEFPreP1B") nom_discr = "VEF";
      type += nom_discr;
      type += typ;

      // MODIF ELI LAUCOIN (10/12/2007) :
      // I keep the normal behavior for VEF and VDF
      if ((nom_discr == "VDF") || (nom_discr == "VEF"))
        {
          Nom nb_inc;
          // Modif Elie Saikali (Nov 2020)
          if (eqn.probleme().que_suis_je().debute_par("Pb_Multiphase") || (diffusivite.nb_comp() == 1 && nom_discr == "VDF")) nb_inc = "_";
          else if (diffusivite.nb_comp() > 1 && diffusivite.le_nom() == "conductivite") nb_inc = "ANISOTROPE_";
          else if (nom_discr == "VDF" && eqn.diffusion_multi_scalaire())
            nb_inc = "_Multi_inco_Multi_scalar_";
          else
            {
              if (nom_discr == "VEF") nb_inc = "_";
              else nb_inc = "_Multi_inco_";
            }

          type += nb_inc;
        }
      else // however, we modify the general case to avoid adding special cases
        type += "_";

      Nom type_inco = eqn.inconnue().que_suis_je();

      type += (type_inco.suffix("Champ_"));
      if (axi == 1) type += "_Axi";

      return type;
    }
  else if (class_operateur == "Operateur_Conv")
    {
      type = "Op_Conv_";
      type += type_operateur;
      Nom tiret = "_";
      type += tiret;
      Nom discr = que_suis_je();

      // the diffusion operators are common to the VEF and VEFP1B discretizations
      if (discr == "VEFPreP1B") discr = "VEF";

      type += discr;
      if (Motcle(type_operateur) == Motcle("ALE")) return type;
      if (type_operateur != "KEps_Comp")
        {
          type += tiret;
          Nom type_inco = eqn.inconnue().que_suis_je();
          if (type_inco == "Champ_Q1NC") type_inco = "Champ_P1NC";
          if (type_inco.debute_par("Champ_P0_VDF")) type_inco = "Champ_P0_VDF";
          if (type_inco.debute_par("Champ_Face")) type_inco = "Champ_Face";

          type += (type_inco.suffix("Champ_"));

          if (axi == 1)
            if (type_operateur == "quick") type += "_Axi";
        }
      return type;
    }
  else if (class_operateur == "Operateur_Evanescence")
    {
      Nom type_inco = eqn.inconnue().que_suis_je();
      if (type_inco == "Champ_Q1NC") type_inco = "Champ_P1NC";
      if (type_inco.debute_par("Champ_P0_VDF")) type_inco = "Champ_P0_VDF";
      if (type_inco.debute_par("Champ_Face")) type_inco = "Champ_Face";

      type_inco.suffix("Champ");
      type = Nom("Op_Evanescence") + (type_operateur != "" ? "_" : "") + type_operateur + "_" + que_suis_je() + type_inco;
    }
  else
    {
      Cerr << class_operateur << " not understood in get_name_of_type_for of  " << que_suis_je() << finl;
      Process::exit();
    }
  return type;
}

int Discretisation_base::verifie_sous_type(Nom& type, const Nom& sous_type, const Motcle& directive) const
{
  const Type_info *base_info = Type_info::type_info_from_name(sous_type);

  if (base_info)
    {
      if (base_info->has_base(type))
        {
          type = sous_type;
          return 0;
        }

      // Elie Saikali: I'll do better the day we write real C++ not like this... MUST TEST IF THE CAST WORKS AND NOT THE NAME !!!!!!!!!!!!!!!!!!
      if (sous_type.debute_par("Champ_Face_dep_expr") && type == "Champ_Face_VDF")
        {
          type = sous_type;
          return 0;
        }

      Cerr << "Error in " << que_suis_je() << " ::discretiser_champ" << finl;
      Cerr << sous_type << " is not a sub type of " << type << " for " << que_suis_je() << " discretization";
      Cerr << "( directive : \"" << directive << "\")" << finl;
      Process::exit();
    }
  else
    {
      Cerr << "Error in " << que_suis_je() << "::discretiser_champ" << finl;
      Cerr << "Unknown class type " << sous_type << finl;
      Process::exit();
    }
  return -1;
}
