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

#include <Equation_rayonnement_base.h>
#include <Modele_rayo_semi_transp.h>
#include <Matrice_Morse_Sym.h>
#include <Champ_Uniforme.h>
#include <Matrice_Bloc.h>
#include <Discret_Thyd.h>
#include <Fluide_base.h>
#include <Param.h>

Implemente_base(Equation_rayonnement_base, "Equation_rayonnement_base", Equation_base);

bool Equation_rayonnement_base::initTimeStep(double dt)
{
  schema_temps().set_dt() = dt;
  return Equation_base::initTimeStep(dt);
}

bool Equation_rayonnement_base::solve()
{
  resoudre(schema_temps().temps_courant() + schema_temps().pas_de_temps());
  return true;
}

Sortie& Equation_rayonnement_base::printOn(Sortie& s) const { return s << que_suis_je() << finl; }

Entree& Equation_rayonnement_base::readOn(Entree& is) { return Equation_base::readOn(is); }

void Equation_rayonnement_base::set_param(Param& param) const
{
  param.ajouter_non_std("conditions_limites|boundary_conditions", (this), Param::REQUIRED);
  param.ajouter_non_std("solveur", (this), Param::REQUIRED);
}

int Equation_rayonnement_base::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  int retval = 1;
  if (mot == "conditions_limites|boundary_conditions")
    {
      lire_cl(is);
      verif_Cl();
    }
  else if (mot == "solveur")
    {
      Cerr << "Reading and typing of the radiation equation solver :" << finl;
      Nom nom_solveur("Solv_");
      Nom type_solv_sys;
      is >> type_solv_sys;
      nom_solveur += type_solv_sys;
      Cerr << "Name of the radiation equation solver : " << nom_solveur << finl;
      solveur.typer(nom_solveur);
      is >> solveur.valeur();
      solveur.nommer("solveur_irradiance");
    }
  else
    retval = -1;

  return retval;
}

void Equation_rayonnement_base::completer()
{
  // typage de l'operateur de diffusion
  if (sub_type(Fluide_base, fluide()))
    if (fluide().is_rayo_semi_transp())
      {
        if (fluide().is_longueur_rayo_discretised())
          {
            const Champ_Don_base& long_rayo = fluide().longueur_rayo();
            terme_diffusif.associer_diffusivite(long_rayo);
          }
        else
          {
            Cerr << "Error while reading the Radiation equation." << finl;
            Cerr << "You may not have discretized the problem of type Pb_Couple_rayo_semi_transp." << finl;
            Cerr << "Remark : The discretization of a coupling problem carries out those of the involved coupled " << finl;
            Cerr << "         problems therefore it is not useful to keep the discretization of these last ones.  " << finl;
            Cerr << finl;
            Process::exit();
          }
      }
    else
      {
        Cerr << "Error : the radiative properties of the incompressible fluid have not" << finl;
        Cerr << "been defined while a semi transparent radiation model is used." << finl;
        Cerr << "The fields kappa and indice which respectively define the absoption coefficient" << finl;
        Cerr << "and the refraction index must be added to the luid properties." << finl;
        Process::exit();
      }
  else
    {
      Cerr << "Error while reading the Radiation equation." << finl;
      Cerr << "Your fluid is of type " << fluide().que_suis_je() << finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered " << finl;
      Cerr << "with the semi transparent radiation model." << finl;
      Process::exit();
    }

  Nom type = "Op_Diff_";
  Nom discr = discretisation().que_suis_je();
  // les operateurs C_D_Turb_T sont communs aux discretisations VEF et VEFP1B
  if (discr == "VEFPreP1B")
    discr = "VEF";
  type += discr;

  Nom nb_inc;
  if (sub_type(Champ_Uniforme, terme_diffusif.diffusivite()))
    nb_inc = "_";
  else
    nb_inc = "_Mult_inco_";
  type += nb_inc;

  Nom type_inco = inconnue().que_suis_je();
  type += (type_inco.suffix("Champ_"));

  if (axi)
    type += "_Axi";

  terme_diffusif.typer(type);
  terme_diffusif.l_op_base().associer_eqn(*this);

  if (sub_type(Fluide_base, fluide()))
    if (fluide().has_kappa())
      {
        const Champ_Don_base& long_rayo = fluide().longueur_rayo();
        terme_diffusif->associer_diffusivite(long_rayo);
        terme_diffusif.completer();
        terme_diffusif->dimensionner(la_matrice);
      }
    else
      {
        Cerr << "Error : the radiative properties of the incompressible fluid have not" << finl;
        Cerr << "been defined while a semi transparent radiation model is used." << finl;
        Cerr << "The fields kappa and indice which respectively define the absoption coefficient" << finl;
        Cerr << "and the refraction index must be added to the luid properties." << finl;
        Process::exit();
      }
  else
    {
      Cerr << "Error while reading the Radiation equation." << finl;
      Cerr << "Your fluid is of type " << fluide().que_suis_je() << finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered " << finl;
      Cerr << "with the semi transparent radiation model." << finl;
      Process::exit();
    }

  Equation_base::completer();
}

/*! @brief Associe un milieu physique a l'equation
 *
 * @param (Milieu_base& un_milieu) le milieu physique a associer a l'equation
 */
void Equation_rayonnement_base::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_base, un_milieu))
    if (fluide().has_kappa())
      {
        const Fluide_base& un_fluide = ref_cast(Fluide_base, un_milieu);
        associer_fluide(un_fluide);
        le_fluide = un_fluide;
      }
    else
      {
        Cerr << "Error : the radiative properties of the incompressible fluid have not" << finl;
        Cerr << "been defined while a semi transparent radiation model is used." << finl;
        Cerr << "The fields kappa and indice which respectively define the absoption coefficient" << finl;
        Cerr << "and the refraction index must be added to the luid properties." << finl;
        Process::exit();
      }
  else
    {
      Cerr << "Error for the method Equation_rayonnement_base::associer_milieu_base" << finl;
      Cerr << "Your fluid is of type " << un_milieu.que_suis_je() << finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered " << finl;
      Cerr << "with the semi transparent radiation model." << finl;
      Process::exit();
    }
}

/*! @brief Renvoie le milieu physique de l'equation (le Fluide_base upcaste en Milieu_base)
 *
 * @return (Milieu_base&) le Fluide_base de l'equation upcaste en Milieu_base
 */
const Milieu_base& Equation_rayonnement_base::milieu() const
{
  if (!le_fluide.non_nul())
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

/*! @brief Renvoie le milieu physique de l'equation (le Fluide_base upcaste en Milieu_base)
 *
 *     (version const)
 *
 * @return (Milieu_base&) le Fluide_base de l'equation upcaste en Milieu_base
 */
Milieu_base& Equation_rayonnement_base::milieu()
{
  if (!le_fluide.non_nul())
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}

/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      exit si i>0
 *     (version const)
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 1 operateur
 */
const Operateur& Equation_rayonnement_base::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    default:
      Cerr << "Error for Equation_rayonnement_base::operateur(int i)" << finl;
      Cerr << "Equation_rayonnement_base has " << nombre_d_operateurs() << " operators " << finl;
      Cerr << "and you are trying to access the " << i << " th one." << finl;
      Process::exit();
    }
  return terme_diffusif;
}

/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      exit si i>0
 *     (version const)
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 1 operateur
 */
Operateur& Equation_rayonnement_base::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    default:
      Cerr << "Error for Equation_rayonnement_base::operateur(int i)" << finl;
      Cerr << "Equation_rayonnement_base has " << nombre_d_operateurs() << " operators " << finl;
      Cerr << "and you are trying to access the " << i << " th one." << finl;
      Process::exit();
    }
  return terme_diffusif;
}

void Equation_rayonnement_base::get_noms_champs_postraitables(Noms& noms, Option opt) const
{
  if (opt == DESCRIPTION)
    Cerr << que_suis_je() << " : " << champs_compris_.liste_noms_compris() << finl;
  else
    noms.add(champs_compris_.liste_noms_compris());
}

void Equation_rayonnement_base::discretiser()
{
  // Discretisation de l'equation de rayonnement
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  Cerr << "Radiation equation discretisation" << finl;
  dis.discretiser_champ("temperature", domaine_dis(), "irradiance", "w/m2", 1, 1 /* une case */, schema_temps().temps_courant(), irradiance_);
  champs_compris_.ajoute_champ(irradiance_);

  Equation_base::discretiser();
}

/*! @brief Renvoie la discretisation associee a l'equation.
 *
 * @return (Discretisation_base&) a discretisation associee a l'equation
 * @throws pas de probleme associe
 */
const Discretisation_base& Equation_rayonnement_base::discretisation() const
{
  return le_modele->discretisation();
}

void Equation_rayonnement_base::associer_pb_base(const Probleme_base& pb)
{
  Equation_base::associer_pb_base(pb);
  le_modele = ref_cast(Modele_rayo_semi_transp, pb);
  associer_sch_tps_base(pb.schema_temps());
}

void Equation_rayonnement_base::Mat_Morse_to_Mat_Bloc(Matrice& matrice_tmp)
{
  int n1 = nb_colonnes_tot();
  int n2 = nb_colonnes();

  Matrice_Bloc& matrice = ref_cast(Matrice_Bloc, matrice_tmp.valeur());
  Matrice_Morse& MBrr = ref_cast(Matrice_Morse, matrice.get_bloc(0, 0).valeur());
  Matrice_Morse& MBrv = ref_cast(Matrice_Morse, matrice.get_bloc(0, 1).valeur());

  IntVect& tab1RR = MBrr.get_set_tab1();
  IntVect& tab2RR = MBrr.get_set_tab2();
  DoubleVect& coeffRR = MBrr.get_set_coeff();
  IntVect& tab1RV = MBrv.get_set_tab1();
  IntVect& tab2RV = MBrv.get_set_tab2();
  DoubleVect& coeffRV = MBrv.get_set_coeff();

  DoubleTab ligne_tmp(n1);
  for (int i = 0; i < n2; i++)
    {
      int k;
      // On recopie le premier bloc de la matrice dans un tableau :
      //      ligne_tmp = 0;
      for (k = la_matrice.get_tab1()(i) - 1; k < la_matrice.get_tab1()(i + 1) - 1; k++)
        ligne_tmp(la_matrice.get_tab2()(k) - 1) = la_matrice.get_coeff()(k);

      // On complete la partie reelle de la matrice
      for (k = tab1RR(i) - 1; k < tab1RR(i + 1) - 1; k++)
        coeffRR[k] = ligne_tmp(tab2RR[k] - 1);

      // On complete la partie virtuelle
      for (k = tab1RV(i) - 1; k < tab1RV(i + 1) - 1; k++)
        coeffRV[k] = ligne_tmp(n2 + tab2RV[k] - 1);
    }
}

void Equation_rayonnement_base::dimensionner_Mat_Bloc_Morse_Sym(Matrice& matrice_tmp)
{
  int n1 = nb_colonnes_tot();
  int n2 = nb_colonnes();
  int iligne;
  const IntVect& tab1 = la_matrice.get_set_tab1();
  const IntVect& tab2 = la_matrice.get_set_tab2();

  matrice_tmp.typer("Matrice_Bloc");
  Matrice_Bloc& matrice = ref_cast(Matrice_Bloc, matrice_tmp.valeur());
  matrice.dimensionner(1, 2);
  matrice.get_bloc(0, 0).typer("Matrice_Morse_Sym");
  matrice.get_bloc(0, 1).typer("Matrice_Morse");

  Matrice_Morse_Sym& MBrr = ref_cast(Matrice_Morse_Sym, matrice.get_bloc(0, 0).valeur());
  Matrice_Morse& MBrv = ref_cast(Matrice_Morse, matrice.get_bloc(0, 1).valeur());
  MBrr.dimensionner(n2, 0);
  MBrv.dimensionner(n2, 0);

  IntVect& tab1RR = MBrr.get_set_tab1();
  IntVect& tab2RR = MBrr.get_set_tab2();
  IntVect& tab1RV = MBrv.get_set_tab1();
  IntVect& tab2RV = MBrv.get_set_tab2();

  IntVect compteur_MBrr(n2);
  IntVect compteur_MBrv(n2);
  compteur_MBrr = 0;
  compteur_MBrv = 0;

  // On parcours les lignes de la_matrice pour compter les elements
  // non nuls de chaque ligne
  int jcolonne;
  for (iligne = 0; iligne < n2; iligne++)
    {
      int k;
      for (k = tab1(iligne) - 1; k < tab1(iligne + 1) - 1; k++)
        {
          jcolonne = tab2(k) - 1;
          if (jcolonne < n2)
            {
              // l'element correspondant est dans la partie RR de la_matrice
              if ((jcolonne >= iligne) && (jcolonne < n2))
                {
                  // l'element correspondant est  situe au dessus de la diagonale de la_matrice
                  compteur_MBrr(iligne)++;
                }
            }
          else
            {
              // l'element correspondant est dans la partie RV de la_matrice
              compteur_MBrv(iligne)++;
            }
        }
    }

  // On remplie tab1RR et tab1RV
  tab1RR(0) = 1;
  tab1RV(0) = 1;
  for (int i = 0; i < n2; i++)
    {
      tab1RR(i + 1) = compteur_MBrr(i) + tab1RR(i);
      tab1RV(i + 1) = compteur_MBrv(i) + tab1RV(i);
    }
  // On dimensionne tab2RR et tab2RV
  MBrr.dimensionner(n2, tab1RR(n2) - 1);
  MBrv.dimensionner(n2, n1 - n2, tab1RV(n2) - 1);

  // On remplit tab2RR et tab2RV
  int compteurRR, compteurRV;
  for (iligne = 0; iligne < n2; iligne++)
    {
      int k;
      compteurRR = tab1RR(iligne) - 1;
      compteurRV = tab1RV(iligne) - 1;
      for (k = tab1(iligne) - 1; k < tab1(iligne + 1) - 1; k++)
        {
          jcolonne = tab2(k) - 1;
          if (jcolonne < n2)
            {
              // l'element correspondant est dans la partie RR de la_matrice
              if ((jcolonne >= iligne) && (jcolonne < n2))
                {
                  // l'element correspondant est  situe au dessus de la diagonale de la_matrice
                  tab2RR(compteurRR) = tab2(k);
                  compteurRR++;
                }
            }
          else
            {
              // l'element correspondant est dans la partie RV de la_matrice
              tab2RV(compteurRV) = tab2(k) - n2;
              compteurRV++;
            }
        }
    }
}
