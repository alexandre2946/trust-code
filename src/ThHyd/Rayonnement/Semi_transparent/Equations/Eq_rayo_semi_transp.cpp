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

#include <Eq_rayo_semi_transp.h>
#include <Pb_rayo_semi_transp.h>
#include <Matrice_Morse_Sym.h>
#include <Champ_Uniforme.h>
#include <Matrice_Bloc.h>
#include <Discret_Thyd.h>
#include <Fluide_base.h>
#include <EChaine.h>
#include <Param.h>

Implemente_instanciable(Eq_rayo_semi_transp, "Eq_rayo_semi_transp", Equation_base);

bool Eq_rayo_semi_transp::initTimeStep(double dt)
{
  schema_temps().set_dt() = dt;
  return Equation_base::initTimeStep(dt);
}

void Eq_rayo_semi_transp::resoudre(double temps)
{
  rayo_solv_->resoudre(temps);
}

bool Eq_rayo_semi_transp::resoudre()
{
  rayo_solv_->resoudre(schema_temps().temps_courant() + schema_temps().pas_de_temps());
  return true;
}

Sortie& Eq_rayo_semi_transp::printOn(Sortie& s) const { return s << que_suis_je() << finl; }

Entree& Eq_rayo_semi_transp::readOn(Entree& is) { return Equation_base::readOn(is); }

void Eq_rayo_semi_transp::set_param(Param& param) const
{
  param.ajouter_non_std("conditions_limites|boundary_conditions", (this), Param::REQUIRED);
  param.ajouter_non_std("solveur", (this), Param::REQUIRED);
}

int Eq_rayo_semi_transp::lire_motcle_non_standard(const Motcle& mot, Entree& is)
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
      solveur_.typer(nom_solveur);
      is >> solveur_.valeur();
      solveur_.nommer("solveur_irradiance");
    }
  else
    retval = -1;

  return retval;
}

void Eq_rayo_semi_transp::completer()
{
  const Domaine_dis_base& dom_dis = domaine_dis();
  for (int i = 0; i < dom_dis.nb_front_Cl(); i++)
    {
      const Frontiere_dis_base& la_fr_dis = dom_dis.frontiere_dis(i);
      le_dom_Cl_dis->les_conditions_limites(i)->associer_fr_dis_base(la_fr_dis);
    }

  // typage de l'operateur de diffusion
  Cerr << "Reading and typing of the diffusion operator of equation " << que_suis_je() << finl;

  if (sub_type(Fluide_base, fluide()))
    if (fluide().is_rayo_semi_transp())
      {
        if (fluide().is_longueur_rayo_discretised())
          terme_diffusif_.associer_diffusivite(fluide().longueur_rayo());
        else
          {
            Cerr << "Error in Eq_rayo_semi_transp::completer." << finl;
            Cerr << "You may not have discretized the problem of type Pb_Couple_rayo_semi_transp." << finl;
            Process::exit();
          }
      }
    else
      {
        Cerr << "Error : the radiative properties of the incompressible fluid have not" << finl;
        Cerr << "been defined while a semi transparent radiation problem is used." << finl;
        Cerr << "The fields kappa and indice which respectively define the absoption coefficient" << finl;
        Cerr << "and the refraction index must be added to the fluid properties." << finl;
        Process::exit();
      }
  else
    {
      Cerr << "Error while reading the Radiation equation. Your fluid is of type " << fluide().que_suis_je() << finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered with the semi transparent radiation model." << finl;
      Process::exit();
    }

  EChaine diff("{ }");
  diff >> terme_diffusif_;
  terme_diffusif_->associer_diffusivite(fluide().longueur_rayo());
  terme_diffusif_.completer();
  terme_diffusif_->dimensionner(la_matrice_);

  Equation_base::completer();

  // On assemble la matrice une fois pour toute au debut du calcul
  // XXX Attention, ceci n'est valable que si kappa est constant au cours du temps
  rayo_solv_->assembler_matrice();
}

/*! @brief Associe un milieu physique a l'equation
 *
 * @param (Milieu_base& un_milieu) le milieu physique a associer a l'equation
 */
void Eq_rayo_semi_transp::associer_milieu_base(const Milieu_base& un_milieu)
{
  if (sub_type(Fluide_base, un_milieu))
    if (fluide().is_rayo_semi_transp())
      {
        const Fluide_base& un_fluide = ref_cast(Fluide_base, un_milieu);
        associer_fluide(un_fluide);
        le_fluide_ = un_fluide;
      }
    else
      {
        Cerr << "Error : the radiative properties of the incompressible fluid have not" << finl;
        Cerr << "been defined while a semi transparent radiation problem is used." << finl;
        Cerr << "The fields kappa and indice which respectively define the absoption coefficient" << finl;
        Cerr << "and the refraction index must be added to the fluid properties." << finl;
        Process::exit();
      }
  else
    {
      Cerr << "Error while reading the Radiation equation. Your fluid is of type " << fluide().que_suis_je() << finl;
      Cerr << "Currently only fluid of type Fluide_base can be considered with the semi transparent radiation model." << finl;
      Process::exit();
    }
}

/*! @brief Renvoie le milieu physique de l'equation (le Fluide_base upcaste en Milieu_base)
 *
 * @return (Milieu_base&) le Fluide_base de l'equation upcaste en Milieu_base
 */
const Milieu_base& Eq_rayo_semi_transp::milieu() const
{
  if (!le_fluide_)
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
}

/*! @brief Renvoie le milieu physique de l'equation (le Fluide_base upcaste en Milieu_base)
 *
 *     (version const)
 *
 * @return (Milieu_base&) le Fluide_base de l'equation upcaste en Milieu_base
 */
Milieu_base& Eq_rayo_semi_transp::milieu()
{
  if (!le_fluide_)
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide_.valeur();
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
const Operateur& Eq_rayo_semi_transp::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif_;
    default:
      Cerr << "Error for Eq_rayo_semi_transp::operateur(int i)" << finl;
      Cerr << "Eq_rayo_semi_transp has " << nombre_d_operateurs() << " operators " << finl;
      Cerr << "and you are trying to access the " << i << " th one." << finl;
      Process::exit();
    }
  return terme_diffusif_;
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
Operateur& Eq_rayo_semi_transp::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif_;
    default:
      Cerr << "Error for Eq_rayo_semi_transp::operateur(int i)" << finl;
      Cerr << "Eq_rayo_semi_transp has " << nombre_d_operateurs() << " operators " << finl;
      Cerr << "and you are trying to access the " << i << " th one." << finl;
      Process::exit();
    }
  return terme_diffusif_;
}

void Eq_rayo_semi_transp::get_noms_champs_postraitables(Noms& noms, Option opt) const
{
  if (opt == DESCRIPTION)
    Cerr << que_suis_je() << " : " << champs_compris_.liste_noms_compris() << finl;
  else
    noms.add(champs_compris_.liste_noms_compris());
}

void Eq_rayo_semi_transp::discretiser()
{
  // Discretisation de l'equation de rayonnement
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  Cerr << "Radiation equation discretisation" << finl;
  dis.discretiser_champ("temperature", domaine_dis(), "irradiance", "w/m2", 1, 1 /* une case */, schema_temps().temps_courant(), irradiance_);
  champs_compris_.ajoute_champ(irradiance_);

  Equation_base::discretiser();

  // typing Rayo_semi_transp_solver_base
  Nom discr = dis.que_suis_je(), type = "Rayo_semi_transp_solver_";
  if (discr == "VEFPreP1B") discr = "VEF";
  type += discr;
  rayo_solv_.typer(type);
  rayo_solv_->associer_equation_rayo(*this);
}

/*! @brief Renvoie la discretisation associee a l'equation.
 *
 * @return (Discretisation_base&) a discretisation associee a l'equation
 * @throws pas de probleme associe
 */
const Discretisation_base& Eq_rayo_semi_transp::discretisation() const
{
  return pb_rayo_semi_transp_->discretisation();
}

void Eq_rayo_semi_transp::associer_pb_base(const Probleme_base& pb)
{
  Equation_base::associer_pb_base(pb);
  pb_rayo_semi_transp_ = ref_cast(Pb_rayo_semi_transp, pb);
  associer_sch_tps_base(pb.schema_temps());
}

void Eq_rayo_semi_transp::Mat_Morse_to_Mat_Bloc(Matrice& matrice_tmp)
{
  const int n1 = rayo_solv_->nb_colonnes_tot();
  const int n2 = rayo_solv_->nb_colonnes();

  Matrice_Bloc& matrice = ref_cast(Matrice_Bloc, matrice_tmp.valeur());
  Matrice_Morse& MBrr = ref_cast(Matrice_Morse, matrice.get_bloc(0, 0).valeur());
  Matrice_Morse& MBrv = ref_cast(Matrice_Morse, matrice.get_bloc(0, 1).valeur());

  auto& tab1RR = MBrr.get_set_tab1();
  auto& tab2RR = MBrr.get_set_tab2();
  auto& coeffRR = MBrr.get_set_coeff();
  auto& tab1RV = MBrv.get_set_tab1();
  auto& tab2RV = MBrv.get_set_tab2();
  auto& coeffRV = MBrv.get_set_coeff();

  DoubleTab ligne_tmp(n1);
  for (int i = 0; i < n2; i++)
    {
      int k;
      // On recopie le premier bloc de la matrice dans un tableau :
      //      ligne_tmp = 0;
      for (k = la_matrice_.get_tab1()(i) - 1; k < la_matrice_.get_tab1()(i + 1) - 1; k++)
        ligne_tmp(la_matrice_.get_tab2()(k) - 1) = la_matrice_.get_coeff()(k);

      // On complete la partie reelle de la matrice
      for (k = tab1RR(i) - 1; k < tab1RR(i + 1) - 1; k++)
        coeffRR[k] = ligne_tmp(tab2RR[k] - 1);

      // On complete la partie virtuelle
      for (k = tab1RV(i) - 1; k < tab1RV(i + 1) - 1; k++)
        coeffRV[k] = ligne_tmp(n2 + tab2RV[k] - 1);
    }
}

void Eq_rayo_semi_transp::dimensionner_Mat_Bloc_Morse_Sym(Matrice& matrice_tmp)
{
  const int n1 = rayo_solv_->nb_colonnes_tot();
  const int n2 = rayo_solv_->nb_colonnes();

  int iligne;
  const auto& tab1 = la_matrice_.get_set_tab1();
  const auto& tab2 = la_matrice_.get_set_tab2();

  matrice_tmp.typer("Matrice_Bloc");
  Matrice_Bloc& matrice = ref_cast(Matrice_Bloc, matrice_tmp.valeur());
  matrice.dimensionner(1, 2);
  matrice.get_bloc(0, 0).typer("Matrice_Morse_Sym");
  matrice.get_bloc(0, 1).typer("Matrice_Morse");

  Matrice_Morse_Sym& MBrr = ref_cast(Matrice_Morse_Sym, matrice.get_bloc(0, 0).valeur());
  Matrice_Morse& MBrv = ref_cast(Matrice_Morse, matrice.get_bloc(0, 1).valeur());
  MBrr.dimensionner(n2, 0);
  MBrv.dimensionner(n2, 0);

  auto& tab1RR = MBrr.get_set_tab1();
  auto& tab2RR = MBrr.get_set_tab2();
  auto& tab1RV = MBrv.get_set_tab1();
  auto& tab2RV = MBrv.get_set_tab2();

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
