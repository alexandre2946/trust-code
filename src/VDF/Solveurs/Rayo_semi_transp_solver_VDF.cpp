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

#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>
#include <Echange_externe_impose_rayo_semi_transp.h>
#include <Echange_global_impose_rayo_semi_transp.h>
#include <Echange_contact_rayo_semi_transp_VDF.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Neumann_paroi_rayo_semi_transp_VDF.h>
#include <Rayo_semi_transp_solver_VDF.h>
#include <Champ_front_uniforme.h>
#include <Eq_rayo_semi_transp.h>
#include <Pb_rayo_semi_transp.h>
#include <Flux_radiatif_VDF.h>
#include <Champ_Uniforme.h>
#include <Domaine_VDF.h>
#include <Domaine_VDF.h>
#include <Fluide_base.h>
#include <Symetrie.h>
#include <Debog.h>

Implemente_instanciable(Rayo_semi_transp_solver_VDF, "Rayo_semi_transp_solver_VDF", Rayo_semi_transp_solver_base);

Sortie& Rayo_semi_transp_solver_VDF::printOn(Sortie& s) const { return s << que_suis_je() << finl; }
Entree& Rayo_semi_transp_solver_VDF::readOn(Entree& is) { return is; }

int Rayo_semi_transp_solver_VDF::nb_colonnes_tot()
{
  Eq_rayo_semi_transp& eq_rayo = eq_rayo_semi_transp_.valeur();
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, eq_rayo.domaine_dis());
  return domaine_VF.domaine().nb_elem_tot();
}

int Rayo_semi_transp_solver_VDF::nb_colonnes()
{
  Eq_rayo_semi_transp& eq_rayo = eq_rayo_semi_transp_.valeur();
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, eq_rayo.domaine_dis());
  return domaine_VF.domaine().nb_elem();
}

/*! @brief Modifies the matrix to account for the presence of radiating faces near boundary elements.
 *
 * @return the modified input stream
 */
void Rayo_semi_transp_solver_VDF::modifier_matrice()
{
  Eq_rayo_semi_transp& eq_rayo = eq_rayo_semi_transp_.valeur();
  Matrice_Morse& matrice = eq_rayo.matrice_rayo();
  Conds_lim& les_cl = eq_rayo.domaine_Cl_dis().les_conditions_limites();

  const auto& fluide = eq_rayo.fluide();
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, eq_rayo.domaine_dis());
  const IntTab& face_voisins = zvdf.face_voisins();
  const DoubleVect& face_surfaces = zvdf.face_surfaces();

  // Loop over boundary conditions associated with the equation
  for (int num_cl = 0; num_cl < les_cl.size(); num_cl++)
    {
      Cond_lim& la_cl = eq_rayo.domaine_Cl_dis().les_conditions_limites(num_cl);
      if (sub_type(Flux_radiatif_VDF, la_cl.valeur()))
        {
          Flux_radiatif_VDF& cl_radiatif = ref_cast(Flux_radiatif_VDF, la_cl.valeur());
          const DoubleTab& epsilon = cl_radiatif.emissivite().valeurs();
          const DoubleTab& long_rayo = fluide.longueur_rayo().valeurs();
          const DoubleTab& kappa = fluide.kappa().valeurs();
          const double A = cl_radiatif.A();

          if (sub_type(Front_VF, la_cl->frontiere_dis()))
            {
              const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
              const int ndeb = le_bord.num_premiere_face();
              const int nfin = ndeb + le_bord.nb_faces();
              assert(cl_radiatif.emissivite().nb_comp() == 1);
              assert(fluide.longueur_rayo().nb_comp() == 1);
              assert(fluide.kappa().nb_comp() == 1);

              for (int face = ndeb; face < nfin; face++)
                {
                  int elem = face_voisins(face, 0);

                  if (elem == -1)
                    elem = face_voisins(face, 1);


                  double eF = zvdf.dist_norm_bord(face);
                  double k = -123., l_r = -123.;

                  if (sub_type(Champ_Uniforme, fluide.kappa()))
                    {
                      k = kappa(0, 0);
                      l_r = long_rayo(0, 0);
                    }
                  else
                    {
                      k = kappa(elem, 0);
                      l_r = long_rayo(elem, 0);
                    }

                  double epsi = -123.;

                  if (sub_type(Champ_front_uniforme, cl_radiatif.emissivite()))
                    epsi = epsilon(0, 0);
                  else
                    epsi = epsilon(face - ndeb, 0);

                  const double numer_coeff = l_r * face_surfaces(face);
                  double denum_coeff = 3 * k * epsi;

                  denum_coeff = 1 / denum_coeff;
                  denum_coeff *= A * (2 - epsi);
                  denum_coeff = denum_coeff + eF;

                  const double coeff = numer_coeff / denum_coeff;

                  // Add this coefficient to the diagonal of the discretization matrix
                  if (epsi < DMINFLOAT) { /* Do nothing */}
                  else
                    matrice(elem, elem) += coeff;
                }
            }
          else
            {
              Cerr << "Error in Rayo_semi_transp_solver_VDF::modifier_matrice()" << finl;
              Cerr << "the boundary associated with la_cl does not derive from Front_VF" << finl;
              Process::exit();
            }
        }
      else if (sub_type(Symetrie, la_cl.valeur()))
        {
          /* Do nothing */
        }
      else
        Process::exit("The boundary condition used is not recognized for the radiation equation!");
    }
}

void Rayo_semi_transp_solver_VDF::assembler_matrice()
{
  Eq_rayo_semi_transp& eq_rayo = eq_rayo_semi_transp_.valeur();
  Matrice_Morse& matrice = eq_rayo.matrice_rayo();
  Operateur_Diff& terme_diffusif = eq_rayo.terme_diffusif_rayo();

  const auto& fluide = eq_rayo.fluide();
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, eq_rayo.domaine_dis());
  const int nb_elem_tot = domaine_VF.nb_elem_tot();

  const DoubleTab& irradi = eq_rayo.inconnue().valeurs();

  matrice.clean();

  // Account for the div((1/3K)grad(irradiance)) part in the discretization matrix.
  terme_diffusif->contribuer_a_avec(irradi, matrice);

  // Modify the matrix to account for the right-hand side term K*irradiance
  const DoubleTab& kappa = fluide.kappa().valeurs();

  // verify
  if (matrice.ordre() != nb_elem_tot)
    Process::exit();

  Cerr << "Matrix order OK" << finl;
  assert(fluide.kappa().nb_comp() == 1);

  double k = -123.;
  for (int i = 0; i < matrice.ordre(); i++)
    {
      if (sub_type(Champ_Uniforme, fluide.kappa()))
        k = kappa(0, 0);
      else
        k = kappa(i, 0);

      const double vol = domaine_VF.volumes(i);

      matrice(i, i) = matrice(i, i) + k * vol;
    }

  // Modify the matrix to account for the effect of radiating walls on boundary elements
  modifier_matrice();
}

void Rayo_semi_transp_solver_VDF::resoudre(double temps)
{
  Eq_rayo_semi_transp& eq_rayo = eq_rayo_semi_transp_.valeur();
  Operateur_Diff& terme_diffusif = eq_rayo.terme_diffusif_rayo();
  Matrice_Morse& matrice = eq_rayo.matrice_rayo();
  SolveurSys& solveur = eq_rayo.solveur_rayo();

  const auto& fluide = eq_rayo.fluide();
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF, eq_rayo.domaine_dis());
  const int nb_elem = domaine_VF.nb_elem();
  const DoubleTab& kappa = fluide.kappa().valeurs();

  /*
   * Note: the matrix assembly is performed once at the start of the
   * computation, which is only valid when kappa does not depend on time.
   */

  //compute the right-hand side
  DoubleTrav secmem(eq_rayo.inconnue().valeurs());
  secmem = 0.;

  // retrieve T from the fluid problem ....
  Probleme_base& pb_fluide = eq_rayo.pb_rayo_semi_transp().probleme_fluide();

  assert(pb_fluide.equation(1).inconnue().le_nom() == "temperature");
  const DoubleTab& temper = pb_fluide.equation(1).inconnue().valeurs();
  const DoubleTab& indice = fluide.indice().valeurs();
  const double sigma = eq_rayo.pb_rayo_semi_transp().valeur_sigma();
  assert(fluide.indice().nb_comp());
  assert(fluide.kappa().nb_comp() == 1);

  double n = -123., k = -123.;
  for (int elem = 0; elem < nb_elem; elem++)
    {
      if (sub_type(Champ_Uniforme, fluide.indice()))
        n = indice(0, 0);
      else
        n = indice(elem, 0);

      if (sub_type(Champ_Uniforme, fluide.kappa()))
        k = kappa(0, 0);
      else
        k = kappa(elem, 0);

      const double vol = domaine_VF.volumes(elem);
      const double T = temper(elem);
      secmem(elem) += +4 * n * n * sigma * pow(T, 4) * k * vol;
    }

  // Update the fields associated with boundary conditions
  // before evaluating their contribution in the discretization matrix
  evaluer_cl_rayonnement(temps);
  terme_diffusif->contribuer_au_second_membre(secmem);

  if (solveur->que_suis_je() == "Solv_GCP")
    if (sub_type(Champ_Uniforme, fluide.kappa()))
      {
        Matrice matrice_tmp;
        eq_rayo.dimensionner_Mat_Bloc_Morse_Sym(matrice_tmp);
        eq_rayo.Mat_Morse_to_Mat_Bloc(matrice_tmp);
        solveur.resoudre_systeme(matrice_tmp.valeur(), secmem, eq_rayo.inconnue().valeurs());
      }
    else
      {
        Cerr << "Error in Rayo_semi_transp_solver_VDF::resoudre()! Cannot solve the equation" << finl;
        Cerr << "of semi-transparent radiation with the solver: " << solveur->que_suis_je() << " because kappa is not constant, therefore the matrix is not symmetric" << finl;
        Process::exit();
      }
  else if (solveur->que_suis_je() == "Solv_Gmres")
    if (Process::nproc() == 1)
      solveur.resoudre_systeme(matrice, secmem, eq_rayo.inconnue().valeurs());
    else
      {
        Cerr << "Error in Rayo_semi_transp_solver_VDF::resoudre()! Cannot solve a problem" << finl;
        Cerr << "of semi-transparent radiation in parallel using the Solv_Gmres solver. If you treat" << finl;
        Cerr << "a problem with constant kappa, you can bypass this limitation by using the GCP solver with GCP preconditioning" << finl;
        Process::exit();
      }
  else
    {
      Cerr << "Error in Rayo_semi_transp_solver_VDF::resoudre()! Cannot use the solver: " << solveur->que_suis_je() << finl;
      Cerr << "to solve the radiation equation in a semi-transparent radiation problem" << finl;
      Process::exit();
    }

  Debog::verifier("Rayo_semi_transp_solver_VEF::resoudre irradiance ", eq_rayo.inconnue().valeurs());
}

void Rayo_semi_transp_solver_VDF::evaluer_cl_rayonnement(double temps)
{
  Eq_rayo_semi_transp& eq_rayo = eq_rayo_semi_transp_.valeur();
  const auto& fluide = eq_rayo.fluide();

  // Search for boundary conditions associated with the temperature equation
  Conds_lim& les_cl_rayo = eq_rayo.domaine_Cl_dis().les_conditions_limites();
  Equation_base& eq_temp = eq_rayo.pb_rayo_semi_transp().probleme_fluide().equation(1);
  assert(eq_temp.inconnue().le_nom() == "temperature");

  // Loop over boundary conditions of the radiation equation
  Conds_lim& les_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites();
  for (int num_cl_rayo = 0; num_cl_rayo < les_cl_rayo.size(); num_cl_rayo++)
    {
      Cond_lim& la_cl_rayo = eq_rayo.domaine_Cl_dis().les_conditions_limites(num_cl_rayo);
      if (sub_type(Flux_radiatif_VDF, la_cl_rayo.valeur()))
        {
          Flux_radiatif_VDF& la_cl_rayon = ref_cast(Flux_radiatif_VDF, la_cl_rayo.valeur());
          // Search for boundary temperatures for this boundary
          Nom nom_cl_rayo = la_cl_rayo->frontiere_dis().le_nom();

          OBS_PTR(Champ_front_base) Tb;
          int test_remplissage_Tb = 0;
          for (int num_cl_temp = 0; num_cl_temp < les_cl_temp.size(); num_cl_temp++)
            {
              Cond_lim& la_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites(num_cl_temp);
              Nom nom_cl_temp = la_cl_temp->frontiere_dis().le_nom();
              if (nom_cl_temp == nom_cl_rayo)
                {
                  if (sub_type(Neumann_paroi_rayo_semi_transp_VDF, la_cl_temp.valeur()))
                    {
                      Neumann_paroi_rayo_semi_transp_VDF& la_cl_temper = ref_cast(Neumann_paroi_rayo_semi_transp_VDF, la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Echange_contact_rayo_semi_transp_VDF, la_cl_temp.valeur()))
                    {
                      Echange_contact_rayo_semi_transp_VDF& la_cl_temper = ref_cast(Echange_contact_rayo_semi_transp_VDF, la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Echange_externe_impose_rayo_semi_transp, la_cl_temp.valeur()))
                    {
                      Echange_externe_impose_rayo_semi_transp& la_cl_temper = ref_cast(Echange_externe_impose_rayo_semi_transp, la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Frontiere_ouverte_temperature_imposee_rayo_semi_transp, la_cl_temp.valeur()))
                    {
                      Frontiere_ouverte_temperature_imposee_rayo_semi_transp& la_cl_temper = ref_cast(Frontiere_ouverte_temperature_imposee_rayo_semi_transp, la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Frontiere_ouverte_rayo_semi_transp, la_cl_temp.valeur()))
                    {
                      Frontiere_ouverte_rayo_semi_transp& la_cl_temper = ref_cast(Frontiere_ouverte_rayo_semi_transp, la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else if (sub_type(Echange_global_impose_rayo_semi_transp, la_cl_temp.valeur()))
                    {
                      Echange_global_impose_rayo_semi_transp& la_cl_temper = ref_cast(Echange_global_impose_rayo_semi_transp, la_cl_temp.valeur());
                      test_remplissage_Tb = 1;
                      la_cl_temper.calculer_temperature_bord(temps);
                      Tb = la_cl_temper.temperature_bord();
                    }
                  else
                    {
                      Cerr << "Error in Rayo_semi_transp_solver_VDF::evaluer_cl_rayonnement! The case of a thermal BC " << la_cl_temp->que_suis_je() << " is not implemented" << finl;
                      Process::exit();
                    }
                }
            }

          // The boundary temperature array was not filled!!!!
          if (test_remplissage_Tb == 0)
            Cerr << "Rayo_semi_transp_solver_VDF::evaluer_cl_rayonnement -- The boundary temperature array was not filled!!!!" << finl;

          const Domaine_VF& zvf = ref_cast(Domaine_VF, eq_rayo.domaine_dis());
          la_cl_rayon.evaluer_cl_rayonnement(Tb.valeur(), fluide.kappa(), fluide.longueur_rayo(), fluide.indice(), zvf, eq_rayo.pb_rayo_semi_transp().valeur_sigma(), temps);
        }
      else if (sub_type(Symetrie, la_cl_rayo.valeur()))
        {
          /* Do nothing */
        }
      else
        {
          Cerr << "The boundary condition " << la_cl_rayo->que_suis_je() << " is not recognized for the radiation equation!" << finl;
          Process::exit();
        }
    }
}
