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
#include <Champ_front_uniforme.h>
#include <Eq_rayo_semi_transp.h>
#include <Pb_rayo_semi_transp.h>
#include <Flux_radiatif_VDF.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <Domaine_VDF.h>
#include <Debog.h>

Implemente_instanciable(Flux_radiatif_VDF, "Flux_radiatif_VDF", Flux_radiatif_base);

Sortie& Flux_radiatif_VDF::printOn(Sortie& s) const { return s << que_suis_je() << finl; }
Entree& Flux_radiatif_VDF::readOn(Entree& s) { return Flux_radiatif_base::readOn(s); }

void Flux_radiatif_VDF::evaluer_cl_rayonnement(Champ_front_base& Tb, const Champ_Don_base& coeff_abs, const Champ_Don_base& longueur_rayo,
                                               const Champ_Don_base& indice, const Domaine_VF& zvf,
                                               const double sigma, double temps)
{
  const DoubleTab& l_rayo = longueur_rayo.valeurs();
  const DoubleTab& n = indice.valeurs();
  const DoubleTab& epsilon = emissivite().valeurs();
  const DoubleTab& kappa = coeff_abs.valeurs();

  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, zvf);
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());
  const IntTab& face_voisins = zvdf.face_voisins();

  // Dimension le_champ_front
  assert(le_champ_front->nb_comp() == 1);
  DoubleTab& tab = le_champ_front->valeurs_au_temps(temps);
  tab.resize(le_bord.nb_faces(), le_champ_front->nb_comp());

  assert(indice.nb_comp() == 1);
  assert(longueur_rayo.nb_comp() == 1);
  assert(coeff_abs.nb_comp() == 1);
  assert(Tb.nb_comp() == 1);
  assert(emissivite().nb_comp() == 1);

  const int ndeb = le_bord.num_premiere_face();
  const int nfin = ndeb + le_bord.nb_faces();
  double nn = -123., k = -123., l_r = -123.;
  double epsi = -123., T = -123.;

  // Loop over the faces of le_bord
  for (int face = ndeb; face < nfin; face++)
    {
      int elem = face_voisins(face, 0);

      const double eF = zvdf.dist_norm_bord(face);

      if (sub_type(Champ_Uniforme, indice))
        nn = n(0, 0);
      else
        nn = n(elem, 0);

      if (sub_type(Champ_Uniforme, longueur_rayo))
        l_r = l_rayo(0, 0);
      else
        l_r = l_rayo(elem, 0);

      if (sub_type(Champ_Uniforme, coeff_abs))
        k = kappa(0, 0);
      else
        k = kappa(elem, 0);

      // Determination of the wall temperature as a function of the face considered
      if (sub_type(Champ_front_uniforme, Tb))
        T = Tb.valeurs()(0, 0);
      else
        T = Tb.valeurs_au_temps(temps)(face - ndeb, 0);

      // Determination of the wall emissivity as a function of the face considered
      if (sub_type(Champ_front_uniforme, emissivite()))
        epsi = epsilon(0, 0);
      else
        epsi = epsilon(face - ndeb, 0);

      // Filling the boundary condition for the radiation equation
      double numer_coeff = l_r;
      numer_coeff *= 4 * nn * nn * sigma * pow(T, 4);

      double denum_coeff = 3 * k * epsi;
      denum_coeff = 1 / denum_coeff;
      denum_coeff *= A_ * (2 - epsi);
      denum_coeff = denum_coeff + eF;

      if (epsi < DMINFLOAT)
        tab(face - ndeb, 0) = 0.;
      else
        tab(face - ndeb, 0) =  numer_coeff / denum_coeff;
    }

  tab.echange_espace_virtuel();
}

void Flux_radiatif_VDF::calculer_flux_radiatif(const Equation_base& eq_temp)
{
  // Retrieve the boundary temperature
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());
  const int nb_faces = le_bord.nb_faces();
  const Conds_lim& les_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites();
  OBS_PTR(Champ_front_base) Tb;

  int test_nom = 0;
  for (int num_cl_temp = 0; num_cl_temp < les_cl_temp.size(); num_cl_temp++)
    {
      const Cond_lim& la_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites(num_cl_temp);
      Nom nom_cl_temp = la_cl_temp->frontiere_dis().le_nom();
      if (nom_cl_temp == frontiere_dis().le_nom())
        {
          test_nom = 1;
          if (sub_type(Neumann_paroi_rayo_semi_transp_VDF, la_cl_temp.valeur()))
            {
              const Neumann_paroi_rayo_semi_transp_VDF& la_cl_temper = ref_cast(Neumann_paroi_rayo_semi_transp_VDF, la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Echange_contact_rayo_semi_transp_VDF, la_cl_temp.valeur()))
            {
              Echange_contact_rayo_semi_transp_VDF& la_cl_temper = ref_cast_non_const(Echange_contact_rayo_semi_transp_VDF, la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Echange_externe_impose_rayo_semi_transp, la_cl_temp.valeur()))
            {
              Echange_externe_impose_rayo_semi_transp& la_cl_temper = ref_cast_non_const(Echange_externe_impose_rayo_semi_transp, la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Frontiere_ouverte_temperature_imposee_rayo_semi_transp, la_cl_temp.valeur()))
            {
              const Frontiere_ouverte_temperature_imposee_rayo_semi_transp& la_cl_temper = ref_cast(Frontiere_ouverte_temperature_imposee_rayo_semi_transp, la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Frontiere_ouverte_rayo_semi_transp, la_cl_temp.valeur()))
            {
              const Frontiere_ouverte_rayo_semi_transp& la_cl_temper = ref_cast(Frontiere_ouverte_rayo_semi_transp, la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Echange_global_impose_rayo_semi_transp, la_cl_temp.valeur()))
            {
              Echange_global_impose_rayo_semi_transp& la_cl_temper = ref_cast_non_const(Echange_global_impose_rayo_semi_transp, la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else
            {
              Cerr << "Implementation needed for other boundary conditions of the temperature equation 1 " << finl;
              Process::exit();
            }
        }
    }
  if (test_nom == 0)
    {
      Cerr << "Error: there is no boundary condition on a boundary named: " << le_nom() << finl;
      Process::exit();
    }

  // Tb contains the boundary temperatures
  // Compute the radiative flux
  DoubleTab& Flux = flux_radiatif_->valeurs();
  Flux.resize(le_bord.nb_faces(), 1);
  const Eq_rayo_semi_transp& eq_rayo = ref_cast(Eq_rayo_semi_transp, domaine_Cl_dis().equation());
  const Fluide_base& fluide = eq_rayo.fluide();
  const DoubleTab& kappa = fluide.kappa().valeurs();
  const DoubleTab& indice = fluide.indice().valeurs();
  const DoubleTab& irradiance = eq_rayo.inconnue().valeurs();

  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_Cl_dis().domaine_dis());
  const IntTab& face_voisins = zvdf.face_voisins();
  const DoubleVect& face_surfaces = zvdf.face_surfaces();
  const double sigma = eq_rayo.pb_rayo_semi_transp().valeur_sigma();

  assert(fluide.kappa().nb_comp() == 1);
  assert(emissivite().nb_comp() == 1);
  assert(fluide.indice().nb_comp() == 1);
  assert(Tb->nb_comp() == 1);

  double bilan_flux = 0.;
  const int ndeb = le_bord.num_premiere_face();
  double kappa_F = -123., epsi = -123.;
  double Tbord = -123., n = -123.;

  // Loop over the faces
  for (int face = 0; face < nb_faces; face++)
    {
      int elem = face_voisins(face + ndeb, 0);
      if (elem < 0)
        elem = face_voisins(face + ndeb, 1);

      if (sub_type(Champ_Uniforme, fluide.kappa()))
        kappa_F = kappa(0, 0);
      else
        kappa_F = kappa(elem, 0);

      if (sub_type(Champ_front_uniforme, emissivite()))
        epsi = emissivite().valeurs()(0, 0);
      else
        epsi = emissivite().valeurs()(face, 0);

      if (sub_type(Champ_Uniforme, fluide.indice()))
        n = indice(0, 0);
      else
        n = indice(elem, 0);

      if (sub_type(Champ_front_uniforme, Tb.valeur()))
        Tbord = Tb->valeurs()(0, 0);
      else
        Tbord = Tb->valeurs()(face, 0);

      const double G_F = irradiance(elem);
      const double eF = zvdf.dist_norm_bord(face + ndeb);

      double denum = A_ * (2 - epsi);
      denum /= 3 * kappa_F * epsi;
      denum += eF;

      const double numer = G_F - 4 * n * n * sigma * pow(Tbord, 4);
      double grad_G = numer / denum;

      if (epsi < DMINFLOAT)
        Flux(face, 0) = 0;
      else
        Flux(face, 0) = -grad_G / (3 * kappa_F);

      bilan_flux += face_surfaces(face + ndeb) * Flux(face, 0);
    }

  Debog::verifier_bord(" Flux_radiatif_VDF::calculer_flux_radiatif_Flux ", Flux, ndeb);

  if (eq_rayo.schema_temps().limpr())
    Cout << "Radiative flux on the boundary " << le_bord.le_nom() << ": " << bilan_flux << finl;
}
