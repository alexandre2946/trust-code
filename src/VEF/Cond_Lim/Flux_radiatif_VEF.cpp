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
#include <Temperature_imposee_paroi_rayo_semi_transp.h>
#include <Neumann_paroi_rayo_semi_transp_VEF.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Champ_front_uniforme.h>
#include <Eq_rayo_semi_transp.h>
#include <Pb_rayo_semi_transp.h>
#include <Flux_radiatif_VEF.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Fluide_base.h>
#include <Domaine_VEF.h>
#include <TRUST_Ref.h>
#include <Debog.h>

Implemente_instanciable(Flux_radiatif_VEF, "Flux_radiatif_VEF", Flux_radiatif_base);

Sortie& Flux_radiatif_VEF::printOn(Sortie& s) const { return s << que_suis_je() << finl; }
Entree& Flux_radiatif_VEF::readOn(Entree& s) { return Flux_radiatif_base::readOn(s); }

void Flux_radiatif_VEF::evaluer_cl_rayonnement(Champ_front_base& Tb, const Champ_Don_base& coeff_abs, const Champ_Don_base& longueur_rayo,
                                               const Champ_Don_base& indice, const Domaine_VF& zvf,
                                               const double sigma, double temps)
{
  const DoubleTab& n = indice.valeurs();
  const DoubleTab& epsilon = emissivite().valeurs();
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());

  // Size the DoubleTab associated with le_champ_front
  assert(le_champ_front->nb_comp() == 1);
  DoubleTab& tab = le_champ_front->valeurs_au_temps(temps);

  assert(emissivite().nb_comp() == 1);
  assert(indice.nb_comp() == 1);
  assert(Tb.nb_comp() == 1);

  const int ndeb = le_bord.num_premiere_face();
  const int nfin = ndeb + le_bord.nb_faces();
  double epsi = -123., nn = -123., T = -123.;

  // Loop over the faces of le_bord
  for (int face = ndeb; face < nfin; face++)
    {
      if (sub_type(Champ_front_uniforme, emissivite()))
        epsi = epsilon(0, 0);
      else
        epsi = epsilon(face - ndeb, 0);

      if (sub_type(Champ_Uniforme, indice))
        nn = n(0, 0);
      else
        nn = n(face, 0);

      if (sub_type(Champ_front_uniforme, Tb))
        T = Tb.valeurs()(0, 0);
      else
        T = Tb.valeurs_au_temps(temps)(face - ndeb, 0);

      const double numer_coeff = 4 * nn * nn * sigma * pow(T, 4) * epsi;
      const double denum_coeff = A_ * (2 - epsi);

      tab(face - ndeb, 0) = numer_coeff / denum_coeff;
    }

  tab.echange_espace_virtuel();
}

void Flux_radiatif_VEF::calculer_flux_radiatif(const Equation_base& eq_temp)
{
  // We need to retrieve the wall temperature
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
          if (sub_type(Neumann_paroi_rayo_semi_transp_VEF, la_cl_temp.valeur()))
            {
              Neumann_paroi_rayo_semi_transp_VEF& la_cl_temper = ref_cast_non_const(Neumann_paroi_rayo_semi_transp_VEF, la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Temperature_imposee_paroi_rayo_semi_transp, la_cl_temp.valeur()))
            {
              Temperature_imposee_paroi_rayo_semi_transp& la_cl_temper = ref_cast_non_const(Temperature_imposee_paroi_rayo_semi_transp, la_cl_temp.valeur());
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
          else
            {
              Cerr << "To be implemented for other boundary conditions of the temperature equation 1 " << finl;
              Process::exit();
            }
        }
    }
  if (test_nom == 0)
    {
      Cerr << "Error: there is no boundary condition on a boundary with the name: " << frontiere_dis().le_nom() << finl;
      Process::exit();
    }

  // Tb contains the wall temperatures
  // Computation of the radiative flux
  DoubleTab& Flux = flux_radiatif().valeurs();
  Flux.resize(le_bord.nb_faces(), 1);
  const Eq_rayo_semi_transp& eq_rayo = ref_cast(Eq_rayo_semi_transp, domaine_Cl_dis().equation());
  const Fluide_base& fluide = eq_rayo.fluide();
  const DoubleTab& indice = fluide.indice().valeurs();
  const DoubleTab& irradiance = eq_rayo.inconnue().valeurs();

  const Domaine_VEF& zvef = ref_cast(Domaine_VEF, domaine_Cl_dis().domaine_dis());
  const DoubleTab& face_normales = zvef.face_normales();
  const double sigma = eq_rayo.pb_rayo_semi_transp().valeur_sigma();

  assert(emissivite().nb_comp() == 1);
  assert(fluide.indice().nb_comp() == 1);
  assert(Tb->nb_comp() == 1);

  double bilan_flux = 0.;
  const int ndeb = le_bord.num_premiere_face();
  double epsi = -123., Tbord = -123., n = -123.;

  // Loop over the faces
  for (int face = 0; face < nb_faces; face++)
    {
      if (sub_type(Champ_front_uniforme, emissivite()))
        epsi = emissivite().valeurs()(0, 0);
      else
        epsi = emissivite().valeurs()(face, 0);

      if (sub_type(Champ_Uniforme, fluide.indice()))
        n = indice(0, 0);
      else
        n = indice(face + ndeb, 0);

      if (sub_type(Champ_front_uniforme, Tb.valeur()))
        Tbord = Tb->valeurs()(0, 0);
      else
        Tbord = Tb->valeurs()(face, 0);

      const double irra = irradiance(face + ndeb);

      const double denum = A_ * (2 - epsi);
      const double numer = epsi * (irra - 4 * n * n * sigma * pow(Tbord, 4));
      Flux(face, 0) = -numer / denum;

      // Compute the balance
      double surface = 0.;
      for (int i = 0; i < dimension; i++)
        surface += (face_normales(face + ndeb, i) * face_normales(face + ndeb, i));

      surface = sqrt(surface);
      bilan_flux += surface * Flux(face, 0);
    }

  Debog::verifier_bord(" Flux_radiatif_VEF::calculer_flux_radiatif_Flux ", Flux, ndeb);

  if (eq_rayo.schema_temps().limpr())
    Cout << "Radiative flux on boundary " << le_bord.le_nom() << " : " << bilan_flux << finl;
}
