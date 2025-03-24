/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Paroi_flux_impose_Rayo_transp_VDF.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Domaine_VDF.h>

Implemente_instanciable(Paroi_flux_impose_Rayo_transp_VDF, "Paroi_flux_impose_Rayo_transp_VDF", Paroi_flux_impose_Rayo_transp);

Sortie& Paroi_flux_impose_Rayo_transp_VDF::printOn(Sortie& s) const { return s; }

Entree& Paroi_flux_impose_Rayo_transp_VDF::readOn(Entree& is) { return Paroi_flux_impose_Rayo_transp::readOn(is); }

void Paroi_flux_impose_Rayo_transp_VDF::completer()
{
  Paroi_flux_impose_Rayo_transp::completer();

  const Domaine_VF& domaine_VDF = ref_cast(Domaine_VDF, domaine_Cl_dis().domaine_dis());
  const DoubleTab& T_f = mon_dom_cl_dis->equation().inconnue().valeurs();
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF, frontiere_dis());
  int ndeb = la_frontiere_VF.num_premiere_face();
  int nb_faces_bord = la_frontiere_VF.nb_faces();

  const IntTab& face_voisins = domaine_VDF.face_voisins();

  for (int numfa = 0; numfa < nb_faces_bord; numfa++)
    {
      int elem = face_voisins(numfa + ndeb, 0);
      if (elem < 0)
        elem = face_voisins(numfa + ndeb, 1);
      teta_i_(numfa) = T_f(elem);
    }
}

void Paroi_flux_impose_Rayo_transp_VDF::calculer_Teta_i()
{
  const Domaine_VDF& le_dom_vdf = ref_cast(Domaine_VDF, domaine_Cl_dis().domaine_dis());
  const Milieu_base& le_milieu = mon_dom_cl_dis->equation().milieu();
  const DoubleTab& T_f = mon_dom_cl_dis->equation().inconnue().valeurs();
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF, frontiere_dis());
  int ndeb = la_frontiere_VF.num_premiere_face();
  int nb_faces_bord = la_frontiere_VF.nb_faces();
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_Cl_dis().domaine_dis());
  const IntTab& face_voisins = zvdf.face_voisins();
  int is_rho_unif = 0;
  int is_conduc_unif = 0;
  int is_Cp_unif = 0;

  double d_rho = 0;
  double d_Lambda = 0;
  double d_Cp = 0;

  const DoubleTab& rho = le_milieu.masse_volumique().valeurs();
  const DoubleTab& Lambda = le_milieu.conductivite().valeurs();
  const DoubleTab& Cp = le_milieu.capacite_calorifique().valeurs();

  if (sub_type(Champ_Uniforme, le_milieu.masse_volumique()))
    {
      is_rho_unif = 1;

      d_rho = rho(0, 0);
    }
  if (sub_type(Champ_Uniforme, le_milieu.conductivite()))
    {
      is_conduc_unif = 1;
      d_Lambda = Lambda(0, 0);
    }

  if (sub_type(Champ_Uniforme, le_milieu.capacite_calorifique()))
    {
      is_Cp_unif = 1;
      d_Cp = Cp(0, 0);
    }

  Schema_Temps_base& sch = mon_dom_cl_dis->equation().probleme().schema_temps();
  double dt = sch.pas_de_temps();

  int is_relax = 1;
  if (le_modele_rayo->relaxation() == 0)
    is_relax = 0;

  for (int numfa = 0; numfa < nb_faces_bord; numfa++)
    {
      // QUI QU'A BU ?????
      // T_f (numfa)!!!! balaise Tf(elem) plus judicieux!!!
      int elem = face_voisins(numfa + ndeb, 0);
      if (elem < 0)
        elem = face_voisins(numfa + ndeb, 1);

      if (!is_conduc_unif)
        d_Lambda = Lambda(elem);
      double omega;
      double e = le_dom_vdf.dist_norm_bord(numfa + ndeb);
      if (is_relax)
        {

          if (!is_rho_unif)
            d_rho = rho(elem);
          if (!is_Cp_unif)
            d_Cp = Cp(elem);

          omega = d_Lambda * dt / (d_Lambda * dt + e * e * d_rho * d_Cp);
        }
      else
        omega = 1.;

      double flux_radia = le_modele_rayo->flux_radiatif(numfa + ndeb);

      if (le_champ_front->valeurs().size() == 1)
        teta_i_(numfa) = omega * ((le_champ_front->valeurs()(0, 0) - flux_radia) / (d_Lambda / e) + T_f(elem)) + (1 - omega) * teta_i_(numfa);
      else if (le_champ_front->valeurs().dimension(1) == 1)
        teta_i_(numfa) = omega * ((le_champ_front->valeurs()(numfa, 0) - flux_radia) / (d_Lambda / e) + T_f(elem)) + (1 - omega) * teta_i_(numfa);
      else
        {
          Cerr << "Paroi_flux_impose_Rayo_transp::calculer_Teta_i() erreur" << finl;
          Process::exit();
        }
    }

  // Impression:
  if (zvdf.domaine().bords_a_imprimer().contient(la_frontiere_VF.le_nom()) && sch.limpr())
    {
      Cout << "Impression des temperatures de paroi sur la frontiere " << la_frontiere_VF.le_nom() << " :" << finl;
      Cout << "---------------------------------------------------------------------" << finl;
      for (int numfa = 0; numfa < nb_faces_bord; numfa++)
        Cout << "T(" << numfa << ") : " << teta_i_(numfa) << " K." << finl;
    }
}
