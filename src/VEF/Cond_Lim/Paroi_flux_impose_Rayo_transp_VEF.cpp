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
#include <Paroi_flux_impose_Rayo_transp_VEF.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Domaine_VF.h>

Implemente_instanciable(Paroi_flux_impose_Rayo_transp_VEF, "Paroi_flux_impose_Rayo_transp_VEF", Paroi_flux_impose_Rayo_transp);

Sortie& Paroi_flux_impose_Rayo_transp_VEF::printOn(Sortie& s) const { return s; }

Entree& Paroi_flux_impose_Rayo_transp_VEF::readOn(Entree& is) { return Paroi_flux_impose_Rayo_transp::readOn(is); }

void Paroi_flux_impose_Rayo_transp_VEF::completer()
{
  Paroi_flux_impose_Rayo_transp::completer();

  const DoubleTab& T_p = mon_dom_cl_dis->equation().inconnue().valeurs();
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF, frontiere_dis());
  int ndeb = la_frontiere_VF.num_premiere_face();
  int nb_faces_bord = la_frontiere_VF.nb_faces();

  for (int numfa = 0; numfa < nb_faces_bord; numfa++)
    teta_i_[numfa] = T_p(numfa + ndeb);
}

void Paroi_flux_impose_Rayo_transp_VEF::calculer_Teta_i()
{
  const DoubleTab& T_p = mon_dom_cl_dis->equation().inconnue().valeurs();
  double Temp;
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF, frontiere_dis());
  int ndeb = la_frontiere_VF.num_premiere_face();
  int nb_faces_bord = la_frontiere_VF.nb_faces();
  int is_relax = 1;
  if (le_modele_rayo->relaxation() == 0)
    is_relax = 0;
  for (int numfa = 0; numfa < nb_faces_bord; numfa++)
    {
      Temp = T_p(numfa + ndeb);
      double omega = 0.8;
      if (is_relax == 0)
        omega = 1.;
      teta_i_[numfa] = omega * (Temp) + (1. - omega) * teta_i_[numfa];
    }
}
