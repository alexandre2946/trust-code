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

#include <Source_rayo_semi_transp_base.h>
#include <Pb_Couple_rayo_semi_transp.h>
#include <Pb_rayo_semi_transp.h>
#include <Fluide_base.h>

Implemente_instanciable(Pb_Couple_rayo_semi_transp, "Pb_Couple_rayo_semi_transp", Probleme_Couple);

Entree& Pb_Couple_rayo_semi_transp::readOn(Entree& is) { return is; }

Sortie& Pb_Couple_rayo_semi_transp::printOn(Sortie& os) const { return Probleme_Couple::printOn(os); }

void Pb_Couple_rayo_semi_transp::initialize()
{
  // 1er truc a faire
  int nb_pb_ray = 0;
  for (int l = 0; l < nb_problemes(); l++)
    {
      Probleme_base& le_pb = ref_cast(Probleme_base, probleme(l));

      if (!sub_type(Pb_rayo_semi_transp, le_pb))
        if (sub_type(Fluide_base, le_pb.milieu()))
          if (ref_cast(Fluide_base, le_pb.milieu()).is_rayo_semi_transp())
            nb_pb_ray++;
    }

  if (nb_pb_ray > 1)
    {
      Cerr << "Pb_Couple_rayo_semi_transp::initialize - We can only treat 1 semi-transparent problem at present. You defined " << nb_pb_ray << " !!!" << finl;
      Process::exit();
    }
  else if (nb_pb_ray == 0)
    Process::exit("Pb_Couple_rayo_semi_transp::initialize - It seems you forgot to define the radiation properties in your medium !!!\n");

  // on associe le pb fluide au pb_rayo
  for (int l = 0; l < nb_problemes(); l++)
    {
      Probleme_base& le_pb = ref_cast(Probleme_base, probleme(l));

      if (!sub_type(Pb_rayo_semi_transp, le_pb))
        if (sub_type(Fluide_base, le_pb.milieu()))
          if (ref_cast(Fluide_base, le_pb.milieu()).is_rayo_semi_transp())
            {
              pb_rayo_semi_transp_->associer_probleme_fluide(le_pb);
              break;
            }
    }

  pb_rayo_semi_transp_->discretise_longueur_rayo();

  Probleme_Couple::initialize();
  Probleme_base& le_pb = pb_rayo_semi_transp_->probleme_fluide();
  // Associer le pb rayo aux sources de rayonnement
  for (int i = 0; i < le_pb.nombre_d_equations(); i++)
    {
      Sources& les_sources = le_pb.equation(i).sources();
      for (int j = 0; j < les_sources.size(); j++)
        {
          Source& la_source = les_sources[j];
          if (sub_type(Source_rayo_semi_transp_base, la_source.valeur()))
            {
              Source_rayo_semi_transp_base& source_rayo = ref_cast(Source_rayo_semi_transp_base, la_source.valeur());
              Cerr << "Association pb rayo semi transp au terme source rayo" << finl;
              source_rayo.associer_pb_rayo_semi_transp(pb_rayo_semi_transp_.valeur());
            }
        }
    }

  pb_rayo_semi_transp_->eq_rayo().resoudre(presentTime());
  pb_rayo_semi_transp_->calculer_flux_radiatif();

  for (int i = 0; i < nb_problemes(); i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base, probleme(i));
      for (int j = 0; j < pb.nombre_d_equations(); j++)
        pb.equation(j).domaine_Cl_dis().calculer_coeffs_echange(presentTime());
    }
}

int Pb_Couple_rayo_semi_transp::associer_(Objet_U& ob)
{
  if (sub_type(Pb_rayo_semi_transp, ob))
    {
      Cerr << "association du pb rayo semi transp au pb couple" << finl;
      if (pb_rayo_semi_transp_.non_nul())
        Process::exit("Attention : on ne peut associer qu'un pb de rayonnement a un Pb_Couple_rayo_semi_transp !!! \n");

      pb_rayo_semi_transp_ = ref_cast(Pb_rayo_semi_transp, ob);
      ajouter(pb_rayo_semi_transp_.valeur());

      return 1;
    }
  else
    return Probleme_Couple::associer_(ob);
}
