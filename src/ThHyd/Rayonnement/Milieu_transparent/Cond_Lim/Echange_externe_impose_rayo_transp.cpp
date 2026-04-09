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

#include <Echange_externe_impose_rayo_transp.h>
#include <Pb_Fluide_base.h>
#include <Front_VF.h>

Implemente_instanciable(Echange_externe_impose_rayo_transp, "Paroi_Echange_externe_impose_rayo_transp", Echange_externe_impose);
// XD paroi_echange_externe_impose_rayo_transp paroi_echange_externe_impose paroi_echange_externe_impose_rayo_transp -1 Radiation External type exchange condition with a heat exchange coefficient and an imposed external temperature.


Sortie& Echange_externe_impose_rayo_transp::printOn(Sortie& is) const { return is; }

Entree& Echange_externe_impose_rayo_transp::readOn(Entree& is)
{
  if (app_domains.size() == 0)
    app_domains = { Motcle("Thermique"), Motcle("indetermine") };

  Motcle motlu;
  Motcles les_motcles(2);
  {
    les_motcles[0] = "h_imp";
    les_motcles[1] = "T_ext";
  }

  int ind = 0;
  while (ind < 2)
    {
      is >> motlu;
      int rang = les_motcles.search(motlu);

      switch(rang)
        {
        case 0:
          {
            is >> h_imp_;
            break;
          }
        case 1:
          {
            is >> le_champ_front;
            break;
          }
        default:
          {
            Cerr << "Erreur a la lecture de la condition aux limites de type " << finl;
            Cerr << "Echange_externe_impose_rayo_transp " << finl;
            Cerr << "On attendait " << les_motcles << "a la place de " << motlu << finl;
            Process::exit();
          }
        }
      ind++;
    }

  if (local_min_vect(h_imp_->valeurs()) < 1.e9)
    {
      Cerr << "Erreur sur l'utilisation de la condition a la limite" << finl;
      Cerr << "Echange_externe_impose_rayo_transp. Celle ci ne peut" << finl;
      Cerr << "etre utilisee pour un probleme de rayonnement que pour " << finl;
      Cerr << "imposer une temperature sur une paroi" << finl;
      Process::exit();
    }

  return is;
}

int Echange_externe_impose_rayo_transp::initialiser(double temps)
{
  assert(le_modele_rayo_.est_nul());

  // on recupere le modele rayo ... !
  const Probleme_base& this_pb = domaine_Cl_dis().equation().probleme();
  if (sub_type(Pb_Fluide_base, this_pb)) // sinon Pb_conduction par exemple ;)
    {
      if (this_pb.milieu().is_rayo_transp())
        {
          le_modele_rayo_ = ref_cast(Pb_Fluide_base, this_pb).get_mod_rayo_transp();

          if (le_modele_rayo_->nom_pb_rayonnant() != this_pb.le_nom())
            error_pb_name(que_suis_je(), this_pb.le_nom(), le_modele_rayo_->nom_pb_rayonnant());
        }
      else
        error_non_rad_bc(que_suis_je(), this_pb.le_nom(), frontiere_dis().frontiere().le_nom(), "paroi_echange_externe_impose");
    }

  return Echange_externe_impose::initialiser(temps);
}

void Echange_externe_impose_rayo_transp::completer()
{
  Echange_externe_impose::completer();
  preparer_surface(frontiere_dis(), domaine_Cl_dis());
}

void Echange_externe_impose_rayo_transp::calculer_Teta_i()
{
  const Front_VF& front_vf = ref_cast(Front_VF, frontiere_dis());
  for (int numfa = 0; numfa < front_vf.nb_faces(); numfa++)
    teta_i_[numfa] = T_ext(numfa);
}

void Echange_externe_impose_rayo_transp::mettre_a_jour(double temps)
{
  Echange_externe_impose::mettre_a_jour(temps);
  calculer_Teta_i();
}
