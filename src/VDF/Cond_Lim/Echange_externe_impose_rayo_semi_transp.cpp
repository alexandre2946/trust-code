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

#include <Echange_externe_impose_rayo_semi_transp.h>
#include <Motcle.h>

Implemente_instanciable(Echange_externe_impose_rayo_semi_transp, "Paroi_Echange_externe_impose_rayo_semi_transp", Echange_externe_impose);

Sortie& Echange_externe_impose_rayo_semi_transp::printOn(Sortie& os) const
{
  return os;
}

Entree& Echange_externe_impose_rayo_semi_transp::readOn(Entree& is)
{
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
            Cerr << "Echange_externe_impose_rayo_semi_transp " << finl;
            Cerr << "On attendait " << les_motcles << "a la place de " << motlu << finl;
            Process::exit();
          }
        }
      ind++;
    }

  if (local_min_vect(h_imp_->valeurs()) < 1.e9)
    {
      Cerr << "Erreur sur l'utilisation de la condition a la limite" << finl;
      Cerr << "Echange_externe_impose_rayo_semi_transp. Celle ci ne peut" << finl;
      Cerr << "etre utilisee pour un probleme de rayonnement que pour " << finl;
      Cerr << "imposer une temperature sur une paroi" << finl;
      Process::exit();
    }

  return is;
}

void Echange_externe_impose_rayo_semi_transp::completer()
{
  Echange_externe_impose::completer();
}

const Cond_lim_base& Echange_externe_impose_rayo_semi_transp::la_cl() const
{
  return (*this);
}

Champ_front_base& Echange_externe_impose_rayo_semi_transp::temperature_bord()
{
  return T_ext();
}
