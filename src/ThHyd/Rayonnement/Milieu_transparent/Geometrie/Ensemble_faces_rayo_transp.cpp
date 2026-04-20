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

#include <Ensemble_faces_rayo_transp.h>
#include <Frontiere_dis_base.h>
#include <TRUSTList.h>
#include <EFichier.h>
#include <Domaine.h>

Implemente_instanciable(Ensemble_faces_rayo_transp, "Ensemble_faces_rayo_transp", Objet_U);

Entree& Ensemble_faces_rayo_transp::readOn(Entree& is) { return is; }

Sortie& Ensemble_faces_rayo_transp::printOn(Sortie& os) const { return os; }

void Ensemble_faces_rayo_transp::lire(const Nom& nom_bord_lu, const Nom& nom_bord, const Domaine& dom)
{
  {
    Nom fic2(dom.le_nom() + "." + nom_bord);
    fic2 += "_xv";
    EFichier fichier(fic2);
    Nom motlu;
    while (motlu != nom_bord_lu)
      {
        fichier >> motlu;
      }
    fichier >> positions_;
    fichier.close();
  }
}

int Ensemble_faces_rayo_transp::contient(int num_face) const
{
  // Dans notre cas une face rayonnantes est exactement une face de bord
  if (num_face_Ensemble_.size() == 0)
    return 1;
  else
    {
      int nb_faces_rayo = num_face_Ensemble_.size();

      for (int i = 0; i < nb_faces_rayo; i++)
        if (num_face == num_face_Ensemble_[i])
          return 1;

      return 0;
    }
}

int Ensemble_faces_rayo_transp::is_ok() const
{
  if (les_cl_base_)
    return 1;
  else
    return 0;
}

void Ensemble_faces_rayo_transp::associer_les_cl(Cond_lim_base& la_cl)
{
  la_cl.is_bc_rayo_milieu_transp(la_cond_lim_rayo_);
  les_cl_base_ = la_cl;
  Frontiere& le_bord = la_cl.frontiere_dis().frontiere();
  nb_faces_bord_ = le_bord.nb_faces();
  // On construit num_face_Ensemble. On a les positions_ on cherche la liste des faces de ce bord
  IntList numface;
  DoubleTab pos;

  int n1 = positions_.dimension(0);
  if (n1 != 0)
    {
      Faces& faces = le_bord.faces();
      // pos contient les centres de gravite des faces du bord
      if (nb_faces_bord_)
        Faces::Calculer_centres_gravite(pos, faces.type_face(), le_bord.domaine().coord_sommets(), faces.les_sommets());
      for (int fac = 0; fac < n1; fac++)
        {
          int marq = 0;
          for (int f2 = 0; f2 < nb_faces_bord_; f2++)
            {
              int ok = 1;
              for (int dir = 0; dir < dimension; dir++)
                {
                  if (!(est_egal(positions_(fac, dir), pos(f2, dir))))
                    ok = 0;
                }
              if (ok == 1)
                {
                  numface.add(f2);
                  marq++;
                  if (marq != 1)
                    {
                      Cerr << "Error in Ensemble_faces_rayo_transp::associer_les_cl" << finl;
                      Cerr << "Contact TRUST support." << finl;
                      Cerr << fac << " face en double " << finl;
                      Process::exit();
                    }
                }
            }
          marq = (int) mp_sum((double) marq);
          if (marq == 0)
            {
              Cerr << "Face " << fac << " du fichier " << la_cl.frontiere_dis().frontiere().le_nom() << "_xv non trouvee !!! positions ";
              for (int dir = 0; dir < dimension; dir++)
                Cerr << positions_(fac, dir) << " ";
              Cerr << finl;
              for (int f2 = 0; f2 < nb_faces_bord_; f2++)
                {
                  Cerr << " face " << f2 << " du bord ";
                  for (int dir = 0; dir < dimension; dir++)
                    Cerr << pos(f2, dir) << " ";
                  Cerr << finl;
                }
              Process::exit();
            }

        }
      int nf = numface.size();
      if (nf == 0)
        nb_faces_bord_ = 0;
      num_face_Ensemble_.resize(nf);
      for (int f3 = 0; f3 < nf; f3++)
        num_face_Ensemble_[f3] = numface[f3];
      positions_.resize(0);
    }
}
