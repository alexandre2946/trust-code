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

#include <Modele_rayo_transp.h>
#include <Domaine_Cl_dis_base.h>
#include <Cond_lim_rayo_milieu_transp.h>
#include <Domaine_VF.h>

void Cond_lim_rayo_milieu_transp::associer_modele_rayo(const Modele_rayo_transp& mod)
{
  le_modele_rayo_ = mod;
}

void Cond_lim_rayo_milieu_transp::completer()
{
  Cerr << "Cond_lim_rayo_milieu_transp::doit etre surchargee" << finl;
  Process::exit();
}

void Cond_lim_rayo_milieu_transp::preparer_surface(const Frontiere_dis_base& fr, const Domaine_Cl_dis_base& zcl)
{
  const Front_VF& la_frontiere_VF = ref_cast(Front_VF, fr);
  const int ndeb = la_frontiere_VF.num_premiere_face();
  const int nb_faces_bord = la_frontiere_VF.nb_faces();

  surf_i_.resize(nb_faces_bord);
  teta_i_.resize(nb_faces_bord);

  // recuperation des surfaces de bords.
  const Domaine_VF& domaine = ref_cast(Domaine_VF, zcl.domaine_dis());
  for (int numfa = 0; numfa < nb_faces_bord; numfa++)
    surf_i_[numfa] = domaine.face_surfaces(numfa + ndeb);
}
