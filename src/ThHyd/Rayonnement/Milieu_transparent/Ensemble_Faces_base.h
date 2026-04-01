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

#ifndef Ensemble_Faces_base_included
#define Ensemble_Faces_base_included

#include <Domaine_forward.h>
#include <Cond_lim_rayo_milieu_transp.h>
#include <TRUST_Ref.h>
#include <Motcle.h>

class Cond_lim_base;

class Ensemble_Faces_base: public Objet_U
{
  Declare_instanciable(Ensemble_Faces_base);
public:

  void associer_les_cl(Cond_lim_base&);
  void lire(const Nom&, const Nom&, const Domaine&);
  int contient(int) const;
  int is_ok() const;

  inline const Cond_lim_rayo_milieu_transp& cond_lim_rayo() const
  {
    assert(la_cond_lim_rayo_ != 0);
    return *la_cond_lim_rayo_;
  }

  inline Cond_lim_rayo_milieu_transp& cond_lim_rayo()
  {
    assert(la_cond_lim_rayo_ != 0);
    return *la_cond_lim_rayo_;
  }

  inline double surface(int numfa) const { return la_cond_lim_rayo_->surface(numfa); }
  inline double teta_i(int numfa) { return la_cond_lim_rayo_->teta_i(numfa); }
  inline const Cond_lim_base& la_cl_base() const { return les_cl_base_; }
  inline Cond_lim_base& la_cl_base() { return les_cl_base_; }
  inline const IntVect& Table_faces() const { return num_face_Ensemble_; }
  inline int nb_faces_bord() const { return nb_faces_bord_; }

protected:
  int nb_faces_bord_ = 0;
  OBS_PTR(Cond_lim_base) les_cl_base_;
  Cond_lim_rayo_milieu_transp *la_cond_lim_rayo_ = nullptr;
  IntVect num_face_Ensemble_; //contient_;
  DoubleTab positions_;
};

#endif /* Ensemble_Faces_base_included */
