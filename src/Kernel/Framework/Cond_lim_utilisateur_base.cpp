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

#include <Cond_lim_utilisateur_base.h>
#include <Entree_complete.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <SFichier.h>

Implemente_base(Cond_lim_utilisateur_base, "Cond_lim_utilisateur_base", Cond_lim_base);

Sortie& Cond_lim_utilisateur_base::printOn(Sortie& s) const { return s << que_suis_je(); }

Entree& Cond_lim_utilisateur_base::readOn(Entree& s) { return s; }

void Cond_lim_utilisateur_base::ecrire(const Nom& ajout)
{
  if (je_suis_maitre())
    {
      SFichier conv("convert_jdd", ios::app);
      conv << (*this) << " # " << ajout << finl;
    }
}

void Cond_lim_utilisateur_base::lire(Entree& s, Equation_base& mon_eq, const Nom& nom_bord)
{
  nom_bord_ = nom_bord;
  mon_equation = mon_eq;
  la_cl_ = new Cond_lim;
  Nom ajout("");
  complement(ajout);
#ifndef NDEBUG
  ecrire(ajout);
#endif
  Entree_complete s_complete(ajout, s);

  s_complete >> *(la_cl_);
  Cerr << "end reading cond lim util" << finl;
}

Cond_lim& Cond_lim_utilisateur_base::la_cl()
{
  return *(la_cl_);
}

void Cond_lim_utilisateur_base::complement(Nom&)
{
  Cerr << "Cond_lim_utilisateur_base::complement(Nom& ) does nothing" << finl;
}

/*! @brief Returns 0 if the problem is not radiating, 1 if it is semi-transparent,
 *
 *                      2 if it is transparent.
 *
 */
int Cond_lim_utilisateur_base::is_pb_rayo()
{
  const Milieu_base& milieu = ref_cast(Milieu_base, mon_equation->probleme().milieu());

  if (milieu.is_rayo_transp())
    return 2;
  else if (milieu.is_rayo_semi_transp())
    return 1;
  else
    return 0;
}
