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

#include <LecFicDiffuse.h>
#include <communications.h>

Implemente_instanciable_sans_constructeur(LecFicDiffuse,"LecFicDiffuse",Lec_Diffuse_base);

Entree& LecFicDiffuse::readOn(Entree& s)
{
  throw;
}

Sortie& LecFicDiffuse::printOn(Sortie& s) const
{
  throw;
}

// LecFicDiffuse::LecFicDiffuse():EFichier()
LecFicDiffuse::LecFicDiffuse(): Lec_Diffuse_base()
{
  set_error_action(ERROR_CONTINUE);
}

/*! @brief ouverture du fichier name.
 *
 * Cette methode doit etre appelee sur tous les processeurs. En cas
 *   d'echec : exit()
 *
 */
LecFicDiffuse::LecFicDiffuse(const char* name,
                             IOS_OPEN_MODE mode_): Lec_Diffuse_base()
{
  set_error_action(ERROR_CONTINUE);
  int ok = ouvrir(name, mode_);
  if (!ok && Process::je_suis_maitre())
    {
      Cerr << "File " << name << " does not exist (LecFicDiffuse)" << finl;
      Process::exit();
    }
}

/*! @brief Ouverture du fichier.
 *
 * Cette methode doit etre appelee par tous les processeurs du groupe.
 *  Valeur de retour: 1 si ok, 0 sinon
 *
 */
int LecFicDiffuse::ouvrir(const char* name,
                          IOS_OPEN_MODE mode_)
{
  int ok = 0;

{
  int x = 2* (Process::me() == 0);
  envoyer_broadcast(x, 0);
  assert(x == 2);
 }

  Process::barrier();
  Process::barrier();
  Process::barrier();

  //if(Process::je_suis_maitre()) {
    Entree_Fichier_base::ouvrir(name, mode_);

{
  int x = 3 * (Process::me() == 0);
  envoyer_broadcast(x, 0);
  assert(x == 3);
 }


	ok = is_open();
  //}

{
  int x = 2 * (Process::me() == 0);
  envoyer_broadcast(x, 0);
  assert(x == 2);
 }

  Process::barrier();
  Process::barrier();
  Process::barrier();

{
  int x = 1 * (Process::me() == 0);
  envoyer_broadcast(x, 0);
  assert(x == 1);
 }

  envoyer_broadcast(ok, 0);

{
  int x = 3 * (Process::me() == 0);
  envoyer_broadcast(x, 0);
  assert(x == 3);
 }

  return ok;
}

/*! @brief
 *
 */
std::istream& LecFicDiffuse::get_istream()
{
  if(!Process::je_suis_maitre())
    {
      std::cerr << "Error get_istream (LecFicDiffuse)" << std::endl;
      Process::exit();
    }
  return Entree::get_istream();
}

/*! @brief
 *
 */
const std::istream& LecFicDiffuse::get_istream() const
{
  if(!Process::je_suis_maitre())
    {
      std::cerr << "Error get_istream (LecFicDiffuse)" << std::endl;
      Process::exit();
    }
  return get_istream();
}
