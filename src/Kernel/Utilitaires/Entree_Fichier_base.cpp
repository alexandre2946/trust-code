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

#include <Entree_Fichier_base.h>
#include <Process.h>
#include <Nom.h>

#ifndef LATATOOLS
#include <EntreeSortie.h>
#endif

Implemente_base_sans_constructeur_ni_destructeur(Entree_Fichier_base,"Entree_Fichier_base",Objet_U);

Entree& Entree_Fichier_base::readOn(Entree& s)
{
  throw;
}

Sortie& Entree_Fichier_base::printOn(Sortie& s) const
{
  throw;
}


std::ifstream& Entree_Fichier_base::get_ifstream()
{
  return *this;
}

Entree_Fichier_base::~Entree_Fichier_base()
{
  //Entree_Fichier_base::close();
}

int Entree_Fichier_base::ouvrir(const char* name, IOS_OPEN_MODE mode_)
{
	//if (is_open()) {
		close();
		clear();
	//}

  	IOS_OPEN_MODE ios_mod = mode_;

	if (is_bin) {
		ios_mod=ios_mod|ios::binary;
	}

	open(name, ios_mod);

	int ok = good();

	if (is_bin) {
      Nom test;
      (*this) >> test;
      if (test == "INT64")
        {
		  set_64_bits(true);
#ifndef INT_is_64_
          Cerr<<"Opening " <<name<< " which is an int64 binary file..."<<finl;
#endif
        }
      else
        {
		  set_64_bits(false);
#ifdef INT_is_64_
          Cerr<<"Opening " <<name<< " which is an int32 binary file..."<<finl;
#endif
          // rewind, to go back at begining of file:
		  close();
		  clear();
		  open(name, ios_mod);
          ok = good();
        }
    }
  return ok;
}
