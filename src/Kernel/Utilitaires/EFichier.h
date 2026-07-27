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

#ifndef EFichier_included
#define EFichier_included


#include <arch.h>
#include <Entree_Fichier_base.h>
#include <Declare_Inst.h>

/*! @brief Fichier en lecture Cette classe est a la classe C++ ifstream ce que la classe Entree est a la
 *
 *     classe C++ istream. Elle redefinit de facon virtuelle les operateurs de lecture dans un fichier.
 *
 */

class EFichier : public Entree_Fichier_base
{

  Declare_instanciable(EFichier);

public:

  using Entree_Fichier_base::Entree_Fichier_base;
  using Input::operator>>;

  EFichier(const char* name,IOS_OPEN_MODE mode=std::ios::in);


    	#pragma GCC diagnostic push
    	#pragma GCC diagnostic ignored "-Wextra"
		EFichier(const EFichier& other): Entree_Fichier_base(dynamic_cast<const Entree_Fichier_base&>(other)) {}
		EFichier(EFichier& other): Entree_Fichier_base(dynamic_cast<Entree_Fichier_base&>(other)) {}
    	#pragma GCC diagnostic pop
};

#endif
