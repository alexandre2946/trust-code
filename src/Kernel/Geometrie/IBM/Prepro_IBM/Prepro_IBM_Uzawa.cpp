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

#include <Prepro_IBM_Uzawa.h>

Implemente_instanciable( Prepro_IBM_Uzawa,"Prepro_IBM_Uzawa|methode_IBM_Uzawa",Prepro_IBM_base );
// XD Prepro_IBM_Uzawa Prepro_IBM_base methode_IBM_Uzawa BRACE To perform the intersection of an IB (Lagrange mesh) in a
// XD_CONT MED-format file .med with the Euler computional mesh.

Sortie& Prepro_IBM_Uzawa::printOn(Sortie& os) const { return Prepro_IBM_base::printOn(os); }

void Prepro_IBM_Uzawa::set_param(Param& param) const
{
  Prepro_IBM_base::set_param(param);
  param.ajouter("choix_de_la_methode_uzawa", &lvl_, Param::OPTIONAL); // XD_ADD_P entier
  // XD_CONT choix de la methode d'Uzawa
}

Entree& Prepro_IBM_Uzawa::readOn(Entree& is)
{
  Prepro_IBM_base::readOn(is);
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);

  if (lvl_ == 1)
    Cout << "Uzawa : methode LVL1" << finl;
  else if (lvl_ == 2)
    Cout << "Uzawa : methode LVL2" << finl;
  else
    {
      Cerr << "Prepro_IBM_Uzawa : choix_de_la_methode_uzawa : invalide argument" << finl;
      Process::exit();
    }

  return is;
}

void Prepro_IBM_Uzawa::associer_pb(const Probleme_base& pb)
{
  Prepro_IBM_base::associer_pb(pb);
  compute_solid_fluid(0);

  // Ecriture eventuelle
  if( save_prepro_ == 1) Save_Med_File();
}

void Prepro_IBM_Uzawa::compute_solid_fluid(int maj_from_ext)
{
}
