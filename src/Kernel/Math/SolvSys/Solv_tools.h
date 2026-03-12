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

#ifndef Solv_tools_included
#define Solv_tools_included

#include <ArrOfBit.h>
#include <TRUSTTab.h>

class Solv_tools
{
public:
  void construit_renum(const DoubleVect&);
  const ArrOfTID& get_ix() const { return ix; }

protected:
  TIDTab renum_;                // Tableau de renumerotation globale lignes matrice TRUST -> matrice CSR
  IntTab index_;                // Tableau de renumerotation locale
  ArrOfBit items_to_keep_;      // Faut t'il conserver dans la matrice CSR la ligne item de la matrice TRUST ?
  ArrOfTID ix;                  // Tableau de travail pour remplissage Vec plus rapide
  int nb_items_to_keep_ = -1;        // Nombre local d'items a conserver
  int nb_rows_ = -1;                 // Nombre de lignes locales de la matrice TRUST
  trustIdType nb_rows_tot_ = -1;             // Nombre de lignes globales de la matrice TRUST
  trustIdType decalage_local_global_ = -1;   // Decalage numerotation local/global pour matrice CSR et vecteur
  int secmem_sz_ = -1;               // (Local) second member size
};

#endif /* Solv_tools_included */
