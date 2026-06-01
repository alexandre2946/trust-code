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
  TIDTab renum_;                // Global row renumbering table: TRUST matrix -> CSR matrix
  IntTab index_;                // Local renumbering table
  ArrOfBit items_to_keep_;      // Whether to keep row item in the CSR matrix from the TRUST matrix
  ArrOfTID ix;                  // Work array for faster Vec filling
  int nb_items_to_keep_ = -1;        // Local number of items to keep
  int nb_rows_ = -1;                 // Number of local rows in the TRUST matrix
  trustIdType nb_rows_tot_ = -1;             // Number of global rows in the TRUST matrix
  trustIdType decalage_local_global_ = -1;   // Local/global index offset for the CSR matrix and vector
  int secmem_sz_ = -1;               // (Local) second member size
};

#endif /* Solv_tools_included */
