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

#ifndef Partitionneur_Fichier_Decoupage_included
#define Partitionneur_Fichier_Decoupage_included

#include <Partitionneur_base.h>
#include <TRUST_Ref.h>

#include <Domaine_forward.h>


/*! @brief Partition of a domain from a disk file containing, for each element, the processor number to which that element is assigned.
 *
 *   The file must contain an array in ArrOfInt format
 *    (number n of values followed by the values).
 *   The size of the array must equal the number of elements
 *   in the domain to be split, and each value indicates
 *   on which processor that element is to be placed.
 *
 */

class Partitionneur_Fichier_Decoupage : public Partitionneur_base
{
  Declare_instanciable_with_param(Partitionneur_Fichier_Decoupage);

public:
  void associer_domaine(const Domaine& domaine) override;
  void initialiser(const char *filename);
  void construire_partition(IntVect& elem_part, int& nb_parts_tot) const override;

protected:
  // Partitioner parameters
  OBS_PTR(Domaine) ref_domaine_;
  Nom filename_ = "";
  bool corriger_partition_ = false;
};
#endif
