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
#ifndef Partitionneur_Sous_Domaines_included
#define Partitionneur_Sous_Domaines_included

#include <Partitionneur_base.h>
#include <TRUST_Ref.h>

#include <Domaine_forward.h>

/*! @brief Domain partitioner based on sub-domains of the domain. See construire_partition().
 *
 */
class Partitionneur_Sous_Domaines : public Partitionneur_base
{
  Declare_instanciable_with_param(Partitionneur_Sous_Domaines);

public:
  void associer_domaine(const Domaine& domaine) override;
  void initialiser(const Noms& noms_sous_domaines);
  void construire_partition(IntVect& elem_part, int& nb_parts_tot) const override;

private:
  // Partitioner parameters
  OBS_PTR(Domaine) ref_domaine_;
  // Names of sub-domains to use
  Noms noms_sous_domaines_;
  // Names of domains to use
  Noms noms_domaines_;
};
#endif
