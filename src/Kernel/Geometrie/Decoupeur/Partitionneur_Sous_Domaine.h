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

#ifndef Partitionneur_Sous_Domaine_included
#define Partitionneur_Sous_Domaine_included

#include <Partitionneur_base.h>
#include <TRUST_Ref.h>

#include <Domaine_forward.h>


/*! @brief Partitioner allowing the splitting of sub-domains (potentially overlapping) created by Create_domain_from_sub_domain in a "conforming" manner: the sub-domain is split in a way
 *
 *   that is "conforming" with the full domain.
 *
 *   Usage:
 *       - create a global domain, and two sub-domains (which partially overlap) for Domaine1 and Domaine2
 *       - split the global domain and write the splitting file
 *       - create Domaine1 and Domaine2 via Create_domain_from_sub_domain
 *       - split D1 and D2 with Partitionneur sous_domaine, using the global domain splitting as the source.
 *
 *   Syntax:
 *     Decouper dom_N
 *     {
 *         partitionneur sous_domaine
 *         {
 *              fichier     decoup/domaine_glob.txt
 *              fichier_ssz sous_domaine_dom_N.file
 *              OR
 *              nom_ssz     sous_domaine
 *         }
 *         Nom_Domaines decoup/dom_N
 *     }
 *
 *
 * @sa Partitionneur_Union Create_domain_from_sub_domain
 */
class Partitionneur_Sous_Domaine : public Partitionneur_base
{
  Declare_instanciable_with_param(Partitionneur_Sous_Domaine);

protected:
  void validate_params() const override;

public:
  void associer_domaine(const Domaine& dom) override { };
  void construire_partition(IntVect& elem_part, int& nb_parts_tot) const override;

protected:
  // Partitioner parameters
  Nom filename_ = "";      ///! Name of the global splitting file
  Nom filename_ssz_ = ""; ///! Name of the sub-domain file
  Nom name_ssz_ = ""; ///! Name of the sub-domain (declared in the data file)
};
#endif
