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

#ifndef Decouper_included
#define Decouper_included

#include <Partitionneur_base.h>
#include <TRUST_Deriv.h>
#include <Interprete_geometrique_base.h>
#include <TRUST_Ref.h>
#include <vector>
#include <Domaine_forward.h>

enum class DomainesFileOutputType { BINARY_MULTIPLE, HDF5_SINGLE };

/*! @brief Interprete Decouper.
 *
 * No algorithm here, only reading of parameters from the .data file and execution of the partitioner
 *   and the cutter. See the interprete() method.
 *
 */
template <typename _SIZE_>
class Decouper_32_64 : public Interprete_geometrique_base_32_64<_SIZE_>
{
  Declare_instanciable_32_64(Decouper_32_64);
public:
  // Those two classes are tightly related and are just easier to code if they are allowed to access members directly:
  friend class Decouper_multi;
  friend class Decouper_etendu;  // For STT

  using int_t = _SIZE_;
  using BigIntVect_t = TRUSTVect<int, _SIZE_>;  // always storing 'int' (=proc num) but might have a lot of entries
  using IntTab_t = IntTab_T<_SIZE_>;
  using Domaine_t = Domaine_32_64<_SIZE_>;
  using Partitionneur_base_t = Partitionneur_base_32_64<_SIZE_>;
  using Static_Int_Lists_t = Static_Int_Lists_32_64<_SIZE_>;

  Entree& lire(Entree& is);                 //reading of parameters
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  Entree& interpreter(Entree& is) override;

  /** Writing of a given elem_part partition.
   *   som_raccord (optional) : som_raccord[s] -> processes to which vertex s is connected by a connector to another domain.
   */
  void ecrire(const Static_Int_Lists_t *som_raccord=nullptr);

  // Whether to be more verbose:
  static int print_more_infos_;

protected:
  // The actual tool used for partitioning (Metis, Tranche, etc.):
  OWN_PTR(Partitionneur_base_t) deriv_partitionneur_;
  int nb_parts_tot_ = -1;
  // Result of the partitionning process:
  BigIntVect_t elem_part_;

  // Splitting parameters (filled by lire()):
  Nom nom_domaine_;
  int epaisseur_joint_ = 1;
  Nom nom_domaines_decoup_ = "?";
  Nom nom_fichier_decoupage_ = "?";
  Nom nom_fichier_decoupage_sommets_ = "?";
  Nom nom_fichier_lata_ = "?";
  Nom nom_fichier_med_ = "?";
  DomainesFileOutputType format_ = DomainesFileOutputType::BINARY_MULTIPLE;
  int reorder_ = 0;

private:
  // Does nothing, we directly override interpreter():
  Entree& interpreter_(Entree& is) override { return is; }

  void lire_partitionneur(Entree& is);
  void ecrire_fichier_decoupage() const;
  void ecrire_fichier_decoupage_som() const;
  void postraiter_decoupage(const Nom& nom_fichier) const;
  void ecrire_sous_domaines(const int nb_parties, const Static_Int_Lists_t* som_raccord) const;

};

using Decouper = Decouper_32_64<int>;
using Decouper_64 = Decouper_32_64<trustIdType>;

#endif
