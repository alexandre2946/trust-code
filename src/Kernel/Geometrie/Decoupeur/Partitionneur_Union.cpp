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
#include <Partitionneur_Union.h>
#include <EFichier.h>
#include <Domaine.h>
#include <Param.h>
#include <EChaine.h>
#include <Sous_Domaine.h>
#include <Interprete.h>

Implemente_instanciable(Partitionneur_Union,"Partitionneur_Union",Partitionneur_base);
// XD partitionneur_union partitionneur_deriv union 1 Let several local domains be generated from a bigger one using the keyword create_domain_from_sub_domain, and let their partitions be generated in the usual way. Provided the list of partition files for each small domain, the keyword 'union' will partition the global domain in a conform fashion with the smaller domains.


Sortie& Partitionneur_Union::printOn(Sortie& os) const
{
  Cerr << "Partitionneur_Union::printOn invalid\n" << finl;
  Process::exit();
  return os;
}


/*! @brief Lecture des parametres du partitionneur sur disque.
 *
 * Format attendu:
 *    {
 *      sous_domaines N ssdom1 ... ssdomN
 *      fichiers_decoupage N file1 ... fileN
 *    }
 *
 */
void Partitionneur_Union::set_param(Param& param) const
{

  param.ajouter("sous_domaines", &sous_domaines_, Param::REQUIRED); // XD_ADD_P listchaine list of sous_domaines names. They must be valid Sous_Domaine that have been declared earlier in the dataset
  param.ajouter("fichiers_decoupage", &fichiers_decoupage_, Param::REQUIRED); // XD_ADD_P listchaine list of files which contain the partittion for the corresponding subdomain
}
void Partitionneur_Union::validate_params() const
{
  Cerr << que_suis_je() << "::validate_params" << finl;
  if (sous_domaines_.size() != fichiers_decoupage_.size())
    {
      Process::exit("Expected same number of elements in sous_domaines and fichiers_decoupage");
    }

  for (int i = 0; i< sous_domaines_.size(); i++)
    {
      if (Interprete::objet_existant(sous_domaines_[i]) == 0)
        {
          Process::exit(sous_domaines_[i] + "is not an existing TRUST object");
        }
      if (not sub_type(Sous_Domaine, Interprete::objet(sous_domaines_[i])))
        {
          Process::exit(sous_domaines_[i] + " is not a Sous_Domaine object");
        }
    }

}


void Partitionneur_Union::associer_domaine(const Domaine& domaine)
{
  ref_domaine_ = domaine;
}

/*! @brief Lit le contenu du fichier "filename_" et stocke le resultat dans elem_part
 *
 */
void Partitionneur_Union::construire_partition(IntVect& elem_part, int& nb_parts_tot) const
{
  elem_part.resize(ref_domaine_->nb_elem());
  elem_part = -1;

  for (int i_dom = 0; i_dom< sous_domaines_.size(); i_dom++)
    {
      //on recupere le sous-domaine par son nom et le decoupage en ouvrant le fichier...
      const Sous_Domaine& ssz = ref_cast(Sous_Domaine, Interprete::objet(sous_domaines_[i_dom]));
      EFichier file;
      file.ouvrir(fichiers_decoupage_[i_dom].getString().c_str());
      IntVect dec_ssz;
      file >> dec_ssz;
      file.close();
      //... et on remplit un morceau de elem_part avec
      if (dec_ssz.size_array() != ssz.nb_elem_tot())
        {
          Cerr << "Partitionneur_Union : incoherent element number for sub-domaine " << sous_domaines_[i_dom] << finl;
          Process::exit();
        }
      for (int i = 0; i < ssz.nb_elem_tot(); i++)
        elem_part[ssz(i)] = dec_ssz(i);
    }
}
