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
// XD partitionneur_union partitionneur_deriv union BRACE Let several local domains be generated from a bigger one using
// XD_CONT the keyword create_domain_from_sub_domain, and let their partitions be generated in the usual way. Provided
// XD_CONT the list of partition files for each small domain, the keyword 'union' will partition the global domain in a
// XD_CONT conform fashion with the smaller domains.

// XD attr sous_domaines bloc_lecture sous_domaines REQ List of the partition files with the following syntaxe:
// XD_CONT {sous_domaine1 decoupage1 ... sous_domaineim decoupageim } where sous_domaine1 ... sous_zomeim are small
// XD_CONT domains names and decoupage1 ... decoupageim are partition files.


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

  param.ajouter("sous_domaines", &fic_ssz, Param::REQUIRED);
}
void Partitionneur_Union::validate_params() const
{
  Cerr << que_suis_je() << "::validate_params" << finl;

  for (const auto& p: fic_ssz)
    {
      // p.first is name of subdomain and p.second is file name
      // waiting for c++17 for syntax: for (const auto& [subdomain, filename]: fic_ssz)
      if (Interprete::objet_existant(p.first) == 0)
        {
          Process::exit(p.first + "is not an existing TRUST object");
        }
      if (not sub_type(Sous_Domaine, Interprete::objet(p.first)))
        {
          Process::exit(p.first + " is not a Sous_Domaine object");
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

  for (const auto& p: fic_ssz)
    {
      Nom subdomain = Nom(p.first);
      const std::string& filename = p.second;
      //on recupere le sous-domaine par son nom et le decoupage en ouvrant le fichier...
      const Objet_U& ref = Interprete::objet(subdomain); // need this or compiler complains of 'possibly' dangling ref
      const Sous_Domaine& ssz = ref_cast(Sous_Domaine, ref);
      EFichier file;
      file.ouvrir(filename.c_str());
      IntVect dec_ssz;
      file >> dec_ssz;
      file.close();
      //... et on remplit un morceau de elem_part avec
      if (dec_ssz.size_array() != ssz.nb_elem_tot())
        {
          Cerr << "Partitionneur_Union : incoherent element number for sub-domaine " <<subdomain << finl;
          Process::exit();
        }
      for (int i = 0; i < ssz.nb_elem_tot(); i++)
        elem_part[ssz(i)] = dec_ssz(i);
    }
}
