/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Reordonner.h>
#include <Domaine.h>

Implemente_instanciable_32_64(Reordonner_32_64,"Reordonner",Interprete_geometrique_base_32_64<_T_>);
// XD resequencing interprete reordonner INHERITS_BRACE The Reordonner interpretor is required sometimes for a VDF mesh
// XD_CONT which is not produced by the internal mesher. Example where this is used: NL2 Read_file dom fichier.geom NL2
// XD_CONT Reordonner dom NL2 Observations: This keyword is redundant when the mesh that is read is correctly sequenced
// XD_CONT in the TRUST sense. This significant mesh operation may take some time... The message returned by TRUST is
// XD_CONT not explicit when the Reordonner (Resequencing) keyword is required but not included in the data set...
// XD attr domain_name ref_domaine domain_name REQ Name of domain to resequence.


/*! @brief Writes the object type to an output stream.
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the modified output stream
 */
template <typename _SIZE_>
Sortie& Reordonner_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  return s << this->que_suis_je() << finl;
}


/*! @brief DOES NOTHING
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the input stream
 */
template <typename _SIZE_>
Entree& Reordonner_32_64<_SIZE_>::readOn(Entree& is )
{
  return is;
}


/*! @brief Main function of the Mailler interpreter. Reorders the nodes of the domain specified by
 *
 *     the directive.
 *        Reordonner_32_64 dom
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 * @throws the object to reorder is not of type Domaine
 */
template <typename _SIZE_>
Entree& Reordonner_32_64<_SIZE_>::interpreter_(Entree& is)
{
  this->associer_domaine(is);
  this->domaine().reordonner();
  return is;
}


template class Reordonner_32_64<int>;
#if INT_is_64_ == 2
template class Reordonner_32_64<trustIdType>;
#endif


