/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Convert_ICoCoTrioField.h>
#include <Champ_Input_Proto.h>
#include <Schema_Temps_base.h>
#include <ICoCoExceptions.h>
#include <Probleme_base.h>
#include <Interprete.h>
#include <TRUSTTab.h>
#include <Domaine.h>

using ICoCo::WrongArgument;

/*! @brief Reading an input field from an input stream Factorization of the readOn methods of the different Champ_input
 *
 *     Format:
 *      { nb_comp number_of_field_components
 *        nom field_name
 *        probleme name_of_associated_problem_base
 *        [sous_domaine name_of_associated_sous_domaine] }
 *
 * @param (Entree& is) input stream
 * @return the modified input stream
 * @throws no opening brace at the start of the format
 * Exception: one of the 3 keywords missing, or another unrecognized keyword
 * @throws no closing brace at the end of the data set
 */
void Champ_Input_Proto::read(Entree& is)
{
  Nom nom_champ;
  Nom nom_pb;
  Nom nom_sous_domaine;
  int nbc=-1;

  // Read nom_champ and nom_pb
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "It was expected { to start to read the Champ_Input" << finl;
      Process::exit();
    }
  is >> motlu;
  while(1)
    {
      if (motlu == "nb_comp")
        {
          is >> nbc;
          set_nb_comp(nbc);
        }
      else if (motlu == "nom")
        {
          is >> nom_champ;
        }
      else if (motlu == "initial_value")
        {
          is >> initial_value_;
        }
      else if (motlu=="probleme")
        {
          is >> nom_pb;
        }
      else if (motlu=="sous_domaine" || motlu=="sous_zone")
        {
          is >> nom_sous_domaine;
        }
      else if (motlu == accolade_fermee)
        break;
      else
        {
          Cerr << "Error in Ch_front_input_uniforme::readOn" << finl;
          Cerr << "We expected a } instead of " << motlu << finl;
          Process::exit();
        }
      is >> motlu;
    }
  if (nom_champ=="??" || nom_pb=="??" || nbc==-1)
    {
      Cerr << "Error in Ch_front_input_uniforme::readOn" << finl;
      Cerr << "It must be define the name, problem and number of components" << motlu << finl;
      Process::exit();
    }
  // Name the field
  set_name(nom_champ);

  // Find the problem and register in its list
  mon_pb=ref_cast(Probleme_base,Interprete::objet(nom_pb));

  mon_sous_domaine.reset();
  if (nom_sous_domaine!="??")
    {
      if (!sous_domaine_ok)
        {
          Cerr << "Error in Ch_front_input_uniforme::readOn" << finl;
          Cerr << "The sous_domaines are not managed by this field" << finl;
          Process::exit();
        }
      mon_sous_domaine=mon_pb->domaine().ss_domaine(nom_sous_domaine);
    }
}

void Champ_Input_Proto::setValueOnTab(const TrioField& afield, DoubleTab& tab)
{
  if (afield._type==0)
    {
      if ((afield._nb_elems!=tab.dimension(0))&&(afield._nb_elems!=tab.dimension_tot(0)))
        {

          throw WrongArgument(mon_pb->le_nom().getChar(),"setInputField",afield.getCharName(),"should have the same _nb_elems as returned by getInputFieldTemplate");
        }
    }
  else
    {
      if (afield._nbnodes!=tab.dimension(0))
        throw WrongArgument(mon_pb->le_nom().getChar(),"setInputField","afield","should have the same _nb_elems as returned by getInputFieldTemplate");
      if (afield._nb_field_components!=tab.dimension(1))

        throw WrongArgument(mon_pb->le_nom().getChar(),"setInputField","afield","should have the same _nb_field_components as returned by getInputFieldTemplate");
    }
  // order is   c1 c2 c3   c1 c2 c3   c1 c2 c3
  memcpy(tab.addr(),afield._field,tab.size()*sizeof(double));
  tab.echange_espace_virtuel();
}

void Champ_Input_Proto::setInitialValue(const ArrOfDouble& vo)
{
  initial_value_=vo;
}

void Champ_Input_Proto::setDoubleValue(const double val)
{
  Cerr<< "You need to override setDoubleValue() in derived class of Champ_Input_Proto"<< finl;
  throw;
}
