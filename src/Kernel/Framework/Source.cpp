/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Discretisation_base.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <TRUSTTabs.h>
#include <Source.h>
#include <Perf_counters.h>

Implemente_instanciable(Source,"Source",OWN_PTR(Source_base));

Sortie& Source::printOn(Sortie& os) const
{
  return OWN_PTR(Source_base)::printOn(os);
}

Entree& Source::readOn(Entree& is)
{
  return OWN_PTR(Source_base)::readOn(is);
}


/*! @brief
 *
 * @param (Nom& typ) the type name to assign to the source
 */
void Source::typer_direct(const Nom& typ)
{
  OWN_PTR(Source_base)::typer(typ);
}


/*! @brief Types the source by computing the required type name from the provided parameters.
 *
 * @param (Nom& typ) the beginning of the type name
 * @param (Equation_base& eqn) the equation associated with the source
 */
void Source::typer(const Nom& typ, const Equation_base& eqn)
{
  Nom type= eqn.discretisation().get_name_of_type_for(que_suis_je(),typ,eqn);
  if (eqn.que_suis_je() == "Transport_Marqueur_FT")
    type = typ;
  if(eqn.que_suis_je() == "Navier_Stokes_standard_sensibility" && typ=="boussinesq_temperature")
    type ="boussinesq_temperature_sensibility_VEFPreP1B_P1NC";


  Cerr << type << finl;
  //Cout << "In source.cpp source type = " << type << finl;
  OWN_PTR(Source_base)::typer(type);
}

/*! @brief Call to the underlying object.
 *
 * Adds the contribution of the source to the array passed as parameter.
 *
 * @param (DoubleTab& xx) the array to which the source contribution is added
 * @return (DoubleTab&) the modified parameter xx
 */
DoubleTab& Source::ajouter(DoubleTab& xx) const
{
  statistics().begin_count(STD_COUNTERS::source_terms,statistics().get_last_opened_counter_level()+1);
  DoubleTab& tmp = valeur().ajouter(xx);
  statistics().end_count(STD_COUNTERS::source_terms);
  return tmp;
}

/*! @brief Call to the underlying object.
 *
 * Assigns the source term to the array passed as parameter.
 *
 * @param (DoubleTab& xx) the array in which the source value is stored
 * @return (DoubleTab&) the modified parameter xx
 */
DoubleTab& Source::calculer(DoubleTab& xx) const
{
  statistics().begin_count(STD_COUNTERS::source_terms,statistics().get_last_opened_counter_level()+1);
  DoubleTab& tmp = valeur().calculer(xx);
  statistics().end_count(STD_COUNTERS::source_terms);
  return tmp;
}
