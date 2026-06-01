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


#ifndef Source_Generique_base_included
#define Source_Generique_base_included

#include <Champ_Generique_base.h>
#include <Source_base.h>


/*! @brief Source_Generique_base This class is the base of the source term hierarchy carrying
 *
 *      an OWN_PTR(Champ_Generique_base) that allows evaluating an expression depending
 *      on fields of the problem.
 *
 * @sa Source_base, Abstract class., Abstract methods:, DoubleTab& ajouter(DoubleTab& ) const, void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base& ) [protected], void associer_pb(const Probleme_base& ) [protected], Nom localisation_source(), Syntax:, Sources { Source_Generique "Champ_Generique { ...} " }, with "Champ_Generique" a generic field to specify., Note: the discretization of the field returned by the generic field, must correspond to that where the source term is evaluated.
 */
class Source_Generique_base : public Source_base
{
  Declare_base(Source_Generique_base);

public :

  DoubleTab& calculer(DoubleTab& ) const override;
  void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base&) override =0;
  void associer_pb(const Probleme_base& ) override;

  virtual Nom localisation_source() =0;
  void completer() override;
  void mettre_a_jour(double temps) override;

protected :

  OWN_PTR(Champ_Generique_base) ch_source_;

};

#endif
