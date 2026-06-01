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

#ifndef Probleme_Couple_included
#define Probleme_Couple_included

#include <TRUSTTabs_forward.h>
#include <Schema_Temps_base.h>
#include <Couplage_U.h>
#include <TRUST_Ref.h>

class Discretisation_base;
class Schema_Temps_base;
class Champ_base;

/*! @brief Probleme_Couple This is the historical coupling class of TRUST.
 *
 *      It is a particular Couplage_U.
 *      Probleme_Couple only couples Probleme_base objects, all
 *      associated with the same time scheme.
 *      It has its own field exchange mechanism via "raccords".
 *      Problems are grouped. Each group performs one iteration, then the
 *      new fields become available for subsequent groups. Within a
 *      group, the old fields are exchanged.
 *
 *
 * @sa Probleme_base Probleme
 */

//  WEC :
//  The time scheme is cloned as many times as there are problems.
//  This cloning is not very clean:
//  * they all have the same name
//  To avoid this:
//  change the syntax of the .data files to have truly one scheme per
//  problem

class Probleme_Couple : public Couplage_U
{
  Declare_instanciable(Probleme_Couple);
public :

  ///////////////////////////////////////////////
  //                                           //
  // Implementation of the Problem interface  //
  //                                           //
  ///////////////////////////////////////////////

  // interface UnsteadyProblem
  bool initTimeStep(double dt) override;
  double computeTimeStep(bool& stop) const override;
  bool solveTimeStep() override;

  // interface IterativeUnsteadyProblem

  bool iterateTimeStep(bool& converged) override;

  ////////////////////////////////////////////////////////
  //                                                    //
  // End of Problem interface implementation  //
  //                                                    //
  ////////////////////////////////////////////////////////

  bool updateGivenFields() override;

  void ajouter(Probleme_base&);
  int associer_(Objet_U&) override;
  virtual void associer_sch_tps_base(Schema_Temps_base&);
  virtual const Schema_Temps_base& schema_temps() const;
  virtual Schema_Temps_base& schema_temps();

  virtual void discretiser(Discretisation_base&);
  inline virtual void mettre_a_jour_modele_rayo(double temps);
  void initialize() override;
  void sauver() const override;

protected:

  // Definitions of problem groups.
  // Integer vector where each entry specifies the size of a group.
  // Empty vector = single group.
  ArrOfInt groupes;

  // List of cloned time schemes
  VECT(OWN_PTR(Schema_Temps_base)) sch_clones;
};

inline void Probleme_Couple::mettre_a_jour_modele_rayo(double temps)
{

  Cerr<<"The method Probleme_Couple::mettre_a_jour_modele_rayo does nothing"<<finl;
  Cerr<<"We should not pass through here"<<finl;
  exit();
}

class Probleme_Couple_Point_Fixe : public Probleme_Couple
{
  Declare_instanciable(Probleme_Couple_Point_Fixe);
public:
  bool solveTimeStep() override;
};

#endif
