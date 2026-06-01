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

#ifndef MorEqn_included
#define MorEqn_included

#include <TRUST_Ref.h>
#include <Motcle.h>

class Champ_base;
class Equation_base;

/*! @brief class MorEqn Class that groups the functionalities of linking with an
 *
 *      Equation. The member MorEqn::mon_equation contains a reference
 *      to an Equation_base object. The classes that are "Equation pieces"
 *      are associated by a 1-1 link to their equation by inheriting from MorEqn,
 *      for example: Operateur_base, Solveur_Masse, Source_base ...
 *
 * @sa Equation_base, This is not a class of the TRUST tree by itself., This class is intended to be a base class of a class, which will also inherit from Objet_U
 */
class MorEqn
{

public:
  void associer_eqn(const Equation_base&);
  virtual void calculer_pour_post(Champ_base& espace_stockage,const Nom& option, int comp) const;
  virtual Motcle get_localisation_pour_post(const Nom& option) const;
  inline const Equation_base& equation() const;
  inline Equation_base& equation();

  inline int mon_equation_non_nul() const;

  /* compatibility with multiphase equations : by default, error message */
  virtual void check_multiphase_compatibility() const;
protected :
  OBS_PTR(Equation_base) mon_equation;
  inline virtual ~MorEqn();
};

MorEqn::~MorEqn()
{}
/*! @brief Returns the reference to the equation pointed to by MorEqn::mon_equation.
 *
 *     (const version)
 *
 * @return (Equation_base&) the equation associated with the object
 * @throws no associated equation
 */
inline const Equation_base& MorEqn::equation() const
{
  if (!mon_equation)
    {
      Cerr << "\nError in MorEqn::equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return mon_equation.valeur();
}

/*! @brief Returns the reference to the equation pointed to by MorEqn::mon_equation.
 *
 * @return (Equation_base&) the equation associated with the object
 */
inline  Equation_base& MorEqn::equation()
{
  if (!mon_equation)
    {
      Cerr << "\nError in MorEqn::equation() : The equation is unknown !" << finl;
      Process::exit();
    }
  return mon_equation.valeur();
}
int MorEqn::mon_equation_non_nul() const
{
  return bool(mon_equation);
}
#endif

