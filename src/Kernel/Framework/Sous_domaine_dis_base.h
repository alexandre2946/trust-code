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

#ifndef Sous_domaine_dis_base_included
#define Sous_domaine_dis_base_included

#include <Domaine_forward.h>
#include <TRUST_Vector.h>
#include <TRUST_Deriv.h>
#include <TRUST_Ref.h>

class Domaine_dis_base;

/*! @brief This class is at the base of the discretized sub-domain hierarchy.
 * At the time of domain discretization, just as a
 * Domaine_dis is created for each domain, a Sous_domaine_dis_base is
 * created for each Sous_Domaine and discretized.
 *
 * Classes inheriting from Sous_domaine_dis_base must override the discretiser() method;
 *
 * The Sous_Domaine defines a subset of the elements of the Domaine.
 * Similarly, the Sous_domaine_dis_base defines a subset of the attributes
 * of Domaine_dis, in particular of faces.
 *
 * @sa Domaine_dis, Sous_Domaine, Domaine_dis_base, Sous_Domaine_VF
*/
class Sous_domaine_dis_base : public Objet_U
{
  Declare_base(Sous_domaine_dis_base);

public:
  void associer_sous_domaine(const Sous_Domaine&);
  void associer_domaine_dis(const Domaine_dis_base&);

  // REF accessor methods
  const Sous_Domaine& sous_domaine() const;
  Sous_Domaine& sous_domaine();
  const Domaine_dis_base& domaine_dis() const;
  Domaine_dis_base& domaine_dis();

  // Specific methods:
  virtual void discretiser() = 0;

protected:
  OBS_PTR(Sous_Domaine) le_sous_domaine;
  OBS_PTR(Domaine_dis_base) le_dom_dis;
};

#endif /* Sous_domaine_dis_base_included */
