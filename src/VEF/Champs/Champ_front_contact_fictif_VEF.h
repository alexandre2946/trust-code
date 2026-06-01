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


#ifndef Champ_front_contact_fictif_VEF_included
#define Champ_front_contact_fictif_VEF_included

#include <Champ_front_contact_VEF.h>

//class Equation_base;
//class Milieu_base;
//class Domaine_dis_base;
//class Domaine_Cl_dis_base;
//class Front_dis_base;

/*! @brief class Champ_front_contact_fictif_VEF Derived class from Champ_front_contact_VEF, which itself derives from
 *
 *         Champ_front_var and represents the
 *      boundary fields obtained by taking the trace
 *      of an object of type OWN_PTR(Champ_Inc_base) (unknown field of an equation)
 *
 * @sa Champ_front_var_instationnaire Champ_Inc
 */
class Champ_front_contact_fictif_VEF : public Champ_front_contact_VEF
{

  Declare_instanciable(Champ_front_contact_fictif_VEF);

public:

  void mettre_a_jour(double temps) override;


protected :

  double conduct_fictif = -100.; // thermal conductivity of the fictitious solid
  double ep_fictif= -100.; // thickness of the fictitious solid

};

#endif
