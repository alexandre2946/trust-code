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


#ifndef Champ_front_var_instationnaire_included
#define Champ_front_var_instationnaire_included

#include <Champ_front_var.h>


/*! @brief class Champ_front_var_instationnaire Derived class from Champ_front_var that represents fields on
 *
 *      boundaries variable in space and in time.
 *      Champ_front_var_instationnaire fields are classified depending on whether
 *      their values depend or not on parameters external to
 *      the equation, as Champ_front_var_instationnaire_indep and
 *      Champ_front_var_instationnaire_dep.
 *      The implementation of the field computation must be done in the
 *      mettre_a_jour method.
 *      In the first case (indep), the initialiser method may
 *      call the mettre_a_jour method, but not in the second case
 *      (dep). In either case, it can use the unknown passed
 *      as a parameter as a first estimate.
 *
 * @sa Champ_front_var
 */
class Champ_front_var_instationnaire : public Champ_front_var
{

  Declare_base(Champ_front_var_instationnaire);

public:
  void fixer_nb_valeurs_temporelles(int nb_cases) override;
  int initialiser(double temps, const Champ_Inc_base& inco) override;
  bool has_valeurs_au_temps(double temps) const override ;
  DoubleTab& valeurs_au_temps(double temps) override;
  const DoubleTab& valeurs_au_temps(double temps) const override;
  int avancer(double temps) override;
  int reculer(double temps) override;

  virtual double valeur_au_temps_et_au_point(double temps,int som,double x,double y, double z,int comp) const
  {
    Cerr<<"valeur_au_temps_et_au_point mut be defined si valeur_au_temps_et_au_point_disponible()==1"<<finl;
    exit();
    return 0;
  };
  inline virtual int valeur_au_temps_et_au_point_disponible() const
  {
    return 0;
  };
protected :
};

#endif
