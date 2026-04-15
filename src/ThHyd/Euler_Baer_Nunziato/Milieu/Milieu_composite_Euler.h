/****************************************************************************
* Copyright (c) 2022, CEA
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

#ifndef Milieu_composite_Euler_included
#define Milieu_composite_Euler_included

#include <Milieu_composite.h>

class Milieu_composite_Euler: public Milieu_composite
{
  Declare_instanciable(Milieu_composite_Euler);
public :
  int initialiser(const double temps) override {  return 1; }
  void discretiser(const Probleme_base& pb, const  Discretisation_base& dis) override;
  void mettre_a_jour(double temps) override { /* Do nothing */ }

  const Interface_base& interface_phase() const { return inter_lu_.valeur(); }
  Interface_base& interface_phase() { return inter_lu_.valeur(); }

  void init_energie_tot(DoubleTab& energie_tot_jdd) const;
  void calculer_pression(DoubleTab& pression_old) const ;
  void calculer_vitesse_son(DoubleTab& c) const;

protected :
  void mettre_a_jour_tabs() override { /* Do nothing */ }
};

#endif /* Milieu_composite_Euler_included */
