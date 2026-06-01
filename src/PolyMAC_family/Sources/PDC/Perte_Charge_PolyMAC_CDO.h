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

#ifndef Perte_Charge_PolyMAC_CDO_included
#define Perte_Charge_PolyMAC_CDO_included

#include <Perte_Charge_Gen.h>
#include <Domaine_PolyMAC_CDO.h>

//! Factorises the features of several pressure drop classes in VEF, face velocity
/**
   Perte_Charge_Isotrope, Perte_Charge_Directionnelle and
   Perte_Charge_Anisotrope inherit from Perte_Charge_PolyMAC_CDO. They
   must override essentially readOn() and perte_charge().
   readOn() is expected to read at least diam_hydr and sous_domaine.

   These classes are intended to replace Perte_Charge_PolyMAC_CDO_Face
   and Perte_Charge_PolyMAC_CDO_P1NC.
*/

class Perte_Charge_PolyMAC_CDO : public Perte_Charge_Gen
{
  Declare_base(Perte_Charge_PolyMAC_CDO);
public:
  DoubleTab& ajouter(DoubleTab& ) const override; //!< Calls perte_charge for each face where it is necessary
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override ;

  const Domaine_PolyMAC_CDO& le_dom_poly()  const  { return ref_cast(Domaine_PolyMAC_CDO, le_dom_vf_.valeur()); }

};

#endif /* Perte_Charge_PolyMAC_CDO_included */
