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

#ifndef Perte_Charge_Circulaire_PolyMAC_CDO_Face_included
#define Perte_Charge_Circulaire_PolyMAC_CDO_Face_included

#include <Perte_Charge_PolyMAC_HFV.h>
#include <PDC_PolyMAC_CDO_impl.h>

//!  Anisotropic pressure drop (along a unit vector v and in the plane orthogonal to this vector)
/**
   Corresponds to the sum of a pressure drop of intensity lambda in
   the direction of v, and a pressure drop of intensity
   lambda_ortho for every direction orthogonal to v.

   Or equivalently, the sum of an isotropic pressure drop of
   intensity lambda_ortho and a directional pressure drop
   of vector v and intensity lambda - lambda_ortho. This is how
   the computation is coded.

   du/dt =
   - lambda_ortho(Re,x,y,z,t) * u * ||u|| / 2 Dh
   - (lambda(Re,x,y,z,t)-lambda_ortho(Re,x,y,z,t))
   * u.v * v * ||u|| / 2 Dh ||v||^2

   Reading arguments:

   Perte_Charge_Circulaire_PolyMAC_CDO_Face diametre_hydraulique {
   lambda expression(Re,x,y,z,t)
   lambda_ortho expression(Re,x,y,z,t)
   diam_hydr champ_don
   direction champ_don
   [sous_domaine nom]
   }

*/

class Perte_Charge_Circulaire_PolyMAC_CDO_Face: public Perte_Charge_PolyMAC_CDO, public PDC_Circulaire_PolyMAC_CDO
{
  Declare_instanciable(Perte_Charge_Circulaire_PolyMAC_CDO_Face);
public:
  void mettre_a_jour(double temps) override
  {
    diam_hydr->mettre_a_jour(temps);
    diam_hydr_ortho->mettre_a_jour(temps);
    v->mettre_a_jour(temps);
  }

protected:
  void set_param(Param& titi) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  //! Implements the effective computation of the pressure drop for a given location
  void coeffs_perte_charge(const DoubleVect&, const DoubleVect&, double, double, double, double, double, double&, double&, double&, DoubleVect&) const override;
};

/////////////////////////////////////////////////

class Perte_Charge_Circulaire_PolyMAC_HFV_Face: public Perte_Charge_PolyMAC_HFV, public PDC_Circulaire_PolyMAC_CDO
{
  Declare_instanciable(Perte_Charge_Circulaire_PolyMAC_HFV_Face);
public:
  void mettre_a_jour(double temps) override
  {
    diam_hydr->mettre_a_jour(temps);
    diam_hydr_ortho->mettre_a_jour(temps);
    v->mettre_a_jour(temps);
  }

protected:
  void set_param(Param& titi) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void coeffs_perte_charge(const DoubleVect&, const DoubleVect&, double, double, double, double, double, double&, double&, double&, DoubleVect&) const override;
};

#endif /* Perte_Charge_Circulaire_PolyMAC_CDO_Face_included */
