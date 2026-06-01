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

#ifndef Perte_Charge_VEF_included
#define Perte_Charge_VEF_included

#include <Terme_Source_Qdm.h>
#include <Source_base.h>
#include <TRUST_Ref.h>
#include <Parser_U.h>
#include <Champ_Don_base.h>

class Sous_domaine_VF;
class Domaine_Cl_VEF;
class Champ_Inc_base;
class Domaine_VEF;
class Fluide_base;
class Param;

//! Factorizes the common functionality of several pressure drop operators in VEF, velocity at faces.
/**
   Perte_Charge_Isotrope, Perte_Charge_Directionnelle and
   Perte_Charge_Anisotrope inherit from Perte_Charge_VEF. They
   must essentially override readOn() and perte_charge().
   readOn() is expected to read at least diam_hydr and sous_domaine.

   These classes are intended to replace Perte_Charge_VEF_Face
   and Perte_Charge_VEF_P1NC.
*/

class Perte_Charge_VEF : public Source_base, public Terme_Source_Qdm
{
  Declare_base(Perte_Charge_VEF);

public:
  DoubleTab& ajouter(DoubleTab& ) const override; //!< Calls perte_charge for each face where needed
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override ;
  DoubleTab& calculer(DoubleTab& ) const override ;
  void associer_pb(const Probleme_base&) override;  //!< Associates le_fluide and la_vitesse
  void completer() override;

protected:
  virtual void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void associer_domaines(const Domaine_dis_base&,const Domaine_Cl_dis_base&) override;   //!< Associates le_dom_VEF and le_dom_Cl_VEF

  //! Called for each face by ajouter()
  /**
     Uses the intermediate quantities: u, norme_u, dh_valeur, reynolds.
     Returns the computed result in p_charge.

     Advantage: good factorization.
     Drawback: cost of a virtual method call inside a loop.

     \param u Velocity at the current face. 2 or 3 components.
     \param pos position of the current face.
     \param t current time.
     \param norme_u Norm of the velocity at the current face.
     \param dh Hydraulic diameter at the current face (derived from diam_hydr).
     \param nu kinematic viscosity.
     \param reynolds Reynolds number at the current face: norme_u * dh_valeur / nu.
     \param coeff_ortho coefficient in the orthogonal direction.
     \param coeff_long coefficient in the longitudinal direction.
     \param u_l velocity in the longitudinal direction.
     \param v_valeur longitudinal direction of p_charge with 2 or 3 components.
     The pressure drop equals -coeff_long*u_l*v_valeur -coeff_ortho(u -u_v*v_valeur).
  */
  virtual void coeffs_perte_charge(const DoubleVect& u, const DoubleVect& pos,
                                   double t, double norme_u, double dh, double nu, double reynolds,
                                   double& coeff_ortho, double& coeff_long,
                                   double& u_l, DoubleVect& v_valeur) const=0;

  //! Hydraulic diameter used in the pressure drop computation
  OWN_PTR(Champ_Don_base) diam_hydr;
  //! Fluid associated with the problem
  OBS_PTR(Fluide_base) le_fluide;
  //! Velocity associated with the solved equation
  OBS_PTR(Champ_Inc_base) la_vitesse;
  //! Domain to which the pressure drop applies.
  OBS_PTR(Domaine_VEF) le_dom_VEF;
  OBS_PTR(Domaine_Cl_VEF) le_dom_Cl_VEF;

  // Sub-domain case
  bool sous_domaine=false; //!< Is the term restricted to a sub-domain?
  Nom nom_sous_domaine; //!< Name of the sub-domain, initialized in readOn()
  OBS_PTR(Sous_domaine_VF) le_sous_domaine_dis; //!< Initialized in completer()
  int implicite_;

  mutable Parser_U lambda;
};

#endif
