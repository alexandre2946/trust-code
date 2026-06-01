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

#ifndef Terme_Boussinesq_base_included
#define Terme_Boussinesq_base_included

#include <Convection_Diffusion_Temperature.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <TRUST_Vector.h>
#include <TRUST_Ref.h>

#include <Parser_U.h>
#include <Domaine.h>

class Convection_Diffusion_std;
class Champ_Don_base;
class Param;

/*! @brief Classe Terme_Boussinesq_base This class represents the gravity term appearing in the momentum equation
 *
 *     divided by the reference density.
 *     We are in the framework of the Boussinesq hypothesis: the density is
 *     assumed constant and equal to its reference value except in the
 *     body force term where a small variation of the density as a function
 *     of one or more scalars transported by the flow (temperature and/or
 *     one or more concentrations) is taken into account.
 *     Special case of Terme_Boussinesq_base for temperature:
 *     the gravity term has the expression: beta*(T-T0) where beta
 *     represents the thermal expansion coefficient and T0 a reference
 *     value for temperature.
 *     Special case of Terme_Boussinesq_base for concentration:
 *     the gravity term has the expression: beta*(C-C0)
 *     or (beta[0]*(C[0]-C0[0]) + .....+ beta[i]*(C[i]-C0[i])) in
 *     the case of a vector of concentrations.
 *     beta represents the variation of density as a function
 *     of the constituent concentration and C0 a reference value
 *     for the concentration.
 *
 *
 */
class Terme_Boussinesq_base : public Source_base
{
  Declare_base(Terme_Boussinesq_base);

public :

  void associer_pb(const Probleme_base& pb) override ;
  inline const Champ_Don_base& gravite() const { return la_gravite_.valeur(); }
  inline int verification() const { return verif_; }
  inline double Scalaire0(int i) const { return Scalaire0_[i]; }
  inline const Champ_Don_base& beta() const { return beta_.valeur(); }
  inline const Convection_Diffusion_std& equation_scalaire() const { return equation_scalaire_.valeur(); }
  DoubleTab& calculer(DoubleTab& resu) const override
  {
    resu=0;
    return ajouter(resu);
  }
  void mettre_a_jour(double temps) override
  {
    for (int i=0; i<Scalaire0_.size_array(); i++)
      {
        fct_Scalaire0_[i].setVar(0,temps);
        Scalaire0_[i] = fct_Scalaire0_[i].eval();
      }
  }
  inline const ArrOfDouble& getScalaire0() const { return Scalaire0_; }

protected :
  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;

  OBS_PTR(Champ_Don_base) la_gravite_;
  int verif_=1;
  ArrOfDouble Scalaire0_; // T0=Scalaire0_(0) ou C0(i)=Scalaire0_(i)
  Nom NomScalaire_; // Temperature ou Concentration
  VECT(Parser_U) fct_Scalaire0_;
  OBS_PTR(Champ_Don_base) beta_;
  OBS_PTR(Convection_Diffusion_std) equation_scalaire_;
  inline void check() const;
};

inline void Terme_Boussinesq_base::check() const
{
  // No verification other than at the first time step, if verification is not disabled and only for temperature
  if (equation_scalaire().probleme().schema_temps().nb_pas_dt()>0 || equation_scalaire().probleme().reprise_effectuee() || verif_==0 || !sub_type(Convection_Diffusion_Temperature,equation_scalaire())) return;

  // New: verify that average(T)==T0 only at the start of the computation
  const double T0 = Scalaire0(0);
  double moyenne_T = mp_moyenne_vect(equation_scalaire().inconnue().valeurs());
  if (inf_ou_egal(moyenne_T,T0-10) || sup_ou_egal(moyenne_T,T0+10))
    {
      Cerr << "New criteria in TRUST for the Boussinesq source term definition:" << finl;
      Cerr << "To avoid an incorrect choice for T0 value" << finl;
      Cerr << "the initial average temperature on the domain " << equation_scalaire().probleme().domaine().le_nom() << finl;
      Cerr << "should be between +/-10 degrees around T0 value." << finl;
      Cerr << "The initial average temperature is : " << moyenne_T << finl;
      Cerr << "T0 value is : " << T0 << finl;
      Cerr << "So, you need to change T0 or the initial temperature field to respect this criteria defined by default." << finl;
      Cerr << "If you want to overcome this criteria, which it is not recommended," << finl;
      Cerr << "you can specify, into the Boussinesq source term definition, the option:" << finl;
      Cerr << "verif_boussinesq 0" << finl;
      Process::exit();
    }
}

// Method to compute the value on a cell-centered field for a uniform or multi-component field
inline double valeur(const DoubleTab& valeurs, const int elem, const int dim)
{
  if(valeurs.nb_dim()==1)
    return valeurs(elem);
  else
    return valeurs(elem,dim);
}

// Method to compute the value on a face bounded by elem1 and elem2 for a uniform or multi-component field
inline double valeur(const DoubleTab& valeurs_champ, int elem1, int elem2, const int compo)
{
  if (valeurs_champ.dimension(0)==1)
    return valeurs_champ(0,compo); // Uniform field
  else
    {
      if (elem2<0) elem2 = elem1; // boundary face
      if (valeurs_champ.nb_dim()==1)
        return 0.5*(valeurs_champ(elem1)+valeurs_champ(elem2));
      else
        return 0.5*(valeurs_champ(elem1,compo)+valeurs_champ(elem2,compo));
    }
}
KOKKOS_INLINE_FUNCTION
double valeur(CDoubleTabView valeurs_champ, int valeurs_champ_dimension0, int nb_dim, int elem1, int elem2, const int compo, int nb_compo)
{
  if (valeurs_champ_dimension0==1)
    return valeurs_champ(compo,0); // Uniform field
  else
    {
      if (elem2 < 0) elem2 = elem1; // boundary face
      if (nb_dim == 1)
        return 0.5*(valeurs_champ(elem1,0)+valeurs_champ(elem2,0));
      else
        return 0.5*(valeurs_champ(elem1,compo)+valeurs_champ(elem2,compo));
    }
}

#endif

