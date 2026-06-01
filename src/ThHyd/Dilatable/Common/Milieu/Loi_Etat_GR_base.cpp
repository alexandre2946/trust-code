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

#include <Fluide_Quasi_Compressible.h>
#include <Loi_Etat_GR_base.h>
#include <Champ_Uniforme.h>
#include <Motcle.h>

Implemente_base_sans_constructeur(Loi_Etat_GR_base,"Loi_Etat_Gaz_Reel_base",Loi_Etat_base);
// XD loi_etat_gaz_reel_base loi_etat_base loi_etat_gaz_reel_base INHERITS_BRACE Basic class for real gases state laws
// XD_CONT used with a dilatable fluid.

Loi_Etat_GR_base::Loi_Etat_GR_base() : MMole_(-1),Cp_(-1),R(-1) { }

Sortie& Loi_Etat_GR_base::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

Entree& Loi_Etat_GR_base::readOn(Entree& is)
{
  return is;
}

/*! @brief Returns the type of fluid associated.
 *
 * @return The fluid type name ("Gaz_Reel").
 */
const Nom Loi_Etat_GR_base::type_fluide() const
{
  return "Gaz_Reel";
}

/*! @brief Initialises the state law by computing Pth.
 *
 */
void Loi_Etat_GR_base::initialiser()
{
  le_fluide->inco_chaleur().nommer("enthalpie");

  const DoubleTab& tab_H = le_fluide->inco_chaleur().valeurs();
  const DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  int i, n = tab_H.dimension(0);
  DoubleTab& tab_T = temperature_->valeurs();
  tab_TempC.resize(n);
  double Pth = le_fluide->pression_th();
  for (i=0 ; i<n ; i++)
    {
      tab_rho_n(i) = tab_rho(i,0);
      tab_TempC(i) = calculer_temperature(Pth,tab_H(i,0));
      tab_T(i) = tab_TempC(i);
    }
  tab_Cp.ref(le_fluide->capacite_calorifique().valeurs());
  calculer_Cp();
}

/*! @brief Initialises the enthalpy unknown.
 *
 */
void Loi_Etat_GR_base::initialiser_inco_ch()
{
  /*
   * XXX XXX XXX
   * The heat unknown inco_chaleur is first filled with the initial temperature.
   * It must be converted to enthalpy;
   * the data are the density and the temperature,
   * so enthalpy and pressure must be computed.
   */

  DoubleTab& tab_TH = le_fluide->inco_chaleur().valeurs();
  double Pth = le_fluide->pression_th();
  DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  tab_rho_n=tab_rho;
  tab_rho_np1=tab_rho;
  int som,n = tab_TH.dimension(0);
  if (le_fluide->inco_chaleur().le_nom() == "enthalpie")
    {
      for (som=0 ; som<n ; som++)
        tab_rho_np1(som) = tab_rho(som,0) = tab_rho_n(som) = calculer_masse_volumique(Pth,tab_TH(som,0));
    }
  else
    {
      for (som=0 ; som<n ; som++)
        {
          tab_TH(som,0) = calculer_H(Pth,tab_TH(som,0));
          tab_rho_np1(som) = tab_rho(som,0) = tab_rho_n(som) = calculer_masse_volumique(Pth,tab_TH(som,0));
        }
    }
  Cerr<<"FIN Loi_Etat_GR_base::initialiser_H Pth = "<<Pth<<"  H = "<<tab_TH(0,0)<<finl;
}

/*! @brief Fills the temperature array from the enthalpy unknown.
 *
 */
void Loi_Etat_GR_base::remplir_T()
{
  const DoubleTab& tab_H = le_fluide->inco_chaleur().valeurs();
  int i, n = tab_TempC.dimension(0);
  DoubleTab& tab_T = temperature_->valeurs();
  double Pth = le_fluide->pression_th();
  for (i=0 ; i<n ; i++)
    {
      tab_TempC(i) = calculer_temperature(Pth,tab_H(i,0));
      tab_T(i,0) = tab_TempC(i);
    }
}

/*! @brief Computes Cp as a function of physical quantities P and h.
 *         Cp = dh/dT = de/dT - 1/rho^2 * drho/dT
 *
 * @param P Pressure.
 * @param h Enthalpy.
 * @return Computed Cp value.
 */
double Loi_Etat_GR_base::Cp_calc(double P, double h) const
{
  double res;
  if (R==-1) res = 1./(DT_DH(P,h));
  else res = Cp_;

  return res;
}

/*! @brief Computes Cp using the PolyCp_ polynomial.
 *
 */
void Loi_Etat_GR_base::calculer_Cp()
{
  double Pth = le_fluide->pression_th();
  const DoubleTab& tab_h = le_fluide->inco_chaleur().valeurs();
  for (int i=0; i<tab_Cp.size(); i++) tab_Cp(i) = Cp_calc(Pth,tab_h(i,0));
}

/*! @brief Computes the conductivity divided by Cp: equivalent to k*dT/dh for using enthalpy in the diffusion operator.
 *
 */
void Loi_Etat_GR_base::calculer_lambda()
{
  const Champ_Don_base& mu = le_fluide->viscosite_dynamique();
  const DoubleTab& tab_mu = mu.valeurs();
  Champ_Don_base& lambda = le_fluide->conductivite();
  DoubleTab& tab_lambda = lambda.valeurs();

  int i, n=tab_lambda.size();
  if (!sub_type(Champ_Uniforme,lambda))
    {
      if (sub_type(Champ_Uniforme,mu))
        {
          for (i=0 ; i<n ; i++) tab_lambda(i,0) = tab_mu(0,0) * tab_Cp(i) / Pr_;
        }
      else
        {
          for (i=0 ; i<n ; i++) tab_lambda(i,0) = tab_mu(i,0) * tab_Cp(i) / Pr_;
        }
    }
}

double Loi_Etat_GR_base::De_DP(double P, double T) const
{
  double res = 0;
  Cerr<<"Loi_Etat_GR_base::De_DP non accede normalement"<<finl;
  abort();
  return res;
}
double Loi_Etat_GR_base::De_DT(double P, double T) const
{
  double res = 0;
  Cerr<<"Loi_Etat_GR_base::De_DP non accede normalement"<<finl;
  abort();
  return res;
}

/*! @brief Computes the thermodynamic pressure from enthalpy and density by Newton iteration.
 *
 * @param H Enthalpy.
 * @param rho Density.
 * @return Thermodynamic pressure Pth.
 */
double Loi_Etat_GR_base::inverser_Pth(double H, double rho)
{
  double P = le_fluide->pression_th();
  double acc = (calculer_masse_volumique(P,H) - rho) / Drho_DP(P,H);
  int i=0;
  while (std::fabs(acc)>1e-8 && i<1000)
    {
      P = P-acc;
      acc = (calculer_masse_volumique(P,H) - rho) / Drho_DP(P,H);
      i++;
    }
  if (std::fabs(acc)>1e-8)
    {
      Cerr<<"Probleme dans l'inversion de la pression : nb_iter="<<i<<finl;
      Cerr<<" Pth="<<P<<" H="<<H<<" rho="<<rho<<finl;
      abort();
    }
  return P;
}

void Loi_Etat_GR_base::calculer_masse_volumique()
{
  Loi_Etat_base::calculer_masse_volumique();
}
