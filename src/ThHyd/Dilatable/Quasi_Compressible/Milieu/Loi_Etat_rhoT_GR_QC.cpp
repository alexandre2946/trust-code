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

#include <Loi_Etat_rhoT_GR_QC.h>
#include <Motcle.h>
#include <Fluide_Quasi_Compressible.h>
#include <Champ_Uniforme.h>

Implemente_instanciable(Loi_Etat_rhoT_GR_QC,"Loi_Etat_rhoT_Gaz_Reel_QC",Loi_Etat_GR_base);
// XD rhoT_gaz_reel_QC loi_etat_gaz_reel_base rhoT_gaz_reel_QC NO_BRACE Class for real gas state law used with a
// XD_CONT quasi-compressible fluid.
// XD attr bloc bloc_lecture bloc REQ Description.


Sortie& Loi_Etat_rhoT_GR_QC::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

Entree& Loi_Etat_rhoT_GR_QC::readOn(Entree& is)
{
  double PolyRho_lu=-1, PolyT_lu=-1, MMole_lu=-1, Pr_lu=-1;

  Motcle accferme="}";
  Motcle accouverte="{";

  Motcle motlu;
  is >> motlu;
  Cerr<<"Reading the real gas state law rhoT"<<finl;
  if (motlu != accouverte)
    {
      Cerr<<" Expected "<<accouverte<<" instead of "<<motlu<<finl;
      abort();
    }
  Motcles les_mots(4);
  {
    les_mots[0] = "Prandtl";
    les_mots[1] = "masse_molaire";
    les_mots[2] = "Poly_rho";
    les_mots[3] = "Poly_T";
    //     les_mots[4] = "Cp";   //for debugging
  }
  is >> motlu;
  while(motlu != accferme )
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :
          {
            is>>Pr_;
            Pr_lu=1;
            break;
          }
        case 1 :
          {
            is>>MMole_;
            MMole_lu=1;
            break;
          }
        case 2 :
          {
            int i,j,im,jm;
            is>>im;
            is>>jm;
            PolyRho_.resize(im,jm);
            for (i=0 ; i<im ; i++)
              for (j=0 ; j<jm ; j++)
                is>>PolyRho_(i,j);
            PolyRho_lu=1;
            break;
          }
        case 3 :
          {
            int i,j,im,jm;
            is>>im;
            is>>jm;
            PolyT_.resize(im,jm);
            for (i=0 ; i<im ; i++)
              for (j=0 ; j<jm ; j++)
                is>>PolyT_(i,j);
            PolyT_lu=1;
            break;
          }
        case 4 :
          {
            is>>Cp_;
            PolyRho_lu=1;
            PolyT_lu=1;
            MMole_lu=1;
            if (Cp_<0)
              {
                Cp_ = -Cp_;
                debug=1;
              }
            else
              debug=0;
            R = Cp_ *(1.-1./1.4);
            Cerr<<"Debogage Gaz Reel : Cp="<<Cp_<<" R="<<R<<finl;
            break;
          }
        default :
          {
            Cerr<<"An equation of state "<<que_suis_je()<<" does not have the property "<<motlu<<finl;
            Cerr<<"Expected a keyword in:"<<finl<<les_mots<<finl;
            abort();
          }
        }
      is >> motlu;
    }

  if (Pr_lu==-1)
    {
      Cerr<<"ERROR: expected the definition of the Prandtl number (Prandtl pr)"<<finl;
      abort();
    }
  if (PolyRho_lu==-1)
    {
      Cerr<<"ERROR: expected the definition of the density polynomial (Poly_Rho nb_coeff1 nb_coeff2 a00 a01 a02 ...)"<<finl;
      abort();
    }
  else
    {
      int i,j;
      Cerr<<"Polynomial for rho:"<<finl;
      for (i=0 ; i<PolyRho_.dimension(0) ; i++)
        {
          for (j=0 ; j<PolyRho_.dimension(1) ; j++)
            {
              Cerr<<"   f("<<i<<","<<j<<")= "<<PolyRho_(i,j);
            }
          Cerr<<finl;
        }
    }
  if (PolyT_lu==-1)
    {
      Cerr<<"ERROR: expected the definition of the temperature polynomial (Poly_T nb_coeff1 nb_coeff2 a00 a01 a02 ...)"<<finl;
      abort();
    }
  else
    {
      int i,j;
      Cerr<<"Polynomial for T:"<<finl;
      for (i=0 ; i<PolyT_.dimension(0) ; i++)
        {
          for (j=0 ; j<PolyT_.dimension(1) ; j++)
            {
              Cerr<<"   f("<<i<<","<<j<<")= "<<PolyT_(i,j);
            }
          Cerr<<finl;
        }
    }
  if (MMole_lu==-1)
    {
      Cerr<<"ERROR: expected the definition of the molar mass (masse_molaire m)"<<finl;
      abort();
    }
  return is;
}

/*! @brief Computes the pointwise density.
 *
 * @param P Pressure.
 * @param h Enthalpy.
 * @return Corresponding density value.
 */
double Loi_Etat_rhoT_GR_QC::calculer_masse_volumique(double P, double h) const
{
  double res = 0;
  if (R==-1)
    {
      double H = h;
      int i,j;
      for (i=0 ; i<PolyRho_.dimension(0) ; i++)
        for (j=0 ; j<PolyRho_.line_size() ; j++)
          res += PolyRho_(i,j) *pow(P,j) *pow(H,i);
    }
  else
    {
      //debug
      res = P*Cp_/(R*h);
    }
  return res;
}

/*! @brief Computes the pointwise temperature.
 *
 * @param P Pressure.
 * @param h Enthalpy.
 * @return Corresponding temperature value.
 */
double Loi_Etat_rhoT_GR_QC::calculer_temperature(double P, double h)
{
  double res = 0;
  if (R==-1)
    {
      int i,j;
      double H = h;
      for (i=0 ; i<PolyT_.dimension(0) ; i++)
        for (j=0 ; j<PolyT_.line_size() ; j++)
          res += PolyT_(i,j) *pow(P,i) *pow(H,j);
    }
  else
    {
      //debug
      res = h/Cp_;
    }
  return res;
}

/*! @brief Real gas case: recomputes enthalpy from pressure and temperature.
 *
 */
double Loi_Etat_rhoT_GR_QC::calculer_H(double Pth_, double T_) const
{
  double res=0;
  if (R==-1)
    {
      //need to solve c0-T_ + c1.H + c2.H^2 +...=0
      int i,j, it,max_iter = 1000;
      double eps = 1.e-6;
      DoubleVect coef(PolyT_.line_size());
      coef = 0;

      for (i=0 ; i<PolyT_.dimension(0) ; i++)
        for (j=0 ; j<PolyT_.line_size() ; j++)
          coef(j) += PolyT_(i,j) * pow(Pth_,i);

      coef(0) -= T_;
      double dH, H = -coef(0)/coef(1), num=coef(0), den=0;
      for (j=1 ; j<PolyT_.line_size() ; j++)
        {
          num += coef(j) * pow(H,j);
          den += j*coef(j) * pow(H,j-1);
        }
      dH = num/den;
      H -= dH;

      it=0;
      while (dH/H>eps && it<max_iter)
        {
          num=coef(0);
          den=0;
          for (j=1 ; j<PolyT_.line_size() ; j++)
            {
              num += coef(j) * pow(H,j);
              den += j*coef(j) * pow(H,j-1);
            }
          dH = num / den;
          H -= dH;
          it++;
        }

      if (dH/H>eps)
        {
          Cerr<<"Problem in Loi_Etat_rhoT_GR_QC::calculer_H: Newton resolution not possible "<<finl;
          Cerr<<"Pth= "<<Pth_<<" t="<<T_<<" H="<<H<<finl;
          abort();
        }
      res = H;
    }
  else
    {
      //debug
      res = Cp_*T_;
    }
  return res;
}

double Loi_Etat_rhoT_GR_QC::Drho_DP(double P, double h) const
{
  double res = 0;
  if (R==-1)
    {
      int i,j;
      double H =h;
      for (i=1 ; i<PolyRho_.dimension(0) ; i++)
        for (j=0 ; j<PolyRho_.line_size() ; j++)
          res += i* PolyRho_(i,j) *pow(P,i-1) *pow(H,j);
    }
  else
    {
      //debug
      res = 1/(R*h/Cp_);
    }
  return res;
}

//actually corresponds to Drho_DH
double Loi_Etat_rhoT_GR_QC::Drho_DT(double P, double h) const
{
  double res = 0;
  if (R==-1)
    {
      int i,j;
      double H = h;
      for (i=1 ; i<PolyRho_.dimension(0) ; i++)
        for (j=0 ; j<PolyRho_.line_size() ; j++)
          res += i* PolyRho_(i,j) *pow(P,j) *pow(H,i-1);
    }
  else
    {
      //debug
      res = -P*Cp_/(R*h*h);
    }
  return res;
}

double Loi_Etat_rhoT_GR_QC::DT_DH(double P, double h) const
{
  double res = 0;
  if (R==-1)
    {
      int i,j;
      double H = h;
      for (i=0 ; i<PolyT_.dimension(0) ; i++)
        for (j=1 ; j<PolyT_.line_size() ; j++)
          res += j* PolyT_(i,j) *pow(P,i) *pow(H,j-1);
    }
  else
    {
      //debug
      res = Cp_-R;
    }
  return res;
}

void Loi_Etat_rhoT_GR_QC::calculer_masse_volumique()
{
  Loi_Etat_GR_base::calculer_masse_volumique();
}
