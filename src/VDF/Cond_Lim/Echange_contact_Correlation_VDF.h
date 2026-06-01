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

#ifndef Echange_contact_Correlation_VDF_included
#define Echange_contact_Correlation_VDF_included

#include <Echange_global_impose.h>
#include <Domaine_forward.h>
#include <Parser_U.h>

class Milieu_base;
class Front_VF;
class Domaine_VDF;
class Param;

class Echange_contact_Correlation_VDF: public Echange_global_impose
{

  Declare_instanciable(Echange_contact_Correlation_VDF);

public:

  void mettre_a_jour(double) override;
  void completer() override;
  virtual void imprimer(double) const;
  virtual int limpr(double, double) const;

  /*
   * Method that can be overridden in derived classes to define a particular Nusselt number.
   * Otherwise, a Nusselt function of local Reynolds and local Prandtl can be entered directly in the data file.
   * It returns the value of the heat exchange coefficient on cell i.
   */
  virtual double calculer_coefficient_echange(int i);

  /*
   * Method that can be overridden in derived classes to define a particular shape.
   * Otherwise, a function of the lateral surface area S and the hydraulic diameter Dh can be entered in the data file.
   * It computes the volume of a fluid slice with lateral surface area s and hydraulic diameter d.
   * This function depends on the geometry for which the correlation is written
   * (parallel plates, cylinders, etc.)
   */
  virtual double volume(double s, double d);

  /*
   * Returns U and T and the physical properties on 1D cell i.
   */
  inline double getU(int i) const { return U(i); }
  inline double getT(int i) const { return T(i); }
  inline double getMu(int i) const { return mu(i); }
  inline double getLambda(int i) const { return lambda(i); }
  inline double getRho(int i) const { return rho(i); }
  inline double getCp() const { return Cp; }
  inline double getDh() const { return diam; }
  inline double getQh() const { return debit; }

protected:
  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;

  void calculer_CL();
  void calculer_Q();
  void calculer_prop_physique();
  void trier_coord();
  void calculer_Vitesse();
  void calculer_Tfluide();
  void calculer_h_mon_pb(DoubleTab&);
  void calculer_h_solide(DoubleTab&, const Equation_base&, const Domaine_VDF&, const Front_VF&, const Milieu_base&);
  void init_tab_echange();

  IntVect correspondance_solide_fluide;
  DoubleTab autre_h;
  bool Reprise_temperature = false;

  DoubleTab tab_ech;

  double Tinf = -100., Tsup = -100.; // Inlet and outlet temperatures
  double T_CL0 = -100., T_CL1 = -100.; // BCs on the domain. Sequential: = Tinf and Tsup; Parallel: = Tvoisin
  int dir = -1; // Direction of the cylinder
  Parser_U lambda_T;
  Parser_U mu_T;
  Parser_U rho_T;
  Parser_U fct_Nu;
  Parser_U fct_vol;
  double dt_impr = -100.;
  double Cp = -100.;
  double debit = -100.;

  double diam = -100.; // hydraulic diameter
  DoubleVect vol; // slice volumes Sdz along the flow direction
  DoubleVect coord; // coordinates of the 1D discretization points

  int N = -1;   // number of 1D mesh points for velocity and temperature
  DoubleVect U, T; // fluid unknowns: velocity and temperature
  DoubleVect Qvol; // volumetric power in the 1D fluid
  DoubleVect rho, mu, lambda; // array of fluid physical properties at a given instant
  DoubleVect h_correlation; // heat transfer coefficient from correlation

};

#endif
