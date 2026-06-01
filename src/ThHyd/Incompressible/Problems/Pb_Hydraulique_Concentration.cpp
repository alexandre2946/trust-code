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

#include <Pb_Hydraulique_Concentration.h>
#include <Fluide_Incompressible.h>
#include <Constituant.h>
#include <Verif_Cl.h>

Implemente_instanciable(Pb_Hydraulique_Concentration,"Pb_Hydraulique_Concentration",Pb_Fluide_base);
// XD pb_hydraulique_concentration Pb_base pb_hydraulique_concentration INHERITS_BRACE Resolution of
// XD_CONT Navier-Stokes/multiple constituent transport equations.
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible REQ The fluid medium associated with the
// XD_CONT problem.
// XD attr constituant constituant constituant OPT Constituents.
// XD attr navier_stokes_standard navier_stokes_standard navier_stokes_standard OPT Navier-Stokes equations.
// XD attr convection_diffusion_concentration convection_diffusion_concentration convection_diffusion_concentration OPT
// XD_CONT Constituent transport vectorial equation (concentration diffusion convection).

/*! @brief Simple call to: Pb_Fluide_base::printOn(Sortie&) Writes the problem to an output stream.
 *
 * @param os an output stream
 * @return the modified output stream
 */
Sortie& Pb_Hydraulique_Concentration::printOn(Sortie& os) const
{
  return Pb_Fluide_base::printOn(os);
}


/*! @brief Simple call to: Pb_Fluide_base::readOn(Entree&) Reads the problem from an input stream.
 *
 * @param is an input stream
 * @return the modified input stream
 */
Entree& Pb_Hydraulique_Concentration::readOn(Entree& is)
{
  return Pb_Fluide_base::readOn(is);
}

/*! @brief Returns the number of equations. Returns 2 since a hydraulic problem with transport has 2 equations:
 *
 *       - the Navier-Stokes equation
 *       - a convection-diffusion equation (possibly vectorial)
 *
 * @return the number of equations
 */
int Pb_Hydraulique_Concentration::nombre_d_equations() const
{
  return 2;
}

/*! @brief Returns the hydraulic equation of type Navier_Stokes_std if i=0, returns the convection-diffusion equation of type
 *
 *     Convection_Diffusion_Concentration if i=1
 *     (the convection-diffusion equation may be vectorial)
 *     (const version)
 *
 * @param i the index of the equation to return
 * @return the equation corresponding to the index
 */
const Equation_base& Pb_Hydraulique_Concentration::equation(int i) const
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Hydraulique_Concentration::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_concentration;
}

/*! @brief Returns the hydraulic equation of type Navier_Stokes_std if i=0, returns the convection-diffusion equation of type
 *
 *     Convection_Diffusion_Concentration if i=1
 *     (the convection-diffusion equation may be vectorial)
 *
 * @param i the index of the equation to return
 * @return the equation corresponding to the index
 */
Equation_base& Pb_Hydraulique_Concentration::equation(int i)
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Hydraulique_Concentration::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_concentration;
}


/*! @brief Associates a medium to the problem. If the medium is of type:
 *
 *       - Fluide_Incompressible, it will be associated with the hydraulic equation
 *       - Constituant, it will be associated with the convection-diffusion equation
 *     Any other medium type causes an error.
 *
 * @param mil the physical medium to associate with the problem
 * @throws wrong type of physical medium
 */
void Pb_Hydraulique_Concentration::associer_milieu_base(const Milieu_base& mil)
{
  if ( sub_type(Fluide_Incompressible,mil) )
    eq_hydraulique.associer_milieu_base(mil);
  else if ( sub_type(Constituant,mil) )
    eq_concentration.associer_milieu_base(mil);
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je() << " ne peut etre associe a " << finl;
      Cerr << "un probleme de type Pb_Hydraulique_Concentration " << finl;
      exit();
    }
}


/*! @brief Teste la compatibilite des equations de convection-diffusion et de l'hydraulique.
 *
 * Le test se fait sur les conditions
 *     aux limites discretisees de chaque equation.
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_concentration(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return (int) code de retour propage
 */
int Pb_Hydraulique_Concentration::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_co = eq_concentration.domaine_Cl_dis();
  return tester_compatibilite_hydr_concentration(domaine_Cl_hydr,domaine_Cl_co);
}
