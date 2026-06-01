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

#include <Assembleur_base.h>

Implemente_base_sans_constructeur(Assembleur_base,"Assembleur_base",Objet_U);

Assembleur_base::Assembleur_base() :
  resoudre_increment_pression_(-1),resoudre_en_u_(-1) // -1 means "not initialized"

{
}

// printOn

Sortie& Assembleur_base::printOn(Sortie& s ) const
{
  return s  << que_suis_je();
}


// readOn

Entree& Assembleur_base::readOn(Entree& is )
{
  return is;
}

/*! @brief Sets the value of the resoudre_increment_pression_ flag.
 *
 * This flag determines whether the pressure solver calculates an increment of
 *  pressure between two time steps or the total pressure (this affects
 *  in particular the calculation of boundary conditions).
 *  flag = 0: pressure resolution
 *  flag = 1: pressure increment resolution
 *
 */
int Assembleur_base::set_resoudre_increment_pression(int flag)
{
  assert(flag == 0 || flag == 1);
  resoudre_increment_pression_ = flag;
  return flag;
}

/*! @brief Returns the value of the resoudre_increment_pression_ flag (0 or 1) Returns -1 if the flag has not been initialized
 *
 */
int Assembleur_base::get_resoudre_increment_pression() const
{
  return resoudre_increment_pression_;
}

/*! @brief Sets the value of the resoudre_en_u__ flag.
 *
 * This flag determines whether the pressure solver solves for u or for rho*u
 *   (this affects
 *  in particular the calculation of Dirichlet boundary conditions).
 *  flag = 1: u resolution
 *  flag = 0: rho_u resolution
 *
 */
int Assembleur_base::set_resoudre_en_u(int flag)
{
  assert(flag == 0 || flag == 1);
  resoudre_en_u_ = flag;
  return flag;
}

/*! @brief Returns the value of the resoudre_en_u_ flag (0 or 1) Returns -1 if the flag has not been initialized
 *
 */
int Assembleur_base::get_resoudre_en_u() const
{
  return resoudre_en_u_;
}

/*! @brief Assembly of the matrix div(porosity/rho * grad P) The type of the "rho" field to be provided depends on the discretization (vdf, vef, p1b, .
 *
 * ..)
 *  For front-tracking for example.
 *
 */
int Assembleur_base::assembler_rho_variable(Matrice& mat, const Champ_Don_base& rho)
{
  Cerr << "Error in Assembleur_base::assembler_rho_variable\n";
  Cerr << "This method has not been implemented for ";
  Cerr << que_suis_je() << finl;
  assert(0);
  exit();
  return 0;
}

int Assembleur_base::assembler_QC(const DoubleTab&, Matrice& matrice)
{
  Cerr << "Stop in Assembleur_base::assembler_QC." << finl;
  assert(0);
  exit();
  return 1;
}

void Assembleur_base::completer(const Equation_base& )
{
  ;
}
