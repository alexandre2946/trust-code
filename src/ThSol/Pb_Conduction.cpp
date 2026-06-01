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

#include <Pb_Conduction.h>

Implemente_instanciable(Pb_Conduction,"Pb_Conduction",Probleme_base);
// XD Pb_Conduction Pb_base Pb_Conduction INHERITS_BRACE Resolution of the heat equation.
// XD attr solide solide solide OPT The medium associated with the problem.
// XD attr Conduction Conduction Conduction OPT Heat equation.

/*! @brief Does nothing.
 *
 * @param s An output stream.
 * @return The output stream.
 */
Sortie& Pb_Conduction::printOn(Sortie& s ) const
{
  return s;
}


/*! @brief Simple call to Probleme_base::readOn(Entree&).
 *
 * @param is An input stream.
 * @return The modified input stream.
 */
Entree& Pb_Conduction::readOn(Entree& is )
{
  return Probleme_base::readOn(is);
}


/*! @brief Returns the number of equations in the problem.
 *
 * Always equal to 1 for a standard conduction problem.
 *
 * @return Number of equations in the problem.
 */
int Pb_Conduction::nombre_d_equations() const
{
  return 1;
}

/*! @brief Returns the Conduction equation when i = 0 (const version).
 *
 * Triggers an error otherwise because the problem has only one equation.
 *
 * @param i Index of the equation to return.
 * @return The Conduction equation.
 */
const Equation_base& Pb_Conduction::equation(int i) const
{
  assert (i==0);
  return eq_conduction;
}

/*! @brief Returns the Conduction equation when i = 0.
 *
 * Triggers an error otherwise because the problem has only one equation.
 *
 * @param i Index of the equation to return.
 * @return The Conduction equation.
 */
Equation_base& Pb_Conduction::equation(int i)
{
  assert (i==0);
  return eq_conduction;
}




