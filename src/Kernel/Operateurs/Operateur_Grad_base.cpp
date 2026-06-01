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

#include <Operateur_Grad_base.h>
#include <TRUSTTrav.h>

Implemente_base(Operateur_Grad_base,"Operateur_Grad_base",Operateur_base);


/*! @brief DOES NOTHING To be overridden in derived classes
 *
 * @param (Sortie& os) an output stream
 * @return (Sortie&) the output stream
 */
Sortie& Operateur_Grad_base::printOn(Sortie& os) const
{
  return os;
}


/*! @brief DOES NOTHING To be overridden in derived classes
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the input stream
 */
Entree& Operateur_Grad_base::readOn(Entree& is)
{
  return is;
}

void Operateur_Grad_base::dimensionner(Matrice_Morse& mat) const
{
  /* try dimensionner_blocs(), but for pressure! */
  if (has_interface_blocs()) dimensionner_blocs({{ "pression", &mat }});
  else Process::exit(que_suis_je() + " : dimensionner() not coded!");
}

DoubleTab& Operateur_Grad_base::ajouter(const DoubleTab& inco, DoubleTab& secmem) const
{
  if (has_interface_blocs())
    secmem *= -1, ajouter_blocs({}, secmem, {{ "pression", inco }}), secmem *= -1; /* to get the correct sign */
  else Process::exit(que_suis_je() + " : ajouter() not coded!");
  return secmem;
}

void Operateur_Grad_base::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  if (has_interface_blocs())
    {
      DoubleTrav secmem(equation().inconnue().valeurs()); //will be discarded
      matrice.get_set_coeff() *= -1, ajouter_blocs({{ "pression", &matrice }}, secmem), matrice.get_set_coeff() *= -1; /* to get the correct sign */
    }
  else Process::exit(que_suis_je() + " : contribuer_a_avec() not coded!");
}

/*! @brief Computation without boundary conditions ?
 *
 */
DoubleVect& Operateur_Grad_base::multvect(const DoubleTab& x, DoubleTab& b) const
{
  //
  //Cerr << que_suis_je()<<"::multvect is not overridden!!" << finl;
  // copy of X with zeros in the vector
  //  (we recover a result taking into account only the boundary conditions)
#if 0
  // Attempt at more efficient coding: does not work!
  // discrepancy on SGE_3D_VEF
  DoubleTrav x0(x); // copy of the structure, initialized to zero
  calculer(x0, b);
  operator_negate(b);
  ajouter(x, b);
#endif
  DoubleTrav x0(x);
  DoubleTrav b0(b);

  calculer(x0, b0);
  calculer(x,b);
  b-=b0;

  // We have subtracted from the result the contribution of the boundary conditions.
  return b;
}
