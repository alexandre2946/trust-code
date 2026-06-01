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

#ifndef Polynome_included
#define Polynome_included

#include <TRUSTTab.h>


/*! @brief Polynomial in n variables, n <= 4. Coefficients are stored using a DoubleTab.
 *
 */
class Polynome : public Objet_U
{
  Declare_instanciable(Polynome);

public :
  inline Polynome(int n1);
  inline Polynome(int n1, int n2);
  inline Polynome(int n1, int n2, int n3);
  inline Polynome(int n1, int n2, int n3, int n4);
  inline Polynome(const DoubleTab& );
  inline int degre(int ) const;
  inline double& coeff(int n1);
  inline double& coeff(int n1, int n2);
  inline double& coeff(int n1, int n2, int n3);
  inline double& coeff(int n1, int n2, int n3, int n4);
  inline double coeff(int n1) const;
  inline double coeff(int n1, int n2) const;
  inline double coeff(int n1, int n2, int n3) const;
  inline double coeff(int n1, int n2, int n3, int n4) const;
  double operator()(double x1) const;
  double operator()(double x1, double x2) const;
  double operator()(double x1, double x2, double x3) const;
  double operator()(double x1, double x2, double x3, double x4) const;

  void derive(int =0);
  void integre(int =0);

  double derive(double x1) const;
  double derive(double x1, double x2) const;
  double derive(double x1, double x2, double x3) const;
  double derive(double x1, double x2, double x3, double x4) const;

  double integre(double x1) const;
  double integre(double x1, double x2) const;
  double integre(double x1, double x2, double x3) const;
  double integre(double x1, double x2, double x3, double x4) const;

  Polynome& operator +=(const Polynome&);
  Polynome& operator -=(const Polynome&);
  Polynome& operator *=(const Polynome&);
  Polynome& operator *=(double);
  Polynome& operator /=(double);

private :
  DoubleTab coeff_;
};



/*! @brief Constructs a polynomial in one variable of degree n1.
 *
 * @param (int n1) degree of the polynomial
 */
inline Polynome::Polynome(int n1) : coeff_(++n1) {}

/*! @brief Constructs a polynomial in 2 variables of degrees n1 and n2.
 *
 * @param (int n1) degree of the first variable
 * @param (int n2) degree of the second variable
 */
inline Polynome::Polynome(int n1, int n2) : coeff_(++n1, ++n2) {}

/*! @brief Constructs a polynomial in 3 variables of degrees n1, n2 and n3.
 *
 * @param (int n1) degree of the first variable
 * @param (int n2) degree of the second variable
 * @param (int n3) degree of the third variable
 */
inline Polynome::Polynome(int n1, int n2, int n3) : coeff_(++n1, ++n2, ++n3) {}

/*! @brief Constructs a polynomial in 4 variables of degrees n1, n2, n3 and n4.
 *
 * @param (int n1) degree of the first variable
 * @param (int n2) degree of the second variable
 * @param (int n3) degree of the third variable
 * @param (int n4) degree of the fourth variable
 */
inline Polynome::Polynome(int n1, int n2, int n3, int n4) : coeff_(++n1, ++n2, ++n3, ++n4) {}

/*! @brief Constructs a polynomial from its coefficient array.
 *
 * @param (const DoubleTab& t) the coefficient array of the polynomial (1 to 4 dimensions)
 */
inline Polynome::Polynome(const DoubleTab& t) : coeff_(t) {}

/*! @brief Returns the degree of the polynomial with respect to the i-th variable.
 *
 * @param (int i) index of the variable
 * @return (int) degree of the polynomial
 */
inline int Polynome::degre(int i) const
{
  return coeff_.dimension(i)-1;
}

/*! @brief Returns the coefficient of the term of degree n1.
 *
 * @param (int n1) degree
 * @return (double&) coefficient
 */
inline double& Polynome::coeff(int n1)
{
  return coeff_(n1);
}

/*! @brief Returns the coefficient of the term of degree (n1, n2).
 *
 * @param (int n1) degree in first variable
 * @param (int n2) degree in second variable
 * @return (double&) coefficient
 */
inline double& Polynome::coeff(int n1, int n2)
{
  return coeff_(n1, n2);
}

/*! @brief Returns the coefficient of the term of degree (n1, n2, n3).
 *
 * @param (int n1) degree in first variable
 * @param (int n2) degree in second variable
 * @param (int n3) degree in third variable
 * @return (double&) coefficient
 */
inline double& Polynome::coeff(int n1, int n2, int n3)
{
  return coeff_(n1, n2, n3);
}

/*! @brief Returns the coefficient of the term of degree (n1, n2, n3, n4).
 *
 * @param (int n1) degree in first variable
 * @param (int n2) degree in second variable
 * @param (int n3) degree in third variable
 * @param (int n4) degree in fourth variable
 * @return (double&) coefficient
 */
inline double& Polynome::coeff(int n1, int n2, int n3, int n4)
{
  return coeff_(n1, n2, n3, n4);
}

/*! @brief Returns the coefficient of the term of degree n1 (const).
 *
 * @param (int n1) degree
 * @return (double) coefficient
 */
inline double Polynome::coeff(int n1) const
{
  return coeff_(n1);
}

/*! @brief Returns the coefficient of the term of degree (n1, n2) (const).
 *
 * @param (int n1) degree in first variable
 * @param (int n2) degree in second variable
 * @return (double) coefficient
 */
inline double Polynome::coeff(int n1, int n2) const
{
  return coeff_(n1, n2);
}

/*! @brief Returns the coefficient of the term of degree (n1, n2, n3) (const).
 *
 * @param (int n1) degree in first variable
 * @param (int n2) degree in second variable
 * @param (int n3) degree in third variable
 * @return (double) coefficient
 */
inline double Polynome::coeff(int n1, int n2, int n3) const
{
  return coeff_(n1, n2, n3);
}

/*! @brief Returns the coefficient of the term of degree (n1, n2, n3, n4) (const).
 *
 * @param (int n1) degree in first variable
 * @param (int n2) degree in second variable
 * @param (int n3) degree in third variable
 * @param (int n4) degree in fourth variable
 * @return (double) coefficient
 */
inline double Polynome::coeff(int n1, int n2, int n3, int n4) const
{
  return coeff_(n1, n2, n3, n4);
}
#endif        //Polynome_H

