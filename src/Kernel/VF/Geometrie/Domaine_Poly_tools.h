/****************************************************************************
* Copyright (c) 2025, CEA
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

#ifndef Domaine_Poly_tools_included
#define Domaine_Poly_tools_included

#include <TRUSTTrav.h>
#include <Lapack.h>
#include <math.h>
#include <vector>

/* matrix product and transpose of DoubleTab */
static inline DoubleTab prod(DoubleTab a, DoubleTab b)
{
  int i, j, k, m = a.dimension(0), n = a.dimension(1), p = b.dimension(1);
  assert(n == b.dimension(0));
  DoubleTab r(m, p);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < p; k++) r(i, k) += a(i, j) * b(j, k);
  return r;
}
static inline DoubleTab transp(DoubleTab a)
{
  int i, j, m = a.dimension(0), n = a.dimension(1);
  DoubleTab r(n, m);
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++) r(j, i) = a(i, j);
  return r;
}

/*! @brief Solves the least squares problem ||M.x - b||_2, placing the kernel of M into P and returning the residual.
 *
 * This function uses Singular Value Decomposition (SVD) to solve the least squares problem for the equation M.x = b.
 * It decomposes the matrix M into its constituent parts using the LAPACK routine DGESVD, which is a Fortran subroutine
 * for computing the SVD of a general rectangular matrix. The kernel of M (null space) is placed into matrix P.
 *
 * @param M The input matrix of the linear system.
 * @param b The right-hand side vector of the linear system.
 * @param eps The threshold for considering singular values as zero.
 * @param P Pointer to a matrix where the kernel of M will be stored.
 * @param x The solution vector that minimizes the least squares problem.
 * @param S Vector to store the singular values of M.
 *
 * @return The residual of the least squares solution, i.e., the Euclidean norm of (M.x - b).
 */
static inline double kersol(const DoubleTab& M, DoubleTab& b, double eps, DoubleTab *P, DoubleTab& x, DoubleTab& S)
{
  int i, nk, m = M.dimension(0), n = M.dimension(1), k = std::min(m, n), l = std::max(m, n), w = 5 * l, info=-1, iP, jP;
  double res2 = 0;
  char a = 'A';
  //lapack en mode Fortran -> on decompose en fait Mt!!
  DoubleTab A = M, U(m, m), Vt(n, n), W(w), iS(n, m);
  S.resize(k);
  F77NAME(dgesvd)(&a, &a, &n, &m, A.addr(), &n, S.addr(), Vt.addr(), &n, U.addr(), &m, W.addr(), &w, &info);
  for (i = 0, nk = n; i < k && S(i) > eps * S(0); i++) nk--;
  if (P) P->resize(n, nk);
  for (i = 0, jP = -1; i < n; i++)
    if (i < k && S(i) > eps * S(0)) iS(i, i) = 1 / S(i); //diagonal term of iS
    else if (P)
      for (iP = 0, jP++; iP < n; iP++) (*P)(iP, jP) = Vt(i, iP); //column of V -> column of P
  x = prod(transp(Vt), prod(iS, prod(transp(U), b)));
  DoubleTab res = prod(M, x);
  for (i = 0; i < m; i++) res2 += std::pow(res(i, 0) - b(i, 0), 2);
  return sqrt(res2);
}

/*! @def CRIMP(a)
 * @brief Compacts a multi-dimensional array by resizing it.
 *
 * This macro adjusts the size of a multi-dimensional array `a` by incrementing and then decrementing its first dimension.
 * It handles arrays with 1, 2, or 3 dimensions.
 *
 * @param a The array to be compacted.
 */
#define CRIMP(a) a.nb_dim() > 2 ? a.resize(a.dimension(0) + 1, a.dimension(1), a.dimension(2)) : a.nb_dim() > 1 ? a.resize(a.dimension(0) + 1, a.dimension(1)) : a.resize(a.dimension(0) + 1), \
        a.nb_dim() > 2 ? a.resize(a.dimension(0) - 1, a.dimension(1), a.dimension(2)) : a.nb_dim() > 1 ? a.resize(a.dimension(0) - 1, a.dimension(1)) : a.resize(a.dimension(0) - 1)
#endif
