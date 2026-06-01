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

#include <Solv_TDMA.h>
#include <TRUSTVect.h>

void Solv_TDMA::resoudre(const DoubleVect& ma, const DoubleVect& mb, const DoubleVect& mc, const DoubleVect& sm, DoubleVect& vi, int M)
{
  DoubleVect malpha, mbeta; //2 intermediate vectors
  DoubleVect my; //solution vector of the system L.y=F
  int i;
  malpha.resize(M);
  mbeta.resize(M);
  my.resize(M);

  //The Thomas algorithm (TDMA) relies on the LU decomposition of the tridiagonal
  //matrix to be inverted. Let A.x = y be the linear system to solve with A = LU.
  //We first solve by forward substitution the system L.z = y with z = U.x, then
  //solve by back substitution the system U.x = z.

  // The 3 diagonals of the matrix are each stored as a vector:
  // ma: main diagonal; mb: lower sub-diagonal
  // mc: upper super-diagonal
  // M is the order of matrix A

  //Matrix L has a main diagonal (malpha) and a lower sub-diagonal (mb).
  //Matrix U has a main diagonal (all elements equal to 1) and an upper
  //super-diagonal (mbeta).

  //Verified solution of an antisymmetric 3x3 linear system (23.05.2003)

  //Fill malpha and mbeta

  malpha(0) = ma(0);
  mbeta(0) = mc(0)/ma(0);

  for (i = 1 ; i<M-1 ; i++)
    {
      malpha(i) = ma(i) - mb(i-1)*mbeta(i-1);
      mbeta(i) = mc(i)/malpha(i);
    }
  malpha(M-1) = ma(M-1) - mb(M-2)*mbeta(M-2);

  //Forward substitution to solve the first system

  my(0) = sm(0)/malpha(0);

  for (i = 1 ; i<M ; i++)
    my(i) = (sm(i)-mb(i-1)*my(i-1))/malpha(i);

  //Back substitution to solve the second system

  vi(M-1) = my(M-1);

  for (i = M-2 ; i>=0 ; i--)
    vi(i) = my(i) - mbeta(i)*vi(i+1);
}
