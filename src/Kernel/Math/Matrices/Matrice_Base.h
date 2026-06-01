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

#ifndef Matrice_Base_included
#define Matrice_Base_included

#include <TRUSTTab.h>
#include <Matrix_tools.h>

/*! @brief Matrice_Base class - Base class of the matrix hierarchy.
 *
 *     This class defines the generic interface of a matrix as used in Trio-U.
 *     Consequently it is not instantiable. All matrix types must derive from
 *     this base class and implement its abstract methods.
 *
 *     In the method comments, A represents a Matrice_base object.
 *
 * @sa Abstract class, Abstract methods:, DoubleVect& multvect_(const DoubleVect&, DoubleVect& ) const, int ordre() const
 */

class Matrice_Base : public Objet_U
{

  Declare_base(Matrice_Base);

public :
  /// If square matrix, returns number of lines, otherwise 0
  virtual int ordre() const=0;
  /// Return local number of lines (=size on the current proc)
  virtual int nb_lignes() const=0;
  /// Return local number of columns (=size on the current proc)
  virtual int nb_colonnes() const=0;

  // Methods for computing r+=Ax, implemented in derived classes
  virtual DoubleVect& ajouter_multvect_(const DoubleVect& x, DoubleVect& r) const =0;
  virtual DoubleVect& ajouter_multvectT_(const DoubleVect& x, DoubleVect& r) const =0;
  virtual DoubleTab& ajouter_multTab_(const DoubleTab& x, DoubleTab& r) const =0;
  // Methods for computing r+=Ax with echange_espace_virtuel(), implemented in Matrice_Base
  virtual inline DoubleVect& ajouter_multvect(const DoubleVect& x, DoubleVect& r) const;
  virtual inline DoubleVect& ajouter_multvectT(const DoubleVect& x, DoubleVect& r) const;
  virtual inline DoubleTab& ajouter_multTab(const DoubleTab& , DoubleTab& r) const;
  // Methods for computing r=Ax, implemented in Matrice_Base
  virtual inline DoubleVect& multvect_(const DoubleVect&, DoubleVect& ) const;
  virtual inline DoubleVect& multvect(const DoubleVect&, DoubleVect& ) const;
  virtual inline DoubleVect& multvectT_(const DoubleVect&, DoubleVect& ) const;
  virtual inline DoubleVect& multvectT(const DoubleVect&, DoubleVect& ) const;
  virtual inline DoubleTab& multTab(const DoubleTab& , DoubleTab& r) const;
  // Matrix printing methods
  virtual inline Sortie& imprimer(Sortie&) const;
  virtual inline Sortie& imprimer_formatte(Sortie& s) const;

  friend DoubleVect operator*(const Matrice_Base&, const DoubleVect&);

  virtual void scale(const double x) =0;

  // Zero out the matrix values
  virtual void clean() { Process::exit("Matrice_base::clean() not implemented.");};

  virtual void get_stencil(Stencil& stencil) const;

  virtual void get_symmetric_stencil(Stencil& stencil) const;

  virtual void get_stencil_and_coefficients(Stencil& stencil, StencilCoeffs& coefficients) const;
  virtual void get_stencil_and_coeff_ptrs(Stencil& stencil, std::vector<const double *>& coeff_ptr) const;

  virtual void get_symmetric_stencil_and_coefficients(Stencil& stencil, StencilCoeffs& coefficients) const;

  int get_stencil_size() const ;
  virtual void build_stencil();

  void set_stencil( const Stencil& stencil );

  bool is_stencil_up_to_date() const ;

protected:
  bool is_stencil_up_to_date_ = false;
  Stencil stencil_ ;
};



/*! @brief Multiplication of a vector by the matrix.
 *
 * Operation: r = A*x
 *
 * @param (DoubleVect& x) the vector to multiply
 * @param (DoubleVect& r) the result vector of the operation
 * @return (DoubleVect&) the result vector of the operation
 */
inline DoubleVect& Matrice_Base::
multvect(const DoubleVect& x, DoubleVect& r) const
{
  multvect_(x,r);
  r.echange_espace_virtuel();
  return r;
}

inline DoubleVect& Matrice_Base::
multvect_(const DoubleVect& x, DoubleVect& r) const
{
  r=0;
  ajouter_multvect_(x,r);
  return r;
}

/*! @brief Multiplication of a vector by the transposed matrix.
 *
 * Operation: r = AT*x
 *
 * @param (DoubleVect& x) the vector to multiply
 * @param (DoubleVect& r) the result vector of the operation
 * @return (DoubleVect&) the result vector of the operation
 */
inline DoubleVect& Matrice_Base::
multvectT(const DoubleVect& x, DoubleVect& r) const
{
  multvectT_(x,r);
  r.echange_espace_virtuel();
  return r;
}

inline DoubleVect& Matrice_Base::
multvectT_(const DoubleVect& x, DoubleVect& r) const
{
  r=0;
  ajouter_multvectT_(x,r);
  return r;
}

/*! @brief NOT IMPLEMENTED Multiplication of a matrix represented by an array by the matrix.
 *
 *     Operation: R = A*X
 *
 * @param (DoubleTab&) the matrix to multiply
 * @param (DoubleTab& r) the result matrix of the operation
 * @return (DoubleTab&) the result matrix of the operation
 * @throws NOT IMPLEMENTED
 */
inline DoubleTab& Matrice_Base::
multTab(const DoubleTab& x, DoubleTab& r) const
{
  r=0;
  ajouter_multTab_(x,r);
  r.echange_espace_virtuel();
  return r;
}


/*! @brief Matrix-vector multiply-accumulate operation (saxpy).
 *
 * Operation: r = r + A*x
 *
 * @param (DoubleVect& x) the vector to multiply
 * @param (DoubleVect& r) the result vector of the operation
 * @return (DoubleVect&) the result vector of the operation
 */
inline DoubleVect& Matrice_Base::
ajouter_multvect(const DoubleVect& x, DoubleVect& r) const
{
  ajouter_multvect_(x,r);
  r.echange_espace_virtuel();
  return r;
}

/*! @brief Matrix-vector multiply-accumulate operation (saxpy).
 *
 * Operation: r = r + A*x
 *
 * @param (DoubleVect& x) the vector to multiply
 * @param (DoubleVect& r) the result vector of the operation
 * @return (DoubleVect&) the result vector of the operation
 */
inline DoubleVect& Matrice_Base::
ajouter_multvectT(const DoubleVect& x, DoubleVect& r) const
{
  ajouter_multvectT_(x,r);
  r.echange_espace_virtuel();
  return r;
}

/*! @brief NOT IMPLEMENTED Matrix-matrix multiply-accumulate operation (saxpy)
 *
 *     (matrix represented by an array)
 *     Operation: R = R + A*X
 *
 * @param (DoubleTab&) the matrix to multiply
 * @param (DoubleTab& r) the result matrix of the operation
 * @return (DoubleTab&) the result matrix of the operation
 * @throws NOT IMPLEMENTED
 */
inline DoubleTab& Matrice_Base::
ajouter_multTab(const DoubleTab& x, DoubleTab& r) const
{
  ajouter_multTab_(x,r);
  r.echange_espace_virtuel();
  return r;
}


/*! @brief Friend function (outside the class) of the Matrice_Base class.
 *
 * Multiplication operator: returns (A*vect)
 *
 * @param (Matrice_Base& A) the multiplying matrix
 * @param (DoubleVect& vect) the vector to multiply
 * @return (DoubleVect) the result vector of the operation
 */
DoubleVect operator * (const Matrice_Base& A, const DoubleVect& vect);


inline Sortie& Matrice_Base::imprimer(Sortie& os) const
{
  return os<<*this;
}

inline Sortie& Matrice_Base::imprimer_formatte(Sortie& os) const
{
  return os<<*this;
}

#endif
