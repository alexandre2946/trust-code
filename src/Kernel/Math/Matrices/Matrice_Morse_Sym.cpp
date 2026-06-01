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

#include <Matrice_Morse_Sym.h>
#include <TRUSTArrays.h>
#include <Array_tools.h>
#include <TRUSTTabs.h>

#include <Sparskit.h>
#include <Noms.h>

Implemente_instanciable_sans_constructeur(Matrice_Morse_Sym,"Matrice_Morse_Sym",Matrice_Morse);

/*! @brief Writes the three arrays of the Morse storage structure to an output stream.
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie& s) the modified output stream
 */
Sortie& Matrice_Morse_Sym::printOn(Sortie& s) const
{
  s<<get_est_definie()<<finl;
  return  Matrice_Morse::printOn(s);
}


/*! @brief NOT IMPLEMENTED
 *
 * @param (Entree& s) an input stream
 * @return (Entree& s) the input stream
 * @throws NOT IMPLEMENTED
 */
Entree& Matrice_Morse_Sym::readOn(Entree& s)
{
  int est_def;
  s>>est_def;
  set_est_definie(est_def);
  return Matrice_Morse::readOn(s) ;
}

/*! @brief Constructor of a Matrice_Morse_Sym with n1 rows and capable of storing at most n2 non-zero elements.
 *
 *     The matrix elements and the connectivity table
 *     are given in the last 3 parameters.
 *
 * @param (int n1) the number of rows of the matrix
 * @param (int n2) the maximum number of non-zero elements storable by the matrix
 * @param (IntLists& voisins) list of neighbours for each row
 * @param (DoubleLists& valeurs) list of values
 * @param (DoubleVect& terme_diag) the vector of diagonal terms
 */
Matrice_Morse_Sym::Matrice_Morse_Sym(int n1, int n2, const IntLists& voisins,
                                     const DoubleLists& valeurs,
                                     const DoubleVect& terme_diag)
  : Matrice_Morse(n1, n2)
{
  symetrique_=1;
  remplir(voisins, valeurs, terme_diag);
  morse_matrix_structure_has_changed_=1;
}


/*! @brief Copy constructor of a Matrice_Morse_Sym from a Matrice_Morse.
 *
 * @param (Matrice_Morse& acopier) the matrix to copy
 */
Matrice_Morse_Sym::Matrice_Morse_Sym(const Matrice_Morse& A)
  :Matrice_Morse(A)
{
  symetrique_=1;
  partie_sup(A);
  compacte();
  morse_matrix_structure_has_changed_=1;
}
Matrice_Morse_Sym::Matrice_Morse_Sym(const Matrice& A)
{
  const Matrice_Base& A_base=A.valeur();
  if(sub_type(Matrice_Morse_Sym, A_base))
    *this=ref_cast(Matrice_Morse_Sym, A_base);
  else if(sub_type(Matrice_Morse, A_base))
    *this=ref_cast(Matrice_Morse, A_base);
  symetrique_=1;
  morse_matrix_structure_has_changed_=1;
}
Matrice_Morse_Sym& Matrice_Morse_Sym::operator=(const Matrice_Morse& A)
{
  tab1_=A.get_tab1();
  tab2_=A.get_tab2();
  coeff_=A.get_coeff();
  m_=A.nb_colonnes();
  symetrique_=A.get_symmetric();
  partie_sup(A);
  compacte();
  morse_matrix_structure_has_changed_=1;
  is_stencil_up_to_date_=A.is_stencil_up_to_date();
  return *this;
}
Matrice_Morse_Sym& Matrice_Morse_Sym::operator=(const Matrice& A)
{
  const Matrice_Base& A_base=A.valeur();
  if(sub_type(Matrice_Morse_Sym, A_base))
    *this=ref_cast(Matrice_Morse_Sym, A_base);
  else if(sub_type(Matrice_Morse, A_base))
    *this=ref_cast(Matrice_Morse, A_base);
  morse_matrix_structure_has_changed_=1;
  return *this;
}

/*! @brief Copy constructor of a Matrice_Morse_Sym from a Matrice_Morse.
 *
 * @param (Matrice_Morse& acopier) the matrix to copy
 */
Matrice_Morse_Sym::Matrice_Morse_Sym(const Matrice_Morse_Sym& acopier) :
  Matrice_Morse(acopier),Matrice_Sym()
{
  set_est_definie(acopier.get_est_definie());
  symetrique_=1;
  morse_matrix_structure_has_changed_=1;
}

/*! @brief Matrix-matrix multiply-accumulate operation (saxpy) where the matrix is represented by an array.
 *
 *     Operation: RESU = RESU + A*X
 *
 * @param (DoubleTab& x) the matrix to multiply
 * @param (DoubleTab& resu) the result matrix of the operation
 * @return (DoubleTab&) the result matrix of the operation
 * @throws result size incompatible with the size of x
 */
DoubleTab& Matrice_Morse_Sym::ajouter_multTab_(const DoubleTab& x, DoubleTab& resu) const
{
  if ( (x.nb_dim() == 1) && (resu.nb_dim() == 1))
    {
      ajouter_multvect(x,resu);
      return resu;
    }

  assert_check_symmetric_morse_matrix_structure( );
  int i,j;
  double aij;
  int comp;
  int nb_com=x.dimension(1);
  assert(resu.dimension(1)==nb_com);
  double* t=new double[nb_com];
  int n=nb_lignes();
  for(i=0; i<n; i++)
    {
      for(comp=0; comp<nb_com; comp++)
        t[comp] = coeff_(tab1_(i)-1)*x(i,comp);
      for (auto k=tab1_(i); k<tab1_(i+1)-1; k++)
        {
          j = tab2_(k)-1;
          aij = coeff_(k);
          for(comp=0; comp<nb_com; comp++)
            {
              t[comp] += aij*x(j,comp);
              resu(j,comp) += aij*x(i,comp);
            }
        }
      for(comp=0; comp<nb_com; comp++)
        resu(i,comp) += t[comp] ;
    }
  delete [] t;
  return resu;
}

double Matrice_Morse_Sym::multvect_et_prodscal(const DoubleVect& x, DoubleVect& resu) const
{
  assert_check_symmetric_morse_matrix_structure( );

  assert(x.size_totale() == nb_lignes() || x.size() == nb_lignes());
  assert(resu.size_totale() == nb_lignes() || resu.size() == nb_lignes());
  assert(m_ <= x.size_totale());

  double prod_scal_local = 0.;
  operator_egal(resu, 0.);

  const int fin = nb_lignes() + 1;
  {
    const double * thecoef = get_coeff().addr() - 1;
    const double * xx = x.addr() - 1;
    double * res = resu.addr() - 1;
    const auto * index1 = get_tab1().addr() - 1;
    const int * index2 = get_tab2().addr() - 1;
#ifdef COMPILER_BLRTS_XLC
#pragma disjoint (*coef, *xx, *res)
#endif
    int i;
    for (i = 1; i < fin; i++)
      {
        auto j = index1[i];
        auto j_next= index1[i+1];
        //int ncoeffs = j_next - j;
        assert(i==index2[j]); // Even a zero diagonal element must be stored in a Mat_Morse_Sym
        double xi = xx[i];
        double resu_tmp = res[i] + thecoef[j] * xi;
        j++;
        for (; j < j_next; j++)
          {
            int nj     = index2[j];
            double coef_j = thecoef[j];
            double xn     = xx[nj];
            resu_tmp   += coef_j * xn;
            res[nj] += coef_j * xi;
          }
        res[i] = resu_tmp;
        prod_scal_local += resu_tmp * xi;
      }
  }

  return prod_scal_local;
}


/*! @brief Matrix-vector multiply-accumulate operation (saxpy).
 *
 * Operation: resu = resu + A*x
 *
 * @param (DoubleVect& x) the vector to multiply
 * @param (DoubleVect& resu) the result vector of the operation
 * @return (DoubleVect&) the result vector of the operation
 */

DoubleVect& Matrice_Morse_Sym::ajouter_multvect_(const DoubleVect& x, DoubleVect& resu) const
{
  assert_check_symmetric_morse_matrix_structure( );

  assert(x.size_totale() == nb_lignes() || x.size() == nb_lignes());
  assert(resu.size_totale() == nb_lignes() || resu.size() == nb_lignes());
  assert(m_ <= x.size_totale());

  const int fin = nb_lignes() + 1;
  {
    const double * thecoef = get_coeff().addr() - 1;
    const double * xx = x.addr() - 1;
    double * res = resu.addr() - 1;
    const auto * index1 = get_tab1().addr() - 1;
    const int * index2 = get_tab2().addr() - 1;
#ifdef COMPILER_BLRTS_XLC
#pragma disjoint (*coef, *xx, *res)
#endif
    int i;
    for (i = 1; i < fin; i++)
      {
        auto j = index1[i];
        auto j_next= index1[i+1];
        assert(i==index2[j]); // Even a zero diagonal element must be stored in a Mat_Morse_Sym
#if 0
        // Note B.M.: sur les procs intel ce deroulage n'a aucun effet benefique.
        int ncoeffs = j_next - j;
        if (ncoeffs == 4)
          {
            int n2 = index2[j+1];
            int n3 = index2[j+2];
            int n4 = index2[j+3];
            double xi = xx[i];
            double res_i  = res[i];
            double res_n2 = res[n2];
            double res_n3 = res[n3];
            double res_n4 = res[n4];
            double a1 = thecoef[j];
            double a2 = thecoef[j+1];
            double a3 = thecoef[j+2];
            double a4 = thecoef[j+3];
            res[i] = res_i + a1 * xi + a2 * xx[n2] + a3 * xx[n3] + a4 * xx[n4];
            res[n2] = res_n2 + a2 * xi;
            res[n3] = res_n3 + a3 * xi;
            res[n4] = res_n4 + a4 * xi;
          }
        else
#endif
          {
            double xi = xx[i];
            double resu_tmp = res[i] + thecoef[j] * xi;
            j++;
            for (; j < j_next; j++)
              {
                int nj     = index2[j];
                double coef_j = thecoef[j];
                double xn     = xx[nj];
                resu_tmp   += coef_j * xn;
                res[nj] += coef_j * xi;
              }
            res[i] = resu_tmp;
          }
      }
  }
#if 0
  int fin2,l;
  m=tab1_(0);
  for(i=debut; i<n; i++)
    {
      xi=x(i);
      l=m;
      m=tab1_(i+1);
      t = coeff_(l-1)*xi;
      assert(i==tab2(l-1)-1); // Even a zero diagonal element must be stored in a Mat_Morse_Sym
      //aij=coeff_(tab1_(i)-1);
      fin2=m-1;
      for (k=l; k<fin2; k++)
        {
          j = tab2_(k)-1;
          assert(j<n);
          aij = coeff_(k);
          t += aij*x(j);
          resu(j) += aij*xi;
        }
      resu(i) += t ;
    }
#endif
  return resu;
}


/*! @brief Transposed-matrix-vector multiply-accumulate operation (saxpy).
 *
 *     Operation: resu = resu + A^{T}*x
 *
 * @param (DoubleVect& x) the vector to multiply
 * @param (DoubleVect& resu) the result vector of the operation
 * @return (DoubleVect&) the result vector of the operation
 */
DoubleVect& Matrice_Morse_Sym::ajouter_multvectT_(const DoubleVect& x,DoubleVect& resu) const
{
  Cerr <<"Matrice_Morse_Sym::ajouter_multvectT_ is not coded" << finl;
  exit();
  return resu;
}


/*! @brief Friend function (non-member) of class Matrice_Morse_Sym. DOES NOTHING: NOT IMPLEMENTED
 *
 * @param (Matrice_Morse_Sym&)
 * @param (Matrice_Morse_Sym&)
 * @return (Matrice_Morse_Sym)
 */
Matrice_Morse_Sym operator+(const Matrice_Morse_Sym& A,const Matrice_Morse_Sym& B)
{
  //The two matrices must obviously have the same size
  if (A.nb_lignes() != B.nb_lignes())
    {
      Cerr << "Error in Matrice_Morse_Sym::operator+" << finl;
      Cerr << "The 2 matrices must have the same dimension" << finl;
      Cerr << "The first matrix is of dimension: " << A.nb_lignes()<< finl;
      Cerr << "The second matrix is of dimension: " << B.nb_lignes()<< finl;
      Process::exit();
    }

  Matrice_Morse_Sym somme(A.nb_lignes(),(int)(A.nb_coeff()+B.nb_coeff()));
  IntList colonnes_A_ligne_i,colonnes_B_ligne_i;
  int min_colonnes_A,min_colonnes_B,non_nuls_ligne_i;
  DoubleList coeffs_A_ligne_i,coeffs_B_ligne_i;
  double coeff_A,coeff_B;

  //To ensure that compacte() will work
  for (auto i=0; i<somme.get_coeff().size(); i++)
    somme.get_set_coeff()(i) = 0.;

  //Fill tab2_ in ascending column order
  somme.get_set_tab1()(0) = 1;
  for (int i=0; i<A.nb_lignes(); i++) //loop over rows
    {
      //To avoid boundary effects
      if (!colonnes_A_ligne_i.est_vide()) colonnes_A_ligne_i.vide();
      if (!colonnes_B_ligne_i.est_vide()) colonnes_B_ligne_i.vide();

      //Initialisation of colonnes_A and colonnes_B
      //BEWARE of C++ indexing
      for (auto j=0; j<A.get_tab1()(i+1)-A.get_tab1()(i); j++)
        {
          colonnes_A_ligne_i.add_if_not(A.get_tab2()(A.get_tab1()(i)+j-1));
          coeffs_A_ligne_i.add(A.get_coeff()(A.get_tab1()(i)+j-1));
        }//fin for

      for (auto j=0; j<B.get_tab1()(i+1)-B.get_tab1()(i); j++)
        {
          colonnes_B_ligne_i.add_if_not(B.get_tab2()(B.get_tab1()(i)+j-1));
          coeffs_B_ligne_i.add(B.get_coeff()(B.get_tab1()(i)+j-1));
        }//fin for

      //Initialize other variables
      non_nuls_ligne_i = 0;

      //Start of the algorithm
      while (!colonnes_A_ligne_i.est_vide() || !colonnes_B_ligne_i.est_vide())
        {
          //If one list is empty but not the other
          if (colonnes_A_ligne_i.est_vide() && !colonnes_B_ligne_i.est_vide())
            {
              somme.get_set_tab2()(somme.get_tab1()(i)+non_nuls_ligne_i-1) =
                colonnes_B_ligne_i[0];

              somme.get_set_coeff()(somme.get_tab1()(i)+non_nuls_ligne_i-1) =
                coeffs_B_ligne_i[0] ;

              colonnes_B_ligne_i.suppr(colonnes_B_ligne_i[0]);
              coeffs_B_ligne_i.suppr(coeffs_B_ligne_i[0]);
              non_nuls_ligne_i++;
            }//fin if

          if (colonnes_B_ligne_i.est_vide() && !colonnes_A_ligne_i.est_vide())
            {
              somme.get_set_tab2()(somme.get_tab1()(i)+non_nuls_ligne_i-1) =
                colonnes_A_ligne_i[0];

              somme.get_set_coeff()(somme.get_tab1()(i)+non_nuls_ligne_i-1) =
                coeffs_A_ligne_i[0];

              colonnes_A_ligne_i.suppr(colonnes_A_ligne_i[0]);
              coeffs_A_ligne_i.suppr(coeffs_A_ligne_i[0]);
              non_nuls_ligne_i++;
            }//fin if

          //Neither list is empty
          min_colonnes_A = colonnes_A_ligne_i[0];
          min_colonnes_B = colonnes_B_ligne_i[0];
          coeff_A = coeffs_A_ligne_i[0];
          coeff_B = coeffs_B_ligne_i[0];

          if (min_colonnes_A < min_colonnes_B)
            {
              somme.get_set_tab2()(somme.get_tab1()(i)+non_nuls_ligne_i-1) = min_colonnes_A;
              somme.get_set_coeff()(somme.get_tab1()(i)+non_nuls_ligne_i-1) = coeff_A;

              colonnes_A_ligne_i.suppr(min_colonnes_A);
              coeffs_A_ligne_i.suppr(coeff_A);
              non_nuls_ligne_i++;
            }//fin if

          if (min_colonnes_B < min_colonnes_A)
            {
              somme.get_set_tab2()(somme.get_tab1()(i)+non_nuls_ligne_i-1) = min_colonnes_B;
              somme.get_set_coeff()(somme.get_tab1()(i)+non_nuls_ligne_i-1) = coeff_B;

              colonnes_B_ligne_i.suppr(min_colonnes_B);
              coeffs_B_ligne_i.suppr(coeff_B);
              non_nuls_ligne_i++;
            }//fin if

          if (min_colonnes_A == min_colonnes_B)
            {
              somme.get_set_tab2()(somme.get_tab1()(i)+non_nuls_ligne_i-1) = min_colonnes_A;
              somme.get_set_coeff()(somme.get_tab1()(i)+non_nuls_ligne_i-1) = coeff_A+coeff_B;

              colonnes_A_ligne_i.suppr(min_colonnes_A);
              colonnes_B_ligne_i.suppr(min_colonnes_B);
              coeffs_A_ligne_i.suppr(coeff_A);
              coeffs_B_ligne_i.suppr(coeff_B);
              non_nuls_ligne_i++;
            }//fin if

        }// fin while

      somme.get_set_tab1()(i+1) = somme.get_tab1()(i)+non_nuls_ligne_i;//?????
    }// fin for

  somme.compacte();
  somme.morse_matrix_structure_has_changed_=1;
  return somme;
}


/*! @brief Friend function (non-member) of class Matrice_Morse_Sym. DOES NOTHING: NOT IMPLEMENTED
 *
 * @param (double)
 * @param (Matrice_Morse_Sym& A)
 * @return (Matrice_Morse_Sym)
 * @throws NOT IMPLEMENTED
 */
Matrice_Morse_Sym operator *(double x, const Matrice_Morse_Sym& A)
{
  Matrice_Morse_Sym mat_res(A);
  mat_res.get_set_coeff()*=x;
  mat_res.morse_matrix_structure_has_changed_=1;
  return(mat_res);
}


/*! @brief Friend function (non-member) of class Matrice_Morse_Sym. Simply calls operator*(double,const Matrice_Morse_Sym&) (which is NOT IMPLEMENTED)
 *
 * @param (Matrice_Morse_Sym& A) the matrix to multiply by x
 * @param (double x) a scalar
 * @return (Matrice_Morse_Sym) the result of the underlying call
 */
Matrice_Morse_Sym operator *(const Matrice_Morse_Sym& A, double x)
{
  return(x*A);
}


/*! @brief DOES NOTHING: NOT IMPLEMENTED
 *
 * @param (double)
 * @return (Matrice_Morse_Sym)
 * @throws NOT IMPLEMENTED
 */
Matrice_Morse_Sym& Matrice_Morse_Sym::operator *=( double x )
{
  scale( x );
  return(*this);
}

void Matrice_Morse_Sym::scale( const double x )
{
  coeff_ *= x;
}

void Matrice_Morse_Sym::get_stencil( Stencil& stencil ) const
{
  assert_check_symmetric_morse_matrix_structure( );

  Stencil symmetric_stencil;
  get_symmetric_stencil( symmetric_stencil );

  Matrice_Sym::unsymmetrize_stencil( nb_lignes( ),
                                     symmetric_stencil,
                                     stencil );
}

void Matrice_Morse_Sym::get_symmetric_stencil( Stencil& stencil ) const
{
  assert_check_symmetric_morse_matrix_structure( );

  stencil.resize( 0, 2 );


  ArrOfInt tmp;


  const int nb_lines = nb_lignes( );
  for ( int i=0; i<nb_lines; ++i )
    {
      auto k0 = tab1_( i ) - 1;
      auto k1 = tab1_( i + 1 ) - 1;
      int size = (int)(k1 - k0);

      tmp.resize_array( 0 );
      tmp.resize_array( size );

      for ( int k=0; k<size; ++k )
        {
          tmp[ k ] = tab2_( k + k0 ) - 1;
        }

      tmp.ordonne_array( );

      for ( int k=0; k<size; ++k )
        {
          stencil.append_line( i, tmp[ k ] );
        }
    }

  const int new_size = stencil.dimension( 0 );

  stencil.resize( new_size, 2 );
}

void Matrice_Morse_Sym::get_stencil_and_coefficients( Stencil&      stencil,
                                                      StencilCoeffs& coefficients ) const
{
  assert_check_symmetric_morse_matrix_structure( );

  Stencil symmetric_stencil;
  StencilCoeffs symmetric_coefficients;
  get_symmetric_stencil_and_coefficients( symmetric_stencil, symmetric_coefficients );

  Matrice_Sym::unsymmetrize_stencil_and_coefficients( nb_lignes( ),
                                                      symmetric_stencil,
                                                      symmetric_coefficients,
                                                      stencil,
                                                      coefficients );
}


void Matrice_Morse_Sym::get_symmetric_stencil_and_coefficients( Stencil&      stencil,
                                                                StencilCoeffs& coefficients ) const
{
  assert_check_symmetric_morse_matrix_structure( );

  stencil.resize( 0, 2 );


  coefficients.resize( 0 );


  IntTab tmp1(0);


  ArrOfDouble tmp2;


  ArrOfInt index;


  const int nb_lines = nb_lignes( );
  for ( int i=0; i<nb_lines; ++i )
    {
      auto k0   = tab1_( i ) - 1;
      auto k1   = tab1_( i + 1 ) - 1;
      int size = (int)(k1 - k0);

      index.resize_array( 0 );

      tmp1.resize( 0 );
      tmp1.resize( size );

      tmp2.resize_array( 0 );
      tmp2.resize_array( size );

      for ( int k=0; k<size; ++k )
        {
          tmp1( k ) = tab2_( k + k0 ) - 1;
          tmp2[ k ] = coeff_( k + k0 );
        }
      tri_lexicographique_tableau_indirect( tmp1, index );

      for ( int k=0; k<size; ++k )
        {
          int l = index[ k ];
          stencil.append_line( i, tmp1[ l ] );
          coefficients.append_array( tmp2[ l ] );
        }
    }

  const int new_size = stencil.dimension( 0 );
  assert( coefficients.size_array( ) == new_size );


  stencil.resize( new_size, 2 );


  coefficients.resize_array( new_size );
}

/*! @brief Unary negation operator, returns the opposite of the matrix: -A.
 *
 *     Calls operator*(const Matrice_Morse_Sym&,double)
 *
 */
Matrice_Morse_Sym Matrice_Morse_Sym::operator -() const
{
  return((*this)*(-1));
}
/*
static int commun(const ArrOfInt& items,
                  int prems,
                  int der,
                  int val,
                  int& item_courant)
{
  if (val==items[item_courant])
    return item_courant;
  if((val<items[0]) || (val>items[der]))
    return item_courant;
  if (val<items[item_courant])
    {
      --item_courant;
      if(val>=items[item_courant]) return item_courant;
      der=item_courant;
      item_courant=(prems+item_courant)/2;
      return commun(items,prems,der,val,item_courant);
    }
  if (val>items[item_courant])
    {
      ++item_courant;
      if(val<items[item_courant]) return (--item_courant);
      if(val==items[item_courant]) return (item_courant);
      prems=item_courant;
      item_courant=(item_courant+der)/2;
      return commun(items,prems,der,val,item_courant);
    }
  return item_courant;
}
*/


/*! @brief Assignment operator of a Matrice_Morse_Sym into a Matrice_Morse_Sym.
 *
 * @param (Matrice_Morse_Sym& a) the right-hand side of the assignment
 * @return (Matrice_Morse_Sym&) the result of the assignment (*this)
 */
Matrice_Morse_Sym& Matrice_Morse_Sym::operator=(const Matrice_Morse_Sym& a )
{
  tab1_  =a.get_tab1();
  tab2_  =a.get_tab2();
  coeff_ =a.get_coeff();
  m_=a.nb_colonnes();
  set_est_definie(a.get_est_definie());
  morse_matrix_structure_has_changed_=1;
  is_stencil_up_to_date_ =  a.is_stencil_up_to_date() ;
  return(*this);
}

int Matrice_Morse_Sym_test()
{
  return 1;
}


int Matrice_Morse_Sym::inverse(const DoubleVect& secmem, DoubleVect& solution,
                               double coeff_seuil) const
{
  Cerr << "Not coded." << finl;
  exit();
  return 0;
}

int Matrice_Morse_Sym::inverse(const DoubleVect& secmem, DoubleVect& solution,
                               double coeff_seuil, int max_iter) const
{
  Cerr << "Not coded." << finl;
  exit();
  return 0;
}

Sortie& Matrice_Morse_Sym::imprimer_formatte(Sortie& s) const
{
  return Matrice_Morse::imprimer_formatte(s, 1);
}

/*! @brief Remove duplicates by sorting tab2.
 *
 *     Check the diagonal: all diagonal elements must be stored.
 *
 */
void Matrice_Morse_Sym::compacte(int elim_coeff_nul)
{
  Matrice_Morse::compacte(elim_coeff_nul);
  // Check whether all diagonal elements are properly stored
  int n=nb_lignes();
  int elements_diagonaux_non_stockes=0;
  auto size_tab2_ = tab2_.size_array();
  for (int i=0; i<n; i++)
    {
      auto k=tab1_(i)-1;
      if (k>=size_tab2_ || i!=tab2_(k)-1)
        {
          elements_diagonaux_non_stockes++;
          // Cerr << "Probleme rencontre dans une Matrice_Morse_Sym:" << finl;
          // Cerr << "La ligne "<<i<<" n'a pas sa diagonale correctement stockee." << finl;
          // Cerr << "Les elements diagonaux meme nuls doivent etre stockes dans une Mat_Morse_Sym TRUST." << finl;
        }
    }
  if (elements_diagonaux_non_stockes)
    {
      // Resize
      auto nnz=size_tab2_+elements_diagonaux_non_stockes;
      tab2_.resize(nnz);
      coeff_.resize(nnz);
      // Modify the matrix starting from the end
      for (int i=n-1; i>=0; i--)
        {
          // Copy with shift
          auto k1=tab1_(i)-1;
          auto k2=tab1_(i+1)-1;
          for (auto j=k2-1; j>=k1; j--)
            {
              tab2_(j+elements_diagonaux_non_stockes)=tab2_(j);
              coeff_(j+elements_diagonaux_non_stockes)=coeff_(j);
            }
          tab1_(i+1)+=elements_diagonaux_non_stockes;

          if (k1>=size_tab2_ || i!=tab2_(k1)-1)
            {
              // Insert the missing zero diagonal element
              elements_diagonaux_non_stockes--;
              tab2_(k1+elements_diagonaux_non_stockes)=i+1;
              coeff_(k1+elements_diagonaux_non_stockes)=0.;
            }
        }
      assert(elements_diagonaux_non_stockes==0);
    }
  morse_matrix_structure_has_changed_=1;
}

/*! @brief Renumbers a matrix to reduce its bandwidth.
 *
 */
void Matrice_Morse_Sym::renumerote() const
{
  Cerr << "Bandwidth of the matrix : " << largeur_de_bande() << finl;
  Cerr << "Renumbering the matrix ..." << finl;
  const Matrice_Morse_Sym& matrice_initial = *this;
  ArrOfInt& tab_iperm = matrice_initial.permutation_inverse();

  // convert a symmetric Morse matrix to a general Morse matrix
  Matrice_Morse matrice2(matrice_initial);
  Matrice_Morse matrice(matrice_initial);

  int mon_ordre = matrice.ordre();
  matrice2.transpose(matrice);
  for (int i=0; i<mon_ordre; i++) matrice2(i, i) = 0.;
  matrice2 += matrice ;

  // compute the permutation to apply
  const int n = mon_ordre;
  matrice2.set_tab1_int32();
  const int* tab1tmp = matrice2.get_tab1_int32().addr();
  const int* tab2tmp = matrice2.get_tab2().addr();
  int init = 1;
  tab_iperm.resize_array(n);
  tab_iperm[0] = 1;

  int* masktmp = new int[n];
  for (int i=0 ; i<n; i++) masktmp[i] = 1;
  const int* mask = (const int*) masktmp;
  const int maskval = 1;
  // GF: changed to n+1 to allow Cholesky on a 1xN mesh
  int* level = new int[n+1];
  int nlev;

  // renumber the nodes
  // subroutine perphn(n,ja,ia,init,iperm,mask,maskval,nlev,riord,levels)
  // SPARSKIT2/ORDERINGS/levset.f
  F77NAME(PERPHN)(&n, tab2tmp, tab1tmp, &init,  mask, &maskval,
                  &nlev, tab_iperm.addr(), level);

  delete []masktmp;
  delete []level;

  ArrOfInt tab_perm(n);
  for (int i=0 ; i<n; i++ ) tab_perm[tab_iperm[i] - 1] = i + 1;

  matrice2 = matrice;

  const double* a = matrice.get_coeff().addr();
  const int* ja = matrice.get_tab2().addr();
  matrice.set_tab1_int32();
  const int* ia = matrice.get_tab1_int32().addr();
  const double* tao = matrice2.get_coeff().addr();
  double* ao = (double*)tao;
  const int* tjao = matrice2.get_tab2().addr();
  int* jao = (int*)tjao;
  matrice2.set_tab1_int32();
  const int* tiao = matrice2.get_tab1_int32().addr();
  int* iao = (int*)tiao;

  const int* perm = tab_perm.addr();
  const int* perm_inv = tab_iperm.addr();//normalement inutile
  const int job = 1;

  // permute matrice2, i.e. compute P A tP
  // only the upper part is permuted here

  // subroutine dperm (nrow,a,ja,ia,ao,jao,iao,perm,qperm,job)
  // SPARSKIT2/FORMATS/unary.f
  F77NAME(DPERM) (&n, a, ja, ia, ao, jao, iao, perm, perm_inv, &job);
  matrice2.set_tab1(matrice2.get_tab1_int32());

  matrice.transpose(matrice2);
  for (int i=0; i<mon_ordre; i++) matrice2(i, i) = 0.;
  matrice2 += matrice;
  matrice.partie_sup(matrice2);
  matrice_renumerotee_.typer("Matrice_Morse_Sym");
  ref_cast(Matrice_Morse_Sym,matrice_renumerotee_.valeur()) = matrice;

  // Build permutation() - provisional, already done with
  int size=permutation_inverse().size_array();
  permutation().resize_array(size);
  for(int i=0; i<size; i++)
    permutation()[permutation_inverse()[i]-1]=i+1;
  Cerr << "Bandwidth of the renumbered matrix : " << ref_cast(Matrice_Morse,matrice_renumerotee_.valeur()).largeur_de_bande() << finl;
  morse_matrix_structure_has_changed_=1;
}

bool Matrice_Morse_Sym::check_symmetric_morse_matrix_structure() const
{
  if ( ! ( check_morse_matrix_structure( ) ) )
    {
      Cerr << "Invalid morse structure" << finl;
      return false;
    }

  const int nb_lines = nb_lignes( );
  for ( int i=0; i<nb_lines; ++i )
    {
      auto k0 = tab1_( i ) - 1;
      auto k1 = tab1_( i + 1 ) - 1;

      for ( auto k=k0; k<k1; ++k )
        {
          int j = tab2_( k ) - 1;
          if  ( j < i )
            {
              Cerr << "( j < i ) : line index : " << i << " and column index " << j << finl;
              return false;
            }
        }
    }
  return true;
}

bool Matrice_Morse_Sym::check_sorted_symmetric_morse_matrix_structure() const
{
  if ( ! ( check_sorted_morse_matrix_structure( ) ) )
    {
      return false;
    }

  const int nb_lines = nb_lignes( );
  for ( int i=0; i<nb_lines; ++i )
    {
      auto k0 = tab1_( i ) - 1;
      auto k1 = tab1_( i + 1 ) - 1;

      for ( auto k=k0; k<k1; ++k )
        {
          int j = tab2_( k ) - 1;
          if  ( j < i )
            {
              Cerr << "( j < i ) : line index : " << i << " and column index " << j << finl;
              return false;
            }
        }
    }
  return true;
}


void Matrice_Morse_Sym::assert_check_symmetric_morse_matrix_structure() const
{
  if (!morse_matrix_structure_has_changed_) return;
#ifndef NDEBUG
  if ( ! ( check_symmetric_morse_matrix_structure( ) ) )
    {
      Cerr << "Error in 'Matrice_Morse_Sym::assert_check_symmetric_morse_matrix_structure( )':" << finl;
      Cerr << "  Exiting..." << finl;
      Process::exit( );
    }
  else
    morse_matrix_structure_has_changed_=0;
#endif
}

void Matrice_Morse_Sym::assert_check_sorted_symmetric_morse_matrix_structure() const
{
  if (!morse_matrix_structure_has_changed_) return;
#ifndef NDEBUG
  if ( ! ( check_symmetric_morse_matrix_structure( ) ) )
    {
      Cerr << "Error in 'Matrice_Morse_Sym::assert_check_sorted_symmetric_morse_matrix_structure( )':" << finl;
      Cerr << "  Exiting..." << finl;
      Process::exit( );
    }
  else
    morse_matrix_structure_has_changed_=0;
#endif
}
