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

#include <Matrice_Morse.h>
#include <Sparskit.h>
#include <unordered_map>
#include <Matrice_Morse_Sym.h>
#include <Check_espace_virtuel.h>
#include <SFichier.h>
#include <Noms.h>
#include <ArrOfBit.h>
#include <Array_tools.h>
#include <TRUSTTrav.h>
#include <TRUSTTrav.h>


Implemente_instanciable_sans_constructeur(Matrice_Morse,"Matrice_Morse",Matrice_Base);

/*! @brief Ecrit les trois tableaux de la structure de stockage Morse sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie& s) le flot de sortie modifie
 */
Sortie& Matrice_Morse::printOn(Sortie& s) const
{
  s << tab1_;
  s << tab2_;
  s << coeff_;
  s << m_ << finl;
  return s;
}

/*! @brief NON CODE
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree& s) le flot d'entree
 * @throws NON CODE
 */
Entree& Matrice_Morse::readOn(Entree& s)
{
  s >> tab1_;
  s >> tab2_;
  s >> coeff_;
  s >> m_;
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  return s;
}

Sortie& Matrice_Morse::imprimer(Sortie& s) const
{
  int n=nb_lignes();
  for(int i=0; i<n; i++)
    {
      s <<i << ": " <<finl;
      s << "--------------------------------" << finl;
      for (auto k=tab1_(i)-1; k<tab1_(i+1)-1; k++)
        {
          s << "("<<(tab2_(k)-1) << "),(" <<coeff_(k)<< ") "
            << " k= " << k << finl;
        }
      s <<finl;
    }
  return s;
}

Sortie& Matrice_Morse::imprimer_formatte(Sortie& s) const
{
  return imprimer_formatte(s,0);
}

Sortie& Matrice_Morse::imprimer_formatte(Sortie& s, int symetrie) const
{
  int numerotation_fortran=(tab1_[0]==1);
  for (int proc=0; proc<Process::nproc(); proc++)
    {
      if (proc==Process::me())
        {
          s << "Matrix morse on the processor " << proc << " : " << finl;
          int n=nb_lignes();
          Noms tab_imp;
          tab_imp.dimensionner(nb_colonnes());
          for(int i=0; i<n; i++)
            {
              for (int k=0; k<nb_colonnes(); k++)
                tab_imp[k]="  .  ";
              if (i<10)
                s <<i << " :" ;
              else
                s <<i << ":" ;

              if (symetrie)
                {
                  for (int j=0; j<i; j++)
                    {
                      for (auto k=tab1_(j)-numerotation_fortran; k<tab1_(j+1)-numerotation_fortran; k++)
                        if (tab2_(k)-numerotation_fortran==i)
                          tab_imp[j] = coeff_(k);
                    }
                  int ligne=tab2_(tab1_(i)-numerotation_fortran)-numerotation_fortran;
                  if (i!=ligne)
                    {
                      Cerr << "Problem detected on this Matrice_Morse_Sym." << finl;
                      Cerr << "The diagonal of the line " << ligne << " must be stored even if it is null." << finl;
                      exit();
                    }
                }
              for (auto k=tab1_(i)-numerotation_fortran; k<tab1_(i+1)-numerotation_fortran; k++)
                if (tab2_(k)+!numerotation_fortran==0)
                  Cerr<<"Line " <<i<< " no coefficient "<<k<<finl;
                else
                  {
                    if (coeff_(k)>=0)
                      tab_imp[tab2_(k)-numerotation_fortran]=" ";
                    else
                      tab_imp[tab2_(k)-numerotation_fortran]="";
                    tab_imp[tab2_(k)-numerotation_fortran] += (Nom)coeff_(k);
                  }
              for(int k=0; k<nb_colonnes(); k++)
                s<<tab_imp[k];
              s<<finl;
            }
        }
      Process::barrier();
    }
  return s;
}

Sortie& Matrice_Morse::imprimer_image(Sortie& s) const
{
  return imprimer_image(s,0);
}

Sortie& Matrice_Morse::imprimer_image(Sortie& s, int symetrie) const
{
  int numerotation_fortran=(tab1_[0]==1);
  for (int proc=0; proc<Process::nproc(); proc++)
    {
      if (proc==Process::me())
        {
          s << "Matrix morse on the processor " << proc << " : " << finl;
          int n=nb_lignes();
          Noms tab_imp;
          tab_imp.dimensionner(nb_colonnes());
          for(int i=0; i<n; i++)
            {
              for (int k=0; k<nb_colonnes(); k++)
                tab_imp[k]="\u2588\u2588";
              if (i<10)
                s <<i << " :" ;
              else
                s <<i << ":" ;

              if (symetrie)
                {
                  for (int j=0; j<i; j++)
                    {
                      for (auto k=tab1_(j)-numerotation_fortran; k<tab1_(j+1)-numerotation_fortran; k++)
                        if (tab2_(k)-numerotation_fortran==i)
                          tab_imp[j] = (std::abs(coeff_(k)) < 1e-20) ? "  " : "\u2592\u2592";
                    }
                  int ligne=tab2_(tab1_(i)-numerotation_fortran)-numerotation_fortran;
                  if (i!=ligne)
                    {
                      Cerr << "Problem detected on this Matrice_Morse_Sym." << finl;
                      Cerr << "The diagonal of the line " << ligne << " must be stored even if it is null." << finl;
                      exit();
                    }
                }
              for (auto k=tab1_(i)-numerotation_fortran; k<tab1_(i+1)-numerotation_fortran; k++)
                if (tab2_(k)+!numerotation_fortran==0)
                  Cerr<<"Line " <<i<< " no coefficient "<<k<<finl;
                else
                  tab_imp[tab2_(k)-numerotation_fortran] = (std::abs(coeff_(k)) < 1e-20) ? "  " : "\u2592\u2592";

              for(int k=0; k<nb_colonnes(); k++)
                s<<tab_imp[k];
              s<<finl;
            }
        }
      Process::barrier();
    }
  return s;
}

void Matrice_Morse::WriteFileMTX(const Nom& name) const
{
  if (Process::is_parallel())
    {
      Cerr << "Warning, matrix market format is not available yet in parallel." << finl;
      return;
    }
  Nom filename(Objet_U::nom_du_cas());
  filename += "_";
  filename += name;
  filename += ".mtx";
  SFichier mtx(filename);
  mtx.precision(14);
  mtx.setf(ios::scientific);
  int rows = nb_lignes();
  Cerr << "Matrix (" << rows << " lines) written into file: " << filename << " ... " << finl;
  mtx << "%%MatrixMarket matrix coordinate real " << (sub_type(Matrice_Morse_Sym, *this) ? "symmetric" : "general") << finl;
  Cerr << "Matrix (" << rows << " lines) written into file: " << filename << finl;
  mtx << "%%matrix" << finl;
  mtx << rows << " " << rows << " " << get_tab1()[rows] << finl;
  for (int row=0; row<rows; row++)
    for (auto j=get_tab1()[row]; j<get_tab1()[row+1]; j++)
      mtx << row+1 << " " << get_tab2()[j-1] << " " << get_coeff()[j-1] << finl;
}

/*! @brief Constructeur par copie d'une Matrice_Morse.
 *
 * Copie de chaque membre donne du paramtre.
 *
 * @param (Matrice_Morse& acopier) la matrice morse a copier
 */
Matrice_Morse::Matrice_Morse(const Matrice_Morse& acopier) :Matrice_Base(),
  tab1_(acopier.tab1_),
  tab2_(acopier.tab2_),
  coeff_(acopier.coeff_),
  m_(acopier.m_),
  symetrique_(0),
  zero_(0)
{
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  is_stencil_up_to_date_ = acopier.is_stencil_up_to_date_ ;
}
/*! @brief Constructeur d'une matrice Morse carree d'ordre n et pouvant stocker au maximum nnz elements non nuls.
 *
 *     Egalement constructeur par defaut car les 2 parametres
 *     ont une valeur par defaut.
 *
 * @param (int n) l'ordre de la matrice carree a construire
 * @param (int nnz) le nombre d'elements non nuls que pourra stocker la matrice.
 */
template<typename _SIZE_>
Matrice_Morse::Matrice_Morse(int n, _SIZE_ nnz) :
  morse_matrix_structure_has_changed_(1), symetrique_(0), zero_(0)
{
  dimensionner(n,nnz), sorted_ = 0;
  is_stencil_up_to_date_ = false ;
}

Matrice_Morse::Matrice_Morse()
{
  dimensionner(0,0);
  morse_matrix_structure_has_changed_=1;
  symetrique_ = 0;
  sorted_ = 0;
  zero_ = 0;
  is_stencil_up_to_date_ = false ;
}


/*! @brief Constructeur d'une matrice Morse avec n lignes et m colonnes pouvant stocker au maximum nnz elements non nuls.
 *
 * @param (int n) le nombre de ligne de la matrice
 * @param (int m) le nombre de colonne de la matrice
 * @param (int nnz) le nombre d'elements non nuls que pourra stocker la matrice.
 */
template<typename _SIZE_>
Matrice_Morse::Matrice_Morse(int n, int m, _SIZE_ nnz):
  morse_matrix_structure_has_changed_(1), symetrique_(0), zero_(0)
{
  dimensionner(n,m,nnz);
  is_stencil_up_to_date_ = false, sorted_ = 0 ;
}


Matrice_Morse::Matrice_Morse(int n, int nnz, const IntLists& voisins,
                             const DoubleLists& valeurs,
                             const DoubleVect& terme_diag)
  :  morse_matrix_structure_has_changed_(1), symetrique_(0) , zero_(0)
{
  dimensionner(n,n,nnz);
  remplir(voisins, valeurs, terme_diag);
  is_stencil_up_to_date_ = false, sorted_ = 0;
}

void Matrice_Morse::set_nb_columns( const int nb_col )
{
  m_ = nb_col;
}

void Matrice_Morse::set_symmetric( const int symmetric )
{
  symetrique_ = symmetric ;
}

/*! @brief Size the matrix with n lines and n columns and nnz zero-values coefficients
 *
 */
template<typename _SIZE_>
void Matrice_Morse::dimensionner(int n, _SIZE_ nnz)
{
  dimensionner(n,n,nnz);
  return ;
}


/*! @brief Redimensionne la matrice creuse en ajoutant eventuellement des coefficients non nuls
 *
 *
 *    Parametre : const IntTab &Ind
 *       Signification : tableau de taille nc * 2
 *                       ou nc est le nombre de couples (i,j)
 *                       pour les indices des nouveaux coefficients
 *
 *
 */
void Matrice_Morse::dimensionner(const IntTab& Ind)
{
  if (Ind.size()==0) return; // On ne fait rien si la structure est vide
  int n_ancien = nb_lignes(), m_ancien = nb_colonnes();

  assert(Ind.nb_dim() == 2);
  assert(Ind.dimension(1) == 2);

  // Calcul du nouveau nombre de lignes
  //   = max (ancien, indices de ligne des nouveaux coeffs)
  //
  // et du nouveau nombre de colonnes
  //   = max (ancien, indices de colonne des nouveaux coeffs)

  int nInd = Ind.dimension(0);
  int n = 0;
  int m = 0;
  for (int i=0; i<nInd; i++)
    {
      if (n < Ind(i,0)) n = Ind(i,0);
      if (m < Ind(i,1)) m = Ind(i,1);
    }
  n++;
  m++;
  if (n < n_ancien) n = n_ancien;
  if (m < m_ancien) m = m_ancien;

  // Copies des anciens tableaux d'indices

  auto tab1_temp(tab1_);

  // Initialisation au nombre de coefficients deja presents a chaque ligne

  tab1_.resize(n+1);
  m_ = m;

  for (int i=1; i<=n_ancien; i++)
    tab1_[i] = tab1_temp[i] - tab1_temp[i-1];
  for (int i=n_ancien+1; i<=n; i++)
    tab1_[i] = 0;

  // Parcourt des indices des nouveaux coeffs pour voir s'ils sont
  // deja presents

  int i_nouveaux = 0;
  for (int i=0; i<nInd; i++)
    {
      int i0 = Ind(i,0);
      int i1 = Ind(i,1) + 1;

      int test_present = 0;

      if (i0 < n_ancien)
        {
          auto kmin = tab1_temp[i0]-1;
          auto kmax = tab1_temp[i0+1]-1;
          for (auto k=kmin; k<kmax; k++)
            if (tab2_[k] == i1)
              {
                test_present = 1;
                break;
              }
        }
      if (!test_present)
        {
          i_nouveaux++;
          tab1_[i0+1] += 1;
        }
    }
  if (i_nouveaux == 0)
    {
      tab1_=tab1_temp;
      return;
    }

  // Nouveau tableau des positions des premiers coeffs de chaque ligne
  tab1_[0] = 1;
  for (int i=1; i<=n; i++)
    tab1_[i] += tab1_[i-1];

  auto nnz_ancien = tab2_.size_array();
  auto nnz = nnz_ancien + i_nouveaux;

  auto tab2_temp(tab2_);
  auto coeff_temp(coeff_);
  tab2_.resize(nnz);
  coeff_.resize(nnz);

  // Recopie des anciens coefficients et de leurs indices
  // de colonne dans les nouveaux tableaux

  tab2_ = -1;
  for (int i=0; i<n_ancien; i++)
    {
      for (auto j1 = tab1_temp[i]-1, j2 = tab1_[i]-1;
           j1 < tab1_temp[i+1]-1;
           j1++, j2++)
        {
          tab2_[j2] = tab2_temp[j1];
          coeff_[j2] = coeff_temp[j1];
        }
    }

  for (int i=0; i<nInd; i++)
    {
      int j0 = Ind(i,0);
      int j1 = Ind(i,1) + 1;
      auto k(tab1_(0));
      for (k=tab1_[j0]-1; tab2_[k] >= 0; k++)
        if (tab2_[k] == j1)
          {
            break;
          }
      if (tab2_[k] < 0)
        {
          tab2_[k] = j1;
          coeff_[k] = 0.0;
        }
    }
  // on remet les coeffs dans l'ordre... pas optimal mais pour voir..
  coeff_=0;
  //
  {
    int nbis=nb_lignes();
    for(int i=0; i<nbis; i++)
      {
        for (auto k=tab1_(i)-1; k<tab1_(i+1)-1; k++)
          {
            for (auto k2=k; k2<tab1_(i+1)-1; k2++)
              {
                int j1=tab2_(k);
                int j2=tab2_(k2);
                if (j1>j2)
                  {
                    tab2_(k)=j2;
                    tab2_(k2)=j1;
                  }
              }
          }
      }
  }
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
}

/*! @brief Size the matrix with n lines, m columns with nnz zero-values coefficients
 *
 */
template<typename _SIZE_>
void Matrice_Morse::dimensionner(int n, int m, _SIZE_ nnz)
{
  tab2_.resize(nnz);
  coeff_.resize(nnz);
  m_=m;

  // on regarde si tab1 a la bonne taille et si tab1[n1]==nnz.
  if ( tab1_.size_array()!=(n+1) || (tab1_[n]-1)!=nnz )
    {
      tab1_.resize(n+1);
      tab1_=1;
    }
  tab1_.resize(n+1);
  tab1_[n]=nnz+1;

  morse_matrix_structure_has_changed_=1, sorted_ = 0;
}

/*! @brief Initialisation a la matrice unite (modif MT)
 *
 */
void Matrice_Morse::unite()
{
  coeff_ = 0.0;
  int i,n  = ordre();
  for (i=0; i<n; i++)
    operator()(i,i) = 1.0;
}

/*! @brief Renvoie l'ordre de la matrice: - le nombre de lignes si la matrice est carree
 *
 *      - 0 sinon
 *
 * @return (int) l'ordre de la matrice
 */
int Matrice_Morse::ordre() const
{
  if(nb_lignes()==nb_colonnes())
    return nb_lignes();
  else
    return 0;
}

/*! @brief Method to check/clean the Matrice_Morse matrix: -Suppress coefficient defined several times
 *
 *  -elim_coeff_nul=0, on ne supprime pas les coefficients nuls de la matrice
 *  -elim_coeff_nul=1, on supprime les coefficients nuls de la matrice
 *  -elim_coeff_nul=2, on supprime les coefficients nuls et quasi-nuls de la matrice
 *
 */
void Matrice_Morse::compacte(int elim_coeff_nul)
{
  int n=nb_lignes();
  int coeff_nuls=0;
  int coeff_quasi_nuls=0;
  auto tab_elim_coeff(tab2_); // Possibly BigArrOfInt
  tab_elim_coeff = 0;
  if (elim_coeff_nul)
    {
      ArrOfDouble tab_coeff_max(n);
      tab_coeff_max = 0.;
      // Recherche des coefficients nuls hors diagonale a supprimer de la matrice morse
      {
        ArrOfInt tab_cnt(1);
        tab_cnt = 0;
        auto tab1 = tab1_.view_ro();
        CDoubleArrView coeff = coeff_.view_ro();
        DoubleArrView coeff_max = tab_coeff_max.view_rw();
        auto elim_coeff = tab_elim_coeff.view_rw();
        IntArrView cnt = tab_cnt.view_rw();
        Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n), KOKKOS_LAMBDA(const int i)
        {
          auto k1 = tab1(i)-1;
          auto k2 = tab1(i+1)-1;
          for (auto k = k1; k < k2; k++)
            {
              double abs_c = Kokkos::fabs(coeff(k));
              if (abs_c > coeff_max(i)) coeff_max(i) = abs_c;
              if (coeff(k) == 0)
                {
                  Kokkos::atomic_add(&cnt(0), 1);
                  elim_coeff(k) = 1;
                }
            }
        });
        end_gpu_timer(__KERNEL_NAME__);
        coeff_nuls = tab_cnt(0);
      }

      if (elim_coeff_nul==2)
        {
          // Recherche des coefficients quasi nuls hors diagonale (1.e-12 plus petit que le coefficient le plus grand de la ligne) a supprimer de la matrice morse
          const double eps = Objet_U::precision_geom;
          ArrOfInt tab_cnt(1);
          tab_cnt = 0;
          auto tab1 = tab1_.view_ro();
          CDoubleArrView coeff = coeff_.view_ro();
          CDoubleArrView coeff_max = tab_coeff_max.view_ro();
          IntArrView elim_coeff = tab_elim_coeff.view_rw();
          IntArrView cnt = tab_cnt.view_rw();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n), KOKKOS_LAMBDA(const int i)
          {
            double cm = coeff_max(i);
            if (!est_egal(cm, 0., eps) && cm < 1e10)
              {
                auto k1 = tab1(i) - 1;
                auto k2 = tab1(i + 1) - 1;
                for (auto k = k1; k < k2; k++)
                  if (coeff(k) != 0 && est_egal(Kokkos::fabs(coeff(k)) / cm, 0., eps))
                    {
                      Kokkos::atomic_add(&cnt(0), 1);
                      elim_coeff(k) = 1;
                    }
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
          coeff_quasi_nuls = tab_cnt(0);
        }
    }
  // Recherche des coefficients doublons
  int nb_doublons=0;
  {
    auto tab1 = tab1_.view_ro();
    CIntArrView tab2 = tab2_.view_ro();
    CDoubleArrView coeff = coeff_.view_ro();
    IntArrView elim_coeff = tab_elim_coeff.view_rw();
    ArrOfInt tab_doublons(1);
    tab_doublons = 0;
    ArrOfInt tab_error(1);
    tab_error = 0;
    IntArrView doublons = tab_doublons.view_rw();
    IntArrView error = tab_error.view_rw();
    Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n), KOKKOS_LAMBDA(const int i)
    {
      auto k1 = tab1(i)-1;
      auto k2 = tab1(i+1)-1;
      int jmax = -1; // Highest column of a coefficient in the line i
      for (auto k = k1; k < k2; k++)
        {
          int j = tab2(k)-1;
          if (j > jmax)
            jmax = j;
          else
            {
              // Found a column j lower than jmax, check if not defined before:
              for (auto kk = k-1; kk >= k1; kk--)
                {
                  int jj = tab2(kk)-1;
                  if (jj == j)
                    {
                      // Already defined!
                      Kokkos::atomic_add(&doublons(0), 1);
                      elim_coeff(k) = 1;
                      // Check if same coefficients:
                      if (coeff(kk) != coeff(k))
                        Kokkos::atomic_add(&error(0), 1);
                      break;
                    }
                }
            }
        }
    });
    end_gpu_timer(__KERNEL_NAME__);
    nb_doublons = tab_doublons(0);
    if (tab_error(0))
      {
        Cerr << "Error in a Matrix Morse: duplicate entries with different values!" << finl;
        exit();
      }
  }

  auto nnz(tab1_(0));
  nnz=0;
  if (nb_doublons || coeff_nuls || coeff_quasi_nuls)
    {
      // Step 1: Count kept entries per row (parallel_for over rows)
      ArrOfInt tab_kept_per_row(n);
      {
        auto tab1 = tab1_.view_ro();
        CIntArrView elim_coeff = tab_elim_coeff.view_ro();
        IntArrView kept_per_row = tab_kept_per_row.view_wo();
        Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n), KOKKOS_LAMBDA(const int i)
        {
          int count = 0;
          auto k1 = tab1(i)-1;
          auto k2 = tab1(i+1)-1;
          for (auto k = k1; k < k2; k++)
            if (!elim_coeff(k)) count++;
          kept_per_row(i) = count;
        });
        end_gpu_timer(__KERNEL_NAME__);
      }

      // Step 2: Save old tab1_ (needed for source offsets in scatter step)
      auto old_tab1(tab1_);

      // Step 3: Update tab1_ via prefix scan (updates tab1_(1..n), tab1_(0)=1 unchanged)
      using tab1_scan_t = decltype(nnz);
      {
        auto tab1 = tab1_.view_rw();
        CIntArrView kept_per_row = tab_kept_per_row.view_ro();
        Kokkos::parallel_scan(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n), KOKKOS_LAMBDA(const int i, tab1_scan_t& update, const bool final)
        {
          update += kept_per_row(i);
          if (final) tab1(i+1) = update + 1;
        });
        end_gpu_timer(__KERNEL_NAME__);
      }

      // Step 4: Out-of-place scatter of coeff_ and tab2_ to new positions (parallel_for over rows)
      // Safe because new_pos(i) <= old_pos(i) always, and rows are processed independently
      nnz = tab1_[n] - 1;
      auto new_coeff(coeff_);
      auto new_tab2(tab2_);
      {
        auto tab1 = tab1_.view_ro();
        auto old_tab1_ro = old_tab1.view_ro();
        CDoubleArrView coeff_src = coeff_.view_ro();
        CIntArrView tab2_src = tab2_.view_ro();
        DoubleArrView coeff_dst = new_coeff.view_wo();
        IntArrView tab2_dst = new_tab2.view_wo();
        CIntArrView elim_coeff = tab_elim_coeff.view_ro();
        Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n), KOKKOS_LAMBDA(const int i)
        {
          auto new_pos = tab1(i) - 1;
          auto k1 = old_tab1_ro(i)-1;
          auto k2 = old_tab1_ro(i+1)-1;
          for (auto k = k1; k < k2; k++)
            if (!elim_coeff(k))
              {
                coeff_dst(new_pos) = coeff_src(k);
                tab2_dst(new_pos) = tab2_src(k);
                new_pos++;
              }
        });
        end_gpu_timer(__KERNEL_NAME__);
      }

      // Step 5: Copy compacted data back
      {
        auto tab2 = tab2_.view_rw();
        auto coeff = coeff_.view_rw();
        CIntArrView new_tab2_ro = new_tab2.view_ro();
        CDoubleArrView new_coeff_ro = new_coeff.view_ro();
        Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, nnz), KOKKOS_LAMBDA(const int i)
        {
          tab2(i) = new_tab2_ro(i);
          coeff(i) = new_coeff_ro(i);
        });
        end_gpu_timer(__KERNEL_NAME__);
      }
    }
  else
    {
      nnz = tab1_[n] - 1;
    }

  // On redimensionne les tableaux
  tab2_.resize(nnz);
  coeff_.resize(nnz);

  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  assert_check_morse_matrix_structure( );
}

/*! @brief Operateur d'affectation d'une Matrice_Morse dans une autre Matrice_Morse.
 *
 * @param (Matrice_Morse& a) la partie droite de l'affectation
 */
Matrice_Morse& Matrice_Morse::operator=(const Matrice_Morse& a )
{
  tab1_ = a.get_tab1();
  tab2_ = a.get_tab2();
  coeff_ = a.get_coeff();
  m_=a.nb_colonnes();
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  is_stencil_up_to_date_=a.is_stencil_up_to_date();
  return(*this);
}

/*! @brief *this = a transposee.
 *
 * @param (Matrice_Morse& a) la matrice a transposee
 */
Matrice_Morse& Matrice_Morse::transpose(const Matrice_Morse& a)
{
  int n=a.nb_lignes();
  int jk=nb_lignes();
  int job=1;
  int ipos=1;
  int m=a.nb_colonnes();
  int l=nb_lignes();
  if(m!=jk)
    {
      Cerr << "Matrice_Morse::transpose bad dimensions" << finl;
      exit();
    }
  m=a.nb_lignes();
  l=nb_colonnes();
  if(m!=l)
    {
      Cerr << "Matrice_Morse::transpose bad dimensions" << finl;
      exit();
    }

  for(int i=0; i<=jk; i++ ) tab1_[i] = 0 ;
  for(int i=0; i<n; i++)
    {
      for(auto k=a.tab1_[i]-1; k<a.tab1_[i+1]-1; k++)
        {
          int j = a.tab2_[k] ;
          tab1_[j] = tab1_[j] +1 ;
        }
    }
  tab1_[0] = ipos ;
  for(int i=0; i<jk; i++) tab1_[i+1] = tab1_[i] + tab1_[i+1] ;
  for(int i=0; i<n; i++)
    {
      for(auto k=a.tab1_[i]-1; k<a.tab1_[i+1]-1; k++)
        {
          int j = a.tab2_[k]-1 ;
          auto next = tab1_[j] ;
          if (job == 1) coeff_[next-1] = a.coeff_[k] ;
          tab2_[next-1] = i+1 ;
          tab1_[j]    = next+1 ;
        }
    }
  for(int i=jk-1; i>=0; i--) tab1_[i+1] = tab1_[i] ;
  tab1_[0] = ipos ;

  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  return(*this);
}


//A=x*A avec x une matrice diagonale stockee dans un vecteur
//la meme methode peut etre utilisee pour stocke le resultat dans
//un autre matrice que la matrice initiale
Matrice_Morse& Matrice_Morse::diagmulmat(const DoubleVect& x)
{
  int m=nb_lignes();
  int l=0;
  int n=x.size_array();
  if(n!=m)
    {
      Cerr << "Matrice_Morse::diagmulmat bad dimensions" << finl;
      exit();
    }
  F77NAME(DIAMUA)(&m ,&l,
                  coeff_.addr(),tab2_.addr(),reinterpret_cast<const int*>(tab1_.addr()),x.addr(),
                  coeff_.addr(),tab2_.addr(),reinterpret_cast<int*>(tab1_.addr()));
  return(*this);
}

//extraction de la partie superieure d'une matrice morse
//la matrice resultat est celle appelante
Matrice_Morse& Matrice_Morse::partie_sup(const Matrice_Morse& a)
{
  int m=nb_lignes();
  int n=a.nb_lignes();
  if(m!=n)
    {
      Cerr << "Matrice_Morse::partie_sup : bad dimensions m!=n." << finl;
      exit();
    }
  double t;
  auto ko(tab1_(0));
  auto kfirst (ko);
  auto kdiag(ko);
  ko = -1;
  for(int i=0; i< n; i++)
    {
      kfirst = ko + 1 ;
      kdiag = -1 ;
      for(auto k = a.tab1_[i]-1; k< a.tab1_[i+1]-1; k++)
        {
          if (a.tab2_[k]-1 >= i)
            {
              ko++ ;
              coeff_[ko] = a.coeff_[k] ;
              tab2_[ko] = a.tab2_[k] ;
              if (a.tab2_[k] == i) kdiag = ko ;
            }
        }
      if (kdiag != -1 && kdiag != kfirst)
        {
          t = coeff_[kdiag] ;
          coeff_[kdiag] = coeff_[kfirst] ;
          coeff_[kfirst] = t ;
          { int ktmp = tab2_[kdiag] ; tab2_[kdiag] = tab2_[kfirst] ; tab2_[kfirst] = ktmp ; }
        }
      tab1_[i] = kfirst+1 ;
    }
  auto nnz = (ko + 1) ;

  tab1_[n] = (nnz) + 1 ;
  tab2_.resize( nnz );
  coeff_.resize( nnz );
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  return(*this);
}

/*! @brief Operation de multiplication-accumulation (saxpy) matrice vecteur.
 *
 * Operation: resu = resu + A*x
 *
 */
DoubleVect& Matrice_Morse::ajouter_multvect_(const DoubleVect& tab_x,DoubleVect& tab_resu) const
{
  assert_check_morse_matrix_structure();
  const int n = tab1_.size_array() - 1;
  assert(tab_x.size_array() == nb_colonnes());
  // Test dans cet ordre car l'attribut size() peut etre invalide:
  assert(tab_resu.size_array() == n || tab_resu.size() == n);
  // If matrix, x, resu are on device, we compute on the device to avoid expensive copy during TRUST GCP:
  if (tab_x.isDataOnDevice() && tab_resu.isDataOnDevice() && coeff_.isDataOnDevice())
    {
      //if (tab_x.line_size()>1) Process::exit("line_size>1 pour x dans Matrice_Morse::ajouter_multvect_");
      // Faster implementation on GPU (ToDo Kokkos: future, use Kokkos kernel?)
      auto tab1 = tab1_.view_ro();
      CIntArrView tab2 = tab2_.view_ro();
      CDoubleArrView coeff = coeff_.view_ro();
      CDoubleArrView x = tab_x.view_ro();
      DoubleArrView resu = tab_resu.view_rw();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__),
                           Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(
                             const int i)
      {
        auto start = tab1(i)-1;
        auto end = tab1(i + 1)-1;
        double tmp {};

        for (auto k = start; k < end; k++)
          {
            int j = tab2(k) - 1;
            tmp+= coeff(k) * x(j);
          }
        resu(i) += tmp;
      });
      end_gpu_timer(__KERNEL_NAME__);
    }
  else
    {
      tab_x.ensureDataOnHost();
      tab_resu.ensureDataOnHost();
      coeff_.ensureDataOnHost();
      // Fast CPU (old) implementation with pointer:
      const DoubleVect& x = tab_x;
      DoubleVect& resu = tab_resu;
      const auto *tab1_ptr = tab1_.addr() + 1;
      const int *tab2_ptr = tab2_.addr();
      const double *coeff_ptr = coeff_.addr();
      const double *x_fortran = x.addr() - 1; // Pour indexer x avec un indice fortran
      auto k_fortran = 1; // indice fortran dans tab2 et coeff
      for (int i = 0; i < n; i++, tab1_ptr++)
        {
          const auto kmax = *tab1_ptr; // tab1_[i+1] = indice fortran dans tab2_
          assert(kmax >= k_fortran && kmax <= tab2_.size_array() + 1);
          double t = resu[i];
          assert(k_fortran == tab1_[i] && tab2_ptr == tab2_.addr() + (k_fortran - 1));
          for (; k_fortran < kmax; k_fortran++, tab2_ptr++, coeff_ptr++)
            {
              int colonne = *tab2_ptr; // indice fortran
              assert(colonne >= 1 && colonne <= nb_colonnes());
              t += (*coeff_ptr) * x_fortran[colonne];
            }
          resu[i] = t;
        }
    }
  return tab_resu;
}

// Multiplication de la matrice par un vecteur x en prenant uniquement les items reels non communs pour x
ArrOfDouble& Matrice_Morse::ajouter_multvect_(const ArrOfDouble& x,ArrOfDouble& resu,ArrOfInt& est_reel_pas_com) const
{
  ToDo_Kokkos("critical ?");
  assert_check_morse_matrix_structure( );
  int n = nb_lignes();

  assert(nb_colonnes()==x.size_array());
  assert(n==resu.size_array());
  for(int i=0; i<n; i++)
    {
      double t = 0.0;
      for (auto k=tab1_(i)-1; k<tab1_(i+1)-1; k++)
        {
          int j=tab2_(k)-1;
          if (est_reel_pas_com[j]) t += coeff_(k)*x[j];
        }
      resu[i] += t ;
    }
  return resu;
}

/*! @brief Operation de multiplication-accumulation (saxpy) matrice matrice (matrice X representee par un tableau)
 *
 *     Operation: RESU = RESU + A*X
 *
 * @param (DoubleTab& x) la matrice a multiplier
 * @param (DoubleTab& resu) la matrice resultat de l'operation
 * @return (DoubleTab&) la matrice resultat de l'operation
 */
DoubleTab& Matrice_Morse::ajouter_multTab_(const DoubleTab& x,DoubleTab& resu) const
{

  if ( (x.nb_dim() == 1) && (resu.nb_dim() == 1))
    {
      ajouter_multvect(x,resu);
      return resu;
    }

  assert_check_morse_matrix_structure( );
  int nb_comp = x.dimension(1);

  assert(resu.dimension(1) == nb_comp);
  double* t= new double[nb_comp];
  int ncomp;
  int n=nb_lignes();
  for(int i=0; i<n; i++)
    {
      for (ncomp=0; ncomp<nb_comp; ncomp++)
        t[ncomp] = 0.0;
      for (auto k=tab1_(i)-1; k<tab1_(i+1)-1; k++)
        for (ncomp=0; ncomp<nb_comp; ncomp++)
          t[ncomp] += coeff_(k)*x(tab2_(k)-1,ncomp);
      for (ncomp=0; ncomp<nb_comp; ncomp++)
        resu(i,ncomp) += t[ncomp] ;
    }
  delete [] t;
  return resu;
}


/*! @brief Operation de multiplication-accumulation (saxpy) matrice vecteur, par la matrice transposee.
 *
 *     Operation: resu = resu + A^{T}*x
 *
 * @param (DoubleVect& x) le vecteur a multiplier
 * @param (DoubleVect& resu) le vecteur resultat de l'operation
 * @return (DoubleVect&) le vecteur resultat de l'operation
 */
DoubleVect& Matrice_Morse::ajouter_multvectT_(const DoubleVect& x,DoubleVect& resu) const
{
  assert_check_morse_matrix_structure( );

  int n=nb_lignes();
  for(int i=0; i<n; i++)
    {
      double xi = x(i);
      for (auto k=tab1_(i)-1; k<tab1_(i+1)-1; k++)
        resu(tab2_(k)-1) += coeff_(k) * xi;
    }
  return resu;
}

// Multiplication de la tranposee de la matrice par un vecteur x en prenant uniquement les items reels non communs
ArrOfDouble& Matrice_Morse::ajouter_multvectT_(const ArrOfDouble& x,ArrOfDouble& resu,ArrOfInt& est_reel_pas_com) const
{
  assert_check_morse_matrix_structure( );
  int n=nb_lignes();

  assert(n==x.size_array());
  assert(nb_colonnes()==resu.size_array());
  for(int i=0; i<n; i++)
    {
      if (est_reel_pas_com[i])
        {
          double xi = x[i];
          for (auto k=tab1_(i)-1; k<tab1_(i+1)-1; k++)
            resu[tab2_(k)-1] += coeff_(k) * xi;
        }
    }
  return resu;
}

/*! @brief Fonction (hors classe) amie de la classe Matrice_Morse Addition de 2 matrices au format Morse.
 *
 *     Operation: renvoie (A+B)
 *
 * @param (Matrice_Morse& A) une matrice au format Morse
 * @param (Matrice_Morse& B) une matrice au format Morse
 * @return (Matrice_Morse) le resultat de l'operation
 */
Matrice_Morse operator+(const Matrice_Morse& A , const Matrice_Morse& B )
{
  int nrow=A.nb_lignes();
  int ncol=A.nb_colonnes();
  Matrice_Morse C;
  // PL: avant de dimensionner a nzmax on verifie si A et B n'ont pas la meme structure par hasard...
  // Cela evite un pic memoire provoque par l'addition de matrices dans Equation_base::dimensionner_matrice
  auto nzmax = A.has_same_morse_matrix_structure(B) ? A.nb_coeff() : A.nb_coeff()+B.nb_coeff();
  C.dimensionner(nrow, ncol, nzmax);
  // Algorithm (per row i):
  //   1. Collect entries from row i of A and B into a small temporary buffer
  //   2. Sort by column index
  //   3. Merge duplicate columns (accumulate coefficients)
  //   4. Write result into C and advance c_tab1
  //
  // Time: O((nnz_A + nnz_B) * log(max_nnz_per_row))  [sort dominates]
  // Space: O(max_nnz_per_row_A + max_nnz_per_row_B)   [reused buffer]
  //
  // ToDo: Kokkos parallel version for GPU once CPU version is validated
  const auto& a_tab1 = A.get_tab1();
  const auto& a_tab2 = A.get_tab2();
  const auto& a_coeff = A.get_coeff();
  const auto& b_tab1 = B.get_tab1();
  const auto& b_tab2 = B.get_tab2();
  const auto& b_coeff = B.get_coeff();
  auto& c_tab1 = C.get_set_tab1();
  auto& c_tab2 = C.get_set_tab2();
  auto& c_coeff = C.get_set_coeff();

  using idx_t = std::remove_reference_t<decltype(c_tab1[0])>;
  idx_t nnz_c = 0; // running count of non-zeros written into C (0-based offset)
  c_tab1[0] = 1;    // 1-based (Morse/Fortran convention)

  std::unordered_map<int, idx_t> col_to_pos;
  col_to_pos.reserve(256);
  for (int i = 0; i < nrow; ++i)
    {
      col_to_pos.clear();

      // Step 1: copy A row i into C, recording each column's position
      for (auto k = a_tab1[i] - 1; k < a_tab1[i + 1] - 1; ++k)
        {
          c_tab2[nnz_c]      = (int)a_tab2[k];
          c_coeff[nnz_c]     = a_coeff[k];
          col_to_pos[(int)a_tab2[k]] = nnz_c;
          ++nnz_c;
        }

      // Step 2: merge B row i — accumulate if column already in C, else append
      for (auto k = b_tab1[i] - 1; k < b_tab1[i + 1] - 1; ++k)
        {
          const int col = (int)b_tab2[k];
          auto it = col_to_pos.find(col);
          if (it != col_to_pos.end())
            c_coeff[it->second] += b_coeff[k]; // column shared with A: accumulate
          else
            {
              c_tab2[nnz_c]  = col;
              c_coeff[nnz_c] = b_coeff[k];
              col_to_pos[col] = nnz_c;
              ++nnz_c;
            }
        }

      c_tab1[i + 1] = nnz_c + 1; // 1-based pointer to start of next row
    }
  const auto nnz = C.tab1_[nrow] - 1;
  C.get_set_tab2().resize(nnz);
  C.get_set_coeff().resize(nnz);
  C.morse_matrix_structure_has_changed_ = 1, C.sorted_ = 0;
  return(C);
}

bool Matrice_Morse::has_same_morse_matrix_structure(const Matrice_Morse& A) const
{
  int nrow = A.nb_lignes();
  for (int i = 0; i < nrow; i++)
    if (tab1_(i) != A.tab1_(i))
      return false;
  auto ncoeff = tab2_.size_array(), ncoeff_A = A.tab2_.size_array();
  if (ncoeff != ncoeff_A) return false;

  for (auto i = 0; i < ncoeff; i++)
    if (tab2_(i) != A.tab2_(i))
      return false;
  return true;
}

/*! @brief Calcule la solution du systeme lineaire: A * solution = secmem.
 *
 * La methode utilisee est GMRES preconditionnee avec ILUT.
 *   ATTENTION: cette methode n'a vraisemblablement jamais ete testee en parallele
 *
 * @param (DoubleVect& secmem) le second membre du systeme lineaire
 * @param (DoubleVect& solution) la solution du systeme
 * @param (double coeff_seuil)
 * @return (int) renvoie toujours 1
 * @throws Erreur dans ilut 'matrix may be wrong' dixit SAAD
 * @throws Erreur dans ilut : debordement dans L
 * @throws Erreur dans ilut : debordement dans U
 * @throws Valeur illegale de lfil : sans doute ecrasement memoire
 * @throws Ligne vide rencontree
 * @throws Pivot nul rencontre ! au pas
 * @throws Il s'est passe quelque chose de bizarre : je prefere tout arreter.
 */
// Delegates to the 4-arg version with max_iter=-1 (retry-on-failure mode, maxits=ordre())
int Matrice_Morse::inverse(const DoubleVect& secmem, DoubleVect& solution,
                           double coeff_seuil) const
{
  return inverse(secmem, solution, coeff_seuil, -1);
}

// Solves A*solution=secmem using ILUT-preconditioned PGMRES.
// max_iter<0: use ordre() as iteration limit and retry with stronger preconditioner on failure (returns 1).
// max_iter>=0: use max_iter as iteration limit and return 0 on failure (used by hyperbolic implicit).
int Matrice_Morse::inverse(const DoubleVect& secmem, DoubleVect& solution,
                           double coeff_seuil, int max_iter) const
{
  if (Process::is_parallel())
    {
      Cerr << "Matrice_Morse::inverse has never been tested in parallel" << finl;
      Cerr << "Try 'Solveur Gmres { diag }' or 'Solveur Petsc Gmres { precond diag { } }'" << finl;
      Cerr << "instead of 'Solveur Gmres { }' which is not parallelized yet." << finl;
      exit();
    }

  const bool retry_on_failure = (max_iter < 0);

  DoubleVect toto(secmem);
  int prems=1;                         // recompute L and U only when prems=1
  int lf_min = 10;
  int lf = std::min(lf_min, ordre()/2); // fill level for ILUT
  int nn = ordre();
  int ima = std::min(lf_min, ordre()/2); // Krylov space dimension
  IntVect ju, jlu;
  DoubleVect alu, vv;
  DoubleVect Resini(toto);

  int ie=1;
  auto n2 = nb_coeff()+(2*lf*nn); // number of non-zeros in LU

  double r, coeff_seuilr;

precond:
  if (prems)
    {
      int iw = (int)(n2 + 2);
      ju.resize(nn);
      jlu.resize(iw);
      alu.resize(iw);
      double to = 1.e-10; // drop tolerance for ILUT
      DoubleVect w(nn+1);
      IntVect jw(2*nn);
      set_tab1_int32();
      F77NAME(ILUT)(&nn, coeff_.addr(), tab2_.addr(), get_tab1_int32().addr(), &lf,
                    &to, alu.addr(), jlu.addr(), ju.addr(),
                    &iw, w.addr(), jw.addr(), &ie);
      switch(ie)
        {
        case  0:
          break;
        case -1:
          Cerr << "Error in ilut 'matrix may be wrong' dixit SAAD" << finl;
          exit();
          break;
        case -2:
          Cerr << "Error in ilut : overflow in L" << finl;
          exit();
          break;
        case -3:
          Cerr << "Error in ilut : overflow in U" << finl;
          exit();
          break;
        case -4:
          Cerr << "Illegal value for lfil : it may be a memory trouble" << finl;
          exit();
          break;
        case -5:
          Cerr << "Empty line met" << finl;
          exit();
          break;
        default:
          Cerr << "Pivot null met ! at step " << ie << finl;
          exit();
          break;
        }
      prems=0;
    }

  vv.resize(nn*(ima+1));
  assert_espace_virtuel_vect(solution);
  multvect(solution, Resini);
  Resini -= toto;
  r = mp_prodscal(Resini, Resini);
  r = sqrt(r);
  Cout << " Initial residu : " << r << finl;
  coeff_seuilr = (r == 0) ? DMAXFLOAT : coeff_seuil/r;
  Resini = toto;
  int minits = 10;
  int maxits = std::max(minits, retry_on_failure ? nn : max_iter);
  int io = 0;
  F77NAME(PGMRES)(&nn, &ima, toto.addr(), solution.addr(), vv.addr(), &coeff_seuilr,
                  &maxits, &io, coeff_.addr(), tab2_.addr(), get_tab1_int32().addr(),
                  alu.addr(), jlu.addr(), ju.addr(), &ie);
  switch(ie)
    {
    case 0:
      Cout << "     ** PGMRES has converged **" << finl;
      break;
    case 1:
      Cout << "     ** No convergence after " << maxits << " iterations **" << finl;
      if (retry_on_failure)
        {
          toto = Resini;
          if (lf < 50)
            {
              lf += 5;
              Cerr << "  The degree of the preconditioning matrix LU is increased: " << lf << finl;
              n2 = (int)tab2_.size_array()+(2*lf*nn);
              prems = 1;
              goto precond;
            }
        }
      else
        return 0;
      break;
    case -1:
      Cerr << "Convergence after 0 iterations !! 'stationnary state may be obtained'" << finl;
      break;
    default:
      Cerr << "Something abnormal has happened : it is preferable to stop." << finl;
      exit();
    }
  return 1;
}


/*! @brief Operateur de multiplication d'une matrice par un vecteur: scaling des lignes de la matrice par les coefficients
 *
 *     correspondants du vecteur passe en parametre.
 *     A *= x, effectue les scaling suivants:
 *       A(i,:) = A(i,:) * x(i), pour toutes les lignes i de A
 *
 * @param (DoubleVect& x) vecteur de scaling
 * @return (Matrice_Morse&) le resultat de l'operation (*this)
 */
Matrice_Morse& Matrice_Morse::operator *=(const DoubleVect& x)
{
  for(int i = 0; i<nb_lignes(); i++)
    for(auto k = tab1_(i)-1; k<tab1_(i+1)-1; k++)
      coeff_(k) *= x(i);
  return *this;
}


/*! @brief Affecte le produit de 2 matrices Morse A et B a l'objet (this).
 *
 * Operation: this = A * B
 *
 * @param (Matrice_Morse& A) une matrice au format Morse
 * @param (Matrice_Morse& B) une matrice au format Morse
 * @return (Matrice_Morse&) le resultat de l'operation (*this)
 */
Matrice_Morse& Matrice_Morse::affecte_prod(const Matrice_Morse& a, const Matrice_Morse& b)
{
  int nrow= a.nb_lignes();                // nb de lignes de A
  int ncol= b.nb_colonnes();                // nb de colonnes de B
  //assert(nrow==ncol);
  // Jloi non?
  assert(a.nb_colonnes()==b.nb_lignes());
  tab1_.resize(nrow+1);
  m_ = ncol;
  int job = 1 ;                      // on recupere tout (tab1, tab2, coeff de matrice_resu)
  auto  nzmax = nb_coeff();                // nb de valeurs maximales de la matrice resultante
  if(nzmax==0)
    {
      nzmax=a.nb_coeff();
      tab2_.resize(nzmax);
      coeff_.resize(nzmax);
      assert(nzmax==nb_coeff());
    }
  IntVect iw(ncol+1);                        // tableau de travail
  double scal=0. ;
  int ii, jj ;
  int values = 0;
  if (job != 0) values = 1 ;
  auto len = -1 ;
  tab1_[0] = 1  ;
  iw = -1 ;
  for(ii=0; ii< nrow; ii++)
    {
      for(auto ka=a.tab1_[ii]-1; ka < a.tab1_[ii+1]-1; ka++)
        {
          if (values == 1) scal = a.coeff_[ka] ;
          jj   = a.tab2_[ka] - 1 ;
          for (auto kb=b.tab1_[jj]-1; kb < b.tab1_[jj+1]-1; kb++)
            {
              int jcol = b.tab2_[kb] -1 ;
              int jpos = iw[jcol]  ;
              if (jpos == -1)
                {
                  len++ ;
                  if (len > nzmax-1)
                    {
                      // Cerr << "Matrice_Morse::affect_prod len > nzmax -1 " << nzmax << finl;
                      nzmax *= 2;
                      coeff_.resize(nzmax);
                      tab2_.resize(nzmax);
                    }
                  tab2_[len] = jcol + 1 ;
                  iw[jcol]= (int)len ;
                  if (values == 1) coeff_[len]  = scal*b.coeff_[kb] ;
                }
              else
                {
                  if (values == 1) coeff_[jpos] += scal*b.coeff_[kb] ;
                }
            }
        }

      for (auto k=tab1_[ii]-1; k < len+1 ; k++) iw[tab2_[k]-1] = -1 ;
      tab1_[ii+1] = (len+1) + 1 ;
    }

  coeff_.resize(tab1_[nrow]-1);
  tab2_.resize(tab1_[nrow]-1);
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  return *this;
}



/*! @brief Fonction (hors classe) amie de la classe Matrice_Morse Scaling de la matrice par un scalaire: multiplie tous
 *
 *     les elements de la matrice par un scalaire.
 *     Operation: renvoie x*A
 *
 * @param (double x) valeur du scaling
 * @param (Matrice_Morse& B) une matrice au format Morse
 * @return (Matrice_Morse) le resultat de l'operation
 */
Matrice_Morse operator *(double x , const Matrice_Morse& A)
{
  Matrice_Morse mat_res(A);
  mat_res.coeff_*=x;
  return(mat_res);
}

/*! @brief Operateur de negation unaire, renvoie l'opposee de la matrice: - A Appelle operator*(double,const Matrice_Morse&)
 *
 * @return (Matrice_Morse) le resultat de l'appel a operator*(double,const Matrice_Morse&)
 */
Matrice_Morse Matrice_Morse::operator -() const
{
  return((-1)*(*this));
}


/*! @brief NE FAIT RIEN
 *
 * @param (Matrice_Morse&) une matrice morse
 * @return (Matrice_Morse&) renvoie toujours *this
 */
Matrice_Morse& Matrice_Morse::operator +=(const Matrice_Morse& A)
{
  // PL: Avant de verifier de faire des operations couteuses en RAM, on verifie
  // si ce n'est pas la meme structure:
  if (has_same_morse_matrix_structure(A))
    {
      auto size = A.nb_coeff();
      const auto& coeff = A.get_coeff();
      for (auto i=0; i<size; i++)
        coeff_(i)+=coeff(i);
    }
  else
    {
      *this = *this + A;
      morse_matrix_structure_has_changed_ = 1, sorted_ = 0;
    }
  return *this;
}


/*! @brief Operateur de multiplication (de tous les elements) d'une matrice par un scalaire.
 *
 *     Operation: A = x * A
 *
 * @param (double x) le parametre de scaling
 * @return (Matrice_Morse&) le resultat de l'operation (*this)
 */
Matrice_Morse& Matrice_Morse::operator *=(double x )
{
  scale( x );
  return(*this);
}

void Matrice_Morse::scale( const double x )
{
  coeff_ *= x;
}

void Matrice_Morse::get_stencil( Stencil& stencil ) const
{
  assert_check_morse_matrix_structure( );

  if( is_stencil_up_to_date_ )
    {
      stencil = stencil_;
      return;
    }


  stencil.resize( 0, 2 );
  auto nnz = tab2_.size_array();
  stencil.resize(nnz, 2);


  ArrOfInt tmp;


  decltype(nnz) compteur = 0;

  const int nb_lines = nb_lignes( );
  for ( int i=0; i<nb_lines; ++i )
    {
      auto k0   = tab1_( i ) - 1;
      auto k1   = tab1_( i + 1 ) - 1;
      const auto size = k1 - k0;
      const int  size_int = (int)size;

      tmp.resize_array( 0 );
      tmp.resize_array( size_int );

      for ( int k=0; k<size_int; ++k )
        {
          tmp[ k ] = tab2_( k + k0 ) - 1;
        }

      tmp.ordonne_array( );

      for ( int k=0; k<size_int; ++k )
        {
          stencil( k+compteur , 0 ) = i;
          stencil( k+compteur , 1 ) =  tmp[ k ];
        }
      compteur += size;
    }


}

// Local template method : copy either value or ptr to value!
namespace
{
template<typename _T_> static inline void _fill_slot(_T_& dest, const double& src);
template<> inline void _fill_slot<double>(double& dest, const double& src)
{
  dest = src;
}
template<> inline void _fill_slot<const double *>(const double*& dest, const double& src)
{
  dest = &src;
}
}


template<typename _TAB_T_, typename _VALUE_T_>
inline void Matrice_Morse::get_stencil_coeff_templ( Stencil& stencil, _TAB_T_& coeffs_span) const
{
  auto nnz = tab2_.size_array();
  coeffs_span.resize(nnz);
  stencil.resize(nnz, 2);
  decltype(nnz) compteur = 0;
  const int nb_lines = nb_lignes( );
  for ( int i=0; i<nb_lines; ++i )
    {
      const auto k0      = tab1_( i ) - 1;
      const auto k1      = tab1_( i + 1 ) - 1;
      const int  size_int = (int)(k1 - k0);
      for ( int k=0; k<size_int; ++k )
        {
          stencil( compteur + k , 0 ) = i;
          stencil( compteur + k , 1 ) = tab2_( k + k0 ) - 1;
          ::_fill_slot<_VALUE_T_>(coeffs_span[ compteur + k ], coeff_(k+k0));
        }
      compteur += size_int;
    }


}


void Matrice_Morse::get_stencil_and_coeff_ptrs(Stencil& stencil,
                                               std::vector<const double *>& coeff_ptr) const
{
  assert_check_morse_matrix_structure( );

  if( is_stencil_up_to_date_ )
    {
      Cerr << "Error in Matrice_Morse::get_symmetric_stencil_and_coeff_ptrs( )"<<finl;
      Cerr << "  stencil up to date - function not impl. in this case."<<finl;
      Cerr << "  Aborting..." << finl;
      Process::abort( );
      return;
    }

  get_stencil_coeff_templ< std::vector<const double *>, const double *>(stencil, coeff_ptr);
  assert( (trustIdType)coeff_ptr.size( ) == stencil.dimension( 0 ));
}


void Matrice_Morse::get_stencil_and_coefficients( Stencil&       stencil,
                                                  StencilCoeffs& coefficients ) const
{
  if( is_stencil_up_to_date_ )
    {
      if( coeff_.size( ) == 0 )
        {
          Cerr << "Error in Matrice_Morse::get_stencil_and_coefficients( )"<<finl;
          Cerr << "  The coefficients are not filled."<<finl;
          Cerr << "  Aborting..." << finl;
          Process::abort( );
        }
      stencil = stencil_ ;
      { const auto sz = coeff_.size_array(); coefficients.resize(sz); for (auto k=sz-sz; k<sz; k++) coefficients[k] = coeff_[k]; }
      return;
    }



  get_stencil_coeff_templ<StencilCoeffs, double>(stencil, coefficients);
  assert( coefficients.size_array( ) == stencil.dimension( 0 ));
}


/*! @brief Operateur de division (de tous les elements) d'une matrice par un scalaire.
 *
 *     Operation: A =  A / x
 *
 * @param (double x) le parametre de scaling
 * @return (Matrice_Morse&) le resultat de l'operation (*this)
 * @throws division par zero impossible
 */
Matrice_Morse& Matrice_Morse::operator /=(double x )
{
  coeff_/=x;
  return(*this);
}

void Matrice_Morse::remplir(const IntLists& voisins,
                            const DoubleLists& valeurs,
                            const DoubleVect& terme_diag)
{

  int num_elem;
  int compteur,rang =0;

  // Remplissage des tableaux tab1, tab2 et coeff_ :
  auto* p_tab1 = tab1_.addr();
  int* p_tab2 = tab2_.addr();
  double* p_coeff = coeff_.addr();

  int* tab2_ptr = p_tab2;
  int n=nb_lignes();

  for (num_elem=0; num_elem<n; num_elem++)
    {

      IntList_Curseur liste_vois(voisins[num_elem]);
      DoubleList_Curseur liste_val(valeurs[num_elem]);
      compteur =0;
      *p_tab1++ = rang;

      *tab2_ptr++=num_elem;
      *p_coeff++ = terme_diag[num_elem];

      while  (liste_vois)
        {
          *tab2_ptr++ = liste_vois.valeur();
          *p_coeff++ = liste_val.valeur();
          ++liste_vois;
          ++liste_val;
          compteur++;
        }
      //       tab2[rang] = compteur;
      rang += (compteur + 1);
    }
  tab1_(num_elem)=rang;
  formeF();
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  is_stencil_up_to_date_=false;
}

void Matrice_Morse::remplir(const IntLists& voisins,
                            const DoubleLists& valeurs)
{

  int num_elem;
  int compteur,rang =0;

  // Remplissage des tableaux tab1, tab2 et coeff_ :
  auto* p_tab1 = tab1_.addr();
  int* p_tab2 = tab2_.addr();
  double* p_coeff = coeff_.addr();

  int* tab2_ptr = p_tab2;
  int n=nb_lignes();

  for (num_elem=0; num_elem<n; num_elem++)
    {

      IntList_Curseur liste_vois(voisins[num_elem]);
      DoubleList_Curseur liste_val(valeurs[num_elem]);
      compteur =0;
      *p_tab1++ = rang;

      while  (liste_vois)
        {
          *tab2_ptr++ = liste_vois.valeur();
          *p_coeff++ = liste_val.valeur();
          ++liste_vois;
          ++liste_val;
          compteur++;
        }
      //       tab2[rang] = compteur;
      rang += (compteur);
    }
  tab1_(num_elem)=rang;
  formeF();
  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  is_stencil_up_to_date_=false;
}

/*! @brief Remplissage d'une matrice morse par une matrice morse plus petite
 *
 */
void Matrice_Morse::remplir(const int ideb, const int jdeb, const int n, const int m, const Matrice_Morse& mat)
{
  // Verification
  assert(ideb<=n);
  assert(jdeb<=m);

  // On va construire une matrice locale
  Matrice_Morse matrice_locale(mat);
  // Cas ou la matrice locale est symetrique
  if (sub_type(Matrice_Morse_Sym,mat))
    {
      // Creation de la partie inferieure L
      Matrice_Morse L(matrice_locale);
      L.transpose(matrice_locale);
      int lordre = L.ordre();
      for (int i=0; i<lordre; i++)
        L(i, i) = 0.;
      // On ajoute M=U+L
      matrice_locale += L;
    }

  // Dimensionnement de la matrice globale
  auto nnz=matrice_locale.nb_coeff();
  dimensionner(n,m,(int)nnz);

  // Remplissage de la matrice globale par la matrice locale:
  // Remplissage de tab1_ avec le decalage par ideb:
  int mon_nb_lignes=matrice_locale.nb_lignes();
  assert(mon_nb_lignes+ideb<=n);
  for (int i=0; i<ideb; i++)
    tab1_(i)=1;
  for (int i=0; i<mon_nb_lignes; i++)
    tab1_(i+ideb)=matrice_locale.tab1_(i);
  for (int i=mon_nb_lignes+ideb; i<n+1; i++)
    tab1_(i)=matrice_locale.tab1_(mon_nb_lignes);

  // Remplissage de tab2_ avec le decalage par jdeb:
  for (auto i=0; i<nnz; i++)
    tab2_(i)=matrice_locale.tab2_(i)+jdeb;

  // Remplissage de coeff_:
  for (auto i=0; i<nnz; i++)
    coeff_(i)=matrice_locale.coeff_(i);

  morse_matrix_structure_has_changed_=1, sorted_ = 0;
  is_stencil_up_to_date_=false;
}

void Matrice_Morse::formeC()
{
  int n=nb_lignes();
  for(int ii=0; ii<=n; ii++)
    tab1_(ii)--;
  for(int ii=0; ii<n; ii++)
    tab2_(tab1_(ii))=nb_vois(ii);
  for(auto k=0; k<nb_coeff(); k++)
    tab2_(k)--;
  morse_matrix_structure_has_changed_=1;
  is_stencil_up_to_date_=false;
}

void Matrice_Morse::formeF()
{
  int n=nb_lignes();
  for(int ii=0; ii<=n; ii++)
    tab1_(ii)++;
  for(auto k=0; k<nb_coeff(); k++)
    tab2_(k)++;
  morse_matrix_structure_has_changed_=1;
  is_stencil_up_to_date_=false;
}




/*! @brief NE FAIT RIEN
 *
 * @return (int) renvoie toujours 1
 */
int Matrice_Morse_test()
{
  return 1;
}

/*! @brief Remplit la matrice avec des zeros.
 *
 */

void Matrice_Morse::clean()
{
  coeff_ = 0;
}

/*! @brief Calcule la largeur de bande d'une matrice morse
 *
 */
int Matrice_Morse::largeur_de_bande() const
{
  int ldist,min = 0;
  const auto* p_tab1_ = get_tab1().addr();
  const int* p_tab2_ = get_tab2().addr();
  int N=ordre();

  for(int i=0; i<N; i++)
    for(auto k = p_tab1_[i]; k < p_tab1_[i+1]; k++)
      {
        if (p_tab2_[k-1]-1<N)
          {
            ldist = p_tab2_[k-1] - i;
            if( min < ldist ) min = ldist;
          }
      };
  return min;
}

bool Matrice_Morse::check_morse_matrix_structure() const
{
  const int nb_lines   = nb_lignes( );
  const int nb_columns = nb_colonnes( );
  const auto nb_coefficients = tab1_( nb_lines ) - 1;

  if ( tab2_.size_array( ) != nb_coefficients )
    {
      Cerr << "invalid tab2 size" << finl;
      return false;
    }

  if ( coeff_.size_array( ) != nb_coefficients )
    {
      Cerr << "invalid coeff size" << finl;
      return false;
    }

  ArrOfBit flags( nb_columns );

  for ( int i=0; i<nb_lines; ++i )
    {
      flags = 0;

      auto k0 = tab1_( i ) - 1;
      auto k1 = tab1_( i + 1 ) - 1;

      for ( auto k=k0; k<k1; ++k )
        {
          int j = tab2_( k ) - 1;

          if ( j < 0 )
            {
              Cerr << "invalid column index (<0): " << j << finl;
              return false;
            }

          if ( j >= nb_columns )
            {
              Cerr << "invalid column index (>nb_cols): " << j << " > " << nb_columns << finl;
              return false;
            }

          if ( flags[ j ] )
            {
              Cerr << "invalid coefficient ( " << i << ", " << j << " ): already defined ( " << k << " )" << finl;
              return false;
            }

          flags.setbit( j );
        }
    }

  return true;
}

bool Matrice_Morse::check_sorted_morse_matrix_structure() const
{
  const int nb_lines   = nb_lignes( );
  const int nb_columns = nb_colonnes( );
  const auto nb_coefficients = tab1_( nb_lines ) - 1;

  if ( tab2_.size_array( ) != nb_coefficients )
    {
      Cerr << "invalid tab2 size" << finl;
      return false;
    }

  if ( coeff_.size_array( ) != nb_coefficients )
    {
      Cerr << "invalid coeff size" << finl;
      return false;
    }

  ArrOfBit flags( nb_columns );

  for ( int i=0; i<nb_lines; ++i )
    {
      flags = 0;

      auto k0 = tab1_( i ) - 1;
      auto k1 = tab1_( i + 1 ) - 1;

      int j0 = tab2_( k0 ) - 1 - 1;

      for ( auto k=k0; k<k1; ++k )
        {
          int j = tab2_( k ) - 1;

          if ( j < 0 )
            {
              Cerr << "invalid column index (<0): " << j << finl;
              return false;
            }

          if ( j >= nb_columns )
            {
              Cerr << "invalid column index (>nb_cols): " << j << " > " << nb_columns << finl;
              return false;
            }

          if ( flags[ j ] )
            {
              Cerr << "invalid coefficient ( " << i << ", " << j << " ): already defined ( " << k << " )" << finl;
              return false;
            }

          if ( j <= j0 )
            {
              Cerr << "unsorted coefficient: ( " << i << ", " << j << " ) after ( " << i << ", " << j0 << " ) " << finl;;
              return false;
            }

          j0 = j;
          flags.setbit( j );
        }
    }

  return true;
}


void Matrice_Morse::assert_check_morse_matrix_structure() const
{
  if (!morse_matrix_structure_has_changed_) return;
#ifndef NDEBUG
  if ( ! ( check_morse_matrix_structure( ) ) )
    {
      Cerr << "Error in 'Matrice_Morse::assert_check_morse_matrix_structure( )':" << finl;
      Cerr << "  Exiting..." << finl;
      Process::exit( );
    }
  else
    morse_matrix_structure_has_changed_=0;
#endif
}

void Matrice_Morse::assert_check_sorted_morse_matrix_structure() const
{
  if (!morse_matrix_structure_has_changed_) return;
#ifndef NDEBUG
  if ( ! ( check_sorted_morse_matrix_structure( ) ) )
    {
      Cerr << "Error in 'Matrice_Morse::assert_check_sorted_morse_matrix_structure( )':" << finl;
      Cerr << "  Exiting..." << finl;
      Process::exit( );
    }
  else
    morse_matrix_structure_has_changed_=0;
#endif
}

// Build a new Morse matrix spanning the rectangular area defined by the two points (nl0, nc0) and (nl1, nc1)
// in the original matrix.
// Indices are provided in C mode (0-based indexing).
void Matrice_Morse::construire_sous_bloc(int nl0, int nc0, int nl1, int nc1, Matrice_Morse& result) const
{
  // count non-zero entries:
  assert(nl0 >= 0);
  assert(nc0 >= 0);
  assert(nl0 <= nl1);
  assert(nc0 <= nc1);

  auto max_nnz = tab1_(nl1+1) - tab1_(nl0); // maximum number of zeros that we will find
  int tot=0;
  IntTab loca((int)max_nnz, 2);
  DoubleTab sub_coeffs((int)max_nnz);
  for (int li=nl0; li <= nl1; li++)
    {
      auto idx_coeff = tab1_(li)-1;
      int nb_coeff_on_line = (int)(tab1_(li+1)-tab1_(li));
      for (int j=0; j < nb_coeff_on_line; j++)
        {
          int col_idx = tab2_(j+idx_coeff)-1;
          if (col_idx >= nc0 && col_idx <= nc1) // is the coeff in the window?
            {
              loca(tot, 0) = li - nl0;
              loca(tot, 1) = col_idx - nc0;
              sub_coeffs(tot) = coeff_(j+idx_coeff);
              tot++;
            }
        }
    }
  loca.resize(tot,2);
  sub_coeffs.resize(tot);

  result.dimensionner(loca);
  // Set coefficient values:
  for (int i =0 ; i < tot; i++)
    {
      int il = loca(i, 0);
      int ic = loca(i, 1);
      result.coef(il, ic) = sub_coeffs(i);
    }
}

void Matrice_Morse::sort_stencil()
{
  if (sorted_) return; //deja fait
  for (int i = 0; i + 1 < tab1_.size_array(); i++) //indice de ligne
    std::sort(tab2_.addr() + tab1_(i) - 1, tab2_.addr() + tab1_(i + 1) - 1);
  morse_matrix_structure_has_changed_ = sorted_ = 1;
}

// Check if the matrix is sorted based on a stencil condition
bool Matrice_Morse::is_sorted_stencil() const
{
  if (!sorted_)
    {
      const int n = nb_lignes();
      for (int i = 0; i < n; i++)
        {
          const auto k0 = tab1_( i ) - 1;
          const auto k1 = tab1_( i + 1 ) - 1;
          for (auto k=k0; k<k1-1; k++)
            if (tab2_(k)>tab2_(k+1))
              return sorted_; // not sorted
        }
      sorted_ = true;
    }
  return sorted_;
}

// Check the matrix is diagonal:
// Faster than using:
// Stencil stencil;
// A.get_stencil(stencil);
// Matrix_tools::is_diagonal_stencil(A.nb_lignes(), A.nb_colonnes(), stencil);
bool Matrice_Morse::is_diagonal()
{
  bool is_diagonal = true;
  const int n = nb_lignes();
  for (int i = 0; i < n; i++)
    {
      const auto k1 = get_tab1()(i) - 1;
      const auto k2 = get_tab1()(i + 1) - 1;
      for (auto k = k1; k < k2; k++)
        {
          if (k2-k1>1 || get_tab2()(k)-1!=i)
            {
              is_diagonal = false;
              break;
            }
        }
    }
  return is_diagonal;
}

// Explicit instantiations for 'auto nnz' abbreviated function templates
template Matrice_Morse::Matrice_Morse(int, int);
template Matrice_Morse::Matrice_Morse(int, int, int);
template void Matrice_Morse::dimensionner(int, int);
template void Matrice_Morse::dimensionner(int, int, int);
#ifdef TRUST_USE_GPU
template Matrice_Morse::Matrice_Morse(int, trustIdType);
template Matrice_Morse::Matrice_Morse(int, int, trustIdType);
template void Matrice_Morse::dimensionner(int, trustIdType);
template void Matrice_Morse::dimensionner(int, int, trustIdType);
#endif
