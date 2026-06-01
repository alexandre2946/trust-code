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
#include <Matrice_Bloc_Sym.h>
#include <MD_Vector_tools.h>
#include <MD_Vector_base.h>
#include <TRUSTTab_parts.h>
#include <Motcle.h>
#include <Param.h>
#include <SSOR.h>

Implemente_instanciable_sans_constructeur(SSOR,"SSOR",Precond_base);
// XD ssor precond_base ssor BRACE Symmetric successive over-relaxation algorithm.

SSOR::SSOR() : omega_(1.6), algo_fortran_(-1), avec_assert_(-1), algo_items_communs_(-1), line_size_(0) { }

Sortie& SSOR::printOn(Sortie& s ) const
{
  s << " { omega  "<<omega_ << " } ";
  return s;
}

Entree& SSOR::readOn(Entree& is )
{
  Param param(que_suis_je());
  param.ajouter("omega", &omega_);  // XD attr omega floattant omega OPT Over-relaxation facteur (between 1 and 2,
  // XD_CONT default value 1.6).
  param.lire_avec_accolades(is);

  if (omega_ <= 0. || omega_ >= 2.)
    {
      Cerr << "SSOR::readOn, omega is not within [0, 2]: SSOR not activated (you should use precond nul instead) "  << finl;
    }
  return is;
}

void SSOR::prepare_(const Matrice_Base& la_matrice, const DoubleVect& secmem)
{
  if (get_status() >= REINIT_ALL)
    {
      md_secmem_ = secmem.get_md_vector();
      line_size_ = secmem.line_size();
      // For the next preconditioning step, check the matrix
      avec_assert_ = 1;

      if (nproc() == 1 || !(md_secmem_)) algo_items_communs_ = 0;
      else
        {
          // Number of sequential items on this proc
          const int sz_tot = secmem.size_totale();
          items_a_traiter_.reset();
          items_a_traiter_.resize(sz_tot / line_size_, line_size_, RESIZE_OPTIONS::NOCOPY_NOINIT);
          items_a_traiter_.set_md_vector(md_secmem_);
          int n = md_secmem_->get_sequential_items_flags(items_a_traiter_, line_size_);
          int sz = md_secmem_->get_nb_items_reels();

          if (sz < 0) // size() is invalid: real items are not grouped at the start!
            algo_items_communs_ = 1;
          else
            {
              assert(sz >= n);
              if (mp_sum(sz) > mp_sum(n)) algo_items_communs_ = 1; // There are shared items among the real items
              else
                {
                  algo_items_communs_ = 0;
                  items_a_traiter_.reset();
                }
            }
        }
    }

  if (get_status() >= REINIT_COEFF) { /* Do nothing */ }
  Precond_base::prepare_(la_matrice, secmem);
}

/*! @brief Calcule la solution du systeme lineaire: A * solution = b avec la methode de relaxation SSOR.
 *
 */
int SSOR::preconditionner_(const Matrice_Base& la_matrice, const DoubleVect& b, DoubleVect& solution)
{
  // for historical compatibility:
  if (omega_ <= 0. || omega_ >= 2)
    {
      operator_egal(solution, b);
      return 1;
    }

  operator_egal(solution, b);

  if (sub_type(Matrice_Morse_Sym, la_matrice))
    {
      const Matrice_Morse_Sym& matrice = ref_cast(Matrice_Morse_Sym, la_matrice);
      ssor(matrice, solution);
    }
  else if (sub_type(Matrice_Bloc_Sym, la_matrice))
    {
      // Matrix corresponding to a multi-location vector (MD_Vector_composite)
      const Matrice_Bloc_Sym& matrice = ref_cast(Matrice_Bloc_Sym, la_matrice);
      ssor(matrice, solution);
    }
  else if (sub_type(Matrice_Bloc, la_matrice))
    {
      // Assume a real-real matrix and a real-virtual matrix
      const Matrice_Bloc& mat = ref_cast(Matrice_Bloc, la_matrice);
      const Matrice_Morse_Sym& matrice = ref_cast(Matrice_Morse_Sym, mat.get_bloc(0, 0).valeur());
      ssor(matrice, solution);
    }
  else
    {
      Cerr << "SSOR::preconditionner not coded for type " << la_matrice.que_suis_je() << finl;
      exit();
    }
  if (echange_ev_solution_) solution.echange_espace_virtuel();
  return 1;
}

void traite_diagonale(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  const int nb_lignes_a_traiter = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  const auto& tab1 = mat.get_tab1();
  const auto& coeff = mat.get_coeff();
  const double psi = (2. - omega) / omega;
  const double *coeff_fortran = coeff.addr() - 1; // indexable by Fortran index
  const auto *tab1_ptr = tab1.addr();
  double *vect_ptr = vecteur.addr();
  for (int i = nb_lignes_a_traiter; i; i--, tab1_ptr++, vect_ptr++)
    {
      const auto j = *tab1_ptr;
      // Diagonal coefficient of row i:
      const double coeff_i_i = coeff_fortran[j];
      (*vect_ptr) *= psi * coeff_i_i;
    }
}

enum class descente_enum { NORMAL , NORMAL_ASSERT , DIAG , DIAG_ASSERT };

template<descente_enum _TYPE_>
void descente_generique(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  static constexpr bool IS_NORMAL_ASSERT = (_TYPE_ == descente_enum::NORMAL_ASSERT), IS_DIAG_ASSERT = (_TYPE_ == descente_enum::DIAG_ASSERT), NOT_DIAG = (_TYPE_ != descente_enum::DIAG && _TYPE_ != descente_enum::DIAG_ASSERT);
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();
  const int nb_lignes_a_traiter = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  // "Fortran" pointer to the solution array (indexable with Fortran index); the pointer is constant, not the pointed data.
  double * const sol_fortran = vecteur.addr() - 1;
  const auto *tab1_ptr = tab1.addr();
  assert(nb_lignes_a_traiter <= tab1.size_array() + 1);
  assert(*tab1_ptr == 1); // otherwise the 2 lines below are wrong.
  const int *tab2_ptr = tab2.addr();
  const double *coeff_ptr = coeff.addr();
  auto last_tab1_de_i = *tab1_ptr;
  tab1_ptr++;
  for (int i = 1; i <= nb_lignes_a_traiter; i++, tab1_ptr++)
    {
      const auto tab1_de_i = *tab1_ptr; // = tab1[i]
      const int nvois = (int)(tab1_de_i - last_tab1_de_i);

      if (IS_NORMAL_ASSERT || IS_DIAG_ASSERT) assert(nvois >= 0 && (tab1_de_i - 1) <= tab2.size_array());

      last_tab1_de_i = tab1_de_i;
      // This check must remain: without it, other rows of the solution vector would be polluted by values depending on shared and virtual items
      {
        double v_i;
        if (NOT_DIAG)
          {
            // The first coefficient must be the diagonal coefficient and must be strictly positive
            if (IS_NORMAL_ASSERT)  assert(nvois >= 1 && (*tab2_ptr) == i && (*coeff_ptr) > 0.);

            const double omega_coeff_i_i = omega / (*coeff_ptr);
            v_i = sol_fortran[i] *= omega_coeff_i_i;
            tab2_ptr++;
            coeff_ptr++;
          }
        else // PRECOND DIAG !!
          {
            // No stored diagonal coefficient (neither in tab2 nor in coeff), the diagonal equals 1:
            v_i = sol_fortran[i] *= omega;
          }

        for (int j = nvois-1; j; j--, tab2_ptr++, coeff_ptr++)
          {
            const int i2 = *tab2_ptr; // column index for the next coefficient
            const double coeff_i_i2 = *coeff_ptr;
            // the matrix only has upper diagonal coefficients and the diagonal has already been handled, so i2 > i. No virtual items allowed!
            if (IS_NORMAL_ASSERT || IS_DIAG_ASSERT) assert(i2 > i && i2 <= nb_lignes_a_traiter);
            // B.M.: ... however, the check on shared items is redundant here because the value will be zeroed out (see **Annulation items communs**);
            // skipping the check gives a performance gain
            sol_fortran[i2] -= coeff_i_i2 * v_i;
          }
      }
    }
}

void descente(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  descente_generique<descente_enum::NORMAL>(omega,mat,vecteur);
}

void descente_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  descente_generique<descente_enum::NORMAL_ASSERT>(omega,mat,vecteur);
}

void descente_diag_ok_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  descente_generique<descente_enum::NORMAL_ASSERT>(omega,mat,vecteur); /* meme qu'avant :D */
}

void descente_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  descente_generique<descente_enum::DIAG>(omega,mat,vecteur);
}

void descente_assert_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  descente_generique<descente_enum::DIAG_ASSERT>(omega,mat,vecteur);
}

// Forward sweep on an off-diagonal block. Method called by SSOR(const Matrice_bloc & ...)
// vecteur: size = "number of rows of the matrix", vecteur2: size = "number of columns"
void descente_bloc_extradiag_assert(const Matrice_Morse& mat, const DoubleVect& vecteur, DoubleVect& vecteur2)
{
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();
  const int nb_lignes = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  const double *vecteur_ptr = vecteur.addr();
  double *vecteur2_fortran_ptr = vecteur2.addr() - 1;
  const int *tab2_fortran_ptr = tab2.addr() - 1;
  const double *coeff_fortran_ptr = coeff.addr() - 1;
  assert(coeff.size_array() == tab2.size_array());
  for (int i_ligne = 0; i_ligne < nb_lignes; i_ligne++, vecteur_ptr++)
    {
      {
        const double v_i = *vecteur_ptr;
        auto index = tab1[i_ligne];
        const auto index_fin = tab1[i_ligne + 1];
        // There may be no coefficient on the row => index_fin == index
        assert(index > 0 && index_fin >= index && index_fin <= tab2.size_array() + 1);
        const int *tab2_ptr = tab2_fortran_ptr + index;
        const double *coeff_ptr = coeff_fortran_ptr + index;
        for (; index < index_fin; index++, tab2_ptr++, coeff_ptr++)
          {
            const int i_colonne = *tab2_ptr;
            assert(i_colonne >= 1 && i_colonne <= vecteur2.size_array());
            const double c = *coeff_ptr;
            // no shared-item check on columns, see **Annulation items communs**
            vecteur2_fortran_ptr[i_colonne] -= c * v_i;
          }
      }
    }
}

template<descente_enum _TYPE_>
void descente_generique(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  static constexpr bool IS_NORMAL_ASSERT = (_TYPE_ == descente_enum::NORMAL_ASSERT), IS_DIAG_ASSERT = (_TYPE_ == descente_enum::DIAG_ASSERT), NOT_DIAG = (_TYPE_ != descente_enum::DIAG && _TYPE_ != descente_enum::DIAG_ASSERT);
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();
  const int nb_lignes_a_traiter = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  // "Fortran" pointer to the solution array (indexable with Fortran index); the pointer is constant, not the pointed data.
  double * const sol_fortran = vecteur.addr() - 1;
  const int *flags_ptr = items_a_traiter.addr();
  assert(nb_lignes_a_traiter <= items_a_traiter.size_array());
  const auto *tab1_ptr = tab1.addr();
  assert(nb_lignes_a_traiter <= tab1.size_array() + 1);
  assert(*tab1_ptr == 1); // otherwise the 2 lines below are wrong.
  const int *tab2_ptr = tab2.addr();
  const double *coeff_ptr = coeff.addr();
  auto last_tab1_de_i = *tab1_ptr;
  tab1_ptr++;
  for (int i = 1; i <= nb_lignes_a_traiter; i++, tab1_ptr++)
    {
      const auto tab1_de_i = *tab1_ptr; // = tab1[i]
      const int nvois = (int)(tab1_de_i - last_tab1_de_i);

      if (IS_NORMAL_ASSERT || IS_DIAG_ASSERT) assert(nvois >= 0 && (tab1_de_i - 1) <= tab2.size_array());

      last_tab1_de_i = tab1_de_i;
      // This check must remain: without it, other rows of the solution vector would be polluted by values depending on shared and virtual items
      if (IS_NORMAL_ASSERT || IS_DIAG_ASSERT) assert((flags_ptr - items_a_traiter.addr()) == (i-1));

      const int item_a_traiter = *(flags_ptr++);
      if (item_a_traiter)
        {
          double v_i;
          if (NOT_DIAG)
            {
              // The first coefficient must be the diagonal coefficient and must be strictly positive
              if (IS_NORMAL_ASSERT) assert(nvois >= 1 && (*tab2_ptr) == i && (*coeff_ptr) > 0.);

              const double omega_coeff_i_i = omega / (*coeff_ptr);
              v_i = sol_fortran[i] *= omega_coeff_i_i;
              tab2_ptr++;
              coeff_ptr++;
            }
          else // PRECOND DIAG !!
            {
              // No stored diagonal coefficient (neither in tab2 nor in coeff), the diagonal equals 1:
              v_i = sol_fortran[i] *= omega;
            }

          for (int j = nvois-1; j; j--, tab2_ptr++, coeff_ptr++)
            {
              const int i2 = *tab2_ptr; // column index for the next coefficient
              const double coeff_i_i2 = *coeff_ptr;
              // the matrix only has upper diagonal coefficients and the diagonal has already been handled, so i2 > i.
              if (IS_NORMAL_ASSERT || IS_DIAG_ASSERT) assert(i2 > i && i2 <= vecteur.size_totale());

              // B.M.: ... however, the check on shared items is redundant here because
              // the value will be zeroed out (see **Annulation items communs**); skipping the check gives a performance gain
              sol_fortran[i2] -= coeff_i_i2 * v_i;
            }
        }
      else
        {
          // **Annulation items communs**
          // this is the last write to sol_fortran[i] because the matrix is upper triangular. To avoid the items_a_traiter_ check
          // in the back-substitution, the solution is zeroed for shared and virtual items (sol_fortran[i] will not be modified again in the forward sweep)
          sol_fortran[i] = 0.;
          coeff_ptr += nvois;
          tab2_ptr += nvois;
        }
    }
  // **Annulation items communs** : For the back-substitution, values in virtual slots must be zeroed
  {
    const int fin = vecteur.size_totale();
    for (int i = nb_lignes_a_traiter+1; i <= fin; i++)
      sol_fortran[i] = 0.;
  }
}

void descente(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  descente_generique<descente_enum::NORMAL>(omega,mat,vecteur,items_a_traiter);
}

void descente_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  descente_generique<descente_enum::NORMAL_ASSERT>(omega,mat,vecteur,items_a_traiter);
}

void descente_diag_ok_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  descente_generique<descente_enum::NORMAL_ASSERT>(omega,mat,vecteur,items_a_traiter); /* meme qu'avant :D */
}

void descente_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  descente_generique<descente_enum::DIAG>(omega,mat,vecteur,items_a_traiter);
}

void descente_assert_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  descente_generique<descente_enum::DIAG_ASSERT>(omega,mat,vecteur,items_a_traiter);
}

// Forward sweep on an off-diagonal block. Method called by SSOR(const Matrice_bloc & ...)
// vecteur: size = "number of rows of the matrix", vecteur2: size = "number of columns"
void descente_bloc_extradiag_assert(const Matrice_Morse& mat, const DoubleVect& vecteur, DoubleVect& vecteur2, const ArrOfInt& items_a_traiter)
{
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();
  const int nb_lignes = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  const int *flags_ptr = items_a_traiter.addr();
  const double *vecteur_ptr = vecteur.addr();
  double *vecteur2_fortran_ptr = vecteur2.addr() - 1;
  const int *tab2_fortran_ptr = tab2.addr() - 1;
  const double *coeff_fortran_ptr = coeff.addr() - 1;
  assert(coeff.size_array() == tab2.size_array());
  for (int i_ligne = 0; i_ligne < nb_lignes; i_ligne++, vecteur_ptr++)
    {
      if (*(flags_ptr++))
        {
          const double v_i = *vecteur_ptr;
          auto index = tab1[i_ligne];
          const auto index_fin = tab1[i_ligne + 1];
          // There may be no coefficient on the row => index_fin == index
          assert(index > 0 && index_fin >= index && index_fin <= tab2.size_array() + 1);
          const int *tab2_ptr = tab2_fortran_ptr + index;
          const double *coeff_ptr = coeff_fortran_ptr + index;
          for (; index < index_fin; index++, tab2_ptr++, coeff_ptr++)
            {
              const int i_colonne = *tab2_ptr;
              assert(i_colonne >= 1 && i_colonne <= vecteur2.size_array());
              const double c = *coeff_ptr;
              // no shared-item check on columns, see **Annulation items communs**
              vecteur2_fortran_ptr[i_colonne] -= c * v_i;
            }
        }
    }
}

enum class remontee_enum { NORMAL , NORMAL_ASSERT , DIAG_OK_ASSERT , DIAG , DIAG_ASSERT };

template<remontee_enum _TYPE_>
void remontee_generique(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  static constexpr bool IS_NORMAL_ASSERT = (_TYPE_ == remontee_enum::NORMAL_ASSERT), IS_DIAG_OK_ASSERT = (_TYPE_ == remontee_enum::DIAG_OK_ASSERT), IS_DIAG_ASSERT = (_TYPE_ == remontee_enum::DIAG_ASSERT),
                        NOT_DIAG = (_TYPE_ != remontee_enum::DIAG && _TYPE_ != remontee_enum::DIAG_ASSERT);
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();

  const int nb_lignes_a_traiter = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  const double psi = IS_DIAG_OK_ASSERT ? -1e10 : (2. - omega) / omega;

  // "fortran" pointer to the solution array (indexable with Fortran index)
  const double *const sol_fortran = vecteur.addr() - 1;
  const auto *tab1_ptr = tab1.addr() + nb_lignes_a_traiter;
  auto last_tab1_de_i = *tab1_ptr;
  tab1_ptr--;
  // We do not go to the end of tab2 because nb_lignes_a_traiter = tab1.size_array() is not guaranteed:
  // -2 because last_tab1_de_i is the Fortran index of the first coefficient of the next row
  //  (-1 for Fortran->C conversion and -1 to move to the last coefficient of the previous row)
  const int *tab2_ptr = tab2.addr() + last_tab1_de_i - 2;
  const double *coeff_ptr = coeff.addr() + last_tab1_de_i - 2;
  double *soli_ptr = vecteur.addr() + nb_lignes_a_traiter - 1;
  for (int i = nb_lignes_a_traiter; i; i--, tab1_ptr--, soli_ptr--)
    {
      const auto tab1_de_i = *tab1_ptr; // = tab1[i]
      const int nvois = (int)(last_tab1_de_i - tab1_de_i);

      if (IS_NORMAL_ASSERT || IS_DIAG_OK_ASSERT || IS_DIAG_ASSERT) assert(nvois >= 0 && tab1_de_i > 0);

      last_tab1_de_i = tab1_de_i;
      // This check must remain: sol[i] must not be modified for shared
      // and virtual items (they are zero and must stay zero to avoid polluting other rows):
      if (1)
        {
          // Operation "diagonale":
          double x = 0.;
          for (int j = nvois - 1; j; j--, tab2_ptr--, coeff_ptr--)
            {
              const int i2 = *tab2_ptr;
              const double coeff_i_i2 = *coeff_ptr;
              // Zero operation for shared and virtual items
              // (their RHS contribution was zeroed during the forward sweep)
              x += coeff_i_i2 * sol_fortran[i2];
            }

          if (NOT_DIAG)
            {
              // here coeff_ptr is the diagonal coefficient
              if (IS_NORMAL_ASSERT || IS_DIAG_OK_ASSERT) assert((*tab2_ptr) == i);

              const double coeff_i_i = *coeff_ptr;
              const double omega_coeff_i_i = omega / coeff_i_i;
              // If IS_DIAG_OK_ASSERT: the diagonal has already been multiplied by psi * coeff_i_i
              *soli_ptr = IS_DIAG_OK_ASSERT ? ((*soli_ptr) - x) * omega_coeff_i_i : (*soli_ptr) * psi * omega - x * omega_coeff_i_i;
              coeff_ptr--;
              tab2_ptr--;
            }
          else // PRECOND DIAG !!
            {
              // Diagonal preconditioning, no diagonal coefficient stored:
              *soli_ptr = ((*soli_ptr) * psi - x) * omega;
            }
        }
      else
        {
          assert((*soli_ptr) == 0.);
          coeff_ptr -= nvois;
          tab2_ptr -= nvois;
        }
    }
}

void remontee(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  remontee_generique<remontee_enum::NORMAL>(omega,mat,vecteur);
}

void remontee_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  remontee_generique<remontee_enum::NORMAL_ASSERT>(omega,mat,vecteur);
}

void remontee_diag_ok_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  remontee_generique<remontee_enum::DIAG_OK_ASSERT>(omega,mat,vecteur);
}

void remontee_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  remontee_generique<remontee_enum::DIAG>(omega,mat,vecteur);
}

void remontee_assert_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur)
{
  remontee_generique<remontee_enum::DIAG_ASSERT>(omega,mat,vecteur);
}

// Back-substitution on an off-diagonal block. Method called by SSOR(const Matrice_bloc & ...)
// vecteur: size = "number of rows of the matrix", vecteur2: size = "number of columns"
void remontee_bloc_extradiag_assert(const Matrice_Morse& mat, DoubleVect& vecteur, const DoubleVect& vecteur2)
{
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();
  const int nb_lignes = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  double *vecteur_ptr = vecteur.addr();
  const double *vecteur2_fortran_ptr = vecteur2.addr() - 1;
  const int *tab2_fortran_ptr = tab2.addr() - 1;
  const double *coeff_fortran_ptr = coeff.addr() - 1;
  assert(coeff.size_array() == tab2.size_array());
  for (int i_ligne = 0; i_ligne < nb_lignes; i_ligne++, vecteur_ptr++)
    {
      {
        double x = *vecteur_ptr;
        auto index = tab1[i_ligne];
        const auto index_fin = tab1[i_ligne + 1];
        // There may be no coefficient on the row => index_fin == index
        assert(index > 0 && index_fin >= index && index_fin <= tab2.size_array() + 1);
        const int *tab2_ptr = tab2_fortran_ptr + index;
        const double *coeff_ptr = coeff_fortran_ptr + index;
        for (; index < index_fin; index++, tab2_ptr++, coeff_ptr++)
          {
            const int i_colonne = *tab2_ptr;
            assert(i_colonne >= 1 && i_colonne <= vecteur2.size_array());
            const double c = *coeff_ptr;
            const double x2 = vecteur2_fortran_ptr[i_colonne];
            // no shared-item check on columns, see **Annulation items communs**
            x -= c * x2;
          }
        *vecteur_ptr = x;
      }
    }
}

template<remontee_enum _TYPE_>
void remontee_generique(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  static constexpr bool IS_NORMAL_ASSERT = (_TYPE_ == remontee_enum::NORMAL_ASSERT), IS_DIAG_OK_ASSERT = (_TYPE_ == remontee_enum::DIAG_OK_ASSERT), IS_DIAG_ASSERT = (_TYPE_ == remontee_enum::DIAG_ASSERT),
                        NOT_DIAG = (_TYPE_ != remontee_enum::DIAG && _TYPE_ != remontee_enum::DIAG_ASSERT);
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();

  const int nb_lignes_a_traiter = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  const int *flags_ptr = items_a_traiter.addr() + nb_lignes_a_traiter - 1;
  assert(nb_lignes_a_traiter <= items_a_traiter.size_array());
  const double psi =  IS_DIAG_OK_ASSERT ? -1e10 : (2. - omega) / omega;

  // "fortran" pointer to the solution array (indexable with Fortran index)
  const double * const sol_fortran = vecteur.addr() - 1;
  const auto *tab1_ptr = tab1.addr() + nb_lignes_a_traiter;
  auto last_tab1_de_i = *tab1_ptr;
  tab1_ptr--;
  // We do not go to the end of tab2 because nb_lignes_a_traiter = tab1.size_array() is not guaranteed:
  // -2 because last_tab1_de_i is the Fortran index of the first coefficient of the next row
  //  (-1 for Fortran->C conversion and -1 to move to the last coefficient of the previous row)
  const int *tab2_ptr = tab2.addr() + last_tab1_de_i - 2;
  const double *coeff_ptr = coeff.addr() + last_tab1_de_i - 2;
  double * soli_ptr = vecteur.addr() + nb_lignes_a_traiter - 1;
  for (int i = nb_lignes_a_traiter; i; i--, tab1_ptr--, soli_ptr--)
    {
      const auto tab1_de_i = *tab1_ptr; // = tab1[i]
      const int nvois = (int)(last_tab1_de_i - tab1_de_i);
      if (IS_NORMAL_ASSERT || IS_DIAG_OK_ASSERT || IS_DIAG_ASSERT) assert(nvois >= 0 && tab1_de_i > 0);

      last_tab1_de_i = tab1_de_i;
      // This check must remain: sol[i] must not be modified for shared
      // and virtual items (they are zero and must stay zero to avoid polluting other rows):
      if (*(flags_ptr--))
        {
          // Operation "diagonale":
          double x = 0.;
          for (int j = nvois-1; j; j--, tab2_ptr--, coeff_ptr--)
            {
              const int i2 = *tab2_ptr;
              const double coeff_i_i2 = *coeff_ptr;
              // Zero operation for shared and virtual items
              // (their RHS contribution was zeroed during the forward sweep)
              x += coeff_i_i2 * sol_fortran[i2];
            }
          if (NOT_DIAG)
            {
              // here coeff_ptr is the diagonal coefficient
              if (IS_NORMAL_ASSERT || IS_DIAG_OK_ASSERT) assert((*tab2_ptr) == i);

              const double coeff_i_i = *coeff_ptr;
              const double omega_coeff_i_i = omega / coeff_i_i;
              // If IS_DIAG_OK_ASSERT: the diagonal has already been multiplied by psi * coeff_i_i
              *soli_ptr = IS_DIAG_OK_ASSERT ? ((*soli_ptr) - x) * omega_coeff_i_i : (*soli_ptr) * psi * omega - x * omega_coeff_i_i;
              coeff_ptr--;
              tab2_ptr--;
            }
          else // PRECOND DIAG !!
            {
              // Diagonal preconditioning, no diagonal coefficient stored:
              *soli_ptr = ((*soli_ptr) * psi - x) * omega;
            }
        }
      else
        {
          assert((*soli_ptr) == 0.);
          coeff_ptr -= nvois;
          tab2_ptr -= nvois;
        }
    }
}

void remontee(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  remontee_generique<remontee_enum::NORMAL>(omega,mat,vecteur,items_a_traiter);
}

void remontee_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  remontee_generique<remontee_enum::NORMAL_ASSERT>(omega,mat,vecteur,items_a_traiter);
}

void remontee_diag_ok_assert(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  remontee_generique<remontee_enum::DIAG_OK_ASSERT>(omega,mat,vecteur,items_a_traiter);
}

void remontee_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  remontee_generique<remontee_enum::DIAG>(omega,mat,vecteur,items_a_traiter);
}

void remontee_assert_precond_diag(const double omega, const Matrice_Morse& mat, DoubleVect& vecteur, const ArrOfInt& items_a_traiter)
{
  remontee_generique<remontee_enum::DIAG_ASSERT>(omega,mat,vecteur,items_a_traiter);
}

// Back-substitution on an off-diagonal block. Method called by SSOR(const Matrice_bloc & ...)
// vecteur: size = "number of rows of the matrix", vecteur2: size = "number of columns"
void remontee_bloc_extradiag_assert(const Matrice_Morse& mat, DoubleVect& vecteur, const DoubleVect& vecteur2, const ArrOfInt& items_a_traiter)
{
  const auto& tab1 = mat.get_tab1();
  const auto& tab2 = mat.get_tab2();
  const auto& coeff = mat.get_coeff();
  const int nb_lignes = vecteur.size_reelle_ok() ? vecteur.size_reelle() : vecteur.size_totale();
  const int *flags_ptr = items_a_traiter.addr();
  double *vecteur_ptr = vecteur.addr();
  const double *vecteur2_fortran_ptr = vecteur2.addr() - 1;
  const int *tab2_fortran_ptr = tab2.addr() - 1;
  const double *coeff_fortran_ptr = coeff.addr() - 1;
  assert(coeff.size_array() == tab2.size_array());
  for (int i_ligne = 0; i_ligne < nb_lignes; i_ligne++, vecteur_ptr++)
    {
      if (*(flags_ptr++))
        {
          double x = *vecteur_ptr;
          auto index = tab1[i_ligne];
          const auto index_fin = tab1[i_ligne + 1];
          // There may be no coefficient on the row => index_fin == index
          assert(index > 0 && index_fin >= index && index_fin <= tab2.size_array() + 1);
          const int *tab2_ptr = tab2_fortran_ptr + index;
          const double *coeff_ptr = coeff_fortran_ptr + index;
          for (; index < index_fin; index++, tab2_ptr++, coeff_ptr++)
            {
              const int i_colonne = *tab2_ptr;
              assert(i_colonne >= 1 && i_colonne <= vecteur2.size_array());
              const double c = *coeff_ptr;
              const double x2 = vecteur2_fortran_ptr[i_colonne];
              // no shared-item check on columns, see **Annulation items communs**
              x -= c * x2;
            }
          *vecteur_ptr = x;
        }
    }
}

// We compute solution = inverse(C)*b with:
//   inverse(C) = inverse((1/w D - E)) * (2-w/w D) * inverse((1/wD -E)t)
//   D : diagonal part of the matrix, E : lower triangular part of the matrix
void SSOR::ssor(const Matrice_Morse_Sym& matrice, DoubleVect& solution)
{
  const auto& tab2 = matrice.get_tab2();
  if (tab2.size_array() > 0 && tab2[0] == 1)
    {
      // The diagonal is present in the matrix
      if (avec_assert_)
        {
          if (algo_items_communs_)
            {
              descente_assert(omega_, matrice, solution, items_a_traiter_);
              remontee_assert(omega_, matrice, solution, items_a_traiter_);
            }
          else
            {
              descente_assert(omega_, matrice, solution);
              remontee_assert(omega_, matrice, solution);
            }
        }
      else
        {
          if (algo_items_communs_)
            {
              descente(omega_, matrice, solution, items_a_traiter_);
              remontee(omega_, matrice, solution, items_a_traiter_);
            }
          else
            {
              descente(omega_, matrice, solution);
              remontee(omega_, matrice, solution);
            }
        }
    }
  else
    {
      // No stored diagonal: assuming all diagonal entries are 1 (diagonal preconditioning)
      if (avec_assert_)
        {
          if (algo_items_communs_)
            {
              descente_assert_precond_diag(omega_, matrice, solution, items_a_traiter_);
              remontee_assert_precond_diag(omega_, matrice, solution, items_a_traiter_);
            }
          else
            {
              descente_assert_precond_diag(omega_, matrice, solution);
              remontee_assert_precond_diag(omega_, matrice, solution);
            }
        }
      else
        {
          if (algo_items_communs_)
            {
              descente_precond_diag(omega_, matrice, solution, items_a_traiter_);
              remontee_precond_diag(omega_, matrice, solution, items_a_traiter_);
            }
          else
            {
              descente_precond_diag(omega_, matrice, solution);
              remontee_precond_diag(omega_, matrice, solution);
            }
        }
    }
  avec_assert_ = 0; // Do not re-check at the next preconditioning step...
}

void SSOR::ssor(const Matrice_Bloc_Sym& matrice, DoubleVect& solution)
{
  DoubleTab_parts s_parts(solution);
  ConstIntTab_parts items_parts;
  if (algo_items_communs_) items_parts.initialize(items_a_traiter_);

  const int nb_parts = s_parts.size();
  assert(s_parts.size() == nb_parts);
  assert(matrice.nb_bloc_lignes() == nb_parts);

  // Forward sweep
  int i_part;
  for (i_part = 0; i_part < nb_parts; i_part++)
    {
      const Matrice_Bloc& matrice0 = ref_cast(Matrice_Bloc, matrice.get_bloc(i_part, i_part).valeur());
      const Matrice_Morse_Sym& MB00 = ref_cast(Matrice_Morse_Sym, matrice0.get_bloc(0, 0).valeur());
      DoubleVect& partie = s_parts[i_part];

      if (algo_items_communs_) descente_diag_ok_assert(omega_, MB00, partie, items_parts[i_part]);
      else descente_diag_ok_assert(omega_, MB00, partie);

      // off-diagonal blocks
      for (int j_part = i_part + 1; j_part < nb_parts; j_part++)
        {
          const Matrice_Bloc& Aij = ref_cast(Matrice_Bloc, matrice.get_bloc(i_part, j_part).valeur());
          const Matrice_Morse& MB00bis = ref_cast(Matrice_Morse, Aij.get_bloc(0, 0).valeur());
          DoubleVect& partie_j = s_parts[j_part];
          // Warning: the shared-item check concerns matrix rows, not columns (we pass items_parts[i_part], not j_part)
          // (note BM: I believe the previous version was buggy but it went unnoticed because the elem-elem matrix comes first and there are no shared items on elements)
          if (algo_items_communs_) descente_bloc_extradiag_assert(MB00bis, partie, partie_j, items_parts[i_part]);
          else descente_bloc_extradiag_assert(MB00bis, partie, partie_j);
        }
    }
  // Diagonal treatment
  for (i_part = 0; i_part < nb_parts; i_part++)
    {
      const Matrice_Bloc& matrice0 = ref_cast(Matrice_Bloc, matrice.get_bloc(i_part, i_part).valeur());
      const Matrice_Morse_Sym& MB00 = ref_cast(Matrice_Morse_Sym, matrice0.get_bloc(0, 0).valeur());
      DoubleVect& partie = s_parts[i_part];
      traite_diagonale(omega_, MB00, partie);
    }

  // Back-substitution
  for (i_part = nb_parts - 1; i_part >= 0; i_part--)
    {
      const Matrice_Bloc& matrice0 = ref_cast(Matrice_Bloc, matrice.get_bloc(i_part, i_part).valeur());
      const Matrice_Morse_Sym& MB00bis = ref_cast(Matrice_Morse_Sym, matrice0.get_bloc(0, 0).valeur());
      DoubleVect& partie = s_parts[i_part];
      if (algo_items_communs_) remontee_diag_ok_assert(omega_, MB00bis, partie, items_parts[i_part]);
      else remontee_diag_ok_assert(omega_, MB00bis, partie);

      // Off-diagonal blocks (horizontal traversal instead of vertical)
      for (int j_part = 0; j_part < i_part; j_part++)
        {
          const Matrice_Bloc& Aij = ref_cast(Matrice_Bloc, matrice.get_bloc(j_part, i_part).valeur());
          const Matrice_Morse& MB00 = ref_cast(Matrice_Morse, Aij.get_bloc(0, 0).valeur());
          DoubleVect& partie_j = s_parts[j_part];
          // Warning: the shared-item check concerns matrix rows, not columns (we pass items_parts[i_part], not j_part)
          if (algo_items_communs_) remontee_bloc_extradiag_assert(MB00, partie_j, partie, items_parts[j_part]);
          else remontee_bloc_extradiag_assert(MB00, partie_j, partie);
        }
    }
}
