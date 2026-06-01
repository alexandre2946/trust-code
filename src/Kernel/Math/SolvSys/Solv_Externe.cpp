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

#include <Solv_Externe.h>
#include <Matrice_Base.h>
#include <Matrice_Morse_Sym.h>
#include <Matrice_Bloc_Sym.h>
#include <Device.h>

Implemente_base_sans_constructeur_ni_destructeur(Solv_Externe,"Solv_Externe",SolveurSys_base);

Sortie& Solv_Externe::printOn(Sortie& s) const
{
  //s << chaine_lue_;
  return s;
}

// readOn
Entree& Solv_Externe::readOn(Entree& is)
{
  //create_solver(is);
  return is;
}

void Solv_Externe::MorseSymToMorse(const Matrice_Morse_Sym& MS, Matrice_Morse& M)
{
  M = MS;
  Matrice_Morse mattmp(MS);
  M.transpose(mattmp);
  int ordre = M.ordre();
  for (int i=0; i<ordre; i++)
    if (M.nb_vois(i))
      M(i, i) = 0.;
  M = mattmp + M;
}

void Solv_Externe::construit_matrice_morse_intermediaire(const Matrice_Base& la_matrice, Matrice_Morse& matrice_morse_intermediaire)
{
  if (sub_type(Matrice_Morse_Sym, la_matrice))
    {
      // Example: pressure matrix in VEFPreP1B
      const Matrice_Morse_Sym& matrice_morse_sym = ref_cast(Matrice_Morse_Sym, la_matrice);
      assert(matrice_morse_sym.get_est_definie());
      MorseSymToMorse(matrice_morse_sym, matrice_morse_intermediaire);
    }
  else if (sub_type(Matrice_Bloc_Sym, la_matrice))
    {
      // Example: pressure matrix in VEF P0+P1+Pa
      const Matrice_Bloc_Sym& matrice = ref_cast(Matrice_Bloc_Sym, la_matrice);
      // Convert the Matrice_Bloc_Sym to Matrice_Morse_Sym format
      Matrice_Morse_Sym matrice_morse_sym;
      matrice.BlocSymToMatMorseSym(matrice_morse_sym);
      MorseSymToMorse(matrice_morse_sym, matrice_morse_intermediaire);
      matrice_morse_sym.dimensionner(0, 0); // Destroy the Morse sym matrix now that it is no longer needed
    }
  else if (sub_type(Matrice_Morse, la_matrice))
    {
      // Example: implicit matrix
      matrice_symetrique_ = false;
    }
  else if (sub_type(Matrice_Bloc, la_matrice))
    {
      // Example: pressure matrix in VDF
      const Matrice_Bloc& matrice_bloc = ref_cast(Matrice_Bloc, la_matrice);
      if (!sub_type(Matrice_Morse_Sym, matrice_bloc.get_bloc(0, 0).valeur()))
        matrice_symetrique_ = false;
      else
        {
          // For a direct solver: fix the matrix if it is not defined (in incompressible VDF, nothing was done before...)
          Matrice_Morse_Sym& mat00 = ref_cast_non_const(Matrice_Morse_Sym, matrice_bloc.get_bloc(0, 0).valeur());
          if (solveur_direct() && mat00.get_est_definie() == 0 && Process::je_suis_maitre())
            mat00(0, 0) *= 2;
        }
      matrice_bloc.BlocToMatMorse(matrice_morse_intermediaire);
    }
  else
    {
      Cerr << "Error, we do not know yet treat a matrix of type " << la_matrice.que_suis_je() << finl;
      exit();
    }
}

const ArrOfInt& Solv_Externe::indice_coeff_to_keep(const Matrice_Morse& matrice_morse)
{
  if (!indice_coeff_to_keep_.size_array())
    {
      // Build indice_coeff_to_keep_
      const auto& tab1 = matrice_morse.get_tab1();
      const int n = tab1.size_array() - 1;
      auto nnz(tab1[0]);
      nnz = 0;
      for (int i = 0; i < n; i++)
        {
          if (items_to_keep_[i])
            nnz += tab1[i + 1] - tab1[i];
        }
      indice_coeff_to_keep_.resize((int)nnz);
      nnz = 0;
      for (int i = 0; i < n; i++)
        {
          if (items_to_keep_[i])
            {
              const auto k0 = tab1[i] - 1;
              const auto k1 = tab1[i + 1] - 1;
              for (auto k = k0; k < k1; k++)
                {
                  indice_coeff_to_keep_[(int)nnz] = (int)k;
                  nnz++;
                }
            }
        }
    }
  return indice_coeff_to_keep_;
}

void Solv_Externe::Create_lhs_rhs_onDevice()
{
  lhs_.resize(nb_rows_);
  rhs_.resize(nb_rows_);
}

template<typename ExecSpace>
void Solv_Externe::Update_lhs_rhs(const DoubleVect& tab_b, DoubleVect& tab_x)
{
  const unsigned int size = tab_b.size_array();
  auto x = tab_x.template view_ro<1, ExecSpace>();
  auto b = tab_b.template view_ro<1, ExecSpace>();
  auto index = static_cast<const ArrOfInt&>(index_).template view_ro<1, ExecSpace>();
  auto lhs = lhs_.template view_wo<1, ExecSpace>();
  auto rhs = rhs_.template view_wo<1, ExecSpace>();
  Kokkos::RangePolicy<ExecSpace> policy(0, size);
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), policy, KOKKOS_LAMBDA(
                         const int i)
  {
    int ind = index[i];
    if (ind != -1)
      {
        lhs[ind] = x[i];
        rhs[ind] = b[i];
      }
  });
  static constexpr bool kernelOnDevice = !std::is_same<ExecSpace, Kokkos::DefaultHostExecutionSpace>::value;
  end_gpu_timer(__KERNEL_NAME__, kernelOnDevice);
}

template<typename ExecSpace>
void Solv_Externe::Update_solution(DoubleVect& tab_x)
{
  const unsigned int size = tab_x.size_array();
  auto index = static_cast<const ArrOfInt&>(index_).template view_ro<1, ExecSpace>();
  auto lhs = lhs_.template view_ro<1, ExecSpace>();
  auto x = tab_x.template view_wo<1, ExecSpace>();
  Kokkos::RangePolicy<ExecSpace> policy(0, size);
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), policy, KOKKOS_LAMBDA(
                         const int i)
  {
    int ind = index[i];
    if (ind != -1)
      x[i] = lhs[ind];
  });
  static constexpr bool kernelOnDevice = !std::is_same<ExecSpace, Kokkos::DefaultHostExecutionSpace>::value;
  end_gpu_timer(__KERNEL_NAME__, kernelOnDevice);
}

#ifdef TRUST_USE_GPU
template void Solv_Externe::Update_lhs_rhs<Kokkos::DefaultExecutionSpace>(const DoubleVect&, DoubleVect&);
template void Solv_Externe::Update_solution<Kokkos::DefaultExecutionSpace>(DoubleVect&);
#endif
template void Solv_Externe::Update_lhs_rhs<Kokkos::DefaultHostExecutionSpace>(const DoubleVect&, DoubleVect&);
template void Solv_Externe::Update_solution<Kokkos::DefaultHostExecutionSpace>(DoubleVect&);

