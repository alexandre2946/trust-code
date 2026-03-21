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

#include <Sources_Multiphase_base.h>

#include <Probleme_base.h>
#include <Matrix_tools.h>

Implemente_base(Sources_Multiphase_base, "Sources_Multiphase_base", Source_base);
Sortie& Sources_Multiphase_base::printOn(Sortie& os) const { return os; }
Entree& Sources_Multiphase_base::readOn(Entree& is) { return is; }

/*! @brief Allocates diagonal stencils in the given matrices for field names
 *         listed in diagonal_fields (and optionally "pression" with its
 *         special index mapping).
 */
void Sources_Multiphase_base::dimensionner_blocs_diagonal(matrices_t matrices, const Probleme_base& pb, int ne, int ne_tot, int Nk,
                                                          const std::set<std::string>& diagonal_fields, bool handle_pression)
{
  for (auto &&n_m : matrices)
    {
      bool is_diag = diagonal_fields.count(n_m.first);
      bool is_pression = handle_pression && (n_m.first == "pression");
      if (!is_diag && !is_pression) continue;

      Matrice_Morse& mat = *n_m.second;
      Matrice_Morse mat2;
      const DoubleTab& dep = pb.get_champ(n_m.first.c_str()).valeurs();
      const int nc = dep.dimension_tot(0);
      const int M = dep.line_size();
      Stencil sten(0, 2);

      if (is_diag)
        for (int e = 0; e < ne; e++)
          for (int n = 0; n < Nk; n++)
            if (n < M)
              sten.append_line(Nk * e + n, M * e + n);

      if (is_pression)
        for (int e = 0; e < ne; e++)
          for (int n = 0, m = 0; n < Nk; n++, m += (M > 1))
            sten.append_line(Nk * e + n, M * e + m);

      Matrix_tools::allocate_morse_matrix(Nk * ne_tot, M * nc, sten, mat2);
      mat.nb_colonnes() ? mat += mat2 : mat = mat2;
    }
}
