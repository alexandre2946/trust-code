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

#ifndef Conservation_Euler_base_included
#define Conservation_Euler_base_included

#include <Convection_Diffusion_std.h>
#include <Operateur_NConserv.h>

class Fluide_base;

class Conservation_Euler_base : public Convection_Diffusion_std
{
  Declare_base(Conservation_Euler_base);
public:
  int nombre_d_operateurs() const override
  {
    Process::exit("Conservation_Euler_base::nombre_d_operateurs !!!  \n");
    return 1;
  }

  const Operateur& operateur(int) const override
  {
    Process::exit("Conservation_Euler_base::operateur !!!  \n");
    return terme_convectif;
  }

  Operateur& operateur(int) override
  {
    Process::exit("Conservation_Euler_base::operateur !!!  \n");
    return terme_convectif;
  }

  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  void completer() override { Equation_base::completer(); }

  inline const Champ_Inc_base& inconnue() const override { return l_inco_ch_; }
  inline Champ_Inc_base& inconnue() override { return l_inco_ch_; }
  inline double calculer_pas_de_temps() const override { return 1e10; }

  void associer_milieu_base(const Milieu_base&) override;
  void associer_fluide(const Fluide_base&);
  const Milieu_base& milieu() const override;
  Milieu_base& milieu() override;
  const Fluide_base& fluide() const;
  Fluide_base& fluide();

  virtual double flux_bord(const double inco_bord, const double vit_n_bord, const double p_bord) const
  {
    Process::exit("Conservation_Euler_base::flux_bord !!!  \n");
    return 0;
  }

  virtual double termes_NonConservatif(const double alpha_bord, const double vitesse_n_inter, const double p_bord) const
  {
    Process::exit("Conservation_Euler_base::termes_NonConservatif !!!  \n");
    return 0;
  }

  void dimensionner_matrice_sans_mem(Matrice_Morse& matrice) override
  {
    Process::exit("Conservation_Euler_base::dimensionner_matrice_sans_mem !!!  \n");
  }

  int has_interface_blocs() const override
  {
    Process::exit("Conservation_Euler_base::has_interface_blocs !!!  \n");
    return -1;
  }

  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = { }) const override
  {
    Process::exit("Conservation_Euler_base::dimensionner_blocs !!!  \n");
  }

  void assembler_blocs_avec_inertie(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = { }) override
  {
    Process::exit("Conservation_Euler_base::assembler_blocs_avec_inertie !!!  \n");
  }

protected :
  Operateur_NConserv terme_nconserv_;
  OWN_PTR(Champ_Inc_base) l_inco_ch_;
  OBS_PTR(Fluide_base) le_fluide_;
};

#endif /* Conservation_Euler_base_included */
