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

#include <Loi_Etat_GP_QC.h>
#include <Fluide_Dilatable_base.h>

Implemente_instanciable(Loi_Etat_GP_QC,"Loi_Etat_Gaz_Parfait_QC",Loi_Etat_Mono_GP_base);
// XD perfect_gaz_QC loi_etat_gaz_parfait_base gaz_parfait_QC BRACE Class for perfect gas state law used with a
// XD_CONT quasi-compressible fluid.
// XD attr Cp double Cp REQ Specific heat at constant pressure (J/kg/K).
// XD attr Cv double Cv OPT Specific heat at constant volume (J/kg/K).
// XD attr gamma double gamma OPT Cp/Cv
// XD attr Prandtl double Prandtl REQ Prandtl number of the gas Pr=mu*Cp/lambda
// XD attr rho_constant_pour_debug field_base rho_constant_pour_debug OPT For developers to debug the code with a
// XD_CONT constant rho.

Sortie& Loi_Etat_GP_QC::printOn(Sortie& os) const
{
  os <<que_suis_je()<< finl;
  return os;
}

Entree& Loi_Etat_GP_QC::readOn(Entree& is)
{
  Cerr << "Lecture de la loi d'etat gaz parfait pour le QC ... " << finl;
  return Loi_Etat_Mono_GP_base::readOn(is);
}

double Loi_Etat_GP_QC::calculer_masse_volumique(double P, double T) const
{
  return rho_constant_pour_debug_ ? rho_constant_pour_debug_->valeurs()(0,0) :
         Loi_Etat_Mono_GP_base::calculer_masse_volumique(P,T);
}

// View version
void Loi_Etat_GP_QC::compute_tab_rho(DoubleTab& tab_rho)
{
  double rho_constant = rho_constant_pour_debug_ ? rho_constant_pour_debug_->valeurs()(0,0) : 0;
  double Pth = le_fluide->pression_th();
  double R = R_;
  CDoubleArrView tab_ICh = static_cast<const ArrOfDouble&>(le_fluide->inco_chaleur().valeurs()).view_ro();
  CDoubleArrView rho_n = static_cast<const ArrOfDouble&>(tab_rho_n).view_ro();
  DoubleArrView rho_np1 = static_cast<ArrOfDouble&>(tab_rho_np1).view_wo();
  DoubleArrView rho = static_cast<ArrOfDouble&>(tab_rho).view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), tab_rho.size(), KOKKOS_LAMBDA(const int som)
  {
    rho_np1(som) = rho_constant ? rho_constant : Loi_Etat_Mono_GP_base::calculer_masse_volumique(Pth, tab_ICh(som), R);
    rho(som) = 0.5 * (rho_n(som) + rho_np1(som));
  });
  end_gpu_timer(__KERNEL_NAME__);
}

/*! @brief Calcule la masse volumique
 *
 */
void Loi_Etat_GP_QC::calculer_masse_volumique()
{
  Loi_Etat_Mono_GP_base::calculer_masse_volumique();
}
