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

#include <Frontiere_ouverte_rho_u_impose.h>
#include <Fluide_Dilatable_base.h>
#include <Frontiere_dis_base.h>
#include <Equation_base.h>

Implemente_instanciable(Frontiere_ouverte_rho_u_impose, "Frontiere_ouverte_rho_u_impose", Entree_fluide_vitesse_imposee_libre);
// XD frontiere_ouverte_rho_u_impose frontiere_ouverte_vitesse_imposee_sortie frontiere_ouverte_rho_u_impose INHERITS_BRACE This keyword is used to designate a condition of imposed mass rate at an open boundary called bord (edge). The imposed mass rate field at the inlet is vectorial and the imposed velocity values are expressed in kg.s-1. This boundary condition can be used only with the Quasi compressible model.


Sortie& Frontiere_ouverte_rho_u_impose::printOn(Sortie& s) const { return s << que_suis_je() << finl; }

Entree& Frontiere_ouverte_rho_u_impose::readOn(Entree& s) { return Entree_fluide_vitesse_imposee_libre::readOn(s); }

void Frontiere_ouverte_rho_u_impose::completer()
{
  le_fluide = ref_cast(Fluide_Dilatable_base, mon_dom_cl_dis->equation().milieu());
}

int Frontiere_ouverte_rho_u_impose::compatible_avec_eqn(const Equation_base& eqn) const
{
  if (!sub_type(Fluide_Dilatable_base, mon_dom_cl_dis->equation().milieu()))
    {
      Cerr << "The boundary condition Frontiere_ouverte_rho_u_impose is only applicable for dialtable fluids and not a fluid of type " <<  mon_dom_cl_dis->equation().milieu().que_suis_je() << finl;
      Process::exit();
    }

  return Cond_lim_base::compatible_avec_eqn(eqn);
}

double Frontiere_ouverte_rho_u_impose::val_imp_au_temps(double temps, int i) const
{
  Cerr << "Access to a boundary condition in rho.u without specifying the component" << finl;
  Process::exit();
  return 0;
}

double Frontiere_ouverte_rho_u_impose::val_imp_au_temps(double temps, int i, int j) const
{
  double rho_u;
  int ndeb = le_champ_front->frontiere_dis().frontiere().num_premiere_face();
  const DoubleTab& tab_rho_u = le_champ_front->valeurs_au_temps(temps);
  assert(tab_rho_u.nb_dim() == 2);

  if (tab_rho_u.dimension(0) == 1)
    rho_u = tab_rho_u(0, j);
  else
    rho_u = tab_rho_u(i, j);

  double rho = le_fluide->rho_face_np1()(i + ndeb);
  return rho_u / rho;
}

const DoubleTab& Frontiere_ouverte_rho_u_impose::tab_val_imp(double temps) const
{
  if (temps==DMAXFLOAT) temps = le_champ_front->get_temps_defaut();
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());
  // ToDo factorize in Champ_front_base::valeurs_face()
  int size = le_champ_front->valeurs().dimension(0) == 1 ? le_bord.nb_faces_tot() : le_champ_front->valeurs().dimension_tot(0);
  if (size>0)
    {
      bool update = le_champ_front->instationnaire();
      if (tab_.dimension(0) != size)
        {
          tab_.resize(size, le_champ_front->valeurs().dimension(1));
          update = true;
        }
      update = true;  // Provisoire
      if (update)
        {
          int ndeb = le_champ_front->frontiere_dis().frontiere().num_premiere_face();
          int nb_comp = tab_.dimension(1);
          CDoubleTabView rho_u = le_champ_front->valeurs_au_temps(temps).view_ro();
          CDoubleArrView rho_face_np1 = static_cast<const ArrOfDouble&>(le_fluide->rho_face_np1()).view_ro();
          DoubleTabView tab = tab_.view_wo();
          bool uniform = (int)rho_u.extent(0) == 1;
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), size, KOKKOS_LAMBDA (const int face)
          {
            for (int comp = 0; comp < nb_comp; comp++)
              tab(face, comp) = rho_u(uniform ? 0 : face, comp) / rho_face_np1(face + ndeb);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
  return tab_;
}

