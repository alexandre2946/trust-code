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

#include <EDO_Pression_th_VEF_Gaz_Parfait.h>
#include <Frontiere_ouverte_rho_u_impose.h>
#include <Fluide_Quasi_Compressible.h>
#include <Neumann_sortie_libre.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>
#include <Loi_Etat_GP_QC.h>
#include <Domaine_VEF.h>
#include <Champ_P1NC.h>
#include <Dirichlet.h>

Implemente_instanciable(EDO_Pression_th_VEF_Gaz_Parfait, "EDO_Pression_th_VEF_Gaz_Parfait", EDO_Pression_th_VEF);

Sortie& EDO_Pression_th_VEF_Gaz_Parfait::printOn(Sortie& os) const { return os << que_suis_je() << finl; }

Entree& EDO_Pression_th_VEF_Gaz_Parfait::readOn(Entree& is) { return is; }

/*! @brief Solves the ODE.
 *
 * @param (double Pth_n) The pressure at the previous time step
 * @return (double) The new pressure value
 */
double EDO_Pression_th_VEF_Gaz_Parfait::resoudre(double Pth_n)
{

  int traitPth = le_fluide_->getTraitementPth();

  if (traitPth == 2)  // Pth constant
    return Pth_n;

  double present = le_fluide_->vitesse().equation().schema_temps().temps_courant();
  double dt = le_fluide_->vitesse().equation().schema_temps().pas_de_temps();
  double futur = present + dt;
  const DoubleTab& tab_tempnp1 = le_fluide_->inco_chaleur().valeurs(futur);    // T(n+1)
  const DoubleTab& tab_tempn = le_fluide_->inco_chaleur().valeurs(present);    // T(n)
  int nb_faces = le_dom->nb_faces();
  double Pth = 0;

  const Domaine_VEF& domaine_vef = ref_cast(Domaine_VEF, le_dom.valeur());

  if (traitPth == 0)   // ODE
    {
      const DoubleTab& tab_rho = le_fluide_->masse_volumique().valeurs();       // n+1/2
      const double rho_moy = Champ_P1NC::calculer_integrale_volumique(domaine_vef, tab_rho, FAUX_EN_PERIO);
      DoubleTrav tab_tmp(tab_rho); // copy of rho
      tab_tmp = tab_rho;
      assert(tab_tmp.size() == nb_faces);
      const double invdt = 1. / dt;
      CDoubleArrView tempn = static_cast<const ArrOfDouble&>(tab_tempn).view_ro();
      CDoubleArrView tempnp1 = static_cast<const ArrOfDouble&>(tab_tempnp1).view_ro();
      DoubleArrView tmp = static_cast<ArrOfDouble&>(tab_tmp).view_wo();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, nb_faces), KOKKOS_LAMBDA(const int i)
      {
        tmp[i] *= (tempnp1[i] - tempn[i]) / tempnp1[i] * invdt;
      });
      end_gpu_timer(__KERNEL_NAME__);
      double S = Champ_P1NC::calculer_integrale_volumique(domaine_vef, tab_tmp, FAUX_EN_PERIO);

      //////////////////
      S /= rho_moy;
      //    double Pth = (Pth_n + S*dt/V);
      //double Pth = Pth_n/(1- S*dt);
      Pth = Pth_n / (1 - S * dt);
      //test to reveal the difference on Pth
      //Pth *= 1e10;
    }

  else if (traitPth == 1)   // Mass conservation (WEC March 2008)
    {

      // We want Mass(n+1) - Mass(n) = dt * Boundary_fluxes
      // where Mass = Pth / R * sum(Vi/Ti) in volume
      // and Flux = Pth / R * sum(Sj.Uj/Tj) at Dirichlet faces

      // So if we denote masse_n = sum(Vi/Ti(n))
      // debit_u_imp = sum(Sj.Uj(imposed)/Tj(n+1)) for imposed velocities
      // and debit_rho_u_imp = sum(Sj.Uj(imposed)/Tj(n+1)) for imposed rho_u
      //                      where U(imposed) = (rho.U)(imposed) / rho(Pth(n),T(n+1))
      // we need Pth(n+1)*masse_np1 - Pth(n)*masse_n = Pth(n+1) * debit_u_imp * dt
      //                                             + Pth(n) * debit_rho_u_imp * dt

      // Warning: porosity should be taken into account in the integrals !!!

      double debit_u_imp = 0, debit_rho_u_imp = 0;

      // We want rho_np1 to be computed with T(n+1) and Pth_n for boundary conditions with imposed rho_u
      le_fluide_->calculer_masse_volumique();

      // Computation of masse_n and masse_np1
      DoubleTrav tab_tmp;
      tab_tmp.copy(tab_tempn, RESIZE_OPTIONS::NOCOPY_NOINIT); // copy the structure only

      CDoubleArrView tempn = static_cast<const ArrOfDouble&>(tab_tempn).view_ro();
      DoubleArrView tmp = static_cast<ArrOfDouble&>(tab_tmp).view_wo();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, nb_faces), KOKKOS_LAMBDA(const int i)
      {
        tmp(i) = 1. / tempn(i);
      });
      end_gpu_timer(__KERNEL_NAME__);
      const double masse_n = Champ_P1NC::calculer_integrale_volumique(domaine_vef, tab_tmp, FAUX_EN_PERIO);

      CDoubleArrView tempnp1 = static_cast<const ArrOfDouble&>(tab_tempnp1).view_ro();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, nb_faces), KOKKOS_LAMBDA(const int i)
      {
        tmp(i) = 1. / tempnp1(i);
      });
      end_gpu_timer(__KERNEL_NAME__);
      const double masse_np1 = Champ_P1NC::calculer_integrale_volumique(domaine_vef, tab_tmp, FAUX_EN_PERIO);

      // Computation of debit_u_imp and debit_rho_u_imp
      for (int n_bord = 0; n_bord < le_dom->nb_front_Cl(); n_bord++)
        {
          const Cond_lim_base& la_cl = le_dom_Cl->les_conditions_limites(n_bord).valeur();
          if (sub_type(Dirichlet, la_cl))
            {
              const Dirichlet& diri = ref_cast(Dirichlet, la_cl);
              const Front_VF& la_front_dis = ref_cast(Front_VF, la_cl.frontiere_dis());
              int ndeb = la_front_dis.num_premiere_face();
              int nfin = ndeb + la_front_dis.nb_faces();
              int dim = Objet_U::dimension;
              CDoubleTabView val_imp = diri.tab_val_imp(futur).view_ro();
              CDoubleTabView face_normales = le_dom->face_normales().view_ro();
              CIntTabView face_voisins = le_dom->face_voisins().view_ro();
              double debit = 0;
              Kokkos::parallel_reduce(start_gpu_timer(__KERNEL_NAME__), range_1D(ndeb, nfin), KOKKOS_LAMBDA(const int face, double& debit_local)
              {
                double debit_v = 0;
                for (int d = 0; d < dim; d++)
                  debit_v += face_normales(face, d) * val_imp(face - ndeb, d) / tempnp1(face);
                int n0 = face_voisins(face, 0);
                if (n0 == -1)
                  debit_v = -debit_v;
                debit_local += debit_v;
              }, debit);
              end_gpu_timer(__KERNEL_NAME__);
              if (sub_type(Frontiere_ouverte_rho_u_impose, la_cl))
                debit_rho_u_imp += debit;
              else
                debit_u_imp += debit;
            }
          else if (sub_type(Neumann_sortie_libre, la_cl))
            {
              Cerr << la_cl.que_suis_je() << " est incompatible avec le traitement conservation_masse." << finl;
              exit();
            }
        }
      // Sum over all procs
      // Optimization: combine 2 mp_sum into 1 collective call
      mp_sum_for_each(debit_u_imp, debit_rho_u_imp);

      // Computation of Pth(n+1)
      Pth = Pth_n * (masse_n - dt * debit_rho_u_imp) / (masse_np1 + dt * debit_u_imp);
    }

  return Pth;
}

void EDO_Pression_th_VEF_Gaz_Parfait::resoudre(DoubleTab& Pth_n)
{

  const int traitPth = le_fluide_->getTraitementPth();
  if (traitPth == 2)
    return; // nothing to do
  else if (traitPth == 0)
    {
      for (int n_bord = 0; n_bord < le_dom->nb_front_Cl(); n_bord++)
        {
          const Cond_lim& la_cl = le_dom_Cl->les_conditions_limites(n_bord);
          if (sub_type(Neumann_sortie_libre, la_cl.valeur()))
            return; // nothing to do
        }

      Cerr << "EDO_Pression_th_VEF_Gaz_Parfait::" << __func__ << " not yet coded ! Call the 911 !!" << finl;
      Process::exit();
    }
  else
    {
      const double present = le_fluide_->vitesse().equation().schema_temps().temps_courant();
      const double dt = le_fluide_->vitesse().equation().schema_temps().pas_de_temps();
      const double futur = present + dt;
      const DoubleTab& tempnp1 = le_fluide_->inco_chaleur().valeurs(futur);    // T(n+1)
      const DoubleTab& tempn = le_fluide_->inco_chaleur().valeurs(present);    // T(n)
      const int nb_faces = le_dom->nb_faces();
      const Domaine_VEF& domaine_vef = ref_cast(Domaine_VEF, le_dom.valeur());

      // We want Mass(n+1) - Mass(n) = dt * Boundary_fluxes
      // where Mass = Pth / R * sum(Vi/Ti) in volume
      // and Flux = Pth / R * sum(Sj.Uj/Tj) at Dirichlet faces

      // So if we denote masse_n = sum(Vi/Ti(n))
      // debit_u_imp = sum(Sj.Uj(imposed)/Tj(n+1)) for imposed velocities
      // and debit_rho_u_imp = sum(Sj.Uj(imposed)/Tj(n+1)) for imposed rho_u
      //                      where U(imposed) = (rho.U)(imposed) / rho(Pth(n),T(n+1))
      // we need Pth(n+1)*masse_np1 - Pth(n)*masse_n = Pth(n+1) * debit_u_imp * dt
      //                                             + Pth(n) * debit_rho_u_imp * dt

      // Warning: porosity should be taken into account in the integrals !!!

      double debit_u_imp = 0, debit_rho_u_imp = 0;

      // We want rho_np1 to be computed with T(n+1) and Pth_n for boundary conditions with imposed rho_u
      le_fluide_->calculer_masse_volumique();

      // Computation of masse_n and masse_np1
      DoubleVect tmp;
      tmp.copy(tempn, RESIZE_OPTIONS::NOCOPY_NOINIT); // copy the structure only
      for (int i = 0; i < nb_faces; i++)
        tmp[i] = 1. / tempn[i];
      const double masse_n = Champ_P1NC::calculer_integrale_volumique(domaine_vef, tmp, FAUX_EN_PERIO);
      for (int i = 0; i < nb_faces; i++)
        tmp[i] = 1. / tempnp1[i];
      const double masse_np1 = Champ_P1NC::calculer_integrale_volumique(domaine_vef, tmp, FAUX_EN_PERIO);

      // Computation of debit_u_imp and debit_rho_u_imp
      for (int n_bord = 0; n_bord < le_dom->nb_front_Cl(); n_bord++)
        {
          const Cond_lim_base& la_cl = le_dom_Cl->les_conditions_limites(n_bord).valeur();
          if (sub_type(Dirichlet, la_cl))
            {
              const Dirichlet& diri = ref_cast(Dirichlet, la_cl);
              const Front_VF& la_front_dis = ref_cast(Front_VF, la_cl.frontiere_dis());
              int ndeb = la_front_dis.num_premiere_face();
              int nfin = ndeb + la_front_dis.nb_faces();
              for (int face = ndeb; face < nfin; face++)
                {
                  double debit_v = 0;
                  for (int d = 0; d < dimension; d++)
                    debit_v += le_dom->face_normales(face, d) * diri.val_imp_au_temps(futur, face - ndeb, d) / tempnp1(face);
                  int n0 = le_dom->face_voisins(face, 0);
                  if (n0 == -1)
                    debit_v = -debit_v;
                  if (sub_type(Frontiere_ouverte_rho_u_impose, la_cl))
                    debit_rho_u_imp += debit_v;
                  else
                    debit_u_imp += debit_v;
                }
            }
          else if (sub_type(Neumann_sortie_libre, la_cl))
            {
              Cerr << la_cl.que_suis_je() << " est incompatible avec le traitement conservation_masse." << finl;
              exit();
            }
        }
      // Sum over all procs
      // Optimization: combine 2 mp_sum into 1 collective call
      mp_sum_for_each(debit_u_imp, debit_rho_u_imp);

      // Computation of Pth(n+1)
      for (int f = 0; f < nb_faces; f++)
        Pth_n(f) = Pth_n(f) * (masse_n - dt * debit_rho_u_imp) / (masse_np1 + dt * debit_u_imp);
    }
}
