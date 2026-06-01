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

#include <Navier_Stokes_Fluide_Dilatable_Proto.h>
#include <Transport_Interfaces_base.h>
#include <Fluide_Dilatable_base.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
#include <Dirichlet.h>
#include <TRUSTTrav.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Debog.h>
#include <kokkos++.h>
#include <Perf_counters.h>

Navier_Stokes_Fluide_Dilatable_Proto::Navier_Stokes_Fluide_Dilatable_Proto() : cumulative_(0) { }

// Multiply density by velocity and return density*velocity (mass flux)
DoubleTab& Navier_Stokes_Fluide_Dilatable_Proto::rho_vitesse_impl(const DoubleTab& tab_rho, const DoubleTab& vit,
                                                                  DoubleTab& rhovitesse) const
{
  const int n = vit.dimension(0), ncomp = vit.line_size();
  CDoubleTabView tab_rho_v = tab_rho.view_ro();
  CDoubleTabView vit_v = vit.view_ro();
  DoubleTabView rhovitesse_v = rhovitesse.view_wo();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), n, KOKKOS_LAMBDA(
                         const int i)
  {
    for (int j=0 ; j<ncomp ; j++)
      rhovitesse_v(i, j) = tab_rho_v(i, 0) * vit_v(i, j);
  });
  end_gpu_timer(__KERNEL_NAME__);

  rhovitesse.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_Fluide_Dilatable_Proto::rho_vitesse : ", rhovitesse);
  return rhovitesse;
}

int Navier_Stokes_Fluide_Dilatable_Proto::impr_impl(const Navier_Stokes_std& eqn,Sortie& os) const
{
  const Fluide_Dilatable_base& fluide_dil = ref_cast(Fluide_Dilatable_base,eqn.fluide());
  const DoubleTab& vit = eqn.vitesse().valeurs(), &rho = fluide_dil.rho_face_np1();
  DoubleTrav mass_flux(vit);
  rho_vitesse_impl(rho,vit,mass_flux);

  DoubleTrav array(eqn.div().valeurs());
  if (tab_W.get_md_vector())
    {
      operator_egal(array, tab_W ); //, VECT_REAL_ITEMS); // initialise
      array*=-1;
    }
  else
    {
      // Note (B.M.): some implementations of this method do not perform virtual space exchange:
      fluide_dil.secmembre_divU_Z(array);
      array*=-1;
    }

  array.echange_espace_virtuel();
  eqn.operateur_divergence().ajouter(mass_flux, array);
  double LocalMassFlowRateError = mp_max_abs_vect(array); // max|sum(rho*u*ndS)|

  os << "-------------------------------------------------------------------"<< finl;
  os << "Cell balance mass flow rate control for the problem " << eqn.probleme().le_nom() << " : " << finl;
  os << "Absolute value : " << LocalMassFlowRateError << " kg/s" << finl; ;

  // Divide array by vol(i)
  eqn.operateur_divergence().volumique(array);

  // Divide by a mean rho
  double rho_moyen = mp_moyenne_vect(rho), dt = eqn.probleme().schema_temps().pas_de_temps();
  double bilan_massique_relatif = mp_max_abs_vect(array) * dt / rho_moyen;
  os << "Relative value : " << bilan_massique_relatif << finl; // =max|LocalMassFlowRateError/(rho_moyen*Vol(i)/dt)|

  // Calculation as OpenFOAM: http://foam.sourceforge.net/docs/cpp/a04190_source.html
  // It is relative errors (normalized by the volume/dt)
  double TotalMass = rho_moyen * eqn.probleme().domaine().volume_total();
  double local = LocalMassFlowRateError / ( TotalMass / dt ), global = mp_somme_vect(array) / ( TotalMass / dt );
  cumulative_ += global;

  os << "time step continuity errors : sum local = " << local << ", global = " << global << ", cumulative = " << cumulative_ << finl;

  if (local > 0.01)
    {
      Cerr << "The mass balance is too bad (relative value > 1%)." << finl;
      Cerr << "Please check and lower the convergence value of the pressure solver." << finl;
      Process::exit();
    }
  return 1;
}

/*! @brief @brief Computes the time derivative of the velocity unknown, i.e. the acceleration dU/dt, and returns it.
 *
 * Calls Equation_base::derivee_en_temps_inco(DoubleTab&) and also computes the pressure.
 *
 * @param vpoint Array of acceleration values dU/dt.
 * @return Array of acceleration values (velocity derivative).
 */
DoubleTab& Navier_Stokes_Fluide_Dilatable_Proto::derivee_en_temps_inco_impl(Navier_Stokes_std& eqn,DoubleTab& vpoint)
{
  const Fluide_Dilatable_base& fluide_dil=ref_cast(Fluide_Dilatable_base,eqn.milieu());
  DoubleTab& press = eqn.pression().valeurs(), &vit = eqn.vitesse().valeurs();
  DoubleTrav secmem(press);
  DoubleTrav inc_pre(press);
  DoubleTrav rhoU(vit);

  if (!tab_W.get_md_vector())
    {
      tab_W.copy(secmem, RESIZE_OPTIONS::NOCOPY_NOINIT); // copy structure
      // initialisation to avoid assert failure when filling in EDO_Pression_th_VEF::secmembre_divU_Z_VEFP1B
      tab_W = 0.;
    }

  // Get champ gradP
  OBS_PTR(Champ_base) gradient_pression;
  eqn.has_champ("gradient_pression", gradient_pression);

  if (!gradient_pression)
    {
      Cerr<<"l'equation ne comprend pas gradient_pression "<<finl;
      Process::exit();
    }

  DoubleTab& gradP = gradient_pression->valeurs();
  DoubleTrav Mmoins1grad(gradP);

  // We use the incremental pressure-projection algorithm (Chorin)

  // Step 1 : prepare operators and solve for a provisional velocity u*
  prepare_and_solve_u_star(eqn,fluide_dil,rhoU, vpoint);

  // Step 2 : solve the poisson equation for the pressure increment (Pi = P^n+1 - P^n)
  // Attention the matrix has a constant coefficient as the variable of NS is rhoU and not U !!!
  solve_pressure_increment(eqn,fluide_dil,rhoU,secmem,inc_pre,vpoint);

  // Step 3 : compute P^n+1 & compute the correct velocity u^n+1
  correct_and_compute_u_np1(eqn,fluide_dil,rhoU,Mmoins1grad,inc_pre,gradP,vpoint);

  return vpoint;
}

void Navier_Stokes_Fluide_Dilatable_Proto::assembler_avec_inertie_impl(const Navier_Stokes_std& eqn, Matrice_Morse& mat_morse,
                                                                       const DoubleTab& present, DoubleTab& tab_secmem)
{
  // ******   before inertia   ******
  // diffusion in div(mu grad u), but we want to implicitize in rho * u => divide contributions by the associated rho_face
  // GF: add after contributing so as to have the correct boundary fluxes
  DoubleTrav rhovitesse(present);

  // Op diff
  eqn.operateur(0).l_op_base().contribuer_a_avec(present,mat_morse);
  eqn.operateur(0).ajouter(tab_secmem);

  const Fluide_Dilatable_base& fluide_dil=ref_cast(Fluide_Dilatable_base,eqn.milieu());
  const DoubleTab& tab_rho_face_np1 = fluide_dil.rho_face_np1(), &tab_rho_face_n=fluide_dil.rho_face_n();
  const int nb_compo = present.line_size();

  auto tab1 = mat_morse.get_tab1().view_ro();
  CIntArrView tab2 = mat_morse.get_tab2().view_ro();
  CDoubleArrView rho_face_np1 = static_cast<const ArrOfDouble&>(tab_rho_face_np1).view_ro();
  DoubleArrView coeff = mat_morse.get_set_coeff().view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), mat_morse.nb_lignes(), KOKKOS_LAMBDA(const int i)
  {
    for (auto k=tab1(i)-1; k<tab1(i+1)-1; k++)
      {
        int j = tab2(k)-1;
        double rapport = rho_face_np1(j/nb_compo);
        coeff(k) /= rapport;
      }
  });
  end_gpu_timer(__KERNEL_NAME__);

  rho_vitesse_impl(tab_rho_face_np1,present,rhovitesse); // rho*U

  // Op conv
  eqn.operateur(1).l_op_base().contribuer_a_avec(rhovitesse,mat_morse);
  eqn.operateur(1).ajouter(rhovitesse,tab_secmem);

  // sources
  eqn.sources().ajouter(tab_secmem);
  eqn.sources().contribuer_a_avec(present,mat_morse);

  // solve in rho*u, so store rho*u in present
  rho_vitesse_impl(tab_rho_face_np1,present,ref_cast_non_const(DoubleTab,present));
  mat_morse.ajouter_multvect(present, tab_secmem);

  /*
   * contribution to the inertia matrix:
   * divide the diagonal by rho^{n+1}_face
   * add inertia in the standard way
   * multiply the diagonal back by rho^{n+1}
   */

  // add inertia
  const double dt=eqn.schema_temps().pas_de_temps();
  eqn.solv_masse().ajouter_masse(dt,mat_morse,0);

  rho_vitesse_impl(tab_rho_face_n,eqn.inconnue().passe(),rhovitesse);
  eqn.solv_masse().ajouter_masse(dt,tab_secmem,rhovitesse,0);

  // boundary condition locking wrong if Dirichlet u!=0 !!!!!! missing multiplication by rho
  for (int op=0; op< eqn.nombre_d_operateurs(); op++) eqn.operateur(op).l_op_base().modifier_pour_Cl(mat_morse,tab_secmem);

  /*
   * final correction for Dirichlet conditions:
   * we must not impose u^{n+1} but rho*u^{n+1} => multiply the result by rho_face_np1
   */
  const Conds_lim& lescl=eqn.domaine_Cl_dis().les_conditions_limites();

  for (auto& itr : lescl)
    {
      const Cond_lim_base& la_cl_base = itr.valeur();
      if (sub_type(Dirichlet,la_cl_base))
        {
          const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl_base.frontiere_dis());
          int ndeb = la_front_dis.num_premiere_face();
          int nfin = ndeb + la_front_dis.nb_faces();
          int dim = present.line_size()==1 ? 1 : Objet_U::dimension;
          DoubleTabView secmem = tab_secmem.view_rw();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(ndeb, nfin), KOKKOS_LAMBDA(const int num_face)
          {
            for (int dir=0; dir<dim; dir++)
              secmem(num_face,dir)*=rho_face_np1(num_face);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
}


void Navier_Stokes_Fluide_Dilatable_Proto::assembler_blocs_avec_inertie(const Navier_Stokes_std& eqn, matrices_t matrices, DoubleTab& tab_secmem, const tabs_t& semi_impl)
{
  statistics().begin_count(STD_COUNTERS::ajouter_blocs,statistics().get_last_opened_counter_level()+1);
  const std::string& nom_inco = eqn.inconnue().le_nom().getString();
  Matrice_Morse *mat = matrices.count(nom_inco)?matrices.at(nom_inco):nullptr;
  const DoubleTab& present = eqn.inconnue().valeurs();

  // ******   before inertia   ******
  // diffusion in div(mu grad u), but we want to implicitize in rho * u => divide contributions by the associated rho_face
  // GF: add after contributing so as to have the correct boundary fluxes
  DoubleTrav rhovitesse(present);

  // Op diff
  eqn.operateur(0).l_op_base().ajouter_blocs(matrices, tab_secmem, semi_impl);

  const Fluide_Dilatable_base& fluide_dil=ref_cast(Fluide_Dilatable_base,eqn.milieu());
  const DoubleTab& tab_rho_face_np1 = fluide_dil.rho_face_np1(), &tab_rho_face_n=fluide_dil.rho_face_n();
  const int nb_compo = present.line_size();

  auto tab1 = mat->get_tab1().view_ro();
  CIntArrView tab2 = mat->get_tab2().view_ro();
  CDoubleArrView rho_face_np1 = static_cast<const ArrOfDouble&>(tab_rho_face_np1).view_ro();
  DoubleArrView coeff = mat->get_set_coeff().view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), mat->nb_lignes(), KOKKOS_LAMBDA(const int i)
  {
    for (auto k=tab1(i)-1; k<tab1(i+1)-1; k++)
      {
        int j = tab2(k)-1;
        double rapport = rho_face_np1(j/nb_compo);
        coeff(k) /= rapport;
      }
  });
  end_gpu_timer(__KERNEL_NAME__);

  rho_vitesse_impl(tab_rho_face_np1,present,rhovitesse); // rho*U

  // Op conv
  eqn.operateur(1).l_op_base().ajouter_blocs(matrices, tab_secmem, {{nom_inco,rhovitesse}});
  statistics().end_count(STD_COUNTERS::ajouter_blocs);

  // sources
  statistics().begin_count(STD_COUNTERS::source_terms,statistics().get_last_opened_counter_level()+1);
  for (int i = 0; i < eqn.sources().size(); i++)
    eqn.sources()(i)->ajouter_blocs(matrices, tab_secmem, semi_impl);
  statistics().end_count(STD_COUNTERS::source_terms);

  statistics().begin_count(STD_COUNTERS::ajouter_blocs,statistics().get_last_opened_counter_level()+1);
  // solve in rho*u, so store rho*u in present
  rho_vitesse_impl(tab_rho_face_np1,present,ref_cast_non_const(DoubleTab,present));
  mat->ajouter_multvect(present,tab_secmem);
  eqn.operateur_gradient()->ajouter_blocs(matrices, tab_secmem, semi_impl);

  /*
   * contribution to the inertia matrix:
   * divide the diagonal by rho^{n+1}_face
   * add inertia in the standard way
   * multiply the diagonal back by rho^{n+1}
   */

  // add inertia
  const double dt=eqn.schema_temps().pas_de_temps();
  eqn.solv_masse().ajouter_masse(dt,*mat,0);
  rho_vitesse_impl(tab_rho_face_n,eqn.inconnue().passe(),rhovitesse);
  eqn.solv_masse().ajouter_masse(dt,tab_secmem,rhovitesse,0);

  // boundary condition locking wrong if Dirichlet u!=0 !!!!!! missing multiplication by rho
  for (int op=0; op< eqn.nombre_d_operateurs(); op++) eqn.operateur(op).l_op_base().modifier_pour_Cl(*mat,tab_secmem);

  /*
   * final correction for Dirichlet conditions:
   * we must not impose u^{n+1} but rho*u^{n+1} => multiply the result by rho_face_np1
   */
  const Conds_lim& lescl=eqn.domaine_Cl_dis().les_conditions_limites();

  for (auto& itr : lescl)
    {
      const Cond_lim_base& la_cl_base = itr.valeur();
      if (sub_type(Dirichlet,la_cl_base))
        {

          const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl_base.frontiere_dis());
          int ndeb = la_front_dis.num_premiere_face();
          int nfin = ndeb + la_front_dis.nb_faces();
          int dim = present.line_size()==1 ? 1 : Objet_U::dimension;
          DoubleTabView secmem = tab_secmem.view_rw();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(ndeb, nfin), KOKKOS_LAMBDA(const int num_face)
          {
            for (int dir=0; dir<dim; dir++)
              secmem(num_face,dir)*=rho_face_np1(num_face);
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
  tab_secmem.echange_espace_virtuel();
  statistics().end_count(STD_COUNTERS::ajouter_blocs);

}


void Navier_Stokes_Fluide_Dilatable_Proto::assembler_impl( Matrice_Morse& mat_morse, const DoubleTab& present, DoubleTab& secmem)
{
  Cerr << "Navier_Stokes_Fluide_Dilatable_Proto::assembler is not coded ! You should use assembler_avec_inertie !" << finl;
  Process::exit();
}

/*
 * ***************
 * Private methods
 * ***************
 */
void Navier_Stokes_Fluide_Dilatable_Proto::prepare_and_solve_u_star(Navier_Stokes_std& eqn,
                                                                    const Fluide_Dilatable_base& fluide_dil,
                                                                    DoubleTab& rhoU, DoubleTab& vpoint)
{
  const DoubleTab& tab_rho_face_n =fluide_dil.rho_face_n(), &tab_rho_face_np1=fluide_dil.rho_face_np1();
  const DoubleTab& tab_rho = fluide_dil.rho_discvit(); // rho with the same discretization as velocity
  const DoubleTab& vit = eqn.vitesse().valeurs();

  fluide_dil.secmembre_divU_Z(tab_W); // Compute W=-dZ/dt, right-hand side of the equation div(rhoU) = W
  vpoint=0;

  // add diffusion (with the dynamic viscosity)
  if (!eqn.schema_temps().diffusion_implicite()) eqn.operateur(0).ajouter(vpoint);

  DoubleTab& rhovitesse = ref_cast_non_const(DoubleTab,eqn.rho_la_vitesse().valeurs());
  rho_vitesse_impl(tab_rho,vit,rhovitesse);

  // add convection using rhovitesse
  if (!eqn.schema_temps().diffusion_implicite()) eqn.operateur(1).ajouter(rhovitesse,vpoint);
  else
    {
      DoubleTrav trav(vpoint);
      eqn.derivee_en_temps_conv(trav,rhovitesse);
      vpoint = trav;
    }

  // add source term
  eqn.sources().ajouter(vpoint);

  // add gradP
  eqn.corriger_derivee_expl(vpoint);

  const Champ_base& rho_vit=eqn.get_champ("rho_comme_v");
  ref_cast_non_const(DoubleTab,rho_vit.valeurs())=tab_rho_face_np1;

  if (eqn.schema_temps().diffusion_implicite())
    {
      DoubleTrav secmemV(vpoint);
      secmemV = vpoint;
      double dt = eqn.schema_temps().pas_de_temps();
      /*
       * secmemV contains M(rho^{n+1} u^{n+1} - rho^n u^n)/dt
       * M^{-1} secmemV*dt + rho^n u^n - rho^{n+1} u^n = rho^{n+1} (u^{n+1} - u^n)
       * dt/rho^{n+1} = (u^{n+1} - u^n)/dt
       * M^{-1} secmemV/rho^{n+1} + (rho^n - rho^{n+1})/rho^{n+1}/dt * u^n
       *
       * modify the mass solver to divide by rho^{n+1}
       * (also useful for implicit diffusion)
       */

      eqn.solv_masse().set_name_of_coefficient_temporel("rho_comme_v");
      eqn.solv_masse().appliquer(secmemV);
      DoubleTrav dr(tab_rho_face_n);

      CDoubleTabView tab_rho_face_n_v = tab_rho_face_n.view_ro();
      CDoubleTabView tab_rho_face_np1_v = tab_rho_face_np1.view_ro();
      DoubleTabView dr_v = dr.view_rw();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), dr.size_totale(),
                           KOKKOS_LAMBDA(
                             const int i)
      {
        dr_v(i, 0) = (tab_rho_face_n_v(i, 0) / tab_rho_face_np1_v(i, 0) - 1.) / dt;
      });
      end_gpu_timer(__KERNEL_NAME__);

      // use vpoint as a temporary array
      rho_vitesse_impl(dr,vit,vpoint);
      secmemV += vpoint;

      DoubleTrav delta_u(eqn.inconnue().futur());
      delta_u = eqn.inconnue().futur();
      eqn.Gradient_conjugue_diff_impl(secmemV, delta_u ) ;

      /*
       * delta_u=unp1 -un => delta_u + un=unp1
       * (delat_u + un)*rhonp1 = rhonp1 * unp1
       * (delat_u + un)*rhonp1  - rhon * un= rhonp1 * unp1 - rhon * un
       */

      delta_u *= dt;
      delta_u += vit;
      rho_vitesse_impl(tab_rho_face_np1,delta_u,vpoint);
      vpoint -= rhovitesse;
      vpoint /= dt;
      eqn.solv_masse().set_name_of_coefficient_temporel("no_coeff");
    }
  else eqn.solv_masse().appliquer(vpoint);

  update_vpoint_on_boundaries(eqn,fluide_dil,vpoint);

} /* END prepare_and_solve_u_star */

void Navier_Stokes_Fluide_Dilatable_Proto::update_vpoint_on_boundaries(const Navier_Stokes_std& eqn,
                                                                       const Fluide_Dilatable_base& fluide_dil,
                                                                       DoubleTab& tab_vpoint)
{
  // add d(rho)/dt at Dirichlet boundaries because the mass solver has set it to zero
  // NOTE: for incompressible flows the term is added by modifier_secmem
  const double dt_ = eqn.schema_temps().pas_de_temps();
  const DoubleTab& tab_rho_face_n = fluide_dil.rho_face_n(), &tab_rho_face_np1=fluide_dil.rho_face_np1();
  const DoubleTab& tab_vit = eqn.vitesse().valeurs();
  const Conds_lim& lescl = eqn.domaine_Cl_dis().les_conditions_limites();
  const IntTab& face_voisins = eqn.domaine_dis().face_voisins();
  const int taille = tab_vpoint.line_size();

  if (taille==1)
    if (orientation_VDF_.size() == 0)
      orientation_VDF_.ref(ref_cast(Domaine_VF,eqn.domaine_dis()).orientation());
  for (auto& itr : lescl)
    {
      const Cond_lim_base& la_cl_base = itr.valeur();
      if (sub_type(Dirichlet,la_cl_base))
        {
          const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl_base.frontiere_dis());
          const Dirichlet& diri=ref_cast(Dirichlet,la_cl_base);
          const int ndeb = la_front_dis.num_premiere_face(), nfin = ndeb + la_front_dis.nb_faces();

          if (taille==1) // VDF //
            {
              ToDo_Kokkos("critical");
              for (int num_face=ndeb; num_face<nfin; num_face++)
                {
                  int n0 = face_voisins(num_face, 0);
                  if (n0 == -1) n0 = face_voisins(num_face, 1);

                  // GF: in case of implicit diffusion, vpoint!=0 so we ignore the old value
                  tab_vpoint(num_face)=(diri.val_imp(num_face-ndeb,orientation_VDF_(num_face))*tab_rho_face_np1(num_face)-
                                        tab_vit(num_face)*tab_rho_face_n(num_face))/dt_;
                }
            }
          else // VEF //
            {
              int dim = Objet_U::dimension;
              CDoubleTabView val_imp = diri.tab_val_imp().view_ro();
              CDoubleArrView rho_face_np1 = static_cast<const ArrOfDouble&>(tab_rho_face_np1).view_ro();
              CDoubleArrView rho_face_n = static_cast<const ArrOfDouble&>(tab_rho_face_n).view_ro();
              CDoubleTabView vit = tab_vit.view_ro();
              DoubleTabView vpoint = tab_vpoint.view_wo();
              Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({ndeb, 0}, {nfin, dim});
              Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), policy, KOKKOS_LAMBDA(const int num_face, const int jj)
              {
                // GF: in case of implicit diffusion, vpoint!=0 so we ignore the old value
                vpoint(num_face,jj)=(rho_face_np1(num_face)*val_imp(num_face-ndeb,jj)
                                     -rho_face_n(num_face)*vit(num_face,jj))/dt_;
              });
              end_gpu_timer(__KERNEL_NAME__);
            }
        }
    }

} /* END update_vpoint_on_boundaries */

void Navier_Stokes_Fluide_Dilatable_Proto::solve_pressure_increment(Navier_Stokes_std& eqn,
                                                                    const Fluide_Dilatable_base& fluide_dil,
                                                                    DoubleTab& rhoU, DoubleTab& secmem,
                                                                    DoubleTab& inc_pre, DoubleTab& vpoint)
{
  const DoubleTab& tab_rho_face_n =fluide_dil.rho_face_n(), &tab_rho_face_np1=fluide_dil.rho_face_np1();
  const DoubleTab& vit = eqn.vitesse().valeurs();
  const double dt_ = eqn.schema_temps().pas_de_temps(), t = eqn.schema_temps().temps_courant();

  // Pressure resolution
  vpoint.echange_espace_virtuel();

  // Compute rhoU(n) :
  rho_vitesse_impl(tab_rho_face_n,vit,rhoU);

  // Add source term to vpoint if interfaces
  Probleme_base& prob=eqn.probleme();
  DoubleTrav vpoint0(vpoint);
  vpoint0 = vpoint;
  for (int i=0; i<prob.nombre_d_equations(); i++)
    if (sub_type(Transport_Interfaces_base,prob.equation(i)))
      {
        Transport_Interfaces_base& eq_transport = ref_cast(Transport_Interfaces_base,prob.equation(i));
        const int nb = vpoint.dimension(0), m = vpoint.line_size();
        DoubleTab source_ibc(nb,m);

        // Add a source term to vpoint to impose the interface velocity on the fluid.
        // source_ibc is local and not post-processable (unlike the FT case where the source term is defined and can be post-processed).
        eq_transport.modifier_vpoint_pour_imposer_vit(rhoU,vpoint0,vpoint,tab_rho_face_np1,source_ibc,t,dt_);
      }

  secmem = tab_W;
  operator_negate(secmem);
  eqn.operateur_divergence().ajouter(rhoU,secmem);
  secmem /= dt_; // (-tabW + Div(rhoU))/dt

  eqn.operateur_divergence().ajouter(vpoint, secmem);
  secmem *= -1;
  secmem.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_Fluide_Dilatable_base::derivee_en_temps_inco, secmem : ", secmem);

  // assembler is called only once during preparer_calcul (instead of assembler_QC)
  // Correction of the right-hand side according to boundary conditions:
  eqn.assembleur_pression()->modifier_secmem(secmem);
  eqn.solveur_pression().resoudre_systeme(eqn.matrice_pression().valeur(),secmem,inc_pre);

} /* END  solve_pressure_increment */

void Navier_Stokes_Fluide_Dilatable_Proto::correct_and_compute_u_np1(Navier_Stokes_std& eqn,
                                                                     const Fluide_Dilatable_base& fluide_dil,
                                                                     DoubleTab& rhoU,DoubleTab& Mmoins1grad,
                                                                     DoubleTab& inc_pre,DoubleTab& gradP,
                                                                     DoubleTab& vpoint)
{
  const DoubleTab& tab_rho_face_np1=fluide_dil.rho_face_np1();
  const DoubleTab& vit = eqn.vitesse().valeurs();
  DoubleTab& press = eqn.pression().valeurs();
  const double dt_ = eqn.schema_temps().pas_de_temps();

  // The virtual space of the pressure is needed to compute the gradient below,
  // and modifier_solution does not always perform the virtual space exchange.
  // We assume that pression and inc_pre have their virtual space up to date.
  // Compute pression += inc_pre:
  operator_add(press, inc_pre, VECT_ALL_ITEMS);
  eqn.assembleur_pression()->modifier_solution(press);

  // Pressure correction of velocity: M^{-1} B^T P
  eqn.solv_masse().appliquer(gradP);
  vpoint += gradP; // M-1 F

  press.echange_espace_virtuel();
  eqn.operateur_gradient().calculer(press, gradP);

  // Save B^T P for the next step.
  Mmoins1grad = gradP;
  eqn.solv_masse().appliquer(Mmoins1grad);

  // Pressure correction
  vpoint -= Mmoins1grad;

  // vpoint = (rhoU(n+1)-rhoU(n))/dt
  vpoint *= dt_;
  vpoint += rhoU; // rhoU(n+1)

  // Compute U(n+1):
  tab_divide_any_shape(vpoint, tab_rho_face_np1);

  // Compute (U(n+1)-U(n))/dt :
  vpoint -= vit;
  vpoint /= dt_;

  vpoint.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_Fluide_Dilatable_base::derivee_en_temps_inco, vpoint : ", vpoint);

} /* END correct_and_compute_u_np1 */

