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

#include <Op_Diff_VEF_Face.h>
#include <Champ_P1NC.h>
#include <Champ_Uniforme.h>
#include <Periodique.h>
#include <Symetrie.h>
#include <Neumann_homogene.h>
#include <Neumann_paroi.h>
#include <Echange_externe_impose.h>
#include <Echange_externe_radiatif.h>
#include <Neumann_sortie_libre.h>
#include <Milieu_base.h>
#include <TRUSTTrav.h>
#include <Probleme_base.h>
#include <Navier_Stokes_std.h>
#include <Porosites_champ.h>
#include <Device.h>
#include <Echange_couplage_thermique.h>
#include <Champ_front_calc_interne.h>
#include <Robin_VEF.h>

Implemente_instanciable_sans_constructeur(Op_Diff_VEF_Face,"Op_Diff_VEF_P1NC",Op_Diff_VEF_base);

Op_Diff_VEF_Face::Op_Diff_VEF_Face()
{
  declare_support_masse_volumique(1);
}

Sortie& Op_Diff_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Entree& Op_Diff_VEF_Face::readOn(Entree& s )
{
  return s ;
}

/*! @brief Associate the diffusivity field.
 *
 */
void Op_Diff_VEF_Face::associer_diffusivite(const Champ_base& diffu)
{
  diffusivite_ = diffu;
}

void Op_Diff_VEF_Face::completer()
{
  Operateur_base::completer();
}

const Champ_base& Op_Diff_VEF_Face::diffusivite() const
{
  return diffusivite_.valeur();
}

int ma_func_qui_renvoie_int()
{
  return 1;
}

void Op_Diff_VEF_Face::ajouter_cas_scalaire(const DoubleTab& tab_inconnue,
                                            DoubleTab& tab_resu, DoubleTab& tab_flux_bords,
                                            DoubleTab& tab_nu,
                                            const Domaine_Cl_VEF& domaine_Cl_VEF,
                                            const Domaine_VEF& domaine_VEF ) const
{
  int nb_faces = domaine_VEF.nb_faces();
  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  int nb_bords=domaine_VEF.nb_front_Cl();
  const int premiere_face_int=domaine_VEF.premiere_face_int();

  {
    CIntTabView elem_faces = domaine_VEF.elem_faces().view_ro();
    CIntTabView face_voisins = domaine_VEF.face_voisins().view_ro();
    CDoubleTabView face_normale = domaine_VEF.face_normales().view_ro();
    CDoubleArrView inverse_volumes = domaine_VEF.inverse_volumes().view_ro();
    CDoubleTabView nu = tab_nu.view_ro();
    CDoubleTabView inconnue = tab_inconnue.view_ro();
    DoubleArrView flux_bords = static_cast<ArrOfDouble&>(tab_flux_bords).view_rw();
    DoubleArrView resu = static_cast<ArrOfDouble&>(tab_resu).view_rw();
    // Process boundary faces
    for (int n_bord = 0; n_bord < nb_bords; n_bord++)
      {
        const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
        const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
        int num1 = 0;
        int num2 = le_bord.nb_faces_tot();
        int nb_faces_bord_reel = le_bord.nb_faces();
        CIntArrView le_bord_num_face = le_bord.num_face().view_ro();
        if (sub_type(Periodique, la_cl.valeur()))
          {
            const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
            CIntArrView face_associee = la_cl_perio.face_associee().view_ro();
            Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__),
                                 Kokkos::RangePolicy<>(num1, nb_faces_bord_reel), KOKKOS_LAMBDA(
                                   const int ind_face)
            {
              int num_face = le_bord_num_face(ind_face);
              int fac_asso = face_associee(ind_face);
              fac_asso = le_bord_num_face(fac_asso);
              for (int kk = 0; kk < 2; kk++)
                {
                  int elem = face_voisins(num_face, kk);
                  for (int i = 0; i < nb_faces_elem; i++)
                    {
                      int j = elem_faces(elem, i);
                      if (j > num_face && j != fac_asso)
                        {
                          double valA = viscA(num_face, j, elem, nu(elem, 0), face_voisins, face_normale,
                                              inverse_volumes);
                          double flux = valA * (inconnue(j, 0) - inconnue(num_face, 0));
                          Kokkos::atomic_add(&resu(num_face), +flux);
                          if (j < nb_faces) // face reelle
                            Kokkos::atomic_add(&resu(j), -0.5 * flux);
                        }
                    }
                }
            });
            end_gpu_timer(__KERNEL_NAME__);
          }
        else     // There is only one component, so we handle
          // a scalar equation (not velocity): no need to use
          // the tangential tau (wall thermal laws do not compute
          // turbulent exchange at the wall for now)
          {
            Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__),
                                 Kokkos::RangePolicy<>(num1, num2), KOKKOS_LAMBDA(
                                   const int ind_face)
            {
              int num_face = le_bord_num_face(ind_face);
              int elem = face_voisins(num_face, 0);
              for (int i = 0; i < nb_faces_elem; i++)
                {
                  int j = elem_faces(elem, i);
                  if (j > num_face || num_face >= nb_faces)
                    {
                      double valA = viscA(num_face, j, elem, nu(elem, 0), face_voisins, face_normale,
                                          inverse_volumes);
                      double flux = valA * (inconnue(j, 0) - inconnue(num_face, 0));
                      if (num_face < nb_faces) // face reelle
                        {
                          Kokkos::atomic_add(&resu(num_face), +flux);
                          Kokkos::atomic_add(&flux_bords(num_face), -flux);
                        }
                      if (j < nb_faces) // face reelle
                        {
                          Kokkos::atomic_add(&resu(j), -flux);
                          if (j < premiere_face_int)
                            Kokkos::atomic_add(&flux_bords(j), +flux);
                        }
                    }
                }
            });
            end_gpu_timer(__KERNEL_NAME__);
          }
      }

    // Internal faces:
    Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__),
                         Kokkos::MDRangePolicy<Kokkos::Rank<2>>({premiere_face_int, 0}, {nb_faces, 2}),
                         KOKKOS_LAMBDA(const int num_face, const int k)
    {
      int elem = face_voisins(num_face, k);
      for (int i = 0; i < nb_faces_elem; i++)
        {
          int j = elem_faces(elem, i);
          if (j > num_face)
            {
              int contrib = 1;

              if (j >= nb_faces) // This is a virtual face
                {
                  int el1 = face_voisins(j, 0);
                  int el2 = face_voisins(j, 1);
                  if ((el1 == -1) || (el2 == -1))
                    contrib = 0;
                }

              if (contrib)
                {
                  double valA = viscA(num_face, j, elem, nu(elem, 0), face_voisins,
                                      face_normale, inverse_volumes);
                  double flux = valA * (inconnue(j, 0) - inconnue(num_face, 0));
                  Kokkos::atomic_add(&resu(num_face), flux);
                  if (j < nb_faces) // Process real faces
                    Kokkos::atomic_add(&resu(j), -flux);
                }
            }
        }
    });
    end_gpu_timer(__KERNEL_NAME__);
  }

  // Neumann :
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();
      if (sub_type(Neumann_paroi,la_cl.valeur()))
        {
          const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
          CDoubleArrView surface = domaine_VEF.face_surfaces().view_ro();
          CDoubleTabView flux_impose = la_cl_paroi.flux_impose().view_ro();
          DoubleArrView flux_bords = static_cast<ArrOfDouble&>(tab_flux_bords).view_rw();
          DoubleArrView resu = static_cast<ArrOfDouble&>(tab_resu).view_rw();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(ndeb, nfin), KOKKOS_LAMBDA(const int face)
          {
            double flux = flux_impose(face-ndeb, 0) * surface(face);
            resu(face) += flux;
            flux_bords(face) = flux;
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else if (sub_type(Echange_externe_impose,la_cl.valeur()))
        {
          const Echange_externe_impose& la_cl_paroi = ref_cast(Echange_externe_impose, la_cl.valeur());
          const double coeff = COEFF_STEFAN_BOLTZMANN;
          const bool has_emissivity = la_cl_paroi.has_emissivite();
          CDoubleArrView surface = domaine_VEF.face_surfaces().view_ro();
          CDoubleArrView text = static_cast<const ArrOfDouble&>(la_cl_paroi.tab_T_ext()).view_ro();
          CDoubleArrView himp = static_cast<const ArrOfDouble&>(la_cl_paroi.tab_h_imp()).view_ro();
          CDoubleArrView eps;
          if (has_emissivity) eps = static_cast<const ArrOfDouble&>(la_cl_paroi.tab_emissivite()).view_ro();
          CDoubleArrView inconnue = static_cast<const ArrOfDouble&>(tab_inconnue).view_ro();
          DoubleArrView resu = static_cast<ArrOfDouble&>(tab_resu).view_rw();
          DoubleArrView flux_bords = static_cast<ArrOfDouble&>(tab_flux_bords).view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(ndeb, nfin), KOKKOS_LAMBDA (const int face)
          {
            int ind_face = face - ndeb;
            double flux = himp(ind_face)*(text(ind_face)-inconnue(face))*surface(face);
            resu[face] += flux;
            flux_bords(face) = flux;

            if (has_emissivity)
              {
                double T = inconnue(face);
                double t_ext = text(ind_face);
                flux = coeff * eps(ind_face) * (t_ext * t_ext * t_ext * t_ext - T * T * T * T) * surface(face);
                resu[face] += flux;
                flux_bords(face) += flux;
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else if (sub_type(Echange_couplage_thermique, la_cl.valeur()))
        {
          ToDo_Kokkos("critical");
          const Echange_couplage_thermique& la_cl_paroi = ref_cast(Echange_couplage_thermique, la_cl.valeur());
          const DoubleVect& surface = domaine_VEF.face_surfaces();
          for (int face=ndeb; face<nfin; face++)
            {
              double h=la_cl_paroi.h_imp(face-ndeb);
              double Text=la_cl_paroi.T_ext(face-ndeb);
              double phiext=la_cl_paroi.flux_exterieur_impose(face-ndeb);
              double flux=(phiext+h*(Text-tab_inconnue(face)))*surface(face);
              tab_resu[face] += flux;
              tab_flux_bords(face) = flux;
            }
        }
      else if (sub_type(Neumann_homogene,la_cl.valeur())
               || sub_type(Symetrie,la_cl.valeur())
               || sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          DoubleArrView flux_bords = static_cast<ArrOfDouble&>(tab_flux_bords).view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(ndeb, nfin), KOKKOS_LAMBDA(const int face)
          {
            flux_bords(face) = 0.;
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
}

void Op_Diff_VEF_Face::ajouter_cas_vectoriel(const DoubleTab& inconnue,
                                             DoubleTab& resu, DoubleTab& tab_flux_bords,
                                             DoubleTab& nu,
                                             const Domaine_Cl_VEF& domaine_Cl_VEF,
                                             const Domaine_VEF& domaine_VEF,
                                             int nb_comp) const
{
  assert(nb_comp==dimension);

  // Build grad_ array if necessary
  if(!grad_.get_md_vector())
    {
      grad_.resize(0, Objet_U::dimension, Objet_U::dimension);
      domaine_VEF.domaine().creer_tableau_elements(grad_);
    }
  Champ_P1NC::calcul_gradient(inconnue,grad_,domaine_Cl_VEF);

  /* ToDo OpenMP: refactor with Op_Dift_VEF_Face.cpp into a template class
  if (le_modele_turbulence->utiliser_loi_paroi())
   {
      Champ_P1NC::calcul_duidxj_paroi(grad_,nu,nu_turb,tau_tan_,domaine_Cl_VEF);
      grad_.echange_espace_virtuel(); // gradient_elem up to date on virtual elements
  }
  DoubleTab Re;
  Re.resize(0, Objet_U::dimension, Objet_U::dimension);
  domaine_VEF.domaine().creer_tableau_elements(Re);
  Re = 0.;
  if (le_modele_turbulence->calcul_tenseur_Re(nu_turb, grad_, Re))
  {
      Cerr << "Using nonlinear turbulent diffusion in NS" << finl;
      for (int elem=0; elem<nb_elem; elem++)
          for (int i=0; i<nbr_comp; i++)
              for (int j=0; j<nbr_comp; j++)
                  Re(elem,i,j) *= nu_turb[elem];
  }
  else
  {
      for (int elem=0; elem<nb_elem; elem++)
          for (int i=0; i<nbr_comp; i++)
              for (int j=0; j<nbr_comp; j++)
                  Re(elem,i,j) = nu_turb[elem]*(grad_(elem,i,j) + grad_(elem,j,i));
  }
  Re.echange_espace_virtuel();
  */

  int nb_faces = domaine_VEF.nb_faces();
  int nb_faces_bord = domaine_VEF.premiere_face_int();
  CIntTabView face_voisins_v = domaine_VEF.face_voisins().view_ro();
  CDoubleTabView face_normales_v = domaine_VEF.face_normales().view_ro();
  CDoubleTabView nu_v = nu.view_ro();
  CDoubleTabView3 grad_v = grad_.view_ro<3>();
  DoubleTabView resu_v = resu.view_rw();
  DoubleTabView tab_flux_bords_v = tab_flux_bords.view_rw();

  auto kern_ajouter = KOKKOS_LAMBDA(int
                                    num_face, int k)
  {
    int elem = face_voisins_v(num_face, k);
    if (elem >= 0)
      {
        int ori = 1 - 2 * k;
        double nu_elem = nu_v(elem, 0);
        for (int i = 0; i < nb_comp; i++)
          for (int j = 0; j < nb_comp; j++)
            {
              double grad_ij = grad_v(elem, i, j);
              double grad_ji = grad_v(elem, j, i);
              double fn = face_normales_v(num_face, j);
              double flux = ori * fn * (nu_elem * grad_ij  /* + Re(elem, i, j) */ );
              Kokkos::atomic_sub(&resu_v(num_face, i), flux);

              if (num_face < nb_faces_bord)
                {
                  double flux_bord = ori * fn * (nu_elem * (grad_ij + grad_ji));
                  Kokkos::atomic_sub(&tab_flux_bords_v(num_face, i), flux_bord);
                }
            }
      }
  };
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {nb_faces,2}) , kern_ajouter);
  end_gpu_timer(__KERNEL_NAME__);


  const int nb_bords=domaine_VEF.nb_front_Cl();
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      // Update flux_bords on symmetry:
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Symetrie,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          DoubleTabView flux_bords = tab_flux_bords.view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(ndeb, nfin), KOKKOS_LAMBDA(const int face)
          {
            flux_bords(face, 0) = 0.;
          });
          end_gpu_timer(__KERNEL_NAME__);
        }

      else if (sub_type(Robin_VEF, la_cl.valeur()))
        {
#ifdef TRUST_USE_GPU
          Cerr << "Warning not tested on GPU" << finl;
#endif
          ToDo_Kokkos("critical");

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          const Robin_VEF& la_cl_robin = ref_cast(Robin_VEF,la_cl.valeur());
          int marq = phi_psi_diffuse(equation());
          const DoubleVect& porosite_face = equation().milieu().porosite_face();
          double scale_factor_is_one = 1.;
          double inv_alpha = 1./la_cl_robin.get_alpha_cl() ;
          double inv_beta =  1./la_cl_robin.get_beta_cl();
          double inv_alpha_minus_inv_beta = 1./la_cl_robin.get_alpha_cl() - 1./la_cl_robin.get_beta_cl();
          DoubleTab normal_vector;
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb +le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              int id_face_bord = face -ndeb;
              double face_surface = domaine_VEF.face_surfaces(face);
              normal_vector = domaine_VEF.normalized_boundaries_outward_vector(face, scale_factor_is_one);

              for (int nc1=0; nc1<nb_comp; nc1++)
                {
                  double flux_robin_uu = 0. ;
                  double flux_robin_rhs = 0.;
                  double flux_tot ;

                  // forward term for velocity
                  for (int nc2 = 0; nc2<nb_comp; nc2++)
                    {
                      const double normal2 = normal_vector(nc1)*normal_vector(nc2);
                      flux_robin_uu += (inv_beta*(nc1==nc2) + inv_alpha_minus_inv_beta*normal2)* (face_surface);
                    }
                  flux_robin_uu *= inconnue(face,nc1);

                  // rhs for robin bc
                  double val;
                  if (dimension == 2)
                    {
                      // add normal component rhs
                      val = inv_alpha * normal_vector(nc1) * la_cl_robin.flux_normal_imp(id_face_bord);

                      // add tangential component rhs
                      double tgte = (2*nc1-1)*normal_vector(1-nc1);
                      val += inv_beta * la_cl_robin.flux_tangentiel_imp(id_face_bord, 0)*tgte;
                    }
                  else
                    {
                      // add normal component rhs
                      val = inv_alpha * normal_vector(nc1) * la_cl_robin.flux_normal_imp(id_face_bord);

                      // add tangential component rhs
                      val += inv_beta * la_cl_robin.flux_tangentiel_imp(id_face_bord, nc1) ;
                    }
                  flux_robin_rhs = val*face_surface;
                  flux_tot = (flux_robin_rhs - flux_robin_uu)* (marq ? porosite_face(face) : 1);
                  resu(face,nc1) +=  flux_tot;
                  tab_flux_bords(face,nc1) +=  flux_tot;


                }
            }
        }
    }
}

void Op_Diff_VEF_Face::ajouter_cas_multi_scalaire(const DoubleTab& inconnue,
                                                  DoubleTab& resu, DoubleTab& tab_flux_bords,
                                                  DoubleTab& nu,
                                                  const Domaine_Cl_VEF& domaine_Cl_VEF,
                                                  const Domaine_VEF& domaine_VEF,
                                                  int nb_comp) const
{
  ToDo_Kokkos("critical");
  const IntTab& elemfaces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  int i0,j,num_face;
  int nb_faces = domaine_VEF.nb_faces();
  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  int n_bord;
  double flux0;
  //DoubleVect n(Objet_U::dimension);
  //DoubleTrav Tgrad(Objet_U::dimension,Objet_U::dimension);

  assert(nb_comp>1);
  int nb_bords=domaine_VEF.nb_front_Cl();
  int ind_face;

  for (n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      // const IntTab& elemfaces = domaine_VEF.elem_faces();
      int num1 = 0;
      int num2 = le_bord.nb_faces_tot();
      int nb_faces_bord_reel = le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          int fac_asso;
          for (ind_face=num1; ind_face<nb_faces_bord_reel; ind_face++)
            {
              fac_asso = la_cl_perio.face_associee(ind_face);
              fac_asso = le_bord.num_face(fac_asso);
              num_face = le_bord.num_face(ind_face);
              for (int kk=0; kk<2; kk++)
                {
                  int elem = face_voisins(num_face, kk);
                  for (i0=0; i0<nb_faces_elem; i0++)
                    {
                      if ( ( (j= elemfaces(elem,i0)) > num_face ) && (j != fac_asso ) )
                        {
                          for (int nc=0; nc<nb_comp; nc++)
                            {
                              double valA = viscA(num_face,j,elem,nu(elem,nc));
                              resu(num_face,nc)+=valA*inconnue(j,nc);
                              resu(num_face,nc)-=valA*inconnue(num_face,nc);
                              if(j<nb_faces) // face reelle
                                {
                                  ////WARNING: NUM_face differs from the reference version
                                  resu(j,nc)+=0.5*valA*inconnue(num_face,nc);
                                  resu(j,nc)-=0.5*valA*inconnue(j,nc);
                                }
                            }
                        }
                    }
                }
            }
        }// end if periodic
      else
        {
          for (ind_face=num1; ind_face<num2; ind_face++)
            {
              num_face = le_bord.num_face(ind_face);
              int elem=face_voisins(num_face,0);

              // Loop over faces:
              for (int i=0; i<nb_faces_elem; i++)
                if (( (j= elemfaces(elem,i)) > num_face ) || (ind_face>=nb_faces_bord_reel))
                  {
                    for (int nc=0; nc<nb_comp; nc++)
                      {
                        double valA = viscA(num_face,j,elem,nu(elem,nc));
                        if (ind_face<nb_faces_bord_reel)
                          {
                            double flux=valA*(inconnue(j,nc)-inconnue(num_face,nc));
                            resu(num_face,nc)+=flux;
                            tab_flux_bords(num_face,nc)-=flux;
                          }

                        if(j<nb_faces) // face reelle
                          {
                            resu(j,nc)+=valA*inconnue(num_face,nc);
                            resu(j,nc)-=valA*inconnue(j,nc);
                          }
                      }
                  }
            }
        }
    }// End for n_bord

  // Process internal faces

  for (num_face=domaine_VEF.premiere_face_int(); num_face<nb_faces; num_face++)
    {
      for (int k=0; k<2; k++)
        {
          int elem = face_voisins(num_face,k);
          for (i0=0; i0<nb_faces_elem; i0++)
            {
              if ( (j= elemfaces(elem,i0)) > num_face )
                {
                  int el1,el2;
                  int contrib=1;
                  if(j>=nb_faces) // This is a virtual face
                    {
                      el1 = face_voisins(j,0);
                      el2 = face_voisins(j,1);
                      if((el1==-1)||(el2==-1))
                        contrib=0;
                    }
                  if(contrib)
                    {
                      for (int nc=0; nc<nb_comp; nc++)
                        {
                          double valA = viscA(num_face,j,elem,nu(elem,nc));
                          resu(num_face,nc)+=valA*inconnue(j,nc);
                          resu(num_face,nc)-=valA*inconnue(num_face,nc);
                          if(j<nb_faces) // Process real faces only
                            {
                              resu(j,nc)+=valA*inconnue(num_face,nc);
                              resu(j,nc)-=valA*inconnue(j,nc);
                            }
                          else
                            {
                              // Face j is virtual
                            }
                        }
                    }
                }
            }
        }
    }// End internal faces


  //Based on what is done for the scalar case
  for (n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_paroi,la_cl.valeur()))
        {
          const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              for (int nc=0; nc<nb_comp; nc++)
                {
                  flux0=la_cl_paroi.flux_impose(face-ndeb,nc)*domaine_VEF.surface(face);
                  resu(face,nc) += flux0;
                  tab_flux_bords(face,nc) = flux0;
                }
            }
        }
      else if (sub_type(Echange_externe_impose,la_cl.valeur()))
        {
          throw;
          const Echange_externe_impose& la_cl_paroi = ref_cast(Echange_externe_impose, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              for (int nc=0; nc<nb_comp; nc++)
                {
                  flux0=la_cl_paroi.h_imp(face-ndeb,nc)*(la_cl_paroi.T_ext(face-ndeb,nc)-inconnue(face,nc))*domaine_VEF.surface(face);
                  resu(face,nc) += flux0;
                  tab_flux_bords(face,nc) = flux0;

                  if (la_cl_paroi.has_emissivite())
                    {
                      const double text = la_cl_paroi.T_ext(face - ndeb, nc), T = inconnue(face, nc);
                      flux0 = COEFF_STEFAN_BOLTZMANN * la_cl_paroi.emissivite(face - ndeb, nc) * (text * text * text * text - T * T * T * T) * domaine_VEF.face_surfaces(face);
                      resu(face, nc) += flux0;
                      tab_flux_bords(face, nc) += flux0;
                    }
                }
            }
        }
      else if (sub_type(Echange_couplage_thermique,la_cl.valeur()))
        {
          const Echange_couplage_thermique& la_cl_paroi = ref_cast(Echange_couplage_thermique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              for (int nc=0; nc<nb_comp; nc++)
                {
                  double phiext = la_cl_paroi.flux_exterieur_impose(face-ndeb,nc);
                  flux0 = (phiext + la_cl_paroi.h_imp(face-ndeb,nc)*(la_cl_paroi.T_ext(face-ndeb,nc)-inconnue(face,nc)))*domaine_VEF.surface(face);
                  resu(face,nc) += flux0;
                  tab_flux_bords(face,nc) = flux0;
                }
            }
        }
      else if (sub_type(Neumann_homogene,la_cl.valeur())
               || sub_type(Symetrie,la_cl.valeur())
               || sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            for (int nc=0; nc<nb_comp; nc++)
              tab_flux_bords(face,nc) = 0.;
        }
    }
}


DoubleTab& Op_Diff_VEF_Face::ajouter(const DoubleTab& inconnue_org, DoubleTab& resu) const
{
  remplir_nu(nu_);
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();

  int nb_comp = 1;
  int nb_dim = resu.nb_dim();
  if(nb_dim==2)
    nb_comp=resu.dimension(1);
  DoubleTab nu;
  DoubleTab tab_inconnue;
  int marq=phi_psi_diffuse(equation());
  const DoubleVect& porosite_face = equation().milieu().porosite_face();
  const DoubleVect& porosite_elem = equation().milieu().porosite_elem();
  // either div(phi nu grad inco)
  // or div(nu grad phi inco)
  // depending on whether phi_psi or psi is diffused
  modif_par_porosite_si_flag(nu_,nu,!marq,porosite_elem);
  const DoubleTab& inconnue=modif_par_porosite_si_flag(inconnue_org,tab_inconnue,marq,porosite_face);

  const Champ_base& inco = equation().inconnue();
  const Nature_du_champ nature_champ = inco.nature_du_champ();

  // Size and initialize the flux balance array:
  if (flux_bords_.size_array()!=domaine_VEF.nb_faces_bord()) flux_bords_.resize(domaine_VEF.nb_faces_bord(),nature_champ==scalaire ? 1 : nb_comp);
  flux_bords_=0.;

  if(nature_champ==scalaire)
    ajouter_cas_scalaire(inconnue, resu, flux_bords_, nu, domaine_Cl_VEF, domaine_VEF);
  else if (nature_champ==vectoriel)
    ajouter_cas_vectoriel(inconnue, resu, flux_bords_, nu, domaine_Cl_VEF, domaine_VEF,nb_comp);
  else if (nature_champ==multi_scalaire)
    ajouter_cas_multi_scalaire(inconnue, resu, flux_bords_, nu, domaine_Cl_VEF, domaine_VEF,nb_comp);
  modifier_flux(*this);

  return resu;
}

DoubleTab& Op_Diff_VEF_Face::calculer(const DoubleTab& inconnue, DoubleTab& resu) const
{
  resu = 0;
  return ajouter(inconnue,resu);
}

void Op_Diff_VEF_Face::ajouter_contribution(const DoubleTab& tab_transporte, Matrice_Morse& tab_matrice) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();

  modifier_matrice_pour_periodique_avant_contribuer(tab_matrice,equation());

  // Fill the nu array because matrix assembly with ajouter_contribution
  // may be performed before the first time step
  remplir_nu(nu_);
  DoubleTrav tab_nu;

  // either div(phi nu grad inco)
  // or div(nu grad phi inco)
  // depending on whether phi_psi or psi is diffused
  int marq = phi_psi_diffuse(equation());
  modif_par_porosite_si_flag(nu_,tab_nu,!marq,equation().milieu().porosite_elem());

  int nb_dim = tab_transporte.nb_dim();
  int nb_comp = (nb_dim==2 ? tab_transporte.dimension(1) : 1);
  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  int nb_bords = domaine_VEF.nb_front_Cl();


  IntTrav tab_face_associee(domaine_VEF.premiere_face_int());
  IntTrav tab_fac2b_idx(domaine_VEF.nb_faces());


  // Retrieve the indices of periodic boundary faces and the
  // associated faces into arrays beforehand to
  // allow structured access in the kernel afterwards
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int num1 = le_bord.num_premiere_face();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());

          int nb_faces = le_bord.nb_faces();
          int num2b = num1 + nb_faces / 2;

          // only iterate over half the periodic faces
          // the result will be copied to the associated face at the end...
          ToDo_Kokkos("critical");
          for (int fac = num1; fac < num2b; fac++)
            {
              int fac_asso = la_cl_perio.face_associee(fac - num1) + num1;

              tab_face_associee(fac) = fac_asso;
              tab_fac2b_idx(fac) = fac;
              tab_fac2b_idx(fac+nb_faces/2) = fac;
            }
        }
    }

  CIntTabView elem_faces = domaine_VEF.elem_faces().view_ro();
  CIntTabView face_voisins = domaine_VEF.face_voisins().view_ro();
  CDoubleArrView porosite_face = equation().milieu().porosite_face().view_ro();
  CDoubleArrView inverse_volumes = domaine_VEF.inverse_volumes().view_ro();
  CDoubleTabView face_normale = domaine_VEF.face_normales().view_ro();
  CDoubleTabView nu = tab_nu.view_ro();
  CIntArrView est_face_bord = domaine_VEF.est_face_bord().view_ro();
  CIntArrView face_associee = static_cast<ArrOfInt&>(tab_face_associee).view_ro();
  CIntArrView fac2b_idx = static_cast<ArrOfInt&>(tab_fac2b_idx).view_ro();
  Matrice_Morse_View matrice;
  matrice.set(tab_matrice);

  auto ajouter_contrib = KOKKOS_LAMBDA(const int fac)
  {
    int type_face = est_face_bord(fac);
    int fac2b = fac2b_idx(fac);

    if (type_face == 2 && fac == fac2b) // faces perio
      {
        int elem1 = face_voisins(fac,0);
        int elem2 = face_voisins(fac,1);

        int fac_asso = face_associee(fac);
        for (int i = 0; i < nb_faces_elem; i++)
          {
            int j = elem_faces(elem1,i);
            if (j > fac)
              {
                double val = viscA(fac,j,elem1,nu(elem1,0), face_voisins, face_normale, inverse_volumes);
                double coeff_face1 = val * (marq ? porosite_face(fac) : 1);
                double coeff_face2 = val * (marq ? porosite_face(j) : 1);

                for (int nc = 0; nc < nb_comp; nc++)
                  {
                    int n0 = fac*nb_comp + nc;
                    int j0 = j*nb_comp + nc;

                    matrice.atomic_add(n0, n0, +coeff_face1);
                    matrice.atomic_add(n0, j0, -coeff_face2);
                    matrice.atomic_add(j0, n0, -coeff_face1);
                    matrice.atomic_add(j0, j0, +coeff_face2);
                  }
              }
            if (elem2 != -1)
              {
                j = elem_faces(elem2,i);
                if (j > fac)
                  {
                    double val = viscA(fac,j,elem2,nu(elem2,0), face_voisins, face_normale, inverse_volumes);
                    double coeff_face1 = val * (marq ? porosite_face(fac) : 1);
                    double coeff_face2 = val * (marq ? porosite_face(j) : 1);

                    for (int nc = 0; nc < nb_comp; nc++)
                      {
                        int n0 = fac*nb_comp + nc;
                        int j0 = j*nb_comp + nc;
                        int n1 = fac_asso*nb_comp+nc;

                        matrice.atomic_add(n0, n0, +coeff_face1);
                        matrice.atomic_add(n0, j0, -coeff_face2);
                        matrice.atomic_add(j0, n1, -coeff_face1);
                        matrice.atomic_add(j0, j0, +coeff_face2);
                      }
                  }
              }
          }
      }
    else if (type_face == 1) // non-periodic boundary faces
      {
        int elem1 = face_voisins(fac,0);
        for (int i = 0; i < nb_faces_elem; i++)
          {
            int j = elem_faces(elem1,i);
            if (j > fac)
              {
                double val = viscA(fac,j,elem1,nu(elem1,0), face_voisins, face_normale, inverse_volumes);
                double coeff_face1 = val * (marq ? porosite_face(fac) : 1);
                double coeff_face2 = val * (marq ? porosite_face(j) : 1);

                for (int nc = 0; nc < nb_comp; nc++)
                  {
                    int n0 = fac*nb_comp + nc;
                    int j0 = j*nb_comp + nc;

                    matrice.atomic_add(n0, n0, +coeff_face1);
                    matrice.atomic_add(n0, j0, -coeff_face2);
                    matrice.atomic_add(j0, n0, -coeff_face1);
                    matrice.atomic_add(j0, j0, +coeff_face2);
                  }
              }
          }
      }
    else if (type_face == 0) // faces internes
      {
        for (int k=0; k<2; k++)
          {
            int elem = face_voisins(fac,k);
            if (elem!=-1)
              {
                for (int i = 0; i < nb_faces_elem; i++)
                  {
                    int j = elem_faces(elem, i);
                    if (j > fac)
                      {
                        double val = viscA(fac, j, elem, nu(elem, 0), face_voisins, face_normale, inverse_volumes);
                        double coeff_face1 = val * (marq ? porosite_face(fac) : 1);
                        double coeff_face2 = val * (marq ? porosite_face(j) : 1);

                        for (int nc = 0; nc < nb_comp; nc++)
                          {
                            int n0 = fac * nb_comp + nc;
                            int j0 = j * nb_comp + nc;

                            matrice.atomic_add(n0, n0, +coeff_face1);
                            matrice.atomic_add(n0, j0, -coeff_face2);
                            matrice.atomic_add(j0, n0, -coeff_face1);
                            matrice.atomic_add(j0, j0, +coeff_face2);
                          }
                      }
                  }
              }
          }
      }
  };
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__),domaine_VEF.nb_faces(), ajouter_contrib);
  end_gpu_timer(__KERNEL_NAME__);

  int premiere_face_int = domaine_VEF.premiere_face_int();
  DoubleTrav tab_h_impose(premiere_face_int);
  DoubleTrav tab_derivee_flux_exterieur_imposee(premiere_face_int);

  // Neumann: fill arrays with boundary conditions for the Kokkos kernel
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Echange_externe_impose,la_cl.valeur()))
        {
          const Echange_externe_impose& la_cl_paroi = ref_cast(Echange_externe_impose, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          const double coeff = COEFF_STEFAN_BOLTZMANN;
          const bool has_emissivity = la_cl_paroi.has_emissivite();
          CDoubleArrView inconnue = static_cast<const ArrOfDouble&>(equation().inconnue().valeurs()).view_ro();
          CDoubleArrView himp = static_cast<const ArrOfDouble&>(la_cl_paroi.tab_h_imp()).view_ro();
          CDoubleArrView eps;
          if (has_emissivity) eps = static_cast<const ArrOfDouble&>(la_cl_paroi.tab_emissivite()).view_ro();
          DoubleArrView h_impose = static_cast<ArrOfDouble&>(tab_h_impose).view_wo();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), Kokkos::RangePolicy<>(ndeb, nfin), KOKKOS_LAMBDA(const int face)
          {
            int ind_face = face - ndeb;
            h_impose(face) = himp(ind_face);
            if (has_emissivity)
              {
                const double T = inconnue(face);
                h_impose(face) = 4 * coeff * eps(ind_face) * T * T * T;
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else if (sub_type(Echange_couplage_thermique, la_cl.valeur()))
        {
          const Echange_couplage_thermique& la_cl_paroi = ref_cast(Echange_couplage_thermique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());

          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          ToDo_Kokkos("critical");
          for (int face = ndeb; face < nfin; face++)
            {
              tab_h_impose(face) = la_cl_paroi.h_imp(face-ndeb);
              tab_derivee_flux_exterieur_imposee(face) = la_cl_paroi.derivee_flux_exterieur_imposee(face-ndeb);
            }
        }
      else if (sub_type(Robin_VEF, la_cl.valeur()))
        {
          ToDo_Kokkos("critical");
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          const Robin_VEF& la_cl_robin = ref_cast(Robin_VEF,la_cl.valeur());
          double inv_alpha_minus_inv_beta = 1./la_cl_robin.get_alpha_cl() - 1./la_cl_robin.get_beta_cl();
          double inv_beta  = 1./la_cl_robin.get_beta_cl();
          double scale_factor_is_one = 1. ;
          DoubleTab normal_vector ;
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          const DoubleVect& tab_porosite_face = equation().milieu().porosite_face();
          for (int face = ndeb; face < nfin; face++)
            {
              double face_surface = domaine_VEF.face_surfaces(face);
              normal_vector = domaine_VEF.normalized_boundaries_outward_vector(face, scale_factor_is_one);
              //int elem  = face_voisins(face, 0) ;
              for (int nc1 = 0; nc1 < nb_comp; nc1++)
                {
                  const int i = face * nb_comp + nc1;

                  // diagonal term
                  double val = (inv_beta + inv_alpha_minus_inv_beta*normal_vector(nc1)*normal_vector(nc1))*face_surface;
                  double robin_contribution=  val * (marq ? tab_porosite_face(face) : 1) ;
                  tab_matrice(i,i) += robin_contribution  ;

                  // extradiagonal term
                  for (int nc2 = 0; nc2<nc1; nc2++)
                    {
                      const int j = face * nb_comp + nc2;
                      const double normal2 = normal_vector(nc1)*normal_vector(nc2);
                      val = inv_alpha_minus_inv_beta*normal2* (face_surface);
                      robin_contribution=  val * (marq ? tab_porosite_face(face) : 1) ;
                      tab_matrice(i, j) += robin_contribution;
                      tab_matrice(j ,i) += robin_contribution;
                    }
                }
            }
        }
    }

  CDoubleArrView h_impose = static_cast<const ArrOfDouble&>(tab_h_impose).view_ro();
  CDoubleArrView derivee_flux_exterieur_imposee = static_cast<const ArrOfDouble&>(tab_derivee_flux_exterieur_imposee).view_ro();
  CDoubleArrView face_surfaces = domaine_VEF.face_surfaces().view_ro();

  // Neumann: compute contributions on boundaries
  auto neumann = KOKKOS_LAMBDA (const int face)
  {
    double h = h_impose(face);
    double dphi_dT = derivee_flux_exterieur_imposee(face);
    matrice.add(face,face, + (h + dphi_dT) * face_surfaces(face));
  };

  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), premiere_face_int, neumann);
  end_gpu_timer(__KERNEL_NAME__);

  modifier_matrice_pour_periodique_apres_contribuer(tab_matrice,equation());
}

void Op_Diff_VEF_Face::ajouter_contribution_multi_scalaire(const DoubleTab& tab_transporte, Matrice_Morse& tab_matrice) const
{
  modifier_matrice_pour_periodique_avant_contribuer(tab_matrice, equation());

  // Fill the nu array because matrix assembly with ajouter_contribution
  // may be performed before the first time step
  remplir_nu(nu_);

  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  const IntTab& tab_elem_faces = domaine_VEF.elem_faces();
  const IntTab& tab_face_voisins = domaine_VEF.face_voisins();
  const ArrOfInt& tab_est_face_bord = domaine_VEF.est_face_bord();

  int nb_dim = tab_transporte.nb_dim();
  int nb_comp = (nb_dim == 2 ? tab_transporte.dimension(1) : 1);

  DoubleTab tab_nu;
  int marq = phi_psi_diffuse(equation());
  const DoubleVect& porosite_elem = equation().milieu().porosite_elem();

  // either div(phi nu grad inco)
  // or div(nu grad phi inco)
  // depending on whether phi_psi or psi is diffused
  modif_par_porosite_si_flag(nu_, tab_nu, !marq, porosite_elem);
  DoubleVect tab_porosite_eventuelle(equation().milieu().porosite_face());
  if (!marq)
    tab_porosite_eventuelle = 1;

  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  int nb_bords = domaine_VEF.nb_front_Cl();
  int nb_faces_tot = domaine_VEF.nb_faces_tot();

  IntVect tab_face_associee(domaine_VEF.premiere_face_int());
  IntVect tab_fac2b_idx(domaine_VEF.nb_faces_tot());

  // Retrieve the indices of periodic boundary faces and the
  // associated faces into arrays beforehand for structured
  // access in the kernel afterwards
  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int num1 = le_bord.num_premiere_face();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());

          int nb_faces = le_bord.nb_faces();
          int num2b = num1 + nb_faces / 2;

          // only iterate over half the periodic faces
          // the result will be copied to the associated face at the end...
          ToDo_Kokkos("critical");
          for (int fac = num1; fac < num2b; fac++)
            {
              int fac_asso = la_cl_perio.face_associee(fac - num1) + num1;

              tab_face_associee(fac) = fac_asso;
              tab_fac2b_idx(fac) = fac;
              tab_fac2b_idx(fac+nb_faces/2) = fac;
            }
        }
    }

  CIntArrView est_face_bord = tab_est_face_bord.view_ro();
  CIntTabView face_voisins = tab_face_voisins.view_ro();
  CDoubleTabView face_normales = domaine_VEF.face_normales().view_ro();
  CDoubleArrView inverse_volumes = domaine_VEF.inverse_volumes().view_ro();
  CIntTabView elem_faces = tab_elem_faces.view_ro();

  CDoubleTabView nu = tab_nu.view_ro();
  CDoubleArrView porosite_eventuelle = tab_porosite_eventuelle.view_ro();

  CIntArrView fac2b_idx = tab_fac2b_idx.view_ro();
  CIntArrView face_associee = tab_face_associee.view_ro();

  Matrice_Morse_View matrice;
  matrice.set(tab_matrice);

  auto kern_elem_faces = KOKKOS_LAMBDA (const int fac)
  {
    int type_face = est_face_bord(fac);
    int fac2b = fac2b_idx(fac);

    // Periodic boundary faces
    if (type_face == 2 && fac == fac2b)
      {
        int fac_asso = face_associee(fac);

        int elem1 = face_voisins(fac, 0);
        int elem2 = face_voisins(fac, 1);

        for (int i = 0; i < nb_faces_elem; i++)
          {
            int j = elem_faces(elem1, i);
            if (j > fac)
              {
                for (int nc = 0; nc < nb_comp; nc++)
                  {
                    double val = viscA(fac, j, elem1, nu(elem1, nc), face_voisins, face_normales, inverse_volumes);

                    int n0 = fac * nb_comp + nc;
                    int j0 = j * nb_comp + nc;

                    matrice.atomic_add(n0, n0, + val * porosite_eventuelle(fac));
                    matrice.atomic_add(n0, j0, - val * porosite_eventuelle(j));
                    matrice.atomic_add(j0, n0, - val * porosite_eventuelle(fac));
                    matrice.atomic_add(j0, j0, + val * porosite_eventuelle(j));
                  }
              }
            if (elem2 != -1)
              {
                j = elem_faces(elem2, i);
                if (j > fac)
                  {
                    for (int nc = 0; nc < nb_comp; nc++)
                      {
                        double val = viscA(fac, j, elem2, nu(elem1, nc), face_voisins, face_normales, inverse_volumes);

                        int n0 = fac * nb_comp + nc;
                        int j0 = j * nb_comp + nc;
                        int n0perio = fac_asso * nb_comp + nc;

                        matrice.atomic_add(n0, n0, + val * porosite_eventuelle(fac));
                        matrice.atomic_add(n0, j0, - val * porosite_eventuelle(j));
                        matrice.atomic_add(j0, n0perio, - val * porosite_eventuelle(fac));
                        matrice.atomic_add(j0, j0, + val * porosite_eventuelle(j));
                      }
                  }
              }
          }
      }
    // Non-periodic boundary faces
    else if (type_face == 1)
      {
        int elem1 = face_voisins(fac, 0);

        for (int i = 0; i < nb_faces_elem; i++)
          {
            int j = elem_faces(elem1, i);
            if (j > fac)
              {
                for (int nc = 0; nc < nb_comp; nc++)
                  {
                    double val = viscA(fac, j, elem1, nu(elem1, nc), face_voisins, face_normales, inverse_volumes);

                    int n0 = fac * nb_comp + nc;
                    int j0 = j * nb_comp + nc;

                    matrice.atomic_add(n0, n0, + val * porosite_eventuelle(fac));
                    matrice.atomic_add(n0, j0, - val * porosite_eventuelle(j));
                    matrice.atomic_add(j0, n0, - val * porosite_eventuelle(fac));
                    matrice.atomic_add(j0, j0, + val * porosite_eventuelle(j));
                  }
              }
          }
      }
    // Faces internes
    else if (type_face == 0)
      {
        int elem1 = face_voisins(fac, 0);
        int elem2 = face_voisins(fac, 1);

        for (int i = 0; i < nb_faces_elem; i++)
          {
            int j = elem_faces(elem1, i);
            if (j > fac)
              {
                for (int nc = 0; nc < nb_comp; nc++)
                  {
                    double val = viscA(fac, j, elem1, nu(elem1, nc), face_voisins, face_normales, inverse_volumes);

                    int n0 = fac * nb_comp + nc;
                    int j0 = j * nb_comp + nc;

                    matrice.atomic_add(n0, n0, + val * porosite_eventuelle(fac));
                    matrice.atomic_add(n0, j0, - val * porosite_eventuelle(j));
                    matrice.atomic_add(j0, n0, - val * porosite_eventuelle(fac));
                    matrice.atomic_add(j0, j0, + val * porosite_eventuelle(j));
                  }
              }

            if (elem2 != -1)
              {
                j = elem_faces(elem2, i);
                if (j > fac)
                  {
                    for (int nc = 0; nc < nb_comp; nc++)
                      {
                        double val = viscA(fac, j, elem2, nu(elem2, nc), face_voisins, face_normales, inverse_volumes);
                        int n0 = fac * nb_comp + nc;
                        int j0 = j * nb_comp + nc;

                        matrice.atomic_add(n0, n0, + val * porosite_eventuelle(fac));
                        matrice.atomic_add(n0, j0, - val * porosite_eventuelle(j));
                        matrice.atomic_add(j0, n0, - val * porosite_eventuelle(fac));
                        matrice.atomic_add(j0, j0, + val * porosite_eventuelle(j));
                      }
                  }
              }
          }
      }
  };

  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), nb_faces_tot, kern_elem_faces);
  end_gpu_timer(__KERNEL_NAME__);

  for (int n_bord = 0; n_bord < nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());

      if (sub_type(Echange_externe_impose, la_cl.valeur()))
        {
          const Echange_externe_impose& la_cl_paroi = ref_cast(Echange_externe_impose, la_cl.valeur());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          ToDo_Kokkos("critical");
          for (int face = ndeb; face < nfin; face++)
            for (int nc = 0; nc < nb_comp; nc++)
              {
                const int i = face * nb_comp + nc;
                tab_matrice(i, i) += la_cl_paroi.h_imp(face - ndeb, nc) * domaine_VEF.surface(face);
              }
        }
      if (sub_type(Echange_externe_radiatif, la_cl.valeur()))
        {
          throw;
        }
    }

  modifier_matrice_pour_periodique_apres_contribuer(tab_matrice, equation());
}

void Op_Diff_VEF_Face::contribue_au_second_membre(DoubleTab& resu ) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  int nb_comp = 1;
  int nb_dim = resu.nb_dim();

  int nb_bords=domaine_VEF.nb_front_Cl();

  if(nb_dim==2)
    nb_comp=resu.dimension(1);

  // Partie imposee :

  if (nb_dim == 1)
    {
      for (int n_bord=0; n_bord<nb_bords; n_bord++)
        {
          const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

          if (sub_type(Neumann_paroi,la_cl.valeur()))
            {
              const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              ToDo_Kokkos("critical");
              for (int face=ndeb; face<nfin; face++)
                resu[face] += la_cl_paroi.flux_impose(face-ndeb)*domaine_VEF.surface(face);
            }
          else if (sub_type(Echange_externe_impose,la_cl.valeur()))
            {
              Cerr << "Non code pour Echange_externe_impose" << finl;
              assert(0);
            }
          else if (sub_type(Echange_externe_radiatif,la_cl.valeur()))
            {
              Cerr << "Non code pour Echange_externe_radiatif" << finl;
              assert(0);
            }

        }
    }
  else
    {
      for (int n_bord=0; n_bord<nb_bords; n_bord++)
        {
          const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

          if (sub_type(Neumann_paroi,la_cl.valeur()))
            {
              const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int ndeb = le_bord.num_premiere_face();
              int nfin = ndeb + le_bord.nb_faces();
              ToDo_Kokkos("critical");
              for (int face=ndeb; face<nfin; face++)
                for (int comp=0; comp<nb_comp; comp++)
                  resu(face,comp) += la_cl_paroi.flux_impose(face-ndeb,comp)*domaine_VEF.surface(face);
            }
        }
    }
}

void Op_Diff_VEF_Face::verifier() const
{
  static int testee=0;
  if(testee)
    return;
  testee=1;
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  //  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  //  const Conds_lim& les_cl = domaine_Cl_VEF.les_conditions_limites();
  const DoubleVect& volumes_entrelaces = domaine_VEF.volumes_entrelaces();

  const DoubleTab& xv=domaine_VEF.xv();
  DoubleTab vit(equation().inconnue().valeurs());
  DoubleTab resu(vit);
  int i, comp;
  if(dimension==2)
    {
      const int nbf = vit.dimension(0);
      Cerr << " Verification of delta(x,0) " << finl;
      for(i=0; i<nbf; i++)
        {
          vit(i,0)=xv(i,0);
          vit(i,1)=0;
        }
      calculer(vit, resu);
      for(i=0; i<nbf; i++)
        for(comp=0; comp<dimension; comp++)
          resu(i,comp)/=(volumes_entrelaces(i));
      for(i=0; i<nbf; i++)
        {
          if(std::fabs(resu(i,0))>1.e-10)
            {
              Cerr << " delta(x,0) ("<<i<<") = "
                   << resu(i,0);
              Cerr << finl;
            }
        }
      Cerr << " Verification of delta(y(1-y),0) " << finl;
      for(i=0; i<nbf; i++)
        {
          vit(i,0)=xv(i,1)*(1-xv(i,1));
          vit(i,1)=0;
        }
      calculer(vit, resu);
      for(i=0; i<nbf; i++)
        for(comp=0; comp<dimension; comp++)
          resu(i,comp)/=(volumes_entrelaces(i));
      for(i=0; i<nbf; i++)
        {
          if(std::fabs(2-resu(i,0))>1.e-10)
            {
              Cerr << " delta(y(1-y),0) ("<<i<<") = "
                   << resu(i,0);
              Cerr << finl;
            }
        }
    }
}



