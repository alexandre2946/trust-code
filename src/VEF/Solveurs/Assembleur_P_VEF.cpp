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

#include <Assembleur_P_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Domaine_VEF.h>
#include <Periodique.h>
#include <Neumann_sortie_libre.h>
#include <Matrice_Bloc.h>
#include <Milieu_base.h>
#include <Robin_VEF.h>


Implemente_instanciable(Assembleur_P_VEF,"Assembleur_P_VEF",Assembleur_base);

Sortie& Assembleur_P_VEF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom() ;
}

Entree& Assembleur_P_VEF::readOn(Entree& s )
{
  return Assembleur_base::readOn(s);
}

int Assembleur_P_VEF::assembler(Matrice& la_matrice)
{
  // If rho is constant, solve with pressure P*=P/rho
  const DoubleVect& volumes_entrelaces_ref=le_dom_VEF->volumes_entrelaces();
  DoubleVect tab_volumes_entrelaces(volumes_entrelaces_ref);
  const DoubleVect& tab_volumes_entrelaces_cl=le_dom_Cl_VEF->volumes_entrelaces_Cl();
  int size=tab_volumes_entrelaces_cl.size();
  {
    CDoubleArrView volumes_entrelaces_cl = static_cast<const ArrOfDouble&>(tab_volumes_entrelaces_cl).view_ro();
    DoubleArrView volumes_entrelaces = static_cast<ArrOfDouble&>(tab_volumes_entrelaces).view_rw();
    Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, size), KOKKOS_LAMBDA(const int f)
    {
      if (volumes_entrelaces_cl(f)!=0)
        volumes_entrelaces(f)=volumes_entrelaces_cl(f);
    });
    end_gpu_timer(__KERNEL_NAME__);
  }

  tab_volumes_entrelaces.echange_espace_virtuel();
  // Assemble the matrix
  return assembler_mat(la_matrice,tab_volumes_entrelaces,1,1);
}



void calculer_inv_volume_special(DoubleTab& tab_inv_volumes_entrelaces, const Domaine_Cl_VEF& domaine_Cl_VEF,const DoubleTab& tab_volumes_entrelaces)
{
  tab_inv_volumes_entrelaces=tab_volumes_entrelaces;
  int taille=tab_volumes_entrelaces.dimension_tot(0);
  {
    CDoubleTabView volumes_entrelaces = tab_volumes_entrelaces.view_ro();
    DoubleTabView inv_volumes_entrelaces = tab_inv_volumes_entrelaces.view_rw();
    const int dimension = Objet_U::dimension;
    Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, taille), KOKKOS_LAMBDA(const int i)
    {
      for (int comp=0; comp<dimension; comp++)
        inv_volumes_entrelaces(i,comp)=1./volumes_entrelaces(i,comp);
    });
    end_gpu_timer(__KERNEL_NAME__);
  }
}
void Assembleur_P_VEF::calculer_inv_volume(DoubleTab& inv_volumes_entrelaces, const Domaine_Cl_VEF& domaine_Cl_VEF,const DoubleVect& volumes_entrelaces)
{
  // the inverse of the volume is now a DoubleTab
  // this is to make Piso work
  const DoubleTab* doubleT = dynamic_cast<const DoubleTab*>(&volumes_entrelaces);
  if (doubleT)
    {
      calculer_inv_volume_special(inv_volumes_entrelaces, domaine_Cl_VEF,*doubleT);
      return;
    }
  int taille=volumes_entrelaces.size_totale();
  inv_volumes_entrelaces.resize(taille,Objet_U::dimension);
  if (0)
    {
      // this way of computing the interlaced volume may be good in
      // the future provided the porosity issue is addressed,
      // and ESPECIALLY the SIMPLER issue where we do not want to go through this
      DoubleTab tmp;
      tmp=(inv_volumes_entrelaces);
      tmp=1;
      domaine_Cl_VEF.equation().solv_masse().appliquer(tmp);
      int sz=inv_volumes_entrelaces.size_totale();
      for (int i=0; i<sz; i++)
        {
          inv_volumes_entrelaces(i)=tmp(i);
          if (!est_egal(inv_volumes_entrelaces(i),1./volumes_entrelaces(i)))
            Cerr<<i<<" "<<inv_volumes_entrelaces(i)-1./volumes_entrelaces(i)<<" "<<inv_volumes_entrelaces(i)<<finl;;
        }
      Process::exit();
    }
  else
    {
      CDoubleArrView porosite_face = static_cast<const ArrOfDouble&>(equation().milieu().porosite_face()).view_ro();
      CDoubleArrView volumes_entrelaces_v = static_cast<const ArrOfDouble&>(volumes_entrelaces).view_ro();
      DoubleTabView inv_volumes_entrelaces_v = inv_volumes_entrelaces.view_rw();
      const int dim = Objet_U::dimension;
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, taille), KOKKOS_LAMBDA(const int i)
      {
        for (int comp=0; comp<dim; comp++)
          inv_volumes_entrelaces_v(i,comp)=1./volumes_entrelaces_v(i)*porosite_face(i);
      });
      end_gpu_timer(__KERNEL_NAME__);
    }
}


int Assembleur_P_VEF::assembler_mat(Matrice& la_matrice, const DoubleVect& volumes_entrelaces, int incr_pression, int resoudre_en_u)
{
  // Set flags of Assembleur_base
  set_resoudre_increment_pression(incr_pression);
  set_resoudre_en_u(resoudre_en_u);
  const Domaine_Cl_VEF& domaine_Cl_VEF = le_dom_Cl_VEF.valeur();
  DoubleTab inverse_quantitee_entrelacee;
  calculer_inv_volume(inverse_quantitee_entrelacee, domaine_Cl_VEF, volumes_entrelaces);
  remplir(la_matrice, inverse_quantitee_entrelacee);
  modifier_matrice(la_matrice);
  return 1;
}

int Assembleur_P_VEF::remplir(Matrice& la_matrice, const DoubleTab& inverse_quantitee_entrelacee)
{
  has_P_ref=0;
  // Pressure matrix: sparse matrix of size nb_poly x nb_poly
  // This function stores the matrix in a Morse matrix structure
  // well suited for sparse matrices.
  // First, compute the sizes of arrays tab1 and tab2
  // (coeff_ has the same size as tab2)
  //   For each polyhedron, associate:
  //   - a list of integers neighbors[i] = {j>i such that Mij is nonzero}
  //   - a list of reals values[i] = {Mij for j in Neighbors[i]}
  //   - a real diagonal term
  // Temporary implementation:
  // Assemble a pressure matrix for a hydraulic equation
  // Inject boundary conditions into this matrix
  // This can be done because the pressure matrix is a priori not
  // shared between several equations on the same domain.

  const Domaine_VEF& le_dom = le_dom_VEF.valeur();
  const Domaine_Cl_VEF& le_dom_cl = le_dom_Cl_VEF.valeur();
  les_coeff_pression.resize(le_dom_cl.nb_faces_Cl());
  int n1 = le_dom.domaine().nb_elem_tot();
  int n2 = le_dom.domaine().nb_elem();


  // Add porosities.

  la_matrice.typer("Matrice_Bloc"); // Actually Matrice_Bloc_Sym?
  Matrice_Bloc& matrice=ref_cast(Matrice_Bloc, la_matrice.valeur());
  matrice.dimensionner(2,2);
  matrice.get_bloc(0,0).typer("Matrice_Morse_Sym");
  matrice.get_bloc(0,1).typer("Matrice_Morse");
  matrice.get_bloc(1,0).typer("Matrice_Morse");
  matrice.get_bloc(1,1).typer("Matrice_Morse"); // Actually Matrice_Morse_Sym ?

  Matrice_Morse_Sym& MBrr = ref_cast(Matrice_Morse_Sym,matrice.get_bloc(0,0).valeur());
  Matrice_Morse& MBrv = ref_cast (Matrice_Morse,matrice.get_bloc(0,1).valeur());
  Matrice_Morse& MBvr = ref_cast (Matrice_Morse, matrice.get_bloc(1,0).valeur());
  Matrice_Morse& MBvv = ref_cast (Matrice_Morse, matrice.get_bloc(1,1).valeur());

  MBrr.dimensionner(n2,0);
  MBrv.dimensionner(n2,0);
  MBvv.dimensionner(n1-n2,0);
  // The sub-block vr is sized and set to zero
  MBvr.dimensionner(n1-n2,n2,0);
  MBvr.get_set_tab1() = 1;

  // Process internal faces:
  int ndeb = le_dom_VEF->premiere_face_int();
  int nfin = le_dom_VEF->nb_faces_tot();
  int nb_faces = le_dom_VEF->nb_faces();

#ifdef TRUST_USE_GPU
  ArrOfTID rang_voisinRR(n2);
  ArrOfTID rang_voisinRV(n2);
  ArrOfTID rang_voisinVV(n1-n2);
#else
  ArrOfInt rang_voisinRR(n2);
  ArrOfInt rang_voisinRV(n2);
  ArrOfInt rang_voisinVV(n1-n2);
#endif
  rang_voisinRR=1; // Diagonale
  rang_voisinRV=0; // No diagonal
  rang_voisinVV=1; // Diagonale

  CIntTabView face_voisins = le_dom.face_voisins().view_ro();
  CIntArrView ind_faces_virt_bord = le_dom_VEF->ind_faces_virt_bord().view_ro();
  auto rang_voisinRR_v = rang_voisinRR.view_rw();
  auto rang_voisinRV_v = rang_voisinRV.view_rw();
  auto rang_voisinVV_v = rang_voisinVV.view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(ndeb, nfin), KOKKOS_LAMBDA(const int num_face)
  {
    int elem1 = face_voisins(num_face, 0);
    int elem2 = face_voisins(num_face, 1);
    const bool is_face_virt_bord = (num_face >= nb_faces) && (ind_faces_virt_bord(num_face - nb_faces) != -1);
    if (!is_face_virt_bord && elem1 != -1 && elem2 != -1)
      {
        if (elem1 > elem2)
          {
            if(elem1 < n2)
              Kokkos::atomic_add(&rang_voisinRR_v(elem2), 1);
            else
              {
                if(elem2 < n2)
                  Kokkos::atomic_add(&rang_voisinRV_v(elem2), 1);
                else
                  Kokkos::atomic_add(&rang_voisinVV_v(elem2 - n2), 1);
              }
          }
        else // elem2 >= elem1
          {
            if(elem2 < n2)
              Kokkos::atomic_add(&rang_voisinRR_v(elem1), 1);
            else
              {
                if(elem1 < n2)
                  Kokkos::atomic_add(&rang_voisinRV_v(elem1), 1);
                else
                  Kokkos::atomic_add(&rang_voisinVV_v(elem1 - n2), 1);
              }
          }
      }
  });
  end_gpu_timer(__KERNEL_NAME__);

  // Account for periodic boundary conditions
  const Conds_lim& les_cl = le_dom_cl.les_conditions_limites();
  for (int i=0; i<les_cl.size(); i++)
    {
      const Cond_lim& la_cl = les_cl[i];

      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int nb_faces_bord_tot = le_bord.nb_faces_tot();
          CIntArrView front_num_face = le_bord.num_face().view_ro();
          CIntArrView face_associee = la_cl_perio.face_associee().view_ro();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, nb_faces_bord_tot), KOKKOS_LAMBDA(const int ind_face)
          {
            if (ind_face < face_associee(ind_face)) // Process each periodic pair once
              {
                int num_face = front_num_face(ind_face);
                int elem1 = face_voisins(num_face, 0);
                int elem2 = face_voisins(num_face, 1);
                if (elem1 != -1 && elem2 != -1)
                  {
                    if (elem1 > elem2)
                      {
                        if(elem1 < n2)
                          Kokkos::atomic_add(&rang_voisinRR_v(elem2), 1);
                        else
                          {
                            if(elem2 < n2)
                              Kokkos::atomic_add(&rang_voisinRV_v(elem2), 1);
                            else
                              Kokkos::atomic_add(&rang_voisinVV_v(elem2 - n2), 1);
                          }
                      }
                    else // elem2 >= elem1
                      {
                        if(elem2 < n2)
                          Kokkos::atomic_add(&rang_voisinRR_v(elem1), 1);
                        else
                          {
                            if(elem1 < n2)
                              Kokkos::atomic_add(&rang_voisinRV_v(elem1), 1);
                            else
                              Kokkos::atomic_add(&rang_voisinVV_v(elem1 - n2), 1);
                          }
                      }
                  }
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }

  auto& tab1RR = MBrr.get_set_tab1();
  auto& tab2RR = MBrr.get_set_tab2();
  auto& tab1RV = MBrv.get_set_tab1();
  auto& tab2RV = MBrv.get_set_tab2();
  auto& tab1VV = MBvv.get_set_tab1();
  auto& tab2VV = MBvv.get_set_tab2();

  tab1RR(0)=1;
  tab1RV(0)=1;
  tab1VV(0)=1;
  auto tab1RR_v = tab1RR.view_rw();
  auto tab1RV_v = tab1RV.view_rw();
  auto tab1VV_v = tab1VV.view_rw();
  using tab1_value_t = typename decltype(tab1RR_v)::value_type;
  Kokkos::parallel_scan(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n2), KOKKOS_LAMBDA(const int i, tab1_value_t& update, const bool final)
  {
    update += rang_voisinRR_v(i);
    if (final) tab1RR_v(i+1) = update + 1;
  });
  end_gpu_timer(__KERNEL_NAME__);
  Kokkos::parallel_scan(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n2), KOKKOS_LAMBDA(const int i, tab1_value_t& update, const bool final)
  {
    update += rang_voisinRV_v(i);
    if (final) tab1RV_v(i+1) = update + 1;
  });
  end_gpu_timer(__KERNEL_NAME__);
  Kokkos::parallel_scan(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n1-n2), KOKKOS_LAMBDA(const int i, tab1_value_t& update, const bool final)
  {
    update += rang_voisinVV_v(i);
    if (final) tab1VV_v(i+1) = update + 1;
  });
  end_gpu_timer(__KERNEL_NAME__);
  MBrr.dimensionner(n2,tab1RR(n2)-1);
  MBrv.dimensionner(n2,n1-n2,tab1RV(n2)-1);
  MBvv.dimensionner(n1-n2,n1-n2,tab1VV(n1-n2)-1);

  auto tab2RR_v = tab2RR.view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n2), KOKKOS_LAMBDA(const int i)
  {
    tab2RR_v(tab1RR_v(i) - 1) = i+1; // Diagonale
    rang_voisinRR_v(i) = tab1RR_v(i);
    rang_voisinRV_v(i) = tab1RV_v(i) - 1;
  });
  end_gpu_timer(__KERNEL_NAME__);
  auto tab2VV_v = tab2VV.view_rw();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, n1-n2), KOKKOS_LAMBDA(const int i)
  {
    tab2VV_v(tab1VV_v(i) - 1) = i+1; // Diagonale
    rang_voisinVV_v(i) = tab1VV_v(i);
  });
  end_gpu_timer(__KERNEL_NAME__);

  int dim = Objet_U::dimension;

  MBrr.dimensionner(n2,tab1RR(n2)-1);
  MBrv.dimensionner(n2,n1-n2,tab1RV(n2)-1);
  MBvv.dimensionner(n1-n2,n1-n2,tab1VV(n1-n2)-1);

  auto tab2RV_v = tab2RV.view_rw();
  DoubleArrView coeffRR = MBrr.get_set_coeff().view_rw();
  DoubleArrView coeffRV = MBrv.get_set_coeff().view_rw();
  DoubleArrView coeffVV = MBvv.get_set_coeff().view_rw();
  CDoubleTabView face_normales = le_dom.face_normales().view_ro();
  CDoubleTabView inverse_quantitee_entrelacee_v = inverse_quantitee_entrelacee.view_ro();
  Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(ndeb, nfin), KOKKOS_LAMBDA(const int num_face)
  {
    int elem1 = face_voisins(num_face, 0);
    int elem2 = face_voisins(num_face, 1);
    const bool is_face_virt_bord = (num_face >= nb_faces) && (ind_faces_virt_bord(num_face - nb_faces) != -1);
    if (!is_face_virt_bord && elem1 != -1 && elem2 != -1)
      {
        double val = 0.;
        for (int d = 0; d < dim; d++)
          val += face_normales(num_face, d) * face_normales(num_face, d) * inverse_quantitee_entrelacee_v(num_face, d);

        // diagonale :
        if (elem1 < n2) Kokkos::atomic_add(&coeffRR(tab1RR_v(elem1) - 1), val);
        else            Kokkos::atomic_add(&coeffVV(tab1VV_v(elem1 - n2) - 1), val);
        if (elem2 < n2) Kokkos::atomic_add(&coeffRR(tab1RR_v(elem2) - 1), val);
        else            Kokkos::atomic_add(&coeffVV(tab1VV_v(elem2 - n2) - 1), val);

        if (elem1 > elem2)
          {
            if (elem1 < n2)
              {
                auto slot = Kokkos::atomic_fetch_add(&rang_voisinRR_v(elem2), 1);
                tab2RR_v(slot) = elem1 + 1;
                coeffRR(slot) -= val;
              }
            else
              {
                if (elem2 < n2)
                  {
                    auto slot = Kokkos::atomic_fetch_add(&rang_voisinRV_v(elem2), 1);
                    tab2RV_v(slot) = (elem1 - n2) + 1;
                    coeffRV(slot) -= val;
                  }
                else
                  {
                    auto slot = Kokkos::atomic_fetch_add(&rang_voisinVV_v(elem2 - n2), 1);
                    tab2VV_v(slot) = (elem1 - n2) + 1;
                    coeffVV(slot) -= val;
                  }
              }
          }
        else
          {
            if (elem2 < n2)
              {
                auto slot = Kokkos::atomic_fetch_add(&rang_voisinRR_v(elem1), 1);
                tab2RR_v(slot) = elem2 + 1;
                coeffRR(slot) -= val;
              }
            else
              {
                if (elem1 < n2)
                  {
                    auto slot = Kokkos::atomic_fetch_add(&rang_voisinRV_v(elem1), 1);
                    tab2RV_v(slot) = (elem2 - n2) + 1;
                    coeffRV(slot) -= val;
                  }
                else
                  {
                    auto slot = Kokkos::atomic_fetch_add(&rang_voisinVV_v(elem1 - n2), 1);
                    tab2VV_v(slot) = (elem2 - n2) + 1;
                    coeffVV(slot) -= val;
                  }
              }
          }
      }
  });
  end_gpu_timer(__KERNEL_NAME__);
  // Process the boundary conditions
  for (int i=0; i<les_cl.size(); i++)
    {

      // Treatment depends on the type of boundary condition:
      //  - For Neumann_sortie_libre or Robin_VEF conditions:
      //    compute the coefficient on the face and account for it
      //    in the diagonal term of the neighboring element.
      // - For a BC face with any other condition: no
      //   contribution to the pressure matrix.

      const Cond_lim& la_cl = les_cl[i];
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      int nb_faces_bord_tot = le_bord.nb_faces_tot();
      if (sub_type(Neumann_sortie_libre,la_cl.valeur()) || sub_type(Robin_VEF,la_cl.valeur())  )
        {
          has_P_ref=1;
          MBrr.set_est_definie(1);
          CIntArrView front_num_face = le_bord.num_face().view_ro();
          DoubleArrView coeff_pression = static_cast<ArrOfDouble&>(les_coeff_pression).view_rw();
          const int coeff_pression_size = les_coeff_pression.size_array();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, nb_faces_bord_tot), KOKKOS_LAMBDA(const int ind_face)
          {
            int num_face = front_num_face(ind_face);
            double val = 0.;
            for (int d = 0; d < dim; d++)
              val += face_normales(num_face, d) * face_normales(num_face, d) * inverse_quantitee_entrelacee_v(num_face, d);

            int elem = face_voisins(num_face, 0);
            if (elem < n2) Kokkos::atomic_add(&coeffRR(tab1RR_v(elem) - 1), val);
            else           Kokkos::atomic_add(&coeffVV(tab1VV_v(elem - n2) - 1), val);
            // Store the pressure coefficients on the real faces
            if (num_face < coeff_pression_size)
              coeff_pression(num_face) = val;
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
      else if (sub_type(Periodique,la_cl.valeur()) )
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          CIntArrView front_num_face_coeff = le_bord.num_face().view_ro();
          CIntArrView face_associee_coeff = la_cl_perio.face_associee().view_ro();
          Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, nb_faces_bord_tot), KOKKOS_LAMBDA(const int ind_face)
          {
            if (ind_face < face_associee_coeff(ind_face)) // Process each periodic pair once
              {
                int num_face = front_num_face_coeff(ind_face);
                int elem1 = face_voisins(num_face, 0);
                int elem2 = face_voisins(num_face, 1);
                double val = 0.;
                for (int d = 0; d < dim; d++)
                  val += face_normales(num_face, d) * face_normales(num_face, d) * inverse_quantitee_entrelacee_v(num_face, d);

                // diagonale :
                if (elem1 < n2) Kokkos::atomic_add(&coeffRR(tab1RR_v(elem1) - 1), val);
                else            Kokkos::atomic_add(&coeffVV(tab1VV_v(elem1 - n2) - 1), val);
                if (elem2 < n2) Kokkos::atomic_add(&coeffRR(tab1RR_v(elem2) - 1), val);
                else            Kokkos::atomic_add(&coeffVV(tab1VV_v(elem2 - n2) - 1), val);

                if (elem1 > elem2)
                  {
                    if (elem1 < n2)
                      {
                        auto slot = Kokkos::atomic_fetch_add(&rang_voisinRR_v(elem2), 1);
                        tab2RR_v(slot) = elem1 + 1;
                        coeffRR(slot) -= val;
                      }
                    else
                      {
                        if (elem2 < n2)
                          {
                            auto slot = Kokkos::atomic_fetch_add(&rang_voisinRV_v(elem2), 1);
                            tab2RV_v(slot) = elem1 - n2 + 1;
                            coeffRV(slot) -= val;
                          }
                        else
                          {
                            auto slot = Kokkos::atomic_fetch_add(&rang_voisinVV_v(elem2 - n2), 1);
                            tab2VV_v(slot) = elem1 - n2 + 1;
                            coeffVV(slot) -= val;
                          }
                      }
                  }
                else
                  {
                    if (elem2 < n2)
                      {
                        auto slot = Kokkos::atomic_fetch_add(&rang_voisinRR_v(elem1), 1);
                        tab2RR_v(slot) = elem2 + 1;
                        coeffRR(slot) -= val;
                      }
                    else
                      {
                        if (elem1 < n2)
                          {
                            auto slot = Kokkos::atomic_fetch_add(&rang_voisinRV_v(elem1), 1);
                            tab2RV_v(slot) = elem2 - n2 + 1;
                            coeffRV(slot) -= val;
                          }
                        else
                          {
                            auto slot = Kokkos::atomic_fetch_add(&rang_voisinVV_v(elem1 - n2), 1);
                            tab2VV_v(slot) = elem2 - n2 + 1;
                            coeffVV(slot) -= val;
                          }
                      }
                  }
              }
          });
          end_gpu_timer(__KERNEL_NAME__);
        }
    }
  has_P_ref = (int)mp_max(has_P_ref);
  return 1;
}

/*! @brief Assembles the pressure matrix for a quasi-compressible fluid: laplacian(P) is replaced by div(grad(P)/rho).
 *
 * @param tab_rho Mass density array.
 * @return Always returns 1.
 */
int Assembleur_P_VEF::assembler_QC(const DoubleTab& tab_rho, Matrice& matrice)
{
  Cerr << "Assembling pressure matrix for Quasi-Compressible in progress..." << finl;
  int stat=assembler(matrice);
  set_resoudre_en_u(0);
  return stat;
}

int Assembleur_P_VEF::modifier_secmem(DoubleTab& secmem)
{
  const Domaine_VEF& le_dom = le_dom_VEF.valeur();
  const Domaine_Cl_VEF& le_dom_cl = le_dom_Cl_VEF.valeur();
  int nb_cond_lim = le_dom_cl.nb_cond_lim();
  const IntTab& face_voisins = le_dom.face_voisins();

  // Modification of the right-hand side:
  for (int i=0; i<nb_cond_lim; i++)
    {
      const Cond_lim_base& la_cl_base = le_dom_cl.les_conditions_limites(i).valeur();
      const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl_base.frontiere_dis());
      const Champ_front_base& champ_front = la_cl_base.champ_front();
      int ndeb = la_front_dis.num_premiere_face();
      int nfin = ndeb + la_front_dis.nb_faces();

      // GF we switched to pressure increment
      if ((sub_type(Neumann_sortie_libre,la_cl_base)) && (!get_resoudre_increment_pression()))
        {
          const Neumann_sortie_libre& la_cl_Neumann = ref_cast(Neumann_sortie_libre, la_cl_base);
          ToDo_Kokkos("critical");
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              double Pimp = la_cl_Neumann.flux_impose(num_face-ndeb);
              double coef = les_coeff_pression[num_face]*Pimp;
              secmem[face_voisins(num_face,0)] += coef;
            }
        }

      /*if ((sub_type(Robin_VEF,la_cl_base)) )
        {
          // [vkr/oswr] we should modify the value of increment_pression_bord in order to respect the oswr algorithm
          double Pstar_OSWR, coef;
          const Robin_VEF& la_cl_robin = ref_cast(Robin_VEF, la_cl_base);
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              Pstar_OSWR = la_cl_robin.increment_pression_bord(num_face);
              coef = les_coeff_pression[num_face]*Pstar_OSWR+100;
              Cerr << "PASSING THROUGH GET RESOUDRE INCR PRESSION" <<finl;
              secmem[face_voisins(num_face,0)] += coef;
            }
        }*/

      else if ( champ_front.instationnaire() && get_resoudre_en_u() )
        {
          const DoubleTab& Gpt = champ_front.derivee_en_temps();
          bool ch_unif = (Gpt.nb_dim()==1);
          ToDo_Kokkos("critical");
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              double Stt = 0.;
              for (int k=0; k<dimension; k++)
                {
                  double Gpoint = ch_unif ? Gpt(k) : Gpt(num_face - ndeb, k);
                  Stt -= Gpoint * le_dom.face_normales(num_face, k);
                }
              secmem(face_voisins(num_face,0)) += Stt;
            }
        }
    }
  secmem.echange_espace_virtuel();
  return 1;
}

int Assembleur_P_VEF::modifier_solution(DoubleTab& pression)
{
  // Projection :
  double press_0;
  if(!has_P_ref)
    {
      // Take the minimum pressure as the reference pressure
      // in order to have the same reference pressure in sequential and parallel
      press_0=DMAXFLOAT;
      int nb_elem=le_dom_VEF->domaine().nb_elem();
      ToDo_Kokkos("critical");
      for(int n=0; n<nb_elem; n++)
        if (pression[n] < press_0)
          press_0 = pression[n];
      press_0 = Process::mp_min(press_0);
      ToDo_Kokkos("critical");
      for(int n=0; n<nb_elem; n++)
        pression[n] -=press_0;

      pression.echange_espace_virtuel();
    }
  return 1;
}

/*! @brief Optionally modifies the matrix to make it definite if it is not.
 *
 * @return 1 if the matrix is modified, 0 otherwise.
 */
int Assembleur_P_VEF::modifier_matrice(Matrice& matrice)
{
  int matrice_modifiee=0;
  Matrice_Bloc& mat_bloc = ref_cast(Matrice_Bloc, matrice.valeur());
  Matrice_Morse_Sym& A00RR = ref_cast(Matrice_Morse_Sym,mat_bloc.get_bloc(0,0).valeur());
  // Find the element on which the reference pressure is imposed
  const bool is_first_proc_with_real_elems = Process::me() == Process::mp_min(le_dom_VEF->nb_elem() ? Process::me() : 1e8);
  if (is_first_proc_with_real_elems && !A00RR.get_est_definie())
    {
      int element_referent=0;
      double distance=DMAXFLOAT;
      const DoubleTab& coord=le_dom_VEF->xp();
      int n = le_dom_VEF->nb_elem();
      ToDo_Kokkos("critical");
      for(int i=0; i<n; i++)
        {
          double tmp=0;
          for (int j=0; j<dimension; j++)
            tmp+=coord(i,j)*coord(i,j);
          if (inf_strict(tmp,distance) && !est_egal(A00RR(i,i),0.))
            {
              distance=tmp;
              element_referent=i;
            }
        }
      Cerr << "Modifying the row (element) " << element_referent << finl;
      A00RR(element_referent,element_referent)*=2;
      matrice_modifiee=1;
    }
  // has_P_ref=1;
  A00RR.set_est_definie(1);
  return matrice_modifiee;
}

const Domaine_dis_base& Assembleur_P_VEF::domaine_dis_base() const
{
  return le_dom_VEF.valeur();
}

const Domaine_Cl_dis_base& Assembleur_P_VEF::domaine_Cl_dis_base() const
{
  return le_dom_Cl_VEF.valeur();
}

void Assembleur_P_VEF::associer_domaine_dis_base(const Domaine_dis_base& le_dom_dis)
{
  le_dom_VEF = ref_cast(Domaine_VEF,le_dom_dis);
}

void Assembleur_P_VEF::associer_domaine_cl_dis_base(const Domaine_Cl_dis_base& le_dom_Cl_dis)
{
  le_dom_Cl_VEF = ref_cast(Domaine_Cl_VEF, le_dom_Cl_dis);
}

void Assembleur_P_VEF::completer(const Equation_base& Eqn)
{
  mon_equation=Eqn;
}
