/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Champ_val_tot_sur_vol_VEF.h>
#include <Equation_base.h>
#include <Domaine_Cl_VEF.h>
#include <Dirichlet.h>
#include <Sous_Domaine.h>
#include <Domaine_VEF.h>

Implemente_instanciable(Champ_val_tot_sur_vol_VEF,"Valeur_totale_sur_volume_VEF",Champ_val_tot_sur_vol_base);


Sortie& Champ_val_tot_sur_vol_VEF::printOn(Sortie& os) const
{
  Champ_val_tot_sur_vol_base::printOn(os);
  return os;
}

Entree& Champ_val_tot_sur_vol_VEF::readOn(Entree& is)
{
  Champ_val_tot_sur_vol_base::readOn(is);
  return is;
}

DoubleVect& Champ_val_tot_sur_vol_VEF::eval_contrib_loc(const Domaine_dis_base& zdis,const Domaine_Cl_dis_base& zcldis,DoubleVect& vol_glob_pond)
{
  const Domaine_VEF& zvef = ref_cast(Domaine_VEF,zdis);
  const Domaine_Cl_VEF& zclvef = ref_cast(Domaine_Cl_VEF,zcldis);
  const int nb_elem = zvef.nb_elem();
  int size_vol = les_sous_domaines.size()+1;
  vol_glob_pond.resize(size_vol);

  const int nb_faces = zvef.nb_faces();
  const int nb_fac_el = zvef.domaine().nb_faces_elem();
  const IntTab& elem_faces = zvef.elem_faces();
  const DoubleVect& vol_entrelaces = zvef.volumes_entrelaces();
  const DoubleVect& vol_entrelaces_Cl =  zclvef.volumes_entrelaces_Cl();
  const ArrOfInt& faces_doubles = zvef.faces_doubles();
  const DoubleVect& por_face = zclvef.equation().milieu().porosite_face();
  int prem_face_std = zvef.premiere_face_std();
  int face_g,face_marq;
  double fac_pond,vol_entrelace;

  int cpt=1;
  IntVect face_fait(nb_faces);

  for (auto& itr : les_sous_domaines)
    {
      const Sous_Domaine& sz = itr.valeur();
      int size_sz = sz.nb_elem_tot();
      face_fait = 0;
      int el,elem0,elem1,elem_test;

      for (int elem=0; elem<size_sz; elem++)
        {
          el = sz(elem);
          // Keep only the real elements
          if (el<nb_elem)
            {
              for (int fac=0; fac<nb_fac_el; fac++)
                {
                  face_g = elem_faces(el,fac);

                  if (!face_fait(face_g))
                    {

                      // Among the two neighboring elements of the processed face, keep
                      // the one that is not the current element in the sub-domain
                      elem0 = zvef.face_voisins(face_g,0);
                      elem1 = zvef.face_voisins(face_g,1);
                      if (elem0==el)
                        elem_test = elem1;
                      else
                        elem_test = elem0;

                      // Four possible situations:
                      // - elem_test=-1: boundary condition
                      // - elem_test is also in the current sub-domain
                      // - elem_test is in another sub-domain
                      // - elem_test is in the default part of the domain (not in a sub-domain)

                      int ok_trouve_loc = 0;
                      // Check if elem_test is in the current sub-domain
                      // (1 if real element, 2 if virtual, 0 otherwise)
                      // ok_trouve_loc set to 1 if elem_test is in the current sub-domain, 0 otherwise
                      for (int poly=0; poly<size_sz; poly++)
                        {
                          if (elem_test==sz(poly))
                            {
                              ok_trouve_loc = 1;
                              if (elem_test>nb_elem-1)
                                ok_trouve_loc = 2;
                              break;
                            }
                        }

                      // elem_test is in the current sub-domain (as a real element) or is a boundary condition
                      if ((ok_trouve_loc==1) || (elem_test==-1))
                        fac_pond = 1.;
                      // elem_test is in another sub-domain or in the default part
                      // or in the current sub-domain but as a virtual element
                      else
                        fac_pond = 0.5;

                      if (face_g<prem_face_std)
                        vol_entrelace = vol_entrelaces_Cl(face_g);
                      else
                        vol_entrelace =  vol_entrelaces(face_g);

                      vol_glob_pond(cpt) += fac_pond*vol_entrelace*por_face(face_g);
                      face_fait(face_g) = 1;
                    }
                }
            }
        }

      // Remove the vol_entrelaces_Cl contributions for Dirichlet faces because the power
      // attributed to these vol_entrelaces_Cl is not effective during the computation
      for (int n_bord=0; n_bord<zvef.nb_front_Cl(); n_bord++)
        {
          const Cond_lim& la_cl = zclvef.les_conditions_limites(n_bord);
          if (sub_type(Dirichlet,la_cl.valeur()))
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
              int face;
              int num1 = 0;
              int num2 = le_bord.nb_faces();
              for (int ind_face=num1; ind_face<num2; ind_face++)
                {
                  face = le_bord.num_face(ind_face);
                  if (face_fait(face)==1)
                    vol_glob_pond(cpt) -= vol_entrelaces_Cl(face)*por_face(face);
                }
            }
        }

      cpt++;
    }

  for (int num_face=0; num_face<prem_face_std; num_face++)
    {
      face_marq = faces_doubles[num_face];
      double contrib_double = double(face_marq);
      vol_glob_pond(0) += (1.-0.5*contrib_double)*vol_entrelaces_Cl(num_face)*por_face(num_face);
    }


  // Remove the vol_entrelaces_Cl contributions for Dirichlet faces because the power
  // attributed to these vol_entrelaces_Cl is not effective during the computation
  for (int n_bord=0; n_bord<zvef.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zclvef.les_conditions_limites(n_bord);
      if (sub_type(Dirichlet,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int face;
          int num1 = 0;
          int num2 = le_bord.nb_faces();
          for (int ind_face=num1; ind_face<num2; ind_face++)
            {
              face = le_bord.num_face(ind_face);
              face_marq = faces_doubles[face];
              double contrib_double = double(face_marq);
              vol_glob_pond(0) -= (1.-0.5*contrib_double)*vol_entrelaces_Cl(face)*por_face(face);
            }
        }
    }

  for (int num_face=prem_face_std; num_face<nb_faces; num_face++)
    {
      face_marq = faces_doubles[num_face];
      double contrib_double = double(face_marq);
      vol_glob_pond(0) += (1.-0.5*contrib_double)*vol_entrelaces(num_face)*por_face(num_face);
    }
  vol_glob_pond(0) = mp_sum(vol_glob_pond(0));

  for (int i=1; i<size_vol; i++)
    {
      vol_glob_pond(i) = mp_sum(vol_glob_pond(i));
      vol_glob_pond(0) -= vol_glob_pond(i);
    }

  return vol_glob_pond;
}
