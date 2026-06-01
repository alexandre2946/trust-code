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

#include <Op_Diff_VDF_Face_Axi_base.h>

Implemente_base(Op_Diff_VDF_Face_Axi_base,"Op_Diff_VDF_Face_Axi_base",Op_Diff_VDF_Face_base);

Sortie& Op_Diff_VDF_Face_Axi_base::printOn(Sortie& s ) const { return s << que_suis_je() ; }
Entree& Op_Diff_VDF_Face_Axi_base::readOn(Entree& s ) { return s ; }

void Op_Diff_VDF_Face_Axi_base::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_cl_dis, const Champ_Inc_base& ch_transporte)
{
  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_dis);
  const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF,domaine_cl_dis);
  const Champ_Face_VDF& inco = ref_cast(Champ_Face_VDF,ch_transporte);
  le_dom_vdf = zvdf;
  la_zcl_vdf = zclvdf;
  inconnue = inco;
  surface.ref(zvdf.face_surfaces());
  volumes_entrelaces.ref(zvdf.volumes_entrelaces());
  orientation.ref(zvdf.orientation());
  porosite.ref(la_zcl_vdf->equation().milieu().porosite_face());
  xp.ref(zvdf.xp());
  xv.ref(zvdf.xv());
  Qdm.ref(zvdf.Qdm());
  face_voisins.ref(zvdf.face_voisins());
  elem_faces.ref(zvdf.elem_faces());
  type_arete_bord.ref(zclvdf.type_arete_bord());
}

double Op_Diff_VDF_Face_Axi_base::calculer_dt_stab() const
{
  return Op_Diff_VDF_base::calculer_dt_stab_(le_dom_vdf.valeur()) ;
}

void Op_Diff_VDF_Face_Axi_base::ajouter_elem(const DoubleTab& inco, DoubleTab& resu) const
{
  if (inco.line_size() > 1) not_implemented(__func__);
  for (int num_elem = 0; num_elem < le_dom_vdf->nb_elem(); num_elem++)
    {
      const int fx0 = elem_faces(num_elem,0), fx1 = elem_faces(num_elem,dimension), fy0 = elem_faces(num_elem,1), fy1 = elem_faces(num_elem,1+dimension);
      // Compute tau11
      const double tau11 = (inco[fx1]-inco[fx0])/(xv(fx1,0) - xv(fx0,0));
      // Compute tau22
      double R = xp(num_elem,0), d_teta = xv(fy1,1) - xv(fy0,1);
      if (d_teta < 0) d_teta += deux_pi;
      double tau22 =  (inco[fy1]-inco[fy0])/(R*d_teta);
      tau22 += 0.5*(inco[fx0]+inco[fx1])/R; // extra terms in axisymmetric

      const double flux_X = tau11*nu_(num_elem)*0.5*(surface(fx0)+surface(fx1)), flux_Y = tau22*nu_(num_elem)*0.5*(surface(fy0)+surface(fy1));
      resu[fx0] += flux_X;
      resu[fx1] -= flux_X;
      resu[fy0] += flux_Y;
      resu[fy1] -= flux_Y;

      // Extra terms in the axisymmetric Laplacian: integrated as source terms
      const double coef_laplacien_axi = +0.5*tau22*nu_(num_elem);
      resu[fx0] -= coef_laplacien_axi*volumes_entrelaces(fx0)*porosite(fx0)/xv(fx0,0);
      resu[fx1] -= coef_laplacien_axi*volumes_entrelaces(fx1)*porosite(fx1)/xv(fx1,0);
    }
}

void Op_Diff_VDF_Face_Axi_base::ajouter_elem_3D(const DoubleTab& inco, DoubleTab& resu) const
{
  for (int num_elem = 0; num_elem < le_dom_vdf->nb_elem(); num_elem++)
    {
      const int fz0 = elem_faces(num_elem,2), fz1 = elem_faces(num_elem,2+dimension);
      // Compute tau33
      const double tau33 = (inco[fz1]-inco[fz0])/(xv(fz1,2) - xv(fz0,2)), flux_Z = tau33*nu_(num_elem)*0.5*(surface(fz0)+surface(fz1));
      resu[fz0] += flux_Z;
      resu[fz1] -= flux_Z;
    }
}

void  Op_Diff_VDF_Face_Axi_base::ajouter_aretes_bords(const DoubleTab& inco, DoubleTab& resu) const
{
  int ndeb = le_dom_vdf->premiere_arete_bord(), nfin = ndeb + le_dom_vdf->nb_aretes_bord();
  for (int n_arete = ndeb; n_arete < nfin; n_arete++)
    {
      const int n_type = type_arete_bord(n_arete-ndeb);

      switch(n_type)
        {
        case TypeAreteBordVDF::PAROI_PAROI: /* fall through */
        case TypeAreteBordVDF::FLUIDE_FLUIDE:
        case TypeAreteBordVDF::PAROI_FLUIDE:
          {
            const int fac1 = Qdm(n_arete,0), fac2 = Qdm(n_arete,1), fac3 = Qdm(n_arete,2), signe  = Qdm(n_arete,3), ori1 = orientation(fac1), ori3 = orientation(fac3);
            const int rang1 = fac1 - le_dom_vdf->premiere_face_bord(), rang2 = fac2 - le_dom_vdf->premiere_face_bord();
            double vit_imp, dist3, tps = inconnue->temps();

            if (n_type == TypeAreteBordVDF::PAROI_FLUIDE) // wall_fluid edge: we must determine which is the fluid face
              {
                if (est_egal(inco[fac1],0)) vit_imp = Champ_Face_get_val_imp_face_bord(tps,rang2,ori3,la_zcl_vdf.valeur());
                else vit_imp = Champ_Face_get_val_imp_face_bord(tps,rang1,ori3,la_zcl_vdf.valeur());
              }
            else vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(tps,rang1,ori3,la_zcl_vdf.valeur())+Champ_Face_get_val_imp_face_bord(tps,rang2,ori3,la_zcl_vdf.valeur()));

            const double db_diffusivite = nu_mean_2_pts_(face_voisins(fac3,0),face_voisins(fac3,1));

            if (ori1 == 0) // boundary with equation R = const
              {
                double flux1;
                if (ori3 == 1)  // flux of tau12 through the boundary
                  {
                    dist3 = xv(fac3,0)-xv(fac1,0);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;

                    const double tau12 = (inco[fac3]-vit_imp)/dist3;
                    flux1 = db_diffusivite*tau12*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                else //if (ori3 == 2)  flux of tau13 through the boundary
                  {
                    assert(ori3 == 2);
                    dist3 = xv(fac3,0)-xv(fac1,0);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau13 = (inco[fac3]-vit_imp)/dist3;
                    flux1 = db_diffusivite*tau13*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                resu[fac3] += signe*flux1;
              }
            else if (ori1 == 1) // boundary with equation theta = const
              {
                double R = xv(fac3,0), d_teta = xv(fac3,1) - xv(fac1,1);
                if (d_teta < 0) d_teta += deux_pi;

                dist3  = R*d_teta;
                if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                if (ori3 == 0) // flux of tau21 through the boundary
                  {
                    double tau21 = (inco[fac3]-vit_imp)/dist3;
                    tau21 -= 0.5*(inco[fac1]+inco[fac2])/R; // Extra term in axi
                    const double flux2 = db_diffusivite*tau21*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                    resu[fac3]+=signe*flux2;

                    // Extra terms in the axisymmetric Laplacian
                    // They are integrated as source terms
                    const double coef_laplacien_axi = 0.5*db_diffusivite*tau21;
                    resu(fac1) += coef_laplacien_axi*volumes_entrelaces(fac1)*porosite(fac1)/xv(fac1,0);
                    resu(fac2) += coef_laplacien_axi*volumes_entrelaces(fac2)*porosite(fac2)/xv(fac2,0);
                  }
                else if (ori3 == 2) // flux of tau23 through the boundary
                  {
                    // XXX : beware if discrepancy: in the constant case it was tau23 = signe*(vit_imp-inco[fac3])/dist3 (normally the same but just in case)
                    const double tau23 = (inco[fac3]-vit_imp)/dist3;
                    const double flux3 = db_diffusivite*tau23*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                    resu[fac3] += signe*flux3;
                  }
              }
            else // (ori1 == 2) boundary with equation Z = const
              {
                double flux4;
                if (ori3 == 0) // flux of tau31 through the boundary
                  {
                    dist3 = xv(fac3,2)-xv(fac1,2);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;

                    const double tau31 = (inco[fac3]-vit_imp)/dist3;
                    flux4 = db_diffusivite*tau31*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                else // if (ori3 == 1)  flux of tau32 through the boundary
                  {
                    assert(ori3 == 1) ;
                    dist3 = xv(fac3,2)-xv(fac1,2);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;

                    const double tau32 = (inco[fac3]-vit_imp)/dist3;
                    flux4 = db_diffusivite*tau32*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                resu[fac3] += signe*flux4;
              }
            break;
          }
        case TypeAreteBordVDF::NAVIER_NAVIER: // pas de flux diffusif calcule
          break;
        default :
          {
            Cerr << "On a rencontre un type d'arete non prevu : [ num arete : " << n_arete << " ], [ type : " << n_type << " ]" << finl;
            Process::exit();
            break;
          }
        }
    }
}

void Op_Diff_VDF_Face_Axi_base::ajouter_aretes_mixtes_internes(const DoubleTab& inco, DoubleTab& resu) const
{
  const int ndeb = le_dom_vdf->premiere_arete_mixte(), nfin = le_dom_vdf->nb_aretes();
  for (int n_arete=ndeb; n_arete<nfin; n_arete++)
    {
      const int fac1 = Qdm(n_arete,0), fac2 = Qdm(n_arete,1), fac3 = Qdm(n_arete,2), fac4 = Qdm(n_arete,3), ori1 = orientation(fac1), ori3 = orientation(fac3);
      const double db_diffusivite = nu_mean_4_pts_(fac3,fac4);

      if (ori1 == 1)  // (only possibility: ori3 =0)  XY edge
        {
          double flux1;
          // Compute tau21
          const double R = xv(fac3,0);
          double d_teta = xv(fac4,1) - xv(fac3,1);
          if (d_teta < 0) d_teta += deux_pi;
          double tau21 = (inco(fac4)-inco(fac3))/(R*d_teta);

          // Extra term in axi
          tau21 -= 0.5*(inco[fac1]+inco[fac2])/R;

          // flux of tau21 on the facet straddling faces fac1 and fac2
          flux1 = db_diffusivite*tau21*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
          resu(fac3) += flux1;
          resu(fac4) -= flux1;

          // Extra terms in the axisymmetric Laplacian: integrated as source terms
          const double coef_laplacien_axi = 0.5*db_diffusivite*tau21;
          resu(fac1) += coef_laplacien_axi*volumes_entrelaces(fac1)*porosite(fac1)/xv(fac1,0);
          resu(fac2) += coef_laplacien_axi*volumes_entrelaces(fac2)*porosite(fac2)/xv(fac2,0);

          // Compute tau12
          const double tau12 = (inco(fac2)-inco(fac1))/(xv(fac2,0) - xv(fac1,0));
          // flux of tau12 on the facet straddling faces fac3 and fac4
          flux1 = db_diffusivite*tau12*0.25*(surface(fac3)+surface(fac4))*(porosite(fac3)+porosite(fac4));
          resu(fac1) += flux1;
          resu(fac2) -= flux1;
        }
      else if (ori3 == 1) // (only possibility: ori1 = 2) YZ edge
        {
          double flux2;
          // Compute tau32
          const double tau32 = (inco(fac4)-inco(fac3))/(xv(fac4,2) - xv(fac3,2));

          // flux of tau32 on the facet straddling faces fac1 and fac2
          flux2 = db_diffusivite*tau32*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
          resu(fac3) += flux2;
          resu(fac4) -= flux2;

          // Compute tau23
          const double R = xv(fac1,0);
          double d_teta = xv(fac2,1) - xv(fac1,1);
          if (d_teta < 0) d_teta += deux_pi;
          const double tau23 = (inco(fac2)-inco(fac1))/(R*d_teta);
          // flux of tau23 on the facet straddling faces fac3 and fac4
          flux2 = db_diffusivite*tau23*0.25*(surface(fac3)+surface(fac4))*(porosite(fac3)+porosite(fac4));
          resu(fac1) += flux2;
          resu(fac2) -= flux2;
        }
      else // only possibility: ori1 = 2 and ori3 = 0: XZ edge
        {
          double flux3;
          // Compute tau31
          const double tau31 = (inco(fac4)-inco(fac3))/(xv(fac4,2) - xv(fac3,2));

          // flux of tau31 on the facet straddling faces fac1 and fac2
          flux3 = db_diffusivite*tau31*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
          resu(fac3) += flux3;
          resu(fac4) -= flux3;

          // Compute tau13
          const double tau13 = (inco(fac2)-inco(fac1))/(xv(fac2,0) - xv(fac1,0));
          // flux of tau13 on the facet straddling faces fac3 and fac4
          flux3 = db_diffusivite*tau13*0.25*(surface(fac3)+surface(fac4))*(porosite(fac3)+porosite(fac4));
          resu(fac1) += flux3;
          resu(fac2) -= flux3;
        }
    }
}

DoubleTab& Op_Diff_VDF_Face_Axi_base::ajouter(const DoubleTab& inco,  DoubleTab& resu) const
{
  ajouter_elem(inco,resu); // Loop over elements
  if (dimension == 3) ajouter_elem_3D(inco,resu); // Loop over extra elements if 3D
  ajouter_aretes_bords(inco, resu); // Loop over boundary edges
  ajouter_aretes_mixtes_internes(inco, resu); // Loop over mixed and internal edges
  return resu;
}

DoubleTab& Op_Diff_VDF_Face_Axi_base::calculer(const DoubleTab& inco, DoubleTab& resu) const
{
  resu = 0;
  return ajouter(inco,resu);
}

void Op_Diff_VDF_Face_Axi_base::fill_coeff_matrice_morse(const int fac1, const int fac2, const double flux, Matrice_Morse& matrice) const
{
  const auto& tab1 = matrice.get_set_tab1();
  const auto& tab2 = matrice.get_set_tab2();
  auto& coeff = matrice.get_set_coeff();
  for (auto k = tab1[fac1]-1; k < tab1[fac1+1]-1; k++)
    {
      if (tab2[k]-1 == fac1) coeff[k] += flux;
      if (tab2[k]-1 == fac2) coeff[k] -= flux;
    }
  for (auto k = tab1[fac2]-1; k < tab1[fac2+1]-1; k++)
    {
      if (tab2[k]-1 == fac1) coeff[k] -= flux;
      if (tab2[k]-1 == fac2) coeff[k] += flux;
    }
}

void Op_Diff_VDF_Face_Axi_base::ajouter_contribution_elem(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  if (inco.line_size() > 1) not_implemented(__func__);

  const auto& tab1 = matrice.get_set_tab1();
  const auto& tab2 = matrice.get_set_tab2();
  auto& coeff = matrice.get_set_coeff();
  for (int num_elem = 0; num_elem < le_dom_vdf->nb_elem(); num_elem++)
    {
      const int fx0 = elem_faces(num_elem,0), fx1 = elem_faces(num_elem,dimension), fy0 = elem_faces(num_elem,1), fy1 = elem_faces(num_elem,1+dimension);
      // Compute tau11
      const double tau11 = 1/(xv(fx1,0) - xv(fx0,0));

      // Compute tau22
      const double R = xp(num_elem,0);
      double d_teta = xv(fy1,1) - xv(fy0,1);
      if (d_teta < 0) d_teta += deux_pi;

      double tau22 = 1/(R*d_teta);
      // extra terms in axi
      tau22 += 0.5/R;
      const double flux_X = tau11*nu_(num_elem)*0.5*(surface(fx0)+surface(fx1));
      const double flux_Y = tau22*nu_(num_elem)*0.5*(surface(fy0)+surface(fy1));
      fill_coeff_matrice_morse(fx0,fx1,flux_X,matrice);
      fill_coeff_matrice_morse(fy0,fy1,flux_Y,matrice);

      // Extra terms in the axisymmetric Laplacian: integrated as source terms
      const double coef_laplacien_axi = +0.5*tau22*nu_(num_elem);

      for (auto k = tab1[fx0]-1; k < tab1[fx0+1]-1; k++)
        if (tab2[k]-1 == fx0) coeff[k] += coef_laplacien_axi*volumes_entrelaces(fx0)*porosite(fx0)/xv(fx0,0);

      for (auto k=tab1[fx1]-1; k<tab1[fx1+1]-1; k++)
        if (tab2[k]-1 == fx1) coeff[k] += coef_laplacien_axi*volumes_entrelaces(fx1)*porosite(fx1)/xv(fx1,0);
    }
}

void Op_Diff_VDF_Face_Axi_base::ajouter_contribution_elem_3D(Matrice_Morse& matrice) const
{
  for (int num_elem = 0; num_elem < le_dom_vdf->nb_elem(); num_elem++)
    {
      const int fz0 = elem_faces(num_elem,2), fz1 = elem_faces(num_elem,2+dimension);
      // Compute tau33
      const double tau33 = 1/(xv(fz1,2) - xv(fz0,2));
      const double flux_Z = tau33*nu_(num_elem)*0.5*(surface(fz0)+surface(fz1));
      fill_coeff_matrice_morse(fz0,fz1,flux_Z,matrice);
    }
}

void Op_Diff_VDF_Face_Axi_base::ajouter_contribution_aretes_bords(Matrice_Morse& matrice) const
{
  const auto& tab1 = matrice.get_set_tab1();
  const auto& tab2 = matrice.get_set_tab2();
  auto& coeff = matrice.get_set_coeff();
  const int ndeb = le_dom_vdf->premiere_arete_bord(), nfin = ndeb + le_dom_vdf->nb_aretes_bord();
  for (int n_arete=ndeb; n_arete<nfin; n_arete++)
    {
      const int n_type = type_arete_bord(n_arete-ndeb);
      switch(n_type)
        {
        case TypeAreteBordVDF::PAROI_PAROI: /* fall through */
        case TypeAreteBordVDF::FLUIDE_FLUIDE:
        case TypeAreteBordVDF::PAROI_FLUIDE:
          {
            const int fac1 = Qdm(n_arete,0), fac2 = Qdm(n_arete,1), fac3 = Qdm(n_arete,2), signe  = Qdm(n_arete,3), ori1 = orientation(fac1), ori3 = orientation(fac3);
            const double db_diffusivite = nu_mean_2_pts_(face_voisins(fac3,0),face_voisins(fac3,1));

            if (ori1 == 0) // boundary with equation R = const
              {
                double flux1;
                if (ori3 == 1)  // flux of tau12 through the boundary
                  {
                    double dist3 = xv(fac3,0)-xv(fac1,0);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau12 = 1/dist3;
                    flux1 = db_diffusivite*tau12*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                else //if (ori3 == 2)  flux of tau13 through the boundary
                  {
                    assert (ori3 == 2);
                    double dist3 = xv(fac3,0)-xv(fac1,0);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau13 =  1/dist3;
                    flux1 = db_diffusivite*tau13*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }

                for (auto k = tab1[fac3]-1; k < tab1[fac3+1]-1; k++)
                  if (tab2[k]-1 == fac3) coeff[k] += signe*flux1;
              }
            else if (ori1 == 1) // boundary with equation theta = const
              {
                const double R = xv(fac3,0);
                double d_teta = xv(fac3,1) - xv(fac1,1);
                if (d_teta < 0) d_teta += deux_pi;

                double dist3  = R*d_teta;
                if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;

                if (ori3 == 0) // flux of tau21 through the boundary
                  {
                    double tau21 = 1/dist3;

                    // Extra term in axi
                    tau21 -= 0.5/R;

                    const double flux2 = db_diffusivite*tau21*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));

                    for (auto k = tab1[fac3]-1; k < tab1[fac3+1]-1; k++)
                      if (tab2[k]-1 == fac3) coeff[k] += signe*flux2;

                    // Extra terms in the axisymmetric Laplacian: integrated as source terms
                    const double coef_laplacien_axi = 0.5*db_diffusivite*tau21;

                    for (auto k = tab1[fac1]-1; k < tab1[fac1+1]-1; k++)
                      if (tab2[k]-1 == fac1) coeff[k] += coef_laplacien_axi*volumes_entrelaces(fac1)*porosite(fac1)/xv(fac1,0);

                    for (auto k = tab1[fac2]-1; k < tab1[fac2+1]-1; k++)
                      if (tab2[k]-1 == fac2) coeff[k] += coef_laplacien_axi*volumes_entrelaces(fac2)*porosite(fac2)/xv(fac2,0);
                  }
                else // if (ori3 == 2) flux of tau23 through the boundary
                  {
                    assert(ori3 == 2);
                    const double tau23 = 1/dist3;
                    const double flux3 = db_diffusivite*tau23*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                    for (auto k = tab1[fac3]-1; k < tab1[fac3+1]-1; k++)
                      if (tab2[k]-1 == fac3) coeff[k] += signe*flux3;
                  }
              }
            else // (ori1 == 2) boundary with equation Z = const
              {
                double flux4;
                if (ori3 == 0) // flux of tau31 through the boundary
                  {
                    double dist3 = xv(fac3,2)-xv(fac1,2);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau31 = 1/dist3;
                    flux4 = db_diffusivite*tau31*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                else //if (ori3 == 1) flux of tau32 through the boundary
                  {
                    assert(ori3 == 1);
                    double dist3 = xv(fac3,2)-xv(fac1,2);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau32 = 1/dist3;
                    flux4 = db_diffusivite*tau32*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }

                for (auto k = tab1[fac3]-1; k < tab1[fac3+1]-1; k++)
                  if (tab2[k]-1 == fac3) coeff[k] += signe*flux4;
              }
            break;
          }
        case TypeAreteBordVDF::NAVIER_NAVIER: // pas de flux diffusif calcule
          break;
        default :
          {
            Cerr << "On a rencontre un type d'arete non prevu : [ num arete : " << n_arete << " ], [ type : " << n_type << " ]" << finl;
            Process::exit();
            break;
          }
        }
    }
}

void Op_Diff_VDF_Face_Axi_base::ajouter_contribution_aretes_mixtes_internes(Matrice_Morse& matrice) const
{
  const auto& tab1 = matrice.get_set_tab1();
  const auto& tab2 = matrice.get_set_tab2();
  auto& coeff = matrice.get_set_coeff();
  const int ndeb = le_dom_vdf->premiere_arete_mixte(), nfin = le_dom_vdf->nb_aretes();
  for (int n_arete = ndeb; n_arete < nfin; n_arete++)
    {
      const int fac1 = Qdm(n_arete,0), fac2 = Qdm(n_arete,1), fac3 = Qdm(n_arete,2), fac4 = Qdm(n_arete,3), ori1 = orientation(fac1), ori3 = orientation(fac3);
      const double db_diffusivite = nu_mean_4_pts_(fac3,fac4);

      if (ori1 == 1)  // (only possibility: ori3 =0)  XY edge
        {
          double flux1;
          // Compute tau21
          const double R = xv(fac3,0);
          double d_teta = xv(fac4,1) - xv(fac3,1);
          if (d_teta < 0) d_teta += deux_pi;
          double tau21 = 1/(R*d_teta);
          // Extra term in axi
          tau21 -= 0.5/R;

          // flux of tau21 on the facet straddling faces fac1 and fac2
          flux1 = db_diffusivite*tau21*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
          fill_coeff_matrice_morse(fac3,fac4,flux1,matrice);

          // Extra terms in the axisymmetric Laplacian: integrated as source terms
          const double coef_laplacien_axi = 0.5*db_diffusivite*tau21;
          for (auto k = tab1[fac1]-1; k < tab1[fac1+1]-1; k++)
            if (tab2[k]-1 == fac1) coeff[k] += coef_laplacien_axi*volumes_entrelaces(fac1)*porosite(fac1)/xv(fac1,0);

          for (auto k = tab1[fac2]-1; k < tab1[fac2+1]-1; k++)
            if (tab2[k]-1 == fac2) coeff[k] += coef_laplacien_axi*volumes_entrelaces(fac2)*porosite(fac2)/xv(fac2,0);

          // Compute tau12
          const double tau12 = 1/(xv(fac2,0) - xv(fac1,0));

          // flux of tau12 on the facet straddling faces fac3 and fac4
          flux1 = db_diffusivite*tau12*0.25*(surface(fac3)+surface(fac4))*(porosite(fac3)+porosite(fac4));
          fill_coeff_matrice_morse(fac1,fac2,flux1,matrice);
        }
      else if (ori3 == 1) // (only possibility: ori1 = 2) YZ edge
        {
          double flux2;
          // Compute tau32
          const double tau32 = 1/(xv(fac4,2) - xv(fac3,2));

          // flux of tau32 on the facet straddling faces fac1 and fac2
          flux2 = db_diffusivite*tau32*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
          fill_coeff_matrice_morse(fac3,fac4,flux2,matrice);

          // Compute tau23
          const double R = xv(fac1,0);
          double d_teta = xv(fac2,1) - xv(fac1,1);
          if (d_teta < 0) d_teta += deux_pi;
          const double tau23 = 1/(R*d_teta);

          // flux of tau23 on the facet straddling faces fac3 and fac4
          flux2 = db_diffusivite*tau23*0.25*(surface(fac3)+surface(fac4))*(porosite(fac3)+porosite(fac4));
          fill_coeff_matrice_morse(fac1,fac2,flux2,matrice);
        }
      else // only possibility: ori1 = 2 and ori3 = 0: XZ edge
        {
          double flux3;
          // Compute tau31
          const double tau31 = 1/(xv(fac4,2) - xv(fac3,2));

          // flux of tau31 on the facet straddling faces fac1 and fac2
          flux3 = db_diffusivite*tau31*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
          fill_coeff_matrice_morse(fac3,fac4,flux3,matrice);

          // Compute tau13
          const double tau13 = 1/(xv(fac2,0) - xv(fac1,0));

          // flux of tau13 on the facet straddling faces fac3 and fac4
          flux3 = db_diffusivite*tau13*0.25*(surface(fac3)+surface(fac4))*(porosite(fac3)+porosite(fac4));
          fill_coeff_matrice_morse(fac1,fac2,flux3,matrice);
        }
    }
}

void Op_Diff_VDF_Face_Axi_base::ajouter_contribution(const DoubleTab& inco, Matrice_Morse& matrice ) const
{
  ajouter_contribution_elem(inco,matrice); // Loop over elements
  if (dimension == 3) ajouter_contribution_elem_3D(matrice); // Loop over extra elements if 3D
  ajouter_contribution_aretes_bords(matrice); // Loop over boundary edges
  ajouter_contribution_aretes_mixtes_internes(matrice); // Loop over mixed and internal edges
}

void Op_Diff_VDF_Face_Axi_base::contribue_au_second_membre(DoubleTab& resu) const
{
  const int ndeb = le_dom_vdf->premiere_arete_bord(), nfin = ndeb + le_dom_vdf->nb_aretes_bord();
  for (int n_arete = ndeb; n_arete < nfin; n_arete++)
    {
      const int n_type = type_arete_bord(n_arete-ndeb);
      switch(n_type)
        {
        case TypeAreteBordVDF::PAROI_PAROI: /* fall through */
        case TypeAreteBordVDF::FLUIDE_FLUIDE:
        case TypeAreteBordVDF::PAROI_FLUIDE:
          {
            const int fac1 = Qdm(n_arete,0), fac2 = Qdm(n_arete,1), fac3 = Qdm(n_arete,2), signe  = Qdm(n_arete,3);
            const int ori1 = orientation(fac1), ori3 = orientation(fac3), rang1 = fac1 - le_dom_vdf->premiere_face_bord(), rang2 = fac2 - le_dom_vdf->premiere_face_bord();
            double vit_imp, tps = inconnue->temps();

            if (n_type == TypeAreteBordVDF::PAROI_FLUIDE) // wall_fluid edge: we must determine which is the fluid face
              {
                if (est_egal(inconnue->valeurs()(fac1), 0))
                  vit_imp = Champ_Face_get_val_imp_face_bord(tps, rang2, ori3, la_zcl_vdf.valeur());
                else
                  vit_imp = Champ_Face_get_val_imp_face_bord(tps, rang1, ori3, la_zcl_vdf.valeur());
              }
            else vit_imp = 0.5*(Champ_Face_get_val_imp_face_bord(tps,rang1,ori3,la_zcl_vdf.valeur())+Champ_Face_get_val_imp_face_bord(tps,rang2,ori3,la_zcl_vdf.valeur()));

            const double db_diffusivite =  nu_mean_2_pts_(face_voisins(fac3,0),face_voisins(fac3,1));

            if (ori1 == 0) // boundary with equation R = const
              {
                double flux1;
                if (ori3 == 1)  // flux of tau12 through the boundary
                  {
                    double dist3 = xv(fac3,0)-xv(fac1,0);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau12 = (-vit_imp)/dist3;
                    flux1 = db_diffusivite*tau12*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                else// if (ori3 == 2)  flux de tau13 a travers le bord
                  {
                    assert(ori3 == 2);
                    double dist3 = xv(fac3,0)-xv(fac1,0);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau13 = (-vit_imp)/dist3;
                    flux1 = db_diffusivite*tau13*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                resu[fac3] += signe*flux1;
              }
            else if (ori1 == 1) // boundary with equation theta = const
              {
                const double R = xv(fac3,0);
                double d_teta = xv(fac3,1) - xv(fac1,1);
                if (d_teta < 0) d_teta += deux_pi;
                double dist3  = R*d_teta;
                if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                if (ori3 == 0) // flux of tau21 through the boundary
                  {
                    // Extra term in axi
                    const double tau21 = (-vit_imp)/dist3;
                    const double flux2 = db_diffusivite*tau21*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                    resu[fac3] += signe*flux2;

                    // Extra terms in the axisymmetric Laplacian: integrated as source terms
                    const double coef_laplacien_axi = 0.5*db_diffusivite*tau21;
                    resu(fac1) += coef_laplacien_axi*volumes_entrelaces(fac1)*porosite(fac1)/xv(fac1,0);
                    resu(fac2) += coef_laplacien_axi*volumes_entrelaces(fac2)*porosite(fac2)/xv(fac2,0);
                  }
                else if (ori3 == 2) // flux of tau23 through the boundary
                  {
                    const double tau23 = (-vit_imp)/dist3;
                    const double flux3 = db_diffusivite*tau23*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                    resu[fac3] += signe*flux3;
                  }
              }
            else // (ori1 == 2) boundary with equation Z = const
              {
                double flux4;
                if (ori3 == 0)  // flux of tau31 through the boundary
                  {
                    double dist3 = xv(fac3,2) - xv(fac1,2);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau31 = (-vit_imp)/dist3;
                    flux4 = db_diffusivite*tau31*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                else //if (ori3 == 1)  flux de tau32 a travers le bord
                  {
                    assert(ori3 == 1);
                    double dist3 = xv(fac3,2)-xv(fac1,2);
                    if (n_type != TypeAreteBordVDF::PAROI_PAROI) dist3 *= 1;
                    const double tau32 = (-vit_imp)/dist3;
                    flux4 = db_diffusivite*tau32*0.25*(surface(fac1)+surface(fac2))*(porosite(fac1)+porosite(fac2));
                  }
                resu[fac3] += signe*flux4;
              }
            break;
          }
        case TypeAreteBordVDF::NAVIER_NAVIER: // pas de flux diffusif calcule
          break;
        default :
          {
            Cerr << "On a rencontre un type d'arete non prevu : [ num arete : " << n_arete << " ], [ type : " << n_type << " ]" << finl;
            Process::exit();
            break;
          }
        }
    }
}
