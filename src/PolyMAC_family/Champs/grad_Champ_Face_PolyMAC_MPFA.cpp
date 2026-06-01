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

#include <grad_Champ_Face_PolyMAC_MPFA.h>
#include <Champ_Fonc_Elem_PolyMAC_CDO.h>
#include <Frottement_impose_base.h>
#include <Champ_Face_PolyMAC_MPFA.h>


#include <Dirichlet.h>
#include <EChaine.h>

Implemente_instanciable(grad_Champ_Face_PolyMAC_MPFA, "grad_Champ_Face_PolyMAC_MPFA", Champ_Fonc_Face_PolyMAC_CDO);

Sortie& grad_Champ_Face_PolyMAC_MPFA::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& grad_Champ_Face_PolyMAC_MPFA::readOn(Entree& s) { return s; }

int grad_Champ_Face_PolyMAC_MPFA::fixer_nb_valeurs_nodales(int n)
{
  const Champ_Fonc_base& self = ref_cast(Champ_Fonc_base, *this);
  const Domaine_PolyMAC_MPFA& domaine = ref_cast(Domaine_PolyMAC_MPFA, self.domaine_dis_base());

  assert(n < 0);

  const MD_Vector& md = domaine.mdv_ch_face;
  // Problem: nb_comp equals 2 but we only want one dimension !!!
  // HACK:
  int old_nb_compo = nb_compo_;
  nb_compo_ /= dimension;
  creer_tableau_distribue(md);
  nb_compo_ = old_nb_compo;
  return n;
}

void grad_Champ_Face_PolyMAC_MPFA::mettre_a_jour(double tps)
{
  if (temps() != tps)
    me_calculer(tps);
  Champ_Fonc_base::mettre_a_jour(tps);
}

void grad_Champ_Face_PolyMAC_MPFA::init_grad()
{
  is_init = 1;
  return;
}

void grad_Champ_Face_PolyMAC_MPFA::me_calculer(double tps)
{
  if (!is_init)
    init_grad();
  update_tab_grad(0); // compute the fgrad coefficients required for the field computation at faces
  calc_gradfve();
  update_ge(); // compute the field at elements from the field at faces
  valeurs().echange_espace_virtuel();
}

void grad_Champ_Face_PolyMAC_MPFA::update_tab_grad(int full_stencil)
{
  const IntTab& f_cl = champ_a_deriver().fcl();
  const Domaine_PolyMAC_MPFA& domaine = ref_cast(Domaine_PolyMAC_MPFA, domaine_vf());
  const Conds_lim& cls = champ_a_deriver().domaine_Cl_dis().les_conditions_limites(); // boundary conditions for the field to differentiate

  domaine.fgrad(champ_a_deriver().valeurs().line_size(), 0, cls, f_cl, nullptr, nullptr, 1, full_stencil, gradve_d, gradve_e, gradve_w);
}

void grad_Champ_Face_PolyMAC_MPFA::calc_gradfve()
{
  const Domaine_PolyMAC_MPFA& domaine = ref_cast(Domaine_PolyMAC_MPFA, domaine_vf());
  const Champ_Face_PolyMAC_MPFA& ch = ref_cast(Champ_Face_PolyMAC_MPFA, champ_a_deriver());

  const IntTab& fcl = champ_a_deriver().fcl();
  const Conds_lim& cls = ch.domaine_Cl_dis().les_conditions_limites(); // boundary conditions for the field to differentiate

  /*  const IntTab&                   fcl = fcl_g;
   const Conds_lim&                cls = cls_g;*/
  int f, n;
  int d_U; //velocity coordinate
  int D = dimension;
  int N = champ_a_deriver().valeurs().line_size(); //number of phases
  int ne_tot = domaine.nb_elem_tot(), nf_tot = domaine.nb_faces_tot(), nf = domaine.nb_faces();

  const DoubleTab& tab_ch = ch.passe();
  DoubleTab& val = valeurs();

  for (f = 0; f < nf; f++)
    {
      for (d_U = 0; d_U < D; d_U++)
        for (n = 0; n < N; n++) // velocity coordinates and phase
          {
            val(f, d_U + n * D) = 0;
            for (int j = gradve_d(f); j < gradve_d(f + 1); j++)
              {
                int e = gradve_e(j);
                int f_bord;
                if (e < ne_tot) //contribution from an element
                  {
                    double val_e = tab_ch(nf_tot + d_U + e * D, n);
                    val(f, d_U + n * D) += gradve_w(j, n) * val_e;
                  }
                else if (fcl(f_bord = e - ne_tot, 0) == 3) //contribution from a boundary: only Dirichlet contributes
                  {
                    double val_f_bord = ref_cast(Dirichlet, cls[fcl(f_bord, 1)].valeur()).val_imp(fcl(f_bord, 2), N * d_U + n);
                    val(f, d_U + n * D) += gradve_w(j, n) * val_f_bord;
                  }
              }
          }
    }

}

void grad_Champ_Face_PolyMAC_MPFA::update_ge()
{
  DoubleTab& val = valeurs();
  const Domaine_PolyMAC_MPFA& domaine = ref_cast(Domaine_PolyMAC_MPFA, domaine_vf());
  const DoubleVect& ve = domaine.volumes(), &fs = domaine.face_surfaces();
  const IntTab& e_f = domaine.elem_faces(), &f_e = domaine.face_voisins();
  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv();
  int e, f, j, d, D = dimension, n, N = val.line_size(), ne_tot = domaine.nb_elem_tot(), nf_tot = domaine.nb_faces_tot();
  double fac;

  for (e = 0; e < ne_tot; e++)
    {
      for (d = 0; d < D; d++)
        for (n = 0; n < N; n++)
          val(nf_tot + D * e + d, n) = 0;
      for (j = 0; j < e_f.dimension(1) && (f = e_f(e, j)) >= 0; j++)
        for (fac = (e == f_e(f, 0) ? 1 : -1) * fs(f) / ve(e), d = 0; d < D; d++)
          for (n = 0; n < N; n++)
            val(nf_tot + D * e + d, n) += fac * (xv(f, d) - xp(e, d)) * val(f, n);
    }
}

void grad_Champ_Face_PolyMAC_MPFA::init_ge2() const
{
}

void grad_Champ_Face_PolyMAC_MPFA::update_ge2(DoubleTab& val, int incr) const
{
}
