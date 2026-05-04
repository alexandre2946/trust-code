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

#include <Op_Conv_Coloc_base.h>
#include <Domaine_Cl_Coloc.h>
#include <Champ_Inc_base.h>
#include <Momentum_Euler.h>
#include <Domaine_Coloc.h>
#include <EcrFicPartage.h>
#include <Pb_Euler.h>

Implemente_base(Op_Conv_Coloc_base,"Op_Conv_Coloc_base",Operateur_Conv_base);

Sortie& Op_Conv_Coloc_base::printOn(Sortie& os) const { return Operateur_Conv_base::printOn(os); }
Entree& Op_Conv_Coloc_base::readOn(Entree& is) { Operateur_Conv_base::readOn(is);  return is; }

void Op_Conv_Coloc_base::completer()
{
  if (!sub_type(Pb_Euler, equation().probleme()))
    {
      Cerr << "WHAT !! Operator " << que_suis_je() << " is only available for Pb_Euler not " << equation().probleme().que_suis_je() << " !! " << finl;
      Process::exit();
    }

  Operateur_base::completer();
  assert(le_dom_coloc_);
}

void Op_Conv_Coloc_base::associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& inc)
{
  le_dom_coloc_ = ref_cast(Domaine_Coloc, domaine_dis);
  le_dcl_coloc_ = ref_cast(Domaine_Cl_Coloc, zcl);
  le_champ_inco = ref_cast(Champ_Inc_base,inc);
}

void Op_Conv_Coloc_base::associer_domaine_cl_dis(const Domaine_Cl_dis_base& zcl)
{
  le_dcl_coloc_ = ref_cast(Domaine_Cl_Coloc, zcl);
}

void Op_Conv_Coloc_base::ajouter_blocs(matrices_t mats, DoubleTab& secmem, const tabs_t& semi_impl) const
{

  const Domaine_Coloc& domaine = ref_cast(Domaine_Coloc, le_dom_coloc_.valeur());
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();
  const int n = secmem.line_size();
  const int nb_faces = domaine.nb_faces();
  DoubleTrav num_flux(nb_faces, n);
  Riemann_solver(num_flux);

  for (int f = 0; f < nb_faces; f++)
    for (int i = 0; i < 2; i++)
      {
        const int e = f_e(f, i);
        if (e >= 0 && e < domaine.nb_elem())
          {
            for (int k = 0; k < n; k++)
              secmem(e, k) -= (i ? -1 : 1) * num_flux(f, k) * fs(f);
          }
      }
}

int Op_Conv_Coloc_base::impr(Sortie& os) const
{
  const Domaine& mon_dom = le_dom_coloc_->domaine();
  const int impr_mom = mon_dom.moments_a_imprimer() && sub_type(Momentum_Euler, equation());
  const int impr_sum = (mon_dom.bords_a_imprimer_sum().est_vide() ? 0 : 1);
  const int impr_bord = (mon_dom.bords_a_imprimer().est_vide() ? 0 : 1);
  const Schema_Temps_base& sch = le_dcl_coloc_->equation().probleme().schema_temps();
  DoubleTab& tab_flux_bords = flux_bords();
  int nb_comp = tab_flux_bords.nb_dim() > 1 ? tab_flux_bords.dimension(1) : 0;
  DoubleVect bilan(nb_comp);
  DoubleTab xgr;
  if (impr_mom)
    xgr = le_dom_coloc_->calculer_xgr();
  if (nb_comp == 0)
    return 1;
  int k, face;
  int nb_front_Cl = le_dom_coloc_->nb_front_Cl();
  DoubleTrav flux_bords2(5, nb_front_Cl, nb_comp);
  flux_bords2 = 0;
  for (int num_cl = 0; num_cl < nb_front_Cl; num_cl++)
    {
      const Cond_lim& la_cl = le_dcl_coloc_->les_conditions_limites(num_cl);
      const Front_VF& frontiere_dis = ref_cast(Front_VF, la_cl->frontiere_dis());
      int ndeb = frontiere_dis.num_premiere_face();
      int nfin = ndeb + frontiere_dis.nb_faces();
      for (face = ndeb; face < nfin; face++)
        {
          for (k = 0; k < nb_comp; k++)
            {
              flux_bords2(0, num_cl, k) += tab_flux_bords(face, k);
              if (mon_dom.bords_a_imprimer_sum().contient(frontiere_dis.le_nom()))
                flux_bords2(3, num_cl, k) += tab_flux_bords(face, k);
            } /* fin for k */
          if (impr_mom)
            {
              if (dimension == 2)
                {
                  flux_bords2(4, num_cl, 0) += tab_flux_bords(face, 1) * xgr(face, 0) - tab_flux_bords(face, 0) * xgr(face, 1);
                }
              else
                {
                  flux_bords2(4, num_cl, 0) += tab_flux_bords(face, 2) * xgr(face, 1) - tab_flux_bords(face, 1) * xgr(face, 2);
                  flux_bords2(4, num_cl, 1) += tab_flux_bords(face, 0) * xgr(face, 2) - tab_flux_bords(face, 2) * xgr(face, 0);
                  flux_bords2(4, num_cl, 2) += tab_flux_bords(face, 1) * xgr(face, 0) - tab_flux_bords(face, 0) * xgr(face, 1);
                }
            }
        } /* fin for face */
    }
  mp_sum_for_each_item(flux_bords2);

  if (je_suis_maitre())
    {
      ouvrir_fichier(Flux, "", 1);
      ouvrir_fichier(Flux_moment, "moment", impr_mom);
      ouvrir_fichier(Flux_sum, "sum", impr_sum);
      Flux.add_col(sch.temps_courant());
      if (impr_mom)
        Flux_moment.add_col(sch.temps_courant());
      if (impr_sum)
        Flux_sum.add_col(sch.temps_courant());
      for (int num_cl = 0; num_cl < nb_front_Cl; num_cl++)
        {
          for (k = 0; k < nb_comp; k++)
            {
              Flux.add_col(flux_bords2(0, num_cl, k));
              if (impr_sum)
                Flux_sum.add_col(flux_bords2(3, num_cl, k));
              bilan(k) += flux_bords2(0, num_cl, k);
            }
          if (dimension == 3)
            {
              for (k = 0; k < nb_comp; k++)
                if (impr_mom)
                  Flux_moment.add_col(flux_bords2(4, num_cl, k));
            }
          else
            {
              if (impr_mom)
                Flux_moment.add_col(flux_bords2(4, num_cl, 0));
            }
        } /* fin for num_cl */
      for (k = 0; k < nb_comp; k++)
        Flux.add_col(bilan(k));
      Flux << finl;
      if (impr_sum)
        Flux_sum << finl;
      if (impr_mom)
        Flux_moment << finl;
    }

  const LIST(Nom) &Liste_bords_a_imprimer = le_dom_coloc_->domaine().bords_a_imprimer();
  if (!Liste_bords_a_imprimer.est_vide())
    {
      EcrFicPartage Flux_face;
      ouvrir_fichier_partage(Flux_face, "", impr_bord);
      for (int num_cl = 0; num_cl < nb_front_Cl; num_cl++)
        {
          const Frontiere_dis_base& la_fr = le_dcl_coloc_->les_conditions_limites(num_cl)->frontiere_dis();
          const Cond_lim& la_cl = le_dcl_coloc_->les_conditions_limites(num_cl);
          const Front_VF& frontiere_dis = ref_cast(Front_VF, la_cl->frontiere_dis());
          int ndeb = frontiere_dis.num_premiere_face();
          int nfin = ndeb + frontiere_dis.nb_faces();
          if (mon_dom.bords_a_imprimer().contient(la_fr.le_nom()))
            {
              if (je_suis_maitre())
                {
                  Flux_face << "# Flux par face sur " << la_fr.le_nom() << " au temps ";
                  sch.imprimer_temps_courant(Flux_face);
                  Flux_face << " : " << finl;
                }
              for (face = ndeb; face < nfin; face++)
                {
                  if (dimension == 2)
                    Flux_face << "# Face a x= " << le_dom_coloc_->xv(face, 0) << " y= " << le_dom_coloc_->xv(face, 1) << " : ";
                  else if (dimension == 3)
                    Flux_face << "# Face a x= " << le_dom_coloc_->xv(face, 0) << " y= " << le_dom_coloc_->xv(face, 1) << " z= " << le_dom_coloc_->xv(face, 2) << " : ";
                  for (k = 0; k < nb_comp; k++)
                    Flux_face << tab_flux_bords(face, k) << " ";
                  Flux_face << finl;
                }
              Flux_face.syncfile();
            }
        }
    }

  return 1;
}
