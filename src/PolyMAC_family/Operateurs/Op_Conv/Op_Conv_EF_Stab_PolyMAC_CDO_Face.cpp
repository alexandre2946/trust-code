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

#include <Op_Conv_EF_Stab_PolyMAC_CDO_Face.h>
#include <Discretisation_base.h>
#include <Dirichlet_homogene.h>
#include <Champ_Face_PolyMAC_CDO.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Schema_Temps_base.h>
#include <Domaine_PolyMAC_CDO.h>
#include <Probleme_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <TRUSTLists.h>
#include <Dirichlet.h>
#include <Param.h>
#include <cmath>

Implemente_instanciable( Op_Conv_EF_Stab_PolyMAC_CDO_Face, "Op_Conv_EF_Stab_PolyMAC_CDO_Face", Op_Conv_PolyMAC_CDO_base );
Implemente_instanciable( Op_Conv_Amont_PolyMAC_CDO_Face, "Op_Conv_Amont_PolyMAC_CDO_Face", Op_Conv_EF_Stab_PolyMAC_CDO_Face );
Implemente_instanciable( Op_Conv_Centre_PolyMAC_CDO_Face, "Op_Conv_Centre_PolyMAC_CDO_Face", Op_Conv_EF_Stab_PolyMAC_CDO_Face );

// XD Op_Conv_EF_Stab_PolyMAC_CDO_Face interprete Op_Conv_EF_Stab_PolyMAC_CDO_Face BRACE Class
// XD_CONT Op_Conv_EF_Stab_PolyMAC_CDO_Face_PolyMAC_CDO

Sortie& Op_Conv_EF_Stab_PolyMAC_CDO_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_CDO_base::printOn(os); }
Sortie& Op_Conv_Amont_PolyMAC_CDO_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_CDO_base::printOn(os); }
Sortie& Op_Conv_Centre_PolyMAC_CDO_Face::printOn(Sortie& os) const { return Op_Conv_PolyMAC_CDO_base::printOn(os); }

Entree& Op_Conv_EF_Stab_PolyMAC_CDO_Face::readOn(Entree& is)
{
  Op_Conv_PolyMAC_CDO_base::readOn(is);
  Param param(que_suis_je());
  param.ajouter("alpha", &alpha_);            // XD_ADD_P double
  // XD_CONT parametre ajustant la stabilisation de 0 (schema centre) a 1 (schema amont)
  param.lire_avec_accolades_depuis(is);
  return is;
}

Entree& Op_Conv_Amont_PolyMAC_CDO_Face::readOn(Entree& is)
{
  alpha_ = 1.0;
  return Op_Conv_PolyMAC_CDO_base::readOn(is);
}

Entree& Op_Conv_Centre_PolyMAC_CDO_Face::readOn(Entree& is)
{
  alpha_ = 0.0;
  return Op_Conv_PolyMAC_CDO_base::readOn(is);
}

double Op_Conv_EF_Stab_PolyMAC_CDO_Face::calculer_dt_stab() const
{

  double dt = 1e10;
  const Domaine_Poly_base& domaine = le_dom_poly_.valeur();
  const DoubleVect& fs = domaine.face_surfaces(), &pf = equation().milieu().porosite_face(), &ve = domaine.volumes(), &pe = equation().milieu().porosite_elem();
  const DoubleTab& vit = vitesse_->valeurs();
  const IntTab& e_f = domaine.elem_faces(), &f_e = domaine.face_voisins();
  const int N = vit.line_size();
  DoubleTrav flux(N); //somme des flux pf * |f| * vf, volume minimal des mailles d'elements/faces affectes par ce flux

  for (int e = 0; e < domaine.nb_elem(); e++)
    {
      // Calcul du volume effectif de l'element
      const double vol = pe(e) * ve(e);
      flux = 0.;

      // Parcourt des faces associees a l'element
      for (int i = 0; i < e_f.dimension(1); i++)
        {
          int f = e_f(e, i);
          if (f < 0) continue; // face in-existante

          for (int n = 0; n < N; n++)
            {
              // Ajout du flux entrant pour la composante n : Seuls les flux entrants comptent
              double flux_f = pf(f) * fs(f) * std::max((e == f_e(f, 1) ? 1 : -1) * vit(f, n), 0.);
              flux(n) += flux_f;
            }
        }

      // Calcul du pas de temps pour chaque composante n
      for (int n = 0; n < N; n++)
        if (std::abs(flux(n)) > 1e-12)
          dt = std::min(dt, vol / flux(n));
    }

  return Process::mp_min(dt);
}

void Op_Conv_EF_Stab_PolyMAC_CDO_Face::completer()
{
  Op_Conv_PolyMAC_CDO_base::completer();

  /* au cas ou... */
  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  if (domaine.domaine().nb_joints() && domaine.domaine().joint(0).epaisseur() < 2)
    {
      Cerr << "Op_Conv_EF_Stab_PolyMAC_CDO_Face : largeur de joint insuffisante (minimum 2)!" << finl;
      Process::exit();
    }
  porosite_f.ref(mon_equation->milieu().porosite_face());
  porosite_e.ref(mon_equation->milieu().porosite_elem());
}

void Op_Conv_EF_Stab_PolyMAC_CDO_Face::dimensionner(Matrice_Morse& mat) const
{
  if (has_interface_blocs())
    {
      Operateur_base::dimensionner(mat);
      return;
    }

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const IntTab& e_f = domaine.elem_faces(), &f_e = domaine.face_voisins(), &equiv = domaine.equiv();
  int i, j, k, l, e, eb, f, fb, fd, fc;

  ch.fcl();

  Stencil stencil(0, 2);

  for (f = 0; f < domaine.nb_faces_tot(); f++)
    if (f_e(f, 0) >= 0 && (f_e(f, 1) >= 0 ||  ch.fcl()(f, 0) == 1 ||  ch.fcl()(f, 0) == 3))
      {
        for (i = 0; i < 2; i++)
          if ((e = f_e(f, i)) >= 0)
            {
              for (k = 0; k < e_f.dimension(1) && (fb = e_f(e, k)) >= 0; k++)
                if (fb < domaine.nb_faces() &&  ch.fcl()(fb, 0) < 2) //partie "faces"
                  {
                    if ((fc = equiv(f, i, k)) >= 0 || f_e(f, 1) < 0)
                      for (j = 0; j < 2; j++) //equivalence : face fd -> face fb
                        {
                          fd = (j == i ? fb : fc); //element/face sources
                          if (fd >= 0) stencil.append_line(fb,fd);
                        }
                    else for (j = 0; j < 2; j++)  //pas d'equivalence : n_f * operateur aux elements
                        {
                          for (eb = f_e(f, j), l = 0; l < e_f.dimension(1) && (fc = e_f(eb, l)) >= 0; l++)
                            stencil.append_line(fb,fc);
                        }
                  }
            }
      }

  // Tri et suppression des doublons
  tableau_trier_retirer_doublons(stencil);

  // Allocation de la matrice
  int taille = domaine.nb_faces_tot() + (dimension < 3 ? domaine.domaine().nb_som_tot() : domaine.domaine().nb_aretes_tot());
  Matrix_tools::allocate_morse_matrix(taille, taille, stencil, mat);
}

// ajoute la contribution de la convection au second membre resu
// renvoie resu

DoubleTab& Op_Conv_EF_Stab_PolyMAC_CDO_Face::ajouter(const DoubleTab& inco, DoubleTab& secmem) const
{
  if (has_interface_blocs())
    return Operateur_base::ajouter(inco, secmem);

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const Conds_lim& cls = la_zcl_poly_->les_conditions_limites();

  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces(), &equiv = domaine.equiv();
  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv(), &vfd = domaine.volumes_entrelaces_dir(), &vit = vitesse_->valeurs();
  const DoubleVect& fs = domaine.face_surfaces(), &ve = domaine.volumes(),  &pf = porosite_f, &pe = porosite_e;
  const DoubleTab& nf = domaine.face_normales();

  int i, j, k, l, e, eb, f, fb, fc, fd, m, n, N = inco.line_size(), d, D = dimension, comp = !incompressible_;
  double mult;

  assert(N == 1);
  DoubleTrav dfac(2, N, N);
  double sum_dfac[2];
  const int nb_faces_tot = domaine.nb_faces_tot();
  const int nb_faces = domaine.nb_faces();
  int nb_face_elem = e_f.dimension(1);
  for (f = 0; f < nb_faces_tot; f++)
    {
      const int elem0 = f_e(f, 0);
      const int elem1 = f_e(f, 1);
      if (elem0 >= 0 && (elem1 >= 0 || ch.fcl()(f, 0) == 1 || ch.fcl()(f, 0) == 3))
        {
          //masse : diagonale + masse ajoutee si correlation
          const double inv_masse = 1.0 / (std::fabs(vit[f]) > 1e-10 ? inco(f) / vit[f] : 1.0);
          for (i = 0, dfac = 0; i < 2; i++)
            {
              //contribution a dfac
              for (eb = f_e(f, i), n = 0; n < N; n++)
                {
                  for (m = 0; m < N; m++)
                    dfac(ch.fcl()(f, 0) == 1 ? 0 : i, n, m) += fs(f) * inco[f] * pe(eb >= 0 ? eb : elem0)
                                                               * (1. + (vit[f] * (i ? -1 : 1) >= 0 ? 1. : vit[f] ? -1.
                                                                        : 0.) *
                                                                  alpha_) * 0.5;
                }
              sum_dfac[i] = 0;
              for (n = 0; n < N; n++)
                for (m = 0; m < N; m++)
                  sum_dfac[i] += dfac(i, n, m);
            }
          for (i = 0; i < 2; i++)
            if ((e = f_e(f, i)) >= 0)
              {
                double inv_ve = 1.0 / ve(e);
                for (k = 0; k < nb_face_elem && (fb = e_f(e, k)) >= 0; k++)
                  if (fb < nb_faces && ch.fcl()(fb, 0) < 2) //partie "faces"
                    {
                      if ((fc = equiv(f, i, k)) >= 0 || elem1 < 0)
                        for (j = 0; j < 2; j++) //equivalence : face fd -> face fb
                          {
                            eb = f_e(f, j), fd = (j == i ? fb : fc); //element/face sources
                            mult = (fd < 0 || domaine.dot(&nf(fb, 0), &nf(fd, 0)) > 0 ? 1 : -1) *
                                   (fd >= 0 ? pf(fd) / pe(eb) : 1); //multiplicateur pour passer de vf a ve
                            for (n = 0; n < N; n++)
                              for (m = 0; m < N; m++)
                                if (dfac(j, n, m))
                                  {
                                    double fac = (e == elem0 ? 1 : -1) * vfd(fb, e != f_e(fb, 0)) *
                                                 dfac(j, n, m) * inv_ve;
                                    if (fd >= 0)
                                      secmem[fb] -= fac * mult * vit[fd]; //autre face calculee
                                    else
                                      {
                                        const Cond_lim_base& my_cl = cls[ch.fcl()(f, 1)].valeur();
                                        if (sub_type(Dirichlet, my_cl)) // sinon : paroi -> pas de contrib
                                          for (d = 0; d < D; d++)  //CL de Dirichlet
                                            secmem[fb] -= fac * nf(fb, d) / fs(fb) *
                                                          ref_cast(Dirichlet, my_cl).val_imp(
                                                            ch.fcl()(f, 2), N * d + m) * inv_masse;
                                      }
                                    if (comp) secmem[fb] += fac * vit[fb]; //partie v div(alpha rho v)
                                  }
                          }
                      else
                        {
                          for (j = 0; j < 2; j++)  //pas d'equivalence : n_f * operateur aux elements
                            {
                              for (eb = f_e(f, j), l = 0; l < nb_face_elem && (fc = e_f(eb, l)) >= 0; l++)
                                {
                                  double num =
                                    (e == f_e(fb, 0) ? 1 : -1) * (e == elem0 ? 1 : -1) * fs(fc) * fs(fb) *
                                    domaine.dot(&xv(fc, 0), &xv(fb, 0), &xp(eb, 0), &xp(e, 0)) *
                                    (eb == f_e(fc, 0) ? 1 : -1);
                                  double den = ve(eb) * ve(e);
                                  if (std::fabs(num) > 1e-9 * den)
                                    {
                                      double num_den = num / den;
                                      secmem[fb] -= sum_dfac[j] * num_den * vit[fc];
                                    }
                                }
                              if (comp)
                                for (l = 0; l < nb_face_elem && (fc = e_f(e, l)) >= 0; l++)
                                  {
                                    double num = (e == f_e(fb, 0) ? 1 : -1) * (e == elem0 ? 1 : -1) * fs(fc) *
                                                 fs(fb) *
                                                 domaine.dot(&xv(fc, 0), &xv(fb, 0), &xp(e, 0), &xp(e, 0)) *
                                                 (e == f_e(fc, 0) ? 1 : -1);
                                    double den = ve(e) * ve(e);
                                    if (std::fabs(num) > 1e-9 * den)
                                      {
                                        double num_den = num / den;
                                        secmem[fb] += sum_dfac[j] * num_den * vit[fc];
                                      }
                                  }
                            }
                        }
                    }
              }
        }
    }

  return secmem;
}

/*! @brief on assemble la matrice.
 *
 */
inline void Op_Conv_EF_Stab_PolyMAC_CDO_Face::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  if (has_interface_blocs())
    {
      Operateur_base::contribuer_a_avec(inco, matrice);
      return;
    }

  const Domaine_PolyMAC_CDO& domaine = le_dom_poly_.valeur();
  const Champ_Face_PolyMAC_CDO& ch = ref_cast(Champ_Face_PolyMAC_CDO, equation().inconnue());
  const IntTab& f_e = domaine.face_voisins(), &e_f = domaine.elem_faces(), &equiv = domaine.equiv();
  const DoubleTab& xp = domaine.xp(), &xv = domaine.xv(), &vfd = domaine.volumes_entrelaces_dir(), &vit = vitesse_->valeurs();
  const DoubleVect& fs = domaine.face_surfaces(), &ve = domaine.volumes(), &pf = porosite_f, &pe = porosite_e;
  const DoubleTab& nf = domaine.face_normales();

  int i, j, k, l, e, eb, f, fb, fc, fd, m, n, N = inco.line_size(), comp = !incompressible_;
  double mult;

  assert(N == 1);
  DoubleTrav dfac(2, N, N), masse(N, N);
  for (f = 0; f < domaine.nb_faces_tot(); f++)
    if (f_e(f, 0) >= 0 && (f_e(f, 1) >= 0 ||  ch.fcl()(f, 0) == 1 ||  ch.fcl()(f, 0) == 3))
      {
        for (i = 0, dfac = 0; i < 2; i++)
          {
            //masse : diagonale + masse ajoutee si correlation
            masse(0, 0) = std::fabs(vit[f]) > 1e-10 ? inco(f) / vit[f] : 1.0;
            //contribution a dfac
            for (eb = f_e(f, i), n = 0; n < N; n++)
              for (m = 0; m < N; m++)
                dfac( ch.fcl()(f, 0) == 1 ? 0 : i, n, m) += fs(f) * inco[f] * pe(eb >= 0 ? eb : f_e(f, 0))
                                                            * (1. + (vit[f] * (i ? -1 : 1) >= 0 ? 1. : vit[f] ? -1. : 0.) * alpha_) / 2;
          }
        for (i = 0; i < 2; i++)
          if ((e = f_e(f, i)) >= 0)
            {
              for (k = 0; k < e_f.dimension(1) && (fb = e_f(e, k)) >= 0; k++)
                if (fb < domaine.nb_faces() &&  ch.fcl()(fb, 0) < 2) //partie "faces"
                  {
                    if ((fc = equiv(f, i, k)) >= 0 || f_e(f, 1) < 0)
                      for (j = 0; j < 2; j++) //equivalence : face fd -> face fb
                        {
                          eb = f_e(f, j), fd = (j == i ? fb : fc); //element/face sources
                          mult = (fd < 0 || domaine.dot(&nf(fb, 0), &nf(fd, 0)) > 0 ? 1 : -1) * (fd >= 0 ? pf(fd) / pe(eb) : 1); //multiplicateur pour passer de vf a ve
                          for (n = 0; n < N; n++)
                            for (m = 0; m < N; m++)
                              if (dfac(j, n, m))
                                {
                                  double fac = (e == f_e(f, 0) ? 1 : -1) * vfd(fb, e != f_e(fb, 0)) * dfac(j, n, m) / ve(e);
                                  if (fd >= 0) matrice(fb,fd) += fac * mult; //autre face calculee
                                  if (comp) matrice(fb,fb) -= fac; //partie v div(alpha rho v)
                                }
                        }
                    else for (j = 0; j < 2; j++)  //pas d'equivalence : n_f * operateur aux elements
                        {
                          for (eb = f_e(f, j), l = 0; l < e_f.dimension(1) && (fc = e_f(eb, l)) >= 0; l++)
                            {
                              double num = (e == f_e(fb, 0) ? 1 : -1) * (e == f_e(f, 0) ? 1 : -1) * fs(fc) * fs(fb) * domaine.dot(&xv(fc, 0), &xv(fb, 0), &xp(eb, 0), &xp(e, 0)) * (eb == f_e(fc, 0) ? 1 : -1);
                              double den = ve(eb) * ve(e);
                              if (std::fabs(num) > 1e-9 * den)
                                {
                                  double num_den = num/den;
                                  for (n = 0; n < N; n++)
                                    for (m = 0; m < N; m++)
                                      if (dfac(j, n, m))
                                        {
                                          double fac = dfac(j, n, m) * num_den;
                                          matrice(fb,fc) += fac;
                                        }
                                }
                            }
                          if (comp)
                            for (l = 0; l < e_f.dimension(1) && (fc = e_f(e, l)) >= 0; l++)
                              {
                                double num = (e == f_e(fb, 0) ? 1 : -1) * (e == f_e(f, 0) ? 1 : -1) * fs(fc) * fs(fb) * domaine.dot(&xv(fc, 0), &xv(fb, 0), &xp(e, 0), &xp(e, 0)) * (e == f_e(fc, 0) ? 1 : -1);
                                double den = ve(e) * ve(e);
                                if (std::fabs(num) > 1e-9 * den)
                                  {
                                    double num_den = num/den;
                                    for (n = 0; n < N; n++)
                                      for (m = 0; m < N; m++)
                                        if (dfac(j, n, m))
                                          {
                                            double fac = dfac(j, n, m) * num_den;
                                            matrice(fb,fc) -= fac;
                                          }
                                  }
                              }
                        }
                  }
            }
      }
}

void Op_Conv_EF_Stab_PolyMAC_CDO_Face::set_incompressible(const int flag)
{
  if (flag == 0)
    {
      Cerr << "Compressible form of operator \"" << que_suis_je() << "\" :" << finl;
      Cerr << "Discretization of \u2207(inco \u2297 v) - v \u2207.(inco)" << finl;
    }
  incompressible_ = flag;
}
