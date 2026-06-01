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

#include <Dirichlet_entree_fluide_leaves.h>
#include <Dirichlet_paroi_defilante.h>
#include <Champ_Uniforme_Morceaux.h>
#include <Dirichlet_paroi_fixe.h>
#include <Modele_turbulence_hyd_base.h>
#include <Champ_Uniforme.h>
#include <Champ_Face_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <distances_VDF.h>
#include <Equation_base.h>
#include <Champ_Don_lu.h>
#include <Fluide_base.h>
#include <Periodique.h>
#include <Option_VDF.h>
#include <Navier.h>

// XXX : Elie Saikali : keeping Champ_Face as well, otherwise typer and reprise would need to be changed in many places ...
Implemente_instanciable(Champ_Face_VDF,"Champ_Face|Champ_Face_VDF",Champ_Face_base);

Sortie& Champ_Face_VDF::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Champ_Face_VDF::readOn(Entree& s)
{
  lire_donnees (s) ;
  return s ;
}

int Champ_Face_VDF::fixer_nb_valeurs_nodales(int nb_noeuds)
{
  assert(nb_noeuds == domaine_vdf().nb_faces());
  const MD_Vector& md = domaine_vdf().md_vector_faces();
  // Problem: nb_comp equals 2 but only one dimension is wanted !!!
  // HACK :
  int old_nb_compo = nb_compo_;
  if(nb_compo_ >= dimension) nb_compo_ /= dimension;
  creer_tableau_distribue(md);
  nb_compo_ = old_nb_compo;
  return nb_noeuds;
}

void Champ_Face_VDF::dimensionner_tenseur_Grad()
{
  tau_diag_.resize(domaine_vdf().nb_elem(), dimension);
  tau_croises_.resize(domaine_vdf().nb_aretes(), 2);
}

Champ_base& Champ_Face_VDF::affecter_(const Champ_base& ch)
{
  DoubleTab& val = valeurs();
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const IntVect& orientation = domaine_VDF.orientation();
  const int N = val.dimension(1), unif = sub_type(Champ_Uniforme, ch), D = dimension;

  if ((sub_type(Champ_Uniforme_Morceaux, ch)) || (sub_type(Champ_Don_lu, ch)))
    {
      const DoubleTab& v = ch.valeurs();
      int ndeb_int = domaine_VDF.premiere_face_int();
      const IntTab& f_e = domaine_VDF.face_voisins();

      for (int f = 0; f < ndeb_int; f++)
        {
          const int ori = orientation(f);
          const int e = f_e(f, f_e(f, 0) == -1);
          for (int n = 0; n < N; n++)
            val(f, n) = v(e, N * ori + n);
        }

      for (int f = ndeb_int; f < domaine_VDF.nb_faces(); f++)
        {
          const int ori = orientation(f);
          for (int n = 0; n < N; n++)
            val(f, n) = 0.5 * (v(f_e(f, 0), N * ori + n) + v(f_e(f, 1), N * ori + n));
        }
    }
  else
    {

      DoubleTab eval;
      if (unif) eval = ch.valeurs();
      else eval.resize(val.dimension(0), N * D), ch.valeur_aux(domaine_VDF.xv(), eval);

      for (int f = 0; f < domaine_VDF.nb_faces(); f++)
        for (int n = 0; n < N; n++)
          val(f, n) = eval(unif ? 0 : f, N * orientation(f) + n);
    }
  val.echange_espace_virtuel();
  return *this;
}

const Champ_Proto& Champ_Face_VDF::affecter(const double x1, const double x2)
{
  const IntVect& orientation = domaine_vdf().orientation();
  DoubleTab& val = valeurs();
  for (int num_face = 0; num_face < val.size(); num_face++)
    {
      int ori = orientation(num_face);
      switch(ori)
        {
        case 0:
          val(num_face) = x1;
          break;
        case 1:
          val(num_face) = x2;
          break;
        }
    }
  return *this;
}

const Champ_Proto& Champ_Face_VDF::affecter(const double x1, const double x2, const double x3)
{
  const IntVect& orientation = domaine_vdf().orientation();
  DoubleTab& val = valeurs();
  for (int num_face = 0; num_face < val.size(); num_face++)
    {
      int ori = orientation(num_face);
      switch(ori)
        {
        case 0:
          val(num_face) = x1;
          break;
        case 1:
          val(num_face) = x2;
          break;
        case 2:
          val(num_face) = x3;
          break;
        }
    }
  return *this;
}

const Champ_Proto& Champ_Face_VDF::affecter(const DoubleTab& v)
{
  const IntVect& orientation = domaine_vdf().orientation();
  DoubleTab& val = valeurs();

  if (v.nb_dim() == 2)
    {
      if (v.dimension(1) == dimension)
        {
          if (v.dimension(0) == val.size())
            for (int num_face = 0; num_face < val.size(); num_face++)
              val(num_face) = v(num_face, orientation(num_face));
          else
            {
              Cerr << "Erreur TRUST dans Champ_Face_VDF::affecter(const DoubleTab& )" << finl;
              Cerr << "Les dimensions du DoubleTab passe en parametre sont incompatibles avec celles du Champ_Face_VDF " << finl;
              Process::exit();
            }
        }
      else
        {
          Cerr << "Erreur TRUST dans Champ_Face_VDF::affecter(const DoubleTab& )" << finl;
          Cerr << "Les dimensions du DoubleTab passe en parametre sont incompatibles avec celles du Champ_Face_VDF " << finl;
          Process::exit();
        }
    }
  return *this;
}

// Cas CL periodique : assure que les valeurs sur des faces periodiques en vis a vis sont identiques. Pour cela on prend la demi somme des deux valeurs.
void Champ_Face_VDF::verifie_valeurs_cl()
{
  const Domaine_Cl_dis_base& zcl = domaine_Cl_dis();
  int nb_cl = zcl.nb_cond_lim();
  DoubleTab& ch_tab = valeurs();
  int ndeb, nfin, num_face;

  for (int i = 0; i < nb_cl; i++)
    {
      const Cond_lim_base& la_cl = zcl.les_conditions_limites(i).valeur();
      if (sub_type(Periodique, la_cl))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl);
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          int voisine;
          double moy;

          for (num_face = ndeb; num_face < nfin; num_face++)
            {
              voisine = la_cl_perio.face_associee(num_face - ndeb) + ndeb;
              if (ch_tab[num_face] != ch_tab[voisine])
                {
                  moy = 0.5 * (ch_tab[num_face] + ch_tab[voisine]);
                  ch_tab[num_face] = moy;
                  ch_tab[voisine] = moy;
                }
            }
        }
    }
  ch_tab.echange_espace_virtuel();
}

/*! @brief Returns the value the field should have on a boundary face, according to the boundary conditions.
 *
 * @brief The index is counted in the list of real boundary faces. The time considered is the
 *  current time of the Champ_Face.
 *  The implementation has changed: these values are no longer stored in the field.
 *
 */
double Champ_Face_VDF::val_imp_face_bord_private(int face, int comp) const
{
  const Domaine_Cl_VDF& zclo = ref_cast(Domaine_Cl_VDF, equation().domaine_Cl_dis());
  return Champ_Face_get_val_imp_face_bord_sym(valeurs(), temps(), face, comp, zclo);
}

// WEC : never called !!
double Champ_Face_VDF::val_imp_face_bord_private(int face, int comp1, int comp2) const
{
  Cerr << "Champ_Face_VDF::val_imp_face_bord(,,) exit" << finl;
  Process::exit();
  return 0; // For compilers
}

// This function returns:
//   1 if the fluid is outgoing on face num_face
//   0 if the face corresponds to a fluid re-entry
int Champ_Face_VDF::compo_normale_sortante(int num_face) const
{
  int signe = 1;
  double vit_norm;
  // signe is -1 if face_voisins(num_face,0) is on the outside
  // signe is  1 if face_voisins(num_face,1) is on the outside
  if (domaine_vdf().face_voisins(num_face, 0) == -1)
    signe = -1;
  vit_norm = valeurs()(num_face) * signe;
  return (vit_norm > 0);
}

DoubleTab& Champ_Face_VDF::trace(const Frontiere_dis_base& fr, DoubleTab& x, double tps, int distant) const
{
  return Champ_Face_VDF_implementation::trace(fr, valeurs(tps), x, distant);
}

void Champ_Face_VDF::mettre_a_jour(double un_temps)
{
  Champ_Inc_base::mettre_a_jour(un_temps);
}

void Champ_Face_VDF::calculer_rotationnel_ordre2_centre_element(DoubleTab& rot) const
{
  const DoubleTab& val = valeurs();
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  int nb_elem = domaine_VDF.nb_elem();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntTab& elem_faces = domaine_VDF.elem_faces();

  if (dimension == 2)
    calrotord2centelemdim2(rot, val, domaine_VDF, nb_elem, face_voisins, elem_faces);
  else if (dimension == 3)
    calrotord2centelemdim3(rot, val, domaine_VDF, nb_elem, face_voisins, elem_faces);
}

int Champ_Face_VDF::imprime(Sortie& os, int ncomp) const
{
  imprime_Face(os, ncomp);
  return 1;
}

void Champ_Face_VDF::calcul_critere_Q(DoubleTab& Q, const Domaine_Cl_VDF& domaine_Cl_VDF)
{
  // Q=0.5*(\Omega_{ij}*\Omega_{ij}-S_{ij}*S_{ij})=-0.25*du_i/dx_j*du_j/dx_i

  const Domaine_VDF& domaine_VDF = domaine_vdf();
  Champ_Face_VDF& vit = *this;
  const DoubleTab& vitesse = valeurs();
  const int nb_elem = domaine_VDF.nb_elem();
  const int nb_elem_tot = domaine_VDF.nb_elem_tot();
  int num_elem, i, j, N = vitesse.line_size();
  double crit, deriv1, deriv2;

  if (N!=1) Process::exit(que_suis_je() + "::calcul_critere_Q : the velocity field must be single phase !!");

  DoubleTab gradient_elem(nb_elem_tot, dimension, dimension, N);
  gradient_elem = 0.;

  vit.calcul_duidxj(vitesse, gradient_elem, domaine_Cl_VDF);

  for (num_elem = 0; num_elem < nb_elem; num_elem++)
    {
      crit = 0.;
      for (i = 0; i < dimension; i++)
        for (j = 0; j < dimension; j++)
          {
            deriv1 = gradient_elem(num_elem, i, j, 0);
            deriv2 = gradient_elem(num_elem, j, i, 0);

            crit += -0.25 * deriv1 * deriv2;
          }
      Q[num_elem] = crit;
    }
}

void Champ_Face_VDF::calcul_y_plus(DoubleTab& y_plus, const Domaine_Cl_VDF& domaine_Cl_VDF)
{
  // Initialize the y_plus field with a negative value,
  // so that when visualizing the field near the wall,
  // one only needs to remove negative values and only
  // wall values remain.

  int ndeb, nfin, elem, ori, l_unif;
  double norm_tau, u_etoile, norm_v = 0, dist, val0, val1, val2, d_visco = 0, visco = 1.;
  y_plus = -1.;

  const Champ_Face_VDF& vit = *this;
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntVect& orientation = domaine_VDF.orientation();
  const Equation_base& eqn_hydr = equation();
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don_base& ch_visco_cin = le_fluide.viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin.valeurs();
  //DoubleTab& tab_visco = ch_visco_cin.valeurs();

  if (sub_type(Champ_Uniforme, ch_visco_cin))
    {
      visco = tab_visco(0, 0);
      l_unif = 1;
    }
  else
    l_unif = 0;

  // Only change values < DMINFLOAT (the old code is not parallel)
  /* GF we should not change tab_visco here!
   if (!l_unif)
   {
   const int n = tab_visco.size_array();
   ArrOfDouble& v = tab_visco;
   for (int i = 0; i < n; i++)
   if (v[i] < DMINFLOAT)
   v[i] = DMINFLOAT;
   }
   */

  DoubleTab yplus_faces(1, 1); // will contain yplus values if available
  int yplus_already_computed = 0; // flag

  const RefObjU& modele_turbulence = eqn_hydr.get_modele(TURBULENCE);
  if (modele_turbulence && sub_type(Modele_turbulence_hyd_base, modele_turbulence.valeur()))
    {
      const Modele_turbulence_hyd_base& mod_turb = ref_cast(Modele_turbulence_hyd_base, modele_turbulence.valeur());
      const Turbulence_paroi_base& loipar = mod_turb.loi_paroi();
      if (loipar.use_shear())
        {
          yplus_faces.resize(domaine_vdf().nb_faces_tot());
          yplus_faces.ref(loipar.tab_d_plus());
          yplus_already_computed = 1;
        }
    }

  for (int n_bord = 0; n_bord < domaine_VDF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VDF.les_conditions_limites(n_bord);

      if (sub_type(Dirichlet_paroi_fixe, la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          for (int num_face = ndeb; num_face < nfin; num_face++)
            {

              if (face_voisins(num_face, 0) != -1)
                elem = face_voisins(num_face, 0);
              else
                elem = face_voisins(num_face, 1);

              if (yplus_already_computed)
                {
                  // y+ is only defined on faces so we take the face value to put in the element
                  y_plus(elem) = yplus_faces(num_face);
                }
              else
                {
                  if (dimension == 2)
                    {
                      ori = orientation(num_face);
                      norm_v = norm_2D_vit(vit.valeurs(), elem, ori, domaine_VDF, val0);
                    }
                  else if (dimension == 3)
                    {
                      ori = orientation(num_face);
                      norm_v = norm_3D_vit(vit.valeurs(), elem, ori, domaine_VDF, val1, val2);
                    } // dim 3

                  if (axi)
                    dist = domaine_VDF.dist_norm_bord_axi(num_face);
                  else
                    dist = domaine_VDF.dist_norm_bord(num_face);
                  if (l_unif)
                    d_visco = visco;
                  else
                    d_visco = tab_visco[elem];

                  // PQ : 01/10/03 : corrections relative to the original version
                  norm_tau = d_visco * norm_v / dist;

                  u_etoile = sqrt(norm_tau);
                  y_plus(elem) = dist * u_etoile / d_visco;

                } // else yplus already computed
            } // loop on faces
        } // End fixed wall
    } // End loop over boundaries
}

/*! @brief Returns gij at elements from the face velocity (gij represents the partial derivative dui/dxj).
 *
 *  @brief From gij, one can compute Sij = 0.5(gij(i,j)+gij(j,i)).
 *
 * @param vitesse Velocity field at faces.
 * @param gij Output velocity gradient tensor at elements.
 * @param domaine_Cl_VDF Boundary condition domain.
 * @return Reference to gij.
 */
DoubleTab& Champ_Face_VDF::calcul_duidxj(const DoubleTab& vitesse, DoubleTab& gij, const Domaine_Cl_VDF& domaine_Cl_VDF) const
{
  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, mon_equation->inconnue());
  const Domaine_Cl_VDF& dclvdf = ref_cast(Domaine_Cl_VDF, vit.domaine_Cl_dis());
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const int nb_elem = domaine_VDF.domaine().nb_elem_tot(), N = vitesse.line_size();
  const IntTab& face_voisins = domaine_VDF.face_voisins(), &elem_faces = domaine_VDF.elem_faces(), &Qdm = domaine_VDF.Qdm();
  const IntVect& orientation = domaine_VDF.orientation();

  const int prem_am = domaine_VDF.premiere_arete_mixte(), dern_am = prem_am + domaine_VDF.nb_aretes_mixtes();
  const int prem_ai = domaine_VDF.premiere_arete_interne(), dern_ai = prem_ai + domaine_VDF.nb_aretes_internes();
  IntVect element(4);
  gij = 0.;

  // Loop over all edges that allow computing the cross terms du_i/dx_j
  // (the non-cross terms are computed by looping over elements)


  // Start with the boundary edges
  int ndeb = domaine_VDF.premiere_arete_bord(), nfin = ndeb + domaine_VDF.nb_aretes_bord();
  for (int num_arete = ndeb; num_arete < nfin; num_arete++)
    for (int n=0; n<N; n++)
      {
        const int n_type = domaine_Cl_VDF.type_arete_bord(num_arete - ndeb);

        if (n_type == 4) // periodic edge type
          {
            const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
            const int i = orientation(num0), j = orientation(num2);

            const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
            const double temp2 = (vitesse(num3, n) - vitesse(num2, n)) / domaine_VDF.dist_face_period(num2, num3, i); // du_j / dx_i

            element(0) = face_voisins(num0, 0);
            element(1) = face_voisins(num0, 1);
            element(2) = face_voisins(num1, 0);
            element(3) = face_voisins(num1, 1);

            for (int k = 0; k < 4; k++)
              {
                // 1) 0.5 : for periodicity, because we will distribute twice over elements touching this edge since it exists twice.
                // 2) 0.25 : distribute the velocity gradient over the 4 surrounding elements.
                gij(element(k), i, j, n) += temp1 * 0.5 * 0.25;
                gij(element(k), j, i, n) += temp2 * 0.5 * 0.25;
              }
          }
        else if (n_type == 3 && Option_VDF::traitement_gradients) /* NAVIER - NAVIER */
          {
            const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), signe = Qdm(num_arete, 3);
            const int i = orientation(num0), j = orientation(num2);

            const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
            const double coeff_frot = (Champ_Face_coeff_frottement_grad_face_bord(num0, n, dclvdf)+Champ_Face_coeff_frottement_grad_face_bord(num1, n, dclvdf))/2.;
            const double temp2 = -signe * coeff_frot * vitesse(num2, n);

            element(0) = face_voisins(num2, 0);
            element(1) = face_voisins(num2, 1);

            for (int k = 0; k < 2; k++)
              {
                // 1) 0.25 : distribute the velocity gradient over the 4 surrounding elements.
                gij(element(k), i, j, n) += temp1 * 0.25;
                gij(element(k), j, i, n) += temp2 * 0.25;
              }
          }
        else if (Option_VDF::traitement_gradients && (n_type == 5 || n_type == 6))
          Process::exit("Issue in Champ_Face_VDF::calcul_duidxj ... This case is not yet considered. Contact the TRUST team.");
        else /* other boundary edges ... */
          {
            const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), signe = Qdm(num_arete, 3);
            const int i = orientation(num0), j = orientation(num2);

            const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
            const double vit_imp = 0.5 * (vit.val_imp_face_bord_private(num0, N*j+n) + vit.val_imp_face_bord_private(num1, N*j+n)); // tangential velocity

            //In this part, we keep the coding of Hyd_SGE_Wale_VDF (num1 not num2) to compute the distance between the cell center and the boundary.
            const double temp2 = -signe * (vitesse(num2, n) - vit_imp) / domaine_VDF.dist_norm_bord(num1);

            element(0) = face_voisins(num2, 0);
            element(1) = face_voisins(num2, 1);

            for (int k = 0; k < 2; k++)
              {
                // 1) 0.25 : distribute the velocity gradient over the 4 surrounding elements.
                gij(element(k), i, j, n) += temp1 * 0.25;
                gij(element(k), j, i, n) += temp2 * 0.25;
              }
          }
      }

  // Continue with corners
  ndeb = domaine_VDF.premiere_arete_coin(), nfin = ndeb + domaine_VDF.nb_aretes_coin();

  for (int num_arete = ndeb; num_arete < nfin; num_arete++)
    for (int n=0; n<N; n++)
      {
        const int n_type = domaine_Cl_VDF.type_arete_coin(num_arete - ndeb);

        if (n_type == 0) // perio-perio edge type
          {
            const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
            const int i = orientation(num0), j = orientation(num2);

            const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
            const double temp2 = (vitesse(num3, n) - vitesse(num2, n)) / domaine_VDF.dist_face_period(num2, num3, i); // du_j / dx_i

            element(0) = face_voisins(num0, 0);
            element(1) = face_voisins(num0, 1);
            element(2) = face_voisins(num1, 0);
            element(3) = face_voisins(num1, 1);

            for (int k = 0; k < 4; k++)
              {
                // 1) 0.5 : for periodicity, because we will distribute twice over elements touching this edge since it exists twice.
                // 2) 0.5 : same as above, because this is a perio-perio corner.
                // 3) 0.25 : distribute the velocity gradient over the 4 surrounding elements.
                gij(element(k), i, j, n) += temp1 * 0.5 * 0.5 * 0.25;
                gij(element(k), j, i, n) += temp2 * 0.5 * 0.5 * 0.25;
              }
          }

        if (n_type == 1) // perio-wall edge type
          {
            const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), signe = Qdm(num_arete, 3);
            const int i = orientation(num1), j = orientation(num2);

            const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
            const double vit_imp = 0.5 * (vit.val_imp_face_bord_private(num0, N*j+n) + vit.val_imp_face_bord_private(num1, N*j+n)); // tangential velocity

            const double temp2 = -signe * (vitesse(num2, n) - vit_imp) / domaine_VDF.dist_norm_bord(num1);

            element(0) = face_voisins(num2, 0);
            element(1) = face_voisins(num2, 1);

            for (int k = 0; k < 2; k++)
              {
                // 1) 0.5 : for periodicity, because we will distribute twice over elements touching this edge since it exists twice.
                // 2) 0.25 : distribute the velocity gradient over the 4 surrounding elements.
                gij(element(k), i, j, n) += temp1 * 0.5 * 0.25;
                gij(element(k), j, i, n) += temp2 * 0.5 * 0.25;
              }
          }

        // XXX : Elie Saikali : adding this for corners only if option_vdf is active for now ...

        if (Option_VDF::traitement_gradients && Option_VDF::traitement_coins)
          {
            if (n_type == 14 || n_type == 15) // fluid-wall or wall-fluid edge type
              {
                const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), signe = Qdm(num_arete, 3);
                const int i = orientation(num1), j = orientation(num2);

                const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
                const double vit_imp = 0.5 * (vit.val_imp_face_bord_private(num0, N*j+n) + vit.val_imp_face_bord_private(num1, N*j+n)); // tangential velocity

                const double temp2 = -signe * (vitesse(num2, n) - vit_imp) / domaine_VDF.dist_norm_bord(num1);

                element(0) = face_voisins(num2, 0);
                element(1) = face_voisins(num2, 1);

                for (int k = 0; k < 2; k++)
                  if (element(k) != -1)
                    {
                      gij(element(k), i, j, n) += temp1 * 0.25;
                      gij(element(k), j, i, n) += temp2 * 0.25;
                    }
              }
            else if (n_type == 3 || n_type == 4 || n_type == 8) // fluid-navier edge type
              {
                const int num0 = Qdm(num_arete, 0),  num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
                const int f1 = num0 > -1 ? num0 : num1, f2 = num2 > -1 ? num2 : num3;
                const int i = orientation(f1), j = orientation(f2);

                const double coeff_frot1 = Champ_Face_coeff_frottement_grad_face_bord(f1, n, dclvdf), coeff_frot2 = Champ_Face_coeff_frottement_grad_face_bord(f2, n, dclvdf);

//                int signe = f2 == num3 ? -1 : 1;
//                const double temp1 = coeff_frot2 * signe * vitesse(f1, n);
//                const double temp2 = coeff_frot1 * signe * vitesse(f2, n);
                const double temp1 = coeff_frot2 * (face_voisins(f2, 0)==-1 ? 1:-1)* vitesse(f1, n);
                const double temp2 = coeff_frot1 * (face_voisins(f1, 0)==-1 ? 1:-1)* vitesse(f2, n);


                element(0) = face_voisins(f1, 0);
                element(1) = face_voisins(f1, 1);

                for (int k = 0; k < 2; k++)
                  if (element(k) != -1)
                    {
                      gij(element(k), i, j, n) += temp1 * 0.25;
                      gij(element(k), j, i, n) += temp2 * 0.25;
                    }
              }
          }
      }

  // Continue with mixed edges

  for (int num_arete = prem_am; num_arete < dern_am; num_arete++)
    for (int n=0; n<N; n++)
      {
        const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
        const int i = orientation(num0), j = orientation(num2);

        const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
        const double temp2 = (vitesse(num3, n) - vitesse(num2, n)) / domaine_VDF.dist_face_period(num2, num3, i); // du_j / dx_i

        element(0) = face_voisins(num0, 0);
        element(1) = face_voisins(num0, 1);
        element(2) = face_voisins(num1, 0);
        element(3) = face_voisins(num1, 1);

        for (int k = 0; k < 4; k++)
          if (element(k) != -1)
            {
              // 1) 0.25 : distribute the velocity gradient over the 3 surrounding elements.
              // We check element(k)!=-1 because in that case it is the "missing cell"!
              gij(element(k), i, j, n) += temp1 * 0.25;
              gij(element(k), j, i, n) += temp2 * 0.25;
            }
      }

  // Continue with internal edges

  for (int num_arete = prem_ai; num_arete < dern_ai; num_arete++)
    for (int n=0; n<N; n++)
      {
        const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
        const int i = orientation(num0), j = orientation(num2);

        const double temp1 = (vitesse(num1, n) - vitesse(num0, n)) / domaine_VDF.dist_face_period(num0, num1, j); // du_i / dx_j
        assert(est_egal(domaine_VDF.dist_face_period(num0, num1, j), domaine_VDF.dist_face(num0, num1, j)));

        const double temp2 = (vitesse(num3, n) - vitesse(num2, n)) / domaine_VDF.dist_face_period(num2, num3, i); // du_j / dx_i
        assert(est_egal(domaine_VDF.dist_face_period(num2, num3, j), domaine_VDF.dist_face(num2, num3, j)));

        element(0) = face_voisins(num0, 0);
        element(1) = face_voisins(num0, 1);
        element(2) = face_voisins(num1, 0);
        element(3) = face_voisins(num1, 1);

        for (int k = 0; k < 4; k++)
          {
            // 1) 0.25 : distribute the velocity gradient over the 4 surrounding elements.
            gij(element(k), i, j, n) += temp1 * 0.25;
            gij(element(k), j, i, n) += temp2 * 0.25;
          }
      }


  // XXX : Elie Saikali : HACK for fluid-fluid corners
  // for this case (admittedly rare), beware of issues with velocity values at corners ... e.g. a champ_fonc_xyz x+y+z does not give the correct value at corners

  // Continue with corners

  ndeb = domaine_VDF.premiere_arete_coin(), nfin = ndeb + domaine_VDF.nb_aretes_coin();

  for (int num_arete = ndeb; num_arete < nfin; num_arete++)
    for (int n=0; n<N; n++)
      {
        const int n_type = domaine_Cl_VDF.type_arete_coin(num_arete - ndeb);

        if (Option_VDF::traitement_gradients && Option_VDF::traitement_coins)
          if (n_type == 16 ) // fluid-fluid edge type
            {
              const int num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2);
              const int i = orientation(num1), j = orientation(num2);

              element(0) = face_voisins(num2, 0);
              element(1) = face_voisins(num2, 1);

              for (int k = 0; k < 2; k++)
                if (element(k) != -1)
                  {
                    // XXX : 1/3 because we want something like: (a+b+c+d)/4 = (a+b+c)/3 => d = (a+b+c)/3
                    gij(element(k), i, j, n) += gij(element(k), i, j, n) / 3.;
                    gij(element(k), j, i, n) += gij(element(k), j, i, n) / 3.;
                  }
            }
      }

  // 2nd part: loop over elements and fill Sij for the non-cross derivatives (du_i / dx_i).
  // In these cases, the gradient is computed directly within the element and not redistributed.

  for (int elem = 0; elem < nb_elem; elem++)
    for (int n=0; n<N; n++)
      for (int i = 0; i < dimension; i++)
        {
          double temp1 = (vitesse(elem_faces(elem, i), n) - vitesse(elem_faces(elem, i + dimension), n)) / domaine_VDF.dim_elem(elem, orientation(elem_faces(elem, i)));
          gij(elem, i, i, n) = -temp1;
        }

  return gij;
}

/*! @brief Returns gij at elements from the element velocity (gij represents the partial derivative dui/dxj).
 *
 *  @brief From gij, one can compute Sij = 0.5(gij(i,j)+gij(j,i)).
 *
 * @param in_vel Input velocity field at elements.
 * @param gij Output velocity gradient tensor at elements.
 * @return Reference to gij.
 */
DoubleTab& Champ_Face_VDF::calcul_duidxj(const DoubleTab& in_vel, DoubleTab& gij) const
{

  const Domaine_VDF& domaine_VDF = domaine_vdf();
  int nb_elem_tot = domaine_VDF.domaine().nb_elem_tot(), N = in_vel.line_size()/dimension;
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntTab& elem_faces = domaine_VDF.elem_faces();

  int num0, num1, num2, num3, num4, num5;
  int f0, f1, f2, f3, f4, f5;

  //
  // Calculate the Sij tensor
  //
  if (dimension == 2)
    {
      for (int element_number = 0; element_number < nb_elem_tot; element_number++)
        for (int n=0; n<N; n++)
          {
            f0 = elem_faces(element_number, 0);
            num0 = face_voisins(f0, 0);
            if (num0 == -1)
              num0 = element_number;
            f1 = elem_faces(element_number, 1);
            num1 = face_voisins(f1, 0);
            if (num1 == -1)
              num1 = element_number;
            f2 = elem_faces(element_number, 2);
            num2 = face_voisins(f2, 1);
            if (num2 == -1)
              num2 = element_number;
            f3 = elem_faces(element_number, 3);
            num3 = face_voisins(f3, 1);
            if (num3 == -1)
              num3 = element_number;

            gij(element_number, 0, 0, n) = 0.5 * ((in_vel(num2, N*0+n) - in_vel(num0, N*0+n)) / domaine_VDF.dim_elem(element_number, 0));
            gij(element_number, 0, 1, n) = 0.5 * ((in_vel(num3, N*0+n) - in_vel(num1, N*0+n)) / domaine_VDF.dim_elem(element_number, 1));
            gij(element_number, 1, 0, n) = 0.5 * ((in_vel(num2, N*1+n) - in_vel(num0, N*1+n)) / domaine_VDF.dim_elem(element_number, 0));
            gij(element_number, 1, 1, n) = 0.5 * ((in_vel(num3, N*1+n) - in_vel(num1, N*1+n)) / domaine_VDF.dim_elem(element_number, 1));
          }
    }
  else
    {
      for (int element_number = 0; element_number < nb_elem_tot; element_number++)
        for (int n=0; n<N; n++)
          {
            f0 = elem_faces(element_number, 0);
            num0 = face_voisins(f0, 0);
            if (num0 == -1)
              num0 = element_number;
            f1 = elem_faces(element_number, 1);
            num1 = face_voisins(f1, 0);
            if (num1 == -1)
              num1 = element_number;
            f2 = elem_faces(element_number, 2);
            num2 = face_voisins(f2, 0);
            if (num2 == -1)
              num2 = element_number;
            f3 = elem_faces(element_number, 3);
            num3 = face_voisins(f3, 1);
            if (num3 == -1)
              num3 = element_number;
            f4 = elem_faces(element_number, 4);
            num4 = face_voisins(f4, 1);
            if (num4 == -1)
              num4 = element_number;
            f5 = elem_faces(element_number, 5);
            num5 = face_voisins(f5, 1);
            if (num5 == -1)
              num5 = element_number;

            gij(element_number, 0, 0, n) = 0.5 * ((in_vel(num3, N*0+n) - in_vel(num0, N*0+n)) / domaine_VDF.dim_elem(element_number, 0));

            gij(element_number, 0, 1, n) = 0.5 * ((in_vel(num4, N*0+n) - in_vel(num1, N*0+n)) / domaine_VDF.dim_elem(element_number, 1));
            gij(element_number, 1, 0, n) = 0.5 * ((in_vel(num3, N*1+n) - in_vel(num0, N*1+n)) / domaine_VDF.dim_elem(element_number, 0));

            gij(element_number, 0, 2, n) = 0.5 * ((in_vel(num5, N*0+n) - in_vel(num2, N*0+n)) / domaine_VDF.dim_elem(element_number, 2));

            gij(element_number, 2, 0, n) = 0.5 * ((in_vel(num3, N*2+n) - in_vel(num0, N*2+n)) / domaine_VDF.dim_elem(element_number, 0));

            gij(element_number, 1, 1, n) = 0.5 * ((in_vel(num4, N*1+n) - in_vel(num1, N*1+n)) / domaine_VDF.dim_elem(element_number, 1));

            gij(element_number, 1, 2, n) = 0.5 * ((in_vel(num5, N*1+n) - in_vel(num2, N*1+n)) / domaine_VDF.dim_elem(element_number, 2));
            gij(element_number, 2, 1, n) = 0.5 * ((in_vel(num4, N*2+n) - in_vel(num1, N*2+n)) / domaine_VDF.dim_elem(element_number, 1));

            gij(element_number, 2, 2, n) = 0.5 * ((in_vel(num5, N*2+n) - in_vel(num2, N*2+n)) / domaine_VDF.dim_elem(element_number, 2));

          }
    }

  return gij;

}

/*! @brief Returns SMA_barre at elements from the face velocity.
 *
 *  @brief SMA_barre = Sij*Sij (summation over indices i and j).
 *
 *  S_barre(num_elem) is computed directly.
 *  The contribution_paroi parameter (here fixed to 0) allows excluding the wall contribution from the product SMA_barre = Sij*Sij.
 *
 */
DoubleVect& Champ_Face_VDF::calcul_S_barre_sans_contrib_paroi(const DoubleTab& vitesse, DoubleVect& SMA_barre, const Domaine_Cl_VDF& domaine_Cl_VDF) const
{
  const int contribution_paroi = 0;

  const Champ_Face_VDF& vit = ref_cast(Champ_Face_VDF, mon_equation->inconnue());
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const IntTab& face_voisins = domaine_VDF.face_voisins(), &elem_faces = domaine_VDF.elem_faces(), &Qdm = domaine_VDF.Qdm();
  const IntVect& orientation = domaine_VDF.orientation();

  const int nb_elem = domaine_VDF.domaine().nb_elem();
  const int prem_am = domaine_VDF.premiere_arete_mixte(), dern_am = prem_am + domaine_VDF.nb_aretes_mixtes();
  const int prem_ai = domaine_VDF.premiere_arete_interne(), dern_ai = prem_ai + domaine_VDF.nb_aretes_internes();

  ArrOfInt element(4);

  int ndeb = domaine_VDF.premiere_arete_bord(), nfin = ndeb + domaine_VDF.nb_aretes_bord();

  for (int num_arete = ndeb; num_arete < nfin; num_arete++)
    {
      int n_type = domaine_Cl_VDF.type_arete_bord(num_arete - ndeb);

      if (n_type == 4) // periodic edge type
        {
          const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
          const int i = orientation(num0), j = orientation(num2);

          const double temp1 = (vitesse[num1] - vitesse[num0]) / domaine_VDF.dist_face_period(num0, num1, j); // dv/dx
          const double temp2 = (vitesse[num3] - vitesse[num2]) / domaine_VDF.dist_face_period(num2, num3, i); // du/dy

          element[0] = face_voisins(num0, 0);
          element[1] = face_voisins(num0, 1);
          element[2] = face_voisins(num1, 0);
          element[3] = face_voisins(num1, 1);

          // compute the sum of cross terms: 2*( (0.5*Sij)^2+(0.5*Sji)^2)
          // Since we are on edges distributing over elements, multiply by 0.25: 0.25*(2*(2*0.5^2))=0.25*4*0.25=0.25
          // The 0.5 factor comes from traversing periodic faces: periodic edges are "the same", so we add twice what is needed -> 0.5
          for (int k = 0; k < 4; k++)
            SMA_barre[element[k]] += 0.5 * (temp1 + temp2) * (temp1 + temp2) * 0.25;
        }
      else
        {
          const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), signe = Qdm(num_arete, 3);
          const int j = orientation(num2);

          const double temp1 = (vitesse[num1] - vitesse[num0]) / domaine_VDF.dist_face_period(num0, num1, j);        // dv/dx
          double vit_imp = 0.5 * (vit.val_imp_face_bord_private(num0, j) + vit.val_imp_face_bord_private(num1, j));                // tangential velocity

          double temp2;

          if (n_type == 0 && contribution_paroi == 0)
            temp2 = 0;
          else
            temp2 = -signe * (vitesse[num2] - vit_imp) / domaine_VDF.dist_norm_bord(num1);

          element[0] = face_voisins(num2, 0);
          element[1] = face_voisins(num2, 1);

          // compute the sum of cross terms: 2*( (0.5*Sij)^2+(0.5*Sji)^2)
          // Since we are on edges distributing over elements, multiply by 0.25: 0.25*(2*(2*0.5^2))=0.25*4*0.25=0.25
          // Account for both symmetric terms: SijSij+SjiSji
          for (int k = 0; k < 2; k++)
            SMA_barre[element[k]] += (temp1 + temp2) * (temp1 + temp2) * 0.25;
        }
    }

  ndeb = domaine_VDF.premiere_arete_coin(), nfin = ndeb + domaine_VDF.nb_aretes_coin();

  for (int num_arete = ndeb; num_arete < nfin; num_arete++)
    {
      int n_type = domaine_Cl_VDF.type_arete_coin(num_arete - ndeb);

      if (n_type == 0) // perio-perio edge type
        {
          const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
          const int i = orientation(num0), j = orientation(num2);

          const double temp1 = (vitesse[num1] - vitesse[num0]) / domaine_VDF.dist_face_period(num0, num1, j);            // dv/dx
          const double temp2 = (vitesse[num3] - vitesse[num2]) / domaine_VDF.dist_face_period(num2, num3, i);            // du/dy

          element[0] = face_voisins(num0, 0);
          element[1] = face_voisins(num0, 1);
          element[2] = face_voisins(num1, 0);
          element[3] = face_voisins(num1, 1);

          // compute the sum of cross terms: 2*( (0.5*Sij)^2+(0.5*Sji)^2)
          // Since we are on edges distributing over elements, multiply by 0.25: 0.25*(2*(2*0.5^2))=0.25*4*0.25=0.25
          // The 0.5 factor comes from traversing periodic faces: periodic edges are "the same", so we add twice what is needed -> 0.5
          // another *0.5 because these are perio-perio edges traversed 4 times
          for (int k = 0; k < 4; k++)
            SMA_barre[element[k]] += 0.5 * 0.5 * (temp1 + temp2) * (temp1 + temp2) * 0.25;
        }

      if (n_type == 1) // perio-wall edge type
        {
          const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), signe = Qdm(num_arete, 3);
          const int j = orientation(num2);

          const double temp1 = (vitesse[num1] - vitesse[num0]) / domaine_VDF.dist_face_period(num0, num1, j);        // dv/dx
          const double vit_imp = 0.5 * (vit.val_imp_face_bord_private(num0, j) + vit.val_imp_face_bord_private(num1, j));             // tangential velocity

          double temp2;

          if (contribution_paroi == 0)
            temp2 = 0;
          else
            temp2 = -signe * (vitesse[num2] - vit_imp) / domaine_VDF.dist_norm_bord(num1);

          element[0] = face_voisins(num2, 0);
          element[1] = face_voisins(num2, 1);

          for (int k = 0; k < 2; k++)
            SMA_barre[element[k]] += 0.5 * (temp1 + temp2) * (temp1 + temp2) * 0.25;
        }

      if (Option_VDF::traitement_gradients && Option_VDF::traitement_coins)
        if (n_type == 14 || n_type == 15 || n_type == 16) // fluid-wall, wall-fluid or fluid-fluid edge type
          {
            const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), signe = Qdm(num_arete, 3);
            const int j = orientation(num2);

            const double temp1 = (vitesse[num1] - vitesse[num0]) / domaine_VDF.dist_face_period(num0, num1, j);        // dv/dx
            const double vit_imp = 0.5 * (vit.val_imp_face_bord_private(num0, j) + vit.val_imp_face_bord_private(num1, j));                // tangential velocity

            double temp2;

            if (n_type == 0 && contribution_paroi == 0)
              temp2 = 0;
            else
              temp2 = -signe * (vitesse[num2] - vit_imp) / domaine_VDF.dist_norm_bord(num1);

            element[0] = face_voisins(num2, 0);
            element[1] = face_voisins(num2, 1);

            for (int k = 0; k < 2; k++)
              if (element[k] != -1)
                SMA_barre[element[k]] += (temp1 + temp2) * (temp1 + temp2) * 0.25;
          }
    }

  for (int num_arete = prem_am; num_arete < dern_am; num_arete++)
    {
      const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
      const int i = orientation(num0), j = orientation(num2);

      const double temp1 = (vitesse[num1] - vitesse[num0]) / domaine_VDF.dist_face_period(num0, num1, j);          // dv/dx
      const double temp2 = (vitesse[num3] - vitesse[num2]) / domaine_VDF.dist_face_period(num2, num3, i);          // du/dy

      element[0] = face_voisins(num0, 0);
      element[1] = face_voisins(num0, 1);
      element[2] = face_voisins(num1, 0);
      element[3] = face_voisins(num1, 1);

      for (int k = 0; k < 4; k++)
        if (element[k] != -1)
          SMA_barre[element[k]] += (temp1 + temp2) * (temp1 + temp2) * 0.25;
    }

  for (int num_arete = prem_ai; num_arete < dern_ai; num_arete++)
    {
      const int num0 = Qdm(num_arete, 0), num1 = Qdm(num_arete, 1), num2 = Qdm(num_arete, 2), num3 = Qdm(num_arete, 3);
      const int i = orientation(num0), j = orientation(num2);

      const double temp1 = (vitesse[num1] - vitesse[num0]) / domaine_VDF.dist_face_period(num0, num1, j);            // dv/dx
      const double temp2 = (vitesse[num3] - vitesse[num2]) / domaine_VDF.dist_face_period(num2, num3, i);            // du/dy

      element[0] = face_voisins(num0, 0);
      element[1] = face_voisins(num0, 1);
      element[2] = face_voisins(num1, 0);
      element[3] = face_voisins(num1, 1);

      for (int k = 0; k < 4; k++)
        SMA_barre[element[k]] += (temp1 + temp2) * (temp1 + temp2) * 0.25;
    }

  // 2nd part: loop over elements and fill Sij for the non-cross derivatives (du/dx and dv/dy)

  for (int elem = 0; elem < nb_elem; elem++)
    {
      for (int i = 0; i < dimension; i++)
        {
          double temp1 = (vitesse[elem_faces(elem, i)] - vitesse[elem_faces(elem, i + dimension)]) / domaine_VDF.dim_elem(elem, orientation(elem_faces(elem, i)));
          SMA_barre(elem) += 2.0 * temp1 * temp1;
        }
    }

  // We take the square root!!!!!  WARNING SMA_barre=squared invariant!!!
  //  racine_carree(SMA_barre)
  return SMA_barre;
}

DoubleVect& Champ_Face_VDF::calcul_S_barre(const DoubleTab& vitesse, DoubleVect& SMA_barre, const Domaine_Cl_VDF& domaine_Cl_VDF) const
{
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const int nb_elem_tot = domaine_VDF.nb_elem_tot();
  const int nb_elem = domaine_VDF.nb_elem(), N = vitesse.line_size();

  int i, j;
  int elem;
  double Sij, temp;

  DoubleTab duidxj(nb_elem_tot, dimension, dimension, N);

  calcul_duidxj(vitesse, duidxj, domaine_Cl_VDF);

  for (elem = 0; elem < nb_elem; elem++)
    for (int n=0; n<N; n++)
      {
        temp = 0.;
        for (i = 0; i < dimension; i++)
          for (j = 0; j < dimension; j++)
            {
              Sij = 0.5 * (duidxj(elem, i, j, n) + duidxj(elem, j, i, n));
              temp += Sij * Sij;
            }
        SMA_barre(elem) = 2. * temp;
      }

  return SMA_barre;

}

DoubleTab& Champ_Face_VDF::calcul_S_barre_Multiphase(const DoubleTab& vitesse, DoubleTab& SMA_barre, const Domaine_Cl_VDF& domaine_Cl_VDF) const
{
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const int nb_elem_tot = domaine_VDF.nb_elem_tot();
  const int nb_elem = domaine_VDF.nb_elem();
  const int N = vitesse.line_size();

  int i, j;
  int elem;
  double Sij, temp;

  DoubleTab duidxj(nb_elem_tot, dimension, dimension, N);

  calcul_duidxj(vitesse, duidxj, domaine_Cl_VDF);

  for (elem = 0; elem < nb_elem; elem++)
    for (int n=0; n<N; n++)
      {
        temp = 0.;
        for (i = 0; i < dimension; i++)
          for (j = 0; j < dimension; j++)
            {
              Sij = 0.5 * (duidxj(elem, i, j, n) + duidxj(elem, j, i, n));
              temp += Sij * Sij;
            }
        SMA_barre(elem,n) = 2. * temp;
      }

  return SMA_barre;

}

void Champ_Face_VDF::calcul_grad_u(const DoubleTab& vitesse, DoubleTab& grad_u, const Domaine_Cl_VDF& domaine_Cl_VDF)
{
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const int nb_elem = domaine_VDF.nb_elem();
  const int nb_elem_tot = domaine_VDF.nb_elem_tot(), N = vitesse.line_size();

  DoubleTab gradient_elem(nb_elem_tot, dimension, dimension, N);
  gradient_elem = 0.;

  calcul_duidxj(vitesse, gradient_elem, domaine_Cl_VDF);

  for (int elem = 0; elem < nb_elem; elem++)
    for (int n=0; n<N; n++)
      for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
          grad_u(elem, N * ( dimension*i+j ) + n) = gradient_elem(elem, i, j, n);
}

void Champ_Face_VDF::calculer_dscald_centre_element(DoubleTab& dscald) const
{
  const DoubleTab& val = valeurs();
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  int nb_elem = domaine_VDF.nb_elem();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const IntTab& elem_faces = domaine_VDF.elem_faces();

  if (dimension == 2)
    caldscaldcentelemdim2(dscald, val, domaine_VDF, nb_elem, face_voisins, elem_faces);
  else if (dimension == 3)
    caldscaldcentelemdim3(dscald, val, domaine_VDF, nb_elem, face_voisins, elem_faces);
}

// Functions to compute the components of the GradU tensor (covariant velocity derivatives)
// in the cylindrical coordinate frame
void Champ_Face_VDF::calculer_dercov_axi(const Domaine_Cl_VDF& domaine_Cl_VDF)
{
  const Domaine_VDF& domaine_VDF = domaine_vdf();
  const DoubleTab& inco = valeurs();
  const IntVect& orientation = domaine_VDF.orientation();
  const IntTab& elem_faces = domaine_VDF.elem_faces();
  const IntTab& Qdm = domaine_VDF.Qdm();
  const DoubleTab& xv = domaine_VDF.xv();
  const DoubleTab& xp = domaine_VDF.xp();
  const IntVect& type_arete_bord = domaine_Cl_VDF.type_arete_bord();

  double d_teta, R;
  double deux_pi = M_PI * 2.0;

  // Fill tau_diag_: diagonal terms of the GradU tensor

  int fx0, fx1, fy0, fy1;
  int num_elem;
  for (num_elem = 0; num_elem < domaine_VDF.nb_elem(); num_elem++)
    {
      fx0 = elem_faces(num_elem, 0);
      fx1 = elem_faces(num_elem, dimension);
      fy0 = elem_faces(num_elem, 1);
      fy1 = elem_faces(num_elem, 1 + dimension);

      // Compute tau11
      tau_diag_(num_elem, 0) = (inco[fx1] - inco[fx0]) / (xv(fx1, 0) - xv(fx0, 0));

      // Compute tau22
      R = xp(num_elem, 0);
      d_teta = xv(fy1, 1) - xv(fy0, 1);
      if (d_teta < 0)
        d_teta += deux_pi;
      tau_diag_(num_elem, 1) = (inco[fy1] - inco[fy0]) / (R * d_teta) + 0.5 * (inco[fx0] + inco[fx1]) / R;
    }

  if (dimension == 3)
    {
      int fz0, fz1;
      for (num_elem = 0; num_elem < domaine_VDF.nb_elem(); num_elem++)
        {
          fz0 = elem_faces(num_elem, 2);
          fz1 = elem_faces(num_elem, 2 + dimension);

          // Compute tau33
          tau_diag_(num_elem, 2) = (inco[fz1] - inco[fz0]) / (xv(fz1, 2) - xv(fz0, 2));

        }
    }

  // Fill tau_croises_: off-diagonal terms of the GradU tensor
  // The cross derivatives of the velocity (off-diagonal terms of the
  // GradU tensor) are computed on edges.
  // There are two derivatives per edge:
  //    For an XY edge: tau12 and tau21
  //    For a YZ edge: tau23 and tau32
  //    For an XZ edge: tau13 and tau31

  // Loop over boundary edges

  int n_arete;
  int ndeb = domaine_VDF.premiere_arete_bord();
  int nfin = ndeb + domaine_VDF.nb_aretes_bord();
  int ori1, ori3;
  int fac1, fac2, fac3, fac4, signe;
  double dist3;

  int n_type;

  for (n_arete = ndeb; n_arete < nfin; n_arete++)
    {
      n_type = type_arete_bord(n_arete - ndeb);

      switch(n_type)
        {
        case TypeAreteBordVDF::PAROI_PAROI:
          // wall-wall
        case TypeAreteBordVDF::FLUIDE_FLUIDE:
          // fluid-fluid
        case TypeAreteBordVDF::PAROI_FLUIDE:
          // wall-fluid
          {
            fac1 = Qdm(n_arete, 0);
            fac2 = Qdm(n_arete, 1);
            fac3 = Qdm(n_arete, 2);
            signe = Qdm(n_arete, 3);
            ori1 = orientation(fac1);
            ori3 = orientation(fac3);
            int rang1 = fac1 - domaine_VDF.premiere_face_bord();
            int rang2 = fac2 - domaine_VDF.premiere_face_bord();
            double vit_imp;

            if (n_type == TypeAreteBordVDF::PAROI_FLUIDE)
              // wall-fluid edge: need to determine which is the fluid face
              {
                if (est_egal(inco[fac1], 0))
                  vit_imp = val_imp_face_bord_private(rang2, ori3);
                else
                  vit_imp = val_imp_face_bord_private(rang1, ori3);
              }
            else
              vit_imp = 0.5 * (val_imp_face_bord_private(rang1, ori3) + val_imp_face_bord_private(rang2, ori3));

            if (ori1 == 0) // boundary with equation R = const
              {
                dist3 = xv(fac3, 0) - xv(fac1, 0);
                if (n_type != TypeAreteBordVDF::PAROI_PAROI)
                  dist3 *= 2;

                if (ori3 == 1)
                  {
                    // compute tau12
                    tau_croises_(n_arete, 0) = signe * (vit_imp - inco[fac3]) / dist3;

                    // compute tau21
                    R = xv(fac1, 0);
                    d_teta = xv(fac2, 1) - xv(fac1, 1);
                    if (d_teta < 0)
                      d_teta += deux_pi;
                    tau_croises_(n_arete, 1) = (inco[fac2] - inco[fac1]) / (R * d_teta);
                  }
                else if (ori3 == 2)
                  {
                    // compute tau13
                    tau_croises_(n_arete, 0) = signe * (vit_imp - inco[fac3]) / dist3;
                    // compute tau31
                    tau_croises_(n_arete, 1) = (inco[fac2] - inco[fac1]) / (xv(fac2, 2) - xv(fac1, 2));
                  }
              }
            else if (ori1 == 1) // boundary with equation theta = const
              {
                R = xv(fac3, 0);
                d_teta = xv(fac3, 1) - xv(fac1, 1);
                if (d_teta < 0)
                  d_teta += deux_pi;
                dist3 = R * d_teta;
                if (n_type != TypeAreteBordVDF::PAROI_PAROI)
                  dist3 *= 2;

                if (ori3 == 0)
                  {
                    // compute tau21
                    tau_croises_(n_arete, 0) = signe * (vit_imp - inco[fac3]) / dist3 - 0.5 * (inco[fac1] + inco[fac2]) / R;
                    // compute tau12
                    tau_croises_(n_arete, 1) = (inco[fac2] - inco[fac1]) / (xv(fac2, 0) - xv(fac1, 0));
                  }
                else if (ori3 == 2)
                  {
                    // compute tau23
                    tau_croises_(n_arete, 0) = signe * (vit_imp - inco[fac3]) / dist3;
                    // compute tau32
                    tau_croises_(n_arete, 1) = (inco[fac2] - inco[fac1]) / (xv(fac2, 2) - xv(fac1, 2));
                  }
              }
            else // (ori1 == 2) boundary with equation Z = const
              {
                dist3 = xv(fac3, 2) - xv(fac1, 2);
                if (n_type != TypeAreteBordVDF::PAROI_PAROI)
                  dist3 *= 2;

                if (ori3 == TypeAreteBordVDF::PAROI_PAROI)
                  {
                    // compute tau31
                    tau_croises_(n_arete, 0) = signe * (vit_imp - inco[fac3]) / dist3;
                    // compute tau13
                    tau_croises_(n_arete, 1) = (inco[fac2] - inco[fac1]) / (xv(fac2, 0) - xv(fac1, 0));
                  }
                else if (ori3 == 1)
                  {
                    // compute tau32
                    tau_croises_(n_arete, 0) = signe * (vit_imp - inco[fac3]) / dist3;

                    // compute tau23
                    R = xv(fac1, 0);
                    d_teta = xv(fac2, 1) - xv(fac1, 1);
                    if (d_teta < 0)
                      d_teta += deux_pi;
                    tau_croises_(n_arete, 1) = (inco[fac2] - inco[fac1]) / (R * d_teta);
                  }
              }
            break;
          }
        case 3:
          {
            // symetrie-symetrie
            // pas de flux diffusif calcule
            break;
          }
        default:
          {
            Cerr << "On a rencontre un type d'arete non prevu\n";
            Cerr << "num arete : " << n_arete;
            Cerr << " type : " << n_type;
            exit();
            break;
          }
        }
    }

  // Loop over mixed and internal edges
  ndeb = domaine_VDF.premiere_arete_mixte();
  nfin = domaine_VDF.nb_aretes();
  for (n_arete = ndeb; n_arete < nfin; n_arete++)
    {
      fac1 = Qdm(n_arete, 0);
      fac2 = Qdm(n_arete, 1);
      fac3 = Qdm(n_arete, 2);
      fac4 = Qdm(n_arete, 3);
      ori1 = orientation(fac1);
      ori3 = orientation(fac3);
      if (ori1 == 1)  // (only possibility: ori3 = 0)  XY edge
        {
          // Compute tau21
          R = xv(fac3, 0);
          d_teta = xv(fac4, 1) - xv(fac3, 1);
          if (d_teta < 0)
            d_teta += deux_pi;
          tau_croises_(n_arete, 1) = (inco(fac4) - inco(fac3)) / (R * d_teta) - 0.5 * (inco[fac1] + inco[fac2]) / R;
          // Compute tau12
          tau_croises_(n_arete, 0) = (inco(fac2) - inco(fac1)) / (xv(fac2, 0) - xv(fac1, 0));
        }
      else if (ori3 == 1) // (only possibility: ori1 = 2) YZ edge
        {
          // Compute tau32
          tau_croises_(n_arete, 1) = (inco(fac4) - inco(fac3)) / (xv(fac4, 2) - xv(fac3, 2));
          // Compute tau23
          R = xv(fac1, 0);
          d_teta = xv(fac2, 1) - xv(fac1, 1);
          if (d_teta < 0)
            d_teta += deux_pi;
          tau_croises_(n_arete, 0) = (inco(fac2) - inco(fac1)) / (R * d_teta);
        }
      else // only possibility: ori1 = 2 and ori3 = 0: XZ edge
        {
          // Compute tau31
          tau_croises_(n_arete, 1) = (inco(fac4) - inco(fac3)) / (xv(fac4, 2) - xv(fac3, 2));
          // Compute tau13
          tau_croises_(n_arete, 0) = (inco(fac2) - inco(fac1)) / (xv(fac2, 0) - xv(fac1, 0));
        }
    }
}

/* ***************************************************** */
/*           USEFUL METHODS OUTSIDE THE CLASS            */
/* ***************************************************** */

double Champ_Face_get_val_imp_face_bord_sym(const DoubleTab& tab_valeurs, const double temp, int face, int comp, const Domaine_Cl_VDF& zclo)
{
  const Domaine_VDF& domaine_vdf = zclo.domaine_VDF();
  int face_locale = -123;
  const int face_globale = face + domaine_vdf.premiere_face_bord(); // Now index in the global face array.
  const Domaine_Cl_dis_base& zcl = zclo; //equation().domaine_Cl_dis();
  // Retrieve the BC associated with the face and the local face index within the boundary.
  //assert(equation().domaine_Cl_dis()==zclo);

  const Cond_lim_base& cl = (face < domaine_vdf.nb_faces()) ? zcl.condition_limite_de_la_face_reelle(face_globale, face_locale) :
                            zcl.condition_limite_de_la_face_virtuelle(face_globale, face_locale);

  const IntTab& face_voisins = domaine_vdf.face_voisins();
  const IntTab& elem_faces = domaine_vdf.elem_faces();
  const DoubleVect& porosite = zclo.equation().milieu().porosite_face();
  const int ori = domaine_vdf.orientation()(face_globale);

  if (sub_type(Navier, cl))
    {
      const int N = tab_valeurs.line_size();
      const int n=comp%N, comploc = (comp-n)/N;
      if (comploc == ori)
        return 0;
      else
        {
          int elem = 0;
          if (face_voisins(face_globale, 0) != -1)
            elem = face_voisins(face_globale, 0);
          else
            elem = face_voisins(face_globale, 1);
          const int comp2 = comploc + Objet_U::dimension;
          return (tab_valeurs(elem_faces(elem, comploc), n) * porosite[elem_faces(elem, comploc)] + tab_valeurs(elem_faces(elem, comp2), n) * porosite[elem_faces(elem, comp2)])
                 / (porosite[elem_faces(elem, comploc)] + porosite[elem_faces(elem, comp2)]);
        }
    }

  if (!cl.champ_front().has_valeurs_au_temps(temp)) // not yet initialized
    return 0.;

  const DoubleTab& vals = cl.champ_front().valeurs_au_temps(temp);
  const int face_de_vals = vals.dimension(0) == 1 ? 0 : face_locale;

  if (sub_type(Dirichlet_entree_fluide, cl))
    return vals(face_de_vals, comp);
  else if (sub_type(Dirichlet_paroi_fixe, cl))
    return 0.;
  else if (sub_type(Dirichlet_paroi_defilante, cl))
    return vals(face_de_vals, comp);

  return 0.; // All other cases
}

double Champ_Face_get_val_imp_face_bord(const double temp, int face, int comp, const Domaine_Cl_VDF& zclo)
{
  const Domaine_VDF& domaine_vdf = zclo.domaine_VDF();
  int face_locale = -123;
  const int face_globale = face + domaine_vdf.premiere_face_bord(); // Now index in the global face array.
  const Domaine_Cl_dis_base& zcl = zclo; //equation().domaine_Cl_dis();
  // Retrieve the BC associated with the face and the local face index within the boundary.
  //assert(equation().domaine_Cl_dis()==zclo);

  const Cond_lim_base& cl = (face < domaine_vdf.nb_faces()) ? zcl.condition_limite_de_la_face_reelle(face_globale, face_locale) :
                            zcl.condition_limite_de_la_face_virtuelle(face_globale, face_locale);
  const int ori = domaine_vdf.orientation()(face_globale);

  if (sub_type(Navier, cl))
    {
      if (comp == ori)
        return 0.;
      else
        {
          Process::exit("You should call Champ_Face_get_val_imp_face_bord_sym and not Champ_Face_get_val_imp_face_bord\n");
          return 1.e9;
        }
    }

  if (!cl.champ_front().has_valeurs_au_temps(temp)) // not yet initialized
    return 0.;

  const DoubleTab& vals = cl.champ_front().valeurs_au_temps(temp);
  int face_de_vals = vals.dimension(0) == 1 ? 0 : face_locale;

  if (sub_type(Dirichlet_entree_fluide, cl))
    return vals(face_de_vals, comp);
  else if (sub_type(Dirichlet_paroi_fixe, cl))
    return 0.;
  else if (sub_type(Dirichlet_paroi_defilante, cl))
    return vals(face_de_vals, comp);

  return 0.; // All other cases
}

double Champ_Face_get_val_imp_face_bord(const double temp, int face, int comp, int comp2, const Domaine_Cl_VDF& zclo)
{
  Process::exit("Champ_Face_VDF::val_imp_face_bord(,,) exit\n");
  return 0.; // For compilers
}

double Champ_Face_coeff_frottement_face_bord(const int f, const int n, const Domaine_Cl_VDF& zclo)
{
  const Domaine_VDF& domaine_vdf = zclo.domaine_VDF();
  const Domaine_Cl_dis_base& zcl = zclo;
  const int face_globale = f + domaine_vdf.premiere_face_bord(); // Now index in the global face array.

  int face_locale = -123;
  const Cond_lim_base& cl = (f < domaine_vdf.nb_faces()) ? zcl.condition_limite_de_la_face_reelle(face_globale, face_locale) :
                            zcl.condition_limite_de_la_face_virtuelle(face_globale, face_locale);

  return sub_type(Navier, cl) ? ref_cast(Navier, cl).coefficient_frottement(face_locale,n) : 0.;
}

double Champ_Face_coeff_frottement_grad_face_bord(const int f, const int n, const Domaine_Cl_VDF& zclo)
{
  const Domaine_VDF& domaine_vdf = zclo.domaine_VDF();
  const Domaine_Cl_dis_base& zcl = zclo;
  const int face_globale = f + domaine_vdf.premiere_face_bord(); // Now index in the global face array.

  int face_locale = -123;
  const Cond_lim_base& cl = (f < domaine_vdf.nb_faces()) ? zcl.condition_limite_de_la_face_reelle(face_globale, face_locale) :
                            zcl.condition_limite_de_la_face_virtuelle(face_globale, face_locale);

  return sub_type(Navier, cl) ? ref_cast(Navier, cl).coefficient_frottement_grad(face_locale,n) : 0.;
}
