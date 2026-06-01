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

#include <Terme_Source_Acceleration_VEF_Face.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Periodique.h>
#include <Navier_Stokes_std.h>
#include <Champ_Fonc_P0_VEF.h>
#include <Milieu_base.h>

Implemente_instanciable(Terme_Source_Acceleration_VEF_Face,"Acceleration_VEF_P1NC",Terme_Source_Acceleration);

Sortie& Terme_Source_Acceleration_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

/*! @brief Call Terme_Source_Acceleration::lire_data.
 *
 */
Entree& Terme_Source_Acceleration_VEF_Face::readOn(Entree& s )
{
  lire_data(s);
  return s;
}

/*! @brief Method called by Source_base::completer() after associer_domaines. Fills the references to the domains and domain_cl.
 *
 */
void Terme_Source_Acceleration_VEF_Face::associer_domaines(const Domaine_dis_base& domaine_dis,
                                                           const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  if (je_suis_maitre())
    Cerr << "Terme_Source_Acceleration_VEF_Face::associer_domaines" << finl;
  le_dom_VEF_    = ref_cast(Domaine_VEF, domaine_dis);
  le_dom_Cl_VEF_ = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis);
}

/*! @brief Helper function for Terme_Source_Acceleration_VEF_Face::ajouter. Adds contributions from a contiguous list of faces of the translation source term:
 *
 *    s_face = source_term * rho
 *    resu  += integral (s_face) over the velocity control volume.
 *   Handles the following cases:
 *     rho = null reference (=> rho = 1.)  otherwise rho != null
 *     boundary faces => free outlet
 *     periodicity
 *     symmetry because in VEF the velocity on a boundary face may not be zero (V_tangential)
 *     internal faces
 *
 */
static void TSAVEF_ajouter_liste_faces(const int premiere_face, const int derniere_face,
                                       const DoubleVect& volumes_entrelaces,
                                       const DoubleVect& volumes_elements,
                                       const DoubleVect& porosite_surf,
                                       const IntTab&      face_voisins,
                                       const OBS_PTR(Champ_base) & ref_rho,
                                       const DoubleTab&   terme_source,
                                       DoubleTab& s_face,
                                       DoubleTab& resu)
{
  int num_face;
  // Constant pointer to a constant array.
  // Null pointer if ref_rho_ is a null reference.
  const DoubleTab * const rho_elem =
    (bool(ref_rho)) ? &(ref_rho->valeurs()) : 0;
  const int dim = Objet_U::dimension;

  for (num_face=premiere_face; num_face<derniere_face; num_face++)
    {
      const double vol = volumes_entrelaces(num_face)*porosite_surf(num_face);
      double src[3] = {0., 0., 0.};
      int j;

      for (j = 0; j < dim; j++)
        src[j] = terme_source(num_face, j);

      double rho = 1.;

      // Compute a mean rho over the velocity control volume
      if (rho_elem)
        {
          const int elem0 = face_voisins(num_face,0);
          const int elem1 = face_voisins(num_face,1);
          double rho0 = 0, rho1 = 0, vol0 = 0, vol1 = 0;
          if (elem0 >= 0)
            {
              rho0 = (*rho_elem)(elem0);
              vol0 = volumes_elements(elem0);
            }
          if (elem1 >= 0)
            {
              rho1 = (*rho_elem)(elem1);
              vol1 = volumes_elements(elem1);
            }
          rho = (rho0 * vol0 + rho1 * vol1) / (vol0 + vol1);
        }

      for (j = 0; j < dim; j++)
        {
          double a = src[j] * rho;
          s_face(num_face, j) = a;
          // Integral over the control volume:
          resu(num_face, j) += a * vol;
        }
    }
}

/*! @brief Adds the term (la_source_ * rho * volume_entrelace) to the resu field.
 *
 * Assumes that resu is discretized like the velocity.
 *   The virtual space of "resu" is NOT updated!
 *
 *   Note on the scheme: the acceleration d/dt(v) is computed at elements,
 *    then multiplied by "rho" at elements, then evaluated at faces via an
 *    average over the neighbouring elements of the face. This is a first attempt,
 *    not necessarily the best approach. Since the acceleration depends on the velocity
 *    which is at faces, two successive interpolations are used.
 *  Side effect:
 *   (la_source_ * rho) is stored in terme_source_post_
 *
 */
DoubleTab& Terme_Source_Acceleration_VEF_Face::ajouter(DoubleTab& resu) const
{
  const Domaine_VF&     domaine               = le_dom_VEF_.valeur();
  const Domaine_Cl_dis_base& domaine_Cl       = le_dom_Cl_VEF_.valeur();
  const IntTab&      face_voisins       = domaine.face_voisins();
  const DoubleVect& porosite_surf      = equation().milieu().porosite_face();
  const DoubleVect& volumes_entrelaces = domaine.volumes_entrelaces();

  DoubleTab& s_face = get_set_terme_source_post().valeurs();
  s_face = 0.;

  // Compute la_source_ from the acceleration fields and the fluid velocity.
  const int dim     = Objet_U::dimension;
  const int nb_faces = resu.dimension(0);
  DoubleTab acceleration_aux_faces(nb_faces, dim);
  calculer_la_source(acceleration_aux_faces);

  // Loop over the boundary conditions to process the boundary faces

  for (int n_bord = 0; n_bord < domaine.nb_front_Cl(); n_bord++)
    {
      // for each boundary condition check its type
      // If Dirichlet face do nothing
      // If Neumann, Periodic or Symmetry face compute the contribution to the source term
      const Cond_lim& la_cl = domaine_Cl.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();

      TSAVEF_ajouter_liste_faces(ndeb, nfin,
                                 volumes_entrelaces,
                                 domaine.volumes(),
                                 porosite_surf,
                                 face_voisins,
                                 ref_rho_,
                                 acceleration_aux_faces,
                                 s_face,
                                 resu);

    }
  // Loop over internal faces
  {
    const int ndeb = domaine.premiere_face_int();
    const int nfin = domaine.nb_faces();
    TSAVEF_ajouter_liste_faces(ndeb, nfin,
                               volumes_entrelaces,
                               domaine.volumes(),
                               porosite_surf,
                               face_voisins,
                               ref_rho_,
                               acceleration_aux_faces,
                               s_face,
                               resu);
  }

  {
    // Enforce periodicity
    int nb_comp=resu.line_size();
    for (int n_bord=0; n_bord<domaine.nb_front_Cl(); n_bord++)
      {
        const Cond_lim& la_cl = domaine_Cl.les_conditions_limites(n_bord);
        if (sub_type(Periodique,la_cl.valeur()))
          {
            const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
            const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
            int nb_faces_bord=le_bord.nb_faces();
            ArrOfInt fait(nb_faces_bord);
            fait = 0;
            for (int ind_face=0; ind_face<nb_faces_bord; ind_face++)
              {
                if (fait[ind_face] == 0)
                  {
                    int ind_face_associee = la_cl_perio.face_associee(ind_face);
                    fait[ind_face] = 1;
                    fait[ind_face_associee] = 1;
                    int face = le_bord.num_face(ind_face);
                    int face_associee = le_bord.num_face(ind_face_associee);
                    for (int comp=0; comp<nb_comp; comp++)
                      {
                        double val = 0.5*(resu(face_associee, comp)+resu(face, comp));
                        resu(face, comp)=resu(face_associee, comp) = val;
                      }
                  }// if fait
              }// for face
          }// sub_type Perio
      }
  }
  return resu;
}

/*! @brief Computes the three-component velocity field at faces from the velocity field at faces of eq_hydraulique_.
 *
 * inconnue().
 *   In VEF: nothing to do. The provided storage is not used; N.S.inconnue() is returned directly.
 *
 * @param v_faces_stockage array to store the result if computations are needed. Not used in VEF.
 * @return reference to the velocity values.
 */
const DoubleTab& Terme_Source_Acceleration_VEF_Face::calculer_vitesse_faces(DoubleTab& v_faces_stockage) const
{
  return get_eq_hydraulique().inconnue().valeurs();
}

/*! @brief Associates the density field. The computed source term will then be homogeneous to d/dt(integral(rho*v)).
 *
 * @param champ_rho a field of type Champ_Fonc_P0_VEF that will be used in calls to "ajouter()" to evaluate the density.
 */

void Terme_Source_Acceleration_VEF_Face::associer_champ_rho(const Champ_base& champ_rho)
{
  // The field must be discretized at elements: possibility
  // to allow other types if needed (Champ_Don, Champ_Inc, etc.)
  // as long as they are P0 fields.
  if (!sub_type(Champ_Fonc_P0_VEF, champ_rho))
    {
      Cerr << "Error in Terme_Source_Acceleration_VEF_Face::associer_champ_rho" << finl;
      Cerr << " The density field must be of type Champ_Fonc_P0_VEF" << finl;
      Cerr << " Type of the associated field: " << champ_rho.que_suis_je() << finl;
      Cerr << " Name of the associated field: " << champ_rho.le_nom() << finl;
      assert(0);
      exit();
    }
  ref_rho_ = champ_rho;
}

