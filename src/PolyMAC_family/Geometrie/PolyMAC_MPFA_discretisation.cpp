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

#include <Champ_Fonc_Elem_PolyMAC_MPFA_TC.h>
#include <Champ_Fonc_Elem_PolyMAC_MPFA_rot.h>
#include <Champ_Fonc_Tabule_Elem_PolyMAC_CDO.h>
#include <grad_Champ_Face_PolyMAC_MPFA.h>
#include <PolyMAC_MPFA_discretisation.h>
#include <Champ_Fonc_Elem_PolyMAC_CDO.h>
#include <Domaine_Cl_PolyMAC_family.h>
#include <Domaine_PolyMAC_MPFA.h>
#include <Champ_Fonc_Tabule.h>
#include <Schema_Temps_base.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <DescStructure.h>
#include <Milieu_base.h>
#include <Motcle.h>

Implemente_instanciable(PolyMAC_MPFA_discretisation, "PolyMAC_MPFA|PolyMAC_P0", PolyMAC_HFV_discretisation);
// XD PolyMAC_MPFA discretisation_base PolyMAC_P0 INHERITS_BRACE PolyMAC_MPFA discretization (previously covimac
// XD_CONT discretization compatible with pb_multi).

Entree& PolyMAC_MPFA_discretisation::readOn(Entree& s) { return s; }

Sortie& PolyMAC_MPFA_discretisation::printOn(Sortie& s) const { return s; }

/**
 * @brief Creates a velocity gradient field for PolyMAC_MPFA discretization
 *
 * This method creates and initializes a gradient field of the velocity field,
 * computing the tensor ∇u where u is the velocity field. The resulting field
 * contains all components of the velocity gradient tensor.
 *
 * @param z           Discretized domain (must be Domaine_PolyMAC_MPFA)
 * @param zcl         Discretized boundary conditions domain (must be Domaine_Cl_PolyMAC_family)
 * @param ch_vitesse  Velocity field (must be Champ_Face_PolyMAC_MPFA)
 * @param[out] ch     Output gradient field (will be typed as grad_Champ_Face_PolyMAC_MPFA)
 *
 * @details The gradient field components are named according to spatial dimension:
 *          - 2D: dU_phase, dV_phase for each phase
 *          - 3D: dU_phase, dV_phase, dW_phase for each phase
 *
 *          Where phase corresponds to each line component of the velocity field.
 *
 * @see grad_Champ_Face_PolyMAC_MPFA for the actual gradient field implementation
 */
void PolyMAC_MPFA_discretisation::grad_u(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  const Champ_Face_PolyMAC_MPFA& vit = ref_cast(Champ_Face_PolyMAC_MPFA, ch_vitesse);
  const Domaine_PolyMAC_MPFA& domaine_poly = ref_cast(Domaine_PolyMAC_MPFA, z);
  const Domaine_Cl_PolyMAC_family& domaine_cl_poly = ref_cast(Domaine_Cl_PolyMAC_family, zcl);

  ch.typer("grad_Champ_Face_PolyMAC_MPFA");

  grad_Champ_Face_PolyMAC_MPFA& ch_grad_u = ref_cast(grad_Champ_Face_PolyMAC_MPFA, ch.valeur()); //

  ch_grad_u.associer_domaine_dis_base(domaine_poly);
  ch_grad_u.associer_domaine_Cl_dis_base(domaine_cl_poly);
  ch_grad_u.associer_champ(vit);
  ch_grad_u.nommer("gradient_vitesse");
  ch_grad_u.fixer_nb_comp(dimension * dimension * vit.valeurs().line_size());

  for (int n = 0; n < ch_grad_u.valeurs().line_size(); n++)
    {
      Nom phase = Nom(n);
      if (dimension == 2)
        {
          ch_grad_u.fixer_nom_compo(dimension * n + 0, Nom("dU_") + phase); // dU
          ch_grad_u.fixer_nom_compo(dimension * n + 1, Nom("dV_") + phase); // dV
        }
      else
        {
          ch_grad_u.fixer_nom_compo(dimension * n + 0, Nom("dU_") + phase); // dU
          ch_grad_u.fixer_nom_compo(dimension * n + 1, Nom("dV_") + phase); // dV
          ch_grad_u.fixer_nom_compo(dimension * n + 2, Nom("dW_") + phase); // dW
        }
    }
  ch_grad_u.fixer_nature_du_champ(multi_scalaire); // tensoriel pour etre precis
  ch_grad_u.fixer_nb_valeurs_nodales(-1);
  ch_grad_u.fixer_unite("s-1");
  ch_grad_u.changer_temps(-1); // so it is calculated at time 0
}

/**
 * @brief Creates a shear rate field for PolyMAC_MPFA discretization
 *
 * This method creates and initializes a field representing the shear rate magnitude,
 * which is derived from the velocity gradient tensor. The shear rate is computed
 * as the magnitude of the strain rate tensor.
 *
 * @param z           Discretized domain (must be Domaine_PolyMAC_MPFA)
 * @param zcl         Discretized boundary conditions domain
 * @param ch_vitesse  Velocity field (must be Champ_Face_PolyMAC_MPFA)
 * @param[out] ch     Output shear rate field (will be typed as Champ_Fonc_Elem_PolyMAC_MPFA_TC)
 *
 * @details The shear rate field contains one scalar component per phase,
 *          named "Taux_cisaillement_phase" where phase is the phase index.
 */
void PolyMAC_MPFA_discretisation::taux_cisaillement(const Domaine_dis_base& z, const Domaine_Cl_dis_base& zcl, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  const Champ_Face_PolyMAC_MPFA& vit = ref_cast(Champ_Face_PolyMAC_MPFA, ch_vitesse);
//  const Domaine_PolyMAC_MPFA&          domaine_poly = ref_cast(Domaine_PolyMAC_MPFA, z);
  const Domaine_PolyMAC_MPFA& domaine = ref_cast(Domaine_PolyMAC_MPFA, vit.domaine_dis_base());

  ch.typer("Champ_Fonc_Elem_PolyMAC_MPFA_TC");
  Champ_Fonc_Elem_PolyMAC_MPFA_TC& ch_grad_u = ref_cast(Champ_Fonc_Elem_PolyMAC_MPFA_TC, ch.valeur()); //

  ch_grad_u.associer_domaine_dis_base(domaine);
  ch_grad_u.associer_champ(vit);
  ch_grad_u.nommer("Taux_cisaillement");
  int N = vit.valeurs().line_size();
  ch_grad_u.fixer_nb_comp(N);
  for (int n = 0; n < N; n++)
    {
      Nom phase = Nom(n);
      ch_grad_u.fixer_nom_compo(n, Nom("Taux_cisaillement_") + phase);
    }
  ch_grad_u.fixer_nature_du_champ(scalaire); // tensoriel pour etre precis
  ch_grad_u.fixer_nb_valeurs_nodales(domaine.nb_elem());
  ch_grad_u.fixer_unite("s-1");
  ch_grad_u.changer_temps(-1); // so it is calculated at time 0
}

/**
 * @brief Creates a vorticity field for PolyMAC_MPFA discretization
 *
 * This method creates and initializes a vorticity field, which represents
 * the curl of the velocity field (∇ × u). The field structure depends on
 * the spatial dimension.
 *
 * @param sch         Time scheme (used for temporal integration context)
 * @param ch_vitesse  Velocity field (must be Champ_Face_PolyMAC_MPFA)
 * @param[out] ch     Output vorticity field (will be typed as Champ_Fonc_Elem_PolyMAC_MPFA_rot)
 *
 * @details Vorticity field structure by dimension:
 *          - 2D: Scalar field with N components (one per phase)
 *                Represents the z-component of curl(u)
 *          - 3D: Vector field with 3×N components
 *                Represents full curl(u) = (ωx, ωy, ωz) for each phase
 *
 * @note Component naming convention:
 *       - Single phase (N=1): "vorticitex", "vorticitey", "vorticitez" for 3D
 *       - Multiple phases: "vorticite_0", "vorticite_1", etc. with offset indexing
 */
void PolyMAC_MPFA_discretisation::creer_champ_vorticite(const Schema_Temps_base& sch, const Champ_Inc_base& ch_vitesse, OWN_PTR(Champ_Fonc_base) &ch) const
{
  const Champ_Face_PolyMAC_MPFA& vit = ref_cast(Champ_Face_PolyMAC_MPFA, ch_vitesse);
  const Domaine_PolyMAC_MPFA& domaine = ref_cast(Domaine_PolyMAC_MPFA, vit.domaine_dis_base());
  int N = vit.valeurs().line_size();

  ch.typer("Champ_Fonc_Elem_PolyMAC_MPFA_rot");
  Champ_Fonc_Elem_PolyMAC_MPFA_rot& ch_rot_u = ref_cast(Champ_Fonc_Elem_PolyMAC_MPFA_rot, ch.valeur());

  ch_rot_u.associer_domaine_dis_base(domaine);
  ch_rot_u.associer_champ(vit);
  ch_rot_u.nommer("vorticite");

  if (dimension == 2)
    {
      ch_rot_u.fixer_nature_du_champ(scalaire);
      ch_rot_u.fixer_nb_comp(N);
    }
  else if (dimension == 3)
    {
      ch_rot_u.fixer_nature_du_champ(vectoriel);
      ch_rot_u.fixer_nb_comp(dimension * N);
      std::string nom_compo = "vorticite";
      if (N == 1)
        {
          std::vector<std::string> dir = { "x", "y", "z" };
          for (int d = 0; d < dimension; d++)
            ch_rot_u.fixer_nom_compo(d, nom_compo + dir[d]);
        }
      else
        {
          for (int n = 0; n < N; n++)
            {
              for (int d = 0; d < dimension; d++)
                {
                  int offset = dimension * n + d;
                  ch_rot_u.fixer_nom_compo(offset, nom_compo + "_" + std::to_string(offset));
                }
            }
        }
    }
  else
    abort();

  ch_rot_u.fixer_nb_valeurs_nodales(domaine.nb_elem());
  ch_rot_u.fixer_unite("s-1");
  ch_rot_u.changer_temps(-1); // so it is calculated at time 0
}

/**
 * @brief Creates a residue field for equation solving diagnostics
 *
 * This method creates a field to store the residue of an equation system,
 * which is useful for monitoring convergence and numerical solution quality.
 * The residue represents the difference between the left and right hand sides
 * of the discretized equation.
 *
 * @param z        Discretized domain
 * @param ch_inco  Unknown field for which the residue is computed
 * @param[out] champ Output residue field
 *
 * @note For non-face fields, the method delegates to the parent class
 *       PolyMAC_HFV_discretisation::residu()
 */
void PolyMAC_MPFA_discretisation::residu(const Domaine_dis_base& z, const Champ_Inc_base& ch_inco, OWN_PTR(Champ_Fonc_base) &champ) const
{
  Nom ch_name(ch_inco.le_nom());
  ch_name += "_residu";
  Cerr << "Discretization of " << ch_name << finl;

  Nom type_ch = ch_inco.que_suis_je();
  if (type_ch.debute_par("Champ_Face"))
    {
      Motcle loc = "champ_face";
      Noms nom(1), unites(1);
      nom[0] = ch_name;
      unites[0] = "units_not_defined";
      int nb_comp = ch_inco.valeurs().line_size() * dimension;

      discretiser_champ(loc, z, vectoriel, nom, unites, nb_comp, ch_inco.temps(), champ);

      Champ_Fonc_base& ch_fonc = ref_cast(Champ_Fonc_base, champ.valeur());
      DoubleTab& tab = ch_fonc.valeurs();
      tab = -10000.0;
      Cerr << "[Information] Discretisation_base::residu : the residue is set to -10000.0 at initial time" << finl;
    }

  else
    PolyMAC_HFV_discretisation::residu(z, ch_inco, champ);
}
