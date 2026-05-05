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

#ifndef Terme_Source_Qdm_Elem_DG_included
#define Terme_Source_Qdm_Elem_DG_included

#include <Terme_Source_Qdm.h>
#include <Source_base.h>
#include <TRUST_Ref.h>
#include <Champ_Fonc_Quad_DG.h>
#include <Probleme_base.h>

class Domaine_Cl_DG;
class Domaine_DG;

/**
 * @brief DG momentum source term for element-based velocity unknowns.
 *
 * This class adds a volumetric body force f(x) to the momentum right-hand side.
 * The contribution to the DG weak form for each element e and each velocity DOF (fb, d) is:
 *   secmem(e, fb + d*nb_bfunc) += integral of f_d(x) * phi_fb(x) dV
 * using a fixed order-5 quadrature rule.
 *
 * The source field is read from the input stream as a Champ_Don_base (la_source),
 * then projected onto the DG quadrature representation (la_source_DG) of type
 * Champ_Fonc_Quad_DG, which stores the field values pre-sampled at quadrature points.
 * Both uniform and spatially varying source fields are supported: for a Champ_Uniforme
 * the same value is used for all elements (index 0 in la_source_DG).
 *
 * The class has no matrix contribution (dimensionner_blocs() is a no-op).
 *
 * @sa Source_base, Terme_Source_Qdm
 */
class Terme_Source_Qdm_Elem_DG : public Source_base, public Terme_Source_Qdm
{
  Declare_instanciable(Terme_Source_Qdm_Elem_DG);
public:
  int has_interface_blocs() const override { return 1; }

  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override { } //rien
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;

  void associer_pb(const Probleme_base& ) override  { }
  void mettre_a_jour(double ) override;

protected:
  OBS_PTR(Domaine_DG) le_dom_DG;
  OBS_PTR(Domaine_Cl_DG) le_dom_Cl_DG;
  void associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base& ) override;

  OWN_PTR(Champ_Don_base) la_source_DG;
};

#endif /* Terme_Source_Qdm_Elem_DG_included */
