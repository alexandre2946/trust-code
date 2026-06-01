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

#include <Sortie_libre_pression_imposee.h>
#include <Fluide_Incompressible.h>
#include <Navier_Stokes_std.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>

Implemente_instanciable_sans_constructeur(Sortie_libre_pression_imposee, "Frontiere_ouverte_pression_imposee", Neumann_sortie_libre);
// XD frontiere_ouverte_pression_imposee neumann frontiere_ouverte_pression_imposee INHERITS_BRACE Imposed pressure
// XD_CONT condition at the open boundary called bord (edge). The imposed pressure field is expressed in Pa.
// XD attr ch front_field_base ch REQ Boundary field type.


Sortie_libre_pression_imposee::Sortie_libre_pression_imposee() : d_rho(-123.) { }

Sortie& Sortie_libre_pression_imposee::printOn(Sortie& s) const { return s << que_suis_je() << finl; }

Entree& Sortie_libre_pression_imposee::readOn(Entree& s)
{
  if (app_domains.size() == 0) app_domains = { Motcle("Hydraulique"), Motcle("indetermine") };

  s >> le_champ_front;
  le_champ_ext.typer("Champ_front_uniforme");
  le_champ_ext->valeurs().resize(1, dimension);
  return s;
}

/*! @brief Completes the boundary conditions.
 *
 * Sets the constant density of the physical medium
 *     of the equation to d_rho.
 *
 */
void Sortie_libre_pression_imposee::completer()
{
  const Milieu_base& mil = mon_dom_cl_dis->equation().milieu();
  if (sub_type(Fluide_Incompressible,mil) && mon_dom_cl_dis->equation().que_suis_je() != "QDM_Multiphase")
    {
      if (sub_type(Champ_Uniforme, mil.masse_volumique()))
        {
          const Champ_Uniforme& rho = ref_cast(Champ_Uniforme, mil.masse_volumique());
          d_rho = rho.valeurs()(0, 0);
        }
      else
        {
          d_rho = -1;

          // GF: in non-constant rho cases (QC, FT) we must not divide by rho, so d_rho=1
          d_rho = 1;
        }
    }
  else
    d_rho = 1;
}

/*! @brief Returns the value of the imposed flux on the i-th component of the field representing the flux at the boundary.
 *
 *     The boundary field is considered constant over all elements of the boundary.
 *     The imposed flux value at the boundary equals the value of the (constant) boundary field divided by d_rho.
 *
 * @param (int i) index along the first dimension of the field
 * @return (double) the value imposed on the specified component of the field
 * @throws second dimension of the boundary field greater than 1
 */
double Sortie_libre_pression_imposee::flux_impose(int i) const
{
  return flux_impose(i,0);
}

/*! @brief Returns the value of the imposed flux on the (i,j)-th component of the field representing the flux at the boundary.
 *
 *     The boundary field is NOT constant over all elements
 *     of the boundary.
 *
 * @param (int i) index along the first dimension of the field
 * @param (int j) index along the second dimension of the field
 * @return (double) the value imposed on the specified component of the field
 */
double Sortie_libre_pression_imposee::flux_impose(int i, int j) const
{
  const Milieu_base& mil = mon_dom_cl_dis->equation().milieu();
  double rho_;
  assert(!est_egal(d_rho, -123.));
  if (d_rho == -1)
    {
      const Champ_base& rho = mil.masse_volumique();
      rho_ = rho.valeurs()(i);
    }
  else
    rho_ = d_rho;

  if (le_champ_front->valeurs().dimension(0) == 1)
    return le_champ_front->valeurs()(0, j) / rho_;
  else if (j < le_champ_front->valeurs().dimension(1))
    return le_champ_front->valeurs()(i, j) / rho_;
  else
    Cerr << "Sortie_libre_pression_imposee::flux_impose error" << finl;
  Process::exit();
  return 0.;
}
