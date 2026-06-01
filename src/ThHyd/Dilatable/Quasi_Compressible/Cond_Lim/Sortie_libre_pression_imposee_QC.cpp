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

#include <Sortie_libre_pression_imposee_QC.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <Motcle.h>

Implemente_instanciable_sans_constructeur(Sortie_libre_pression_imposee_QC, "Frontiere_ouverte_pression_totale_imposee", Neumann_sortie_libre);

Sortie_libre_pression_imposee_QC::Sortie_libre_pression_imposee_QC() : Pthn(1.e5), d_rho(-1) { }

Sortie& Sortie_libre_pression_imposee_QC::printOn(Sortie& s) const { return s << que_suis_je() << finl; }

Entree& Sortie_libre_pression_imposee_QC::readOn(Entree& s)
{
  if (app_domains.size() == 0) app_domains = { Motcle("Hydraulique"), Motcle("indetermine") };

  s >> le_champ_front;
  le_champ_ext.typer("Champ_front_uniforme");
  le_champ_ext->valeurs().resize(1, dimension);
  return s;
}

/*! @brief Completes the boundary conditions.
 *
 * Sets d_rho to the constant density of the physical medium of the equation.
 *
 */
void Sortie_libre_pression_imposee_QC::completer()
{
  const Milieu_base& mil = mon_dom_cl_dis->equation().milieu();
  if (sub_type(Champ_Uniforme, mil.masse_volumique()))
    {
      const Champ_Uniforme& rho = ref_cast(Champ_Uniforme, mil.masse_volumique());
      d_rho = rho.valeurs()(0, 0);
    }
  else
    d_rho = -1;
}

/*! @brief Returns the imposed flux value for the i-th component of the field representing the flux at the boundary.
 *
 * The boundary field is considered constant over all elements of the boundary.
 * The imposed flux value equals the (constant) boundary field value divided by d_rho.
 *
 * @param i Index along the first dimension of the field.
 * @return Imposed value for the specified field component.
 * @throws Second dimension of the boundary field greater than 1.
 */
double Sortie_libre_pression_imposee_QC::flux_impose(int i) const
{
  const Milieu_base& mil = mon_dom_cl_dis->equation().milieu();
  const Champ_base& rho = mil.masse_volumique();
  double rho_;
  if (d_rho == -1)
    rho_ = rho.valeurs()(i);
  else
    rho_ = d_rho;

  if (le_champ_front->valeurs().size() == 1)
    return (le_champ_front->valeurs()(0, 0) - Pthn) / rho_;
  else if (le_champ_front->valeurs().dimension(1) == 1)
    return (le_champ_front->valeurs()(i, 0) - Pthn) / rho_;
  else
    Cerr << "Neumann::flux_impose error" << finl;

  Process::exit();
  return 0.;
}

/*! @brief Returns the imposed flux value for the (i,j)-th component of the field representing the flux at the boundary.
 *
 * The boundary field is NOT constant over all elements of the boundary.
 *
 * @param i Index along the first dimension of the field.
 * @param j Index along the second dimension of the field.
 * @return Imposed value for the specified field component.
 */
double Sortie_libre_pression_imposee_QC::flux_impose(int i, int j) const
{
  const Milieu_base& mil = mon_dom_cl_dis->equation().milieu();
  const Champ_base& rho = mil.masse_volumique();
  double rho_;

  if (d_rho == -1)
    rho_ = rho.valeurs()(i);
  else
    rho_ = d_rho;

  if (le_champ_front->valeurs().dimension(0) == 1)
    return (le_champ_front->valeurs()(0, j) - Pthn) / rho_;
  else
    return (le_champ_front->valeurs()(i, j) - Pthn) / rho_;
}
