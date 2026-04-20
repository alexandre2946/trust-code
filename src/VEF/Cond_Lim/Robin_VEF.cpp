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

#include <Frontiere_dis_base.h>
#include <Navier_Stokes_std.h>
#include <Equation_base.h>
#include <Robin_VEF.h>
#include <Front_VF.h>
#include <Param.h>

Implemente_instanciable( Robin_VEF, "Robin_VEF", Cond_lim_base ) ;
// XD robin_vef condlim_base robin_vef 1 Robin condition at the boundary (edge)

Sortie& Robin_VEF::printOn( Sortie& os ) const
{
  return Cond_lim_base::printOn( os );
}

Entree& Robin_VEF::readOn( Entree& is )
{
  Param param(que_suis_je());

  param.ajouter("alpha", 		&alpha_robin_cl_, 		Param::REQUIRED); // XD_ADD_P floattant Robin coefficient for the normal field
  param.ajouter("beta" , 		&beta_robin_cl_ , 		Param::REQUIRED); // XD_ADD_P floattant Robin coefficient for the tangent field
  param.ajouter("champ_front_normal_et_tangentiel_robin",&le_champ_front,Param::REQUIRED); // XD_ADD_P front_field_base The boundary field
  param.lire_avec_accolades_depuis(is);


  Cerr << "Reading Robin boundary condition, the value of the normal Robin coefficient is  "
       << alpha_robin_cl_ << " and the tangential Robin coefficient is  " << beta_robin_cl_ << finl;

  if (alpha_robin_cl_<=0 || beta_robin_cl_<=0)
    Process::exit("Error of Robin boundary conditions : Alpha and Beta parameters should be positive. \n ") ;

  if (supp_discs.size() == 0)
    supp_discs = { Nom("VEF") , Nom("VEFPreP1B") };

  if (app_domains.size() == 0)
    app_domains = { Motcle("Hydraulique") };

  return is; //Cond_lim_base::readOn(is);
}

void Robin_VEF::completer()
{
  Cond_lim_base::completer();

  assert( mon_dom_cl_dis);

  if (!sub_type(Navier_Stokes_std, mon_dom_cl_dis->equation()))
    Process::exit("Robin_VEF is not yet coded on scalar equations ... \n");

  // pour le moment ....
  if (mon_dom_cl_dis->equation().operateur(0).l_op_base().que_suis_je() != "Op_Diff_VEF_P1NC")
    Process::exit("Robin_VEF is now only available for the laminar Op_Diff_VEF_P1NC operator ... \n");
}

/*! @brief Returns the value of the imposed normal and tangential Robin flux on the (i,j)-th component of the field representing the flux at the boundary.
 *
 * @param (int i) index along the first dimension of the field
 * @param (int j) index along the second dimension of the field (j=1 returns the normal component, and for j>1, returns the (j-1)-th tangential component)
 * @return (double) the imposed value on the specified field component
 */
double Robin_VEF::flux_robin_normal_et_trangentiel_imp(int i, int j) const
{
  // if dim = 2 : le_champ_front has 1+2 dimensions  : 1 column for the localisation of the field + 1 column for the normal field value, and 1 column for the tangential field value
  // if dim = 3 : le_champ_front has 1+4 dimensions  : 1 column for the localisation of the field + 1 column for the normal field value, and 3 column for the tangential field value

  if ((le_champ_front->valeurs().dimension(1)==2 && dimension == 2) || (le_champ_front->valeurs().dimension(1)==4 && dimension == 3))
    {
      if (le_champ_front->valeurs().dimension(0)==1) // a une valeur identique sur toutes les mailles du bord
        return le_champ_front->valeurs()(0,j);
      else
        return le_champ_front->valeurs()(i,j);
    }
  else
    Cerr << "Robin::flux_robin_normal_et_tangentiel error" << finl;
  Process::exit();
  return 0.;
}

/*! @brief Returns the value of the imposed normal flux on the i-th component of the field representing the flux at the boundary.
 *
 * @param (int i) index along the first dimension of the field
 * @return (double) the imposed value on the specified field component
 */

double Robin_VEF::flux_normal_imp(int i) const {	return flux_robin_normal_et_trangentiel_imp(i, 0);}

/*! @brief Returns the value of the imposed tangential flux on the (i,j)-th component of the field representing the tangential flux at the boundary.
 *
 * @param (int i) index along the first dimension of the field
 * @param (int j) index along the second dimension of the field
 * @return (double) the imposed value on the specified field component
 */
double Robin_VEF::flux_tangentiel_imp(int i, int j ) const
{
  if ((dimension==2 && j<2) || (dimension ==3 && j<4))
    return flux_robin_normal_et_trangentiel_imp(i, j+1);
  else
    Cerr << " Robin::flux_tangentiel_imp error, dimension does not match" << finl;
  Process::exit();
  return 0.;
}

void Robin_VEF::mettre_a_jour(double temps)
{
  le_champ_front->mettre_a_jour(temps);
}

// for domain decomposition, not implemented yet
double Robin_VEF::flux_robin_imp_au_temps(double temps, int i) const
{
  return 0.;

}

// for domain decomposition, not implemented yet
double Robin_VEF::flux_robin_imp_au_temps(double temps, int i, int j ) const
{
  return 0.;
}

/*! @brief Updates and returns the imposed flux array.
 *
 * This function checks the total number of faces in the boundary and resizes the
 * `flux_normal_impose_` array if necessary. It updates the values based on the imposed
 * flux from the front field if the field is unsteady or if the dimensions have changed.
 *
 * @return const DoubleTab& Reference to the updated imposed flux array.
 */
const DoubleTab& Robin_VEF::flux_robin_normal_imp() const
{
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());
  int nb_faces_tot = le_bord.nb_faces_tot();
  if (nb_faces_tot>0)
    {
      if (flux_normal_impose_.dimension(0) != nb_faces_tot)
        flux_normal_impose_.resize(nb_faces_tot, 1);
      int size = flux_normal_impose_.dimension(0);
      //int nb_comp = flux_robin_impose_.dimension(1);
      for (int i = 0; i < size; i++)
        //for (int j = 0; j < nb_comp; j++)
        flux_normal_impose_(i, 0) = flux_normal_imp(i);
    }
  return flux_normal_impose_;
}


/*! @brief Updates and returns the imposed flux array.
 *
 * This function checks the total number of faces in the boundary and resizes the
 * `flux_tangentiel_impose_` array if necessary. It updates the values based on the imposed
 * flux from the front field if the field is unsteady or if the dimensions have changed.
 *
 * @return const DoubleTab& Reference to the updated imposed flux array.
 */
const DoubleTab& Robin_VEF::flux_robin_tangentiel_imp() const
{
  const Front_VF& le_bord = ref_cast(Front_VF, frontiere_dis());
  int nb_faces_tot = le_bord.nb_faces_tot();
  if (nb_faces_tot>0)
    {
      if (flux_tangentiel_impose_.dimension(0) != nb_faces_tot)
        flux_tangentiel_impose_.resize(nb_faces_tot, le_champ_front->valeurs().dimension(1));
      int size = flux_tangentiel_impose_.dimension(0);
      int nb_comp = flux_tangentiel_impose_.dimension(1);
      for (int i = 0; i < size; i++)
        for (int j = 0; j < nb_comp; j++)
          flux_tangentiel_impose_(i, 0) = flux_tangentiel_imp(i,j);
    }
  return flux_tangentiel_impose_;
}

