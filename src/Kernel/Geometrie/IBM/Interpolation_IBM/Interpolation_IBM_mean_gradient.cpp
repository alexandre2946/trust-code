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

#include <Interpolation_IBM_mean_gradient.h>
#include <TRUSTTrav.h>

#include <Domaine.h>
#include <Param.h>
#include <Source_PDF_base.h>

Implemente_instanciable( Interpolation_IBM_mean_gradient, "Interpolation_IBM_gradient_moyen|IBM_gradient_moyen", Interpolation_IBM_base ) ;
// XD interpolation_ibm_mean_gradient interpolation_ibm_base ibm_gradient_moyen 1 Immersed Boundary Method (IBM): mean gradient interpolation.

Sortie& Interpolation_IBM_mean_gradient::printOn( Sortie& os ) const
{
  return Objet_U::printOn( os );
}

Entree& Interpolation_IBM_mean_gradient::readOn( Entree& is )
{
  Param param(que_suis_je());
  param.ajouter("points_solides",&solid_points_lu_,Param::OPTIONAL);  // XD_ADD_P field_base Node field giving the projection of the node on the immersed boundary
  param.ajouter("est_dirichlet",&is_dirichlet_lu_,Param::OPTIONAL);   // XD_ADD_P field_base Node field of booleans indicating whether the node belong to an element where the interface is
  param.ajouter("correspondance_elements",&corresp_elems_lu_,Param::OPTIONAL); // XD_ADD_P field_base Cell field giving the SALOME cell number
  param.ajouter("elements_solides",&solid_elems_lu_,Param::OPTIONAL); // XD_ADD_P field_base Node field giving the element number containing the solid point
  param.ajouter_flag("get_solid_points_from_prepro", &solid_points_from_prepro_); // XD_ADD_P rien get IBM solid points from prepro.
  param.ajouter_flag("get_solid_elems_from_prepro", &solid_elems_from_prepro_); // XD_ADD_P rien get IBM solid elems from prepro.
  param.ajouter_flag("get_is_dirichlet_from_prepro", &is_dirichlet_from_prepro_); // XD_ADD_P rien get IBM is_dirichlet from prepro.
  param.ajouter_flag("get_corresp_elems_from_prepro", &corresp_elems_from_prepro_); // XD_ADD_P rien get IBM corresp_elems from prepro.
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Interpolation_IBM_mean_gradient::discretise(const Discretisation_base& dis, Domaine_dis_base& le_dom_)
{
  Interpolation_IBM_base::discretise(dis, le_dom_);
  int nb_comp = Objet_U::dimension;
  Noms units(nb_comp);
  Noms c_nam(nb_comp);

  if (corresp_elems_.non_nul()) has_corresp_ = true;

  dis.discretiser_champ("champ_sommets",le_dom_,"solid_elems","none",1,0., solid_elems_);
  if (solid_elems_from_prepro_)
    {
      OBS_PTR(Prepro_IBM_base) my_prep =  my_source_->getpreproLu();
      if ((&my_prep)->non_nul())
        {
          DoubleTab& the_values = ref_cast_non_const(DoubleTab, my_prep->get_champ_solid_elems());
          solid_elems_->valeurs() = the_values;
        }
    }
  else
    {
      if (solid_elems_lu_.non_nul()) solid_elems_->affecter(solid_elems_lu_);
    }

  if(is_dirichlet_.est_nul())
    {
      Cerr<<"Interpolation_IBM_mean_gradient: field est_dirichlet is required. exit()"<<endl;
      Process::exit();
    }
  else
    my_is_dirichlet_ = is_dirichlet_;
  computeSommetsVoisins(le_dom_, solid_points_, corresp_elems_, has_corresp_);
}

void Interpolation_IBM_mean_gradient::set_fields_from_prepro_to_interp(Prepro_IBM_base& un_prepro)
{
  Interpolation_IBM_base::set_fields_from_prepro_to_interp(un_prepro);
  solid_elems_->valeurs() = un_prepro.get_champ_solid_elems();
}
