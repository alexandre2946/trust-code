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

#include <Champ_front_contact_rayo_transp_VEF.h>
#include <Schema_Temps_base.h>
#include <Pb_Fluide_base.h>
#include <Champ_Inc_base.h>
#include <Pb_Conduction.h>
#include <Equation_base.h>
#include <Domaine_VEF.h>

Implemente_instanciable(Champ_front_contact_rayo_transp_VEF, "Champ_front_contact_rayo_transp_VEF", Champ_front_contact_VEF);

Sortie& Champ_front_contact_rayo_transp_VEF::printOn(Sortie& os) const { return os; }

Entree& Champ_front_contact_rayo_transp_VEF::readOn(Entree& is)
{
  fixer_nb_comp(1);
  return Champ_front_contact_VEF::readOn(is);
}

int Champ_front_contact_rayo_transp_VEF::initialiser(double temps, const Champ_Inc_base& inco)
{
  assert (!le_modele_rayo_);

  // retrieve the rayo model ... but we need the right problem !
  // XXX not yet entered Champ_front_contact_VEF::initialiser ... so some things must be done manually here ...
  const Probleme_base& pb1 = ref_cast(Probleme_base, Interprete::objet(nom_pb1));
  if (pb1.milieu().is_rayo_transp()) // good, we are on the fluid side !!
    {
      assert (sub_type(Pb_Fluide_base, pb1));
      le_modele_rayo_ =ref_cast(Pb_Fluide_base, pb1).get_mod_rayo_transp();
    }
  else // the other problem...
    {
      const Probleme_base& pb2 = ref_cast(Probleme_base, Interprete::objet(nom_pb2));
      if (pb2.milieu().is_rayo_transp()) // pb fluide trouve !!
        {
          assert (sub_type(Pb_Fluide_base, pb2));
          le_modele_rayo_ =ref_cast(Pb_Fluide_base, pb2).get_mod_rayo_transp();
        }
      else
        {
          Cerr << finl << "Big issue in Champ_front_contact_rayo_transp_VEF::initialiser !!!" << finl;
          Cerr << "It seems that you defined a radiation contact BC between the boundaries " << nom_bord1 << " and " << nom_bord2 << finl;
          Cerr << "of problems " << nom_pb1 << " and " << nom_pb2 << " , but neither is a fluid radiation problem !!!" << finl;
          Process::exit("Please fix your data file and use a classical paroi_contact BC for these boundaries ... \n");
        }
    }

  int nb_faces = frontiere_dis().frontiere().nb_faces();
  flux_radiatif_.resize(nb_faces);
  return Champ_front_contact_VEF::initialiser(temps, inco);
}

Champ_front_base& Champ_front_contact_rayo_transp_VEF::affecter_(const Champ_front_base& ch)
{
  return *this;
}

void Champ_front_contact_rayo_transp_VEF::mettre_a_jour_flux_radiatif()
{
  if (is_conduction)   // The model is known by the other problem
    {
      // Computation of the radiative flux on the boundary of the other problem
      DoubleTab flux_radiatif_autre_pb;
      int nb_faces = fr_vf_autre_pb->frontiere().nb_faces();
      int ndeb = fr_vf_autre_pb->frontiere().num_premiere_face();
      flux_radiatif_autre_pb.resize(nb_faces);
      for (int fac_front = 0; fac_front < nb_faces; fac_front++)
        flux_radiatif_autre_pb(fac_front) = le_modele_rayo_->flux_radiatif(fac_front + ndeb);
      // Fetch it from the other problem
      trace_face_raccord(fr_vf_autre_pb.valeur(), flux_radiatif_autre_pb, flux_radiatif_);
    }
  else
    {
      int nb_faces = frontiere_dis().frontiere().nb_faces();
      int ndeb = frontiere_dis().frontiere().num_premiere_face();
      flux_radiatif_.resize(nb_faces);
      for (int fac_front = 0; fac_front < nb_faces; fac_front++)
        flux_radiatif_(fac_front) = le_modele_rayo_->flux_radiatif(fac_front + ndeb);
    }
}

void Champ_front_contact_rayo_transp_VEF::mettre_a_jour(double temps)
{
  Champ_front_contact_VEF::mettre_a_jour(temps);
}

void Champ_front_contact_rayo_transp_VEF::calculer_coeffs_echange(double temps)
{
  Champ_front_contact_VEF::calculer_coeffs_echange(temps);

  // verify that the coupling is indeed radiative
  const Nom& nom_pb_rayonnant = le_modele_rayo_->nom_pb_rayonnant();

  // the radiating problem is neither of the 2 problems: the coupling is non-radiative
  if ((nom_pb_rayonnant != nom_pb1) && (nom_pb_rayonnant != nom_pb2))
    Process::exit("Error in Champ_front_contact_rayo_transp_VEF::calculer_coeffs_echange !! You should not be here since the problem is not a radiation problem !!!");

  // Computation and storage of the radiative flux
  mettre_a_jour_flux_radiatif();
  // Modify the gradients to account for the radiative flux
  DoubleVect gradient_num_transf_autre_pb(gradient_num_transf);
  modifie_gradients_pour_rayonnement(gradient_num_transf, gradient_num_transf_autre_pb);
}

void Champ_front_contact_rayo_transp_VEF::calculer_temperature_bord(double temps)
{
  // coding close to Champ_front_contact_VEF::mettre_a_jour
  const Frontiere& la_front = la_frontiere_dis->frontiere();
  int nb_faces = la_front.nb_faces();

  // Retrieve the coefficients gradient_num_transf and gradient_fro_transf from the other problem
  DoubleVect gradient_num_transf_autre_pb(nb_faces);
  DoubleVect gradient_fro_transf_autre_pb(nb_faces);
  trace_face_raccord(fr_vf_autre_pb.valeur(), ch_fr_autre_pb->get_gradient_num_transf(), gradient_num_transf_autre_pb);
  trace_face_raccord(fr_vf_autre_pb.valeur(), ch_fr_autre_pb->get_gradient_fro_transf(), gradient_fro_transf_autre_pb);

  // Computation and storage of the radiative flux

  // Retrieve the arrays used to compute omega, the damping factor
  DoubleVect coeff_amort_num_autre_pb(nb_faces);
  DoubleVect coeff_amort_denum_autre_pb(nb_faces);
  trace_face_raccord(fr_vf_autre_pb.valeur(), ch_fr_autre_pb->get_coeff_amort_num(), coeff_amort_num_autre_pb);
  trace_face_raccord(fr_vf_autre_pb.valeur(), ch_fr_autre_pb->get_coeff_amort_denum(), coeff_amort_denum_autre_pb);

  // Computation of the wall temperature
  DoubleTab& tab = valeurs_au_temps(temps);
  for (int fac_front = 0; fac_front < nb_faces; fac_front++)
    {
      // COMPUTATION OF THE DAMPING TERM
      Schema_Temps_base& sch = l_inconnue->equation().probleme().schema_temps();
      double dt = sch.pas_de_temps();
      double e = coeff_amort_num_autre_pb(fac_front) + coeff_amort_num(fac_front);
      double omega = dt / (dt + e / (coeff_amort_denum_autre_pb(fac_front) + coeff_amort_denum(fac_front)));
      // For now, omega is hard-coded to 1
      omega = 1;
      // END OF DAMPING TERM COMPUTATION

      double tab_past = tab(fac_front, 0);

      tab(fac_front, 0) = -(gradient_num_local(fac_front) + gradient_num_transf_autre_pb(fac_front)) / (gradient_fro_transf_autre_pb(fac_front) + gradient_fro_local(fac_front));

      tab(fac_front, 0) = omega * tab(fac_front, 0) + (1 - omega) * tab_past;
    }
}

void Champ_front_contact_rayo_transp_VEF::calcul_grads_locaux(double temps)
{
  Champ_front_contact_VEF::calcul_grads_locaux(temps);

  // Computation of the damping coefficients
  calcul_coeff_amort();
}

void Champ_front_contact_rayo_transp_VEF::modifie_gradients_pour_rayonnement(DoubleVect& tab_gradient_num_transf, DoubleVect& gradient_num_transf_autre_pb)
{
  const Frontiere& la_front = la_frontiere_dis->frontiere();
  const int nb_faces = la_front.nb_faces();

  if (is_conduction)
    {
      for (int fac_front = 0; fac_front < nb_faces; fac_front++)
        gradient_num_local(fac_front) = gradient_num_local(fac_front) + flux_radiatif_(fac_front);
    }
  else
    {
      for (int fac_front = 0; fac_front < nb_faces; fac_front++)
        gradient_num_transf_autre_pb(fac_front) = gradient_num_transf_autre_pb(fac_front) - flux_radiatif_(fac_front);
    }
}
