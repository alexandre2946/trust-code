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

#include <Modele_turbulence_hyd_2_eq_base.h>
#include <TRUST_2_PDI.h>
#include <TRUSTTrav.h>
#include <Param.h>

Implemente_base(Modele_turbulence_hyd_2_eq_base, "Modele_turbulence_hyd_2_eq_base", Modele_turbulence_hyd_base);

Sortie& Modele_turbulence_hyd_2_eq_base::printOn(Sortie& is) const { return Modele_turbulence_hyd_base::printOn(is); }
Entree& Modele_turbulence_hyd_2_eq_base::readOn(Entree& is) { return Modele_turbulence_hyd_base::readOn(is); }

void Modele_turbulence_hyd_2_eq_base::set_param(Param& param) const
{
  Modele_turbulence_hyd_base::set_param(param);
  param.ajouter_non_std("Transport_equation", (this)); // cannot be REQUIRED because of Bicephale models
  param.ajouter_non_std("Transport_K_Epsilon", (this));
  param.ajouter_non_std("Transport_K_Omega", (this));
  param.ajouter_non_std("Transport_K_Epsilon_Realisable", (this));
  param.ajouter_non_std("transport_k", (this));
  param.ajouter_non_std("transport_epsilon", (this));
  param.ajouter("k_min", &K_MIN_);
  param.ajouter_flag("quiet", &lquiet_);
}

int Modele_turbulence_hyd_2_eq_base::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "transport_equation")
    {
      Motcle typ_eq;
      is>> typ_eq;
      Nom name_transport_eq=typ_eq;
      ptr_eq_transport_.typer(name_transport_eq);
      get_set_eq_transport().associer_modele_turbulence(*this);
      is >> get_set_eq_transport();
    }
  else if (mot == "Transport_K_Epsilon" || mot == "Transport_K_Omega" || mot == "Transport_K_Epsilon_Realisable")
    {
      Cerr << "Error: You are using an obsolete syntaxe for " << mot << " in your datafile!!!!!!!!!!!!!!!!" << finl;
      Cerr << "         Since v1.9.7, you should replace:" << finl;
      Cerr << "             " << mot << finl;
      Cerr << "         by:" << finl;
      Cerr << "             Transport_equation " << mot << finl;
      Cerr << "Please update your datafile" << finl;
      Process::exit();
    }
  else if (mot == "transport_k" || mot == "transport_epsilon")
    {
      Cerr << "Error: You are using an obsolete syntaxe for " << mot << " in your datafile!!!!!!!!!!!!!!!!" << finl;
      Cerr << "       Since v1.9.7, you should use a block like:" << finl;
      Cerr << "            List_transport_equations {" << finl;
      Cerr << "                                      K_equation transport_k_ou_eps { .... }" << finl;
      Cerr << "                                      Eps_equation transport_k_ou_eps { .... }" << finl;
      Cerr << "                                     }" << finl;
      Cerr << "Please update your datafile" << finl;
      Process::exit();
    }
  else
    return Modele_turbulence_hyd_base::lire_motcle_non_standard(mot, is);
  return 1;
}

void Modele_turbulence_hyd_2_eq_base::verifie_loi_paroi()
{
  Nom lp = loipar_->que_suis_je();
  if (lp == "negligeable_VEF" || lp == "negligeable_VDF")
    {
      Cerr << "The turbulence model of type " << que_suis_je() << finl;
      Cerr << "must not be used with a wall law of type negligeable." << finl;
      Cerr << "Another wall law must be selected with this kind of turbulence model." << finl;
      Process::exit();
    }
}

int Modele_turbulence_hyd_2_eq_base::reprendre_generique(Entree& is)
{
  if(TRUST_2_PDI::is_PDI_restart())
    {
      Cerr <<"Problem in the resumption of Modele_turbulence_hyd_2_eq_base" << finl;
      Cerr << "PDI format does not require to navigate through file..." << finl;
      Process::exit();
    }

  double dbidon;
  Nom bidon;
  DoubleTrav tab_bidon;
  is >> bidon >> bidon;
  is >> dbidon;
  tab_bidon.jump(is);
  return 1;
}

const Transport_2eq_base& Modele_turbulence_hyd_2_eq_base::get_eq_transport() const
{
  const Transport_2eq_base& eq_transport = ptr_eq_transport_.valeur();
  return eq_transport;
}

Champ_Inc_base& Modele_turbulence_hyd_2_eq_base::get_set_unknown()
{
  return get_set_eq_transport().inconnue();
}

const Champ_Inc_base& Modele_turbulence_hyd_2_eq_base::get_unknown() const
{
  const Champ_Inc_base& unkwown = get_eq_transport().inconnue();
  return unkwown;
}
