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

#include <Terme_Derivee_Forme_base.h>
#include <Probleme_base.h>
#include <Milieu_base.h>
#include <Equation_base.h>
#include <Champ_Uniforme.h>
#include <Champ_Fonc_Tabule.h>
#include <Discretisation_base.h>
#include <Champ_val_tot_sur_vol_base.h>

Implemente_base(Terme_Derivee_Forme_base,"Terme_Derivee_Forme_base",Source_base);
// XD  derivee_forme  source_base derivee_forme 0 Class to define a source term corresponding to the shape derivative in the projection equation.

Entree& Terme_Derivee_Forme_base::readOn(Entree& s)
{
  Cerr << "Lecture Derivee de forme" << finl;

  equation().discretisation().discretiser_champ(equation().discretisation().is_poly_family() ? "temperature" : "champ_elem", equation().domaine_dis(), "derivee_forme", ",",1,0., source_derivee_forme);
  source_derivee_forme->nommer("derivee_forme");
  source_derivee_forme->valeurs() = 0.;
  champs_compris_.ajoute_champ(source_derivee_forme);
  return s;
}

Sortie& Terme_Derivee_Forme_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

void Terme_Derivee_Forme_base::modify_name_file(Nom& fichier) const
{
  if (source_derivee_forme->le_nom()!="derivee_forme")
    {
      fichier+="_";
      fichier+=source_derivee_forme->le_nom();
    }
}

void Terme_Derivee_Forme_base::set_source_derivee_forme(DoubleTab& chp) const
{
  DoubleTab& chpval = ref_cast_non_const(DoubleTab, source_derivee_forme->valeurs());
  assert(chpval.dimension(0) == chp.dimension(0));
  assert(chpval.dimension(1) == chp.dimension(1));
  for (int i=0; i<chpval.dimension(0); i++)
    {
      for (int j=0; j<chpval.dimension(1); j++)
        {
          chpval(i,j) = chp(i,j);
        }
    }
}

void Terme_Derivee_Forme_base::mettre_a_jour(double temps)
{
  source_derivee_forme->mettre_a_jour(temps);
}

DoubleTab& Terme_Derivee_Forme_base::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

void Terme_Derivee_Forme_base::creer_champ(const Motcle& motlu)
{
  const Domaine_dis_base& le_dom_dis = equation().domaine_dis();
  const Probleme_base& pb = equation().probleme();
  if (motlu=="source_from_shape_deriv" && !champ_derivee_forme_.non_nul())
    {
      Noms noms(1);
      noms[0]="source_from_shape_deriv";
      Noms unites(1);
      unites[0] = "-";
      pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,scalaire,noms,unites,1,0.,champ_derivee_forme_);
      champs_compris_.ajoute_champ(champ_derivee_forme_);
    }
}

void Terme_Derivee_Forme_base::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("source_from_shape_deriv");
  if (opt==DESCRIPTION)
    Cerr<<" Terme_Derivee_Forme_base : "<< noms_compris <<finl;
  else
    nom.add(noms_compris);
}

bool Terme_Derivee_Forme_base::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  if (nom == "source_from_shape_deriv" && champ_derivee_forme_.non_nul())
    {
      ref_champ = get_champ(nom);
      return true;
    }
  else if (champs_compris_.has_champ(nom, ref_champ))
    return true;
  else
    return false;
}

bool Terme_Derivee_Forme_base::has_champ(const Motcle& nom) const
{
  if (nom == "source_from_shape_deriv" && champ_derivee_forme_.non_nul())
    return true;
  else
    return champs_compris_.has_champ(nom);
}

const Champ_base& Terme_Derivee_Forme_base::get_champ(const Motcle& nom) const
{
  if (nom=="source_from_shape_deriv")
    {
      if (champ_derivee_forme_.est_nul())
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      // Initialisation a 0 du champ_derivee_forme_
      DoubleTab& valeurs = champ_derivee_forme_->valeurs();
      valeurs=0.;
      const DoubleTab& derivee_forme_Array = source_derivee_forme->valeurs();
      if (derivee_forme_Array.size_array()>0)
        {
          int nb_d0 = derivee_forme_Array.dimension(0);
          int nb_d1 = derivee_forme_Array.dimension(1);
          for (int num0=0; num0<nb_d0; num0++)
            for (int num1=0; num1<nb_d1; num1++)
              valeurs(num0,num1)=derivee_forme_Array(num0,num1);
        }
      valeurs.echange_espace_virtuel();
      champ_derivee_forme_ ->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else
    return champs_compris_.get_champ(nom);
}
