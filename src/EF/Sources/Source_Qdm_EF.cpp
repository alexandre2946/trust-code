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

#include <Domaine_Cl_dis_base.h>
#include <Champ_Uniforme.h>

#include <Source_Qdm_EF.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Domaine_Cl_EF.h>
#include <Domaine_EF.h>
#include <Domaine.h>

Implemente_instanciable(Source_Qdm_EF,"Source_Qdm_EF",Source_base);

Sortie& Source_Qdm_EF::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

Entree& Source_Qdm_EF::readOn(Entree& s )
{
  Motcle type;
  s >> type;
  la_source_lu.typer(type);
  Champ_Don_base& ch_source_lu = ref_cast(Champ_Don_base,la_source_lu.valeur());
  s >> ch_source_lu;
  const int nb_comp = ch_source_lu.nb_comp();

  if (ch_source_lu.nb_comp() != dimension)
    {
      Cerr << "Error reading the source term of type " << que_suis_je() << finl;
      Cerr << "the source field must have " << dimension << " components" << finl;
      exit();
    }

  equation().probleme().discretisation().discretiser_champ("champ_elem", equation().domaine_dis(), "Source_Qdm", "N/m3",nb_comp,0., la_source);
  la_source_lu->fixer_nb_comp(nb_comp);
  if (ch_source_lu.le_nom()=="anonyme") ch_source_lu.nommer("Source_Qdm");

  for (int n = 0; n < nb_comp; n++) la_source_lu->fixer_nom_compo(n, ch_source_lu.le_nom() + (nb_comp > 1 ? Nom(n) :""));
  for (int n = 0; n < nb_comp; n++) la_source->fixer_nom_compo(n, ch_source_lu.le_nom() + (nb_comp > 1 ? Nom(n) :""));
  // PL: nommer_completer_champ_physique must be called for both fields (crash otherwise for a Champ_fonc_tabule source type)
  equation().discretisation().nommer_completer_champ_physique(equation().domaine_dis(),ch_source_lu.le_nom(),"N/m3",la_source_lu,equation().probleme());
  equation().discretisation().nommer_completer_champ_physique(equation().domaine_dis(),ch_source_lu.le_nom(),"N/m3",la_source,equation().probleme());
  la_source->valeurs() = 0.;
  la_source->affecter(ch_source_lu);

  return s ;
}

void Source_Qdm_EF::associer_pb(const Probleme_base& )
{
  ;
}

void Source_Qdm_EF::associer_domaines(const Domaine_dis_base& domaine_dis,
                                      const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_EF = ref_cast(Domaine_EF, domaine_dis);
  le_dom_Cl_EF = ref_cast(Domaine_Cl_EF, domaine_Cl_dis);
}


DoubleTab& Source_Qdm_EF::ajouter(DoubleTab& resu) const
{
  const Domaine_EF& domaine_EF = le_dom_EF.valeur();
  int ncomp=dimension;
  const IntTab& elems= domaine_EF.domaine().les_elems() ;
  int nb_som_elem=domaine_EF.domaine().nb_som_elem();
  int nb_elems=domaine_EF.domaine().nb_elem_tot();
  const DoubleTab& IPhi_thilde=domaine_EF.IPhi_thilde();
  int is_source_unif=0;

  if (sub_type(Champ_Uniforme,la_source.valeur()))
    is_source_unif=1;
  const DoubleTab& tab_source=la_source->valeurs();
  for (int num_elem=0; num_elem<nb_elems; num_elem++)
    for (int comp=0; comp<ncomp; comp++)
      {
        double sourcel ;
        int cc=0;
        if (!is_source_unif) cc=num_elem;
        sourcel=tab_source(cc,comp);
        for (int i=0; i<nb_som_elem; i++)
          {
            // assert is wrong if porosities are present
            //	assert(est_egal(domaine_EF.volumes(num_elem)/nb_som_elem,IPhi_thilde(num_elem,i)));
            resu(elems(num_elem,i),comp)+=sourcel*IPhi_thilde(num_elem,i);
          }
      }

  return resu;
}

DoubleTab& Source_Qdm_EF::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

void Source_Qdm_EF::mettre_a_jour(double temps)
{
  la_source_lu->mettre_a_jour(temps);
  la_source->affecter(la_source_lu.valeur());
}


/*##################################################################################################
####################################################################################################
################################# POSTRAITEMENT ####################################################
####################################################################################################
##################################################################################################*/

void Source_Qdm_EF::creer_champ(const Motcle& motlu)
{
  if (motlu=="source_qdm" && !champ_source_qdm_)
    {
      Motcle directive("champ_elem");
      int nb_comp = dimension;
      Noms noms(nb_comp);
      noms[0]="source_qdm";
      Noms unites(nb_comp);
      unites[0] = "N/m3";
      double temps=0.;
      equation().discretisation().discretiser_champ(directive,equation().domaine_dis(),vectoriel,noms,unites,nb_comp,temps,champ_source_qdm_);
      champs_compris_.ajoute_champ(champ_source_qdm_);
    }
}

bool Source_Qdm_EF::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  if (nom == "source_qdm" && champ_source_qdm_)
    {
      ref_champ = Source_Qdm_EF::get_champ(nom);
      return true;
    }
  else
    return false; /* nothing found */
}

bool Source_Qdm_EF::has_champ(const Motcle& nom) const
{
  if (Source_base::has_champ(nom)) return true;

  if (nom == "source_qdm" && champ_source_qdm_)
    return true;
  else
    return false; /* nothing found */
}

const Champ_base& Source_Qdm_EF::get_champ(const Motcle& nom) const
{
  if (nom=="source_qdm")
    {
      if (!champ_source_qdm_)
        throw std::runtime_error(std::string("Field ") + nom.getString() + std::string(" not found !"));

      int is_source_unif=0;
      if (sub_type(Champ_Uniforme,la_source.valeur()))
        is_source_unif=1;

      // Initialize the post-processing field to 0
      DoubleTab& valeurs = champ_source_qdm_->valeurs();
      valeurs=0.;
      const DoubleTab& tab_source=la_source->valeurs();
      if (tab_source.size_array()>0)
        {
          int nb_elem=tab_source.dimension(0);
          for (int num_el=0; num_el<nb_elem; num_el++)
            for (int k=0; k<dimension; k++)
              {
                int cc=0;
                if (!is_source_unif) cc=num_el;
                valeurs(num_el,k)=tab_source(cc,k);
              }
        }
      valeurs.echange_espace_virtuel();
      champ_source_qdm_->mettre_a_jour(equation().probleme().schema_temps().temps_courant());
      return champs_compris_.get_champ(nom);
    }
  else
    return champs_compris_.get_champ(nom);
}

void Source_Qdm_EF::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Source_base::get_noms_champs_postraitables(nom,opt);

  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("source_qdm");
  if (opt==DESCRIPTION)
    Cerr<<" Source_Qdm_EF : "<< noms_compris <<finl;
  else
    nom.add(noms_compris);
}
