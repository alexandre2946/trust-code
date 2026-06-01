/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Terme_Derivee_Forme_EF.h>

Implemente_instanciable(Terme_Derivee_Forme_EF,"Derivee_Forme_EF",Terme_Derivee_Forme_base);

// ################################# READON AND PRINTON ###############################################
Entree& Terme_Derivee_Forme_EF::readOn(Entree& s)
{
  Terme_Derivee_Forme_base::readOn(s);
  return s;
}


Sortie& Terme_Derivee_Forme_EF::printOn(Sortie& s ) const
{
  Terme_Derivee_Forme_base::printOn(s);
  return s ;
}
// ##################################################################################################

void Terme_Derivee_Forme_EF::associer_domaines(const Domaine_dis_base& domaine_dis,
                                               const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  le_dom_EF = ref_cast(Domaine_EF, domaine_dis);
}

void Terme_Derivee_Forme_EF::associer_pb(const Probleme_base& pb ) { }

DoubleTab& Terme_Derivee_Forme_EF::ajouter(DoubleTab& resu) const
{
  const Domaine_EF& domaine_EF = le_dom_EF.valeur();
  int ncomp=equation().inconnue().nb_comp();
  const IntTab& elems= domaine_EF.domaine().les_elems() ;
  int nb_som_elem=domaine_EF.domaine().nb_som_elem();
  int nb_elems=domaine_EF.domaine().nb_elem_tot();

  // Note: computation is IPhi *S(e)
  //  const DoubleTab& IPhi_thilde=domaine_EF.IPhi_thilde();

  const DoubleTab& IPhi=domaine_EF.IPhi();

  int is_source_unif=0;

  if (sub_type(Champ_Uniforme,source_derivee_forme.valeur()))
    is_source_unif=1;
  const DoubleTab& tab_source=source_derivee_forme->valeurs();
  if (ncomp>1)
    for (int num_elem=0; num_elem<nb_elems; num_elem++)
      for (int comp=0; comp<ncomp; comp++)
        {
          double source ;
          int cc=0;
          if (!is_source_unif) cc=num_elem;
          source=tab_source(cc,comp);
          for (int i=0; i<nb_som_elem; i++)
            resu(elems(num_elem,i),comp)+=source*IPhi(num_elem,i);
        }
  else
    {
      for (int num_elem=0; num_elem<nb_elems; num_elem++)
        {
          double source ;

          if (!is_source_unif)
            source=tab_source(num_elem);
          else
            source=tab_source(0,0);
          for (int i=0; i<nb_som_elem; i++)
            resu(elems(num_elem,i))+=source*IPhi(num_elem,i);

        }

    }

  return resu;
}

DoubleTab& Terme_Derivee_Forme_EF::calculer(DoubleTab& resu) const
{
  return Terme_Derivee_Forme_base::calculer(resu);
}

void Terme_Derivee_Forme_EF::mettre_a_jour(double temps)
{
  Terme_Derivee_Forme_base::mettre_a_jour(temps);
}
