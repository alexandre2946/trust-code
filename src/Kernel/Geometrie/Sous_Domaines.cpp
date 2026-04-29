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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Sous_Domaines.cpp
// Directory:   $TRUST_ROOT/src/Kernel/Utilitaires
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#include <Sous_Domaines.h>
#include <Param.h>
#include <Domaine.h>
#include <EChaine.h>
#include <Synonyme_info.h>
#include <Interprete.h>

Implemente_instanciable(Sous_Domaines,"Sous_Domaines",VECT(Sous_Domaine));
Add_synonym(Sous_Domaines, "Sous_Zones");

Sortie& Sous_Domaines::printOn(Sortie& s) const
{
  return VECT(Sous_Domaine)::printOn(s);
}

// Syntaxe identique au decoupeur tranches mais en specifiant un domaine pour la bounding box: { bounding_box domaine_name tranches Nx Ny [Nz] }
Entree& Sous_Domaines::readOn(Entree& s)
{
  Motcle ouvert,ferme,bounding_box,tranches;
  ArrOfInt n(dimension);
  n=1;
  Nom nom_domaine;
  s >> ouvert >> bounding_box >> nom_domaine >> tranches;
  for (int dir=0; dir<dimension; dir++) s >> n[dir];
  s >> ferme;
  if (ouvert!="{" || bounding_box!="BOUNDING_BOX" || tranches!="TRANCHES" || ferme!="}")
    Process::exit("Syntax is: Lire|Read Sous_Domaines { bounding_box domaine_name tranches Nx Ny [Nz] }");
  DoubleTab BB = ref_cast(Domaine, Interprete::objet(nom_domaine)).getBoundingBox();
  ArrOfDouble cotes(dimension);
  double alpha = 1.0001; // Astuce pour pouvoir decouper par faisceau ou demi-faisceau si le canal de croix est maille (mais alors decoupage non equilibre)
  for (int dir=0; dir<dimension; dir++)
    {
      BB(dir, 0) = mp_min(BB(dir, 0));
      BB(dir, 1) = mp_max(BB(dir, 1));
      cotes(dir) = (BB(dir, 1) - BB(dir,0)) / n[dir] * alpha;
    }
  Domaine& dom = ref_cast(Sous_Domaine, front()).domaine();
  dimensionner_force(0);
  int canal=0;
  int nb_elem_ssz = 0;
  for (int i=0; i<n[0]; i++)
    for (int j=0; j<n[1]; j++)
      for (int k=0; k<n[2]; k++)
        {
          // Creation d'une Sous_Domaine de type Rectangle/Boite
          Nom tmp("{ ");
          tmp+=(dimension==2 ? "Rectangle " : "Boite ");
          tmp+=" Origine";
          for (int dir=0; dir<dimension; dir++)
            {
              tmp += " ";
              tmp += (Nom)(BB(dir, 0) + (dir==0 ? i : (dir==1 ? j : k))*cotes(dir));
            }
          tmp+=" Cotes";
          for (int dir=0; dir<dimension; dir++)
            {
              tmp += " ";
              tmp += (Nom)cotes(dir);
            }
          tmp+=" }";
          //Cerr << "Echaine: " << tmp << finl;
          Sous_Domaine& ssz = add(Sous_Domaine());
          ssz.associer_domaine(dom);
          Nom name(dom.le_nom());
          name+="_Sous_Domaine_";
          name+=(Nom)canal;
          ssz.nommer(name);
          EChaine is(tmp);
          is >> ssz;
          if (ssz.nb_elem_tot()!=0)
            {
              dom.add(ssz);
              nb_elem_ssz+=ssz.nb_elem_tot();
              canal++;
            }
        }
  // Add check:
  if (nb_elem_ssz!=dom.nb_elem_tot())
    {
      Cerr << "Error, all the domain is not covered by subdomaines:" << finl;
      Cerr << "There is " << dom.nb_elem_tot() - nb_elem_ssz << " cells not marked." << finl;
      Process::exit();
    }
  return s;
}

int Sous_Domaines::associer_(Objet_U& ob)
{
  if (sub_type(Domaine, ob))
    {
      // Creation d'une premiere sous_domaine pour stocker reference vers domaine/domaine
      Sous_Domaine& ssz = add(Sous_Domaine());
      ssz.associer_domaine(ref_cast(Domaine, ob));
      return 1;
    }
  else
    Process::exit("You need to associate a Domaine to Sous_Domaines.");
  return 0;
}
