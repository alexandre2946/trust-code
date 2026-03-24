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

#include <LecFicDiffuse.h>
#include <VerifierCoin.h>
#include <Scatter.h>
#include <Domaine.h>
#include <Param.h>
#include <time.h>

Implemente_instanciable(VerifierCoin,"VerifierCoin",Interprete_geometrique_base);

/*! @brief Simple appel a: Interprete::printOn(Sortie&) */
Sortie& VerifierCoin::printOn(Sortie& os) const { return Interprete::printOn(os); }

/*! @brief Simple appel a: Interprete::readOn(Entree&) */
Entree& VerifierCoin::readOn(Entree& is) { return Interprete::readOn(is); }

// Split element `elem` into (dimension+1) sub-elements by inserting a centroid node.
// The centroid coordinates are taken from xp(elem,:).
// New elements are appended to les_elems; the sub-domains index is updated accordingly.
void VerifierCoin::cut_elem(int elem, const DoubleTab& xp)
{
  Domaine& dom=domaine();
  DoubleTab& sommets = dom.les_sommets();
  int nouveau_sommet = sommets.dimension(0);
  IntTab& les_elems=dom.les_elems();
  sommets.resize(nouveau_sommet+1, dimension);
  for (int j = 0; j < dimension; j++)
    sommets(nouveau_sommet,j) = xp(elem,j);

  int oldsz = les_elems.dimension(0);
  les_elems.resize(oldsz+dimension, dimension+1);
  Cerr << "-> The element number " << elem << " is cut in " << dimension+1 << " elements." << finl;

  int i0 = les_elems(elem,0);
  int i1 = les_elems(elem,1);
  int i2 = les_elems(elem,2);

  les_elems(elem,0) = i0;
  les_elems(elem,1) = i1;
  les_elems(elem,2) = nouveau_sommet;

  les_elems(oldsz,0) = i1;
  les_elems(oldsz,1) = i2;
  les_elems(oldsz,2) = nouveau_sommet;

  les_elems(oldsz+1,0) = i0;
  les_elems(oldsz+1,1) = i2;
  les_elems(oldsz+1,2) = nouveau_sommet;

  if (dimension == 3)
    {
      int i3 = les_elems(elem,3);

      les_elems(elem,3)   = i3;
      les_elems(oldsz,3)  = i3;
      les_elems(oldsz+1,3)= i3;

      les_elems(oldsz+2,0) = i0;
      les_elems(oldsz+2,1) = i1;
      les_elems(oldsz+2,2) = i2;
      les_elems(oldsz+2,3) = nouveau_sommet;
    }

  mettre_a_jour_sous_domaine(dom, elem, oldsz, dimension);
}

/*! @brief Fonction principale de l'interprete: resoudre un probleme
 *
 *     On cherche dynamiquement le type du probleme a resoudre
 *     on resoud le probleme et on effectue les postraitements.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 * @throws type de probleme inconnu
 */
Entree& VerifierCoin::interpreter_(Entree& is)
{
// XD verifiercoin interprete verifiercoin -1 This keyword subdivides inconsistent 2D/3D cells used with VEFPreP1B discretization. Must be used before the mesh is discretized. NL1 The Read_file option can be used only if the file.decoupage_som was previously created by TRUST. This option, only in 2D, reverses the common face at two cells (at least one is inconsistent), through the nodes opposed. In 3D, the option has no effect.
// XD  attr domain_name  ref_domaine dom 0 Name of the domaine
// XD  attr bloc verifiercoin_bloc bloc 0 not_set
// XD  verifiercoin_bloc objet_lecture nul 1 not_set
// XD  attr Read_file|Lire_fichier chaine filename 1 name of the *.decoupage_som file


  associer_domaine(is);
  Domaine& dom=domaine();
  if (dom.type_elem()->que_suis_je() != "Triangle" && dom.type_elem()->que_suis_je() != "Tetraedre")
    {
      Cerr << "Error for "<<que_suis_je() <<" interpreter : it can be applied only for triangular or tetraedral meshing." << finl;
      exit();
    }

  Nom decoup_som("");
  Param param(que_suis_je());
  param.ajouter("Lire_fichier|Read_file",&decoup_som);
  bool expert_only_obsolete = false;
  param.ajouter_flag("expert_only",&expert_only_obsolete);
  param.lire_avec_accolades_depuis(is);
  if (expert_only_obsolete) Process::exit("Error: expert_only is not supported anymore for VerifierCoin");
  int lecture_decoupage_som=(decoup_som!=""?1:0);

  if (Process::is_parallel())
    {
      Cerr << que_suis_je() << " interpreter can be used only for sequential calculation." << finl;
      Cerr << "For parallel calculation, please use "<< que_suis_je() <<" interpreter" << finl;
      Cerr << "during the domain partitioning step." << finl;
      exit();
    }
  Scatter::uninit_sequential_domain(dom);
  DoubleTab xp;
  dom.calculer_centres_gravite(xp);
  DoubleTab& sommets = dom.les_sommets();
  int nbsom=sommets.dimension(0);

  IntTab& les_elems=dom.les_elems();
  int nbelem=dom.nb_elem();

  // On compte les elements attaches a chaque sommet:
  ArrOfInt nb_elem_per_som(nbsom);
  nb_elem_per_som = 0;
  for (int ne = 0; ne < nbelem; ne++)
    for (int ns = 0; ns < dimension+1; ns++)
      nb_elem_per_som(les_elems(ne,ns))++;

  //On decoupe les elements pour le sommet qui pose probleme
  // PQ : 25/05/07
  // - soit de maniere automatique (tout sommet rattache qu'a un seul element)
  // - soit a partir d'une liste de sommets (liste generee lors d'un pre-calcul et tenant compte du type de CL)
  //   avec comme option de decoupage dans ce cas (lue dans le fichier "decoupage_som") :
  //
  // 0 : decoupage traditionnel (centre de gravite)
  // 1 : decoupage en passant par le sommet oppose de l'element voisin

  int option_decoupage=-1;

  LecFicDiffuse fic;
  if (lecture_decoupage_som)
    {
      Cerr<<"The file" << decoup_som <<" is checked before reading."<<finl;
      fic.ouvrir(decoup_som);
      if(fic.good())
        {
          int dim_cas, nbsom_cas;
          fic >> option_decoupage;
          fic >> dim_cas;
          fic >> nbsom_cas;
          if(dim_cas!=dimension)
            {
              Cerr << "Error for " <<que_suis_je()<<"::interpreter_" << finl;
              Cerr << "The file " << decoup_som << " is planned for a case of dimension " << dim_cas << "D" << finl;
              Cerr << "while the considered domain " << dom.le_nom() << " has a dimension " << dimension << "D." << finl;
              Process::exit();
            }
          if(nbsom!=nbsom_cas)
            {
              Cerr << "Error for " <<que_suis_je()<<"::interpreter_" << finl;
              Cerr << "The nodes number (" <<nbsom_cas<< ") readen in the file " << decoup_som << finl;
              Cerr << "dot not correspond to those of the meshing (" << nbsom << ")" << finl;
              Cerr << "Please check that the file you have specified " << decoup_som << finl;
              Cerr << "is the one to consider for this meshing." << finl;
              Process::exit();
            }
        }
      else
        {
          Cerr << "Error for " <<que_suis_je()<<"::interpreter_" << finl;
          Cerr << "Problem while trying to open the file " << decoup_som << finl;
          Process::exit();
        }
      Cerr << "option_decoupage " << option_decoupage << finl;

      for (int somm=0; somm<nbsom; somm++)
        {
          if (nb_elem_per_som(somm) != 1) continue;
          //On decoupe les elements pour le sommet qui pose probleme
          int somm_lu = -1, elem_opp, somm_opp, somm1, somm2, elem;
          fic >> somm_lu >> somm_opp >> somm1 >> somm2;
          if (dimension==3)
            {
              int somm3;
              fic >> somm3;
            }
          fic >> elem >> elem_opp;
          if(somm_lu!=-1) somm = somm_lu; // -1 indice de fin de fichier
          if(nbsom<=somm)
            {
              Cerr << "Error in VerifierCoin::interpreter" << finl;
              Cerr << "The node " << somm << " is not found." << finl;
              Cerr << "Check the .Zones files are up to date with your mesh file." << finl;
              Process::exit();
            }
          Cerr<<"-> VerifierCoin is applied on the node "<<somm<< " of coordinates: ";
          for(int dir=0; dir<dimension; dir++) Cerr<<sommets(somm,dir)<<" ";
          Cerr<<"..."<<finl;

          if (option_decoupage==1 && dimension==2 && somm_lu!=-1)  // inversion des sommets
            {
              les_elems(elem,0) = somm ;
              les_elems(elem,1) = somm_opp ;
              les_elems(elem,2) = somm1 ;

              les_elems(elem_opp,0) = somm ;
              les_elems(elem_opp,1) = somm_opp ;
              les_elems(elem_opp,2) = somm2 ;
            }
          else // creation d'un nouveau sommet au centre de gravite de l'element
            cut_elem(elem, xp);
        }
    }
  else
    {
      std::map<int,int> som_elem; // Pour trier les sommets comme avant
      for (int elem = 0; elem < nbelem; elem++)
        for (int ns = 0; ns < dimension+1; ns++)
          {
            int somm = les_elems(elem, ns);
            if (nb_elem_per_som(somm) == 1) som_elem.insert({somm, elem});
          }
      for (auto pair : som_elem)
        {
          int somm = pair.first;
          int elem = pair.second;
          nb_elem_per_som(somm) = 0;  // mark: vertex processed, avoid double split
          Cerr << "-> VerifierCoin is applied on the node " << somm << " of coordinates: ";
          for (int dir = 0; dir < dimension; dir++) Cerr << sommets(somm, dir) << " ";
          Cerr << "..." << finl;
          cut_elem(elem, xp);
        }
    }

  Scatter::init_sequential_domain(dom);

  return is;
}
