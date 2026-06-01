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

#include <Extraire_domaine.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <NettoieNoeuds.h>
#include <Sous_Domaine.h>
#include <Parser_U.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Scatter.h>
#include <Param.h>
#include <ParserView.h>

Implemente_instanciable(Extraire_domaine,"Extraire_Domaine",Interprete_geometrique_base);
// XD extraire_domaine interprete extraire_domaine BRACE Keyword to create a new domain built with the domain elements
// XD_CONT of the pb_name problem verifying the two conditions given by Condition_elements. The problem pb_name should
// XD_CONT have been discretized.

Sortie& Extraire_domaine::printOn(Sortie& os) const { return Interprete::printOn(os); }
Entree& Extraire_domaine::readOn(Entree& is) { return Interprete::readOn(is); }

Entree& Extraire_domaine::interpreter_(Entree& is)
{
  Nom nom_pb;
  Nom nom_dom;
  Nom nom_sous_domaine;
  Parser_U condition_elements;
  Nom expr_elements("1");
  condition_elements.setNbVar(dimension);
  condition_elements.addVar("x");
  condition_elements.addVar("y");
  if (dimension==3)
    condition_elements.addVar("z");

  Param param(que_suis_je());
  param.ajouter("domaine",&nom_dom,Param::REQUIRED); // XD_ADD_P ref_domaine
  // XD_CONT Domain in which faces are saved
  param.ajouter("probleme",&nom_pb,Param::REQUIRED); // XD_ADD_P ref_Pb_base
  // XD_CONT Problem from which faces should be extracted
  param.ajouter("condition_elements",&expr_elements); // XD_ADD_P chaine
  // XD_CONT not_set
  param.ajouter("sous_domaine|sous_zone",&nom_sous_domaine); // XD_ADD_P ref_sous_zone
  // XD_CONT not_set
  param.lire_avec_accolades_depuis(is);

  condition_elements.setString(expr_elements);
  condition_elements.parseString();

  associer_domaine(nom_dom);
  Domaine& dom=domaine();
  // retrieve the problem
  if(! sub_type(Probleme_base, objet(nom_pb)))
    {
      Cerr << nom_pb << " is of type " << objet(nom_pb).que_suis_je() << finl;
      Cerr << "and not of type Probleme_base" << finl;
      exit();
    }
  Probleme_base& pb=ref_cast(Probleme_base, objet(nom_pb));
  const Domaine_VF& domaine_vf=ref_cast(Domaine_VF,pb.domaine_dis());
  dom.les_sommets()=domaine_vf.domaine().les_sommets();
  dom.typer(domaine_vf.domaine().type_elem()->que_suis_je());

  int nb_elem=domaine_vf.nb_elem();
  IntTrav tab_marq_elem;
  domaine_vf.domaine().creer_tableau_elements(tab_marq_elem);
  IntArrView marq_elem = static_cast<ArrOfInt&>(tab_marq_elem).view_wo();
  int nb_elem_m=0;
  // mark the elements of interest
  if (nom_sous_domaine== Nom())
    {
      int dim = Objet_U::dimension;
      ParserView parser_condition_elements(condition_elements);
      parser_condition_elements.parseString();
      CDoubleTabView xp = domaine_vf.xp().view_ro();
      Kokkos::parallel_reduce(start_gpu_timer(__KERNEL_NAME__), nb_elem, KOKKOS_LAMBDA(const int elem, int& local_nb_elem_m)
      {
        int threadId = parser_condition_elements.acquire();
        parser_condition_elements.setVar(0,xp(elem,0),threadId);
        parser_condition_elements.setVar(1,xp(elem,1),threadId);
        if (dim==3)
          parser_condition_elements.setVar(2,xp(elem,2),threadId);
        double res=parser_condition_elements.eval(threadId);
        parser_condition_elements.release(threadId);
        if (std::fabs(res)>1e-5)
          {
            marq_elem(elem)=1;
            local_nb_elem_m++;
          }
        else
          marq_elem(elem)=0;
      }, nb_elem_m);
      end_gpu_timer(__KERNEL_NAME__);
    }
  else
    {
      const Sous_Domaine& ssz= ref_cast(Sous_Domaine,objet(nom_sous_domaine));
      CIntArrView les_elems = ssz.les_elems().view_ro();
      Kokkos::parallel_reduce(start_gpu_timer(__KERNEL_NAME__), ssz.nb_elem_tot(), KOKKOS_LAMBDA(const int i, int& local_nb_elem_m)
      {
        int elem = les_elems(i);
        if (elem<nb_elem)
          {
            marq_elem(elem)=1;
            local_nb_elem_m++;
          }
      }, nb_elem_m);
      end_gpu_timer(__KERNEL_NAME__);
    }
  // Note: trick - we do not exchange virtual spaces so that the joint becomes a boundary
  // marq_elem.echange_espace_virtuel();
  int nb_faces=domaine_vf.nb_faces();
  // find the faces on the boundary of the domain (joints included)
  int nb_t=0;
  ArrOfInt tab_marq(nb_faces);
  IntArrView marq = tab_marq.view_rw();
  CIntTabView face_voisin = domaine_vf.face_voisins().view_ro();
  Kokkos::parallel_reduce(start_gpu_timer(__KERNEL_NAME__), nb_faces, KOKKOS_LAMBDA(const int fac, int& local_nb_t)
  {
    int elem0 = face_voisin(fac, 0);
    int elem1 = face_voisin(fac, 1);
    int val0 = -1;
    if (elem0 != -1) val0 = marq_elem(elem0);
    int val1 = -1;
    if (elem1 != -1) val1 = marq_elem(elem1);

    if (val0 != val1)
      {
        if ((val0 == 1) || (val1 == 1))
          {
            if (marq[fac] != -1) //not a joint
              {
                marq[fac] = 1;
                local_nb_t++;
              }
          }
      }
  }, nb_t);
  end_gpu_timer(__KERNEL_NAME__);

  IntTab& les_elems=dom.les_elems();
  const IntTab& les_elems_old=domaine_vf.domaine().les_elems();
  int nb_som_elem=les_elems_old.dimension(1);
  les_elems.resize(nb_elem_m,nb_som_elem);
  const IntTab& face_sommets=domaine_vf.face_sommets();
  int nb=0;
  ToDo_Kokkos("critical");
  for (int elem=0; elem<nb_elem; elem++)
    {
      if (tab_marq_elem(elem)==1)
        {
          for (int k=0; k<nb_som_elem; k++)
            les_elems(nb,k)=les_elems_old(elem,k);
          nb++;
        }
    }
  assert(nb==nb_elem_m);
  Bord faces;
  faces.nommer("Bord");
  faces.typer_faces(domaine_vf.domaine().type_elem()->type_face());
  Faces& les_faces=faces.faces();
  int nb_som_face=face_sommets.dimension(1);
  IntTab& indfaces=les_faces.les_sommets();
  indfaces.resize(nb_t,nb_som_face);
  IntTab& facesv=les_faces.voisins();
  facesv.resize(nb_t,2);
  facesv=-1;
  nb=0;
  ToDo_Kokkos("critical");
  for (int fac=0; fac<nb_faces; fac++)
    if (tab_marq[fac]==1)
      {
        for (int s=0; s<nb_som_face; s++)
          indfaces(nb,s)=face_sommets(fac,s);
        // compute the normal
        nb++;
      }
  Cerr<<finl;;
  dom.faces_bord().add(faces);
  dom.faces_bord().associer_domaine(dom);
  dom.type_elem()->associer_domaine(dom);
  dom.fixer_premieres_faces_frontiere();

  assert(nb==nb_t);
  NettoieNoeuds::nettoie(dom);

  return is;
}
