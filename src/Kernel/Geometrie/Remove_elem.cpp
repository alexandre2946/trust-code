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

#include <NettoieNoeuds.h>
#include <Remove_elem.h>
#include <TRUSTLists.h>
#include <TRUSTArray.h>
#include <Scatter.h>
#include <Param.h>

Implemente_instanciable_32_64(Remove_elem_32_64, "Remove_elem", Interprete_geometrique_base_32_64<_T_>);
// XD remove_elem interprete remove_elem INHERITS_BRACE Keyword to remove element from a VDF mesh (named domaine_name),
// XD_CONT either from an explicit list of elements or from a geometric condition defined by a condition f(x,y)>0 in 2D
// XD_CONT and f(x,y,z)>0 in 3D. All the new borders generated are gathered in one boundary called : newBord (to rename
// XD_CONT it, use RegroupeBord keyword. To split it to different boundaries, use decoupebord keyword). Example of a
// XD_CONT removed zone of radius 0.2 centered at (x,y)=(0.5,0.5): NL2 Remove_elem dom { fonction
// XD_CONT $0.2*0.2-(x-0.5)^2-(y-0.5)^2>0$ } NL2 Warning : the thickness of removed zone has to be large enough to avoid
// XD_CONT singular nodes as decribed below : \includeimage{{removeelem.jpeg}}
// XD attr domaine ref_domaine domain REQ Name of domain
// XD attr bloc remove_elem_bloc bloc REQ not_set
template <typename _SIZE_>
Sortie& Remove_elem_32_64<_SIZE_>::printOn(Sortie& os) const { return Interprete::printOn(os); }

template <typename _SIZE_>
Entree& Remove_elem_32_64<_SIZE_>::readOn(Entree& is) { return Interprete::readOn(is); }

template <typename _SIZE_>
int Remove_elem_32_64<_SIZE_>::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  int retval = 1;
  if (mot == "liste")
    {
      int_t nb_elem, elem;
      is >> nb_elem;
      for (int_t i = 0; i < nb_elem; i++)
        {
          is >> elem;
          listelem.add(elem);
        }
    }
  else
    retval = -1;

  return retval;
}

// XD remove_elem_bloc objet_lecture nul BRACE not_set
template <typename _SIZE_>
Entree& Remove_elem_32_64<_SIZE_>::interpreter_(Entree& is)
{
  this->associer_domaine(is);
  Param param(this->que_suis_je());
  Nom fonction;
  param.ajouter_non_std("liste", this); // XD_ADD_P listentier
  // XD_CONT not_set
  param.ajouter("fonction", &fonction); // XD_ADD_P chaine
  // XD_CONT not_set
  param.lire_avec_accolades_depuis(is);
  if (fonction == "??")
    f_ok = 0;
  else
    {
      Cerr << "Reading and interpretation of the function " << fonction << " ... ";
      f.setNbVar(3);
      f.setString(fonction);
      f.addVar("x");
      f.addVar("y");
      f.addVar("z");
      f.parseString();
      f_ok = 1;
      Cerr << " Ok" << finl;
    }

  Scatter::uninit_sequential_domain(this->domaine());
  remove_elem_(this->domaine());
  Scatter::init_sequential_domain(this->domaine());
  NettoieNoeuds_32_64<_SIZE_>::nettoie(this->domaine());
  Cerr << "Refinement... OK" << finl;
  return is;
}

template <typename _SIZE_>
void Remove_elem_32_64<_SIZE_>::recreer_faces(Domaine_t& domaine, Faces_t& faces, IntTab_t& som_face) const
{
  IntTab_t& sommets = faces.les_sommets();
  int_t nb_faces = sommets.dimension(0);
  int nbs = (Objet_U::dimension == 2) ? 2 : 4;  // number of vertices per face
  IntTab_t faces_recreees(nb_faces, nbs);

  int_t ii = 0;

  for (int_t i = 0; i < nb_faces; i++)
    {
      ArrOfInt_t ind(4);
      ind[0] = sommets(i, 0);
      ind[1] = sommets(i, 1);
      ind[2] = (Objet_U::dimension == 3) ? sommets(i, 2) : -1;
      ind[3] = (Objet_U::dimension == 3) ? sommets(i, 3) : -1;
      ind.ordonne_array(); // sort indices to avoid unexpected index ordering surprises

      int trouve = 0;

      int j = 0;
      for (; j < som_face.dimension(2); j++)
        {
          if (som_face(ind[3], 0, j) == ind[2] && som_face(ind[3], 1, j) == ind[1] && som_face(ind[3], 2, j) == ind[0])
            {
              trouve = 1;
              break;
            }
        }

      if (trouve == 1) // "remove" the face by resetting the indices to -1
        {
          som_face(ind[3], 0, j) = -1;
          som_face(ind[3], 1, j) = -1;
          som_face(ind[3], 2, j) = -1;
        }
      else
        {
          faces_recreees(ii, 0) = sommets(i, 0);
          faces_recreees(ii, 1) = sommets(i, 1);
          if (Objet_U::dimension == 3)
            {
              faces_recreees(ii, 2) = sommets(i, 2);
              faces_recreees(ii, 3) = sommets(i, 3);
            }
          ii++;
        }
    }

  faces_recreees.resize(ii, nbs);

  sommets.reset();
  sommets.ref(faces_recreees);
}

template <typename _SIZE_>
void Remove_elem_32_64<_SIZE_>::creer_faces(Domaine_t& dom, Faces_t& faces, IntTab_t& som_face) const
{
  faces.dimensionner(1);
  IntTab_t& sommets = faces.les_sommets();
  int_t nbsom = this->domaine().les_sommets().dimension(0);
  int nbs = (Objet_U::dimension == 2) ? 2 : 4;  // number of vertices per face
  IntTab_t faces_recreees(1, nbs);

  int ii = 0;

  for (int i = 0; i < nbsom; i++)
    {
      for (int j = 0; j < som_face.dimension(2); j++)
        {
          if (som_face(i, 0, j) != -1)
            {
              faces_recreees.resize(ii + 1, nbs);
              faces_recreees(ii, 0) = i;
              faces_recreees(ii, 1) = som_face(i, 0, j);
              if (Objet_U::dimension == 3)
                {
                  faces_recreees(ii, 2) = som_face(i, 1, j);
                  faces_recreees(ii, 3) = som_face(i, 2, j);
                }
              ii++;
            }
        }
    }

  faces.dimensionner(ii);
  sommets.ref(faces_recreees);
}

template <typename _SIZE_>
void Remove_elem_32_64<_SIZE_>::remplir_liste(IntTab_t& som_face, int_t ind1, int_t ind2, int_t ind3, int_t ind4) const
{
  ArrOfInt_t ind(4);
  ind[0] = ind1;
  ind[1] = ind2;
  ind[2] = ind3;
  ind[3] = ind4;
  ind.ordonne_array(); // sort indices to avoid unexpected index ordering surprises (previously encountered)

  int trouve = 0;

  int j = 0;
  for (; j < som_face.dimension(2); j++)
    {
      if (som_face(ind[3], 0, j) == ind[2] && som_face(ind[3], 1, j) == ind[1] && som_face(ind[3], 2, j) == ind[0])
        {
          trouve = 1;
          break;
        }
    }

  if (trouve == 1) // "remove" the face by resetting the relevant indices to -1
    {
      som_face(ind[3], 0, j) = -1;
      som_face(ind[3], 1, j) = -1;
      som_face(ind[3], 2, j) = -1;
    }
  else // "add" the face by assigning indices to the first "free" slots (=-1)
    {
      for (j = 0; j < som_face.dimension(2); j++)
        {
          if (som_face(ind[3], 0, j) == -1)
            break;
        }
      som_face(ind[3], 0, j) = ind[2];
      som_face(ind[3], 1, j) = ind[1];
      som_face(ind[3], 2, j) = ind[0];
    }
}

template <typename _SIZE_>
void Remove_elem_32_64<_SIZE_>::remove_elem_(Domaine_t& dom)
{
  if (dom.type_elem()->que_suis_je() == "Rectangle" || dom.type_elem()->que_suis_je() == "Rectangle_64"
      || dom.type_elem()->que_suis_je() == "Hexaedre" || dom.type_elem()->que_suis_je() == "Hexaedre_64")
    {

      IntTab_t& les_elems = dom.les_elems();
      int_t oldsz = les_elems.dimension(0);
      ArrOfInt_t marq_remove(oldsz);
      int_t nbsom = this->domaine().les_sommets().dimension(0);

      int nbs = (Objet_U::dimension == 2) ? 4 : 8;  // number of vertices per element
      IntTab_t new_elems(oldsz, nbs);

      int nbfacesom = (Objet_U::dimension == 2) ? 4 : 4 * 3; // number of faces connected to a vertex
      IntTab_t som_face(nbsom, 3, nbfacesom);
      som_face = -1;

      if (f_ok)
        {
          DoubleTab_t xg(oldsz, Objet_U::dimension);
          dom.type_elem()->calculer_centres_gravite(xg);
          for (int_t i = 0; i < oldsz; i++)
            {
              f.setVar(0, xg(i, 0));
              f.setVar(1, xg(i, 1));
              if (Objet_U::dimension == 3)
                f.setVar(2, xg(i, 2));
              //if(f.eval()) listelem.add(i);
              if ((int) (f.eval() + 0.5))
                marq_remove[i] = 1; //listelem.add(i); // to be consistent with what is done in DecoupeBord
            }
        }
      else
        {
          for (int i = 0; i < listelem.size(); i++)
            {
              marq_remove[listelem[i]] = 1;
            }

        }

      int_t j = 0;
      Cerr << "-> " << listelem.size() << " elements will be removed from the domain " << this->domaine().le_nom() << finl;
      /*    if (listelem.size()==0)
       {
       Cerr << "May be an error when applying Remove_elem : no elements found." << finl;
       Process::exit();
       } */
      for (int_t i = 0; i < oldsz; i++)
        {
          if (marq_remove[i] == 0)
            {
              for (int k = 0; k < nbs; k++)
                new_elems(j, k) = les_elems(i, k);
              j++;
            }
          else
            {
              if (Objet_U::dimension == 2)
                {
                  int_t i0 = les_elems(i, 0);
                  int_t i1 = les_elems(i, 1);
                  int_t i2 = les_elems(i, 2);
                  int_t i3 = les_elems(i, 3);

                  this->remplir_liste(som_face, i0, i1, -1, -1);
                  this->remplir_liste(som_face, i0, i2, -1, -1);
                  this->remplir_liste(som_face, i1, i3, -1, -1);
                  this->remplir_liste(som_face, i2, i3, -1, -1);
                }
              else
                {
                  int_t i0 = les_elems(i, 0);
                  int_t i1 = les_elems(i, 1);
                  int_t i2 = les_elems(i, 2);
                  int_t i3 = les_elems(i, 3);
                  int_t i4 = les_elems(i, 4);
                  int_t i5 = les_elems(i, 5);
                  int_t i6 = les_elems(i, 6);
                  int_t i7 = les_elems(i, 7);

                  this->remplir_liste(som_face, i0, i1, i2, i3);
                  this->remplir_liste(som_face, i0, i1, i4, i5);
                  this->remplir_liste(som_face, i0, i4, i2, i6);
                  this->remplir_liste(som_face, i1, i5, i3, i7);
                  this->remplir_liste(som_face, i2, i3, i6, i7);
                  this->remplir_liste(som_face, i4, i5, i6, i7);
                }
            }
        }

      new_elems.resize(j, nbs);
      les_elems.ref(new_elems);

      // Rebuild the octree
      dom.invalide_octree();
      dom.construit_octree();
      dom.reordonner();

      {
        Cerr << " Regeneration of boundaries" << finl;
        for (auto &itr : dom.faces_bord())
          {
            Faces_t& les_faces = itr.faces();
            if (Objet_U::dimension == 2)
              les_faces.typer(Type_Face::segment_2D);
            else
              les_faces.typer(Type_Face::quadrangle_3D);
            this->recreer_faces(dom, les_faces, som_face);
          }
        Cerr << " addition of a new boundary issued from removed elements" << finl;
        Bord_32_64<_SIZE_>& new_bord = dom.faces_bord().add(Bord_32_64<_SIZE_>());
        new_bord.nommer("newBord");
        if (Objet_U::dimension == 2)
          new_bord.typer_faces(Type_Face::segment_2D);
        else
          new_bord.typer_faces(Type_Face::quadrangle_3D);
        Faces_t& les_faces = new_bord.faces();
        creer_faces(dom, les_faces, som_face);
      }

      {
        // Internal boundaries
        Cerr << "Regeneration of internal faces" << finl;
        for (auto &itr : dom.bords_int())
          {
            Faces_t& les_faces = itr.faces();
            if (Objet_U::dimension == 2)
              les_faces.typer(Type_Face::segment_2D);
            else
              les_faces.typer(Type_Face::quadrangle_3D);
            recreer_faces(dom, les_faces, som_face);
          }
      }

      Cerr << "END of Remove_elem..." << finl;
      Cerr << "  1 NbElem=" << dom.les_elems().dimension(0) << "  NbNod=" << dom.nb_som() << finl;
    }

  else

    {
      Cerr << "We do not yet know how to Remove_elem the " << dom.type_elem()->que_suis_je() << "s" << finl;
      Process::exit();
    }
}

template class Remove_elem_32_64<int>;
#if INT_is_64_ == 2
template class Remove_elem_32_64<trustIdType>;
#endif



