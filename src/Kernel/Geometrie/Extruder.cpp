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

#include <Frontiere.h>
#include <Connectivite_som_elem.h>
#include <Static_Int_Lists.h>
#include <Faces_builder.h>
#include <Extruder.h>
#include <Domaine.h>
#include <Scatter.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur_32_64(Extruder_32_64, "Extruder", Interprete_geometrique_base_32_64<_T_>);
// XD extruder interprete extruder BRACE Class to create a 3D tetrahedral/hexahedral mesh (a prism is cut in 14) from a
// XD_CONT 2D triangular/quadrangular mesh.

template <typename _SIZE_>
Extruder_32_64<_SIZE_>::Extruder_32_64() { direction.resize(3, RESIZE_OPTIONS::NOCOPY_NOINIT); }

template <typename _SIZE_>
Sortie&  Extruder_32_64<_SIZE_>::printOn(Sortie& os) const { return Interprete::printOn(os); }

template <typename _SIZE_>
Entree&  Extruder_32_64<_SIZE_>::readOn(Entree& is) { return Interprete::readOn(is); }

/*! @brief Main function of the Extruder interpreter. Extrudes the domain
 *
 *     specified by the directive, element by element.
 *     The domain is extruded using the method:
 *       void Extruder_32_64<_SIZE_>::extruder(Domaine_t& domaine) const
 *     Extruding here means transforming geometric elements of a domain
 *     into 3D elements by translation.
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the input stream
 * @throws the object to mesh is not of type Domaine
 */
template <typename _SIZE_>
Entree&  Extruder_32_64<_SIZE_>::interpreter_(Entree& is)
{
  Nom nom_dom;
  Param param(this->que_suis_je());
  param.ajouter("domaine",&nom_dom,Param::REQUIRED);  // XD attr domaine ref_domaine domain_name REQ Name of the domain.
  param.ajouter("nb_tranches",&NZ,Param::REQUIRED);   // XD attr nb_tranches entier nb_tranches REQ Number of elements
  // XD_CONT in the extrusion direction.
  param.ajouter_arr_size_predefinie("direction",&direction,Param::REQUIRED); // XD attr direction troisf direction REQ
  // XD_CONT Direction of the extrude operation.
  param.lire_avec_accolades_depuis(is);
  this->associer_domaine(nom_dom);
  Scatter::uninit_sequential_domain(this->domaine());
  extruder(this->domaine());
  Scatter::init_sequential_domain(this->domaine());
  return is;
}

inline void check_boundary_name(const Nom& name)
{
// Check that the boundary is not already named "devant" or "derriere"
  if (name=="devant" || name=="derriere")
    {
      Cerr << "Problem : you must change the name of the boundary  : " <<name<<finl;
      Cerr << "because the extrusion keyword is going to create new boundaries named derriere and devant" << finl;
      Cerr << "which will create a conflict." << finl;
      Process::exit();
    }
}
/*! @brief Extrudes all elements of a domain: transforms the geometric elements of the domain into 3D elements.
 *
 * Currently only Rectangles/Quadrangles and Triangles can be extruded.
 *
 * @param dom The domain whose elements are to be extruded.
 */
template <typename _SIZE_>
void Extruder_32_64<_SIZE_>::extruder(Domaine_t& dom)
{


  if(dom.type_elem()->que_suis_je() == "Rectangle" || dom.type_elem()->que_suis_je() ==  "Quadrangle"
      || dom.type_elem()->que_suis_je() == "Rectangle_64" || dom.type_elem()->que_suis_je() ==  "Quadrangle_64" )
    {
      extruder_hexa(dom);
    }
  else if( dom.type_elem()->que_suis_je() == "Triangle" || dom.type_elem()->que_suis_je() == "Triangle_64")
    {
      int_t oldnbsom = dom.nb_som();
      IntTab_t& les_elems=dom.les_elems();
      int_t oldsz=les_elems.dimension(0);
      double dx = direction[0]/NZ;
      double dy = direction[1]/NZ;
      double dz = direction[2]/NZ;

      Faces_t les_faces;
      //domaine.creer_faces(les_faces);
      {
        // block to be factored out with Domaine_VF.cpp:
        Type_Face type_face = dom.type_elem()->type_face(0);
        les_faces.typer(type_face);
        les_faces.associer_domaine(dom);

        Static_Int_Lists_t connectivite_som_elem;
        const int_t     nb_sommets_tot = dom.nb_som_tot();
        const IntTab_t&   elements       = dom.les_elems();

        construire_connectivite_som_elem(nb_sommets_tot,
                                         elements,
                                         connectivite_som_elem,
                                         1 /* include virtual elements */);

        Faces_builder_t faces_builder;
        IntTab_t elem_faces; // Array that will not be needed
        faces_builder.creer_faces_reeles(dom,
                                         connectivite_som_elem,
                                         les_faces,
                                         elem_faces);
      }
      const int_t nbfaces2D = les_faces.nb_faces();


      int_t newnbsom = (oldnbsom*(NZ+1)+NZ*oldsz+nbfaces2D*NZ);
      DoubleTab_t new_soms(newnbsom, 3);
      DoubleTab_t& coord_sommets=dom.les_sommets();
      Objet_U::dimension=3;


      // vertices of the 2D mesh are translated first
      for (int i=0; i<oldnbsom; i++)
        {
          double x = coord_sommets(i,0);
          double y = coord_sommets(i,1);
          double z=0.;
          if (coord_sommets.dimension(1)>2)
            z=coord_sommets(i,2);

          for (int k=0; k<=NZ; k++)
            {
              new_soms(k*oldnbsom+i,0)=x;
              new_soms(k*oldnbsom+i,1)=y;
              new_soms(k*oldnbsom+i,2)=z;

              x += dx;
              y += dy;
              z += dz;
            }
        }


      // then create the centroids of the 2D elements and translate these points
      for (int_t i=0; i<oldsz; i++)
        {
          int_t i0=les_elems(i,0);
          int_t i1=les_elems(i,1);
          int_t i2=les_elems(i,2);

          double xg = 1./3.*(coord_sommets(i0,0)+coord_sommets(i1,0)+coord_sommets(i2,0))+0.5*dx;
          double yg = 1./3.*(coord_sommets(i0,1)+coord_sommets(i1,1)+coord_sommets(i2,1))+0.5*dy;
          double z = 0.5*dz;
          if (coord_sommets.dimension(1)>2)
            z = 1./3.*(coord_sommets(i0,2)+coord_sommets(i1,2)+coord_sommets(i2,2))+0.5*dz;
          for (int k=0; k<NZ; k++)
            {

              new_soms((oldnbsom*(NZ+1)+k*oldsz+i),0)=xg;
              new_soms((oldnbsom*(NZ+1)+k*oldsz+i),1)=yg;
              new_soms((oldnbsom*(NZ+1)+k*oldsz+i),2)=z;

              xg += dx;
              yg += dy;
              z += dz;
            }
        }


      // finally, create the face centers of the 2D mesh and translate these points
      for (int i=0; i<nbfaces2D; i++)
        {
          int_t i0=les_faces.sommet(i,0);
          int_t i1=les_faces.sommet(i,1);

          double x01 = 0.5*(coord_sommets(i0,0)+coord_sommets(i1,0))+0.5*dx;
          double y01 = 0.5*(coord_sommets(i0,1)+coord_sommets(i1,1))+0.5*dy;
          double z = 0.5*dz;
          if (coord_sommets.dimension(1)>2)
            z = 0.5*(coord_sommets(i0,2)+coord_sommets(i1,2))+0.5*dz;

          for (int k=0; k<NZ; k++)
            {
              new_soms((oldnbsom*(NZ+1)+NZ*oldsz+k*nbfaces2D+i),0)=x01;
              new_soms((oldnbsom*(NZ+1)+NZ*oldsz+k*nbfaces2D+i),1)=y01;
              new_soms((oldnbsom*(NZ+1)+NZ*oldsz+k*nbfaces2D+i),2)=z;
              x01 += dx;
              y01 += dy;
              z += dz;
            }
        }


      coord_sommets.resize(0);
      dom.ajouter(new_soms);

      int_t newnbelem = 14*NZ*oldsz;
      IntTab_t new_elems(newnbelem, 4); // the new elements
      int_t cpt=0;


      // first, store the top and bottom tetrahedra: NZ*nb_triangle*2 tetrahedra
      for (int_t i=0; i<oldsz; i++)
        {
          int_t i0=les_elems(i,0);
          int_t i1=les_elems(i,1);
          int_t i2=les_elems(i,2);

          int_t ig=oldnbsom*(NZ+1)+i;

          for (int k=0; k<NZ; k++)
            {
              new_elems(2*k*oldsz+2*i,0) = i0;
              new_elems(2*k*oldsz+2*i,1) = i1;
              new_elems(2*k*oldsz+2*i,2) = i2;
              new_elems(2*k*oldsz+2*i,3) = ig;
              cpt++;

              new_elems(2*k*oldsz+2*i+1,0) = i0+oldnbsom;
              new_elems(2*k*oldsz+2*i+1,1) = i1+oldnbsom;
              new_elems(2*k*oldsz+2*i+1,2) = i2+oldnbsom;
              new_elems(2*k*oldsz+2*i+1,3) = ig;
              cpt++;

              this->mettre_a_jour_sous_domaine(dom,i,(2*k*oldsz+2*i),2);

              i0+=oldnbsom;
              i1+=oldnbsom;
              i2+=oldnbsom;
              ig+=oldsz;
            }
        }



      // then the remaining tetrahedra
      for (int_t i=0; i<nbfaces2D; i++)
        {
          for (int ivois=0; ivois<2; ivois++)
            {
              int_t elem = les_faces.voisin(i,ivois);

              if (elem>=0)
                {
                  int_t i0=les_faces.sommet(i,0);
                  int_t i1=les_faces.sommet(i,1);
                  int_t i01=oldnbsom*(NZ+1)+NZ*oldsz+i;

                  for (int_t k=0; k<NZ; k++)
                    {
                      int_t ig=oldnbsom*(NZ+1)+k*oldsz+elem;

                      new_elems(cpt,0) = i0;
                      new_elems(cpt,1) = i1;
                      new_elems(cpt,2) = i01;
                      new_elems(cpt++,3) = ig;

                      new_elems(cpt,0) = i0+oldnbsom;
                      new_elems(cpt,1) = i1+oldnbsom;
                      new_elems(cpt,2) = i01;
                      new_elems(cpt++,3) = ig;

                      new_elems(cpt,0) = i1;
                      new_elems(cpt,1) = i1+oldnbsom;
                      new_elems(cpt,2) = i01;
                      new_elems(cpt++,3) = ig;

                      new_elems(cpt,0) = i0;
                      new_elems(cpt,1) = i0+oldnbsom;
                      new_elems(cpt,2) = i01;
                      new_elems(cpt++,3) = ig;


                      i0+=oldnbsom;
                      i1+=oldnbsom;
                      i01+=nbfaces2D;
                      ig+=oldsz;
                    }
                }
            }
        }

      les_elems.ref(new_elems);

      // Rebuild the octree
      dom.invalide_octree();
      dom.typer("Tetraedre");

      extruder_dvt(dom, les_faces,oldnbsom, oldsz);

    }
  else
    {
      Cerr << "It is not known yet how to extrude "
           << dom.type_elem()->que_suis_je() <<"s"<<finl;
      this->exit();
    }
}


template <typename _SIZE_>
void Extruder_32_64<_SIZE_>::traiter_faces_dvt(Faces_t& les_faces_bord, Faces_t& les_faces, int_t oldnbsom, int_t oldsz, int_t nbfaces2D )
{
  int_t size_2D = les_faces_bord.nb_faces();

  IntTab_t les_sommets(4*size_2D*NZ, 3);

  for (int_t i=0; i<size_2D; i++)
    {
      int_t i0=les_faces_bord.sommet(i,0);
      int_t i1=les_faces_bord.sommet(i,1);

      //double x01 = 0.5*(coord_sommets(i0,0)+coord_sommets(i1,0));
      //double y01 = 0.5*(coord_sommets(i0,1)+coord_sommets(i1,1));

      // find the index of this boundary face: not optimal!
      int_t jface=-1;
      for (int_t iface=0; iface<nbfaces2D; iface++)
        {
          int_t j0=les_faces.sommet(iface,0);
          int_t j1=les_faces.sommet(iface,1);

          if (((i0==j0) &&(i1==j1)) || ((i0==j1) &&(i1==j0)))
            {
              jface=iface;
              break;
            }
        }
      assert(jface>=0);

      for (int_t k=0; k<NZ; k++)
        {
          //double z = (k+0.5)*dz;
          //int i01 = domaine.chercher_sommets(x01, y01, z);
          int_t j01 = oldnbsom*(NZ+1)+NZ*oldsz+k*nbfaces2D+jface;

          les_sommets(k*4*size_2D+4*i,0) = i0;
          les_sommets(k*4*size_2D+4*i,1) = i1;
          les_sommets(k*4*size_2D+4*i,2) = j01;

          les_sommets(k*4*size_2D+4*i+1,0) = j01;
          les_sommets(k*4*size_2D+4*i+1,1) = i1;
          les_sommets(k*4*size_2D+4*i+1,2) = i1+oldnbsom;


          les_sommets(k*4*size_2D+4*i+2,0) = j01;
          les_sommets(k*4*size_2D+4*i+2,1) = i0+oldnbsom;
          les_sommets(k*4*size_2D+4*i+2,2) = i1+oldnbsom;


          les_sommets(k*4*size_2D+4*i+3,0) = j01;
          les_sommets(k*4*size_2D+4*i+3,1) = i0+oldnbsom;
          les_sommets(k*4*size_2D+4*i+3,2) = i0;

          i0+=oldnbsom;
          i1+=oldnbsom;
        }
    }

  les_faces_bord.typer(Type_Face::triangle_3D);
  les_faces_bord.les_sommets().ref(les_sommets);
  les_faces_bord.voisins().resize(4*size_2D*NZ, 2);
  les_faces_bord.voisins()=-1;
}


template <typename _SIZE_>
void Extruder_32_64<_SIZE_>::extruder_dvt(Domaine_t& dom, Faces_t& les_faces, int_t oldnbsom, int_t oldsz)
{

  const int_t nbfaces2D = les_faces.nb_faces();
  IntTab_t& les_elems=dom.les_elems();

  for (auto &itr : dom.faces_bord())
    {
      check_boundary_name(itr.le_nom());
      Faces_t& les_faces_bord = itr.faces();
      traiter_faces_dvt(les_faces_bord, les_faces, oldnbsom, oldsz, nbfaces2D);
    }

  for (auto &itr : dom.faces_raccord())
    {
      check_boundary_name(itr->le_nom());
      Faces_t& les_faces_bord = itr->faces();
      traiter_faces_dvt(les_faces_bord, les_faces, oldnbsom, oldsz, nbfaces2D);
    }

  Bord_t& devant = dom.faces_bord().add(Bord_t());
  devant.nommer("devant");
  Faces_t& les_faces_dvt=devant.faces();
  les_faces_dvt.typer(Type_Face::triangle_3D);

  IntTab_t som_dvt(oldsz, 3);
  les_faces_dvt.voisins().resize(oldsz, 2);
  les_faces_dvt.voisins()=-1;

  Bord_t& derriere = dom.faces_bord().add(Bord_t());
  derriere.nommer("derriere");
  Faces_t& les_faces_der=derriere.faces();
  les_faces_der.typer(Type_Face::triangle_3D);

  IntTab_t som_der(oldsz, 3);
  les_faces_der.voisins().resize(oldsz, 2);
  les_faces_der.voisins()=-1;

  for (int_t i=0; i<oldsz; i++)
    {
      int_t i0=les_elems(2*i,0);
      int_t i1=les_elems(2*i,1);
      int_t i2=les_elems(2*i,2);

      som_dvt(i,0) = i0;
      som_dvt(i,1) = i1;
      som_dvt(i,2) = i2;

      som_der(i,0) = i0+oldnbsom*NZ;
      som_der(i,1) = i1+oldnbsom*NZ;
      som_der(i,2) = i2+oldnbsom*NZ;

    }


  les_faces_dvt.les_sommets().ref(som_dvt);
  les_faces_der.les_sommets().ref(som_der);

}
template <typename _SIZE_>
void Extruder_32_64<_SIZE_>::extruder_hexa(Domaine_t& dom)
{

  int_t oldnbsom = dom.nb_som();
  IntTab_t& les_elems=dom.les_elems();
  int_t oldsz=les_elems.dimension(0);
  double dx = direction[0]/NZ;
  double dy = direction[1]/NZ;
  double dz = direction[2]/NZ;

  Faces_t les_faces;
  {
    // block to be factored out with Domaine_VF.cpp :
    Type_Face type_face = dom.type_elem()->type_face(0);
    les_faces.typer(type_face);
    les_faces.associer_domaine(dom);

    Static_Int_Lists_t connectivite_som_elem;
    const int_t     nb_sommets_tot = dom.nb_som_tot();
    const IntTab_t&    elements       = dom.les_elems();

    construire_connectivite_som_elem(nb_sommets_tot,
                                     elements,
                                     connectivite_som_elem,
                                     1 /* include virtual elements */);

    Faces_builder_t faces_builder;
    IntTab_t elem_faces; // Array that will not be needed
    faces_builder.creer_faces_reeles(dom,
                                     connectivite_som_elem,
                                     les_faces,
                                     elem_faces);
  }

  int_t newnbsom = oldnbsom*(NZ+1);
  DoubleTab_t new_soms(newnbsom, 3);
  DoubleTab_t& coord_sommets=dom.les_sommets();
  Objet_U::dimension=3;

  int_t i;
  // vertices of the 2D mesh are translated
  for (i=0; i<oldnbsom; i++)
    {
      double x = coord_sommets(i,0);
      double y = coord_sommets(i,1);
      double z=0.;
      if (coord_sommets.dimension(1)>2)
        z=coord_sommets(i,2);
      for (int k=0; k<=NZ; k++)
        {
          new_soms(k*oldnbsom+i,0)=x;
          new_soms(k*oldnbsom+i,1)=y;
          new_soms(k*oldnbsom+i,2)=z;

          x += dx;
          y += dy;
          z += dz;
        }
    }


  coord_sommets.resize(0);
  dom.ajouter(new_soms);

  int_t newnbelem = NZ*oldsz;
  IntTab_t new_elems(newnbelem, 8); // the new elements


  // define the new hexahedra
  for (i=0; i<oldsz; i++)
    {
      int_t i0=les_elems(i,0);
      int_t i1=les_elems(i,1);
      int_t i2=les_elems(i,2);
      int_t i3=les_elems(i,3);


      for (int_t k=0; k<NZ; k++)
        {
          new_elems(k*oldsz+i,0) = i0;
          new_elems(k*oldsz+i,1) = i1;
          new_elems(k*oldsz+i,2) = i2;
          new_elems(k*oldsz+i,3) = i3;
          new_elems(k*oldsz+i,4) = i0+oldnbsom;
          new_elems(k*oldsz+i,5) = i1+oldnbsom;
          new_elems(k*oldsz+i,6) = i2+oldnbsom;
          new_elems(k*oldsz+i,7) = i3+oldnbsom;


          i0+=oldnbsom;
          i1+=oldnbsom;
          i2+=oldnbsom;
          i3+=oldnbsom;
        }
    }

  les_elems.ref(new_elems);

  // Rebuild the octree
  dom.invalide_octree();
  if ((dom.type_elem()->que_suis_je()) ==  "Quadrangle")
    dom.typer("Hexaedre_VEF");
  else
    dom.typer("Hexaedre");

  extruder_dvt_hexa(dom, les_faces,oldnbsom, oldsz);
}
template <typename _SIZE_>
void Extruder_32_64<_SIZE_>::traiter_faces_dvt_hexa(Faces_t& les_faces_bord, int_t oldnbsom)
{
  int_t size_2D = les_faces_bord.nb_faces();

  IntTab_t les_sommets(size_2D*NZ, 4);

  for (int_t i=0; i<size_2D; i++)
    {
      int_t i0=les_faces_bord.sommet(i,0);
      int_t i1=les_faces_bord.sommet(i,1);


      for (int_t k=0; k<NZ; k++)
        {
          les_sommets(k*size_2D+i,0) = i0;
          les_sommets(k*size_2D+i,1) = i1;
          les_sommets(k*size_2D+i,2) = i0+oldnbsom;
          les_sommets(k*size_2D+i,3) = i1+oldnbsom;

          i0+=oldnbsom;
          i1+=oldnbsom;
        }
    }

  les_faces_bord.typer(Type_Face::quadrangle_3D);
  les_faces_bord.les_sommets().ref(les_sommets);
  les_faces_bord.voisins().resize(size_2D*NZ, 2);
  les_faces_bord.voisins()=-1;
}
template <typename _SIZE_>
void Extruder_32_64<_SIZE_>::extruder_dvt_hexa(Domaine_t& dom, Faces_t& les_faces, int_t oldnbsom, int_t oldsz)
{

  IntTab_t& les_elems=dom.les_elems();

  for (auto &itr : dom.faces_bord())
    {
      check_boundary_name(itr.le_nom());
      Faces_t& les_faces_bord = itr.faces();
      traiter_faces_dvt_hexa(les_faces_bord, oldnbsom);
    }

  for (auto &itr : dom.faces_raccord())
    {
      check_boundary_name(itr->le_nom());
      Frontiere_32_64<_SIZE_>& f = dynamic_cast<Frontiere_32_64<_SIZE_>&>(itr);
      Faces_t& les_faces_bord = f.faces();
      traiter_faces_dvt_hexa(les_faces_bord, oldnbsom);
    }

  Bord_t& devant = dom.faces_bord().add(Bord_t());
  devant.nommer("devant");
  Faces_t& les_faces_dvt=devant.faces();
  les_faces_dvt.typer(Type_Face::quadrangle_3D);

  IntTab_t som_dvt(oldsz, 4);
  les_faces_dvt.voisins().resize(oldsz, 2);
  les_faces_dvt.voisins()=-1;

  Bord_t& derriere = dom.faces_bord().add(Bord_t());
  derriere.nommer("derriere");
  Faces_t& les_faces_der=derriere.faces();
  les_faces_der.typer(Type_Face::quadrangle_3D);

  IntTab_t som_der(oldsz, 4);
  les_faces_der.voisins().resize(oldsz, 2);
  les_faces_der.voisins()=-1;

  for (int_t i=0; i<oldsz; i++)
    {
      int_t i0=les_elems(i,0);
      int_t i1=les_elems(i,1);
      int_t i2=les_elems(i,2);
      int_t i3=les_elems(i,3);

      som_dvt(i,0) = i0;
      som_dvt(i,1) = i1;
      som_dvt(i,2) = i2;
      som_dvt(i,3) = i3;

      som_der(i,0) = i0+oldnbsom*NZ;
      som_der(i,1) = i1+oldnbsom*NZ;
      som_der(i,2) = i2+oldnbsom*NZ;
      som_der(i,3) = i3+oldnbsom*NZ;

    }


  les_faces_dvt.les_sommets().ref(som_dvt);
  les_faces_der.les_sommets().ref(som_der);

}

template class Extruder_32_64<int>;
#if INT_is_64_ == 2
template class Extruder_32_64<trustIdType>;
#endif
