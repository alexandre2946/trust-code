/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Tetra_EF.h>
#include <Domaine.h>
#include <Domaine_EF.h>
#include <Champ_P1_EF.h>
#include <Equation_base.h>
#include <Milieu_base.h>

Implemente_instanciable_sans_constructeur(Tetra_EF,"Tetra_EF",Elem_EF_base);

// printOn and readOn


Sortie& Tetra_EF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

Entree& Tetra_EF::readOn(Entree& s )
{
  return s ;
}
/*! @brief Returns for facet fa7: for j=0, j=1: the local indices of the 2 faces surrounding fa7.
 *
 * @brief For j=2, j=3: the local indices of the tetrahedron vertices belonging to fa7.
 *
 */
Tetra_EF::Tetra_EF()
{
}

void Tetra_EF::normale(int num_Face,DoubleTab& Face_normales,
                       const IntTab& Face_sommets,
                       const IntTab& Face_voisins,
                       const IntTab& elem_faces,
                       const Domaine& domaine_geom) const
{

  //Cerr << " num_Face " << num_Face << finl;
  const DoubleTab& les_coords = domaine_geom.coord_sommets();

  // Cerr << "les face sommet " << Face_sommets << finl;
  double x1,y1,z1,x2,y2,z2;
  double nx,ny,nz;
  int f0,no4;

  int n0 = Face_sommets(num_Face,0);
  int n1 = Face_sommets(num_Face,1);
  int n2 = Face_sommets(num_Face,2);


  x1 = les_coords(n0,0) - les_coords(n1,0);
  y1 = les_coords(n0,1) - les_coords(n1,1);
  z1 = les_coords(n0,2) - les_coords(n1,2);

  x2 = les_coords(n2,0) - les_coords(n1,0);
  y2 = les_coords(n2,1) - les_coords(n1,1);
  z2 = les_coords(n2,2) - les_coords(n1,2);

  nx = (y1*z2 - y2*z1)/2;
  ny = (-x1*z2 + x2*z1)/2;
  nz = (x1*y2 - x2*y1)/2;
  // Cerr << "nx " << nx << " ny " << ny << " nz " << nz << finl;

  // Orient the normal from elem1 toward elem2
  // by searching for the vertex of elem1 that is not on the Face
  int elem1 = Face_voisins(num_Face,0);
  if ( (f0 = elem_faces(elem1,0)) == num_Face )
    f0 = elem_faces(elem1,1);

  if ( (no4 = Face_sommets(f0,0)) != n0    &&   no4 != n1
       &&   no4 != n2)
    { /* Do nothing */}
  else if ( (no4 = Face_sommets(f0,1)) != n0 && no4 != n1
            && no4 != n2 )
    { /* Do nothing */}
  else
    no4 = Face_sommets(f0,2);

  x1 = les_coords(no4,0) - les_coords(n0,0);
  y1 = les_coords(no4,1) - les_coords(n0,1);
  z1 = les_coords(no4,2) - les_coords(n0,2);

  if ( (nx*x1+ny*y1+nz*z1) > 0 )
    {
      Face_normales(num_Face,0) = - nx;
      Face_normales(num_Face,1) = - ny;
      Face_normales(num_Face,2) = - nz;
    }
  else
    {
      Face_normales(num_Face,0) = nx;
      Face_normales(num_Face,1) = ny;
      Face_normales(num_Face,2) = nz;
    }

  // Cerr << "Face_normales " << Face_normales << finl;

}

void Tetra_EF::calcul_vc(const ArrOfInt& Face,ArrOfDouble& vc,
                         const ArrOfDouble& vs,const DoubleTab& vsom,
                         const Champ_Inc_base& vitesse,int type_cl) const
{
  int comp;
  const DoubleVect& porosite_face = vitesse.equation().milieu().porosite_face();
  switch(type_cl)
    {

    case 0: // the tetrahedron has no Dirichlet face
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.25*vs[comp];
        break;
      }

    case 1: // the tetrahedron has one Dirichlet face: KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse.valeurs()(Face[3],comp)*porosite_face[Face[3]];
        break;
      }

    case 2: // the tetrahedron has one Dirichlet face: KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse.valeurs()(Face[2],comp)*porosite_face[Face[2]];
        break;
      }

    case 4: // the tetrahedron has one Dirichlet face: KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse.valeurs()(Face[1],comp)*porosite_face[Face[1]];
        break;
      }

    case 8: // the tetrahedron has one Dirichlet face: KEL0
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse.valeurs()(Face[0],comp)*porosite_face[Face[0]];
        break;
      }

    case 3: // the tetrahedron has two Dirichlet faces: KEL3 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom(0,comp) + vsom(1,comp));
        break;
      }

    case 5: // the tetrahedron has two Dirichlet faces: KEL3 and KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom(0,comp) + vsom(2,comp));
        break;
      }

    case 6: // the tetrahedron has two Dirichlet faces: KEL1 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom(0,comp) + vsom(3,comp));
        break;
      }

    case 9: // the tetrahedron has two Dirichlet faces: KEL0 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom(1,comp) + vsom(2,comp));
        break;
      }

    case 10: // the tetrahedron has two Dirichlet faces: KEL0 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom(1,comp) + vsom(3,comp));
        break;
      }

    case 12: // the tetrahedron has two Dirichlet faces: KEL0 and KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5*(vsom(2,comp) + vsom(3,comp));
        break;
      }

    case 7: // the tetrahedron has three Dirichlet faces: KEL1, KEL2 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom(0,comp);
        break;
      }

    case 11: // the tetrahedron has three Dirichlet faces: KEL0, KEL2 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom(1,comp);
        break;
      }

    case 13: // the tetrahedron has three Dirichlet faces: KEL0, KEL1 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom(2,comp);
        break;
      }

    case 14: // the tetrahedron has three Dirichlet faces: KEL0, KEL1 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom(3,comp);
        break;
      }

    } // end of switch
}

/*! @brief Computes the coordinates xg of the centre of a non-standard element.
 *
 * @brief Also computes idirichlet = number of Dirichlet faces of the element.
 * @param xg Output centre coordinates.
 * @param x Vertex coordinate table for the element.
 * @param type_elem_Cl Element boundary condition type.
 * @param idirichlet Output number of Dirichlet faces.
 * @param n1 Output first null facet index (when idirichlet >= 2).
 * @param n2 Output second null facet index (when idirichlet >= 3).
 * @param n3 Output third null facet index (when idirichlet == 3).
 */
void Tetra_EF::calcul_xg(DoubleVect& xg,const DoubleTab& x, const int type_elem_Cl,
                         int& idirichlet,int& n1,int& n2,int& n3) const
{
  int j,dim=xg.size();

  switch(type_elem_Cl)
    {

    case 0:  // the tetrahedron has no Dirichlet face; it has 6 facets
      {
        for (j=0; j<dim; j++)
          xg[j]=0.25*(x(0,j)+x(1,j)+x(2,j)+x(3,j));

        idirichlet=0;
        break;
      }

    case 1:  // the tetrahedron has one Dirichlet face. The 'centre'
      // of the tetrahedron is at the midpoint of face 3 (vertices 0, 1, 2).
      // It has 3 real facets: 0 at nodes 2 3 xg
      //                       1 at nodes 1 3 xg
      //                       3 at nodes 3 0 xg
      // the 3 other facets lie on face 3

      {
        for (j=0; j<dim; j++)
          xg[j]=(x(0,j)+x(1,j)+x(2,j))/3.;

        idirichlet=1;
        break;

      }

    case 2:  // the tetrahedron has one Dirichlet face. The 'centre'
      // of the tetrahedron is at the midpoint of face 2 (vertices 0, 1, 3).
      // It has 3 real facets: 0 at nodes 2 3 xg
      //                       2 at nodes 1 2 xg
      //                       4 at nodes 2 0 xg

      {
        for (j=0; j<dim; j++)
          xg[j]=(x(0,j)+x(1,j)+x(3,j))/3.;

        idirichlet=1;
        break;
      }

    case 4:  // the tetrahedron has one Dirichlet face. The 'centre'
      // of the tetrahedron is at the midpoint of face 1 (vertices 0, 2, 3).
      // It has 3 real facets: 1 at nodes 1 3 xg
      //                       2 at nodes 1 2 xg
      //                       5 at nodes 1 0 xg

      {
        for (j=0; j<dim; j++)
          xg[j]=(x(0,j)+x(2,j)+x(3,j))/3.;

        idirichlet=1;
        break;
      }

    case 8:  // the tetrahedron has one Dirichlet face. The 'centre'
      // of the tetrahedron is at the midpoint of face 0 (vertices 1, 2, 3).
      // It has 3 real facets: 3 at nodes 3 0 xg
      //                       4 at nodes 2 0 xg
      //                       5 at nodes 1 0 xg

      {
        for (j=0; j<dim; j++)
          xg[j]=(x(1,j)+x(2,j)+x(3,j))/3.;

        idirichlet=1;
        break;
      }

    case 3:  // the tetrahedron has two Dirichlet faces 2 and 3. The 'centre'
      // is at the midpoint of the edge with endpoints 0 and 1.
      // It has 1 null facet: 5

      {
        for (j=0; j<dim; j++)
          xg[j]= 0.5*(x(0,j)+x(1,j));

        n1=5;
        idirichlet=2;
        break;
      }


    case 5:  // the tetrahedron has two Dirichlet faces 3 and 1. The 'centre'
      // is at the midpoint of the edge with endpoints 0 and 2.
      // It has 1 null facet: 4

      {
        for (j=0; j<dim; j++)
          xg[j]= 0.5*(x(0,j)+x(2,j));

        n1=4;
        idirichlet=2;
        break;
      }

    case 6:  // the tetrahedron has two Dirichlet faces 1 and 2. The 'centre'
      // is at the midpoint of the edge with endpoints 0 and 3.
      // It has 1 null facet: 3

      {
        for (j=0; j<dim; j++)
          xg[j]= 0.5*(x(0,j)+x(3,j));

        n1=3;
        idirichlet=2;
        break;
      }

    case 9:  // the tetrahedron has two Dirichlet faces 0 and 3. The 'centre'
      // is at the midpoint of the edge with endpoints 1 and 2.
      // It has 1 null facet: 2

      {
        for (j=0; j<dim; j++)
          xg[j]= 0.5*(x(1,j)+x(2,j));

        n1=2;
        idirichlet=2;
        break;
      }

    case 10:  // the tetrahedron has two Dirichlet faces 0 and 2. The 'centre'
      // is at the midpoint of the edge with endpoints 1 and 3.
      // It has 1 null facet: 1

      {
        for (j=0; j<dim; j++)
          xg[j]= 0.5*(x(1,j)+x(3,j));

        n1=1;
        idirichlet=2;
        break;
      }


    case 12:  // the tetrahedron has two Dirichlet faces 0 and 1. The 'centre'
      // is at the midpoint of the edge with vertices 2 and 3.
      // It has 1 null facet

      {
        for (j=0; j<dim; j++)
          xg[j]= 0.5*(x(2,j)+x(3,j));

        n1=0;
        idirichlet=2;
        break;
      }

    case 7:  // three Dirichlet faces: 1, 2, 3. The centre is at vertex 0.
      // There are 3 null facets: 3, 4, 5

      {
        for (j=0; j<dim; j++)
          xg[j]= x(0,j);

        n1=3;
        n2=4;
        n3=5;
        idirichlet=3;
        break;

      }

    case 11:  // three Dirichlet faces: 0, 2, 3. The centre is at vertex 1.
      // There are 3 null facets: 1, 2, 5

      {
        for (j=0; j<dim; j++)
          xg[j]= x(1,j);

        n1=1;
        n2=2;
        n3=5;
        idirichlet=3;
        break;

      }

    case 13:  // three Dirichlet faces: 0, 1, 3. The centre is at vertex 2.
      // There are 3 null facets: 0, 2, 4

      {
        for (j=0; j<dim; j++)
          xg[j]= x(2,j);

        n1=0;
        n2=2;
        n3=4;
        idirichlet=3;
        break;

      }
    case 14:  // three Dirichlet faces: 0, 1, 2. The centre is at vertex 3.
      // There are 3 null facets: 0, 1, 3

      {
        for (j=0; j<dim; j++)
          xg[j]= x(3,j);

        n1=0;
        n2=1;
        n3=3;
        idirichlet=3;
        break;

      }
    }
}

