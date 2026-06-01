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

#include <Quadri_EF.h>
#include <Domaine.h>

Implemente_instanciable_sans_constructeur(Quadri_EF,"Quadri_EF",Elem_EF_base);

// printOn and readOn

Sortie& Quadri_EF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

Entree& Quadri_EF::readOn(Entree& s )
{
  return s ;
}

/*! @brief KEL_(0,fa7), KEL_(1,fa7) are the local indices of the 2 faces surrounding facet with local index fa7.
 *
 * @brief The local index of fa7 is that of the vertex that carries it.
 *
 */
Quadri_EF::Quadri_EF()
{
}

/*! @brief Fills the face_normales array in the Domaine_EF.
 *
 * @param num_Face Local face index.
 * @param Face_normales Array of face normals to fill.
 * @param Face_sommets Face-to-vertex connectivity table.
 * @param Face_voisins Face neighbour element table.
 * @param elem_faces Element-to-face connectivity table.
 * @param domaine_geom Geometric domain.
 */
void Quadri_EF::normale(int num_Face,DoubleTab& Face_normales,
                        const  IntTab& Face_sommets,
                        const IntTab& Face_voisins,
                        const IntTab& elem_faces,
                        const Domaine& domaine_geom) const
{
  const DoubleTab& les_coords = domaine_geom.coord_sommets();
  double x1,y1;
  double nx,ny;
  double x1g=0,y1g=0;
  double x2g=0,y2g=0;
  double grx,gry,psc;
  int sign=1,i;
  int n0 = Face_sommets(num_Face,0);
  int n1 = Face_sommets(num_Face,1);
  x1 = les_coords(n0,0)-les_coords(n1,0);
  y1 = les_coords(n0,1)-les_coords(n1,1);
  nx = -y1;
  ny = x1;
  int elem1=Face_voisins(num_Face,0);
  int elem2=Face_voisins(num_Face,1);

  // Orient the normal toward the element with the highest index.
  // First check whether we are on a boundary.
  if (elem2!=-1)
    {
      // orient from the centre of gravity
      // compute the centre of gravity of each element
      for(i=0; i<4; i++)
        {
          x1g+=les_coords(Face_sommets(elem_faces(elem1,i),0),0);
          x1g+=les_coords(Face_sommets(elem_faces(elem1,i),1),0);
          y1g+=les_coords(Face_sommets(elem_faces(elem1,i),0),1);
          y1g+=les_coords(Face_sommets(elem_faces(elem1,i),1),1);
          x2g+=les_coords(Face_sommets(elem_faces(elem2,i),0),0);
          x2g+=les_coords(Face_sommets(elem_faces(elem2,i),1),0);
          y2g+=les_coords(Face_sommets(elem_faces(elem2,i),0),1);
          y2g+=les_coords(Face_sommets(elem_faces(elem2,i),1),1);
        }

      grx=(x2g-x1g)*0.125;
      gry=(y2g-y1g)*0.125;

      // check the sign of the dot product
      psc=grx*nx+gry*ny;
      if(psc<0)
        {
          if(elem1<elem2)
            sign=-1;
        }
      else if(elem2<elem1)
        sign=-1;
    }
  else
    {
      // orient from the centre of gravity and the midpoint of the
      // current face

      for(i=0; i<4; i++)
        {
          x1g+=les_coords(Face_sommets(elem_faces(elem1,i),0),0);
          x1g+=les_coords(Face_sommets(elem_faces(elem1,i),1),0);
          y1g+=les_coords(Face_sommets(elem_faces(elem1,i),0),1);
          y1g+=les_coords(Face_sommets(elem_faces(elem1,i),1),1);
        }
      // Cerr << "xg et yg de Face_normales: " << x1g << " " << y1g << finl;

      x2g = les_coords(n0,0)+les_coords(n1,0);
      y2g = les_coords(n0,1)+les_coords(n1,1);
      grx=x2g*0.5-x1g*0.125;
      gry=y2g*0.5-y1g*0.125;

      //   Cerr << "grx et gry : " << grx << " " << gry << finl;
      // check the sign of the dot product
      psc=grx*nx+gry*ny;
      if(psc<0)
        sign=-1;
    }
  double scale = 1.0;
  if (bidim_axi)
    {
      const double r0 = les_coords(n0, 0);
      const double r1 = les_coords(n1, 0);
      const double r_bar = 0.5 * (r0 + r1);
      scale = 2.0 * M_PI * ((r_bar <=1e-10) ? x1g / 8.0 : r_bar);
    }
  Face_normales(num_Face, 0) = sign * nx * scale;
  Face_normales(num_Face, 1) = sign * ny * scale;
}

/*! @brief
 *
 */
void Quadri_EF::calcul_vc(const ArrOfInt& Face,ArrOfDouble& vc,
                          const ArrOfDouble& vs,const DoubleTab& vsom,
                          const Champ_Inc_base& vitesse,int type_cl) const
{
  //Cerr << " DANS Quadri_EF::calcul_vc , type_cl = " << type_cl << finl;
  //Cerr << "vs " << vs << " et vsom " << vsom << " et vitesse " << vitesse << finl;

//   switch(type_cl) {
//   case 0: //  no Dirichlet face
//     {
  vc[0] = vs[0]*0.25;
  vc[1] = vs[1]*0.25;
//       break;
//     }

//   case 1: // one Dirichlet face: Face 3
//     {
//       vc[0] = vitesse.valeurs()(Face[3],0);
//       vc[1] = vitesse.valeurs()(Face[3],1);
//       break;
//     }

//   case 3: // one Dirichlet face: Face 2
//     {
//       vc[0] = vitesse.valeurs()(Face[2],0);
//       vc[1] = vitesse.valeurs()(Face[2],1);
//       break;
//     }

//   case 9: // one Dirichlet face: Face 1
//     {
//       vc[0] = vitesse.valeurs()(Face[1],0);
//       vc[1] = vitesse.valeurs()(Face[1],1);
//       break;
//     }

//   case 27: // one Dirichlet face: Face 0
//     {
//       vc[0] = vitesse.valeurs()(Face[0],0);
//       vc[1] = vitesse.valeurs()(Face[0],1);
//       break;
//     }

//   case 4: // two Dirichlet faces: Faces 2,3
//     {
//       vc[0]= vsom(3,0);
//       vc[1]= vsom(3,1);
//       break;
//     }

//   case 28: // two Dirichlet faces: Faces 0,3
//     {
//       vc[0]= vsom(2,0);
//       vc[1]= vsom(2,1);
//       break;
//     }

//   case 12: // two Dirichlet faces: Faces 1,2
//     {
//       vc[0]= vsom(1,0);
//       vc[1]= vsom(1,1);
//       break;
//     }

//   case 36: // two Dirichlet faces: Faces 0,1
//     {
//       vc[0]= vsom(0,0);
//       vc[1]= vsom(0,1);
//       break;
//     }


//   case 10: // two Dirichlet faces: Faces 1,3
//     {
//       vc[0] = vs[0]*0.25;
//       vc[1] = vs[1]*0.25;
//       break;
//     }

//   case 30: // two Dirichlet faces: Faces 0,2
//     {
//       vc[0] = vs[0]*0.25;
//       vc[1] = vs[1]*0.25;
//       break;
//     }

//   case 13: //three Dirichlet faces: Faces 3,2,1
//     {
//       vc[0]= vitesse.valeurs()(Face[2],0);
//       vc[1]= vitesse.valeurs()(Face[2],1);
//       break;
//     }

//   case 31: //three Dirichlet faces: Faces 0,3,2
//     {
//       vc[0]= vitesse.valeurs()(Face[3],0);
//       vc[1]= vitesse.valeurs()(Face[3],1);
//       break;
//     }

//   case 37: //three Dirichlet faces: Faces 1,0,3
//     {
//       vc[0]= vitesse.valeurs()(Face[0],0);
//       vc[1]= vitesse.valeurs()(Face[0],1);
//       break;
//     }

//   case 39: //three Dirichlet faces: Faces 2,1,0
//     {
//       vc[0]= vitesse.valeurs()(Face[1],0);
//       vc[1]= vitesse.valeurs()(Face[1],1);
//       break;
//     }

//   default :
//     {
//       Cerr << "\n  type inconnu : " << type_cl ;
//       exit();
//     }

//   } // end of switch

}

/*! @brief Computes the coordinates xg of the centre of a non-standard element.
 *
 * @brief Also computes idirichlet = number of Dirichlet faces of the element.
 *  If idirichlet=2, n1 is the index of the vertex coinciding with G.
 * @param xg Output centre coordinates.
 * @param x Vertex coordinate table for the element.
 * @param type_elem_Cl Element boundary condition type.
 * @param idirichlet Output number of Dirichlet faces.
 * @param n1 Output index of the vertex coinciding with G (when idirichlet=2).
 */
void Quadri_EF::calcul_xg(DoubleVect& xg, const DoubleTab& x,
                          const int type_elem_Cl,int& idirichlet,int& n1,int& ,int& ) const
{
  int j,dim=xg.size();
//   switch(type_elem_Cl) {

//   case 0:  //  no Dirichlet face: it has 4 facets
//     //  point G is the barycenter of the element vertices
//     {
  for (j=0; j<dim; j++)
    xg[j]=(x(0,j)+x(1,j)+x(2,j)+x(3,j))*0.25;
  idirichlet=0;
//       break;
//     }

//   case 1: // one Dirichlet face: Face 3
//     // point G is the center of face 3
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(2,j)+x(3,j))*0.5;
//       idirichlet=1;
//       break;
//     }

//   case 3: // one Dirichlet face: Face 2
//     // point G is the center of face 2
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(1,j)+x(3,j))*0.5;
//       idirichlet=1;
//       break;
//     }

//   case 9: // one Dirichlet face: Face 1
//     // point G is the center of face 1
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(0,j)+x(1,j))*0.5;
//       idirichlet=1;
//       break;
//     }

//   case 27: // one Dirichlet face: Face 0
//     // point G is the center of face 0
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(0,j)+x(2,j))*0.5;
//       idirichlet=1;
//       break;
//     }

//   case 4: // two Dirichlet faces: Faces 2,3
//     // point G is the vertex common to the two Dirichlet faces
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=x(3,j);
//       idirichlet=2;
//       break;
//     }

//   case 28: // two Dirichlet faces: Faces 0,3
//     // point G is the vertex common to the two Dirichlet faces
//      {
//       for (j=0; j<dim; j++)
// 	xg[j]=x(2,j);
//       idirichlet=2;
//       break;
//     }

//   case 12: // two Dirichlet faces: Faces 1,2
//     // point G is the vertex common to the two Dirichlet faces
//      {
//       for (j=0; j<dim; j++)
// 	xg[j]=x(1,j);
//       idirichlet=2;
//       break;
//     }

//   case 36: // two Dirichlet faces: Faces 0,1
//     // point G is the vertex common to the two Dirichlet faces
//      {
//       for (j=0; j<dim; j++)
// 	xg[j]=x(0,j);
//       idirichlet=2;
//       break;
//     }

//   case 10: // two Dirichlet faces: Faces 1,3
//      // keep the same control volumes as for internal faces
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(0,j)+x(1,j)+x(2,j)+x(3,j))*0.25;
//       idirichlet=0;
//       break;
//     }

//   case 30: // two Dirichlet faces: Faces 0,2
//      // keep the same control volumes as for internal faces
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(0,j)+x(1,j)+x(2,j)+x(3,j))*0.25;
//       idirichlet=0;
//       break;
//     }

//   case 13: //three Dirichlet faces: Faces 3,2,1
//     // point G is the center of the Dirichlet face opposite to the non-Dirichlet face
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(1,j)+x(3,j))*0.5;
//       idirichlet=3;
//       break;
//     }

//   case 31: //three Dirichlet faces: Faces 0,3,2
//     // point G is the center of the Dirichlet face opposite to the non-Dirichlet face
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(2,j)+x(3,j))*0.5;
//       idirichlet=3;
//       break;
//     }

//   case 37: //three Dirichlet faces: Faces 1,0,3
//     // point G is the center of the Dirichlet face opposite to the non-Dirichlet face
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(0,j)+x(2,j))*0.5;
//       idirichlet=3;
//       break;
//     }

//   case 39: //three Dirichlet faces: Faces 2,1,0
//     // point G is the center of the Dirichlet face opposite to the non-Dirichlet face
//     {
//       for (j=0; j<dim; j++)
// 	xg[j]=(x(0,j)+x(1,j))*0.5;
//       idirichlet=3;
//       break;
//     }

//   } // end of switch

}
