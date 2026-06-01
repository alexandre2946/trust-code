/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Hexa_poly.h>
#include <Domaine.h>

Implemente_instanciable_sans_constructeur(Hexa_poly,"Hexa_poly",Elem_poly_base);

// printOn and readOn

Sortie& Hexa_poly::printOn(Sortie& s ) const
{
  return s << que_suis_je() << finl;
}

Entree& Hexa_poly::readOn(Entree& s )
{
  return s ;
}

/*! @brief KEL_(0,fa7), KEL_(1,fa7) are the local indices of the 2 faces surrounding the facet with local index fa7.
 *
 *  The local index of fa7 is the index of the vertex that carries it.
 *
 */
Hexa_poly::Hexa_poly()
{
}

/*! @brief Fills in the face_normales array in the Domaine_poly.
 *
 */

void Hexa_poly::normale(int num_Face,DoubleTab& Face_normales,
                        const  IntTab& Face_sommets,
                        const IntTab& Face_voisins,
                        const IntTab& elem_faces,
                        const Domaine& domaine_geom) const
{
  const DoubleTab& les_coords = domaine_geom.coord_sommets();

  double x1,y1,z1,x2,y2,z2;
  double nx,ny,nz;
  int f0,no4;
  int elem1;

  // we have the 4 vertices of the face
  int n0 = Face_sommets(num_Face,0);
  int n1 = Face_sommets(num_Face,1);
  int n2 = Face_sommets(num_Face,2);
  int n3 = Face_sommets(num_Face,3);

  // split the face into two triangles
  // initialise Face_normales(num_face,dimension) to 0
  // because we will sum the two face_normales of the triangles
  Face_normales(num_Face,0) =  0;
  Face_normales(num_Face,1) =  0;
  Face_normales(num_Face,2) =  0;

  // take vertices S1, S2 and S4

  x1 = les_coords(n0,0) - les_coords(n1,0);
  y1 = les_coords(n0,1) - les_coords(n1,1);
  z1 = les_coords(n0,2) - les_coords(n1,2);

  x2 = les_coords(n3,0) - les_coords(n1,0);
  y2 = les_coords(n3,1) - les_coords(n1,1);
  z2 = les_coords(n3,2) - les_coords(n1,2);

  nx = (y1*z2 - y2*z1)/2;
  ny = (-x1*z2 + x2*z1)/2;
  nz = (x1*y2 - x2*y1)/2;

  // Orient the normal from elem1 to elem2
  // by searching for the vertex of elem1 that is not on the Face

  elem1=Face_voisins(num_Face,0);

  if ( (f0 = elem_faces(elem1,0)) == num_Face )
    f0 = elem_faces(elem1,1);

  if ( (no4 = Face_sommets(f0,0)) != n0    &&   no4 != n1
       &&   no4 != n2 &&   no4 != n3)
    { /* Do nothing */}
  else if ( (no4 = Face_sommets(f0,1)) != n0 && no4 != n1
            && no4 != n2 &&   no4 != n3)
    { /* Do nothing */}
  else if ( (no4 = Face_sommets(f0,2)) != n0 && no4 != n1
            && no4 != n2 &&   no4 != n3)
    { /* Do nothing */}
  else
    no4 = Face_sommets(f0,3);

  x1 = les_coords(no4,0) - les_coords(n0,0);
  y1 = les_coords(no4,1) - les_coords(n0,1);
  z1 = les_coords(no4,2) - les_coords(n0,2);

  if ( (nx*x1+ny*y1+nz*z1) > 0 )
    {
      Face_normales(num_Face,0) +=  -nx;
      Face_normales(num_Face,1) +=  -ny;
      Face_normales(num_Face,2) +=  -nz;
    }
  else
    {
      Face_normales(num_Face,0) += nx;
      Face_normales(num_Face,1) += ny;
      Face_normales(num_Face,2) += nz;
    }

  // take vertices S1, S4 and S3

  x1 = les_coords(n0,0) - les_coords(n2,0);
  y1 = les_coords(n0,1) - les_coords(n2,1);
  z1 = les_coords(n0,2) - les_coords(n2,2);

  x2 = les_coords(n3,0) - les_coords(n2,0);
  y2 = les_coords(n3,1) - les_coords(n2,1);
  z2 = les_coords(n3,2) - les_coords(n2,2);

  nx = (y1*z2 - y2*z1)/2;
  ny = (-x1*z2 + x2*z1)/2;
  nz = (x1*y2 - x2*y1)/2;

  // Orient the normal from elem1 to elem2
  // by searching for the vertex of elem1 that is not on the Face

  elem1=Face_voisins(num_Face,0);

  if ( (f0 = elem_faces(elem1,0)) == num_Face )
    f0 = elem_faces(elem1,1);

  if ( (no4 = Face_sommets(f0,0)) != n0    &&   no4 != n1
       &&   no4 != n2 &&   no4 != n3)
    { /* Do nothing */}
  else if ( (no4 = Face_sommets(f0,1)) != n0 && no4 != n1
            && no4 != n2 &&   no4 != n3)
    { /* Do nothing */}
  else if ( (no4 = Face_sommets(f0,2)) != n0 && no4 != n1
            && no4 != n2 &&   no4 != n3)
    { /* Do nothing */}
  else
    no4 = Face_sommets(f0,3);

  x1 = les_coords(no4,0) - les_coords(n0,0);
  y1 = les_coords(no4,1) - les_coords(n0,1);
  z1 = les_coords(no4,2) - les_coords(n0,2);

  if ( (nx*x1+ny*y1+nz*z1) > 0 )
    {
      Face_normales(num_Face,0) +=  -nx;
      Face_normales(num_Face,1) +=  -ny;
      Face_normales(num_Face,2) +=  -nz;
    }
  else
    {
      Face_normales(num_Face,0) += nx;
      Face_normales(num_Face,1) += ny;
      Face_normales(num_Face,2) += nz;
    }
}



/*! @brief
 *
 */
void Hexa_poly::calcul_vc(const ArrOfInt& Face,ArrOfDouble& vc,
                          const ArrOfDouble& vs,const DoubleTab& vsom,
                          const Champ_Inc_base& vitesse,int type_cl) const
{

  vc[0] = vs[0]/6;
  vc[1] = vs[1]/6;
  vc[2] = vs[2]/6;
}

/*! @brief Computes the coordinates xg of the centre of a non-standard element;
 *   also computes idirichlet = number of Dirichlet faces of the element.
 *
 *  If idirichlet=2, n1 is the index of the vertex coinciding with G.
 *
 */
void Hexa_poly::calcul_xg(DoubleVect& xg, const DoubleTab& x,
                          const int type_elem_Cl,int& idirichlet,int& n1,int& ,int& ) const
{
  int j,dim=xg.size();
  for (j=0; j<dim; j++)
    xg[j]=(x(0,j)+x(1,j)+x(2,j)+x(3,j)+x(4,j)+x(5,j)+x(6,j)+x(7,j))*0.125;
  idirichlet=0;
}


