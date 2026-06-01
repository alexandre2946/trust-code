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

#include <Faces.h>
#include <DomaineAxi1d.h>
#include <communications.h>
#include <Linear_algebra_tools_impl.h>
#include <TRUSTLists.h>

Implemente_instanciable_32_64(Faces_32_64,"Faces",Objet_U);

/*! @brief Writes the faces to an output stream.
 *
 * The following are written:
 *       - the face type
 *       - the vertices
 *       - the neighbors
 *     Only the type Type_Face::vide_0D is written if the number of faces is zero.
 *
 * @param s An output stream.
 * @return The modified output stream.
 */
template <typename _SIZE_>
Sortie& Faces_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  if(nb_faces()==0)
    s << "vide_0D" << finl;
  else
    {
      s << (type(type_face_)) << finl;
      s <<  sommets;
      s << faces_voisins;
    }
  return s ;
}

/*! @brief Writes the faces to an output stream (binary format).
 *
 * The following are written:
 *       - the face type
 *       - the vertices
 *       - the neighbors
 *     Only the type Type_Face::vide_0D is written if the number of faces is zero.
 *
 * @param s An output stream.
 * @return The modified output stream.
 */
template <typename _SIZE_>
Sortie& Faces_32_64<_SIZE_>::ecrit(Sortie& s ) const
{
  if(nb_faces()==0)
    s << "vide_0D" << finl;
  else
    {
      s << (type(type_face_)) << finl;
      sommets.ecrit(s);
      faces_voisins.ecrit(s);
    }
  return s ;
}

/*! @brief Reads the specifications of a face object from an input stream.
 *
 *     The following are read:
 *       - the face type
 *       - the vertices
 *       - the neighbors
 *     Creates a Faces_32_64 object with 0 faces if the type is "vide_0D".
 *
 * @param s An input stream.
 * @return The modified input stream.
 */
template <typename _SIZE_>
Entree& Faces_32_64<_SIZE_>::readOn(Entree& s)
{
  Motcle typ;
  s >> typ;
  typer(typ);
  if (typ != "vide_0D")
    {
      // Vertices and neighbors are only written if the type is not empty:
      s >> sommets;
      s >> faces_voisins;
    }
  else
    {
      // to avoid a crash in the subsequent code that assumes nb_dim()==2 and line_size_ > 0
      sommets.resize(0,1);
      faces_voisins.resize(0,2);
    }
  return s;
}


/*! @brief Reads the specifications of a face object from an input stream (binary format).
 *
 *     The following are read:
 *       - the face type
 *       - the vertices
 *       - the neighbors
 *     Vertices are reordered after reading. Does nothing if the type read is "vide_0D".
 *
 * @param s An input stream.
 * @return The modified input stream.
 */
template <typename _SIZE_>
Entree& Faces_32_64<_SIZE_>::lit(Entree& s)
{
  Motcle typ;
  s >> typ;
  if(typ!="vide_0D")
    {
      typer(typ);
      sommets.lit(s);
      faces_voisins.lit(s);
      reordonner();
    }
  else
    {
      typer(typ);
      sommets.resize(0,0);
      faces_voisins.resize(0,2);
    }
  return s ;
}


/*! @brief Returns the index of the smallest vertex (in terms of vertex numbering) of a face.
 *
 * @param (int face) the index of the face for which the smallest vertex is sought
 * @param (int nb_som) the number of vertices per face
 * @return (int) the index of the smallest vertex of the specified face
 */
template <typename _SIZE_>
typename Faces_32_64<_SIZE_>::int_t Faces_32_64<_SIZE_>::ppsf(int_t face, int nb_som) const
{
  int_t som = sommet(face,0);
  for(int imin=1; imin < nb_som; imin++)
    som=std::min(som,sommet(face,imin));
  return som;
}

/*! @brief Compares 2 faces from 2 face sets f1 and f2 and returns true if they are equal, false otherwise.
 *
 * @param (int f1) the index of the face in the first face set represented by *this
 * @param (Faces& faces2) the second face set
 * @param (int f2) the index of the face in the second face set
 * @param (int nb_som) the number of vertices per face
 * @return (bool) true if the 2 faces are the same, false otherwise
 */
template <typename _SIZE_>
bool Faces_32_64<_SIZE_>::same_face(int_t f1, const Faces_32_64& faces2, int_t f2, int nb_som) const
{
  // [ABN] ??? I don't understand this ugliness - should compare std::set<>?
  bool ok=true;
  for(int i=0; i<nb_som && ok ; i++)
    {
      bool nok=true;
      for(int j=0; j<nb_som && nok ; j++)
        {
          if(sommet(f1, i) == faces2.sommet(f2, j))
            nok=false;           // same vertices
        }
      if(nok) ok=false;          // vertex not found => not the same face
    }
  return ok;
}

/*! @brief Returns a Type_Face object representing the type with the specified name.
 *
 *     Recognized type names are:
 *              "vide_0D"
 *       "point_1D"
 *       "point_1D_axi"
 *       "segment_2D"
 *       "triangle_3D"
 *       "quadrangle_3D"
 *       "segment_2D_axi"
 *       "quadrangle_3D_axi"
 *       "quadrilatere_2D_axi"
 *
 * @param (Motcle& mot) the keyword representing a face type
 * @return (Type_Face) the face type corresponding to the parameter
 * @throws unrecognized face type
 */
template <typename _SIZE_>
Type_Face Faces_32_64<_SIZE_>::type(const Motcle& mot) const
{
  Motcles les_mots(10);
  {
    les_mots[0]="vide_0D";
    les_mots[1]="point_1D";
    les_mots[2]="point_1D_axi";
    les_mots[3]="segment_2D";
    les_mots[4]="triangle_3D";
    les_mots[5]="quadrangle_3D";
    les_mots[6]="segment_2D_axi";
    les_mots[7]="quadrangle_3D_axi";
    les_mots[8]="quadrilatere_2D_axi";
    les_mots[9]="polygone_3d";
  }
  int rang=les_mots.search(mot);
  switch(rang)
    {
    case  0 :
      return Type_Face::vide_0D;
    case  1 :
      return Type_Face::point_1D;
    case  2 :
      return Type_Face::point_1D_axi;
    case  3 :
      return Type_Face::segment_2D;
    case  4 :
      return Type_Face::triangle_3D;
    case  5 :
      return Type_Face::quadrangle_3D;
    case  6 :
      return Type_Face::segment_2D_axi;
    case  7 :
      return Type_Face::quadrangle_3D_axi;
    case  8 :
      return Type_Face::quadrilatere_2D_axi;
    case 9 :
      return Type_Face::polygone_3D;
    default :
      {
        Cerr << "In the area number " << Process::me() << " " << mot << " is not a type of face." << finl;
        Cerr << "Check your splitting." << finl;
        // If mot is empty, an attempt is being made to re-read Domains created with a version prior to 1.5.1
        if (mot=="")
          {
            Cerr << "Your splitting seems to have been done with a TRUST version 1.5 or earlier. Since" << finl;
            Cerr << "the 1.5.1, the files format of .Zones containing the splitting of your mesh has evolved." << finl;
            Cerr << "You must therefore rebuild the splitting of your mesh with TRUST version 1.5.1 or newer." <<
                 finl;
          }
        exit();
      }
    }
  // for the compiler:
  return Type_Face::point_1D;
}

/*! @brief Returns the name associated with a face type.
 *
 * (inverse of Type_Face Faces_32_64<_SIZE_>::type(const Motcle& ) const)
 *
 * @param (Type_Face& typ) a face type
 * @return (Motcle&) the name corresponding to the specified type
 * @throws unrecognized face type
 */
template <typename _SIZE_>
Motcle& Faces_32_64<_SIZE_>::type(const Type_Face& typ) const
{
  static Motcle mot;
  switch(typ)
    {
    case  Type_Face::vide_0D :
      mot="vide_0D";
      break;
    case  Type_Face::point_1D :
      mot="point_1D";
      break;
    case  Type_Face::point_1D_axi :
      mot="point_1D_axi";
      break;
    case  Type_Face::segment_2D :
      mot="segment_2D";
      break;
    case  Type_Face::segment_2D_axi :
      mot="segment_2D_axi";
      break;
    case  Type_Face::triangle_3D :
      mot="triangle_3D";
      break;
    case  Type_Face::quadrangle_3D :
      mot="quadrangle_3D";
      break;
    case  Type_Face::quadrangle_3D_axi :
      mot="quadrangle_3D_axi";
      break;
    case Type_Face::quadrilatere_2D_axi :
      mot="quadrilatere_2D_axi";
      break;
    case Type_Face::polygone_3D:
      mot="polygone_3D";
      break;
    default :
      {
        Cerr << "Error TRUST in Motcle& Faces_32_64<_SIZE_>::type(const Type_Face& typ)" << finl;
        exit();
      }
    }
  return mot;
}

/*! @brief Adds faces.
 *
 * Adding faces is equivalent to adding vertices.
 *
 * @param (IntTab& sommet) the array of vertices to add
 * @throws the specified vertices do not have the correct dimension (1D,2D,3D)
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::ajouter(const IntTab_t& tab_sommet)
{
  assert(sommets.dimension(1)==tab_sommet.dimension(1));
  int_t oldsz=sommets.dimension(0);
  int_t addsz=tab_sommet.dimension(0);
  dimensionner(oldsz+addsz);
  for(int_t i=0; i<addsz; i++)
    for(int j=0; j<tab_sommet.dimension(1); j++)
      sommets(oldsz+i,j)=tab_sommet(i,j);
}

/*! @brief (Re-)sizes the faces. The neighbors are resized accordingly.
 *
 *     Implicitly added vertices are initialized to -1.
 *
 * @param (int i) the new number of faces
 * @return (int) the new number of faces
 */
template <typename _SIZE_>
typename Faces_32_64<_SIZE_>::int_t Faces_32_64<_SIZE_>::dimensionner(int_t i)
{
  int_t oldsz=nb_faces();
  int nombre_som_faces=nb_som_faces();
  if (sommets.size()!=0) sommets.detach_vect(); // Fixed bug: To have the possibility of merging borders after discretization
  sommets.resize(i,nombre_som_faces);
  faces_voisins.resize(i, 2);
  for(int_t j=oldsz; j<i; j++)
    {
      for(int k=0; k<nombre_som_faces; k++)
        sommets(j,k)=-1;
      faces_voisins(j,0)=faces_voisins(j,1)=-1;
    }
  return i;
}

/*! @brief Initializes the neighbors of joint faces to -1.
 *
 * @param (int nb_faces_joint) the number of faces representing Joints
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::initialiser_faces_joint(int_t nb_faces_joint)
{
  for(int_t i=0; i<nb_faces_joint; i++)
    //      faces_voisins(i,0)=faces_voisins(i,1)=-2;
    faces_voisins(i,0)=faces_voisins(i,1)=-1;
}

/*! @brief Initializes the vertices of joint faces to -1.
 *
 * @param (int nb_faces_joint) the number of faces representing Joints
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::initialiser_sommets_faces_joint(int_t nb_faces_joint)
{
  for(int_t i=0; i<nb_faces_joint; i++)
    for(int k=0; k<nb_som_faces(); k++)
      //         sommets(i,k)=-2;
      sommets(i,k)=-1;
}

/*! @brief Sets the type of the faces.
 *
 * @param (Motcle& typ) the type to assign to the faces
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::typer(const Motcle& typ)
{
  typer(type(typ));
}

/*! @brief Sets the type of the faces.
 *
 * @param (Type_Face& typ) the type to assign to the faces
 * @throws unknown face type
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::typer(const Type_Face& typ)
{
  type_face_ = typ;
  // int axi_ = 0;
  int max_dim = 3;
  switch(typ)
    {
    case  Type_Face::vide_0D :
      nb_som_face=0;
      break;
    case  Type_Face::point_1D :
    case  Type_Face::point_1D_axi :
      nb_som_face=1;
      break;
    case  Type_Face::segment_2D :
      nb_som_face=2;
      break;
    case  Type_Face::segment_2D_axi :
      nb_som_face=2; //axi_=1;
      break;
    case  Type_Face::triangle_3D :
      nb_som_face=3;
      break;
    case  Type_Face::quadrangle_3D :
      nb_som_face=4;
      break;
    case  Type_Face::quadrangle_3D_axi :
      nb_som_face=4; //axi_=1;
      break;
    case  Type_Face::quadrilatere_2D_axi :
      nb_som_face=2; //axi_=1;
      break;
    case Type_Face::polygone_3D:
      nb_som_face=-1; // will be set later
      break;
    default :
      {
        Cerr << "Error TRUST in Faces_32_64<_SIZE_>::typer(const Motcle& typ)" << finl;
        exit();
      }
    }
  /* TroisDto2D marche pas
     if (axi_==1 && bidim_axi!=1 && axi!=1)
     {
     Cerr << "The mesh contains faces that show that is a revolution mesh." << finl;
     Cerr << "Add Bidim_axi or Axi in your data file." << finl;
     exit();
     } */
  if (dimension<std::min(max_dim,nb_som_face))
    {
      Cerr << "You are in dimension " << dimension << finl;
      Cerr << "and the mesh contains faces with " << nb_som_face << " nodes." << finl;
      Cerr << "and this seems to be a mesh of dimension " << std::min(max_dim,nb_som_face) << " ..." << finl;
      exit();
    }
}


/*! @brief Completes the specified face: updates its neighbors.
 *
 * @param (int face) the index of the face to complete
 * @param (int num_elem) the index of the neighboring element of the face
 * @throws face already complete
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::completer(int_t face, int_t num_elem)
{
  if ( voisin(face,0) == -1)
    voisin(face, 0) = num_elem;
  else if( voisin(face,1) == -1)
    voisin(face, 1) = num_elem;
  else
    {
      Cerr << finl;
      Cerr << "Problem for the face number " << face << " and the cell " << num_elem << finl;
      Cerr << "Indeed, the face already belongs to the cells " << voisin(face, 0) << " and " << voisin(face, 1) << finl;
      Cerr << "The nodes of this face are:" << finl;
      for (int i=0; i<nb_som_face; i++)
        Cerr << "Node " << sommet(face,i) << finl;
      Cerr << "Check your mesh. Perhaps some cells " << finl;
      Cerr << "are defined 2 times." << finl;
      exit();
    }
}

/*! @brief Computes the surface area of the faces.
 *
 * @param (DoubleVect& surfaces) the vector containing the surface area of each face.
 * @throws unrecognized face type
 * @throws surface area computation not implemented for this face type
 * @throws erroneous surface computation (surface <= 0)
 * @throws face type not consistent with the space dimension
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::calculer_surfaces(DoubleVect_t& surfaces) const
{
  surfaces.resize(nb_faces_tot());
  const Domaine_t& dom=domaine();
  // Check that in cylindrical coordinates, r is positive (with numerical tolerance)
  if (axi || bidim_axi)
    {
      const int_t nb_som = dom.les_sommets().dimension(0);
      // Radius scale to set a robust tolerance
      double rmax = 0.;
      for (int_t i = 0; i < nb_som; i++)
        rmax = std::max(rmax, std::fabs(dom.coord(i,0)));
      const double tol = std::max(Objet_U::precision_geom * std::max(1.0, rmax), 1e-300);
      for (int_t i = 0; i < nb_som; i++)
        {
          const double r = dom.coord(i,0);
          if (r < -tol)
            {
              Cerr << "In axisymmetric, the coordinates of the mesh according to radius" << finl;
              Cerr << "ie along X, cannot be negative beyond tolerance." << finl;
              Cerr << "The revolution axis must be at x=0." << finl;
              Cerr << "we got x = " << r << "  (tolerance = " << tol << ")" << finl;
              exit();
            }
        }
    }
  switch(type_face_)
    {
    case Type_Face::segment_2D :
      {
        assert(dimension==2);
        for (int_t face = 0; face < nb_faces_tot(); face++)
          {
            const double x0 = dom.coord(sommet(face,0), 0);
            const double y0 = dom.coord(sommet(face,0), 1);
            const double x1 = dom.coord(sommet(face,1), 0);
            const double y1 = dom.coord(sommet(face,1), 1);
            const double dx = x0 - x1;
            const double dy = y0 - y1;
            const double L  = std::sqrt(dx*dx + dy*dy);

            if (!Objet_U::bidim_axi)
              {
                // pure 2D Cartesian: surface = length
                surfaces(face) = L;
                if(surfaces(face) == 0.)
                  {
                    Cerr << "area("<<face<<")=0 ! Check your mesh." << finl;
                    exit();
                  }
              }
            else
              {
                // RZ (bidim_axi): surface of the torus generated by the segment
                // S = Δθ * r_bar * L, with r_bar = (r0 + r1)/2 and r ≡ x
                const double r0 = x0;
                const double r1 = x1;
                const double rbar = 0.5 * (r0 + r1);
                surfaces(face) = 2.0 * M_PI * L * rbar;
              }
          }
        break;
      }
    case Type_Face::quadrilatere_2D_axi :
      {
        assert(dimension==2);
        for (int_t face = 0; face < nb_faces_tot(); face++)
          {
            const double r0 = dom.coord(sommet(face,0), 0);
            const double z0 = dom.coord(sommet(face,0), 1);
            const double r1 = dom.coord(sommet(face,1), 0);
            const double z1 = dom.coord(sommet(face,1), 1);
            const double dr = r1 - r0;
            const double dz = z1 - z0;
            const double L  = std::sqrt(dr*dr + dz*dz);       // segment length
            const double rbar = 0.5*(r0 + r1);                 // linear average of r(s) over the segment
            surfaces(face) = 2.0 * M_PI * rbar * L; // S = Δθ ∫_Γ r ds = Δθ * r̄ * L
          }
        break;
      }

    case  Type_Face::segment_2D_axi :
      {
        assert(dimension==2);
        double r0,r1,teta0,teta1,d_teta;
        for(int_t face=0; face <nb_faces_tot(); face++)
          {
            r0 = dom.coord(sommet(face ,0), 0);
            r1 = dom.coord(sommet(face ,1), 0);
            if ( est_egal(r0,r1) ) // surface = r0*abs(teta1-teta0)
              {
                teta0 = dom.coord(sommet(face ,0), 1);
                teta1 = dom.coord(sommet(face ,1), 1);
                d_teta = std::fabs(teta1-teta0);
                if(d_teta > M_PI) d_teta=2.*M_PI-d_teta;
                surfaces(face)=r0*d_teta;
              }
            else //  surface = r1-r2
              surfaces(face)=std::fabs(r1-r0);
          }
        break;
      }
    case  Type_Face::triangle_3D :
      {
        assert(dimension==3);

        double delta0, delta1, delta2;
        double longueur0, longueur1;
        double prod,sa;

        for(int_t face=0; face <nb_faces_tot(); face++)
          {
            prod=0;
            delta0=(dom.coord(sommet(face ,1), 0) - dom.coord(sommet(face ,0), 0));
            delta1=(dom.coord(sommet(face ,1), 1) - dom.coord(sommet(face ,0), 1));
            delta2=(dom.coord(sommet(face ,1), 2) - dom.coord(sommet(face ,0), 2));
            longueur0=(delta0*delta0+delta1*delta1+delta2*delta2);
            sa=delta0;
            delta0=(dom.coord(sommet(face ,2), 0) - dom.coord(sommet(face ,0), 0));
            prod=sa*delta0;
            sa=delta1;
            delta1=(dom.coord(sommet(face ,2), 1) - dom.coord(sommet(face ,0), 1));
            prod+=sa*delta1;
            sa=delta2;
            delta2=(dom.coord(sommet(face ,2), 2) - dom.coord(sommet(face ,0), 2));
            prod+=sa*delta2;
            longueur1=(delta0*delta0+delta1*delta1+delta2*delta2);
            surfaces(face)=0.5*(sqrt(longueur0*longueur1-prod*prod));
            if(surfaces(face)==0.)
              {
                Cerr << "area("<<face<<")=0 ! Check your mesh." << finl;
                exit();
              }
          }
        break;
      }
    case  Type_Face::quadrangle_3D :
      {
        // Based on Hexa_VEF::normale():
        for(int_t face=0; face <nb_faces_tot(); face++)
          {
            int_t n0 = sommet(face, 0),
                  n1 = sommet(face, 1),
                  n2 = sommet(face, 2),
                  n3 = sommet(face, 3);
            // NB: In a prism, there are also triangles as faces... Therefore:
            if (n3<0) n3 = n2;

            double x1 = dom.coord(n0, 0) - dom.coord(n1, 0);
            double y1 = dom.coord(n0, 1) - dom.coord(n1, 1);
            double z1 = dom.coord(n0, 2) - dom.coord(n1, 2);

            double x2 = dom.coord(n3, 0) - dom.coord(n1, 0);
            double y2 = dom.coord(n3, 1) - dom.coord(n1, 1);
            double z2 = dom.coord(n3, 2) - dom.coord(n1, 2);

            double nx = (y1*z2 - y2*z1)/2;
            double ny = (-x1*z2 + x2*z1)/2;
            double nz = (x1*y2 - x2*y1)/2;

            x1 = dom.coord(n0,0) - dom.coord(n2,0);
            y1 = dom.coord(n0,1) - dom.coord(n2,1);
            z1 = dom.coord(n0,2) - dom.coord(n2,2);

            x2 = dom.coord(n3,0) - dom.coord(n2,0);
            y2 = dom.coord(n3,1) - dom.coord(n2,1);
            z2 = dom.coord(n3,2) - dom.coord(n2,2);

            nx -= (y1*z2 - y2*z1)/2;
            ny -= (-x1*z2 + x2*z1)/2;
            nz -= (x1*y2 - x2*y1)/2;

            surfaces(face)=sqrt(nx*nx+ny*ny+nz*nz);
          }
        break;
      }
    case Type_Face::quadrangle_3D_axi :
      {
        assert(dimension==3);
        double r0,r1,teta0,teta1,teta2,z0,z2,d_teta,delta_r;
        for(int_t face=0; face <nb_faces_tot(); face++)
          {
            r0 = dom.coord(sommet(face ,0), 0);
            r1 = dom.coord(sommet(face ,1), 0);
            delta_r = std::fabs(r1 - r0);
            if ( est_egal(r0,r1) )
              {
                teta0 = dom.coord(sommet(face ,0), 1);
                teta1 = dom.coord(sommet(face ,1), 1);
                z0 = dom.coord(sommet(face ,0), 2);
                z2 = dom.coord(sommet(face ,2), 2);
                d_teta = std::fabs(teta1-teta0);
                if(d_teta > M_PI) d_teta=2.*M_PI-d_teta;
                surfaces(face)=r0*d_teta*std::fabs(z2-z0);
              }
            else
              {
                teta0 = dom.coord(sommet(face ,0), 1);
                teta2 = dom.coord(sommet(face ,2), 1);
                if (teta0 == teta2)
                  {
                    z0 = dom.coord(sommet(face ,0), 2);
                    z2 = dom.coord(sommet(face ,2), 2);
                    surfaces(face) = delta_r*std::fabs(z2-z0);
                  }
                else
                  {
                    d_teta=std::fabs(teta2-teta0);
                    if(d_teta > M_PI) d_teta=2.*M_PI-d_teta;
                    surfaces(face) = 0.5*(r0+r1)*d_teta*delta_r;
                  }
              }
          }
        break;
      }
    case Type_Face::point_1D:
    case Type_Face::vide_0D :
      {
        for(int_t face=0; face <nb_faces_tot(); face++)
          surfaces(face) = 1.0;

        break;
      }
    case Type_Face::point_1D_axi:
      {
        assert(dimension==3);

        const DomaineAxi1d_t& domax = ref_cast(DomaineAxi1d_t,dom);

        for(int_t face=0; face <nb_faces_tot(); face++)
          {
            int_t elem = voisin(face,0)==-1 ? voisin(face,1) : voisin(face,0);

            double x0 = domax.origine_repere(elem,0);
            double y0 = domax.origine_repere(elem,1);
            double x = dom.coord(sommet(face ,0), 0);
            double y = dom.coord(sommet(face ,0), 1);

            double r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
            surfaces(face) = 2.*M_PI*r;
          }
        break;
      }
    case Type_Face::polygone_3D:
      {
        const DoubleTab_t& coord=dom.coord_sommets();
        int nmax=les_sommets().dimension_int(1);
        for(int_t face=0; face <nb_faces_tot(); face++)
          {
            double n0=0,n1=0,n2=0;
            int n=nmax-1;
            while (n >= 0 && sommet(face,n)==-1) n--;
            for (int i0=0; i0<=n; i0++)
              {

                int ip1_0=0;
                if (i0<n) ip1_0=i0+1;
                int_t i=sommet(face,i0);
                int_t ip1=sommet(face,ip1_0);
                n0+=coord(i,1)*coord(ip1,2)-coord(i,2)*coord(ip1,1);
                n1+=coord(i,2)*coord(ip1,0)-coord(i,0)*coord(ip1,2);
                n2+=coord(i,0)*coord(ip1,1)-coord(i,1)*coord(ip1,0);
              }
            surfaces(face)=sqrt(n0*n0+n1*n1+n2*n2)/2.;
          }
        break;
      }
    default :
      {
        Cerr << "Error TRUST in type of Faces_32_64 not recognized " << finl;
        exit();
      }
    }
}


/*! @brief Computes the centers of gravity of each face.
 *
 * @param (DoubleTab& xv) array containing the coordinates of the centers of gravity of each face. xv(i,j) contains the j-th coordinate of the center of gravity of the i-th face. The array xv is resized and assigned the parallel descriptor of the vertex array before being filled.
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::calculer_centres_gravite(DoubleTab_t& xv) const
{
  // The array xv is sized in ::calculer_centres_gravite
  const Domaine_t& dom=domaine();
  const DoubleTab_t& coord=dom.coord_sommets();
  Faces_32_64::Calculer_centres_gravite(xv, type_face_, coord, sommets);
}


template <typename _SIZE_>
void Faces_32_64<_SIZE_>::Calculer_centres_gravite(DoubleTab_t& xv, Type_Face type_face_, const DoubleTab_t& coord,  const IntTab_t& sommet)
{
  int_t nb_faces_tot = sommet.dimension_tot(0);
  int dim = coord.dimension_int(1);
  xv.resize(nb_faces_tot, dim);
  // Copy the parallel descriptor of the vertex array:
  xv.set_md_vector(sommet.get_md_vector());

  if(nb_faces_tot!=0)
    {
      int nb_som_faces=sommet.dimension_int(1);
      double inv=1./nb_som_faces;
      xv = 0;
      if(Objet_U::axi)
        {
          switch(type_face_)
            {
            case  Type_Face::segment_2D_axi :
              {
                for (int_t face=0; face<nb_faces_tot; face++)
                  {
                    for(int som=0; som < nb_som_faces; som++)
                      xv(face, 0)+=coord(sommet(face ,som), 0);
                    double teta_0=coord(sommet(face ,0), 1);
                    double teta_1=coord(sommet(face ,1), 1);
                    double teta_min=std::min(teta_1, teta_0);
                    double teta_max=std::max(teta_1, teta_0);
                    double d_teta=teta_max-teta_min;
                    if( (teta_min==0.) && (d_teta>M_PI) )
                      {
                        teta_min+=2.*M_PI;
                      }
                    {
                      xv(face, 1)+=(teta_min+teta_max);
                    }
                  }
                break;
              }
            case Type_Face::quadrangle_3D_axi :
              {
                for (int_t face=0; face<nb_faces_tot; face++)
                  {
                    for(int som=0; som < nb_som_faces; som++)
                      xv(face, 0)+=coord(sommet(face ,som), 0);
                    for(int som=0; som < nb_som_faces; som++)
                      xv(face, 2)+=coord(sommet(face ,som), 2);
                    double teta_0=coord(sommet(face ,0), 1);
                    double teta_1=coord(sommet(face ,1), 1);
                    double teta_2=coord(sommet(face ,2), 1);
                    double teta_3=coord(sommet(face ,3), 1);
                    double teta_min=std::min(teta_1, teta_0);
                    teta_min=std::min(teta_min, teta_2);
                    teta_min=std::min(teta_min, teta_3);
                    double teta_max=std::max(teta_1, teta_0);
                    teta_max=std::max(teta_min, teta_2);
                    teta_max=std::max(teta_min, teta_3);
                    double d_teta=teta_max-teta_min;
                    if( (teta_min==0.) && (d_teta>M_PI) )
                      {
                        if(teta_0==0.) teta_0+=2.*M_PI;
                        if(teta_1==0.) teta_1+=2.*M_PI;
                        if(teta_2==0.) teta_2+=2.*M_PI;
                        if(teta_3==0.) teta_3+=2.*M_PI;
                      }
                    {
                      xv(face, 1)+=(teta_0+teta_1+teta_2+teta_3) ;
                    }
                  }
                break;
              }
            default:
              {
                Cerr << "Face type number " << (int)type_face_ << " not provided in Faces_32_64<_SIZE_>::calculer_centres_gravite" << finl;
                break;
              }
            }
          xv*=inv;
        }
      else
        {
          for (int_t face=0; face<nb_faces_tot; face++)
            {
              int nb_som=nb_som_faces;
              while (sommet(face,nb_som-1)==-1) nb_som--;

              inv=1./nb_som;
              for(int k=0; k<dim; k++)
                for(int som=0; som < nb_som; som++)
                  xv(face, k)+=coord(sommet(face ,som), k);
              for(int k=0; k<dim; k++)
                xv(face, k)*=inv;
            }
        }
    }
}

/*! @brief Reorders the faces.
 *
 * (only quadrangles are reordered)
 *
 * @throws unrecognized face type
 */
template <typename _SIZE_>
void Faces_32_64<_SIZE_>::reordonner()
{
  //Cerr << "Faces_32_64<_SIZE_>::reordonner()" << finl;
  switch(type_face_)
    {
    case  Type_Face::point_1D :
    case  Type_Face::point_1D_axi :
      {
        break;
      }
    case  Type_Face::segment_2D :
      {
        // one can have boundary faces from a triangular mesh in 3D
        assert(dimension>=2);
        break;
      }
    case Type_Face::quadrilatere_2D_axi :
      {
        assert(dimension==2);
        const Domaine_t& dom=domaine();
        SmallArrOfTID_t S(2);
        SmallArrOfTID_t NS(2);
        int i;
        const int_t nombre_faces=nb_faces();
        for(int_t face=0; face < nombre_faces; face++)
          {
            NS=-1;
            for(i=0; i<2; i++)
              S[i] = sommet(face, i);
            assert( S[0] >=0 );
            assert( S[1] >=0 );
            if (dom.coord(S[0], 0) > dom.coord(S[1], 0))
              {
                NS[0]=S[1];
                NS[1]=S[0];
              }
            else
              {
                NS[0]=S[0];
                NS[1]=S[1];
              }

            for(i=0; i<2; i++)
              {
                assert(NS[i]!=-1);
                sommet(face, i)=NS[i];
              }
          }
        break;
      }
    case  Type_Face::segment_2D_axi :
      {
        assert(dimension==2);
        break;
      }
    case  Type_Face::triangle_3D :
      {
        assert(dimension==3);
        break;
      }
    case  Type_Face::quadrangle_3D :
      {
        assert(dimension==3);
        const Domaine_t& dom=domaine();
        SmallArrOfTID_t S(4);
        SmallArrOfTID_t NS(4);
        int i;
        double xmin, ymin, zmin;
        double xmax, ymax, zmax;
        const int_t nombre_faces=nb_faces();
        for(int_t face=0; face < nombre_faces; face++)
          {
            NS=-1;
            for(i=0; i<4; i++)
              {
                S[i] = sommet(face, i);
                assert( S[i] >=0 );
              }
            if (1)
              {
                // adapted from Hexaedre_VEF::reordonne
                const DoubleTab_t& coord=dom.les_sommets();
                int_t s0=S[0], s1=S[1], s2=S[2], s3=S[3];
                double x03=coord(s3,0)-coord(s0,0);
                double y03=coord(s3,1)-coord(s0,1);
                double z03=coord(s3,2)-coord(s0,2);
                double x02=coord(s2,0)-coord(s0,0);
                double y02=coord(s2,1)-coord(s0,1);
                double z02=coord(s2,2)-coord(s0,2);

                double x01=coord(s1,0)-coord(s0,0);
                double y01=coord(s1,1)-coord(s0,1);
                double z01=coord(s1,2)-coord(s0,2);

                Vecteur3 OA(x01,y01,z01);
                Vecteur3 OB(x03,y03,z03);
                Vecteur3 OC(x02,y02,z02);

                Vecteur3 nplan, nmedian;
                // rectangle C B
                //           O A
                // compute n = OA ^ OC
                Vecteur3::produit_vectoriel(OA,OC,nplan);
                // compute the normal nmedian of the plane (OB,n)
                Vecteur3::produit_vectoriel(OB,nplan,nmedian);

                // check that points A and C are on opposite sides of the plane
                double psC=Vecteur3::produit_scalaire(nmedian,OA);
                double psA=Vecteur3::produit_scalaire(nmedian,OC);

                if (psA*psC>0)
                  {
                    // swap vertices 2 and 3
                    sommet(face,2)=S[3];
                    sommet(face,3)=S[2];
                    Cerr << "Permutation of local nodes 2 and 3 on the face " <<face<<finl;

                  }
              }
            else
              {

                xmin=std::min(dom.coord(S[0], 0), dom.coord(S[1], 0));
                xmin=std::min(xmin, dom.coord(S[2], 0));
                xmax=std::max(dom.coord(S[0], 0), dom.coord(S[1], 0));
                xmax=std::max(xmax, dom.coord(S[2], 0));
                ymin=std::min(dom.coord(S[0], 1), dom.coord(S[1], 1));
                ymin=std::min(ymin, dom.coord(S[2], 1));
                ymax=std::max(dom.coord(S[0], 1), dom.coord(S[1], 1));
                ymax=std::max(ymax, dom.coord(S[2], 1));
                zmin=std::min(dom.coord(S[0], 2), dom.coord(S[1], 2));
                zmin=std::min(zmin, dom.coord(S[2], 2));
                zmax=std::max(dom.coord(S[0], 2), dom.coord(S[1], 2));
                zmax=std::max(zmax, dom.coord(S[2], 2));

                if(est_egal(zmin, zmax))
                  {
                    for(i=0; i<4; i++)
                      if ( est_egal(dom.coord(S[i], 0),xmin) && est_egal(dom.coord(S[i], 1),ymin))
                        NS[0]=S[i];
                    for(i=0; i<4; i++)
                      if ( !est_egal(dom.coord(S[i], 0),xmin) && est_egal(dom.coord(S[i], 1),ymin))
                        NS[1]=S[i];
                    for(i=0; i<4; i++)
                      if ( est_egal(dom.coord(S[i], 0),xmin) && !est_egal(dom.coord(S[i], 1),ymin))
                        NS[2]=S[i];
                    for(i=0; i<4; i++)
                      if ( !est_egal(dom.coord(S[i], 0),xmin) && !est_egal(dom.coord(S[i], 1),ymin))
                        NS[3]=S[i];
                  }
                if(est_egal(ymin, ymax))
                  {
                    for(i=0; i<4; i++)
                      if ( est_egal(dom.coord(S[i], 0),xmin) && est_egal(dom.coord(S[i], 2),zmin))
                        NS[0]=S[i];
                    for(i=0; i<4; i++)
                      if ( !est_egal(dom.coord(S[i], 0),xmin) && est_egal(dom.coord(S[i], 2),zmin))
                        NS[1]=S[i];
                    for(i=0; i<4; i++)
                      if ( est_egal(dom.coord(S[i], 0),xmin) && !est_egal(dom.coord(S[i], 2),zmin))
                        NS[2]=S[i];
                    for(i=0; i<4; i++)
                      if ( !est_egal(dom.coord(S[i], 0),xmin) && !est_egal(dom.coord(S[i], 2),zmin))
                        NS[3]=S[i];
                  }
                if(est_egal(xmin, xmax))
                  {
                    for(i=0; i<4; i++)
                      if ( est_egal(dom.coord(S[i], 1),ymin) && est_egal(dom.coord(S[i], 2),zmin))
                        NS[0]=S[i];
                    for(i=0; i<4; i++)
                      if ( !est_egal(dom.coord(S[i], 1),ymin) && est_egal(dom.coord(S[i], 2),zmin))
                        NS[1]=S[i];
                    for(i=0; i<4; i++)
                      if ( est_egal(dom.coord(S[i], 1),ymin) && !est_egal(dom.coord(S[i], 2),zmin))
                        NS[2]=S[i];
                    for(i=0; i<4; i++)
                      if ( !est_egal(dom.coord(S[i], 1),ymin) && !est_egal(dom.coord(S[i], 2),zmin))
                        NS[3]=S[i];
                  }

                for(i=0; i<4; i++)
                  {
                    assert(NS[i]!=-1);
                    sommet(face, i)=NS[i];
                  }
              }
          }
        break;
      }
    case Type_Face::quadrangle_3D_axi :
      {
        assert(dimension==3);
        break;
      }
    case Type_Face::polygone_3D:
      assert(dimension==3);
      break;
    default :
      {
        Cerr << "Error TRUST in type of Faces_32_64 not recognized " << finl;
        exit();
      }
    }
}

/*! @brief Compares the Faces_32_64 object passed as parameter with this Faces_32_64 object.
 *
 * @param (Faces_32_64& faces) the faces to compare against
 * @param (IntVect& renum) the renumbering vector - renum has size 1 and contains -1 if the number of faces differs; otherwise has as many elements as faces and renum(i) = the index of face i from the parameter faces in the current object
 * @return (IntVect&)
 */
template <typename _SIZE_>
typename Faces_32_64<_SIZE_>::IntVect_t& Faces_32_64<_SIZE_>::compare(const Faces_32_64& faces, IntVect_t& renum)
{
  const Domaine_t& domaine = this->domaine();
  const Domaine_t& son_domaine=faces.domaine();
  if ( (nb_faces() != faces.nb_faces()) || ( type_face_ != faces.type_face_) )
    {
      renum.resize(1);
      renum=-1;
      return renum;
    }
  if(domaine.le_nom() == son_domaine.le_nom())
    {
      // just compare the vertex indices:
      TRUSTLists<_SIZE_> listes;
      int_t premier=0;
      int_t numerol, face;
      int numeroli;
      int ok;
      int_t rang;
      for(face=0; face<nb_faces(); face++)
        {
          numerol=this->ppsf(face, nb_som_face);
          assert(numerol < std::numeric_limits<int>::max());
          numeroli = (int)numerol;
          listes[numeroli].add(premier++);
        }
      for(face=0; face<nb_faces(); face++)
        {
          numerol=faces.ppsf(face, nb_som_face);
          assert(numerol < std::numeric_limits<int>::max());
          numeroli = (int)numerol;
          TRUSTList_Curseur<_SIZE_> curseur(listes[numeroli]);
          ok=1;
          while (curseur && ok)
            {
              rang=curseur.valeur();
              ok=!faces.same_face(face, *this, rang, nb_som_face);
              if(!ok)
                {
                  // "face==rang"
                  renum(face)=rang;
                }
              else
                ++curseur;
            }
          if(!curseur)
            {
              renum.resize(1);
              renum=-1;
              return renum;
            }
        }
    }
  else
    {
      // we need to compare the vertex coordinates:
    }
  return renum;
}


template class Faces_32_64<int>;
#if INT_is_64_ == 2
template class Faces_32_64<trustIdType>;
#endif
