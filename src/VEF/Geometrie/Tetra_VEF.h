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

#ifndef Tetra_VEF_included
#define Tetra_VEF_included

#include <Elem_VEF_base.h>

class Tetra_VEF : public Elem_VEF_base
{

  Declare_instanciable_sans_constructeur(Tetra_VEF);

public:
  Tetra_VEF();
  inline int nb_facette() const override
  {
    return 6;
  };
  void creer_facette_normales(const Domaine_VEF&, const IntVect& ) const override;
  void creer_normales_facettes_Cl(DoubleTab&, int ,int ,
                                  const DoubleTab& ,const DoubleVect& , const Domaine&) const override ;
  void modif_volumes_entrelaces(int ,int ,const Domaine_VEF& ,DoubleVect& ,int ) const override ;
  void modif_volumes_entrelaces_faces_joints(int ,int ,const Domaine_VEF& ,DoubleVect& ,int ) const override ;
  void modif_normales_facettes_Cl(DoubleTab& ,int ,int ,int ,int, int ,int ) const override ;
  void calcul_vc(const ArrOfInt& ,ArrOfDouble& ,const ArrOfDouble& ,
                 const DoubleTab& ,const Champ_Inc_base& ,int, const DoubleVect& ) const override ;
  void calcul_xg(DoubleVect& ,const DoubleTab& ,const int ,int& ,
                 int& ,int& ,int& ) const override ;
  void creer_face_normales(DoubleTab&, const IntTab& ,const IntTab&,
                           const IntTab& ,const Domaine& )  const override ;
};

KOKKOS_INLINE_FUNCTION void calcul_vc_tetra(const int* Face, double *vc, const double * vs, const double * vsom,
                                            const double* vitesse,int type_cl, const double* poro)
{
  // Casting (justified by size) type_cl and comp to True_int to avoid nvc++ bug
  int comp;
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
          vc[comp] = vitesse[9+comp] * poro[3];
        break;
      }

    case 2: // the tetrahedron has one Dirichlet face: KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse[6+comp]  *  poro[2];
        break;
      }

    case 4: // the tetrahedron has one Dirichlet face: KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse[3+comp]  *  poro[1];
        break;
      }

    case 8: // the tetrahedron has one Dirichlet face: KEL0
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse[comp] *  poro[0];
        break;
      }

    case 3: // the tetrahedron has two Dirichlet faces: KEL3 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[comp] + vsom[3+comp]);
        break;
      }

    case 5: // the tetrahedron has two Dirichlet faces: KEL3 and KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[comp] + vsom[6+comp]);
        break;
      }

    case 6: // the tetrahedron has two Dirichlet faces: KEL1 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[comp] + vsom[9+comp]);
        break;
      }

    case 9: // the tetrahedron has two Dirichlet faces: KEL0 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[3+comp] + vsom[6+comp]);
        break;
      }

    case 10: // the tetrahedron has two Dirichlet faces: KEL0 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[3+comp] + vsom[9+comp]);
        break;
      }

    case 12: // the tetrahedron has two Dirichlet faces: KEL0 and KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5*(vsom[6+comp] + vsom[9+comp]);
        break;
      }

    case 7: // the tetrahedron has three Dirichlet faces: KEL1, KEL2 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[comp];
        break;
      }

    case 11: // the tetrahedron has three Dirichlet faces: KEL0, KEL2 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[3+comp];
        break;
      }

    case 13: // the tetrahedron has three Dirichlet faces: KEL0, KEL1 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[6+comp];
        break;
      }

    case 14: // the tetrahedron has three Dirichlet faces: KEL0, KEL1 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[9+comp];
        break;
      }

    } // end of switch

}

KOKKOS_INLINE_FUNCTION void calcul_vc_tetra_views(const int* Face, double *vc, const double * vs, const double * vsom,
                                                  CDoubleTabView vitesse,int type_cl, CDoubleArrView porosite_face)
{
  // Casting (justified by size) type_cl and comp to True_int to avoid nvc++ bug
  int comp;
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
          vc[comp] = vitesse(Face[3], comp) * porosite_face(Face[3]);
        break;
      }

    case 2: // the tetrahedron has one Dirichlet face: KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse(Face[2], comp) * porosite_face(Face[2]);
        break;
      }

    case 4: // the tetrahedron has one Dirichlet face: KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse(Face[1], comp) * porosite_face(Face[1]);
        break;
      }

    case 8: // the tetrahedron has one Dirichlet face: KEL0
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vitesse(Face[0], comp) * porosite_face(Face[0]);
        break;
      }

    case 3: // the tetrahedron has two Dirichlet faces: KEL3 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[comp] + vsom[3+comp]);
        break;
      }

    case 5: // the tetrahedron has two Dirichlet faces: KEL3 and KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[comp] + vsom[6+comp]);
        break;
      }

    case 6: // the tetrahedron has two Dirichlet faces: KEL1 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[comp] + vsom[9+comp]);
        break;
      }

    case 9: // the tetrahedron has two Dirichlet faces: KEL0 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[3+comp] + vsom[6+comp]);
        break;
      }

    case 10: // the tetrahedron has two Dirichlet faces: KEL0 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5* (vsom[3+comp] + vsom[9+comp]);
        break;
      }

    case 12: // the tetrahedron has two Dirichlet faces: KEL0 and KEL1
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = 0.5*(vsom[6+comp] + vsom[9+comp]);
        break;
      }

    case 7: // the tetrahedron has three Dirichlet faces: KEL1, KEL2 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[comp];
        break;
      }

    case 11: // the tetrahedron has three Dirichlet faces: KEL0, KEL2 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[3+comp];
        break;
      }

    case 13: // the tetrahedron has three Dirichlet faces: KEL0, KEL1 and KEL3
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[6+comp];
        break;
      }

    case 14: // the tetrahedron has three Dirichlet faces: KEL0, KEL1 and KEL2
      {
        for (comp=0; comp<3; comp++)
          vc[comp] = vsom[9+comp];
        break;
      }

    } // end of switch

}

/*! @brief Computes the coordinates xg of the center of a non-standard element.
 * Also computes idirichlet = number of Dirichlet faces of the element.
 *
 */
KOKKOS_INLINE_FUNCTION
void calcul_xg_tetra(double * xg, const double *x, const int type_elem_Cl, int& idirichlet,int& n1,int& n2,int& n3)
{
  // Casting (justified by size) type_elem_cl and comp to True_int to avoid nvc++ bug
  int dim = 3;
  switch(type_elem_Cl)
    {
    case 0:  // the tetrahedron has no Dirichlet face. It has 6 facets
      {
        for (int j=0; j<dim; j++)
          xg[j]=0.25*(x[j]+x[dim+j]+x[2*dim+j]+x[3*dim+j]);

        idirichlet=0;
        break;
      }

    case 1:  // the tetrahedron has one Dirichlet face. The 'center'
      // of the tetrahedron is at the midpoint of face 3 with vertices 0, 1, 2.
      // It has 3 real facets: 0   at nodes 2 3 xg
      //                       1   at nodes 1 3 xg
      //                       3   at nodes 3 0 xg
      // the 3 other facets are on face 3

      {
        for (int j=0; j<dim; j++)
          xg[j]=(x[j]+x[dim+j]+x[2*dim+j])/3.;

        idirichlet=1;
        break;

      }

    case 2:  // the tetrahedron has one Dirichlet face. The 'center'
      // of the tetrahedron is at the midpoint of face 2 with vertices 0, 1, 3.
      // It has 3 real facets: 0   at nodes 2 3 xg
      //                       2   at nodes 1 2 xg
      //                       4   at nodes 2 0 xg

      {
        for (int j=0; j<dim; j++)
          xg[j]=(x[j]+x[dim+j]+x[3*dim+j])/3.;

        idirichlet=1;
        break;
      }

    case 4:  // the tetrahedron has one Dirichlet face. The 'center'
      // of the tetrahedron is at the midpoint of face 1 with vertices 0, 2, 3.
      // It has 3 real facets: 1   at nodes 1 3 xg
      //                       2   at nodes 1 2 xg
      //                       5   at nodes 1 0 xg

      {
        for (int j=0; j<dim; j++)
          xg[j]=(x[j]+x[2*dim+j]+x[3*dim+j])/3.;

        idirichlet=1;
        break;
      }

    case 8:  // the tetrahedron has one Dirichlet face. The 'center'
      // of the tetrahedron is at the midpoint of face 0 with vertices 1, 2, 3.
      // It has 3 real facets: 3   at nodes 3 0 xg
      //                       4   at nodes 2 0 xg
      //                       5   at nodes 1 0 xg

      {
        for (int j=0; j<dim; j++)
          xg[j]=(x[dim+j]+x[2*dim+j]+x[3*dim+j])/3.;

        idirichlet=1;
        break;
      }

    case 3:  // the tetrahedron has two Dirichlet faces 2 and 3. The 'center'
      // is at the midpoint of the edge with endpoints 0, 1.
      // It has 1 null facet: 5

      {
        for (int j=0; j<dim; j++)
          xg[j]= 0.5*(x[j]+x[dim+j]);

        n1=5;
        idirichlet=2;
        break;
      }


    case 5:  // the tetrahedron has two Dirichlet faces 3 and 1. The 'center'
      // is at the midpoint of the edge with endpoints 0, 2.
      // It has 1 null facet: 4

      {
        for (int j=0; j<dim; j++)
          xg[j]= 0.5*(x[j]+x[2*dim+j]);

        n1=4;
        idirichlet=2;
        break;
      }

    case 6:  // the tetrahedron has two Dirichlet faces 1 and 2. The 'center'
      // is at the midpoint of the edge with endpoints 0, 3.
      // It has 1 null facet: 3

      {
        for (int j=0; j<dim; j++)
          xg[j]= 0.5*(x[j]+x[3*dim+j]);

        n1=3;
        idirichlet=2;
        break;
      }

    case 9:  // the tetrahedron has two Dirichlet faces 0 and 3. The 'center'
      // is at the midpoint of the edge with endpoints 1, 2.
      // It has 1 null facet: 2

      {
        for (int j=0; j<dim; j++)
          xg[j]= 0.5*(x[dim+j]+x[2*dim+j]);

        n1=2;
        idirichlet=2;
        break;
      }

    case 10:  // the tetrahedron has two Dirichlet faces 0 and 2. The 'center'
      // is at the midpoint of the edge with endpoints 1, 3.
      // It has 1 null facet: 1

      {
        for (int j=0; j<dim; j++)
          xg[j]= 0.5*(x[dim+j]+x[3*dim+j]);

        n1=1;
        idirichlet=2;
        break;
      }


    case 12:  // the tetrahedron has two Dirichlet faces 0 and 1. The 'center'
      // is at the midpoint of the edge with vertices 2, 3.
      // It has 1 null facet

      {
        for (int j=0; j<dim; j++)
          xg[j]= 0.5*(x[2*dim+j]+x[3*dim+j]);

        n1=0;
        idirichlet=2;
        break;
      }

    case 7:  // three Dirichlet faces: 1, 2, 3. The center is at vertex 0.
      // There are 3 null facets: 3, 4, 5

      {
        for (int j=0; j<dim; j++)
          xg[j]= x[j];

        n1=3;
        n2=4;
        n3=5;
        idirichlet=3;
        break;

      }

    case 11:  // three Dirichlet faces: 0, 2, 3. The center is at vertex 1.
      // There are 3 null facets: 1, 2, 5

      {
        for (int j=0; j<dim; j++)
          xg[j]= x[dim+j];

        n1=1;
        n2=2;
        n3=5;
        idirichlet=3;
        break;

      }

    case 13:  // three Dirichlet faces: 0, 1, 3. The center is at vertex 2.
      // There are 3 null facets: 0, 2, 4

      {
        for (int j=0; j<dim; j++)
          xg[j]= x[2*dim+j];

        n1=0;
        n2=2;
        n3=4;
        idirichlet=3;
        break;

      }
    case 14:  // three Dirichlet faces: 0, 1, 2. The center is at vertex 3.
      // There are 3 null facets: 0, 1, 3

      {
        for (int j=0; j<dim; j++)
          xg[j]= x[3*dim+j];

        n1=0;
        n2=1;
        n3=3;
        idirichlet=3;
        break;

      }
    }
}
#endif



