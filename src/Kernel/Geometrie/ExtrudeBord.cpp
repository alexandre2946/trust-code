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

#include <ExtrudeBord.h>
#include <EFichier.h>
#include <TroisDto2D.h>
#include <Transformer.h>
#include <SFichier.h>
#include <Extruder.h>
#include <Extruder_en3.h>
#include <Extruder_en20.h>
#include <Scatter.h>
#include <Param.h>
#include <TRUSTVects.h>

Implemente_instanciable(ExtrudeBord,"ExtrudeBord",Interprete_geometrique_base);
// XD extrudebord interprete extrudebord BRACE Class to generate an extruded mesh from a boundary of a tetrahedral or an
// XD_CONT hexahedral mesh. NL2 Warning: If the initial domain is a tetrahedral mesh, the boundary will be moved in the
// XD_CONT XY plane then extrusion will be applied (you should maybe use the Transformer keyword on the final domain to
// XD_CONT have the domain you really want). You can use the keyword Postraiter_domaine to generate a lata|med|... file
// XD_CONT to visualize your initial and final meshes. NL2 This keyword can be used for example to create a periodic box
// XD_CONT extracted from a boundary of a tetrahedral or a hexaedral mesh. This periodic box may be used then to
// XD_CONT engender turbulent inlet flow condition for the main domain.NL2 Note that ExtrudeBord in VEF generates 3 or
// XD_CONT 14 tetrahedra from extruded prisms.

Sortie& ExtrudeBord::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}

Entree& ExtrudeBord::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

Entree& ExtrudeBord::interpreter_(Entree& is)
{
  Nom nom_dom_volumique;
  Nom nom_front,nom_dom_surfacique;
  DoubleVect vect_dir(3);
  int nbpas = -1;

  if(dimension!=3)
    {
      Cerr << "Interpreter "<<que_suis_je()<<" can be used only in 3 dimensions." <<finl;
      exit();
    }
  Param param(que_suis_je());
  param.ajouter("domaine_init",&nom_dom_volumique,Param::REQUIRED);// XD_ADD_P ref_domaine
  // XD_CONT Initial domain with hexaedras or tetrahedras.
  param.ajouter_arr_size_predefinie("direction",&vect_dir,Param::REQUIRED);// XD_ADD_P listf
  // XD_CONT Directions for the extrusion.
  param.ajouter("nb_tranches",&nbpas,Param::REQUIRED);// XD_ADD_P entier
  // XD_CONT Number of elements in the extrusion direction.
  param.ajouter("domaine_final",&nom_dom_surfacique,Param::REQUIRED);// XD_ADD_P chaine
  // XD_CONT Extruded domain.
  param.ajouter("nom_bord",&nom_front,Param::REQUIRED);// XD_ADD_P chaine
  // XD_CONT Name of the boundary of the initial domain where extrusion will be applied.
  param.ajouter_flag("hexa_old",&hexa_old);// XD_ADD_P rien
  // XD_CONT Old algorithm for boundary extrusion from a hexahedral mesh.
  param.ajouter_flag("Trois_Tetra",&Trois_Tetra);// XD_ADD_P rien
  // XD_CONT To extrude in 3 tetrahedras instead of 14 tetrahedras.
  param.ajouter_flag("Vingt_Tetra",&Vingt_Tetra);// XD_ADD_P rien
  // XD_CONT To extrude in 20 tetrahedras instead of 14 tetrahedras.
  param.ajouter("sans_passer_par_le2D",&en3D_);// XD_ADD_P entier
  // XD_CONT Only for non-regression
  param.lire_avec_accolades_depuis(is);

  associer_domaine(nom_dom_volumique);

  const Domaine& dom=domaine();

  if (dom.nb_som_elem()==8 && (hexa_old))
    extruder_hexa_old(nom_front, nom_dom_surfacique, vect_dir, nbpas);
  else if (dom.nb_som_elem()==4 || dom.nb_som_elem()==8)
    extruder_bord(nom_front, nom_dom_surfacique, vect_dir, nbpas);

  Domaine& dom_surfacique=ref_cast(Domaine, objet(nom_dom_surfacique));
  Scatter::init_sequential_domain(dom_surfacique);
  return is;
}

void ExtrudeBord::extruder_bord(Nom& nom_front, Nom& nom_dom_surfacique, DoubleVect& vect_dir, int nbpas)
{
  const Domaine& dom=domaine();

  const Bord& front=dom.bord(nom_front);

  Domaine& dom_surfacique=ref_cast(Domaine, objet(nom_dom_surfacique));

  // Extract from the volume domain dom the boundary front and fill the dom_surfacique domain:
  TroisDto2D tD2dD;
  int& coupe=tD2dD.coupe();
  coupe=0;
  if (en3D_)
    coupe=2;
  tD2dD.extraire_2D(dom,dom_surfacique,front,nom_front,0);

  double xa = tD2dD.getXa();
  double ya = tD2dD.getYa();
  double za = tD2dD.getZa();
  double Ix = tD2dD.getIx();
  double Iy = tD2dD.getIy();
  double Iz = tD2dD.getIz();
  double Jx = tD2dD.getJx();
  double Jy = tD2dD.getJy();
  double Jz = tD2dD.getJz();
  double Kx = tD2dD.getKx();
  double Ky = tD2dD.getKy();
  double Kz = tD2dD.getKz();

  if (!en3D_)
    {
      dimension=2;
    }
  else
    {
      Ix=1;
      Iy=0;
      Iz=0;
      Jx=0;
      Jy=1;
      Jz=0;
      Kx=0;
      Ky=0;
      Kz=1;

    }

  if (Trois_Tetra)  // Evolution of the ExtrudeBord method adapted for 3tetra
    {
      Extruder_en3 extr3;
      extr3.setNbTranches(nbpas);
      extr3.setDirection(vect_dir(0)*Ix+vect_dir(1)*Iy+vect_dir(2)*Iz,
                         vect_dir(0)*Jx+vect_dir(1)*Jy+vect_dir(2)*Jz,
                         vect_dir(0)*Kx+vect_dir(1)*Ky+vect_dir(2)*Kz);

      IntVect num;
      dom_surfacique.ajouter(dom_surfacique.coord_sommets(), num);
      extr3.extruder(dom_surfacique,num);

      if ((dom_surfacique.type_elem()->que_suis_je())== "Rectangle")
        {
          Cerr << "It is not possible to apply ExtrudeBord with 3Tetra option to : "   <<  dom_surfacique.type_elem()->que_suis_je() << " mesh cells" << finl;
          // exit();
        }
    }
  else if (Vingt_Tetra)
    {
      Extruder_en20 extr;
      extr.setNbTranches(nbpas);
      extr.setDirection(vect_dir(0)*Ix+vect_dir(1)*Iy+vect_dir(2)*Iz,
                        vect_dir(0)*Jx+vect_dir(1)*Jy+vect_dir(2)*Jz,
                        vect_dir(0)*Kx+vect_dir(1)*Ky+vect_dir(2)*Kz);
      extr.extruder(dom_surfacique);
    }
  else
    {
      Extruder extr;
      extr.setNbTranches(nbpas);
      extr.setDirection(vect_dir(0)*Ix+vect_dir(1)*Iy+vect_dir(2)*Iz,
                        vect_dir(0)*Jx+vect_dir(1)*Jy+vect_dir(2)*Jz,
                        vect_dir(0)*Kx+vect_dir(1)*Ky+vect_dir(2)*Kz);
      extr.extruder(dom_surfacique);
    }
  if (en3D_)
    {
      // in 3D we do not rotate the boundary
      return;
    }
  // put the new 3D domain back into the coordinate system of the initial 3D domain
  dimension = 3;

  Noms les_fcts(3);
  les_fcts[0] = Nom(xa,"%.16e")+Nom("+x*(")+Nom(Ix,"%.16e")+Nom(")+y*(")+Nom(Jx,"%.16e")+Nom(")+z*(")+Nom(Kx,"%.16e")+Nom(")");
  les_fcts[1] = Nom(ya,"%.16e")+Nom("+x*(")+Nom(Iy,"%.16e")+Nom(")+y*(")+Nom(Jy,"%.16e")+Nom(")+z*(")+Nom(Ky,"%.16e")+Nom(")");
  les_fcts[2] = Nom(za,"%.16e")+Nom("+x*(")+Nom(Iz,"%.16e")+Nom(")+y*(")+Nom(Jz,"%.16e")+Nom(")+z*(")+Nom(Kz,"%.16e")+Nom(")");

  Transformer transf;
  transf.transformer(dom_surfacique,les_fcts);
  // To avoid precision problems later when manipulating
  // the extruded domain, we impose the coordinates of the nom_front boundary of the
  // initial domain onto the coordinates of the "front" boundary of the extruded domain.
  // We also need to properly recompute the "back" boundary.
  const Faces& faces=dom.frontiere(dom.rang_frontiere(nom_front)).faces();
  int nb_faces=faces.nb_faces();
  int nb_som_faces=faces.nb_som_faces();

  const Faces& faces2=dom_surfacique.frontiere(dom_surfacique.rang_frontiere("devant")).faces();
  const Faces& faces3=dom_surfacique.frontiere(dom_surfacique.rang_frontiere("derriere")).faces();
  for (int j=0; j<nb_faces; j++)
    for (int k=0; k<nb_som_faces; k++)
      {
        int som = faces.sommet(j,k);
        int som2 = faces2.sommet(j,k);
        int som3 = faces3.sommet(j,k);
        for (int i=0; i<3; i++)
          {
            dom_surfacique.coord(som2,i) = dom.coord(som,i);
            dom_surfacique.coord(som3,i) = dom.coord(som,i) + vect_dir(i);
          }
      }
  Cerr << "ExtrudeBord does not create any more PERIO boundary." << finl;
  Cerr << "Use RegroupeBord if you want to create periodic boundary condition..." << finl;
}

void ExtrudeBord::extruder_hexa_old(Nom& nom_front, Nom& nom_dom_surfacique, DoubleVect& vect_dir, int nbpas)
{

  const Domaine& dom=domaine();

  int nbfr=dom.nb_front_Cl();


  Domaine& dom_surfacique=ref_cast(Domaine, objet(nom_dom_surfacique));

  dom_surfacique.typer("Hexaedre");
  DoubleTab& sommets2=dom_surfacique.les_sommets();
  IntTab& les_elems2=dom_surfacique.les_elems();
  int i;

  for (int l=0; l<nbfr; l++)
    {
      const Frontiere& fr=dom.frontiere(l);
      const Nom& nomfr=fr.le_nom();

      if(nomfr==nom_front)
        {
          int nbfaces=fr.nb_faces();
          const IntTab& som=fr.les_sommets_des_faces();
          IntTab newsom(4*nbfaces);
          IntTab compt(4*nbfaces);
          newsom=-100;
          compt=0;

          int nbsombord,trouve;
          nbsombord = 0;

          // ******** Correspondence table of boundary vertices: dom_volumique - dom_surfacique

          for (i=0; i<nbfaces; i++)
            for (int j=0; j<4; j++)
              {
                trouve = 0;
                int sommet=som(i,j);

                for (int k=0; k<nbsombord+1; k++)
                  {
                    if(sommet==newsom(k))
                      {
                        k=nbsombord+1;
                        trouve = 1;
                        break;
                      }
                  }

                if (trouve==0)
                  {
                    newsom[nbsombord]=sommet;
                    nbsombord++;
                  }
              }

          newsom.resize(nbsombord);
          compt.resize(nbsombord);
          sommets2.resize((nbpas+1)*nbsombord,3);

          // ******** Creation of vertices derived by translation

          for (int j=0; j<nbpas+1; j++)
            for (i=0; i<nbsombord; i++)
              {
                sommets2(i+j*nbsombord,0)=dom.coord(newsom(i),0)+j*vect_dir(0)/nbpas;
                sommets2(i+j*nbsombord,1)=dom.coord(newsom(i),1)+j*vect_dir(1)/nbpas;
                sommets2(i+j*nbsombord,2)=dom.coord(newsom(i),2)+j*vect_dir(2)/nbpas;
              }


          // ******** Definition of element-to-vertex correspondence arrays

          les_elems2.resize(nbpas*nbfaces,8);

          int  sommet0,sommet1,sommet2,sommet3;
          int  newsom0=-1,newsom1=-1,newsom2=-1,newsom3=-1;

          for (i=0; i<nbfaces; i++)
            {
              sommet0=som(i,0);
              sommet1=som(i,1);
              sommet2=som(i,2);
              sommet3=som(i,3);


              for(int j=0; j<nbsombord; j++)
                {
                  if(newsom[j]==sommet0) newsom0=j;
                  if(newsom[j]==sommet1) newsom1=j;
                  if(newsom[j]==sommet2) newsom2=j;
                  if(newsom[j]==sommet3) newsom3=j;
                }

              for(int k=0; k<nbpas; k++)
                {
                  les_elems2(i+k*nbfaces,0)=newsom0+k*nbsombord;
                  les_elems2(i+k*nbfaces,1)=newsom1+k*nbsombord;
                  les_elems2(i+k*nbfaces,2)=newsom2+k*nbsombord;
                  les_elems2(i+k*nbfaces,3)=newsom3+k*nbsombord;
                  les_elems2(i+k*nbfaces,4)=newsom0+(k+1)*nbsombord;
                  les_elems2(i+k*nbfaces,5)=newsom1+(k+1)*nbsombord;
                  les_elems2(i+k*nbfaces,6)=newsom2+(k+1)*nbsombord;
                  les_elems2(i+k*nbfaces,7)=newsom3+(k+1)*nbsombord;
                }
            }



          // ******** Definition of periodic boundaries


          Bord& bordperio=dom_surfacique.faces_bord().add(Bord());
          bordperio.nommer("PERIO");
          bordperio.typer_faces("QUADRANGLE_3D");
          bordperio.dimensionner(0);
          IntTab faces_perio(2*nbfaces,4);

          for (i=0; i<nbfaces; i++)
            {
              sommet0=som(i,0);
              sommet1=som(i,1);
              sommet2=som(i,2);
              sommet3=som(i,3);


              for(int j=0; j<nbsombord; j++)
                {
                  if(newsom[j]==sommet0) newsom0=j;
                  if(newsom[j]==sommet1) newsom1=j;
                  if(newsom[j]==sommet2) newsom2=j;
                  if(newsom[j]==sommet3) newsom3=j;
                }

              faces_perio(i,0)=newsom0;
              faces_perio(i,1)=newsom1;
              faces_perio(i,2)=newsom2;
              faces_perio(i,3)=newsom3;

              faces_perio(i+nbfaces,0)=newsom0+nbpas*nbsombord;
              faces_perio(i+nbfaces,1)=newsom1+nbpas*nbsombord;
              faces_perio(i+nbfaces,2)=newsom2+nbpas*nbsombord;
              faces_perio(i+nbfaces,3)=newsom3+nbpas*nbsombord;
            }


          bordperio.ajouter_faces(faces_perio);




          // ******** Definition of wall boundaries


          Bord& bordparoi=dom_surfacique.faces_bord().add(Bord());
          bordparoi.nommer("PAROI");
          bordparoi.typer_faces("QUADRANGLE_3D");
          bordparoi.dimensionner(0);


          int nbsombordparoi=0;

          // ***** search for vertices belonging to the contour of the boundary face *****

          for (i=0; i<nbfaces; i++)
            for (int j=0; j<4; j++)
              for (int k=0; k<nbsombord; k++)
                if(som(i,j)==newsom(k))  compt(k)++; // number of faces referencing each vertex

          for (int k=0; k<nbsombord; k++)
            if(compt(k)<3) nbsombordparoi++; // number of vertices to extrude to form the wall faces

          // PQ : 05/05/04 : WARNING: internal vertices can be connected to only
          //                        3 faces (e.g. "O"-type mesh of a cylinder)
          //                        -> we therefore consider only the case where boundary vertices are connected to
          //                        2 elements (excluding VDF-type stair-step boundaries).


          // **** Extension of edges to form the faces of the wall boundary

          IntTab faces_paroi(nbsombordparoi*nbpas,4);
          int ip=0;

          for (i=0; i<nbfaces; i++)
            {
              sommet0=som(i,0);
              sommet1=som(i,1);
              sommet2=som(i,2);
              sommet3=som(i,3);

              int compt0=4,compt1=4,compt2=4,compt3=4;

              for(int j=0; j<nbsombord; j++)
                {
                  if(newsom[j]==sommet0)
                    {
                      newsom0=j;
                      compt0=compt(j);
                    }
                  if(newsom[j]==sommet1)
                    {
                      newsom1=j;
                      compt1=compt(j);
                    }
                  if(newsom[j]==sommet2)
                    {
                      newsom2=j;
                      compt2=compt(j);
                    }
                  if(newsom[j]==sommet3)
                    {
                      newsom3=j;
                      compt3=compt(j);
                    }
                }


              if( (compt0<3) && (compt1<3) ) // edge[0-1]
                {
                  for(int k=0; k<nbpas; k++)
                    {
                      faces_paroi(ip+k,0)=newsom0+k*nbsombord;
                      faces_paroi(ip+k,1)=newsom0+(k+1)*nbsombord;
                      faces_paroi(ip+k,2)=newsom1+k*nbsombord;
                      faces_paroi(ip+k,3)=newsom1+(k+1)*nbsombord;
                    }
                  ip+=nbpas;
                }


              if( (compt1<3) && (compt3<3) ) // edge[1-3]
                {
                  for(int k=0; k<nbpas; k++)
                    {
                      faces_paroi(ip+k,0)=newsom1+k*nbsombord;
                      faces_paroi(ip+k,1)=newsom1+(k+1)*nbsombord;
                      faces_paroi(ip+k,2)=newsom3+k*nbsombord;
                      faces_paroi(ip+k,3)=newsom3+(k+1)*nbsombord;
                    }
                  ip+=nbpas;
                }


              if( (compt3<3) && (compt2<3) ) // edge[3-2]
                {
                  for(int k=0; k<nbpas; k++)
                    {
                      faces_paroi(ip+k,0)=newsom3+k*nbsombord;
                      faces_paroi(ip+k,1)=newsom3+(k+1)*nbsombord;
                      faces_paroi(ip+k,2)=newsom2+k*nbsombord;
                      faces_paroi(ip+k,3)=newsom2+(k+1)*nbsombord;
                    }
                  ip+=nbpas;
                }


              if( (compt2<3) && (compt0<3) ) // edge[2-0]
                {
                  for(int k=0; k<nbpas; k++)
                    {
                      faces_paroi(ip+k,0)=newsom2+k*nbsombord;
                      faces_paroi(ip+k,1)=newsom2+(k+1)*nbsombord;
                      faces_paroi(ip+k,2)=newsom0+k*nbsombord;
                      faces_paroi(ip+k,3)=newsom0+(k+1)*nbsombord;
                    }
                  ip+=nbpas;
                }
            }

          bordparoi.ajouter_faces(faces_paroi);

        }//end of condition on the boundary to extrude
    }//end of loop over all boundaries
  Nom mfile=nom_dom_surfacique;
  mfile+=".geom";

  {
    SFichier file(mfile);
    file<<dom_surfacique;
    file.close();
  }
  {
    EFichier file(mfile);
    file>>dom_surfacique;
    file.close();
  }
}
