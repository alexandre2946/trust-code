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

#include <Prepro_IBM_Ponderation.h>
#include <Connectivite_som_elem.h>

Implemente_instanciable( Prepro_IBM_Ponderation,"Prepro_IBM_Ponderation|methode_IBM_ponderation",Prepro_IBM_base );

Sortie& Prepro_IBM_Ponderation::printOn(Sortie& os) const { return Prepro_IBM_base::printOn(os); }
// XD Prepro_IBM_Ponderation Prepro_IBM_base methode_IBM_ponderation BRACE To perform the intersection of an IB
// XD_CONT (Lagrange mesh) in a MED-format file .med with the Euler computional mesh.

void Prepro_IBM_Ponderation::set_param(Param& param) const
{
  Prepro_IBM_base::set_param(param);
  param.ajouter("type_de_ponderation",&pond_,Param::OPTIONAL); // XD_ADD_P entier
  // XD_CONT choice of the weighting method
}

Entree& Prepro_IBM_Ponderation::readOn(Entree& is)
{
  Prepro_IBM_base::readOn(is);
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);

  if(pond_==1)
    Cout<<"Weighting method = arimethic weight"<<finl;
  else if(pond_==2)
    Cout<<"Weighting method = area weight"<<finl;
  else if(pond_==3)
    Cout<<"Weighting method =  inverse distance weight"<<finl;
  else if(pond_==4)
    Cout<<"Weighting method = area and inverse distance weight"<<finl;
  else
    {
      Cerr<<"Prepro_IBM_ponderation : Type_de_ponderation : invalid argument = "<<pond_<<finl;
      Process::exit();
    }
  return is;
}

void Prepro_IBM_Ponderation::associer_pb(const Probleme_base& pb)
{
  Prepro_IBM_base::associer_pb(pb);
  compute_solid_fluid(0);

  // Optional output
  if( save_prepro_ == 1) Save_Med_File();
}

void Prepro_IBM_Ponderation::compute_solid_fluid(int maj_from_ext)
{
  if (verbose_) Cerr<<"Prepro_IBM_Ponderation:: Computing arrays..."<<finl;

  int dim_esp = Objet_U::dimension;
  DoubleTab& normalArray = champ_normal_->valeurs();
  DoubleTab& aireArray = champ_aire_->valeurs();
  DoubleTab& isNodeDirichletArray = isNodeDirichlet_->valeurs();
  DoubleTab& rotationArray = champ_rotation_->valeurs();
  int nbElemVol = rotationArray.dimension(0); // Number of cells in the volumetric mesh

  //optional update of the normal per element and Dirichlet nodes from champ_aire_ and champ_rotation_
  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();
  const IntTab& elems = le_dom_dis.domaine().les_elems();
  if (maj_from_ext == 1)
    {
      isNodeDirichletArray = -1.;
      normalArray = 0.;
      for (int elem = 0; elem <nbElemVol; elem++)
        {
          if (aireArray(elem) > 0.)
            {
              for (int k = 0; k <dim_esp; k++) normalArray(elem,k) = rotationArray(elem,3*k+(dim_esp-1));
              for (int l = 0; l < le_dom_dis.domaine().nb_som_elem(); l++)
                {
                  isNodeDirichletArray(elems(elem,l)) = 1.0 ;
                }
            }
        }
      normalArray.echange_espace_virtuel();
      isNodeDirichletArray.echange_espace_virtuel();
    }

  // compute neighboring nodes of a given node
  int nb_niveau = 1;
  compute_NeighNode(nb_niveau);

  // compute characteristic distances
  compute_h_max_elem();

  //Solid projection: compute solid projected points and weighted contributions
  projectSolidPoints();

  //Fluid projection: using solid projections and normals, compute the fluid points
  projectFluidPoints();

//     //Update the elementary rotation matrix with nodal normals
  // for (int elem = 0; elem <nbElemVol; elem++)
  //   {
  // 	computeLocalFrame(normalArray, t1Arr, t2Arr, elem);
  // 	computeMatRot(normalArray, t1Arr, t2Arr, elem);
  //   }
}

void Prepro_IBM_Ponderation::projectSolidPoints()
{
  const DoubleTab& normalArray = champ_normal_->valeurs();  // Face normals
  const DoubleTab& aireArray = champ_aire_->valeurs();
  const DoubleTab& baryArray = champ_bary_->valeurs();      // Barycenters of surface faces
  DoubleTab& solideArray = solid_points_->valeurs();  // Array of solid projections
  const DoubleTab& isNodeDirichletArray = isNodeDirichlet_->valeurs();
  DoubleTab& solid_elemsArray = solid_elems_->valeurs();
  solideArray=0.0;
  solid_elemsArray = -2.;

  Cout << "============================================================================"<<finl;
  Cout<<"prepro_IBM::projectSolidPoints - Weighting method = "<<pond_<<finl;

  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();// contains TRUST mesh
  const Domaine& le_dom =  le_dom_dis.domaine();
  const int nb_elem = le_dom.nb_elem();
  const int nb_som = le_dom.nb_som();
  const int nb_som_elem =le_dom.nb_som_elem();
  const IntTab& connectDom = le_dom.les_elems() ;// element-node connectivity
  const DoubleTab coordsDom3D = le_dom.coord_sommets(); // coordinates of TRUST nodes
  const int dim_esp=Objet_U::dimension;
  assert(dim_esp == 3);

  DoubleTab x(dim_esp);
  DoubleTrav sumContribProjSolidArray(nb_som);
  sumContribProjSolidArray = 0.;
  IntLists elemPerNeighbor(nb_som);
  for (int e=0; e<nb_elem; e++)
    {
      if (aireArray(e) > eps_effec_)
        {
          //////////////////////////////////////////
          //compute solid projection for e
          //////////////////////////////////////////

          //Traverse the nodes of the element
          for (int noeud=0; noeud<nb_som_elem; noeud++)
            {
              int sommet = connectDom(e,noeud);
              if(dimTab_(0))
                x(0) = coordsDom3D(sommet,0) - baryArray(e,0);
              else
                x(0) = 0.0;
              if(dimTab_(1))
                x(1) = coordsDom3D(sommet,1) - baryArray(e,1);
              else
                x(1) = 0.0;
              if(dimTab_(2))
                x(2) = coordsDom3D(sommet,2) - baryArray(e,2);
              else
                x(2) = 0.0;
              // for (int d=0; d<dim_esp; d++) Cerr<<x(d)<<" ";
              // Cerr<<finl;
              // for (int d=0; d<dim_esp; d++) Cerr<<coordsDom3D(sommet,d)<<" ";
              // Cerr<<finl;

              double distance=0.0 , distance_normale=0.0;
              for (int d=0; d<dim_esp; d++)
                {
                  distance_normale += x(d)*normalArray(e,d);
                  distance += x(d)*x(d);
                }
              distance = sqrt(distance);
              // Cerr<<"distance_normale : "<<distance_normale<<" distance : "<<distance<<finl;

              // Weighting type 1 to 4
              double pond;
              switch (pond_)
                {
                case 1:
                  pond = 1.0;
                  break;
                case 2:
                  pond = aireArray(e);
                  break;
                case 3:
                  pond = (distance>eps_effec_)? (1.0/distance) : 1.0;
                  break;
                case 4:
                  pond = (distance>eps_effec_)? (aireArray(e)/distance) : aireArray(e);
                  break;
                default:
                  pond = 1.0;
                  break;
                }

              //compute solid projection for the nodes of element e
              for (int d=0; d<dim_esp; d++) solideArray(sommet,d) += pond * (coordsDom3D(sommet,d)- distance_normale*normalArray(e,d));
              sumContribProjSolidArray(sommet) += pond;

              /////////////////////////////////////////////////////////
              //compute solid projection for the neighbors of e
              /////////////////////////////////////////////////////////
              //compute solid projection for the neighbors of the nodes of element e
              IntList& voisins = sommets_voisins_[sommet];
              int taille = 0;
              if (!(voisins.est_vide())) taille = voisins.size();
              // Cerr<<"taille : "<<taille<<finl;
              // Cerr<<"isNodeDirichletArray : "<<isNodeDirichletArray(voisins[0])<<" elemPerNeighbor[num_som_v]: "<<elemPerNeighbor[voisins[0]].contient(e)<<finl;
              for (int l=0; l<taille; l++)
                {
                  int num_som_v= voisins[l];
                  if ( isNodeDirichletArray(num_som_v)<0 && !(elemPerNeighbor[num_som_v].contient(e)) )
                    {
                      // Add element e for this neighbor
                      // Cerr<<"isNodeDirichletArray = "<<isNodeDirichletArray(num_som_v)<<" elemPerNeighbor[num_som_v]: "<<elemPerNeighbor[num_som_v].contient(e)<<finl;
                      elemPerNeighbor[num_som_v].add(e);
                      // Compute the vector neighbor -> barycenter
                      double dx = coordsDom3D(num_som_v, 0) - baryArray(e, 0);
                      double dy = coordsDom3D(num_som_v, 1) - baryArray(e, 1);
                      double dz = coordsDom3D(num_som_v, 2) - baryArray(e, 2);

                      double distance_normale_v = dx * normalArray(e, 0) + dy * normalArray(e, 1) + dz * normalArray(e, 2);
                      double d1 = sqrt(dx * dx + dy * dy + dz * dz);

                      // weighting
                      double pond_v;
                      if (d1 > eps_effec_)
                        {
                          switch (pond_)
                            {
                            case 3:
                              pond_v = fabs(1.0 / d1);
                              break;
                            case 4:
                              pond_v = fabs(aireArray(e) / d1);
                              break;
                            default:
                              pond_v = 1.0;
                              break;
                            }
                        }
                      else
                        {
                          Cerr << "DIST_PROB : d1 ~ 0 for neighbor" << num_som_v << " of the element " << e << finl;
                          pond_v = (pond_ == 4) ? aireArray(e) : 1.0;
                        }

                      // Weighted contribution
                      solideArray(num_som_v, 0) += pond_v * (coordsDom3D(num_som_v, 0) - distance_normale_v * normalArray(e, 0));
                      solideArray(num_som_v, 1) += pond_v * (coordsDom3D(num_som_v, 1) - distance_normale_v * normalArray(e, 1));
                      solideArray(num_som_v, 2) += pond_v * (coordsDom3D(num_som_v, 2) - distance_normale_v * normalArray(e, 2));

                      sumContribProjSolidArray(num_som_v) += pond_v;
                    }
                }
            }
        }
    }

  for(int sommet=0; sommet<nb_som; sommet++)
    {
      if (sumContribProjSolidArray(sommet))
        {
          for (int d=0; d<dim_esp; d++) solideArray(sommet,d) /= sumContribProjSolidArray(sommet);
          int interSoElem = le_dom.chercher_elements(solideArray(sommet,0),solideArray(sommet,1),solideArray(sommet,2));
          if (interSoElem == -1)
            solid_elemsArray(sommet) = -1.;
          else
            solid_elemsArray(sommet) = float(interSoElem);

          // Cerr<<"solid point projection vertex: "<<sommet<<" : ";
          // for (int d=0; d<dim_esp; d++) Cerr<<solideArray(sommet,d)<<" ";
          // Cerr<<";solid element = "<<solid_elemsArray(sommet)<<" ";
          // Cerr<<finl;
        }
    }
}

void Prepro_IBM_Ponderation::projectFluidPoints()
{
  const DoubleTab& aire = champ_aire_->valeurs();             // Aire
  DoubleTab& normalArray = champ_normal_->valeurs();    // Face normals
  DoubleTab& fluide = fluid_points_->valeurs();         // Array of fluid projections
  const DoubleTab& dirichlet = isNodeDirichlet_->valeurs();
  const DoubleTab& hmax_node = h_max_node_->valeurs();
  DoubleVect& elem_fluide = fluid_elems_->valeurs();;


  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();
  const Domaine& le_dom = le_dom_dis.domaine();
  const DoubleTab& coordsDom3D = le_dom.coord_sommets();
  const IntTab& connect = le_dom.les_elems();
  const int nb_som_elem = le_dom.nb_som_elem();
  const int nb_som = le_dom.nb_som();
  const int nb_sommets_tot = le_dom.nb_som_tot();
  const int nb_elem = le_dom.nb_elem();
  int dim_esp = Objet_U::dimension;
  assert(dim_esp == 3);

  DoubleTab normale_som(nb_som,dim_esp);    // Normales nodales
  double fact = 1.1;
  double x,  y, z, xp, yp, zp, norme;

  int interFlElem;

  for (int som = 0; som < nb_som; som++)
    {
      if (dirichlet(som) == 1.)
        {
          //////////////////////////////////////////
          //compute fluid projection for som
          //////////////////////////////////////////
          x = coordsDom3D(som, 0);
          y = coordsDom3D(som, 1);
          z = coordsDom3D(som, 2);
          xp = solid_points_->valeurs()(som, 0);
          yp = solid_points_->valeurs()(som, 1);
          zp = solid_points_->valeurs()(som, 2);
          norme = sqrt((x - xp)*(x - xp) + (y - yp)*(y - yp) + (z - zp)*(z - zp));

          if (norme < eps_effec_)
            {
              elem_fluide(som) = -3; //O(h) because the projection and node are located at the same position (i.e. the node is already on the immersed boundary)
              for (int d = 0; d < dim_esp; d++)
                {
                  fluide(som, d) = coordsDom3D(som, d);
                  normale_som(som, d) = 0.0;
                }
            }
          else
            {
              for (int d = 0; d < dim_esp; d++)
                {
                  normale_som(som, d) = (coordsDom3D(som, d) - solid_points_->valeurs()(som, d)) / norme;
                  fluide(som, d) = dimTab_(d) ? coordsDom3D(som, d) + fact * hmax_node(som) * normale_som(som, d) : coordsDom3D(som, d);
                }
              interFlElem = le_dom.chercher_elements(fluide(som,0),fluide(som,1),fluide(som,2));
              if (interFlElem == -1) // element not found
                elem_fluide(som) = -3.;
              else
                {
                  elem_fluide(som) = float(interFlElem);
                  if (aire(interFlElem) > 0.0)
                    elem_fluide(som) = -1.; // O(h) element traversed by IBC
                  else
                    {
                      bool flag = false;
                      for (int j = 0; j < nb_som_elem; j++)
                        {
                          int s = connect(interFlElem, j);
                          if (dirichlet(s) == 1.0) flag = true; // element with a Dirichlet vertex
                        }
                      if (flag) elem_fluide(som) = -1.0; // for Mean Gradient
                    }
                }
            }

          // If no fluid element was found, extend the zone to neighboring elements taking into account the c_prepro_ coefficient
          if (c_prepro_ > 0 && elem_fluide(som) == -1)
            {
              fact  += c_prepro_;
              for (int d = 0; d < dim_esp; d++)
                fluide(som, d) = dimTab_(d) ? coordsDom3D(som, d) + fact * hmax_node(som) * normale_som(som, d) : coordsDom3D(som, d);
              interFlElem = le_dom.chercher_elements(fluide(som,0),fluide(som,1),fluide(som,2));
              if (interFlElem == -1) // element not found
                elem_fluide(som) = -3.;
              else
                {
                  elem_fluide(som) = float(interFlElem);
                  if (aire(interFlElem) > 0.0)
                    elem_fluide(som) = -1.; // O(h) element traversed by IBC
                  else
                    {
                      bool flag = false;
                      for (int j = 0; j < nb_som_elem; j++)
                        {
                          int s = connect(interFlElem, j);
                          if (dirichlet(s) == 1.0) flag = true; // element with a Dirichlet vertex
                        }
                      if (flag) elem_fluide(som) = -1.0; // for Mean Gradient
                    }
                }
              fact  -= c_prepro_;
            }

          ///////////////////////////////////////////////////////
          //compute fluid projection for the neighbors of som
          ///////////////////////////////////////////////////////
          IntList& voisins = sommets_voisins_[som];
          int taille = 0;
          if (!(voisins.est_vide())) taille = voisins.size();
          // Cerr<<"taille : "<<taille<<finl;
          // Cerr<<"dirichlet : "<<dirichlet(voisins[0])<<finl;
          for (int l=0; l<taille; l++)
            {
              int num_som_v= voisins[l];
              if ( dirichlet(num_som_v)<0 )
                {
                  x = coordsDom3D(num_som_v, 0);
                  y = coordsDom3D(num_som_v, 1);
                  z = coordsDom3D(num_som_v, 2);
                  xp = solid_points_->valeurs()(num_som_v, 0);
                  yp = solid_points_->valeurs()(num_som_v, 1);
                  zp = solid_points_->valeurs()(num_som_v, 2);
                  norme = sqrt((x - xp)*(x - xp) + (y - yp)*(y - yp) + (z - zp)*(z - zp));

                  if (norme < eps_effec_)
                    {
                      elem_fluide(num_som_v) = -3; //O(h) because the projection and node are located at the same position (i.e. the node is already on the immersed boundary)
                      for (int d = 0; d < dim_esp; d++)
                        {
                          fluide(num_som_v, d) = coordsDom3D(som, d);
                          normale_som(num_som_v, d) = 0.0;
                        }
                    }
                  else
                    {
                      for (int d = 0; d < dim_esp; d++)
                        {
                          normale_som(num_som_v, d) = (coordsDom3D(num_som_v, d) - solid_points_->valeurs()(num_som_v, d)) / norme;
                          fluide(num_som_v, d) = dimTab_(d) ? coordsDom3D(num_som_v, d) + fact * hmax_node(num_som_v) * normale_som(num_som_v, d) : coordsDom3D(num_som_v, d);
                        }
                      interFlElem = le_dom.chercher_elements(fluide(num_som_v,0),fluide(num_som_v,1),fluide(num_som_v,2));
                      if (interFlElem == -1) // element not found
                        elem_fluide(num_som_v) = -3.;
                      else
                        {
                          elem_fluide(num_som_v) = float(interFlElem);
                          if (aire(interFlElem) > 0.0)
                            elem_fluide(num_som_v) = -1.; // O(h) element traversed by IBC
                          else
                            {
                              bool flag = false;
                              for (int j = 0; j < nb_som_elem; j++)
                                {
                                  int s = connect(interFlElem, j);
                                  if (dirichlet(s) == 1.0) flag = true; // element with a Dirichlet vertex
                                }
                              if (flag) elem_fluide(num_som_v) = -1.0; // for Mean Gradient
                            }
                        }
                    }

                  // If no fluid element was found, extend the zone to neighboring elements taking into account the c_prepro_ coefficient
                  if (c_prepro_ > 0 && elem_fluide(num_som_v) == -1)
                    {
                      fact  += c_prepro_;
                      for (int d = 0; d < dim_esp; d++)
                        fluide(num_som_v, d) = dimTab_(d) ? coordsDom3D(num_som_v, d) + fact * hmax_node(num_som_v) * normale_som(num_som_v, d) : coordsDom3D(num_som_v, d);
                      interFlElem = le_dom.chercher_elements(fluide(num_som_v,0),fluide(num_som_v,1),fluide(num_som_v,2));
                      if (interFlElem == -1) // element not found
                        elem_fluide(num_som_v) = -3.;
                      else
                        {
                          elem_fluide(num_som_v) = float(interFlElem);
                          if (aire(interFlElem) > 0.0)
                            elem_fluide(num_som_v) = -1.; // O(h) element traversed by IBC
                          else
                            {
                              bool flag = false;
                              for (int j = 0; j < nb_som_elem; j++)
                                {
                                  int s = connect(interFlElem, j);
                                  if (dirichlet(s) == 1.0) flag = true; // element with a Dirichlet vertex
                                }
                              if (flag) elem_fluide(num_som_v) = -1.0; // for Mean Gradient
                            }
                        }
                      fact  -= c_prepro_;
                    }
                }
            }
        }
    }

  // Update the IBC rotation matrix using nodal normals
  // ROTATION MATRIX

  Static_Int_Lists connectivite_som_elem;
  construire_connectivite_som_elem(nb_sommets_tot, connect, connectivite_som_elem, 1 /* include virtual elements */);

  DoubleTrav normalWeighArray(nb_elem);
  for (int som = 0; som < nb_som; som++)
    {
      double n0 = normale_som(som,0);
      double n1 = normale_som(som,1);
      double n2 = normale_som(som,2);
      norme = sqrt(n0*n0+n1*n1+n2*n2);
      if (norme > eps_effec_)
        {
          int sz = connectivite_som_elem.get_list_size(som);
          for (int j = 0; j <sz ; j++)
            {
              int numElem = connectivite_som_elem(som,j);

              if (aire(numElem)  < eps_effec_)
                {
                  normalArray(numElem,0) += n0;
                  normalArray(numElem,1) += n1;
                  normalArray(numElem,2) += n2;
                  normalWeighArray(numElem) += 1.0;
                }
            }
        }
    }

  DoubleTrav t1Arr(nb_elem, dim_esp), t2Arr(nb_elem, dim_esp);
  for (int e=0; e<nb_elem; e++)
    {
      if (normalWeighArray(e) > eps_effec_)
        {
          norme = sqrt(normalArray(e,0)*normalArray(e,0)+normalArray(e,1)*normalArray(e,1)+normalArray(e,2)*normalArray(e,2));
          if (norme > eps_effec_)
            {
              for (int d = 0; d < dim_esp; d++) normalArray(e,d) /= norme;
              // compute the rotation matrix
              computeLocalFrame(normalArray, t1Arr, t2Arr, e);
              computeMatRot(normalArray, t1Arr, t2Arr, e);
            }
        }
    }

}

