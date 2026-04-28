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

#include <Prepro_IBM_base.h>
#include <Sous_Domaine.h>
#include <Domaine_VF.h>
#include <Ecrire_MED.h>
#include <LireMED.h>
#include <Domaine.h>
#include <stdexcept>
#include <Faces.h>

Implemente_base(Prepro_IBM_base, "Prepro_IBM_base", Objet_U);
// XD Prepro_IBM_base Objet_U Prepro_IBM_base -1 To perform the intersection of an IB (Lagrange mesh) in a MED-format file .med with the Euler computional mesh.

Entree& Prepro_IBM_base::readOn(Entree& s)
{
  Cout << "Prepro_IBM => " << finl;
  Param param(que_suis_je());
  set_param(param);
  Cout<<"constante_prepro_c_IBM = "<<c_prepro_<<endl;
  return s;
}

void Prepro_IBM_base::set_param(Param& param) const
{
  param.ajouter("epsilon_prepro_IBM",&eps_,Param::OPTIONAL); // XD_ADD_P double geometric precision (<<1)
  param.ajouter("constant_c_IBM",&c_prepro_,Param::OPTIONAL);  // XD_ADD_P double additive coefficient to search the purely fluid point (cf. publications G Billo)
  param.ajouter_non_std("directions_pt_fluid",(this),Param::OPTIONAL); // XD_ADD_P listentier corresponding to each direction to search the purely fluid point
  // Exemples :
  //  *  [True, True, False] if the case is 2D in plane XY
  //  *  [False, True, True] if the case is 2D in plane YZ
  //  *  [True, True, True] if the case is 3D
  param.ajouter_non_std("MESH_Lagrange_file",(this), Param::OPTIONAL); // XD_ADD_P chaine Name of the .med file including the IB surfacic mesh.
  param.ajouter("MESH_Lagrange_name",&nom_maillage_IB_, Param::OPTIONAL); // XD_ADD_P chaine Name of the IB surfacic mesh.
  param.ajouter_flag("write_results_prepro",&save_prepro_); // XD_ADD_P flag to save output from prepro IBM
  param.ajouter("Out_MED_file_name",&nom_fichier_med_Out_, Param::OPTIONAL); // XD_ADD_P chaine Name of the .med file to save output from prepro IBM

  param.ajouter("prepro_IBM_verbose",&verbose_,Param::OPTIONAL); // XD_ADD_P entier to get verbose version
  param.ajouter_flag("verify_results_prepro",&verify_results_prepro_); // XD_ADD_P flag to save output from prepro IBM
}

int Prepro_IBM_base::lire_motcle_non_standard(const Motcle& un_mot, Entree& is)
{
  if (un_mot == "directions_pt_fluid")
    {
      Cout << "reading search directions for reference fluid pt ... " << finl;
      int dim_lu;
      is >> dim_lu;
      int dim_geom = Objet_U::dimension;
      if(dim_lu != dim_geom)
        {
          Cerr<<"Prepro_IBM : dim_lu <> dim_geom = "<<dim_geom<<endl;
          Process::exit();
        }
      dimTab_.resize(dim_lu);
      for (int i = 0; i < dim_lu; i++) is >> dimTab_(i);
      for (int i = 0; i < dim_lu; i++)
        {
          if (i == 0) Cout<<"directions : "<<dimTab_(0);
          else Cout<<" "<<dimTab_(i);
          if (i == (dim_lu-1)) Cout<<endl;
        }
    }

  if (un_mot == "MESH_Lagrange_file")
    {
      is >> nom_fichier_med_IB_;
      Cout << "reading MED IB mesh ... " <<nom_fichier_med_IB_ <<finl;
      if (nom_maillage_IB_ == "??")
        {
          Cerr<<"Name of the IB Lagrange mesh is required. Please, give it using the keyword MESH_Lagrange_name option before MESH_Lagrange_file."<< endl;
          Process::exit();
        }

      LireMED liremed(nom_fichier_med_IB_, nom_maillage_IB_);
      Nom nom_dom_ = "dom_IB";
      dom_med_IB_.nommer(nom_dom_);
      liremed.associer_domaine(dom_med_IB_);
      liremed.retrieve_MC_objects();
      aSkinUMesh_ = liremed.get_mc_mesh();

      int space_dim = aSkinUMesh_->getSpaceDimension();
      if (space_dim != Objet_U::dimension)
        {
          Cerr<<"Prepro_IBM: MC space dimension = " << space_dim << " is different from Trust space dimension = "<< Objet_U::dimension << endl;
          Process::exit();
        }
      int nbElemSur = int(aSkinUMesh_->getNumberOfCells());
      int nbNodeSur = int(aSkinUMesh_->getNumberOfNodes());
      Cerr << "Detecting " << nbNodeSur << " nodes and " << nbElemSur << " cells." << finl;
      const double *coord = aSkinUMesh_->getCoords()->begin();
      coordsSur3D_.resize(nbNodeSur, space_dim);
      std::copy(coord, coord+coordsSur3D_.size_array(), coordsSur3D_.addr());

      const double *normal = aSkinUMesh_->buildOrthogonalField()->getArray()->begin();
      normalArr_.resize(nbElemSur, space_dim);
      std::copy(normal, normal+normalArr_.size_array(), normalArr_.addr());

      const double *bary = aSkinUMesh_->computeIsoBarycenterOfNodesPerCell()->begin();
      barySurf_.resize(nbElemSur, space_dim);
      std::copy(bary, bary+barySurf_.size_array(), barySurf_.addr());

      if (verbose_)
        {
          Cout << "=== Test : Affichage des 8 premiers sommets ===" << finl;
          Cout << "--- Coordonnées des sommets (max 8) ---" << finl;
          for (int i = 0; i < std::min(8, coordsSur3D_.dimension(0)); i++)
            {
              Cout << "Sommet " << i << " : ("
                   << coordsSur3D_(i, 0) << ", " << coordsSur3D_(i, 1);
              if (space_dim == 3) Cout << ", " << coordsSur3D_(i, 2);
              Cout << ")" << finl;
            }
          Cout << "=== Test : Affichage des 2 premiers barycentres/normales ===" << finl;
          Cout << "Nombre d'éléments surfaciques : " << nbElemSur << finl;
          Cout << "--- Coordonnées des barycentres/normales (max 2) ---" << finl;
          for (int i = 0; i < std::min(2, barySurf_.dimension(0)); i++)
            {
              Cout << "Barycentres " << i << " : ("
                   << barySurf_(i, 0) << ", " << barySurf_(i, 1);
              if (space_dim == 3) Cout << ", " << barySurf_(i, 2);
              Cout << ")" << finl;
            }
          for (int i = 0; i < std::min(2, normalArr_.dimension(0)); i++)
            {
              Cout << "Normales " << i << " : ("
                   << normalArr_(i, 0) << ", " << normalArr_(i, 1);
              if (space_dim == 3) Cout << ", " << normalArr_(i, 2);
              Cout << ")" << finl;
            }
        }
    }

  return 1;
}

Sortie& Prepro_IBM_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

void Prepro_IBM_base::associer_pb(const Probleme_base& pb)
{
  mon_pb_=pb;
  if(!mon_pb_.non_nul())
    {
      Cerr << "Error : " << que_suis_je() << "not associated to a problem ! " << finl;
      Process::exit();
    }
  discretiser();
}

void Prepro_IBM_base::discretiser()
{
  const Probleme_base& pb = mon_pb_.valeur();
  const Domaine_dis_base& le_dom_dis = pb.domaine_dis();
  int dim_esp = Objet_U::dimension;

  // Rotation
  assert(dim_esp==3);
  int nb_comp=dim_esp*dim_esp;
  Noms nom_c(nb_comp);
  Noms unites(nb_comp);
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nom_c,unites,nb_comp,0.,champ_rotation_);
  DoubleTab& rotationArray = champ_rotation_->valeurs();
  rotationArray = 0.;

  // Aire
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,"aire","m^2",1,0., champ_aire_);
  DoubleTab& aireArray = champ_aire_->valeurs();
  aireArray = 0.;

  // Barycentres
  nb_comp=dim_esp;
  Noms nom_bary(nb_comp);
  Noms unites_bary(nb_comp);
  nom_bary[0] = "X";
  nom_bary[1] = "Y";
  nom_bary[2] = "Z";
  unites_bary[0] = "m";
  unites_bary[1] = "m";
  unites_bary[2] = "m";
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nom_bary,unites_bary,nb_comp,0., champ_bary_);
  DoubleTab& baryArray = champ_bary_->valeurs();
  baryArray = 0.;

  // Normales
  nb_comp=dim_esp;
  Noms nom_c3(nb_comp);
  Noms unites3(nb_comp);
  nom_c3[0] = "nx";
  nom_c3[1] = "ny";
  nom_c3[2] = "nz";
  unites3[0] = "m";
  unites3[1] = "m";
  unites3[2] = "m";
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,vectoriel,nom_c3,unites3,nb_comp,0., champ_normal_);
  DoubleTab& normalArray = champ_normal_->valeurs();
  normalArray = 0.;

  // IsNodeDirichlet
  pb.discretisation().discretiser_champ("champ_sommets",le_dom_dis,"is_node_dirichlet"," ",1,0., isNodeDirichlet_);
  DoubleTab& isNodeDirichletArray = isNodeDirichlet_->valeurs();
  isNodeDirichletArray = -1.;

  // h_max_elem & h_max_node
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,"h_max_elem"," ",1,0., h_max_elem_);
  DoubleTab& h_max_elem = h_max_elem_->valeurs();
  h_max_elem = 0.;
  pb.discretisation().discretiser_champ("champ_sommets",le_dom_dis,"h_max_node"," ",1,0., h_max_node_);
  DoubleTab& h_max_node = h_max_node_->valeurs();
  h_max_node = 0.;

  // coord projection solide
  nb_comp=dim_esp;
  Noms nom_c4(nb_comp);
  nom_c4[0] = "X";
  nom_c4[1] = "Y";
  nom_c4[2] = "Z";
  pb.discretisation().discretiser_champ("champ_sommets", le_dom_dis,vectoriel,nom_c4,unites3,nb_comp,0., solid_points_);
  DoubleTab&  solideArray = solid_points_->valeurs();
  solideArray = 0.;

  // elem number projection solide
  pb.discretisation().discretiser_champ("champ_sommets",le_dom_dis,"so_elem_number","",1,0., solid_elems_);
  DoubleTab& solid_elemsArray = solid_elems_->valeurs();
  solid_elemsArray = -2.0;

  // coord projection fluide
  nb_comp=dim_esp;
  Noms nom_c5(nb_comp);
  nom_c5[0] = "X";
  nom_c5[1] = "Y";
  nom_c5[2] = "Z";
  pb.discretisation().discretiser_champ("champ_sommets", le_dom_dis,vectoriel,nom_c5,unites3,nb_comp,0., fluid_points_);
  DoubleTab&  fluidArray = fluid_points_->valeurs();
  fluidArray = 0.;

  // elem number projection fluide
  pb.discretisation().discretiser_champ("champ_sommets",le_dom_dis,"fl_elem_number","",1,0., fluid_elems_);
  DoubleTab& fluid_elemsArray = fluid_elems_->valeurs();
  fluid_elemsArray = -2.0;

  // correspondance elems
  pb.discretisation().discretiser_champ("champ_elem",le_dom_dis,"TRUST_SALOME_correspondance","",1,0., corresp_elems_);
  DoubleTab& corresp_elemsArray = corresp_elems_->valeurs();
  int nbElemVol = corresp_elemsArray.dimension(0);  // Nombre de cellules du maillage volumique
  for (int elem = 0; elem <nbElemVol; elem++) corresp_elemsArray(elem) = elem;

  // computing effective error :
  compute_effective_error();

  // Intersection Euler/Lagrange si lecture du maillage Lagrangien
  if (barySurf_.dimension(0)) computeAire2();
}

void Prepro_IBM_base::computeLocalFrame(const DoubleTab& normal, DoubleTab& t1, DoubleTab& t2, int e)
{
  int dim_esp = Objet_U::dimension;
  assert(dim_esp==dimTab_.dimension(0));
  assert((t1.dimension(1)==dim_esp)&&(t2.dimension(1)==dim_esp));
  if (dimTab_(0) && dimTab_(1) && dimTab_(2))
    {
      t1(e,0) = normal(e,1)+normal(e,2);
      t1(e,1) = -normal(e,0)+normal(e,2);
      t1(e,2) = -(normal(e,0)+normal(e,1));
    }
  else if (dimTab_(0) && dimTab_(1) && !(dimTab_(2)))
    {
      t1(e,0) = normal(e,1);
      t1(e,1) = -normal(e,0);
      t1(e,2) = 0.;
    }
  else if (dimTab_(0) && !(dimTab_(1)) && dimTab_(2))
    {
      t1(e,0) = normal(e,2);
      t1(e,1) = 0.;
      t1(e,2) = -normal(e,0);
    }
  else if (!(dimTab_(0)) && dimTab_(1) && dimTab_(2))
    {
      t1(e,0) = 0.;
      t1(e,1) = normal(e,2);
      t1(e,2) = -normal(e,1);
    }
  else if (dimTab_(0) && !(dimTab_(1)) && !(dimTab_(2)))
    {
      t1(e,0) = 0.;
      t1(e,1) = 1.;
      t1(e,2) = 0.;
    }
  else if (!(dimTab_(0)) && dimTab_(1) && !(dimTab_(2)))
    {
      t1(e,0) = 0.;
      t1(e,1) = 0.;
      t1(e,2) = 1.;
    }
  else if (!(dimTab_(0)) && !(dimTab_(1)) && dimTab_(2))
    {
      t1(e,0) = 1.;
      t1(e,1) = 0.;
      t1(e,2) = 0.;
    }
  double nt = 0.;
  for (int c= 0; c < dim_esp; c++) { nt += t1(e,c)*t1(e,c); }
  if (sqrt(nt) > 1.0e-12) t1 /= sqrt(nt) ;
  t2(e,0) = t1(e,1)*normal(e,2)-t1(e,2)*normal(e,1);
  t2(e,1) = t1(e,2)*normal(e,0)-t1(e,0)*normal(e,2);
  t2(e,2) = t1(e,0)*normal(e,1)-t1(e,1)*normal(e,0);
  nt = 0.;
  for (int c= 0; c < dim_esp; c++) { nt += t2(e,c)*t2(e,c); }
  if (sqrt(nt) > 1.0e-12) t2 /= sqrt(nt) ;
}

void Prepro_IBM_base::computeMatRot(const DoubleTab& normal, DoubleTab& t1, DoubleTab& t2, int elem) // calcul de la matrice de rotation
{
  DoubleTab& rotationArray = champ_rotation_->valeurs();
  assert(rotationArray.dimension(1) == Objet_U::dimension * Objet_U::dimension);

  rotationArray(elem, 0) = t2(elem,0);
  rotationArray(elem, 1) = t1(elem,0);
  rotationArray(elem, 2) = normal(elem,0) ;
  rotationArray(elem, 3) = t2(elem,1);
  rotationArray(elem, 4) = t1(elem,1);
  rotationArray(elem, 5) = normal(elem,1);
  rotationArray(elem, 6) = t2(elem,2);
  rotationArray(elem, 7) = t1(elem,2);
  rotationArray(elem, 8) = normal(elem,2);
}

void Prepro_IBM_base::computeAire2()
{
  // PARAM
  int idebug = verbose_;

  const Domaine_dis_base& le_dom_dis = mon_pb_.valeur().domaine_dis();
  const Domaine& le_dom = mon_pb_.valeur().domaine();
  const MEDCouplingUMesh* le_mc_mesh = le_dom.get_mc_mesh();

  DoubleTab& normalArray = champ_normal_->valeurs();
  DoubleTab& aireArray = champ_aire_->valeurs();
  DoubleTab& isNodeDirichletArray = isNodeDirichlet_->valeurs();
  DoubleTab& baryArray = champ_bary_->valeurs();

  int nbPtsDom =  le_dom_dis.domaine().nb_som();// Nombre de noeuds du maillage volumique
  int nbElemVol = aireArray.dimension(0);  // Nombre de cellules du maillage volumique
  assert(nbElemVol==le_mc_mesh->getNumberOfCells());
  int nbElemSur = normalArr_.dimension(0); // Nombre de cellules (faces) du maillage surfacique
  int dim_esp = Objet_U::dimension; // dimension Euler
  assert(normalArr_.dimension(1)==dim_esp);

  DoubleTrav t1Arr(nbElemSur, dim_esp), t2Arr(nbElemSur, dim_esp);
  Sous_Domaine sdom_mesh3D_loc;
  sdom_mesh3D_loc.associer_domaine(le_dom);
  std::vector<mcIdType> numNodes, numNodes2D, numNodes3D;
  ArrOfDouble mesh2DBBox(2*dim_esp);
  Octree_Double octree_mesh3D;
  octree_mesh3D.build_nodes(le_dom.les_sommets(), 0, eps_); //ne pas inclure les sommets virtuels
  MEDCouplingFieldDouble *measure_aSkinUMesh = aSkinUMesh_->getMeasureField(true);
  assert(measure_aSkinUMesh->getNumberOfComponents()==1);
  double mesure_tot_aSkinUMesh = measure_aSkinUMesh->accumulate(0);

  Cerr << "<<<<<<<<<<<<<<<< Prepro_IBM: BEGINNING OF AREA COMPUTATION..." << finl;
  Cerr << "Prepro_IBM: Computing intersection... with " << nbElemSur << " elements (Lagrangien) and " << nbElemVol << " / " << nbPtsDom << " elements/nodes (Eulerian)"  << finl;
  if (idebug) Cerr << "Domaine 3D : "<< le_dom.nb_som() << " nodes and "<< le_dom.nb_som_elem() <<" par element (Eulerian)"  << finl;

  Cerr<<"dimTab_ ("<<dimTab_(0)<<" , "<<dimTab_(1)<<" , "<<dimTab_(2)<<")"<<finl;

  aireArray = 0.;
  normalArray = 0.;
  baryArray = 0.;
  isNodeDirichletArray = -1.;
  double aire_tot_calcul = 0.;
  double mesure_sig_aSkinUMesh = 0.;
  /////////////////////////////////////////////////////////////////////
  // Boucle sur les elements du maillage surfacique dans l'espace 3D
  /////////////////////////////////////////////////////////////////////
  for (int e = 0; e < nbElemSur; e++)
    {
      Cerr<<" "<<finl;
      Cerr<<">>>> element de surface numero : "<< e <<finl;
      numNodes.clear();
      aSkinUMesh_->getNodeIdsOfCell(e, numNodes);
      int nbNodesSur = int(numNodes.size());
      Cerr<<"nbNodesSur = "<< nbNodesSur <<finl;

      // definition du repere local (normalArr_, t1Arr, t2Arr)
      Cerr << "normalArr_ (" << normalArr_(e,0) ;
      for (int cc = 1; cc <dim_esp; cc++) Cerr << " ,"<<normalArr_(e,cc);
      Cerr << ")" << finl;
      computeLocalFrame(normalArr_, t1Arr, t2Arr, e);
      Cerr << "t1Arr[e]_ (" << t1Arr(e,0);
      for (int cc = 1; cc <dim_esp; cc++) Cerr << " ," << t1Arr(e,cc);
      Cerr << " )" << finl;
      Cerr << "t2Arr[e]_ (" << t2Arr(e,0);
      for (int cc = 1; cc <dim_esp; cc++) Cerr << " ," << t2Arr(e,cc);
      Cerr << ")" << finl;
      if (idebug) Cerr<<"barySurf[e]_ ("<<barySurf_(e, 0)<<" , "<<barySurf_(e, 1)<<" , "<<barySurf_(e, 2)<<")"<<finl;

      // Mesh 2D de la facette 3D numero e
      MEDCouplingUMesh * mesh2DSurf = MEDCouplingUMesh::New("surf_"+std::to_string(e),2);
      DoubleTrav surfCoords2D(nbNodesSur, 2);
      for (int node =0; node < nbNodesSur; node++)
        {
          int numNode = int(numNodes[node]);
          surfCoords2D(node,0) = (coordsSur3D_(numNode,0)-barySurf_(e,0))*t1Arr(e,0) + (coordsSur3D_(numNode,1)-barySurf_(e,1))*t1Arr(e,1) + (coordsSur3D_(numNode,2)-barySurf_(e,2))*t1Arr(e,2);
          surfCoords2D(node,1) = (coordsSur3D_(numNode,0)-barySurf_(e,0))*t2Arr(e,0) + (coordsSur3D_(numNode,1)-barySurf_(e,1))*t2Arr(e,1) + (coordsSur3D_(numNode,2)-barySurf_(e,2))*t2Arr(e,2);
        }
      MCAuto<MEDCoupling::DataArrayDouble> array(MEDCoupling::DataArrayDouble::New());
      array->useArray(surfCoords2D.addr(), false, MEDCoupling::DeallocType::CPP_DEALLOC, nbNodesSur, 2);
      mesh2DSurf->setCoords(array);
      array->decrRef();
      mesh2DSurf->allocateCells(1);
      IntTrav mesh2DSurf_conect(1, nbNodesSur);
      for (int node =0; node < nbNodesSur; node++) mesh2DSurf_conect(0,node) = node;
      auto ptr = mesh2DSurf_conect.addr();
      std::vector<mcIdType> tmp(ptr, ptr+nbNodesSur);
      INTERP_KERNEL::NormalizedCellType typ_fac_mc = INTERP_KERNEL::NORM_POLYGON;
      mesh2DSurf->insertNextCell(typ_fac_mc, nbNodesSur, tmp.data());

      if (idebug)
        {
          numNodes2D.clear();
          mesh2DSurf->getNodeIdsOfCell(0, numNodes2D);
          const double *mesh2DSurf_coords = mesh2DSurf->getCoords()->begin();
          DoubleTab Surf2DCoordonnes;
          Surf2DCoordonnes.resize(int(mesh2DSurf->getNumberOfNodes()), 2);
          std::copy(mesh2DSurf_coords, mesh2DSurf_coords+Surf2DCoordonnes.size_array(), Surf2DCoordonnes.addr());
          Cerr << "mesh2DSurf : "<<mesh2DSurf->getNumberOfCells()<<" cells and "<<mesh2DSurf->getNumberOfNodes()<< " nodes in 2D space"<<finl;
          for (int n = 0; n <mesh2DSurf->getNumberOfNodes(); n++)
            {
              for (int cc = 0; cc <2; cc++) Cerr <<Surf2DCoordonnes(n, cc) <<" ";
              Cerr << endl;
            }
          Cerr << "mesh2DSurf : nb nodes cell #1 = "<<int(numNodes2D.size())<<finl;
          for (int n = 0; n <int(numNodes2D.size()); n++)
            {
              Cerr <<int(numNodes2D[n]) <<" ";
            }
          Cerr << endl;
        }

      // Calcul de la Bounding Box 3D du polygone 2D
      for (int c =0; c < dim_esp; c++)
        {
          int index0 = 2*c;
          int index1 = index0 + 1;
          mesh2DBBox(index0) = 1.e+30; // min
          mesh2DBBox(index1) = -1.e+30; // max
          for (int k =0; k < nbNodesSur; k++)
            {
              int numNode = int(numNodes[k]);
              if (coordsSur3D_(numNode,c) <= mesh2DBBox(index0)) mesh2DBBox(index0) = coordsSur3D_(numNode,c);
              if (coordsSur3D_(numNode,c) >= mesh2DBBox(index1)) mesh2DBBox(index1) = coordsSur3D_(numNode,c);
              // Cerr<<"coordsSur3D_("<<numNode<<","<<c<<") = "<<coordsSur3D_(numNode,c)<<finl;
            }
          mesh2DBBox(index0) -= eps_effec_;
          mesh2DBBox(index1) += eps_effec_;
        }
      Cerr<<"mesh2DBBox = ";
      for (int cc = 0; cc <dim_esp; cc++) Cerr<<"["<<mesh2DBBox(2*cc)<<", "<<mesh2DBBox(2*cc+1)<<"] ";
      Cerr<<finl;

      // Détection des éléments volumiques dans la Bounding Box
      const double bbox[] = {mesh2DBBox(0), mesh2DBBox(1), mesh2DBBox(2), mesh2DBBox(3), mesh2DBBox(4), mesh2DBBox(5)};
      MEDCoupling::DataArrayIdType * cellIdsArr = le_mc_mesh->getCellsInBoundingBox(bbox, 0.); //precision = 0 pour ne pas avoir d'element 3D non coupe par la frontiere
      const mcIdType *daP = cellIdsArr->begin();
      int NbCellsInBB = int(cellIdsArr->getNumberOfTuples());
      int NbCompsInBB = int(cellIdsArr->getNumberOfComponents());
      Cerr<<"getCellsInBoundingBox cellIdsArr : Nb Cells, Nb comp = "<< NbCellsInBB << " "<< NbCompsInBB<< finl;
      Cerr<<" " ;
      for (int k = 0; k <NbCellsInBB ; k++) Cerr<< (int)daP[k] <<" " ;
      Cerr<< finl;


      double aire_elem_calcul = 0.;
      ////////////////////////////////////////////////////////////////////////////
      // Calcul de l'intersection volume-surface pour chaque element intercepte
      ////////////////////////////////////////////////////////////////////////////
      for (int k = 0; k < NbCellsInBB; k++)
        {
          int my_elem = (int)daP[k];
          if (idebug) Cerr<<"**** num elem = "<<my_elem<<" ****"<<finl ;

          // maillage 3D Euler MC 1 cell
          const mcIdType cellIds[1]= {my_elem};
          MEDCouplingUMesh *mesh3D = le_mc_mesh->buildPartOfMySelf(cellIds,cellIds+1,false);
          assert(mesh3D->getNumberOfCells()==1);
          int NumberOfNodes3D = int(mesh3D->getNumberOfNodes());

          // verification element mesh3D intercepte par plan infini
          const double *crd3D = mesh3D->getCoords()->begin();
          DoubleTab coord3D;
          coord3D.resize(NumberOfNodes3D, dim_esp);
          std::copy(crd3D, crd3D+coord3D.size_array(), coord3D.addr());
          double signmult0;
          int istop = 1;
          for (int node3D = 0; node3D <NumberOfNodes3D ; node3D++)
            {
              double ps3D = 0.;
              for (int cc = 0; cc <dim_esp; cc++) ps3D += (coord3D(node3D,cc) - barySurf_(e,cc)) * normalArr_(e,cc);
              if (node3D == 0) {signmult0 = ps3D;}
              else if (signmult0 * ps3D < 0.) istop = 0;
            }
          if (istop)
            {
              if (idebug) Cerr<<"Element not intercepted " <<finl;
              continue;
            }

          //  on calcule l'intersection entre l'element volumique mesh3D et l'element surfacique "e"
          const double origin[] = {barySurf_(e,0), barySurf_(e,1), barySurf_(e,2)};
          const double vec[] = {normalArr_(e,0), normalArr_(e,1), normalArr_(e,2)};
          if (idebug)
            {
              Cerr<<"origin element surfacique "<<e<<" = ";
              for (int cc = 0; cc <dim_esp; cc++) Cerr<< origin[cc]<<" ";
              Cerr << endl;
              Cerr<<"vec normal element surfacique "<<e<<" = ";
              for (int cc = 0; cc <dim_esp; cc++) Cerr<< vec[cc]<<" ";
              Cerr << endl;
            }
          MEDCoupling::DataArrayIdType * polycellIds ;
          try
            {
              // Intersection par un plan infini
              MEDCouplingUMesh * interPoly3D = mesh3D->buildSlice3D(origin, vec, 0., polycellIds);
              polycellIds->decrRef();
              int nbCellsPoly3D = int(interPoly3D->getNumberOfCells());
              int NumberOfNodes = int(interPoly3D->getNumberOfNodes());
              const double *interCoords3D = interPoly3D->getCoords()->begin();
              DoubleTrav interCoordonnes3D;
              interCoordonnes3D.resize(NumberOfNodes, dim_esp);
              std::copy(interCoords3D, interCoords3D+interCoordonnes3D.size_array(), interCoordonnes3D.addr());

              // Mesh 2D du polygone d'intersection en 3D interPoly3D
              numNodes3D.clear();
              interPoly3D->getNodeIdsOfCell(0, numNodes3D);
              int nbNodesCell1 = int(numNodes3D.size());
              if (nbCellsPoly3D != 1)
                {
                  Cerr<<"<Nb cells> interPoly3D = "<<nbCellsPoly3D<<" <> 1. Exit"<<finl;
                  Process::exit();
                }
              // A priori, les noeuds de mesh3D sont repris dans interPoly3D (en debut de la liste de noeuds)
              if (idebug)
                {
                  Cerr << "<Nb cells> interPoly3D = "<<nbCellsPoly3D<<endl;
                  Cerr << "<Nb nodes> mesh3D, interPoly3D and Cell#1 = "<<NumberOfNodes3D<<" "<<NumberOfNodes<<" "<<nbNodesCell1<<endl;
                }
              assert((NumberOfNodes - NumberOfNodes3D) == nbNodesCell1);
              MEDCouplingUMesh * interPoly2D = MEDCouplingUMesh::New("interPoly2D_"+std::to_string(my_elem),2);
              DoubleTrav interPolyCoords2D(nbNodesCell1, (dim_esp - 1));
              for (int node = 0; node <nbNodesCell1 ; node++)
                {
                  int numNode = int(numNodes3D[node]);
                  interPolyCoords2D(node,0) = (interCoordonnes3D(numNode,0)-barySurf_(e,0))*t1Arr(e,0) + (interCoordonnes3D(numNode,1)-barySurf_(e,1))*t1Arr(e,1) + (interCoordonnes3D(numNode,2)-barySurf_(e,2))*t1Arr(e,2);
                  interPolyCoords2D(node,1) = (interCoordonnes3D(numNode,0)-barySurf_(e,0))*t2Arr(e,0) + (interCoordonnes3D(numNode,1)-barySurf_(e,1))*t2Arr(e,1) + (interCoordonnes3D(numNode,2)-barySurf_(e,2))*t2Arr(e,2);
                }
              MCAuto<MEDCoupling::DataArrayDouble> interPolyarray2D(MEDCoupling::DataArrayDouble::New());
              interPolyarray2D->useArray(interPolyCoords2D.addr(), false, MEDCoupling::DeallocType::CPP_DEALLOC, nbNodesCell1, 2);
              interPoly2D->setCoords(interPolyarray2D);
              interPolyarray2D->decrRef();
              interPoly2D->allocateCells(1);
              IntTrav Poly2D_conect(1, nbNodesCell1);
              for (int node =0; node < nbNodesCell1; node++) Poly2D_conect(0,node) = node;
              auto Poly2D_ptr = Poly2D_conect.addr();
              std::vector<mcIdType> Poly2D_tmp(Poly2D_ptr, Poly2D_ptr+ nbNodesCell1);
              typ_fac_mc = INTERP_KERNEL::NORM_POLYGON;
              interPoly2D->insertNextCell(typ_fac_mc, nbNodesCell1, Poly2D_tmp.data());

              if (idebug)
                {
                  numNodes2D.clear();
                  interPoly2D->getNodeIdsOfCell(0, numNodes2D);
                  const double *interPoly2D_coords = interPoly2D->getCoords()->begin();
                  DoubleTab interPoly2DCoordonnes;
                  interPoly2DCoordonnes.resize(int(interPoly2D->getNumberOfNodes()), (dim_esp-1));
                  std::copy(interPoly2D_coords, interPoly2D_coords+interPoly2DCoordonnes.size_array(), interPoly2DCoordonnes.addr());
                  Cerr << "interPoly2D : "<<interPoly2D->getNumberOfCells()<<" cells and "<<interPoly2D->getNumberOfNodes()<< " nodes in 2D space"<<finl;
                  for (int n = 0; n <interPoly2D->getNumberOfNodes(); n++)
                    {
                      for (int cc = 0; cc <(dim_esp-1); cc++) Cerr <<interPoly2DCoordonnes(n, cc) <<" ";
                      Cerr << endl;
                    }
                  Cerr << "interPoly2D : nb nodes cell #1 = "<<int(numNodes2D.size())<<finl;
                  for (int n = 0; n <int(numNodes2D.size()); n++)
                    {
                      Cerr <<int(numNodes2D[n]) <<" ";
                    }
                  Cerr << endl;
                }

              ////////////////////////////////////////////////////////////////////////////////////
              // Intersection des polygones 2D de la facette e et de la coupe de l'element my_elem
              ////////////////////////////////////////////////////////////////////////////////////
              int status_polypoly = 0;
              DoubleTab finalcoords;
              IntTab Tab_conect;
              intersectPolyPoly2D(interPoly2D, mesh2DSurf, finalcoords, Tab_conect, eps_effec_, status_polypoly);
              if (idebug && (status_polypoly != 0)) Cerr<<"status intersectPolyPoly2D : "<<status_polypoly<<finl;
              if (status_polypoly != 0) continue; // next cell

              MEDCouplingUMesh * finalMesh = MEDCouplingUMesh::New("output_mesh",(dim_esp-1));
              finalMesh->allocateCells(1);

              MCAuto<MEDCoupling::DataArrayDouble> outputCoords(MEDCoupling::DataArrayDouble::New());
              outputCoords->useArray(finalcoords.addr(), false, MEDCoupling::DeallocType::CPP_DEALLOC, finalcoords.dimension(0), finalcoords.dimension(1));
              finalMesh->setCoords(outputCoords);
              auto connect_ptr = Tab_conect.addr();
              std::vector<mcIdType> finalConnect(connect_ptr, connect_ptr+Tab_conect.dimension(1) );
              finalMesh->insertNextCell(typ_fac_mc, Tab_conect.dimension(1), finalConnect.data());
              finalMesh->finishInsertingCells();


              // Definition 3D des noeuds de finalmesh 2D
              DoubleTab interCoordonnes;
              interCoordonnes.resize(int(finalMesh->getNumberOfNodes()), dim_esp);
              for (int node =0; node < interCoordonnes.dimension(0); node++)
                for (int kd=0; kd<dim_esp; kd++)
                  interCoordonnes(node,kd) = barySurf_(e,kd) + finalcoords(node,0)*t1Arr(e,kd) + finalcoords(node,1)*t2Arr(e,kd);

              if (idebug)
                {
                  Cerr << "finalMesh : "<<int(finalMesh->getNumberOfCells())<<" cells and "<<int(finalMesh->getNumberOfNodes())<< " node in 2D space : "<<" ";
                  for (int cell=0; cell < finalMesh->getNumberOfCells(); cell++)
                    {
                      numNodes2D.clear();
                      finalMesh->getNodeIdsOfCell(cell, numNodes2D);
                      for (int n = 0; n <int(numNodes2D.size()); n++)
                        {
                          Cerr <<int(numNodes2D[n]) <<" ";
                        }
                      Cerr << endl;
                    }
                  for (int n = 0; n <finalMesh->getNumberOfNodes(); n++)
                    {
                      for (int cc = 0; cc <2; cc++) Cerr <<finalcoords(n, cc) <<" ";
                      Cerr<<" ( ";
                      for (int cc = 0; cc <dim_esp; cc++) Cerr <<interCoordonnes(n, cc) <<" ";
                      Cerr <<" ) "<< endl;
                    }
                }

              /////////////////////////////////////////////////////////////////////////
              // Exploitation des results d'intersection pour calcul barycentre et aire
              /////////////////////////////////////////////////////////////////////////
              int nbpts = interCoordonnes.dimension(0);
              double inv_nbpts =1.0 /( double(nbpts));
              DoubleTab baryMean(dim_esp);
              baryMean= 0.;
              for (int cc = 0; cc <dim_esp; cc++)
                for (int n = 0; n <nbpts; n++) baryMean(cc) += interCoordonnes(n, cc);
              baryMean *= inv_nbpts;
              if (idebug)
                {
                  Cerr << "finalMesh : mean barycenter = ";
                  for (int cc = 0; cc <dim_esp; cc++) Cerr << baryMean(cc) <<" ";
                  Cerr << finl;
                }
              MEDCouplingFieldDouble * measure = finalMesh->getMeasureField(true);
              assert(measure->getNumberOfValues()==1);
              double aire = measure->accumulate(0);
              if (idebug) Cerr<<"Element measure finalMesh = "<<aire<<finl;
              aire_elem_calcul += aire;

              aireArray(my_elem) += aire;
              for (int cc = 0; cc <dim_esp; cc++) baryArray(my_elem, cc) += aire * baryMean(cc);
              for (int cc = 0; cc <dim_esp; cc++) normalArray(my_elem, cc) += aire*normalArr_(e,cc);
              numNodes.clear();
              le_mc_mesh->getNodeIdsOfCell(my_elem, numNodes);
              int nbNodesElem = int(numNodes.size()), the_node = 0;
              for (int l = 0; l <nbNodesElem; l++)
                {
                  the_node = int(numNodes[l]);
                  isNodeDirichletArray(the_node) = 1.0 ;
                }

              if (idebug)
                {
                  Cerr<<"TRUST element aireArray measure = "<<aireArray(my_elem)<<finl;
                  Cerr<<"TRUST element aire weighted baryArray = ";
                  for (int cc = 0; cc <dim_esp; cc++) Cerr<<baryArray(my_elem,cc)<<" ";
                  Cerr<<finl;
                  Cerr<<"TRUST element aire weighted normalArray = ";
                  for (int cc = 0; cc <dim_esp; cc++) Cerr<<normalArray(my_elem,cc)<<" ";
                  Cerr<<finl;
                  Cerr<<"TRUST isNodeDirichletArray= 1 for nodes = ";
                  for (int l = 0; l <nbNodesElem; l++) Cerr<< numNodes[l] <<" ";
                  Cerr<<finl;
                }
            }
          catch (const std::runtime_error& err)
            {
              std::cerr << "Erreur : " << err.what() << std::endl;
            }
        }

      // Bilan AIre element e
      aire_tot_calcul += aire_elem_calcul ;
      double val1;
      measure_aSkinUMesh->getArray()->getTuple(e, &val1);
      Cerr<<"///////////// Measure of surface for element "<<e<<" = "<<val1<<finl;
      Cerr<<"///////////// Measure of computed surface for element "<<e<<" = "<<aire_elem_calcul <<finl;

      sdom_mesh3D_loc.les_elems().reset();
      mesure_sig_aSkinUMesh += val1;
    }

  // Barycentres et Normales moyens ponderes par les aires
  DoubleTrav t1EulerArr(nbElemVol, dim_esp), t2EulerArr(nbElemVol, dim_esp);
  for (int elem = 0; elem <nbElemVol; elem++)
    {
      if (aireArray(elem) >0.)
        {
          double norm = 0.;
          for (int cc = 0; cc <dim_esp; cc++)
            {
              baryArray(elem, cc) /= aireArray(elem) ;
              normalArray(elem, cc) /= aireArray(elem) ;
              norm += normalArray(elem, cc)*normalArray(elem, cc);
            }
          norm = sqrt(norm);
          for (int cc = 0; cc <dim_esp; cc++) normalArray(elem, cc) /= norm; // normalisation;
          // Matrice rotation
          computeLocalFrame(normalArray, t1EulerArr, t2EulerArr, elem);
          computeMatRot(normalArray, t1EulerArr, t2EulerArr, elem);
        }
    }

  // Bilan AIre totale
  Cerr<<" "<<finl;
  Cerr<<"Measure of total surface for Immersed Boundary = "<<mesure_tot_aSkinUMesh<<finl;
  Cerr<<"Measure of sum of elementary surfaces for Immersed Boundary = "<<mesure_sig_aSkinUMesh<<finl;
  Cerr<<"Measure of total computed surface = "<<aire_tot_calcul<<" >>>>>>>>>>>>>>>>"<<finl;

  aireArray.echange_espace_virtuel();
  baryArray.echange_espace_virtuel();
  normalArray.echange_espace_virtuel();
  DoubleTab& rotationArray = champ_rotation_->valeurs();
  rotationArray.echange_espace_virtuel();
  isNodeDirichletArray.echange_espace_virtuel();
}

void Prepro_IBM_base::intersectPolyPoly2D(MEDCouplingUMesh * polyMesh1, MEDCouplingUMesh * polyMesh2, DoubleTab& finalcoords, IntTab& Tab_conect, double eps, int& status)
{
  int idebug = verbose_;
  if (status != 0) status = 0;
  int dim_esp = Objet_U::dimension;
  MCAuto<MEDCoupling::DataArrayDouble> outputCoords(MEDCoupling::DataArrayDouble::New());
  MCAuto<MEDCoupling::DataArrayDouble> MC_p_test(MEDCoupling::DataArrayDouble::New());
  MCAuto<MEDCoupling::DataArrayDouble> dv_outputCoords(MEDCoupling::DataArrayDouble::New());

  // Edges polyMesh1
  const double *mc_coords1 = polyMesh1->getCoords()->begin();
  int nbNodes1 = int(polyMesh1->getNumberOfNodes());
  DoubleTrav coords1;
  coords1.resize(nbNodes1, (dim_esp-1));
  std::copy(mc_coords1, mc_coords1+coords1.size_array(), coords1.addr());
  MEDCouplingUMesh * polyEdgeMesh1 = polyMesh1->buildBoundaryMesh(false);
  int nbEdge1 = int(polyEdgeMesh1->getNumberOfCells());
  if (idebug)
    {
      Cerr<<">>> Prepro_IBM_base::intersectPolyPoly2D: Cell and node nb polyEdgeMesh1 = "<<int(polyEdgeMesh1->getNumberOfCells())<<" "<<int(polyEdgeMesh1->getNumberOfNodes())<<finl;
    }

  // Edges polyMesh2
  const double *mc_coords2 = polyMesh2->getCoords()->begin();
  int nbNodes2 = int(polyMesh2->getNumberOfNodes());
  DoubleTrav coords2;
  coords2.resize(nbNodes2, (dim_esp-1));
  std::copy(mc_coords2, mc_coords2+coords2.size_array(), coords2.addr());
  MEDCouplingUMesh * polyEdgeMesh2 = polyMesh2->buildBoundaryMesh(false);
  if (idebug)
    {
      Cerr<<">>> Prepro_IBM_base::intersectPolyPoly2D: Cell and node nb polyEdgeMesh2 = "<<int(polyEdgeMesh2->getNumberOfCells())<<" "<<int(polyEdgeMesh2->getNumberOfNodes())<<finl;
    }

  // Prendre tous les noeuds de polyMesh1 inclus dans polyMesh2
  DoubleTrav p_test(dim_esp-1);
  for (int i=0; i<nbNodes1; i++)
    {
      for (int k=0; k<(dim_esp-1); k++) p_test(k) = coords1(i,k);
      if (polyMesh2->getCellContainingPoint(p_test.addr(), eps) != -1)
        {
          MC_p_test->useArray(p_test.addr(), false, MEDCoupling::DeallocType::CPP_DEALLOC, 1, (dim_esp-1));
          if (outputCoords->isAllocated())
            {
              outputCoords->aggregate(MC_p_test);
            }
          else
            {
              outputCoords->alloc(1,(dim_esp-1));
              double * tmp = outputCoords->getPointer();
              std::copy(p_test.addr(), p_test.addr() + (dim_esp-1), tmp);
              outputCoords->declareAsNew();  // you have modified data pointed by internal pointer notify object
            }
        }
    }
  // Prendre tous les noeuds de polyMesh2 inclus dans polyMesh1
  for (int i=0; i<nbNodes2; i++)
    {
      for (int k=0; k<(dim_esp-1); k++) p_test(k) = coords2(i,k);
      if (polyMesh1->getCellContainingPoint(p_test.addr(), eps) != -1)
        {
          MC_p_test->useArray(p_test.addr(), false, MEDCoupling::DeallocType::CPP_DEALLOC, 1, (dim_esp-1));
          if (outputCoords->isAllocated())
            {
              outputCoords->aggregate(MC_p_test);
            }
          else
            {
              outputCoords->alloc(1,(dim_esp-1));
              double * tmp = outputCoords->getPointer();
              std::copy(p_test.addr(), p_test.addr() + (dim_esp-1), tmp);
              outputCoords->declareAsNew();  // you have modified data pointed by internal pointer notify object
            }
        }
    }
  if (!(outputCoords->isAllocated()))
    {
      status = -2;
      return;
    }

  // Prendre les points d'intersection
  for (int i=0; i<nbEdge1; i++)
    {
      std::vector< mcIdType > numNodes;
      polyEdgeMesh1->getNodeIdsOfCell(i, numNodes);
      if (numNodes.size() > 2)
        {
          Cerr<<">>>Prepro_IBM_base::intersectPolyPoly2D(): node nb = "<<numNodes.size()<<" for edge "<<i<<" => CHELOU..."<<finl;
          status = -1;
          return;
        }
      DoubleTrav point0(dim_esp-1), point1(dim_esp-1);
      for (int k=0; k<(dim_esp-1); k++) point0(k) = coords1(int(numNodes[0]),k);
      for (int k=0; k<(dim_esp-1); k++) point1(k) = coords1(int(numNodes[1]),k);
      if (idebug)
        {
          Cerr<<">>> Prepro_IBM_base::intersectPolyPoly2D(): polyEdgeMesh1: edge "<<i<<" nodes ";
          for (int ii=0; ii<int(numNodes.size()); ii++) Cerr<<numNodes[ii]<<" ";
          Cerr<<finl;
          Cerr<<"             point0: "<<point0(0)<<" "<<point0(1)<<finl;
          Cerr<<"             point1: "<<point1(0)<<" "<<point1(1)<<finl;
        }
      int status_segpoly = 0;
      intersectSegPoly2D(polyEdgeMesh2, point0, point1, eps, outputCoords, status_segpoly);
    }
  if (idebug) Cerr<<">>> Prepro_IBM_base::intersectPolyPoly2D(): outputCoords nodes after intersectSegPoly2D = "<<int(outputCoords->getNumberOfTuples())<<finl;

  MEDCouplingUMesh * outputMesh = MEDCouplingUMesh::New("output_mesh",(dim_esp-1));
  outputMesh->allocateCells(1);

  MEDCoupling::DataArrayIdType *connect = 0, *connInd = 0;
  outputCoords->findCommonTuples(eps, -1, connect, connInd);
  if (idebug) Cerr<<"             connect and connInd tuples after findCommonTuples = "<<int(connect->getNumberOfTuples())<<" "<<int(connInd->getNumberOfTuples())<<finl;
  if(int(connect->getNumberOfTuples()) > 0)
    {
      dv_outputCoords = outputCoords->getDifferentValues(eps);
      if (idebug) Cerr<<"             nodes nb after getDifferentValues = "<<int(dv_outputCoords->getNumberOfTuples())<<finl;
      outputMesh->setCoords(dv_outputCoords);
    }
  else outputMesh->setCoords(outputCoords);

  MCAuto<MEDCoupling::DataArrayDouble> finalCoords = outputMesh->getCoords();
  int finalCoords_len = int(finalCoords->getNumberOfTuples());
  if (idebug) Cerr<<"             Nb of tuples in finalCoords = "<<finalCoords_len<<finl;
  if (finalCoords_len < 3)
    {
      status = -3; //seulement une droite intersection
      return;
    }

  Tab_conect.resize(1, finalCoords_len);
  for (int node =0; node < finalCoords_len; node++) Tab_conect(0,node) = node;
  auto connect_ptr = Tab_conect.addr();
  std::vector<mcIdType> Connect(connect_ptr, connect_ptr+Tab_conect.dimension(1));
  INTERP_KERNEL::NormalizedCellType typ_fac_mc = INTERP_KERNEL::NORM_POLYGON;
  outputMesh->insertNextCell(typ_fac_mc, Tab_conect.dimension(1), Connect.data());

  // Butterfly ?
  // std::vector<mcIdType> cells;
  // cells.clear();
  // outputMesh->checkButterflyCells(cells, eps);
  // if( cells.size() > 0)
  //   {
  //     MEDCoupling::DataArrayIdType * modif_cell = outputMesh->convexEnvelop2D();
  //     if (modif_cell != nullptr)
  //       {
  //         int modNodNb = int(modif_cell->getNumberOfComponents());
  //         Cerr<<"             Modified cell connectivity by convexEnvelop2D for "<<modNodNb<<" nodes"<<finl;
  //       }
  //   }

  // Connectivite 2D
  std::vector<mcIdType> numNodes2D;
  numNodes2D.clear();

  // Test aire < eps ;
  MEDCouplingFieldDouble * measure = outputMesh->getMeasureField(true);
  double aire = measure->accumulate(0);
  if (aire<eps)
    {
      DoubleTab coords;
      coords.resize(int(outputMesh->getNumberOfNodes()), outputMesh->getMeshDimension());
      const double *mc_Coords = outputMesh->getCoords()->begin();
      std::copy(mc_Coords, mc_Coords+coords.size_array(), coords.addr());
      Cerr << "intersectPolyPoly2D: Element measure outputMesh = "<<aire<<" with tuples : " << endl;
      for (int n = 0; n <int(outputMesh->getNumberOfNodes()); n++)
        {
          for (int cc = 0; cc < outputMesh->getMeshDimension(); cc++) Cerr <<coords(n, cc) <<" ";
          Cerr << endl;
        }
      outputMesh->getNodeIdsOfCell(0, numNodes2D);
      for (int node =0; node <  int(numNodes2D.size()); node++) Cerr <<int(numNodes2D[node]) <<" ";
      Cerr << endl;

      if (int(outputMesh->getCoords()->getNumberOfTuples()) == 4)
        {
          // exchange conectivity 3<->2 of initial mesh
          const mcIdType old2newIds[] =  {0, 1, 3, 2};
          outputMesh->renumberNodesInConn(old2newIds);
          numNodes2D.clear();
          outputMesh->getNodeIdsOfCell(0, numNodes2D);
          aire = (outputMesh->getMeasureField(true))->accumulate(0);
          Cerr<<"             Element measure outputMesh = "<<aire<<" with modified cell connectivity => ";
          for (int node =0; node <  int(numNodes2D.size()); node++) Cerr <<int(numNodes2D[node]) <<" ";
          Cerr << endl;

          if (aire<eps)
            {
              // exchange conectivity 3<->0 of initial mesh
              const mcIdType old2newIds_bis[] =  {0, 1, 3, 2};
              outputMesh->renumberNodesInConn(old2newIds_bis);
              const mcIdType old2newIds_ter[] =  {3, 1, 2, 0};
              outputMesh->renumberNodesInConn(old2newIds_ter);
              numNodes2D.clear();
              outputMesh->getNodeIdsOfCell(0, numNodes2D);
              aire = (outputMesh->getMeasureField(true))->accumulate(0);
              Cerr<<"             Element measure outputMesh = "<<aire<<" with modified cell connectivity => ";
              for (int node =0; node <  int(numNodes2D.size()); node++) Cerr <<int(numNodes2D[node]) <<" ";
              Cerr << endl;

              if (aire<eps)
                {
                  Cerr<<"Element measure outputMesh = "<<aire;
                  Cerr<<" nothing to do :-( "<<finl;
                  status = -4;
                  // Process::exit();
                  return;
                }
            }
        }
      else
        {
          Cerr<<"             Element measure outputMesh = 0. Wrong number of nodes <> 4 : "<<int(outputMesh->getCoords()->getNumberOfTuples())<<finl;
          status = -4;
          // Process::exit();
          return;
        }
    }

  // Tables de sortie
  const double *mc_finalCoords = outputMesh->getCoords()->begin();
  finalcoords.resize(int(outputMesh->getNumberOfNodes()), outputMesh->getMeshDimension());
  std::copy(mc_finalCoords, mc_finalCoords+finalcoords.size_array(), finalcoords.addr());

  numNodes2D.clear();
  outputMesh->getNodeIdsOfCell(0, numNodes2D);
  for (int node =0; node <  int(numNodes2D.size()); node++) Tab_conect(0,node) = int(numNodes2D[node]);

  return;
}

// intersection d'un segment avec un polygone
void Prepro_IBM_base::intersectSegPoly2D(MEDCouplingUMesh * polyEdgeMesh, DoubleTab& p1, DoubleTab& p2, double eps, MCAuto<MEDCoupling::DataArrayDouble> outputCoords, int& status)
{
  int idebug = verbose_;
  int dim_esp = Objet_U::dimension;
  if (status != 0) status = 0 ;

  int nbEdge = int(polyEdgeMesh->getNumberOfCells());
  int nbNodes = int(polyEdgeMesh->getNumberOfNodes());
  const double *mc_coords = polyEdgeMesh->getCoords()->begin();
  DoubleTrav coords;
  coords.resize(nbNodes, (dim_esp-1));
  std::copy(mc_coords, mc_coords+coords.size_array(), coords.addr());
  std::vector<mcIdType> numNodes;

  MCAuto<MEDCoupling::DataArrayDouble> MC_p(MEDCoupling::DataArrayDouble::New());
  DoubleTrav pe1(dim_esp-1);
  DoubleTrav pe2(dim_esp-1);
  if (idebug) Cerr<<">> Prepro_IBM_base::intersectSegPoly2D: Edge nb "<<nbEdge<<finl;
  for (int i=0; i<nbEdge; i++)
    {
      if (int(polyEdgeMesh->getNumberOfNodesInCell(i)) > 2)
        {
          Cerr<<">> Prepro_IBM_base::intersectSegPoly2D(): >getNumberOfNodesInCell("<<i<<") = "<<int(polyEdgeMesh->getNumberOfNodesInCell(i))<<" CHELOU... Input problem"<<finl;
          status = -1;
          return;
        }
      numNodes.clear();
      polyEdgeMesh->getNodeIdsOfCell(i, numNodes);

      for (int k=0; k<(dim_esp-1); k++) pe1(k) = coords(int(numNodes[0]),k);
      for (int k=0; k<(dim_esp-1); k++) pe2(k) = coords(int(numNodes[1]),k);
      if (idebug)
        {
          Cerr<<">> Prepro_IBM_base::intersectSegPoly2D: Edge "<<i<<" index[] = "<<int(numNodes[0])<<" "<<int(numNodes[1])<<finl;
          Cerr<<">> Prepro_IBM_base::intersectSegPoly2D: pe1 : "<<pe1(0)<<" "<<pe1(1)<<" pe2 : "<<pe2(0)<<" "<<pe2(1)<<finl;
        }
      DoubleTab p(dim_esp - 1);
      int status_segseg = 0;
      intersectSegSeg2D(pe1, pe2, p1, p2, p, eps, status_segseg);
      if (status_segseg == 0)
        {
          MC_p->useArray(p.addr(), false, MEDCoupling::DeallocType::CPP_DEALLOC, 1, (dim_esp-1));
          MC_p->rearrange((dim_esp-1));
          outputCoords->aggregate(MC_p);
        }
    }
  MC_p->decrRef();

  if(outputCoords->getNumberOfTuples () <= 0)
    {
      status = -2;
      return;
    }
}

void Prepro_IBM_base::intersectSegSeg2D(DoubleTab& p11, DoubleTab& p12, DoubleTab& p21, DoubleTab& p22, DoubleTab& p, double eps, int& status)
{
  int idebug = verbose_;
  int dim_esp = Objet_U::dimension;
  if (status != 0) status = 0 ;
  p = 0.;

  assert(p11.size() == (dim_esp-1)); // p11 : sommet 1 du segment 1
  assert(p11.size() == p12.size());
  assert(p12.size() == p21.size()); // p12 : sommet 2 du segment 1
  assert(p21.size() == p22.size()); // p21 : sommet 1 du segment 2
  assert(p22.size() == p.size()); // p22 : sommet 2 du segment 2

  DoubleTrav n1(dim_esp-1); //vecteur directeur du premier segment
  n1(0) = p12(0)-p11(0);
  n1(1) = p12(1)-p11(1);
  DoubleTrav n2(dim_esp-1); //vecteur directeur du second segment
  n2(0) = p22(0)-p21(0);
  n2(1) = p22(1)-p21(1);
  double absprodvect = abs(n1(0)*n2(1) - n2(0)*n1(1)); // n1 ^ n2

  if(absprodvect < eps) // n1 et n2 colineaires
    {
      status = -2;
      if (idebug) Cerr<<"> Prepro_IBM_base::intersectSegSeg2D status: "<<status<<" ; p: "<<p(0)<<" "<<p(1)<<" n1 et n2 colineaires"<< finl;
      return;
    }
  else
    {
      double a = n1(1); // representation en equation de la droite 1
      double b = -n1(0);
      double c = - a*p11(0) - b*p11(1);
      double t1 = - (c + a*p21(0) + b*p21(1)) / (a*n2(0)+b*n2(1)); // representation parametrique de la droite 2
      a = n2(1); // representation en equation de la droite 2
      b = -n2(0);
      c = - a*p21(0) - b*p21(1);
      double t2 = - (c + a*p11(0) + b*p11(1)) / (a*n1(0)+b*n1(1)); // representation parametrique de la droite 1
      if ((t1 < -eps) || (t1 > 1.0 + eps) || (t2 < -eps) || (t2 > 1.0 + eps)) //si le point est en dehors des deux segments
        {
          status = -1;
          return;
        }
      else
        {
          p(0) = p21(0) + n2(0)*t1;
          p(1) = p21(1) + n2(1)*t1;
        }

      if (idebug) Cerr<<"> Prepro_IBM_base::intersectSegSeg2D status: "<<status<<" ; p: "<<p(0)<<" "<<p(1)<<finl;
    }
}

void Prepro_IBM_base::Save_Med_File()
{
  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();
  const Domaine& le_dom = mon_pb_.valeur().domaine();
  int nb_som = le_dom.nb_som();
  const Nom& type_elem = le_dom.type_elem()->que_suis_je();
  int dim_esp = Objet_U::dimension;
  assert(dim_esp==3);
  Noms nom_rot(dim_esp*dim_esp);
  Noms unites_rot(dim_esp*dim_esp);
  Noms nom_vect(dim_esp );
  Noms unites_vect(dim_esp );
  unites_vect[0] = "m";
  unites_vect[1] = "m";
  unites_vect[2] = "m";

  Ecrire_MED ecr_med(nom_fichier_med_Out_, le_dom);
  ecr_med.ecrire_domaine(false);

  ecr_med.ecrire_champ("CHAMPMAILLE", "aire", champ_aire_->valeurs(), champ_aire_->unites(), champ_aire_->noms_compo(), type_elem, 0.);
  ecr_med.ecrire_champ("CHAMPMAILLE", "rotation", champ_rotation_->valeurs(), unites_rot, nom_rot, type_elem, 0.);
  nom_vect[0] = "X";
  nom_vect[1] = "Y";
  nom_vect[2] = "Z";
  ecr_med.ecrire_champ("CHAMPMAILLE", "barycentre", champ_bary_->valeurs(), unites_vect, nom_vect, type_elem, 0.);
  nom_vect[0] = "nx";
  nom_vect[1] = "ny";
  nom_vect[2] = "nz";
  ecr_med.ecrire_champ("CHAMPMAILLE", "normale", champ_normal_->valeurs(), unites_vect, nom_vect, type_elem, 0.);
  ecr_med.ecrire_champ("CHAMPMAILLE", "corresp_elems", corresp_elems_->valeurs(), corresp_elems_->unites(), corresp_elems_->noms_compo(), type_elem, 0.);

  ecr_med.ecrire_champ("CHAMPPOINT", "is_node_dirichlet", isNodeDirichlet_->valeurs(), isNodeDirichlet_->unites(), isNodeDirichlet_->noms_compo(), type_elem, 0.);

  nom_vect[0] = "X";
  nom_vect[1] = "Y";
  nom_vect[2] = "Z";
  ecr_med.ecrire_champ("CHAMPPOINT", "projection_solide", solid_points_->valeurs(), unites_vect, nom_vect, type_elem, 0.);
  ecr_med.ecrire_champ("CHAMPPOINT", "element_projection_solide", solid_elems_->valeurs(), unites_vect, nom_vect, type_elem, 0.);

  ecr_med.ecrire_champ("CHAMPPOINT", "projection_fluide", fluid_points_->valeurs(), unites_vect, nom_vect, type_elem, 0.);
  ecr_med.ecrire_champ("CHAMPPOINT", "element_projection_fluide", fluid_elems_->valeurs(), fluid_elems_->unites(), fluid_elems_->noms_compo(), type_elem, 0.);

  ecr_med.ecrire_champ("CHAMPMAILLE", "h_max_elem", h_max_elem_->valeurs(), h_max_elem_->unites(), h_max_elem_->noms_compo(), type_elem, 0.);
  ecr_med.ecrire_champ("CHAMPPOINT", "h_max_node", h_max_node_->valeurs(), h_max_node_->unites(), h_max_node_->noms_compo(), type_elem, 0.);

  OWN_PTR(Champ_Don_base) sommets_voisins;
  mon_pb_.valeur().discretisation().discretiser_champ("champ_sommets",le_dom_dis,"sommets_voisins","",1,0., sommets_voisins);
  DoubleTab& sommets_voisinsArray = sommets_voisins->valeurs();
  sommets_voisinsArray = 0.0;
  for (int node=0; node<nb_som; node++)
    {
      IntList& the_list = sommets_voisins_[node];
      if (!(the_list.est_vide()))
        {
          for (int index=0; index < the_list.size(); index++) sommets_voisinsArray(the_list[index]) += 1.;
        }
    }
  ecr_med.ecrire_champ("CHAMPPOINT", "sommets_voisins", sommets_voisins->valeurs(), sommets_voisins->unites(), sommets_voisins->noms_compo(), type_elem, 0.);
}

void Prepro_IBM_base::compute_effective_error()
{
  assert(Objet_U::dimension == 3);

  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();
  const Domaine& le_dom= le_dom_dis.domaine();
  const DoubleTab coordsDom3D=le_dom.les_sommets();
  int nb_som_elem =le_dom.nb_som_elem();
  const IntTab& elems = le_dom.les_elems() ;
  int nbElemVol = le_dom.nb_elem();

  // computing effective error (eps_ multiply by min geometric scale)
  double hmin= 1.e30;
  for (int my_elem = 0; my_elem < nbElemVol; my_elem++)
    {
      double xmin = 1.e30, ymin = 1.e30, zmin = 1.e30;
      double xmax = -1.e30, ymax = -1.e30, zmax = -1.e30;
      for (int id_loc=0; id_loc < nb_som_elem; id_loc++)
        {
          int id = elems(my_elem,id_loc);
          xmin = std::min(xmin, coordsDom3D(id , 0));
          ymin = std::min(ymin, coordsDom3D(id , 1));
          zmin = std::min(zmin, coordsDom3D(id , 2));

          xmax = std::max(xmax, coordsDom3D(id , 0));
          ymax = std::max(ymax, coordsDom3D(id , 1));
          zmax = std::max(zmax, coordsDom3D(id , 2));
        }
      double hinter = std::min({xmax - xmin, ymax - ymin, zmax - zmin});
      // dimension geometrique min local elementaire
      hmin = std::min(hmin, hinter);
    }
  eps_effec_ = std::max(hmin*eps_, 1e-16);
  Cerr<<"Prepro_IBM_base:: Eps and effective eps = "<<eps_<<" "<<eps_effec_<<finl;
}

void Prepro_IBM_base::compute_h_max_elem()
{
  DoubleTab& hmaxArray = h_max_elem_->valeurs();
  hmaxArray = 0.;
  DoubleTab& hmaxNodeArray = h_max_node_->valeurs();
  hmaxNodeArray = 0.;

  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();// contient maillage TRUST
  const Domaine& le_dom= le_dom_dis.domaine();
  int nb_som_elem =le_dom.nb_som_elem(); // number of nodes per elems ( numerotation local)
  const IntTab& elems = le_dom.les_elems() ; // numerotation global des noeuds d'un element ( les_elems.dimension(0)=>nb_elems, les_elems.dimension(1)=> nb_som_elem)
  const DoubleTab coordsDom3D=le_dom.les_sommets(); // coordonnées des noeuds TRUST
  int nbElemVol = elems.dimension(0);

  assert(Objet_U::dimension == 3);
  assert(Objet_U::dimension == dimTab_.dimension(0));
  double hmax1, hmax2, hmax3;

  for (int my_elem = 0; my_elem < nbElemVol; my_elem++)
    {
      if (dimTab_(0))
        {
          double xmin = 1.e30;
          double xmax = -1.e30;
          for (int id_loc=0; id_loc < nb_som_elem; id_loc++)
            {
              int id = elems(my_elem,id_loc); // global number of node
              xmin = std::min(xmin, coordsDom3D(id , 0));
              xmax = std::max(xmax, coordsDom3D(id , 0));
            }
          hmax1 = xmax-xmin;
        }
      else hmax1 = 0.0;

      if (dimTab_(1))
        {
          double ymin = 1.e30;
          double ymax = -1.e30;
          for (int id_loc=0; id_loc < nb_som_elem; id_loc++)
            {
              int id = elems(my_elem,id_loc); // global number of node
              ymin = std::min(ymin, coordsDom3D(id , 1));
              ymax = std::max(ymax, coordsDom3D(id , 1));
            }
          hmax2 = ymax-ymin;
        }
      else hmax2 = 0.0;

      if (dimTab_(2))
        {
          double zmin = 1.e30;
          double zmax = -1.e30;
          for (int id_loc=0; id_loc < nb_som_elem; id_loc++)
            {
              int id = elems(my_elem,id_loc); // global number of node
              zmin = std::min(zmin, coordsDom3D(id , 2));
              zmax = std::max(zmax, coordsDom3D(id , 2));
            }
          hmax3 = zmax-zmin;
        }
      else hmax3 = 0.0;

      hmaxArray(my_elem) = sqrt(hmax1*hmax1+hmax2*hmax2+hmax3*hmax3) ;

      for (int id_loc=0; id_loc < nb_som_elem; id_loc++)
        {
          int id = elems(my_elem,id_loc);
          hmaxNodeArray(id) = max(hmaxNodeArray(id), hmaxArray(my_elem)) ;
        }
    }
  hmaxArray.echange_espace_virtuel();
  hmaxNodeArray.echange_espace_virtuel();
}

void Prepro_IBM_base::compute_NeighNode(int nb_niveau)
{
  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();
  const Domaine_VF& the_dom_VF = ref_cast(Domaine_VF,le_dom_dis);
  const Domaine& le_dom =  le_dom_dis.domaine();
  int nb_elem_tot = le_dom.nb_elem_tot();
  int nb_som_tot = le_dom.nb_som_tot();
  int nb_som_elem = le_dom.nb_som_elem();
  int nb_faces_elem = le_dom.nb_faces_elem();
  int nb_som_face = the_dom_VF.nb_som_face();
  const IntTab& face_voisins = the_dom_VF.face_voisins();
  const IntTab& elems = le_dom.les_elems() ;
  const IntTab& elem_face = the_dom_VF.elem_faces();
  const IntTab& faces_som = the_dom_VF.face_sommets();

  DoubleTab& aireArray = champ_aire_->valeurs();

  if (nb_niveau != 1)
    {
      Cerr<<"Prepro_IBM_base::compute_NeighNode :  nb_niveau != 1"<<endl;
      Process::exit();
    }

  sommets_voisins_ = IntLists(nb_som_tot);
  int nb_nodes_in_list = 0;
  for (int num_elem = 0; num_elem < nb_elem_tot; num_elem++) // inclus les elem ghost
    {
      if (aireArray(num_elem) > eps_effec_)
        {
          for (int fac=0; fac<nb_faces_elem; fac++)
            {
              int num_fac = elem_face(num_elem, fac);
              for (int voisin=0; voisin<2; voisin++)
                {
                  int num_elem_v = face_voisins(num_fac,voisin);
                  if ( (num_elem_v!=-1) && (num_elem_v!=num_elem) )
                    {
                      // element num_elem_v, voisin de num_elem par une face
                      for (int nodf=0; nodf<nb_som_face; nodf++)
                        {
                          // noeud num_nod, dont on veut trouve les noeuds voisins
                          int num_nod = faces_som(num_fac, nodf);
                          // tous les noeuds de num_elem_v sont des neouds voisins de num_nod
                          for (int nod_v=0; nod_v<nb_som_elem; nod_v++)
                            {
                              int num_nod_v = elems(num_elem_v, nod_v);
                              if (num_nod_v != num_nod)
                                {
                                  sommets_voisins_[num_nod].add_if_not(num_nod_v);
                                  nb_nodes_in_list += 1;
                                }
                            }
                          // element num_elem_v_v, voisin de num_elem_v par une face
                          for (int fac_v=0; fac_v<nb_faces_elem; fac_v++)
                            {
                              int num_fac_v = elem_face(num_elem_v, fac_v);
                              for (int voisin_v=0; voisin_v<2; voisin_v++)
                                {
                                  int num_elem_v_v = face_voisins(num_fac_v,voisin_v);
                                  // boucle sur les noeuds de l'element num_elem_v_v pour savoir
                                  // s'il contient num_nod. si oui, tous ses noeuds sont des
                                  // voisins de num_nod
                                  if ( (num_elem_v_v!=-1) && (num_elem_v_v!=num_elem) )
                                    {
                                      int iok = 0 ;
                                      for (int nod_v_v=0; nod_v_v<nb_som_elem; nod_v_v++)
                                        {
                                          int num_nod_v_v = elems(num_elem_v_v, nod_v_v);
                                          if (num_nod_v_v == num_nod) iok = 1;
                                        }
                                      if (iok == 1)
                                        {
                                          for (int nod_v_v=0; nod_v_v<nb_som_elem; nod_v_v++)
                                            {
                                              int num_nod_v_v = elems(num_elem_v_v, nod_v_v);
                                              if (num_nod_v_v != num_nod)
                                                {
                                                  sommets_voisins_[num_nod].add_if_not(num_nod_v_v);
                                                  nb_nodes_in_list += 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  Cerr<<"Prepro_IBM_base::compute_NeighNode : nb nodes for neighborhood = " <<nb_nodes_in_list<<finl;
}

void Prepro_IBM_base::calculer_normal_proj_solid(DoubleTab& nor, DoubleTab& dist)
{
  const DoubleTab& solid_points = solid_points_->valeurs();
  const DoubleTab& aire = champ_aire_->valeurs();
  int nb_elem = aire.dimension(0);
  const Domaine_dis_base& le_dom_dis = mon_pb_->domaine_dis();
  const IntTab& elems = le_dom_dis.domaine().les_elems();
  const Domaine& dom = le_dom_dis.domaine();
  int nb_som_elem=dom.nb_som_elem();

  assert(nor.dimension(0) == solid_points.dimension(0));
  int dim_esp =  solid_points.dimension(1);
  assert(nor.dimension(1) == dim_esp);
  nor = 0.;
  double eps = 1.0e-10;
  for (int e=0; e<nb_elem; e++)
    {
      if (aire(e)>0.)
        {
          for (int il=0; il<nb_som_elem; il++)
            {
              int i = elems(e,il);
              double d1 = 0.0;
              for (int k=0; k<dim_esp; k++)
                {
                  double xk = dom.coord(i,k);
                  double xks = solid_points(i,k);
                  nor(i, k) = xk - xks;
                  d1 += (xk-xks)*(xk-xks);
                }
              d1 = sqrt(d1);
              if (d1  > eps )
                for (int k=0; k<dim_esp; k++) nor(i, k) /= d1 ;
              dist(i) = d1;
            }
        }
    }
}

void Prepro_IBM_base::compute_solid_fluid(int maj_from_ext)
{
}
