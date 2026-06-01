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

#include <Linear_algebra_tools_impl.h>
#include <TRUSTList.h>
#include <Polyedre.h>
#include <Domaine.h>
#include <algorithm>

using std::swap;

Implemente_instanciable_sans_constructeur_32_64(Polyedre_32_64,"Polyedre",Poly_geom_base_32_64<_T_>);

template <typename _SIZE_>
Polyedre_32_64<_SIZE_>::Polyedre_32_64(): PolyhedronIndex_(1)
{
  PolyhedronIndex_[0]=0;
  FacesIndex_.resize(1);
  FacesIndex_[0]=0;
  nb_som_elem_max_=-1;
  nb_face_elem_max_=0;
  nb_som_face_max_=0;
}

template <typename _SIZE_>
Sortie& Polyedre_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  s<< Nodes_        <<finl;
  s<< FacesIndex_     <<finl;
  s<< PolyhedronIndex_ <<finl;
  s<< nb_som_elem_max_ <<finl;
  s<< nb_face_elem_max_ <<finl;
  s<< nb_som_face_max_ <<finl;;
  WARN;
  return s;
}

template <typename _SIZE_>
Entree& Polyedre_32_64<_SIZE_>::readOn(Entree& s )
{
  s>> Nodes_;
  s>>FacesIndex_;
  s>>PolyhedronIndex_;
  s>>nb_som_elem_max_;
  s>>nb_face_elem_max_;
  s>>nb_som_face_max_;;
  return s;
}

template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::calculer_centres_gravite(DoubleTab_t& xp) const
{
  const Domaine_t& domaine=mon_dom.valeur();
  const IntTab_t& elem=domaine.les_elems();
  const DoubleTab_t& coord=domaine.coord_sommets();

  int_t nb_elem;
  if(xp.dimension(0)==0)
    {
      nb_elem = domaine.nb_elem_tot();
      xp.resize(nb_elem,dimension);
    }
  else
    nb_elem=xp.dimension(0);


  for (int_t num_poly=0; num_poly<nb_elem; num_poly++)
    {
      double volume=0;
      Vecteur3 xg(0,0,0);
      int nb_som_max=elem.dimension_int(1);
      int nb_som_reel;
      for (nb_som_reel=0; nb_som_reel<nb_som_max; nb_som_reel++)
        {
          int_t n=elem(num_poly,nb_som_reel);
          if (n<0)
            break;
          Vecteur3 S(coord(n,0),coord(n,1),coord(n,2));
          xg+=S;
        }
      xg*=1./nb_som_reel;
      Vecteur3 moinsS0(xg);
      moinsS0*=-1;

      Vecteur3 vraixg(0,0,0);

      for (int_t f=PolyhedronIndex_[num_poly]; f<PolyhedronIndex_[num_poly+1]; f++)
        {
          int somm_loc3=Nodes_[FacesIndex_[f+1]-1];
          int_t n3=elem(num_poly,somm_loc3);
          Vecteur3 S3(coord(n3,0),coord(n3,1),coord(n3,2));
          Vecteur3 S3sa(S3);
          S3 += moinsS0;
          for (int_t s=FacesIndex_[f]; s<FacesIndex_[f+1]-2; s++)
            {
              int somm_loc1=Nodes_[s];
              int_t n1=elem(num_poly,somm_loc1);
              int somm_loc2=Nodes_[s+1];
              int_t n2=elem(num_poly,somm_loc2);
              Vecteur3 S1(coord(n1,0),coord(n1,1),coord(n1,2));
              Vecteur3 S2(coord(n2,0),coord(n2,1),coord(n2,2));

              Vecteur3 xgl(xg);
              xgl+=S3sa;
              xgl+=S1;
              xgl+=S2;
              xgl*=0.25;
              S1+=  moinsS0;
              S2+=  moinsS0;
              double vol_l= std::fabs(
                              S1[0] * ( S2[1] * S3[2] - S3[1] * S2[2] )
                              + S2[0] * ( S3[1] * S1[2] - S1[1] * S3[2] )
                              + S3[0] * ( S1[1] * S2[2] - S2[1] * S1[2] ) );
              volume+=vol_l;
              xgl*=vol_l;
              vraixg+=xgl;
            }
        }
      vraixg*=1./volume;
      xp(num_poly,0)=vraixg[0];
      xp(num_poly,1)=vraixg[1];
      xp(num_poly,2)=vraixg[2];
    }
}

template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::calculer_un_centre_gravite(const int_t num_elem,DoubleVect& xp) const
{
  const IntTab_t& les_Polys = mon_dom->les_elems();
  const Domaine_t& le_domaine = mon_dom.valeur();

  xp.resize(dimension);
  int nb_som_reel=nb_som();
  while (les_Polys(num_elem,nb_som_reel-1)==-1)  nb_som_reel--;
  for(int s=0; s<nb_som_reel; s++)
    {
      int_t num_som = les_Polys(num_elem,s);
      for(int i=0; i<dimension; i++)
        xp(i) += le_domaine.coord(num_som,i)/nb_som_reel;
    }

  WARN;
}

/*! @brief Returns the LML name of a polyhedron = "POLYEDRE_" + nb_som_max.
 *
 * @return LML name string for this polyhedron type.
 */
template <typename _SIZE_>
const Nom& Polyedre_32_64<_SIZE_>::nom_lml() const
{
  static Nom nom;
  nom="POLYEDRE_";
  Nom n(this->get_nb_som_elem_max());
  nom+=n;
  return nom;
}


/*! @brief Returns 1 if element "num_poly" of the associated domain contains the point with coordinates given by "pos". Returns 0 otherwise.
 *
 * @param pos Coordinates of the point to locate.
 * @param num_poly Index of the domain element in which to search for the point.
 * @return 1 if the point belongs to element "num_poly", 0 otherwise.
 */
template <typename _SIZE_>
int Polyedre_32_64<_SIZE_>::contient(const ArrOfDouble& pos, int_t num_poly ) const
{
  // check if point P is on the same side as xg for each face.
  const Domaine_t& domaine = mon_dom.valeur();
  const IntTab_t& elem=domaine.les_elems();
  const DoubleTab_t& coord=domaine.coord_sommets();
  Vecteur3 P(pos[0],pos[1],pos[2]);
  Vecteur3 xg(0,0,0);
  int nb_som_max= elem.dimension_int(1); // yes, cast, this is a higher dimension.
  int nb_som_reel;
  for (nb_som_reel=0; nb_som_reel<nb_som_max; nb_som_reel++)
    {
      int_t n=elem(num_poly,nb_som_reel);
      if (n<0)
        break;
      Vecteur3 S(coord(n,0),coord(n,1),coord(n,2));
      xg+=S;
    }
  xg*=1./nb_som_reel;

  for (int_t f=PolyhedronIndex_[num_poly]; f<PolyhedronIndex_[num_poly+1]; f++)
    {
      Vecteur3 n(0,0,0);
      int somm_loc3=Nodes_[FacesIndex_[f+1]-1];
      int_t n3=elem(num_poly,somm_loc3);
      Vecteur3 moinsS3(-coord(n3,0),-coord(n3,1),-coord(n3,2));

      for (int_t s=FacesIndex_[f]; s<FacesIndex_[f+1]-2; s++)
        {
          int somm_loc1=Nodes_[s];
          int_t n1=elem(num_poly,somm_loc1);
          int somm_loc2=Nodes_[s+1];
          int_t n2=elem(num_poly,somm_loc2);
          Vecteur3 S1(coord(n1,0),coord(n1,1),coord(n1,2));
          Vecteur3 S2(coord(n2,0),coord(n2,1),coord(n2,2));
          S1+=moinsS3;
          S2+=moinsS3;
          Vecteur3 nTr;
          Vecteur3::produit_vectoriel(S1,S2,nTr);
          n+=nTr;
        }

      Vecteur3 S3G(xg);
      S3G+=moinsS3;
      Vecteur3 S3P(P);
      S3P+=moinsS3;
      double prod_scal1=Vecteur3::produit_scalaire(n,S3G);
      double prod_scal2=Vecteur3::produit_scalaire(n,S3P);
      if (prod_scal1*prod_scal2<-Objet_U::precision_geom*prod_scal1*prod_scal1)
        return 0;
    }

  return 1;
}


/*! @brief Not yet implemented — always returns 0. Returns 1 if the vertices specified by "pos" are those of element "num_poly" in the associated domain.
 *
 * @param pos Vertex indices to compare.
 * @param num_poly Index of the domain element whose vertices are to be compared.
 * @return 1 if the vertices match, 0 otherwise.
 */
template <typename _SIZE_>
int Polyedre_32_64<_SIZE_>::contient(const SmallArrOfTID_t& pos, int_t num_poly ) const
{
  BLOQUE;
  return 0;
}


/*! @brief Computes the volumes of the elements of the associated domain.
 *
 * @param volumes Vector to fill with the volumes of domain elements.
 */
template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::calculer_volumes(DoubleVect_t& volumes) const
{
  const Domaine_t& domaine=mon_dom.valeur();
  const IntTab_t& elem=domaine.les_elems();
  const DoubleTab_t& coord=domaine.coord_sommets();

  int_t size_tot = domaine.nb_elem_tot();
  assert(volumes.size_totale()==size_tot);
  for (int_t num_poly=0; num_poly<size_tot; num_poly++)
    {
      double volume=0;
      Vecteur3 xg(0,0,0);
      int nb_som_max = elem.dimension_int(1);
      int nb_som_reel;
      for (nb_som_reel=0; nb_som_reel<nb_som_max; nb_som_reel++)
        {
          int_t n=elem(num_poly,nb_som_reel);
          if (n<0)
            break;
          Vecteur3 S(coord(n,0),coord(n,1),coord(n,2));
          xg+=S;
        }
      xg*=1./nb_som_reel;
      Vecteur3 moinsS0(xg);
      moinsS0*=-1;
      for (int_t f=PolyhedronIndex_[num_poly]; f<PolyhedronIndex_[num_poly+1]; f++)
        {
          int somm_loc3=Nodes_[FacesIndex_[f+1]-1];
          int_t n3=elem(num_poly,somm_loc3);
          Vecteur3 S3(coord(n3,0),coord(n3,1),coord(n3,2));
          S3+=  moinsS0;
          for (int_t s=FacesIndex_[f]; s<FacesIndex_[f+1]-2; s++)
            {
              int somm_loc1=Nodes_[s];
              int_t n1=elem(num_poly,somm_loc1);
              int somm_loc2=Nodes_[s+1];
              int_t n2=elem(num_poly,somm_loc2);
              Vecteur3 S1(coord(n1,0),coord(n1,1),coord(n1,2));
              Vecteur3 S2(coord(n2,0),coord(n2,1),coord(n2,2));
              S1+=  moinsS0;
              S2+=  moinsS0;
              volume += std::fabs(
                          S1[0] * ( S2[1] * S3[2] - S3[1] * S2[2] )
                          + S2[0] * ( S3[1] * S1[2] - S1[1] * S3[2] )
                          + S3[0] * ( S1[1] * S2[2] - S2[1] * S1[2] ) );
            }
        }
      volumes(num_poly)=volume/6.;
    }
}

/*! @brief Fills the table faces_som_local(i,j), which gives for 0 <= i < nb_faces() and 0 <= j < nb_som_face(i) the local vertex index on the element.
 *
 * We have 0 <= faces_sommets_locaux(i,j) < nb_som().
 * If faces do not all have the same number of vertices, the number of columns is
 * the maximum, and unused entries are set to -1.
 * Returns 1 if all faces have the same number of vertices, 0 otherwise.
 *
 * @param faces_som_local Table to fill with local face-vertex indices.
 * @return 1 if all faces have the same vertex count, 0 otherwise.
 */
template <typename _SIZE_>
int Polyedre_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local) const
{
  return 0;
}

template <typename _SIZE_>
int Polyedre_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local,int_t ele) const
{
  faces_som_local.resize(nb_face_elem_max_,nb_som_face_max_);
  faces_som_local=-1;
  if (ele == 0 && PolyhedronIndex_.size_array() == 1)
    return 1; //no elements!

  // look for the faces of the element
  int fl=0;
  for (int_t f=PolyhedronIndex_[ele]; f<PolyhedronIndex_[ele+1]; f++) // for all faces of the element 'ele'
    {
      int sl=0;
      for (int_t s=FacesIndex_[f]; s<FacesIndex_[f+1]; s++)  // for all (global) indices of vertex constituting the face
        {
          int somm_loc=Nodes_[s];   // its local numbering
          faces_som_local(fl,sl)=somm_loc;
          sl++;
        }
      fl++;
    }

  return 1;
}

// From the indirection arrays FacesIndex, PolyhedronIndex and Nodes,
// compute les_elems, nb_som_face_max_, nb_face_elem_max_, nb_som_elem_max_,
// and local Nodes...
template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::affecte_connectivite_numero_global(const ArrOfInt_t& Nodes, const ArrOfInt_t& FacesIndex,const ArrOfInt_t& PolyhedronIndex,IntTab_t& les_elems)
{
  nb_som_elem_max_=0;
  // determine the maximum number of vertices per element
  TRUSTList<_SIZE_> prov;
  nb_face_elem_max_=0;
  nb_som_face_max_=0;
  int_t nelem=PolyhedronIndex.size_array()-1;
  for (int_t ele=0; ele<nelem; ele++)
    {
      prov.vide();
      int nbf=(int)(PolyhedronIndex[ele+1]-PolyhedronIndex[ele]); // num of faces always int
      if (nbf>nb_face_elem_max_) nb_face_elem_max_=nbf;
      for (int_t f=PolyhedronIndex[ele]; f<PolyhedronIndex[ele+1]; f++)
        {
          //Cerr<<" ici "<<ele << " " <<f <<" "<<FacesIndex(f+1)-FacesIndex(f)<<finl;
          int nbs=(int)(FacesIndex[f+1]-FacesIndex[f]);
          if (nbs>nb_som_face_max_) nb_som_face_max_=nbs;
          for (int_t s=FacesIndex[f]; s<FacesIndex[f+1]; s++)
            prov.add_if_not(Nodes[s]);
        }
      int nbsom=prov.size();
      if (nbsom>nb_som_elem_max_) nb_som_elem_max_=nbsom;
    }
  nb_som_elem_max_ = Process::mp_max(nb_som_elem_max_);
  nb_som_face_max_ = Process::mp_max(nb_som_face_max_);
  nb_face_elem_max_ = Process::mp_max(nb_face_elem_max_);
  Cerr<<" Polyhedron information nb_som_elem_max "<< nb_som_elem_max_<<" nb_som_face_max "<<nb_som_face_max_<<" nb_face_elem_max "<<nb_face_elem_max_<<finl;
  les_elems.resize(nelem,nb_som_elem_max_);
  les_elems=-1;
  // second pass to determine les_elems
  for (int_t ele=0; ele<nelem; ele++)
    {
      prov.vide();
      for (int_t f=PolyhedronIndex[ele]; f<PolyhedronIndex[ele+1]; f++)
        for (int_t s=FacesIndex[f]; s<FacesIndex[f+1]; s++)
          prov.add_if_not(Nodes[s]);
      int nbsom=prov.size();
      // sort prov in ascending order
      // not strictly necessary but ensures a consistent elem table even after face permutations
      bool perm=1;
      while (perm)
        {
          perm=false;
          for (int i=0; i<nbsom-1; i++)
            if (prov[i]>prov[i+1])
              {
                perm=true;
                int_t sa=prov[i];
                prov[i]=prov[i+1];
                prov[i+1]=sa;
              }
        }
      for (int s=0; s<nbsom; s++)
        les_elems(ele,s)=prov[s];
    }
  FacesIndex_=FacesIndex;
  PolyhedronIndex_=PolyhedronIndex;

  // Determine Nodes_...
  Nodes_.resize_array(Nodes.size_array());
  Nodes_=-2;
  for (int_t ele=0; ele<nelem; ele++)
    for (int_t f=PolyhedronIndex[ele]; f<PolyhedronIndex[ele+1]; f++)
      for (int_t s=FacesIndex[f]; s<FacesIndex[f+1]; s++)
        {
          int_t somm_glob=Nodes[s];
          for (int sl=0; sl<nb_som_elem_max_; sl++)
            if (les_elems(ele,sl)==somm_glob)
              {
                Nodes_[s]=sl;
                break;
              }
        }
  assert(min_array(Nodes_)>-1);
}

template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::remplir_Nodes_glob(ArrOfInt_t& Nodes_glob,const IntTab_t& les_elems) const
{
  Nodes_glob.resize_array(Nodes_.size_array());
  int_t nelem=les_elems.dimension_tot(0);
  for (int_t ele=0; ele<nelem; ele++)
    for (int_t f=PolyhedronIndex_[ele]; f<PolyhedronIndex_[ele+1]; f++)
      for (int_t s=FacesIndex_[f]; s<FacesIndex_[f+1]; s++)
        {
          int somm_loc=Nodes_[s];
          Nodes_glob[s]=les_elems(ele,somm_loc);
        }
}

/*! @brief Appends elements of type new_elem to those already present in les_elems and new_elems.
 *
 * @param type_elem Geometric element type for the new elements.
 * @param new_elems Connectivity table of the new elements to append.
 * @param les_elems Connectivity table to update in place.
 */
template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::ajouter_elements(const Elem_geom_base_32_64<_SIZE_>& type_elem, const IntTab_t& new_elems, IntTab_t& les_elems)
{
  // Append new_elems to les_elems
  int_t nb_old_elem=les_elems.dimension(0);
  int_t nb_new_elem=new_elems.dimension(0);
  int nb_som_old_elem=les_elems.dimension_int(1);
  int nb_som_new_elem=new_elems.dimension_int(1);
  nb_som_elem_max_ = std::max(nb_som_old_elem,nb_som_new_elem);

  les_elems.resize(nb_old_elem+nb_new_elem,nb_som_elem_max_, RESIZE_OPTIONS::COPY_NOINIT);
  for (int_t el=0; el<nb_new_elem; el++)
    {
      for (int s=0; s<nb_som_new_elem; s++)
        les_elems(nb_old_elem+el,s)=new_elems(el,s);

      for (int s=nb_som_new_elem; s<nb_som_old_elem; s++)
        les_elems(nb_old_elem+el,s)=-1;
    }
  // Add the faces: retrieve the face creation array
  IntTab faces_som_local;
  type_elem.get_tab_faces_sommets_locaux(faces_som_local);
  int nb_face_new_elem=faces_som_local.dimension(0);
  int nb_som_face_new_elem=faces_som_local.dimension(1);
  nb_face_elem_max_=std::max(nb_face_elem_max_,nb_face_new_elem);
  nb_som_face_max_=std::max(nb_som_face_max_,nb_som_face_new_elem);

  PolyhedronIndex_.resize_array(1+nb_old_elem+nb_new_elem);

  int_t old_faces_index= FacesIndex_.size_array();
  FacesIndex_.resize_array(old_faces_index+nb_new_elem*nb_face_new_elem);

  int_t old_nodes_index= Nodes_.size_array();
  Nodes_.resize_array(old_nodes_index+nb_new_elem*nb_face_new_elem*nb_som_face_new_elem);

  int new_s=0;
  for (int_t el=0; el<nb_new_elem; el++)
    {
      PolyhedronIndex_[nb_old_elem+el+1]=PolyhedronIndex_[nb_old_elem+el]+nb_face_new_elem;
      for (int f=0; f<nb_face_new_elem; f++)
        {
          int nb_som_face_this_elem=0;
          for (int s=0; s<nb_som_face_new_elem; s++)
            if (faces_som_local(f,s)!=-1)
              {
                Nodes_[old_nodes_index+new_s++]=faces_som_local(f,s);
                nb_som_face_this_elem++;
              }
          if (nb_som_face_this_elem==4)
            {
              // swap 3 and 4 in case of quadrangle
              int_t last=old_nodes_index+new_s-1;
              swap(Nodes_[last],Nodes_[last-1]);
            }
          FacesIndex_[old_faces_index+(el*nb_face_new_elem)+f]=
            FacesIndex_[old_faces_index+(el*nb_face_new_elem)+f-1]+nb_som_face_this_elem;
        }
    }
  Nodes_.resize_array(old_nodes_index+new_s);
}

/* Build a reduced version of the polytope connectivity when splitting domains in DomainCutter - this always produce
 * a 32b object:
 */
template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::build_reduced(OWN_PTR(Elem_geom_base_32_64<int>)& type_elem, const ArrOfInt_t& elems_sous_part) const
{
  type_elem.typer("Polyedre");
  Polyedre_32_64<int>& reduced = ref_cast(Polyedre_32_64<int>, type_elem.valeur());
  reduced.nb_som_elem_max_  = nb_som_elem_max_;
  reduced.nb_face_elem_max_ = nb_face_elem_max_;
  reduced.nb_som_face_max_  = nb_som_face_max_;

  ArrOfInt& Pi = reduced.PolyhedronIndex_, &Fi = reduced.FacesIndex_;
  ArrOfInt& N = reduced.Nodes_;

  for (int_t i = 0; i < elems_sous_part.size_array(); i++)
    {
      int_t e = elems_sous_part[i];
      for (int_t f = PolyhedronIndex_[e]; f < PolyhedronIndex_[e + 1]; f++)
        {
          for (int_t s = FacesIndex_[f]; s < FacesIndex_[f + 1]; s++)
            N.append_array(Nodes_[s]);
          Fi.append_array(N.size_array());
        }
      Pi.append_array(Fi.size_array() - 1);
    }
}

template <typename _SIZE_>
void Polyedre_32_64<_SIZE_>::compute_virtual_index()
{
  using BigIntTab_t = TRUSTTab<int, _SIZE_>;   // A big tab but only containing ints.

  // Brute-force approach to begin with...
  BigIntTab_t faces_som(0,nb_face_elem_max_,nb_som_face_max_);
  mon_dom->creer_tableau_elements(faces_som);

  IntTab faces_som_local;
  int_t nb_elem=mon_dom->nb_elem();
  int_t nb_elem_tot=mon_dom->nb_elem_tot();
  for (int_t ele=0; ele<nb_elem; ele++)
    {
      get_tab_faces_sommets_locaux(faces_som_local,ele);
      for (int k=0; k<nb_face_elem_max_; k++)
        for (int l=0; l<nb_som_face_max_; l++)
          faces_som(ele,k,l)=faces_som_local(k,l);
    }
  faces_som.echange_espace_virtuel();
  int_t nbs=0;

  PolyhedronIndex_.resize(nb_elem_tot+1);
  //Cerr<<"uuu "<< PolyhedronIndex_<<finl;

  for (int_t ele=nb_elem; ele<nb_elem_tot; ele++)
    {
      int nbf=0;
      for (int k=0; k<nb_face_elem_max_; k++)
        {
          if (faces_som(ele,k,0)!=-1)
            nbf++;
          for (int l=0; l<nb_som_face_max_; l++)
            {
              if (faces_som(ele,k,l)!=-1)
                nbs++;
              // Cerr <<" INFO " << ele<< " k" <<k << " l "<< l<< " iiii "<<faces_som(ele,k,l)<<finl;
            }
        }
      PolyhedronIndex_[ele+1]=PolyhedronIndex_[ele]+nbf;
    }
  FacesIndex_.resize(PolyhedronIndex_[nb_elem_tot]+1);
  int_t nbs_old=Nodes_.size_array();
  Nodes_.resize(nbs_old+nbs);
  int_t nbft=PolyhedronIndex_[nb_elem];
  nbs=nbs_old;

  for (int_t ele=nb_elem; ele<nb_elem_tot; ele++)
    {
      for (int k=0; k<nb_face_elem_max_; k++)
        {
          for (int l=0; l<nb_som_face_max_; l++)
            if (faces_som(ele,k,l)!=-1)
              {
                Nodes_[nbs]=faces_som(ele,k,l);
                nbs++;
              }
          if (faces_som(ele,k,0)!=-1)
            {
              FacesIndex_[nbft+1]=nbs;
              nbft++;
            }
        }
    }

  for (int_t ele=nb_elem; ele<nb_elem_tot; ele++)
    {
      get_tab_faces_sommets_locaux(faces_som_local,ele);
      for (int k=0; k<nb_face_elem_max_; k++)
        for (int l=0; l<nb_som_face_max_; l++)
          {
            int_t ind1=faces_som(ele,k,l),
                  ind2=faces_som_local(k,l);
            if (ind1!=ind2)
              {
                Cerr << "PPPPB "<< ele<< " k " <<k << " l "<< l<< " iiii "<<ind1<<" "<<ind2<<finl;
                abort();
              }
          }
    }
}

template <typename _SIZE_>
typename Polyedre_32_64<_SIZE_>::int_t Polyedre_32_64<_SIZE_>::get_somme_nb_faces_elem() const
{
  int_t titi= PolyhedronIndex_[mon_dom->nb_elem()];
  return titi;
}


template class Polyedre_32_64<int>;
#if INT_is_64_ == 2
template class Polyedre_32_64<trustIdType>;
#endif
