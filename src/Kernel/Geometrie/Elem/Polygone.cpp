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

#include <TRUSTList.h>
#include <Polygone.h>
#include <Triangle.h>
#include <Domaine.h>
#include <Polygon_geom_tools.h>
#include <algorithm>

Implemente_instanciable_sans_constructeur_32_64(Polygone_32_64,"Polygone",Poly_geom_base_32_64<_T_>);

template <typename _SIZE_>
Polygone_32_64<_SIZE_>::Polygone_32_64(): PolygonIndex_(1)
{
  PolygonIndex_[0]=0;
  FacesIndex_.resize(1);
  FacesIndex_[0]=0;
  nb_som_elem_max_=-1;
  nb_face_elem_max_=0;
}

template <typename _SIZE_>
Sortie& Polygone_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  s<< FacesIndex_     <<finl;
  s<< PolygonIndex_ <<finl;
  s<< nb_som_elem_max_ <<finl;
  s<< nb_face_elem_max_ <<finl;
  WARN;
  return s;
}

template <typename _SIZE_>
Entree& Polygone_32_64<_SIZE_>::readOn(Entree& s )
{
  s>>FacesIndex_;
  s>>PolygonIndex_;
  s>>nb_som_elem_max_;
  s>>nb_face_elem_max_;
  return s;
}


template <typename _SIZE_>
void Polygone_32_64<_SIZE_>::rebuild_index()
{
  const IntTab_t& les_elems = mon_dom->les_elems();
  int_t nb_elem = les_elems.dimension_tot(0);
  ArrOfInt_t PolygonIndex_OK(nb_elem+1);
  PolygonIndex_OK[0]=0;
  for (int_t ele=0; ele<nb_elem; ele++)
    {
      int nbf=get_nb_som_elem_max();
      while (les_elems(ele,nbf-1)<0)
        nbf--;
      PolygonIndex_OK[ele+1]= PolygonIndex_OK[ele]+nbf;
    }
  ArrOfInt_t FacesIndex_OK(PolygonIndex_OK[nb_elem]);
  int_t f=0;
  for (int_t ele=0; ele<nb_elem; ele++)
    for (int ss=0; ss<(int)(PolygonIndex_OK[ele+1]-PolygonIndex_OK[ele]); ss++)  // yes, difference of long giving an int, hence the cast -> face size
      FacesIndex_OK[f++]= les_elems(ele,ss);

  assert(f==PolygonIndex_OK[nb_elem]);

  FacesIndex_=FacesIndex_OK;
  PolygonIndex_=PolygonIndex_OK;
}

/* Build a reduced version of the polytope connectivity when splitting domains in DomainCutter - this always produce
 * a 32b object:
 */
template <typename _SIZE_>
void Polygone_32_64<_SIZE_>::build_reduced(OWN_PTR(Elem_geom_base_32_64<int>)& type_elem, const ArrOfInt_t& elems_sous_part) const
{
  type_elem.typer("Polygone");
  Polygone_32_64<int>& reduced = ref_cast(Polygone_32_64<int>, type_elem.valeur());
  reduced.nb_som_elem_max_  = nb_som_elem_max_;
  reduced.nb_face_elem_max_ = nb_face_elem_max_;

  const IntTab_t& les_elems = mon_dom->les_elems();
  ArrOfInt& Pi = reduced.PolygonIndex_, &Fi = reduced.FacesIndex_;
  Fi.resize(0);

  for (int_t i = 0; i < elems_sous_part.size_array(); i++)
    {
      int_t e = elems_sous_part[i];
      for (int_t f = PolygonIndex_[e]; f < PolygonIndex_[e + 1]; f++)
        {
          int nf = static_cast<int>(f - PolygonIndex_[e]);  // num of faces always an int
          // The below is not necessary (contrarly to what's done for polyedrons) since get_tab_faces_sommets_locaux() below only uses Pi
//          Fi.append_array(les_elems(e, nf));
          Fi.append_array(les_elems(e, nf) > 0 ? 1 : -1);
        }
      Pi.append_array(Fi.size_array());  // this is what will be used by get_tab_faces_sommets_locaux()
    }
}


template <typename _SIZE_>
void Polygone_32_64<_SIZE_>::compute_virtual_index()
{
  rebuild_index();
}


template <typename _SIZE_>
_SIZE_ Polygone_32_64<_SIZE_>::get_somme_nb_faces_elem() const
{
  return PolygonIndex_[mon_dom->nb_elem()];
}

template <typename _SIZE_>
int Polygone_32_64<_SIZE_>::get_nb_som_elem_max() const
{
  if (nb_som_elem_max_>-1)
    return nb_som_elem_max_ ;
  else
    return mon_dom->les_elems().dimension_int(1);
}

/*! @brief Returns the LML name of a polygon = "POLYEDRE_" + 2*nb_som_max (or "POLYGONE_" + nb_som_max in 3D).
 *
 * @return LML name string for this polygon type.
 */
template <typename _SIZE_>
const Nom& Polygone_32_64<_SIZE_>::nom_lml() const
{
  static Nom nom;
  nom="POLYEDRE_";
  Nom n(2*get_nb_som_elem_max());
  if (dimension==3) nom="POLYGONE_";
  if (dimension==3) n=Nom(get_nb_som_elem_max());
  nom+=n;
  return nom;
}


// ToDo move to Triangle
template <typename _SIZE_>
int contient_triangle(const ArrOfDouble& pos, _SIZE_ som0, _SIZE_ som1, _SIZE_ som2, const TRUSTTab<double, _SIZE_>& coord)
{
  double prod,p0,p1,p2;

  // Determine the orientation (counter-clockwise or clockwise) for the vertex numbering:
  // Compute prod = 01 cross 02 along z
  // prod > 0 : counter-clockwise
  // prod < 0 : clockwise
  prod = (coord(som1,0)-coord(som0,0))*(coord(som2,1)-coord(som0,1))
         - (coord(som1,1)-coord(som0,1))*(coord(som2,0)-coord(som0,0));
  double signe;
  if (prod >= 0)
    signe = 1;
  else
    signe = -1;
  // Compute p0 = 0M cross 1M along z
  p0 = (pos[0]-coord(som0,0))*(pos[1]-coord(som1,1))
       - (pos[1]-coord(som0,1))*(pos[0]-coord(som1,0));
  p0 *= signe;
  // Compute p1 = 1M cross 2M along z
  p1 = (pos[0]-coord(som1,0))*(pos[1]-coord(som2,1))
       - (pos[1]-coord(som1,1))*(pos[0]-coord(som2,0));
  p1 *= signe;
  // Compute p2 = 2M cross 0M along z
  p2 = (pos[0]-coord(som2,0))*(pos[1]-coord(som0,1))
       - (pos[1]-coord(som2,1))*(pos[0]-coord(som0,0));
  p2 *= signe;
  double epsilon=std::fabs(prod)*Objet_U::precision_geom;
  if ((p0>-epsilon) && (p1>-epsilon) && (p2>-epsilon))
    return 1;
  else
    return 0;
}

/*! @brief Returns 1 if element "num_poly" of the domain associated with this geometric element contains the point with coordinates "pos_r".
 *
 * Returns 0 otherwise. Implemented by decomposing the polygon into triangles.
 *
 * @param pos_r Coordinates of the point to locate.
 * @param num_poly Index of the domain element in which to search for the point.
 * @return 1 if the point belongs to element "num_poly", 0 otherwise.
 */
template <typename _SIZE_>
int Polygone_32_64<_SIZE_>::contient(const ArrOfDouble& pos_r, int_t num_poly ) const
{
  const Domaine_t& domaine=mon_dom.valeur();
  const IntTab_t& elem=domaine.les_elems();
  const DoubleTab_t& coord=domaine.coord_sommets();
  //DoubleTab pos(3,dimension);
  // decompose the polygon into triangles all sharing vertex 0.

  int_t s0=elem(num_poly,0);
  for (int s=1; s<nb_som_elem_max_-1 ; s++)
    {
      int_t s1=elem(num_poly,s);
      int_t s2=elem(num_poly,s+1);
      if (s2<0)
        break;

      if (contient_triangle(pos_r,s0,s1,s2,coord))
        return 1;
    }

  return 0;
}


/*! @brief Not yet implemented — always returns 0. Returns 1 if the vertices specified by "pos" are those of element "element" in the associated domain.
 *
 * @param pos Vertex indices to compare.
 * @param element Index of the domain element whose vertices are to be compared.
 * @return 1 if the vertices match, 0 otherwise.
 */
template <typename _SIZE_>
int Polygone_32_64<_SIZE_>::contient(const SmallArrOfTID_t& pos, int_t element ) const
{
  BLOQUE;
  return 0;
}


/*! @brief Computes the volumes (areas) of the elements of the associated domain.
 *
 * @param volumes Vector to fill with the volumes of domain elements.
 */
template <typename _SIZE_>
void Polygone_32_64<_SIZE_>::calculer_volumes(DoubleVect_t& volumes) const
{
  const Domaine_t& domaine = mon_dom.valeur();
  const IntTab_t& elem = domaine.les_elems();
  const DoubleTab_t& coord = domaine.coord_sommets();
  int_t size = domaine.nb_elem();

  assert(volumes.size_totale()==domaine.nb_elem_tot());

  for (int_t num_poly = 0; num_poly < size; num_poly++)
    {
      // Determine the actual number of vertices for this polygon (terminated by -1)
      int nbsom = 0;
      const int nbsom_max = get_nb_som_elem_max();
      while (nbsom < nbsom_max && elem(num_poly, nbsom) >= 0) nbsom++;
      if (nbsom < 3)
        {
          volumes(num_poly) = 0.;
          continue;
        }

      const auto index_of = [&](int i) -> int_t { return elem(num_poly, i); };
      const Polygon_geom_data geom = compute_polygon_geom(coord, dimension, nbsom, index_of, Objet_U::bidim_axi);

      if (!Objet_U::bidim_axi)
        volumes(num_poly) = geom.area_;
      else
        volumes(num_poly) = 2.0 * M_PI * std::fabs(geom.moment_r_);
    }

  volumes.echange_espace_virtuel();
  return;
}

/*! @brief Fills faces_som_local(i,j) giving for 0 <= i < nb_faces() and 0 <= j < nb_som_face(i) the local vertex index on the element.
 *
 * We have 0 <= faces_sommets_locaux(i,j) < nb_som().
 * If faces do not all have the same number of vertices, the number of columns
 * equals the maximum, and unused entries are set to -1.
 * Returns 1 if all faces have the same vertex count, 0 otherwise.
 *
 * @param faces_som_local Table to fill with local face-vertex indices.
 * @return 1 if all faces have the same vertex count, 0 otherwise.
 */
template <typename _SIZE_>
int Polygone_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local) const
{
  return 0;
}

template <typename _SIZE_>
int Polygone_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local,int_t ele) const
{
  faces_som_local.resize(nb_face_elem_max_,nb_som_face());
  faces_som_local=-1;

  // look for the faces of the element
  int nb_face = static_cast<int>(PolygonIndex_[ele+1]-PolygonIndex_[ele]); // always within int

  // [ABN] Duh?! always assume consecutive connectivity??
  for (int fl=0; fl<nb_face-1; fl++)
    {
      faces_som_local(fl,0)=fl;
      faces_som_local(fl,1)=fl+1;
    }

  // Last face:
  int fl=nb_face-1;
  faces_som_local(fl,0)=fl;
  faces_som_local(fl,1)=0;

  return 1;
}

// From the indirection arrays FacesIndex and PolygonIndex,
// compute les_elems, nb_som_face_max_, nb_face_elem_max_, nb_som_elem_max_.
template <typename _SIZE_>
void Polygone_32_64<_SIZE_>::affecte_connectivite_numero_global(const ArrOfInt_t& FacesIndex,const ArrOfInt_t& PolygonIndex,IntTab_t& les_elems)
{
  nb_som_elem_max_=0;
  // determine the maximum number of vertices per element
  TRUSTList<_SIZE_> prov;
  nb_face_elem_max_=0;
  int_t nelem=PolygonIndex.size_array()-1;
  for (int_t ele=0; ele<nelem; ele++)
    {
      prov.vide();
      int_t nbf=PolygonIndex[ele+1]-PolygonIndex[ele];
      if (nbf>nb_face_elem_max_) nb_face_elem_max_=(int)nbf;
      for (int_t f=PolygonIndex[ele]; f<PolygonIndex[ele+1]; f++)
        prov.add_if_not(FacesIndex[f]);
      int nbsom=prov.size();
      if (nbsom>nb_som_elem_max_) nb_som_elem_max_=nbsom;
    }
  nb_som_elem_max_ = Process::mp_max(nb_som_elem_max_);
  nb_face_elem_max_ = Process::mp_max(nb_face_elem_max_);
  Cerr<<" Polygon information nb_som_elem_max "<< nb_som_elem_max_<<" nb_face_elem_max "<<nb_face_elem_max_<<finl;
  les_elems.resize(nelem,nb_som_elem_max_);
  les_elems=-1;
  // second pass to determine les_elems
  for (int_t ele=0; ele<nelem; ele++)
    {
      prov.vide();
      for (int_t f=PolygonIndex[ele]; f<PolygonIndex[ele+1]; f++)
        prov.add_if_not(FacesIndex[f]);
      int nbsom=prov.size();
      for (int s=0; s<nbsom; s++)
        les_elems(ele,s)=prov[s];
    }
  FacesIndex_=FacesIndex;
  PolygonIndex_=PolygonIndex;
  assert(nb_face_elem_max_==nb_som_elem_max_);
}


template <typename _SIZE_>
void Polygone_32_64<_SIZE_>::calculer_centres_gravite(DoubleTab_t& xp) const
{
  const Domaine_t& domaine=mon_dom.valeur();
  const IntTab_t& elem=domaine.les_elems();
  const DoubleTab_t& coord=domaine.coord_sommets();
  int_t nb_elem;
  if(xp.dimension(0)==0)
    {
      nb_elem = mon_dom->nb_elem_tot();
      xp.resize(nb_elem,dimension);
    }
  else
    nb_elem=xp.dimension(0);

  xp=0;
  DoubleTab pos(3,dimension);
  ArrOfDouble xpl(dimension);
  for (int_t num_poly=0; num_poly<nb_elem; num_poly++)
    {
      double aire=0;
      xpl=0;
      int_t s0=elem(num_poly,0);
      for (int d=0; d<dimension; d++)
        pos(0,d)=coord(s0,d);
      for (int s=1; s<get_nb_som_elem_max()-1 ; s++)
        {
          int_t s1=elem(num_poly,s);
          int_t s2=elem(num_poly,s+1);
          if (s2<0)
            break;
          for (int d=0; d<dimension; d++)
            {
              pos(1,d)=coord(s1,d);
              pos(2,d)=coord(s2,d);
            }
          double airel = aire_triangle(pos);
          for (int d=0; d<dimension; d++)
            xpl[d]+=airel*(pos(0,d)+pos(1,d)+pos(2,d));
          aire+=airel;
        }
      aire*=3.;
      for (int d=0; d<dimension; d++)
        xp(num_poly,d)=xpl[d]/(aire);
    }
}

template <typename _SIZE_>
void Polygone_32_64<_SIZE_>::calculer_un_centre_gravite(const int_t num_poly,DoubleVect& xp) const
{
  const Domaine_t& domaine=mon_dom.valeur();
  const IntTab_t& elem=domaine.les_elems();
  const DoubleTab_t& coord=domaine.coord_sommets();
  xp.resize(dimension);

  xp=0;
  DoubleTab pos(3,dimension);
  ArrOfDouble xpl(dimension);
  {
    double aire=0;
    xpl=0;
    int_t s0=elem(num_poly,0);
    for (int d=0; d<dimension; d++)
      pos(0,d)=coord(s0,d);
    for (int s=1; s<nb_som_elem_max_-1 ; s++)
      {
        int_t s1=elem(num_poly,s);
        int_t s2=elem(num_poly,s+1);
        if (s2<0)
          break;
        for (int d=0; d<dimension; d++)
          {
            pos(1,d)=coord(s1,d);
            pos(2,d)=coord(s2,d);
          }
        double airel = aire_triangle(pos);
        for (int d=0; d<dimension; d++)
          xpl[d]+=airel*(pos(0,d)+pos(1,d)+pos(2,d));
        aire+=airel;
      }
    aire*=3.;
    for (int d=0; d<dimension; d++)
      xp(d)=xpl[d]/(aire);
  }
}



template class Polygone_32_64<int>;
#if INT_is_64_ == 2
template class Polygone_32_64<trustIdType>;
#endif
