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

#include <Point.h>
#include <Domaine.h>

Implemente_instanciable_32_64(Point_32_64,"Point",Elem_geom_base_32_64<_T_>);
// XD point points point NO_BRACE Point as class-daughter of Points.

/*! @brief Does nothing.
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the output stream
 */
template<typename _SIZE_>
Sortie& Point_32_64<_SIZE_>::printOn(Sortie& s ) const
{
  return s;
}


/*! @brief Does nothing.
 *
 * @param (Entree& s) an input stream
 * @return (Entree&) the input stream
 */
template<typename _SIZE_>
Entree& Point_32_64<_SIZE_>::readOn(Entree& s )
{
  return s;
}

/*! @brief Returns the LML name of a point element = "VOXEL8".
 *
 * @return (Nom&) always equal to "VOXEL8"
 */
template<typename _SIZE_>
const Nom& Point_32_64<_SIZE_>::nom_lml() const
{
  static Nom nom="POINT";
  return nom;
}


/*! @brief Returns 1 if element ielem of the domain associated with the geometric element contains the point
 *
 *               with coordinates specified by parameter "pos".
 *     Returns 0 otherwise.
 *
 * @param (DoubleVect& pos) coordinates of the point to locate
 * @param (int element) the index of the domain element in which the point is searched.
 * @return (int) 1 if the point belongs to element "element", 0 otherwise
 */
template<typename _SIZE_>
int Point_32_64<_SIZE_>::contient(const ArrOfDouble& pos, int_t element ) const
{
  assert(pos.size_array()==this->dimension);

  const Domaine_t& dom=this->mon_dom.valeur();
  const IntTab_t& elem=dom.les_elems();
  int ok=1;
  for (int d=0; (d<3)&&(ok==1); d++)
    {
      double ps = dom.coord(elem(element,0), d);
      double pv = pos[d];
      if( !est_egal(ps,pv)) ok=0;
    }
  return ok;
}

/*! @brief Returns 1 if the vertices specified by parameter "pos" are the vertices of element "element" of the domain associated with
 *
 *     the geometric element.
 *
 * @param (IntVect& pos) the vertex indices to compare with those of element "element"
 * @param (int element) the index of the domain element whose vertices are to be compared
 * @return (int) 1 if the vertices passed as parameter are those of the specified element, 0 otherwise
 */
template<typename _SIZE_>
int Point_32_64<_SIZE_>::contient(const SmallArrOfTID_t& pos, int_t element ) const
{
  abort();

  assert(pos.size_array()==1);
  const Domaine_t& dom=this->mon_dom.valeur();
  if((dom.sommet_elem(element,0)==pos[0])&&
      (dom.sommet_elem(element,1)==pos[1]))
    return 1;
  else
    return 0;
}


/*! @brief voir ElemGeomBase::get_tab_faces_sommets_locaux
 */
template<typename _SIZE_>
int Point_32_64<_SIZE_>::get_tab_faces_sommets_locaux(IntTab& faces_som_local) const
{
  faces_som_local.resize(1,1);
  faces_som_local[0]=0;
  return 1;
}



template class Point_32_64<int>;
#if INT_is_64_ == 2
template class Point_32_64<trustIdType>;
#endif

