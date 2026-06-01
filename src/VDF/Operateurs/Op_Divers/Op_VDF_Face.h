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

#ifndef Op_VDF_Face_included
#define Op_VDF_Face_included

#include <Domaine_VDF.h>
#include <TRUSTTab.h>
#include <TRUSTTabs_forward.h>

class Matrice_Morse;
class Domaine_Cl_VDF;

class Op_VDF_Face
{
public :
  void dimensionner(const Domaine_VDF&, const Domaine_Cl_VDF&, Matrice_Morse&) const;
  void modifier_pour_Cl(const Domaine_VDF&, const Domaine_Cl_VDF&, Matrice_Morse&, DoubleTab&) const;

private:
  void modifier_pour_Cl_(const int , const int , const int ,  Matrice_Morse& ) const;
};

// Internal method for the Op_VDF_Face class!
inline int face_bord_amont2(const Domaine_VDF& le_dom , const int num_face , const int k , const int i)
{
  const int ori = le_dom.orientation(num_face);
  int elem = le_dom.face_voisins(num_face,0);
  if(elem != -1)
    {
      const int face = le_dom.elem_faces(elem, k+i*Objet_U::dimension), elem_bis = le_dom.face_voisins(face,i);
      elem = (elem_bis != -1) ? le_dom.elem_faces(elem_bis, ori+Objet_U::dimension) : -1;
    }
  else Process::exit();

  return elem;
}

#endif /* Op_VDF_Face_included */
