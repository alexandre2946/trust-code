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

#ifndef ExtrudeBord_included
#define ExtrudeBord_included

#include <Interprete_geometrique_base.h>
#include <TRUSTTabs_forward.h>


/*! @brief class ExtrudeBord
 *
 */
/**
 * Performs an extrusion of a boundary.
 * ExtrudeBord can extrude a boundary meshed either with triangles or quadrangles:
 *  - if the boundary is meshed with quadrangles, the extruded domain will be composed of hexahedra
 *  - if the boundary is meshed with triangles, the extruded domain will be composed of tetrahedra
 * The syntax is as follows:
 *         ExtrudeBord
 *                 {
 *                 domaine_init         NAME_OF_SOURCE_DOMAIN
 *                 direction         X Y Z
 *                 nb_tranches         N
 *                 domaine_final         NAME_OF_EXTRUDED_DOMAIN
 *                 nom_bord         NAME_OF_BOUNDARY_TO_EXTRUDE
 *                hexa_old  //keyword to use the old version of hexa extrusion
 *                 }
 *
 * (MODIF OC 12/2004)
 */
class ExtrudeBord : public Interprete_geometrique_base
{
  Declare_instanciable(ExtrudeBord);

public :

  Entree& interpreter_(Entree&) override;

private:

  void extruder_bord(Nom& nom_front, Nom& nom_dom2, DoubleVect& vect_dir, int nbpas);
  void extruder_hexa_old(Nom& nom_front, Nom& nom_dom2, DoubleVect& vect_dir, int nbpas);
  bool hexa_old=false;    // flag for old version of hexa extrusion: 0 = old version
  bool Trois_Tetra=false; // flag for extrusion into three tetrahedra rather than 14 (default option)
  bool Vingt_Tetra=false; // flag for extrusion into twenty tetrahedra rather than 14 (default option)
  int en3D_=1;
};


#endif

