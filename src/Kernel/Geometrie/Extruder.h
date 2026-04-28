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

#ifndef Extruder_included
#define Extruder_included

#include <Faces_builder.h>
#include <Interprete_geometrique_base.h>
#include <Static_Int_Lists.h>
#include <Domaine_forward.h>
#include <Domaine_forward.h>
#include <Bord.h>

/*! @brief Classe Extruder Cette classe est un interprete qui sert a lire et executer
 *
 *     la directive Extruder:
 *         Extruder nom_domaine
 *     Cette directive est a utiliser en discretisation VEF 2D pour obtenir
 *     un maillage 3D par extrusion puis decoupage.
 *
 * @sa Interprete Extruder, Cette classe est utilisable en 3D
 */

template <typename _SIZE_>
class Extruder_32_64 : public Interprete_geometrique_base_32_64<_SIZE_>
{
  Declare_instanciable_sans_constructeur_32_64(Extruder_32_64);

public :
  using int_t = _SIZE_;
  using Domaine_t = Domaine_32_64<_SIZE_>;
  using Faces_t = Faces_32_64<_SIZE_>;
  using Bord_t =  Bord_32_64<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;
  using Faces_builder_t = Faces_builder_32_64<_SIZE_>;
  using Static_Int_Lists_t = Static_Int_Lists_32_64<_SIZE_>;

  Extruder_32_64();

  Entree& interpreter_(Entree&) override;
  inline void setDirection(double lx, double ly, double lz)
  {
    direction[0]=lx;
    direction[1]=ly;
    direction[2]=lz;
  }
  inline void setNbTranches(int n) { NZ = n; }

  void extruder(Domaine_t&) ;


protected:
  void extruder_hexa(Domaine_t&);
  virtual void extruder_dvt(Domaine_t&, Faces_t&, int_t, int_t ) ;
  virtual void extruder_dvt_hexa(Domaine_t&, Faces_t&, int_t , int_t ) ;

  ArrOfDouble direction;
  int NZ = -10;

private:
  void traiter_faces_dvt_hexa(Faces_t&, int_t);
  void traiter_faces_dvt(Faces_t&, Faces_t&, int_t, int_t, int_t);
};

using Extruder = Extruder_32_64<int>;
using Extruder_64 = Extruder_32_64<trustIdType>;

#endif
