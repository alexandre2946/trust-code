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

#ifndef TroisDto2D_included
#define TroisDto2D_included

#include <Interprete_geometrique_base.h>
#include <Domaine.h>


/*! @brief Class TroisDto2D x->alpha x
 *
 * @sa Interprete Pave, Currently the only type of object recognized by TRUST for meshing a domain is the Pave.
 */
template <typename _SIZE_>
class TroisDto2D_32_64 : public Interprete_geometrique_base_32_64<_SIZE_>
{
  Declare_instanciable_32_64(TroisDto2D_32_64);

public :
  using int_t = _SIZE_;
  using ArrOfInt_t = ArrOfInt_T<_SIZE_>;
  using IntTab_t = IntTab_T<_SIZE_>;
  using DoubleTab_t = DoubleTab_T<_SIZE_>;

  using Domaine_t = Domaine_32_64<_SIZE_>;
  using Bord_t = Bord_32_64<_SIZE_>;
  using Faces_t = Faces_32_64<_SIZE_>;
  using Frontiere_t = Frontiere_32_64<_SIZE_>;
  using Raccord_t = OWN_PTR(Raccord_base_32_64<_SIZE_>);

  Entree& interpreter_(Entree&) override;
  void extraire_2D(const Domaine_t&, Domaine_t&, const Bord_t&,const Nom& , int);

  // OC, 28/10/2004, Modification to be able to extract a boundary of any orientation.
  // Returns the X component of the origin A of the local coordinate system at the boundary, calculated by the class.
  double getXa() { return xa; }
  // Returns the Y component of the origin A of the local coordinate system at the boundary, calculated by the class.
  double getYa() { return ya; }
  // Returns the Z component of the origin A of the local coordinate system at the boundary, calculated by the class.
  double getZa() { return za; }
  // Returns the X component of the first basis vector I of the local coordinate system at the boundary, calculated by the class.
  double getIx() { return Ix; }
  // Returns the Y component of the first basis vector I of the local coordinate system at the boundary, calculated by the class.
  double getIy() { return Iy; }
  // Returns the Z component of the first basis vector I of the local coordinate system at the boundary, calculated by the class.
  double getIz() { return Iz; }
  // Returns the X component of the basis vector J of the local coordinate system at the boundary, calculated by the class.
  double getJx() { return Jx; }
  // Returns the Y component of the basis vector J of the local coordinate system at the boundary, calculated by the class.
  double getJy() { return Jy; }
  // Returns the Z component of the basis vector J of the local coordinate system at the boundary, calculated by the class.
  double getJz() { return Jz; }
  // Returns the X component of the basis vector K of the local coordinate system at the boundary, calculated by the class.
  double getKx() { return Kx; }
  // Returns the Y component of the basis vector K of the local coordinate system at the boundary, calculated by the class.
  double getKy() { return Ky; }
  // Returns the Z component of the basis vector K of the local coordinate system at the boundary, calculated by the class.
  double getKz() { return Kz; }

  inline int& coupe() { return coupe_; }

protected :

  int coupe_ = 1;

private :
  double xa = -100., ya = -100., za = -100.;
  double Ix = -100., Iy = -100., Iz = -100.;
  double Jx = -100., Jy = -100., Jz = -100.;
  double Kx = -100., Ky = -100., Kz = -100.;
};

using TroisDto2D = TroisDto2D_32_64<int>;
using TroisDto2D_64 = TroisDto2D_32_64<trustIdType>;


#endif

