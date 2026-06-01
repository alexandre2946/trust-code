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

#ifndef Robin_VEF_included
#define Robin_VEF_included

#include <Cond_lim_base.h>
#include <Champ_front_txyz.h>

/*! @brief Class Robin_VEF for Robin boundary conditions
 *
 *    Robin boundary conditions functions are set as it follows :
 *      g = \alpha( \nu  \dn u \cdot n - p ) + u \cdot n
 *      xi = \beta( \nu  \dn u \times n)     + u \times n
 *
 *      with \dn = \nabla \cdot n
 *
 *      The parameters \alpha, \beta, and the functions g and xi comes from the data file e.g.
 *
 *      Conditions_limites {
 *                          Robin_boudary_name {
 *                                                  alpha val
 *                                                  beta val
 *                                                  champ_front_normal_et_tangentiel_robin champ_front_fonc_xyz 4 (in 3D, 2 in 2D) normal_field1D, tangential_field_X, tangential_field_Y, tangential_field_Z (tangential_fiel_1D for 2D)
 *                                              }
 *                          }
 *
 * @sa Cond_lim_base Robin_VEF
 */

class Robin_VEF : public Cond_lim_base// , Champ_front_txyz
{

  Declare_instanciable( Robin_VEF ) ;

public :

  void completer() override;

  inline  double get_alpha_cl() const { return alpha_robin_cl_; } ;
  inline double get_beta_cl() const { return beta_robin_cl_; };

  // get the normal and the tangential flux
  double flux_robin_normal_et_trangentiel_imp(int i, int j) const;

  // get the normal flux
  double flux_normal_imp(int i) const;

  // get the tangential flux
  double flux_tangentiel_imp(int i, int j ) const;

  // udapte time data
  void mettre_a_jour(double temps) override;

  // for domain decomposition, not implemented yet
  double flux_robin_imp_au_temps(double temps, int i) const ;
  double flux_robin_imp_au_temps(double temps, int i, int j ) const;

  const DoubleTab& flux_robin_normal_imp() const;
  const DoubleTab& flux_robin_tangentiel_imp() const;

  inline double increment_pression_bord(int face) const { return 0.; } ;

protected :

  double alpha_robin_cl_ = -123.  ;
  double beta_robin_cl_ = -123.   ;
  // Champ_front_normal_robin champ_normal_robin_cl_;
  // Champ_front_tangentiel_robin champ_tangent_robin_cl_;
  mutable DoubleTab flux_normal_impose_; // Stores all flux values on all boundary faces (no assumption of a uniform field). Useful for GPU.
  mutable DoubleTab flux_tangentiel_impose_; // Stores all flux values on all boundary faces (no assumption of a uniform field). Useful for GPU.
};

#endif /* Robin_VEF_included */
