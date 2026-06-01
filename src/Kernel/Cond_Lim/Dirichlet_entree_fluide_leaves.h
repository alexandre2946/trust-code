/****************************************************************************
* Copyright (c) 2023, CEA
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

#ifndef Dirichlet_entree_fluide_leaves_included
#define Dirichlet_entree_fluide_leaves_included

#include <Dirichlet.h>

/// \cond DO_NOT_DOCUMENT
class Dirichlet_entree_fluide_leaves
{ };
/// \endcond

/*! @brief Dirichlet_entree_fluide This class represents a boundary condition imposing a quantity
 *
 *     at the fluid inlet. Derived classes will specialise the imposed quantity: velocity, temperature, concentration, fraction_massique ...
 *
 * @sa Dirichlet Entree_fluide_vitesse_imposee Entree_fluide_vitesse_imposee_libre, Entree_fluide_temperature_imposee Entree_fluide_T_h_imposee, Entree_fluide_Fluctu_Temperature_imposee Entree_fluide_Flux_Chaleur_Turbulente_imposee, Entree_fluide_K_Eps_impose Entree_fluide_V2_impose, Entree_fluide_concentration_imposee Entree_fluide_fraction_massique_imposee
 */
class Dirichlet_entree_fluide: public Dirichlet
{
  Declare_base(Dirichlet_entree_fluide);
};

/* ========================================================================================= */

/*! @brief Entree_fluide_vitesse_imposee Special case of the class Dirichlet_entree_fluide
 *
 *     for imposed velocity: imposes the inlet velocity of the fluid in an equation of type Navier_Stokes
 *
 * @sa Dirichlet_entree_fluide Navier_Stokes_std
 */
class Entree_fluide_vitesse_imposee: public Dirichlet_entree_fluide
{
  Declare_instanciable(Entree_fluide_vitesse_imposee);
};

/* ========================================================================================= */

/*! @brief Entree_fluide_vitesse_imposee_libre Special case of the class Entree_fluide_vitesse_imposee for imposed velocity:
 *
 *   imposes the fluid velocity in an equation of type Navier_Stokes while leaving the other fields free on this open boundary
 *
 * @sa Dirichlet_entree_fluide Navier_Stokes_std
 */
class Entree_fluide_vitesse_imposee_libre: public Entree_fluide_vitesse_imposee
{
  Declare_instanciable(Entree_fluide_vitesse_imposee_libre);
};

/* ========================================================================================= */

/*! @brief Entree_fluide_concentration_imposee Special case of the class Dirichlet_entree_fluide
 *
 *     for imposed concentration: imposes the inlet concentration of the
 *     fluid in an equation of type Convection_Diffusion_Concentration
 *
 * @sa Dirichlet_entree_fluide Convection_Diffusion_Concentration, CLASS: Entree_fluide_Fluctu_temperature_imposee :, Special case of the class Dirichlet_entree_fluide, for the imposed dissipation rate and temperature variance
 */
class Entree_fluide_Fluctu_Temperature_imposee: public Dirichlet_entree_fluide
{
  Declare_instanciable(Entree_fluide_Fluctu_Temperature_imposee);
};

/* ========================================================================================= */

/*! @brief Entree_fluide_temperature_imposee Special case of the class Dirichlet_entree_fluide for imposed temperature
 *
 *   imposes the inlet temperature of the fluid in an equation of type Convection_Diffusion_Temperature
 *
 * @sa Dirichlet_entree_fluide Convection_Diffusion_Temperature
 */
class Entree_fluide_temperature_imposee: public Dirichlet_entree_fluide
{
  Declare_instanciable(Entree_fluide_temperature_imposee);
};

/* ========================================================================================= */
/*! @brief Entree_fluide_V2_impose Special case of the class Dirichlet_entree_fluide for the imposed velocity fluctuations of the K_Eps_V2 model.
 *
 *     This is the same type of class as Entree_fluide_concentration_imposee, imposing turbulent quantities.
 *     Imposes the V2 inlet values of the fluid in an equation of type Transport_V2
 *
 * @sa Dirichlet_entree_fluide Entree_fluide_concentration_imposee, Transport_V2
 */
class Entree_fluide_V2_impose: public Dirichlet_entree_fluide
{
  Declare_instanciable(Entree_fluide_V2_impose);
};

/* ========================================================================================= */

class Entree_fluide_fraction_massique_imposee: public Dirichlet_entree_fluide
{
  Declare_instanciable(Entree_fluide_fraction_massique_imposee);
};

/* ========================================================================================= */

class Entree_fluide_alpha_impose: public Dirichlet_entree_fluide
{
  Declare_instanciable(Entree_fluide_alpha_impose);
};

/* ========================================================================================= */

class Entree_fluide_Flux_Chaleur_Turbulente_imposee: public Dirichlet_entree_fluide
{
  Declare_instanciable(Entree_fluide_Flux_Chaleur_Turbulente_imposee);
};

/* ========================================================================================= */

#endif /* Dirichlet_entree_fluide_leaves_included */
