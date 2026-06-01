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

#ifndef Op_Diff_negligeable_included
#define Op_Diff_negligeable_included

#include <Operateur_negligeable.h>
#include <Operateur_Diff_base.h>
#include <TRUST_Ref.h>

class Champ_base;

/*! @brief Classe Op_Diff_negligeable This class represents a negligible diffusion operator.
 *
 *     When an operator of this type is used in an equation,
 *     it amounts to neglecting the diffusion term.
 *     The methods for modifying and participating in a computation of
 *     the operator are in fact calls to the same methods of
 *     Operateur_negligeable which do nothing.
 *
 *
 * @sa Operateur_Diff_base Operateur_negligeable
 */

class Op_Diff_negligeable: public Operateur_negligeable,
  public Operateur_Diff_base
{
  Declare_instanciable(Op_Diff_negligeable);

public :

  inline DoubleTab& ajouter(const DoubleTab&, DoubleTab& ) const override;
  inline DoubleTab& calculer(const DoubleTab&, DoubleTab& ) const override;
  inline void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const override;
  inline void contribuer_au_second_membre(DoubleTab& ) const override;
  inline void modifier_pour_Cl(Matrice_Morse&, DoubleTab&) const override;
  inline void dimensionner(Matrice_Morse& ) const override;
  inline void mettre_a_jour(double) override;
  void associer_diffusivite(const Champ_base& ) override ;
  const Champ_base& diffusivite() const override;
  inline void associer_champ_masse_volumique(const Champ_base&) override;
  void calculer_pour_post(Champ_base& espace_stockage,const Nom& option, int comp) const override;
  Motcle get_localisation_pour_post(const Nom& option) const override;

  void ajouter_flux(const DoubleTab& inconnue, DoubleTab& contribution) const override;
  void calculer_flux(const DoubleTab& inconnue, DoubleTab& flux) const override;

  int has_interface_blocs() const override
  {
    return 1;
  };
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = { }) const override { };
  void ajouter_blocs(matrices_t matrices, DoubleTab& resu, const tabs_t& semi_impl = { }) const override { };
  void check_multiphase_compatibility() const override { };

protected :

  OBS_PTR(Champ_base) la_diffusivite;
  inline void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&, const Champ_Inc_base& ) override ;
};

class Op_Dift_negligeable: public Op_Diff_negligeable
{
  Declare_instanciable(Op_Dift_negligeable);
};


/*! @brief Adds the contribution of the operator to an array passed as parameter.
 *
 * Simple call to Operateur_negligeable::ajouter(const DoubleTab&,DoubleTab&)
 *
 * @param (DoubleTab& x) the array on which the operator is applied
 * @param (DoubleTab& y) array to which the operator contribution is added
 * @return (DoubleTab&) the input parameter y unchanged
 */
inline DoubleTab&
Op_Diff_negligeable::ajouter(const DoubleTab& x, DoubleTab& y) const
{
  return Operateur_negligeable::ajouter(x,y);
}


/*! @brief Initializes the array parameter with the contribution of the negligible operator: initializes the array to ZERO.
 *
 *     Simple call to Operateur_negligeable::(calculer(const DoubleTab&, DoubleTab&)
 *
 * @param (DoubleTab& x) the array on which the operator is applied
 * @param (DoubleTab& y) array in which the operator contribution is stored
 * @return (DoubleTab&) the input array y set to zero
 */
inline DoubleTab&
Op_Diff_negligeable::calculer(const DoubleTab& x, DoubleTab& y) const
{
  return Operateur_negligeable::calculer(x,y);
}

/*! @brief Time update of a negligible operator: DOES NOTHING. Simple call to Operateur_negligeable::mettre_a_jour(double)
 *
 * @param (double temps) the current time
 */

/*! @brief Assembles the matrix.
 *
 */

inline void Op_Diff_negligeable::contribuer_a_avec(const DoubleTab& inco,
                                                   Matrice_Morse& amatrice) const
{
  ;
}

/*! @brief Adds the contribution to the right-hand side.
 *
 */

inline void Op_Diff_negligeable::contribuer_au_second_membre(DoubleTab& resu) const
{
  ;
}

// Modification of boundary conditions
inline void  Op_Diff_negligeable::modifier_pour_Cl(Matrice_Morse& amatrice, DoubleTab& resu) const
{
  ;
}

inline void Op_Diff_negligeable::dimensionner(Matrice_Morse& amatrice) const
{
  ;
}

inline void Op_Diff_negligeable::mettre_a_jour(double temps)
{
  Operateur_negligeable::mettre_a_jour(temps);
}


/*! @brief Associates various objects to a negligible operator: DOES NOTHING Simple call to Operateur_negligeable::associer(const Domaine_dis_base&,
 *
 *                                                      const Domaine_Cl_dis_base&,
 *                                                      const Champ_Inc_base&)
 *
 * @param (Domaine_dis_base& z)
 * @param (Domaine_Cl_dis_base& zcl)
 * @param (Champ_Inc_base& ch)
 */
inline void Op_Diff_negligeable::associer(const Domaine_dis_base& z,
                                          const Domaine_Cl_dis_base& zcl,
                                          const Champ_Inc_base& ch)
{
  Operateur_negligeable::associer(z, zcl, ch);
}

// Override so that the method does nothing
inline void Op_Diff_negligeable::associer_champ_masse_volumique(const Champ_base&)
{

}

#endif

