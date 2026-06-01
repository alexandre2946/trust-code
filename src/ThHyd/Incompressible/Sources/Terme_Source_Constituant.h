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

#ifndef Terme_Source_Constituant_included
#define Terme_Source_Constituant_included

#include <Equation_base.h>
#include <TRUST_Ref.h>

#include <SFichier.h>

class Champ_base;
class Champ_Don_base;

/*! @brief Classe Terme_Source_Constituant This class represents a source term of the constituent transport equation
 *
 *     of the volumetric thermal power release type.!!!!TO MODIFY
 *     A Terme_Source_Constituant object contains the power (OWN_PTR(Champ_base) given
 *     by the user) and references to the density (rho) and the specific heat (Cp).
 *
 * @sa Non-instantiable class., The implementation of the terms will depend on their discretization.
 */
class Terme_Source_Constituant
{

public :

  inline void associer_champs(const Champ_base&);
  void lire_donnees(Entree& );
  void ouvrir_fichier(const Equation_base& eq, const Nom& out, const Nom& qsj, const Nom& description, SFichier& os,const Nom& type, const int flag) const;
  int completer(const Champ_Inc_base& inco);

  inline const Champ_Don_base& get_source() const
  {
    return la_source_constituant ;
  };

  void mettre_a_jour(double temps)
  {
    la_source_constituant->mettre_a_jour(temps);
  };

protected:
  int colw_;
  OBS_PTR(Champ_base) rho_ref;
  OWN_PTR(Champ_Don_base) la_source_constituant;

};


/*! @brief Associates the given fields rho (density) and Cp (specific heat) to the object.
 *
 * @param (Champ_Don_base& rho) given field representing density
 * @param (Champ_Don_base& cp) given field representing specific heat
 */
inline void Terme_Source_Constituant::associer_champs(const Champ_base& rho)
{
  rho_ref=rho;
}

#endif
