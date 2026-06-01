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

#ifndef Domaine_Cl_EF_included
#define Domaine_Cl_EF_included



/*! @brief class Domaine_Cl_EF
 *
 * @brief This class holds the arrays used to apply boundary conditions
 *   in the EF formulation.
 *
 *
 * @sa Domaine_Cl_dis_base
 */
#include <Domaine_Cl_dis_base.h>

class Champ_Don_base;
class Domaine_EF;
class Matrice_Morse;

class Domaine_Cl_EF : public Domaine_Cl_dis_base
{

  Declare_instanciable(Domaine_Cl_EF);

public :

  void completer(const Domaine_dis_base& ) override;
//  void mettre_a_jour(double );
  int initialiser(double temps) override;
  void imposer_cond_lim(Champ_Inc_base&, double) override;

  void imposer_symetrie(DoubleTab&,int tous_les_sommets_sym=0) const;
  void imposer_symetrie_partiellement(DoubleTab&,const Noms&) const;
  void imposer_symetrie_matrice_secmem(Matrice_Morse& la_matrice, DoubleTab& secmem) const;

  void modifie_gradient(ArrOfDouble& grad_mod,const ArrOfDouble& grad,int num_som) const;


  int nb_faces_sortie_libre() const;
  Domaine_EF& domaine_EF();
  const Domaine_EF& domaine_EF() const;

  int nb_bord_periodicite() const;
  inline  const ArrOfInt& get_type_sommet() const
  {
    return type_sommet_ ;
  };
protected:

  // Attributes:

  int modif_perio_fait_=0;

  // Functions for creating private domain members:

  void remplir_type_elem_Cl(const Domaine_EF& );
  ArrOfInt type_sommet_;  // -1 internal, 0 Neumann, 1 Symmetry, >2 Dirichlet
  // A vertex is first Dirichlet, then Symmetry, then Neumann, then internal
  OWN_PTR(Champ_Don_base) normales_symetrie_,normales_symetrie_bis_,normales_symetrie_ter_;
};

//
// Inline functions of class Domaine_Cl_EF
//





#endif
