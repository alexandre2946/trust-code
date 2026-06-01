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
#ifndef Ensemble_Lagrange_base_included
#define Ensemble_Lagrange_base_included

#include <Domaine.h>

class Equation_base;

/*! @brief class Ensemble_Lagrange_base Base class for classes representing a geometric structure made up
 *
 *      of a set of points tracked using Lagrangian tracking.
 *      Currently only one instantiable derived class: Maillage_FT_Disc.
 *      A Lagrangian set is characterized by the coordinates of its points.
 *
 * @sa Abstract class.
 * @sa Abstract methods:
 *       void associer_equation_transport)
 *       Equation_base& equation_associee()
 */

class Ensemble_Lagrange_base : public Objet_U
{
  Declare_base_sans_constructeur(Ensemble_Lagrange_base);

public :
  Ensemble_Lagrange_base();
  virtual void associer_equation_transport(const Equation_base&
                                           equation) = 0;
  virtual const Equation_base& equation_associee() const = 0;
  void associer_domaine(const Domaine& domaine);
  virtual inline const Domaine& domaine() const;
  void remplir_sommets_tmp(DoubleTab& soms_tmp);
  void generer_marqueurs_sz(DoubleTab& soms_tmp);
  inline const IntVect& nb_marqs_par_sz() const;
  inline const DoubleTab& sommets_lu() const;
  inline DoubleTab& sommets_lu();

protected :

  Noms nom_sz;                //name of the sub-domains where particles are generated
  IntVect nb_marqs_sz;        //number of markers per sub-domain
  IntTab nb_marqs_par_dir;        //number of markers in each direction of a sub-domain
  //if uniform distribution over the sub-domain
  OBS_PTR(Domaine) mon_dom_;                //REF to the Eulerian mesh domain

  DoubleTab sommets_lu_;      //coordinates of vertices read in the case of reading from a file

private :

};

inline const Domaine& Ensemble_Lagrange_base::domaine() const
{
  return mon_dom_.valeur();
}

inline const IntVect& Ensemble_Lagrange_base::nb_marqs_par_sz() const
{
  return nb_marqs_sz;
}

inline const DoubleTab& Ensemble_Lagrange_base::sommets_lu() const
{
  return sommets_lu_;
}

inline DoubleTab& Ensemble_Lagrange_base::sommets_lu()
{
  return sommets_lu_;
}

#endif
