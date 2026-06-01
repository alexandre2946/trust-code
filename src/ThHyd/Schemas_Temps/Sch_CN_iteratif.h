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

#ifndef Sch_CN_iteratif_included
#define Sch_CN_iteratif_included

#include <Schema_Temps_base.h>
#include <TRUSTTabs_forward.h>

/*! @brief Time scheme alternating a half implicit Euler step and a half LeapFrog step.
 *
 *      The implicit solve is iterative (fixed-point).
 *      The time step is computed as the product of the explicit stability time step by a facsec.
 *      The facsec is adjusted automatically so that the solve converges in a predefined number of iterations.
 *      Characteristics of each iteration are written to the file dt_CN.
 *
 *      The solve is governed by 4 parameters (default values in parentheses):
 *      * seuil (1e-3) : convergence threshold. The lower, the more accurate the solve.
 *      * facsec_max (2) : maximum facsec value (to avoid instabilities and capture physical phenomena).
 *      * niter_min (2) : minimum number of iterations. Below this, iteration continues even if convergence seems reached.
 *      * niter_avg (3) : target number of iterations to reach convergence.
 *      * niter_max (6) : number of iterations beyond which a smaller facsec is retried.
 *
 *      Advice for choosing facsec adjustment parameters:
 *      * Choose seuil based on desired precision.
 *      * Choose niter_min: 2 guarantees a second-order time scheme.
 *      * If seeking a steady state, choose seuil_statio >= seuil.
 *      * Choose facsec_max based on the physical phenomena to capture.
 *      * Start by testing with a large value of niter_avg. Observe the number of iterations.
 *        It will plateau at a maximum before dropping back.
 *      * Choose niter_avg at 2/3 of this maximum, and niter_max at approximately 4/3 or double.
 *
 *
 * @sa Schema_Temps_base
 */
class Sch_CN_iteratif : public Schema_Temps_base
{
  Declare_instanciable(Sch_CN_iteratif);

public :

  ////////////////////////////////
  //                            //
  // Scheme characteristics     //
  //                            //
  ////////////////////////////////

  int nb_valeurs_temporelles() const override;
  int nb_valeurs_futures() const override;
  double temps_futur(int i) const override;
  double temps_defaut() const override;

  /////////////////////////////////////////
  //                                     //
  // End of scheme characteristics       //
  //                                     //
  /////////////////////////////////////////

  bool initTimeStep(double dt) override;
  bool iterateTimeStep(bool& converged) override;

  int faire_un_pas_de_temps_eqn_base(Equation_base&) override;

  void completer() override {}
  void set_param(Param& titi) const override;

protected :

  enum type_convergence {DIVERGENCE, NON_CONVERGENCE, CONVERGENCE_LENTE, CONVERGENCE_RAPIDE, CONVERGENCE_OK}; // slow/fast convergence

  virtual bool convergence(const DoubleTab& u0, const DoubleTab& up1, const DoubleTab& delta, int p) const;
  virtual bool divergence(const DoubleTab& u0, const DoubleTab& up1, const DoubleTab& delta, int p) const;
  virtual void ajuster_facsec(type_convergence cv);

  virtual bool iterateTimeStepOnEquation(int i,bool& converged);

  double seuil=1.e-3; // To determine convergence
  int niter_min=2; // Minimum number of iterations (before, continue to iterate)
  int niter_max=6; // Maximum number of iterations (after, considered as not convergent)
  int niter_avg=3; // Average number of iterations wanted (facsec adjusted to fit that number)
  double facsec_max=2; // Maximum facsec (not to miss the physics)

  // Used internally.
  int iteration = -1; // Number of iterations done for the current time step.
  double last_facsec = -100.; // facsec at the beginning of time step resolution,
  // to avoid changing facsec several times for the same time step.
};

#endif
