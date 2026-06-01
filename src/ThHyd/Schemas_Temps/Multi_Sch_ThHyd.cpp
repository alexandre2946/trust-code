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

#include <Multi_Sch_ThHyd.h>
#include <Probleme_base.h>
#include <Navier_Stokes_std.h>
#include <Param.h>

Implemente_instanciable(Multi_Sch_ThHyd,"Multi_Schema_ThHyd",Schema_Temps_base);

Sortie& Multi_Sch_ThHyd::printOn(Sortie& s) const
{
  return  Schema_Temps_base::printOn(s);
}

Entree& Multi_Sch_ThHyd::readOn(Entree& s)
{
  return Schema_Temps_base::readOn(s) ;
}

/*! @brief Returns the number of time values to keep.
 *
 * Here: the max of the two schemes used.
 *
 */
int Multi_Sch_ThHyd::nb_valeurs_temporelles() const
{
  return std::max(sch_ns_->nb_valeurs_temporelles(),sch_scalaires_->nb_valeurs_temporelles());
}

/*! @brief Returns the number of future time values.
 *
 * Here: the value common to the two schemes used.
 *
 */
int Multi_Sch_ThHyd::nb_valeurs_futures() const
{
  int n=sch_ns_->nb_valeurs_futures();
  assert (n==sch_scalaires_->nb_valeurs_futures());
  return n;
}

/*! @brief Returns the time at the i-th future value.
 *
 * Here: the value common to the two schemes used.
 *
 */
double Multi_Sch_ThHyd::temps_futur(int i) const
{
  double t=sch_ns_->temps_futur(i);
  assert(t==sch_scalaires_->temps_futur(i));
  return t;
}

/*! @brief Returns the time that fields should return when calling valeurs().
 *
 *     Here: the value common to the two schemes used.
 *
 */
double Multi_Sch_ThHyd::temps_defaut() const
{
  double t=sch_ns_->temps_defaut();
  assert(t==sch_scalaires_->temps_defaut());
  return t;
}

/////////////////////////////////////////
//                                     //
// End of scheme characteristics       //
//                                     //
/////////////////////////////////////////


void Multi_Sch_ThHyd::completer()
{
  // OC: I don't understand these lines:
  /*  double dts=dt_;
      dt_=0;
      mettre_a_jour();
      nb_pas_dt_--;
      dt_=dts;
  */

  // OC: I would rather have a "completer" method of this kind:

  Schema_Temps_base& le_sch_ns = sch_ns_.valeur();
  le_sch_ns.set_temps_init()=temps_init();
  le_sch_ns.set_temps_max()=temps_max();
  le_sch_ns.set_temps_courant()=temps_courant();
  le_sch_ns.set_nb_pas_dt()=nb_pas_dt();
  le_sch_ns.set_nb_pas_dt_max()=nb_pas_dt_max();
  le_sch_ns.set_dt_min()=pas_temps_min();
  le_sch_ns.set_dt_max()=pas_temps_max();
  le_sch_ns.set_dt_sauv()=temps_sauv();
  le_sch_ns.set_dt_impr()=temps_impr();
  le_sch_ns.set_facsec()=facteur_securite_pas();
  le_sch_ns.set_seuil_statio()=seuil_statio();
  le_sch_ns.set_stationnaire_atteint()=isStationary();
  le_sch_ns.set_diffusion_implicite()=diffusion_implicite();
  le_sch_ns.set_seuil_diffusion_implicite()=seuil_diffusion_implicite();
  le_sch_ns.set_niter_max_diffusion_implicite()=niter_max_diffusion_implicite();
  le_sch_ns.set_dt()=pas_de_temps();
  le_sch_ns.set_mode_dt_start()=mode_dt_start();
  le_sch_ns.set_indice_tps_final_atteint()=indice_tps_final_atteint();
  le_sch_ns.set_indice_nb_pas_dt_max_atteint()=indice_nb_pas_dt_max_atteint();
  le_sch_ns.set_lu()=lu();

  Schema_Temps_base& le_sch_scalaires = sch_scalaires_.valeur();
  le_sch_scalaires.set_temps_init()=temps_init();
  le_sch_scalaires.set_temps_max()=temps_max();
  le_sch_scalaires.set_temps_courant()=temps_courant();
  le_sch_scalaires.set_nb_pas_dt()=nb_pas_dt();
  le_sch_scalaires.set_nb_pas_dt_max()=nb_pas_dt_max();
  le_sch_scalaires.set_dt_min()=pas_temps_min();
  le_sch_scalaires.set_dt_max()=pas_temps_max();
  le_sch_scalaires.set_dt_sauv()=temps_sauv();
  le_sch_scalaires.set_dt_impr()=temps_impr();
  le_sch_scalaires.set_facsec()=facteur_securite_pas();
  le_sch_scalaires.set_seuil_statio()=seuil_statio();
  le_sch_scalaires.set_stationnaire_atteint()=isStationary();
  le_sch_scalaires.set_diffusion_implicite()=diffusion_implicite();
  le_sch_scalaires.set_seuil_diffusion_implicite()=seuil_diffusion_implicite();
  le_sch_scalaires.set_niter_max_diffusion_implicite()=niter_max_diffusion_implicite();
  le_sch_scalaires.set_dt()=pas_de_temps();
  le_sch_scalaires.set_mode_dt_start()=mode_dt_start();
  le_sch_scalaires.set_indice_tps_final_atteint()=indice_tps_final_atteint();
  le_sch_scalaires.set_indice_nb_pas_dt_max_atteint()=indice_nb_pas_dt_max_atteint();
  le_sch_scalaires.set_lu()=lu();
}

int Multi_Sch_ThHyd::mettre_a_jour()
{
  Schema_Temps_base::mettre_a_jour();
  sch_ns_->mettre_a_jour();
  sch_scalaires_->mettre_a_jour();
  return 1;
}
void Multi_Sch_ThHyd::set_param(Param& param) const
{
  param.ajouter("nb_ss_pas_dt",&nb_ss_pas_dt_);
  param.ajouter("Schema_Temps_NS",&sch_ns_);
  param.ajouter("Schema_Temps_scalaires",&sch_scalaires_);
  Schema_Temps_base::set_param(param);
}


/*! @brief Performs one explicit Euler time step on the given equation.
 *
 * @param eqn The equation to advance by one time step.
 * @return Always returns 1.
 */
int Multi_Sch_ThHyd::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{
  //  double dtok=dt_;
  if(sub_type(Navier_Stokes_std, eqn))
    {
      sch_ns_->set_dt()=pas_de_temps();
      sch_ns_->faire_un_pas_de_temps_eqn_base(eqn);
      facteur_securite_pas()=sch_ns_->facteur_securite_pas();
    }
  else
    {
      //sch_scalaires.preparer_pas_temps();
      sch_scalaires_->set_dt()=pas_de_temps();
      sch_scalaires_->faire_un_pas_de_temps_eqn_base(eqn);
    }
  set_stationnaire_atteint()=sch_ns_->isStationary() && sch_scalaires_->isStationary() ;
  return 1;
}

bool Multi_Sch_ThHyd::iterateTimeStep(bool& converged)
{
  Probleme_base& prob=pb_base();
  double temps=temps_courant_+dt_;
  int nb_eqn=pb_base().nombre_d_equations();
  for(int i=nb_eqn-1; i>-1; i--)
    {
      Equation_base& eqn_i=prob.equation(i);
      if (eqn_i.equation_non_resolue())
        {
          Cout<< "====================================================" << finl;
          Cout<< eqn_i.que_suis_je()<<" equation is not solved."<<finl;
          Cout<< "====================================================" << finl;
          // Compute the derivative once to obtain the boundary fluxes
          if (eqn_i.schema_temps().nb_pas_dt()==0)
            {
              DoubleTab inconnue_valeurs(eqn_i.inconnue().valeurs());
              eqn_i.derivee_en_temps_inco(inconnue_valeurs);
            }
        }
      else
        {
          eqn_i.domaine_Cl_dis().mettre_a_jour(temps);
          faire_un_pas_de_temps_eqn_base(eqn_i);
        }
    }
  converged=true;
  return true;

}

/*! @brief Corrects the time step so that dt_min <= dt <= dt_max.
 *
 * @return Returns the result of corriger_pas_temps from the base class.
 */
bool Multi_Sch_ThHyd::corriger_dt_calcule(double& dt) const
{
  bool ok=sch_ns_->corriger_dt_calcule(dt);
  ok = ok && sch_scalaires_->corriger_dt_calcule(dt);
  ok = ok && Schema_Temps_base::corriger_dt_calcule(dt);
  return ok;
}


/*! @brief Call to the underlying object. Changes the current time.
 *
 * @param t The new value of the current time.
 */
void Multi_Sch_ThHyd::changer_temps_courant(const double t)
{
  sch_ns_->changer_temps_courant(t);
  sch_scalaires_->changer_temps_courant(t);

  Schema_Temps_base::changer_temps_courant(t);
}

/*! @brief Call to the underlying object. Returns 1 if the calculation should be stopped for various reasons:
 *
 *         - the final time has been reached
 *         - the maximum number of time steps has been exceeded
 *         - the stationary state has been reached
 *         - a file-based stop indicator
 *     Returns 0 otherwise.
 *
 * @return 1 if the calculation should stop, 0 otherwise.
 */
int Multi_Sch_ThHyd::stop() const
{
  int ls2 = sch_ns_->stop();
  int ls3 = sch_scalaires_->stop();

  return (ls2 | ls3 | Schema_Temps_base::stop());
}

/*! @brief Call to the underlying object. Prints the time scheme to an output stream (if applicable).
 *
 * @param os The output stream.
 */
void Multi_Sch_ThHyd::imprimer(Sortie& os) const
{
  sch_ns_->imprimer(os);
  sch_scalaires_->imprimer(os);

  Schema_Temps_base::imprimer(os);
}
