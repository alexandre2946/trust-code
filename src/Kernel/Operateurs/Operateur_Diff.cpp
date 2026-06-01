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

#include <Discretisation_base.h>
#include <Operateur_Diff.h>
#include <Champ_base.h>

Implemente_instanciable(Operateur_Diff,"Operateur_Diff",OWN_PTR(Operateur_Diff_base));
// XD diffusion_deriv objet_lecture diffusion_deriv NO_BRACE not_set
// XD bloc_diffusion objet_lecture nul NO_BRACE not_set
// XD attr aco chaine(into=["{"]) aco REQ Opening curly bracket.
// XD attr operateur diffusion_deriv operateur OPT if none is specified, the diffusive scheme used is a 2nd-order
// XD_CONT scheme.
// XD attr op_implicite op_implicite op_implicite OPT To have diffusive implicitation, it use Uzawa algorithm. Very
// XD_CONT useful when viscosity has large variations.
// XD attr acof chaine(into=["}"]) acof REQ Closing curly bracket.

/* Here we put the doc for all shared (VDF, VEF,...) keywords: */

// XD diffusion_negligeable diffusion_deriv negligeable NO_BRACE the diffusivity will not taken in count

// XD diffusion_option diffusion_deriv option NO_BRACE not_set
// XD attr bloc_lecture bloc_lecture bloc_lecture REQ not_set

Sortie& Operateur_Diff::printOn(Sortie& os) const
{
  return Operateur::ecrire(os);
}

Entree& Operateur_Diff::readOn(Entree& is)
{
  return Operateur::lire(is);
}

/*! @brief Types the operator: types as "Op_Diff_"+discretisation() + ("_" or "_Multi_inco_") + inconnue().suffix
 *     Associates the diffusivity to the base operator.
 *
 */
void Operateur_Diff::typer()
{
  Cerr << "Operateur_Diff::typer("<<typ<<")" << finl;
  if (Motcle(typ)==Motcle("negligeable"))
    {
      OWN_PTR(Operateur_Diff_base)::typer("Op_Diff_negligeable");
      valeur().associer_diffusivite(diffusivite());
    }
  else
    {
      assert(la_diffusivite);
      Equation_base& eqn=mon_equation.valeur();
      Nom nom_type= eqn.discretisation().get_name_of_type_for(que_suis_je(),typ,eqn,diffusivite());
      OWN_PTR(Operateur_Diff_base)::typer(nom_type);
      valeur().associer_diffusivite(diffusivite());
    }
  Cerr << valeur().que_suis_je() << finl;
}

/*! @brief Returns the underlying object upcast to Operateur_base
 *
 * @return (Operateur_base&) the underlying object upcast to Operateur_base
 */
Operateur_base& Operateur_Diff::l_op_base()
{
  return valeur();
}
/*! @brief Returns the underlying object upcast to Operateur_base (const version)
 *
 * @return (Operateur_base&) the underlying object upcast to Operateur_base
 */
const Operateur_base& Operateur_Diff::l_op_base() const
{
  return valeur();
}

/*! @brief Call to the underlying object.
 *
 * Adds the contribution of the operator to the array
 *     passed as parameter
 *
 * @param (DoubleTab& donnee) array containing the data on which the operator is applied.
 * @param (DoubleTab& resu) array to which the contribution of the operator is added
 * @return (DoubleTab&) the array containing the result
 */
DoubleTab& Operateur_Diff::ajouter(const DoubleTab& donnee,
                                   DoubleTab& resu) const
{
  statistics().begin_count(STD_COUNTERS::diffusion,statistics().get_last_opened_counter_level()+1);
  DoubleTab& tmp = valeur().ajouter(donnee, resu);
  statistics().end_count(STD_COUNTERS::diffusion);
  return tmp;
}

/*! @brief Call to the underlying object.
 *
 * Initializes the array passed as parameter with the contribution
 *     of the operator.
 *
 * @param (DoubleTab& donnee) array containing the data on which the operator is applied.
 * @param (DoubleTab& resu) array in which the contribution of the operator is stored
 * @return (DoubleTab&) the array containing the result
 */
DoubleTab& Operateur_Diff::calculer(const DoubleTab& donnee,
                                    DoubleTab& resu) const
{
  statistics().begin_count(STD_COUNTERS::diffusion,statistics().get_last_opened_counter_level()+1);
  DoubleTab& tmp = valeur().calculer(donnee, resu);
  statistics().end_count(STD_COUNTERS::diffusion);
  return tmp;
}


/*! @brief Returns the field representing the diffusivity.
 *
 * @return (Champ_Don_base&) the field representing the diffusivity
 */
const Champ_base& Operateur_Diff::diffusivite() const
{
  return la_diffusivite.valeur();
}


/*! @brief Associates the diffusivity to the operator.
 *
 * @param (Champ_Don_base& nu) the field representing the diffusivity
 * @return the field representing the diffusivity
 */
void Operateur_Diff::associer_diffusivite(const Champ_base& nu)
{
  la_diffusivite=nu;
}

void  Operateur_Diff::associer_diffusivite_pour_pas_de_temps(const Champ_base& nu)
{
  valeur().associer_diffusivite_pour_pas_de_temps(nu);
}

void Operateur_Diff::associer_diffusivite_volumique(const Champ_base& champ)
{
  valeur().associer_diffusivite_volumique(champ);
}

/*! @brief Types the operator.
 *
 * @param (Nom& typ) the name representing the type of the operator
 */
void Operateur_Diff::typer(const Nom& un_type)
{
  OWN_PTR(Operateur_Diff_base)::typer(un_type);
}
