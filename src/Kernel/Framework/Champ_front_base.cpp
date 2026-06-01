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

#include <Champ_front_base.h>
#include <Frontiere_dis_base.h>
#include <Frontiere.h>

Implemente_base_sans_constructeur(Champ_front_base,"Champ_front_base",Field_base);
// XD front_field_base objet_u champ_front_base INHERITS_BRACE Basic class for fields at domain boundaries.

Champ_front_base::Champ_front_base() { temps_defaut = -1.; }
/*! @brief Prints the field name to an output stream
 *
 * @param (Sortie& s) an output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Champ_front_base::printOn(Sortie& s ) const { return s << que_suis_je() << " " << le_nom(); }

/*! @brief DOES NOTHING - to override in derived classes.
 *
 * @param (Entree& s) an input stream
 * @return (Entree&) the input stream
 */
Entree& Champ_front_base::readOn(Entree& s ) { return s ; }


/*! @brief Initialization at the beginning of calculation.
 *
 * Imperatively this method must not use data
 *     external to the equation (coupling). If mettre_a_jour does,
 *     then initializer must not call mettre_a_jour.
 *
 * @return (0 in case of error, 1 otherwise.)
 */
int Champ_front_base::initialiser(double temps, const Champ_Inc_base& inco)
{
  return 1;
}

/*! @brief Associates a discretized boundary with the field.
 *
 * @param (Frontiere_dis_base& fr) the discretized boundary to associate with the field
 */
void Champ_front_base::associer_fr_dis_base(const Frontiere_dis_base& fr)
{
  la_frontiere_dis = fr;
}


/*! @brief DOES NOTHING, to override.
 *
 * This method is called at the beginning of each time step or
 *    sub-time-step, it may possibly use data
 *    external to the equation. It is up to the algorithm to ensure
 *    that this data is relevant...
 *    Calculates the value of the boundary condition at the requested time.
 *
 * @param (double)
 */
void Champ_front_base::mettre_a_jour(double temps)
{
}

/*! @brief Called by Conds_lim::completer. By default does nothing.
 *
 *      To override in unsteady front fields.
 *
 */
void Champ_front_base::fixer_nb_valeurs_temporelles(int nb_cases)
{
}

/*! @brief DOES NOTHING, to override. This method can calculate and store useful data for the
 *
 *    BC, and depending only on the unknown on which it acts
 *    this BC (not from outside). cf Champ_front_contact_VEF for example.
 *    It is called when the unknown is modified.
 *
 * @param (double)
 */
void Champ_front_base::calculer_coeffs_echange(double temps)
{
}

/*! @brief Returns the vector of field values for the given face.
 *
 * By default for func fields, we assume that the values array
 *     holds nb_faces * nb_compo_ values.
 *     Special case example: champ_front_uniforme::valeurs_face
 *
 * @param (num_face) the index of a face on the boundary 0 <= num_face < frontiere_dis().frontiere().nb_faces()
 * @param (val) Resize this array and fill it.
 */
void Champ_front_base::valeurs_face(int num_face, DoubleVect& val) const
{
#ifndef NDEBUG
  const int nb_faces = frontiere_dis().frontiere().nb_faces();
  // If crash here, it means the values array does not contain a
  // value for each face, you need to reimimplement the method in the
  // derived class...
  assert(num_face >= 0 && num_face < nb_faces);
  assert(valeurs().dimension(0) == nb_faces);
#endif
  const int n = nb_compo_;
  val.resize(n);
  const DoubleTab& valeurs_a_copier=valeurs();
  for (int i = 0; i < n; i++)
    val[i] = valeurs_a_copier(num_face, i);
}

const Domaine_dis_base& Champ_front_base::domaine_dis() const
{
  return frontiere_dis().domaine_dis();
}


/*! @brief To implement in derived classes.
 *
 * Advances in time: the new current time will be the time passed
 *      as a parameter.
 *
 * @return (int) 1 if OK, 0 otherwise
 */
int Champ_front_base::avancer(double temps)
{
  Cerr << "Champ_front_base::avancer(double temps) should be overloaded" << finl;
  Process::exit();
  return 0;
}

/*! @brief To implement in derived classes.
 *
 * Rewinds in time: the new current time will be the time passed
 *      as a parameter.
 *
 * @return (int) 1 if OK, 0 otherwise
 */
int Champ_front_base::reculer(double temps)
{
  Cerr << "Champ_front_base::reculer(double temps) should be overloaded " << "by " << que_suis_je() << finl;
  Process::exit();
  return 0;
}

/*! @brief Changes the time value for the i-th temporal value after the present
 *
 */
void Champ_front_base::changer_temps_futur(double temps,int i)
{
  les_valeurs->futur(i).changer_temps(temps);
}

/*! @brief Computes the rate of change of the field between t1 and t2 and stores it in Gpoint_
 *
 */
void Champ_front_base::calculer_derivee_en_temps(double t1, double t2)
{
  if (std::abs(t2-t1) < DMINFLOAT)
    {
      Gpoint_ = 0;
    }
  else
    {
      const DoubleTab& v1 = valeurs_au_temps(t1);
      const DoubleTab& v2 = valeurs_au_temps(t2);
      if (!Gpoint_.get_md_vector() && v1.dimension(0) == 1)
        {
          // Uniform unsteady field
          int dim = v1.dimension(1);
          Gpoint_.resize(dim);
          for (int i = 0; i < dim; i++)
            Gpoint_(i) = (v2(0, i) - v1(0, i)) / (t2 - t1);
        }
      else
        {
          // Variable unsteady field
          Gpoint_ = v1;
          Gpoint_ *= -1;
          Gpoint_ += v2;
          Gpoint_ /= (t2 - t1);
        }
    }
}

