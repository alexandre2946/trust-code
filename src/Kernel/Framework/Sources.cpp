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

#include <Matrice_Morse.h>
#include <TRUSTTab.h>
#include <Sources.h>

Implemente_instanciable(Sources, "Sources", LIST(Source));
// XD sources listobj sources INHERITS_BRACE source_base INHERITS_COMMA The sources.

Sortie& Sources::printOn(Sortie& os) const { return LIST(Source)::printOn(os); }

/*! @brief Reading of a list of sources from an input stream.
 *
 * Reads a list of sources separated by commas
 *     and adds it to the list.
 *     format:
 *     {
 *      reading block of a source
 *      [, reading block of a source]
 *      ...
 *     }
 *
 * @param (Entree& is) the input stream
 * @return (Entree&) the modified input stream
 * @throws opening brace expected
 * @throws closing brace or comma expected
 */
Entree& Sources::readOn(Entree& is)
{
  Nom accouverte="{";
  Nom accfermee="}";
  Nom virgule=",";
  Nom typ;
  is >> typ;
  if (typ!=accouverte)
    {
      Cerr << "We were waiting for { before " << typ << finl;
      exit();
    }
  Source t;
  is >> typ;
  if(typ==accfermee)
    return is;
  while(1)
    {
      Source& so=add(t);
      Cerr << "Reading of term " << typ << finl;
      Cerr << "and typing: ";
      // Cout << ">>>>>>>>>>>>>>>>>>>> Type de Source = " << typ << finl;
      so.typer(typ,mon_equation.valeur());
      so->associer_eqn(mon_equation.valeur());
      is >> so.valeur();
      is >> typ;
      if(typ==accfermee)
        return is;
      if(typ!=virgule)
        {
          Cerr << typ << " : we expected a ',' or a '}'" << finl;
          exit();
        }
      assert (typ==virgule);
      is >> typ;
    }
}


/*! @brief Adds the contribution of all sources in the list to the array passed as parameter, and returns this array.
 *
 * @param (DoubleTab& xx) the array in which the contribution of the source terms should be accumulated
 * @return (DoubleTab&) the modified parameter xx
 */
DoubleTab& Sources::ajouter(DoubleTab& xx) const
{
  for (const auto& itr : *this) itr.ajouter(xx);
  return xx;
}

/*! @brief Calculates the contribution of all sources in the list and stores the result in the array passed as parameter,
 *
 *     and returns this array.
 *
 * @param (DoubleTab& xx) the array in which the sum of source contributions should be stored
 * @return (DoubleTab&) the modified parameter xx
 */
DoubleTab& Sources::calculer(DoubleTab& xx) const
{
  xx = 0.;
  for (const auto& itr : *this) itr.ajouter(xx);
  return xx;
}

/*! @brief Time update of all sources in the list
 *
 * @param (double temps) the time step for update
 */
void Sources::mettre_a_jour(double temps)
{
  for (auto& itr : *this) itr->mettre_a_jour(temps);
}

/*! @brief Rest all sources to a given time
 *  See ProblemeTrio::resetTime()
 */
void Sources::resetTime(double temps)
{
  for (auto& itr : *this) itr->resetTime(temps);
}

/*! @brief Calls Source::completer() on all sources in the list.
 *
 * See Source_base::completer().
 *
 */
void Sources::completer()
{
  for (auto& itr : *this) itr->completer();
}

/*! @brief For each source in the list, calls the associer_champ_rho method of the source.
 *
 *   If the density is variable, the density field must be declared
 *   to the sources using this method (front-tracking)
 *   Otherwise, by default, calculations are performed with rho=1
 *
 */
void Sources::associer_champ_rho(const Champ_base& champ_rho)
{
  for (auto& itr : *this)
    {
      Source& src = itr;
      Source_base& src_base = src.valeur();
      src_base.associer_champ_rho(champ_rho);
    }
}

/*! @brief For each source in the list, calls a_pour_Champ_Fonc(mot,ch_ref).
 *
 * This method is called by Equation_base::a_pour_Champ_Fonc.
 *
 */
int Sources::a_pour_Champ_Fonc(const Motcle& mot,
                               OBS_PTR(Champ_base)& ch_ref) const
{
  int ok = 0;
  for (const auto& itr : *this)
    {
      const Source& src = itr;
      const Source_base& src_base = src.valeur();
      if (src_base.a_pour_Champ_Fonc(mot, ch_ref))
        {
          ok = 1;
          break;
        }
    }
  return ok;
}

/*! @brief Calls Source::impr() on all sources in the list.
 *
 * See Source_base::impr().
 *
 */
int Sources::impr(Sortie& os) const
{
  for (const auto& itr : *this) itr->impr(os);
  return 1;
}

/*! @brief Sizing of the implicit matrix of source terms.
 *
 * Traverses all sources in the list to size the matrix.
 *
 */
void Sources::dimensionner(Matrice_Morse& matrice) const
{
  for (const auto& itr : *this)
    {
      const Source& src = itr;
      const Source_base& src_base = src.valeur();
      Matrice_Morse mat;
      src_base.dimensionner(mat);
      if (mat.nb_colonnes()) matrice += mat;
    }
}
/*! @brief Contribution to the implicit matrix of source terms. By default, no contribution.
 *
 */
void Sources::contribuer_a_avec(const DoubleTab& a, Matrice_Morse& matrice) const
{
  for (const auto& itr : *this)
    {
      const Source& src = itr;
      const Source_base& src_base = src.valeur();
      src_base.contribuer_a_avec(a,matrice);
    }
}

void Sources::contribuer_jacobienne(Matrice_Bloc& matrice, int n) const
{
  for (const auto& itr : *this)
    {
      const Source& src = itr;
      const Source_base& src_base = src.valeur();
      src_base.contribuer_jacobienne(matrice, n);
    }
}

/*! @brief Calls Source::initialiser(temps) on all sources in the list.
 *
 * See Source_base::initialiser(double temps).
 *
 */
int Sources::initialiser(double temps)
{
  int ok=1;
  for (auto& itr : *this)
    ok = ok && itr->initialiser(temps);
  return ok;
}

void Sources::check_multiphase_compatibility() const
{
  for (const auto& itr : *this)
    {
      const Source& src = itr;
      const Source_base& src_base = src.valeur();
      src_base.check_multiphase_compatibility();
    }
}
