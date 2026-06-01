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

#include <Schema_Temps_base.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Process.h>
#include <Source_base.h>
#include <TRUSTTrav.h>
#include <TRUSTTabs.h>
#include <SFichier.h>
#include <map>

Implemente_base(Source_base,"Source_base",Objet_U);
// XD source_base objet_u source_base INHERITS_BRACE Basic class of source terms introduced in the equation.

/*! @brief DOES NOTHING - to override in derived classes.
 *
 *     Prints the source to an output stream.
 *
 * @param (Sortie& os) the output stream
 * @return (Sortie&) the modified output stream
 */
Sortie& Source_base::printOn(Sortie& os) const
{
  return os;
}

/*! @brief DOES NOTHING - to override in derived classes.
 *
 *     Reading of a source term from an input stream.
 *
 * @param (Entree& is) the input stream
 * @return (Entree&) the modified input stream
 */
Entree& Source_base::readOn(Entree& is)
{
  return is;
}

/*! @brief DOES NOTHING - to override in derived classes.
 *
 * Time update of the source term.
 * @param (double) the time step of the update
 */
void Source_base::mettre_a_jour(double )
{
  Cerr << finl;
  Cerr << "You must to overload the method Source_base::mettre_a_jour() because there's a big" << finl;
  Cerr << "chance that the source term " <<  que_suis_je() << " must be updated," << finl;
  Cerr << "especially if your source has a field, you must update it..." << finl << finl;
  exit();
}

/*! @brief DOES NOTHING - to override in derived classes.
 *
 * Reset current time.
 * @param (double) new current time to be set.
 */
void Source_base::resetTime(double t)
{
  Cerr << finl;
  Cerr << "You must to overload the method Source_base::resetTime() because there's big" << finl;
  Cerr << "chance that the source term " <<  que_suis_je() << " must be reset," << finl;
  Cerr << "especially if your source has a field, you must reset its current time..." << finl << finl;
  exit();
}

/*! @brief Updates internal references of the Source_base object.
 *
 * Calls 2 protected pure virtual methods:
 *        Source_base::associer_domaines(const Domaine_dis_base& ,const Domaine_Cl_dis_base&)
 *        Source_base::associer_pb(const Probleme_base&)
 *
 */
void Source_base::completer()
{
  const Equation_base& eqn = equation();
  const Domaine_dis_base& zdis= eqn.domaine_dis();
  const Domaine_Cl_dis_base& zcldis = eqn.domaine_Cl_dis();
  associer_domaines(zdis, zcldis);
  associer_pb(eqn.probleme());
  // Initialize the bilan_ array:
  bilan_.resize(eqn.inconnue().nb_comp());
  bilan_=0;
  for (auto& itr : col_names_)
    col_width_ = std::max(col_width_, itr.longueur());
}

/*! @brief This method (or the method of the derived class) is called by Sources::associer_champ_rho for each source in the list
 *
 *   (for example, at the initialization of a front-tracking calculation).
 *   This method must be re-implemented in the derived classes
 *   used in problems with variable rho.
 *
 *   The ajouter method calculates the following term:
 *     INTEGRAL            (source term)
 *     over interleaved volume
 *
 *   In problems where rho is variable, "source term" homogeneous to rho*v.
 *   Otherwise, "source term" is homogeneous to v.
 *
 */
void Source_base::associer_champ_rho(const Champ_base& champ_rho)
{
  Cerr << "In Source_base::associer_champ_rho" << finl;
  Cerr << " field : " << champ_rho.le_nom() << finl;
  Cerr << " Source of type : " << que_suis_je() << finl;
  Cerr << " This source does not support the association of a field rho." << finl;
  Cerr << " (method associer_champ_rho must be coded for this source term)." << finl;
  assert(0);
  exit();
}

/*! @brief If the source understands the keyword "mot", it fills the reference to ch_ref and returns 1, otherwise returns 0 (see Source_Translation for example)
 *
 * )
 *
 */
int Source_base::a_pour_Champ_Fonc(const Motcle& mot, OBS_PTR(Champ_base) &ch_ref) const
{
  // The base class understands no keywords
  return 0;
}

const Champ_base& Source_base::get_champ(const Motcle& nom) const
{
  return champs_compris_.get_champ(nom);
}

bool Source_base::has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const
{
  return champs_compris_.has_champ(nom, ref_champ);
}

bool Source_base::has_champ(const Motcle& nom) const
{
  return champs_compris_.has_champ(nom);
}

void Source_base::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  if (opt == DESCRIPTION)
    Cerr << que_suis_je() << " : " << champs_compris_.liste_noms_compris() << finl;
  else
    nom.add(champs_compris_.liste_noms_compris());
}

/*! @brief DOES NOTHING - to override in derived classes.
 *
 *     Time update of the source term.
 *
 * @param (double) the time step of the update
 */
static std::map<std::string, int> counters; // Todo provisional
int Source_base::impr(Sortie& os) const
{
  if (out_=="??")
    {
      //Cerr << "No balance printed for " << que_suis_je() << finl;
    }
  else
    {
      int nb_compo=bilan_.size();
      if (nb_compo==0)
        {
          Cerr << "No balance printed for " << que_suis_je() << finl;
          Cerr << "cause bilan_ array is not filled." << finl;
        }
      else
        {
          int flag=je_suis_maitre();
          ouvrir_fichier(Flux,"",flag);
          const Probleme_base& pb=equation().probleme();
          const Schema_Temps_base& sch=pb.schema_temps();
          double temps=sch.temps_courant();

          if(Process::je_suis_maitre())
            Flux.add_col(temps);
          DoubleVect bilan_p(bilan_);
          mp_sum_for_each_item(bilan_p);
          /*
          for(int k=0; k<nb_compo; k++)
               bilan_(k)=Process::mp_sum(bilan_(k)); */
          // Replaced by:
          // mp_sum_for_each_item(bilan_); // Fixed bug: double line!

          if (Process::je_suis_maitre())
            {
              for(int k=0; k<nb_compo; k++)
                Flux.add_col(bilan_p(k));
              Flux << finl;
            }
        }
    }
  return 1;
}
/*! @brief Sizing of the implicit matrix of source terms.
 *
 * By default does nothing.
 *
 */
void Source_base::dimensionner(Matrice_Morse& mat) const
{
  if (has_interface_blocs()) dimensionner_blocs({{ equation().inconnue().le_nom().getString(), &mat }});
}

void Source_base::dimensionner_bloc_vitesse(Matrice_Morse& mat) const
{
  if (has_interface_blocs()) dimensionner_blocs({{ "vitesse", &mat }});
}

DoubleTab& Source_base::ajouter(DoubleTab& secmem) const
{
  if (has_interface_blocs()) ajouter_blocs({}, secmem, {});
  else Process::exit(que_suis_je() + " : ajouter() not coded!");
  return secmem;
}

DoubleTab& Source_base::calculer(DoubleTab& secmem) const
{
  secmem = 0;
  return ajouter(secmem);
}

/*! @brief contribution to the implicit matrix of source terms, by default no contribution
 *
 */
void Source_base::contribuer_a_avec(const DoubleTab&, Matrice_Morse& mat) const
{
  if (!has_interface_blocs()) return;
  DoubleTrav secmem(equation().inconnue().valeurs()); //sera jete
  ajouter_blocs({{ equation().inconnue().le_nom().getString(), &mat}}, secmem, {});
}

/*! @brief contribution to the right-hand side of source terms implicitly, by default error
 *
 *  method presented for consistency with Operateur_base
 *
 */
void Source_base::contribuer_au_second_membre(DoubleTab& ) const
{
  Cerr<<"Source_base::contribuer_au_second_membre(DoubleTab& ) const uncoded"<<finl;
  exit();
}

/* by default error */
void Source_base::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  Process::exit(que_suis_je() + " : dimensionner_blocs() not coded!");
}
void Source_base::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  Process::exit(que_suis_je() + " : ajouter_blocs() not coded!");
}

/*! @brief Unlike the mettre_a_jour methods, the initializer methods of sources cannot depend on the outside
 *
 *     (itself may not be initialized)
 *     By default, mettre_a_jour
 *
 */
int Source_base::initialiser(double temps)
{
  mettre_a_jour(temps);
  return 1;
}


/*! @brief Opening/creation of a file for printing a source term. To override in derived classes.
 *
 * @throws method to override
 */
void Source_base::ouvrir_fichier(SFichier& os,const Nom& type, const int flag) const
{
  // null flag means we don't open the file
  if (flag==0)
    return ;
  // Todo provisional:
  counters[type.getString()]++;
  if (counters[type.getString()]>1 && type!="")
    {
      Cerr << "Code should be rewritten to have only one call to Source_base::ouvrir_fichier for " << type << " source and not " << counters[type.getString()] << " times." << finl;
      Process::exit();
    }
  const Probleme_base& pb=equation().probleme();
  const Schema_Temps_base& sch=pb.schema_temps();
  const int precision = sch.precision_impr(), wcol = std::max(col_width_, sch.wcol()), gnuplot_header = sch.gnuplot_header();
  os.set_col_width(wcol);

  Nom nomfichier(out_);
  if (type!="") nomfichier+=(Nom)"_"+type;
  nomfichier+=".out";

  // Create the file on the first printing with the header
  if (sch.nb_impr()==1 && !pb.reprise_effectuee())
    {
      os.ouvrir(nomfichier);
      SFichier& fic=os;
      //Nom espace="\t\t";
      fic << (Nom)"# Printing of the source term "+que_suis_je()+" of the equation "+equation().que_suis_je()+" of the problem "+equation().probleme().le_nom() << finl;
      fic << "# " << description() << finl;
      if (!gnuplot_header) fic << "#";
      os.set_col_width(wcol - !gnuplot_header);
      fic.add_col("Time");
      os.set_col_width(wcol);
      for (int i = 0; i < col_names_.size(); i++)
        fic.add_col(col_names_[i]);
      fic << finl;
    }
  // Otherwise open it
  else
    {
      os.ouvrir(nomfichier,ios::app);
    }
  os.precision(precision);
  os.setf(ios::scientific);
}

void Source_base::set_fichier(const Nom& nom)
{
  out_=Objet_U::nom_du_cas();
  out_+="_";
  out_+=equation().probleme().le_nom()+"_"+nom;
}

