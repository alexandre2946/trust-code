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

#include <Schema_Euler_Implicite.h>
#include <Probleme_Couple.h>
#include <EcrFicCollecte.h>
#include <Postraitement.h>
#include <Equation_base.h>
#include <Milieu_base.h>
#include <TRUSTTabs.h>
#include <Domaine_VF.h>
#include <TRUST_2_PDI.h>
#include <Ecrire_YAML.h>

Implemente_instanciable(Probleme_Couple,"Probleme_Couple",Couplage_U);
Implemente_instanciable(Probleme_Couple_Point_Fixe,"Probleme_Couple_Point_Fixe",Probleme_Couple);
// XD coupled_problem pb_gen_base probleme_couple INHERITS_BRACE This instruction causes a probleme_couple type object
// XD_CONT to be created. This type of object has an associated problem list, that is, the coupling of n problems among
// XD_CONT them may be processed. Coupling between these problems is carried out explicitly via conditions at particular
// XD_CONT contact limits. Each problem may be associated either with the Associate keyword or with the Read/groupes
// XD_CONT keywords. The difference is that in the first case, the four problems exchange values then calculate their
// XD_CONT timestep, rather in the second case, the same strategy is used for all the problems listed inside one group,
// XD_CONT but the second group of problem exchange values with the first group of problems after the first group did
// XD_CONT its timestep. So, the first case may then also be written like this: NL2 Probleme_Couple pbc NL2 Read pbc {
// XD_CONT groupes { { pb1 , pb2 , pb3 , pb4 } } } NL2 There is a physical environment per problem (however, the same
// XD_CONT physical environment could be common to several problems). NL2 Each problem is resolved in a domain. NL2
// XD_CONT Warning : Presently, coupling requires coincident meshes. In case of non-coincident meshes, boundary
// XD_CONT condition \'paroi_contact\' in VEF returns error message (see paroi_contact for correcting procedure).
// XD attr groupes list_list_nom groupes OPT { groupes { { pb1 , pb2 } , { pb3 , pb4 } } }
// XD ref domaine_2 domaine
// XD ref pb_1 Pb_base
// XD ref pb_2 Pb_base
// XD ref pb_3 Pb_base
// XD ref pb_4 Pb_base
// XD ref scheme_2 schema_temps_base

///////////////////////////////////////////////
//                                           //
// Implementation of the Problem interface  //
//                                           //
///////////////////////////////////////////////

bool Probleme_Couple::initTimeStep(double dt)
{
  /*
  double& residu_max = schema_temps().residu();
  double test=residu_max;
  for (int i=1; i<nb_problemes(); i++)
    residu_max = std::max(residu_max,sch_clones[i]->residu());
  if (test!=residu_max) abort();
  */
  bool ok =  Couplage_U::initTimeStep(dt);
  return ok;
}

double Probleme_Couple::computeTimeStep(bool& stop) const
{
  double dt_=DMAXFLOAT;

  // We take the minimum dt.
  // We stop if one of the problems wants to stop.
  for (int i=0; i<nb_problemes(); i++)
    {
      double dt1=probleme(i).computeTimeStep(stop);
      if (stop)
        return 0;
      if (dt1<dt_)
        dt_=dt1;
    }
  return dt_;
}

bool Probleme_Couple::solveTimeStep()
{
  // Trigger domain-specific time step logic at most once per distinct domain,
  // even if several problems share it.
  std::set<const Domaine*> processed_domains;
  for (int i = 0; i < nb_problemes(); i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base, probleme(i));
      Domaine& dom = pb.domaine();
      if (processed_domains.insert(&dom).second)
        dom.mettre_a_jour(schema_temps().temps_courant(), pb.domaine_dis(), pb);
    }

  // WEC: To be changed!!!!
  if (sch_clones.size())
    if (sub_type(Schema_Euler_Implicite,schema_temps()))
      {
        Schema_Euler_Implicite& sch_eul_imp=ref_cast(Schema_Euler_Implicite,schema_temps());
        int ok = 1, cv;
        cv = sch_eul_imp.faire_un_pas_de_temps_pb_couple(*this, ok);
        if (!ok) return false;
        // we propagate a certain number of things to the clones
        for (int i=1; i<sch_clones.size(); i++)
          {
            sch_clones[i]->facteur_securite_pas()=schema_temps().facteur_securite_pas();
            sch_clones[i]->set_stationnaire_atteint()=schema_temps().stationnaire_atteint();
            sch_clones[i]->residu()=schema_temps().residu();
          }
        return cv;
      }

  bool ok=Couplage_U::solveTimeStep();
  double& residu_max = schema_temps().residu();
  for (int i=1; i<nb_problemes(); i++)
    residu_max = std::max(residu_max,ref_cast(Probleme_base,probleme(i)).schema_temps().residu());
  return ok;
}

bool Probleme_Couple_Point_Fixe::solveTimeStep()
{
  std::vector<Schema_Euler_Implicite*> per_pb_schemas;
  per_pb_schemas.reserve(nb_problemes());
  for (int i = 0; i < nb_problemes(); i++)
    {
      Schema_Temps_base& sch_i = ref_cast(Probleme_base, probleme(i)).schema_temps();
      if (!sub_type(Schema_Euler_Implicite, sch_i))
        {
          per_pb_schemas.clear();
          break;
        }
      per_pb_schemas.push_back(&ref_cast(Schema_Euler_Implicite, sch_i));
    }

  if (int(per_pb_schemas.size()) != nb_problemes()) Process::exit("Probleme_Couple_Point_Fixe::solveTimeStep only works with problems using Schema_Euler_Implicite");

  for (int i = 0; i < nb_problemes(); i++)
    per_pb_schemas[i]->Initialiser_Champs(ref_cast(Probleme_base, probleme(i)));

  const int max_fp_iter = 1000;
  bool converged = false;
  int compteur = 0;
  int ok = 1;
  int nb_iter_min = 0; // minimum number of iterations before testing convergence

  while (compteur < nb_iter_min || (!converged && ok && compteur < max_fp_iter))
    {
      compteur++;
      converged = true;

      // 2. maillage
      std::set<const Domaine*> processed_domains;
      for (int i = 0; i < nb_problemes(); i++)
        {
          Probleme_base& pb = ref_cast(Probleme_base, probleme(i));
          Domaine& dom = pb.domaine();
          if (processed_domains.insert(&dom).second)
            dom.mettre_a_jour(schema_temps().temps_courant(), pb.domaine_dis(), pb);
        }

      for (int i = 0; i < nb_problemes(); i++)
        {
          Probleme_base& pb = ref_cast(Probleme_base, probleme(i));
          const double temps = per_pb_schemas[i]->temps_courant() + per_pb_schemas[i]->pas_de_temps();
          for (int eq = 0; eq < pb.nombre_d_equations(); eq++)
            pb.equation(eq).domaine_Cl_dis().calculer_coeffs_echange(temps);
        }

      // 1. problems
      for (int i = 0; ok && i < nb_problemes(); i++)
        {
          Probleme_base& pb = ref_cast(Probleme_base, probleme(i));
          pb.updateGivenFields();
          const int cv_pb = per_pb_schemas[i]->Iterer_Pb(pb, compteur, ok);
          converged = converged && cv_pb;
        }
    }

  if (!ok || !converged)
    {
      for (int i = 0; i < nb_problemes(); i++)
        per_pb_schemas[i]->notify_failed_timestep();
      if (limpr())
        Cerr << le_nom() << " : Failure in Probleme_Couple_Point_Fixe::solveTimeStep after " << compteur << " iterations." << finl;
      return false;
    }

  Cout << "Fixed-point convergence at t = " << schema_temps().temps_courant() << " in " << compteur << " iterations." << finl;
  for (int i = 0; i < nb_problemes(); i++)
    per_pb_schemas[i]->test_stationnaire(ref_cast(Probleme_base, probleme(i)));

  double& residu_max = schema_temps().residu();
  for (int i = 1; i < nb_problemes(); i++)
    residu_max = std::max(residu_max, ref_cast(Probleme_base, probleme(i)).schema_temps().residu());

  return ok;
}

bool Probleme_Couple::iterateTimeStep(bool& converged)
{
  bool ok=true;
  converged=true;

  int debut_gr=0;
  int fin_gr=0;
  int gr=0;
  while(ok && fin_gr<nb_problemes())
    {

      if (gr<groupes.size_array())
        fin_gr += groupes[gr++];
      else
        fin_gr=nb_problemes();

      for(int i=debut_gr; ok && i<fin_gr; i++)
        ok &= probleme(i).updateGivenFields();

      for(int i=debut_gr; ok && i<fin_gr; i++)
        {
          bool cv;
          ok &= probleme(i).iterateTimeStep(cv);
          converged = converged && cv;
        }

      debut_gr=fin_gr;
    }

  return ok;
}

////////////////////////////////////////////////////////
//                                                    //
// End of the implementation of the Problem interface //
//                                                    //
////////////////////////////////////////////////////////

Entree& Probleme_Couple::readOn(Entree& is)
{

  Cerr << "Reading of Probleme_Couple " << le_nom() << finl;

  Motcle motlu;
  is >> motlu;
  if (motlu != Motcle("{"))
    {
      Cerr << "We expected { to start to read the Probleme_Couple" << finl;
      exit();
    }

  is >> motlu;
  while (motlu!=Motcle("}"))   // end of readOn
    {

      if (motlu != Motcle("groupes"))
        {
          Cerr << "The keyword " << motlu << " is not understood" << finl;
          exit();
        }

      if (nb_problemes())
        {
          Cerr << "We can associate problems to Probleme_Couple" << finl;
          Cerr << "* either by \"associer prob_couple pb\" (in which case they are all in the same group)" << finl;
          Cerr << "* either by the keyword \"groupes\" while reading the object Probleme_Couple" << finl;
          Cerr << "but not both!" << finl;
        }
      assert(nb_problemes()==0);

      LIST(LIST(Nom)) les_noms;
      is >> les_noms;

      groupes.resize_array(les_noms.size());
      for (int i=0; i<les_noms.size(); i++)
        {
          groupes[i]=les_noms[i].size();
          for (int j=0; j<les_noms[i].size(); j++)
            {
              Nom nom_pb=les_noms[i][j];
              Objet_U& ob=Interprete::objet(nom_pb);
              Probleme_base& pb=ref_cast(Probleme_base,ob);
              ajouter(pb);
            }
        }

      is >> motlu;
    }

  return is;
}

/*! @brief Overrides Objet_U::printOn(Sortie&): prints the coupled problems to an output stream.
 *
 * @param (Sortie& os) the output stream to print to
 * @return (Sortie&) the modified output stream
 */
Sortie& Probleme_Couple::printOn(Sortie& os) const
{
  for(int i=1; i< nb_problemes(); i++)
    os << probleme(i) << finl;

  return os;
}

Entree& Probleme_Couple_Point_Fixe::readOn(Entree& is) { return Probleme_Couple::readOn(is); }
Sortie& Probleme_Couple_Point_Fixe::printOn(Sortie& os) const { return Probleme_Couple::printOn(os); }

bool Probleme_Couple::updateGivenFields()
{
  // No fields coming from outside the coupling.
  // Internal exchanges are done during iterateTimeStep.
  return true;
}

/*! @brief Adds a problem to the list of coupled problems.
 *
 * Sets up the reference from the problem back to this.
 *     Verifies that the conduction, thHyd order is respected.
 *
 * @param (Probleme_base& pb) the problem to add to the coupling
 */
void Probleme_Couple::ajouter(Probleme_base& pb)
{
  addProblem(pb);
  pb.associer_pb_couple(*this);

  int nb_pb=nb_problemes();


  Nom nom_pb=probleme(nb_pb-1).que_suis_je();
  if (nb_pb!=1
      && probleme(0).que_suis_je()=="Pb_Conduction"
      && nom_pb.prefix("Turbulent") != probleme(nb_pb-1).que_suis_je())
    {
      Cerr << finl;
      Cerr << "Warning !" << finl;
      Cerr << "You define a conduction problem named "<<probleme(0).le_nom()<<" before your turbulent thermalhydraulic problem named "<<probleme(nb_pb-1).le_nom()<<"." << finl;
      Cerr << "Please, define first your turbulent thermalhydraulic problem in your data file like:" << finl;
      Cerr << finl;
      Cerr << "Probleme_Couple "<<le_nom()<< finl;
      for (int i=nb_pb-1; i>=0; i--)
        Cerr << "Associer "<<le_nom()<<" "<<probleme(i).le_nom()<<finl;
      Cerr << finl;
      exit();
    }
}
void Probleme_Couple::initialize()
{
  // Assign value 1 to schema_impr_ for the time scheme of problem 0 and 0 for the others. Only one scheme should print.
  for (int i = 0; i < nb_problemes(); i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base, probleme(i));
      pb.schema_temps().schema_impr() = (i == 0);
    }

  Couplage_U::initialize();
}

/*! @brief Overrides Objet_U::associer_(Objet_U&): associates an object with the coupled problem, checking the type dynamically.
 *
 *     The object can be:
 *       - a time scheme (Schema_Temps_base), associated with the problems
 *       - a problem (Probleme_base), added to the list
 *
 * @param (Objet_U& ob) the object to associate
 * @return (int) 1 if the association succeeded, 0 otherwise
 * @throws if the object is not of an expected type
 */
int Probleme_Couple::associer_(Objet_U& ob)
{
  if( sub_type(Schema_Temps_base, ob))
    {
      associer_sch_tps_base(ref_cast(Schema_Temps_base, ob));
      return 1;
    }
  else if( sub_type(Probleme_base, ob))
    {
      Probleme_base& pb = ref_cast(Probleme_base, ob);
      ajouter(pb);
      return 1;
    }
  else
    return 0;
}

/*! @brief Associates a copy of the time scheme with each problem of the coupled problem.
 *
 * @param (Schema_Temps_base& sch) the time scheme to associate
 */
void Probleme_Couple::associer_sch_tps_base(Schema_Temps_base& sch)
{
  sch_clones.dimensionner(nb_problemes());
  for (int i=0; i<nb_problemes(); i++)
    {
      sch_clones[i]=sch; // Clone the scheme
      Probleme_base& pb=ref_cast(Probleme_base,probleme(i));
      pb.associer_sch_tps_base(sch_clones[i].valeur()); // association
      //We assign the value 1 to schema_impr_ for problem 0's scheme
      //and 0 for the others. Only one scheme should print.
      if (i!=0) pb.schema_temps().schema_impr()=0;
    }
}
/*! @brief Returns the time scheme associated with the coupled problems (const version).
 *
 * @return (Schema_Temps_base&) the associated time scheme
 */
const Schema_Temps_base& Probleme_Couple::schema_temps() const
{
  if (nb_problemes()==0)
    {
      Cerr << "You forgot to associate problems to the coupled problem named " << le_nom() << finl;
      Process::exit();
    }
  const Probleme_base& pb=ref_cast(Probleme_base,probleme(0));
  return pb.schema_temps();
}

/*! @brief Returns the time scheme associated with the coupled problems.
 *
 * @return (Schema_Temps_base&) the associated time scheme
 */
Schema_Temps_base& Probleme_Couple::schema_temps()
{
  if (nb_problemes()==0)
    {
      Cerr << "You forgot to associate problems to the coupled problem named " << le_nom() << finl;
      Process::exit();
    }
  Probleme_base& pb=ref_cast(Probleme_base,probleme(0));
  return pb.schema_temps();
}

/*! @brief Associates a discretization with all problems of the coupled problem.
 *
 *     Calls Probleme_Base::discretiser(const Discretisation_base&)
 *     on each of the sub-problems of the coupled problem.
 *     see Probleme_Base::discretiser(const Discretisation_base&)
 *
 * @param (Discretisation_base& dis) a discretization for all problems
 */
void Probleme_Couple::discretiser(Discretisation_base& dis)
{
  // Loop to discretize each problem of the coupled problem:
  for(int i=0; i< nb_problemes(); i++)
    {
      Cerr<<"The problem that we are starting to discretize is "<<probleme(i).le_nom()<<finl;
      Probleme_base& pb=ref_cast(Probleme_base,probleme(i));
      pb.discretiser(dis);
    }
}

void Probleme_Couple::sauver() const
{
  Ecrire_YAML yaml_file;
  bool pdi_format = false;
  Nom yaml_fname;
  for (int i=0; i<nb_problemes(); i++)
    {
      const Probleme_base& pb=ref_cast(Probleme_base,probleme(i));
      const Nom& format = pb.checkpoint_format();
      if(Motcle(format) == "pdi")
        {
          if(i>0 && pb.yaml_filename() != yaml_fname)
            {
              Cerr << "Probleme_Couple::sauver() Error! You have provided different yaml files for each of your problems to initialize PDI. It has to be the same. " << finl;
              Process::exit();
            }
          yaml_fname = pb.yaml_filename();
          const Nom& fname = pb.checkpoint_filename();
          yaml_file.add_pb_base(pb, fname);
          pdi_format = true;
        }
    }
  // we need to initialize PDI with a yaml file that contains the information of all the problems with a PDI checkpoint format
  // (allows to initialize PDI once for all checkpoints and not multiple times)
  if(pdi_format && !TRUST_2_PDI::is_PDI_initialized())
    {
      if(yaml_fname == "??")
        {
          yaml_fname = Nom("save_") + le_nom() + Nom(".yml");
          yaml_file.write_checkpoint_file(yaml_fname.getString());
        }
      TRUST_2_PDI::init(yaml_fname.getString());
    }

  for(int i=0; i<nb_problemes(); i++)
    ref_cast(Probleme_base,probleme(i)).sauver();

}
