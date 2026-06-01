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

#include <Cond_lim_utilisateur_base.h>
#include <Domaine_Cl_dis_base.h>
#include <Frontiere_dis_base.h>
#include <Schema_Temps_base.h>
#include <Equation_base.h>
#include <Probleme_base.h>

Implemente_base(Domaine_Cl_dis_base,"Domaine_Cl_dis_base",Objet_U);

Sortie& Domaine_Cl_dis_base::printOn(Sortie& os) const
{
  return os;
}

/*! @brief Overrides Objet_U::readOn(Sortie&) Reads the discretized boundary conditions from an input stream
 *
 *     The expected format is as follows:
 *     {
 *      Name Cond_lim [REPEAT READ AS MANY TIMES AS NECESSARY]
 *     }
 *
 * @param (Entree& is) an input stream
 * @return (Entree&) the modified input stream
 * @throws opening brace expected
 * @throws invalid number of boundary conditions read
 */
Entree& Domaine_Cl_dis_base::readOn(Entree& is)
{
  assert(mon_equation);
  const Domaine& ledomaine=equation().domaine_dis().domaine();
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Nom nomlu;
  Motcle motlu;
  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Error while reading the boundaries conditions\n";
      Cerr << "We expected a " << accolade_ouverte << " instead of \n"
           << motlu;
      exit();
    }

  int n = ledomaine.nb_front_Cl();
  IntTab front_deja_lu(n);
  front_deja_lu = 0;
  les_conditions_limites().dimensionner(n);
  int nb_clim=0;

  while(1)
    {

      // Reading a boundary name or }
      is >> nomlu;
      motlu=nomlu;
      if (motlu == accolade_fermee)
        break;

      Cerr << "Reading the " << nomlu << " boundary condition of the " << equation().que_suis_je() << " equation from the problem " << equation().probleme().le_nom() << finl;
      Journal()<< "Reading the boundary condition " << nomlu <<finl;

      int rang=ledomaine.rang_frontiere(nomlu);

      // Additional test on boundary conditions that have been read:
      // we test if two boundaries have the same name
      if (front_deja_lu(rang) == 0)
        front_deja_lu(rang) = 1;
      else
        {
          Cerr <<finl;
          Cerr <<"Error in the data set "<<finl;
          Cerr <<"the boundary condition associated"<<finl;
          Cerr <<"to the boundary "<< nomlu <<"is read twice !!"<<finl;
          exit();
        }
      //      const Frontiere& frontiere=domaine.frontiere(rang);
      is >> les_conditions_limites(rang);
      // If we had a condition_limite_utilisateur
      // We must retrieve the bc stored in it
      // copy it and destroy the one contained in the cond_utilisateur
      if (sub_type(Cond_lim_utilisateur_base,les_conditions_limites(rang).valeur()))
        {
          Cond_lim_utilisateur_base& la_cl=ref_cast(Cond_lim_utilisateur_base,les_conditions_limites(rang).valeur());
          la_cl.lire(is,equation(),nomlu);
          Cond_lim* sa=&(la_cl.la_cl());

          // WEC: The new formulation avoids copying the Cond_lim_base object.
          // This copy can cause problems for fields that register themselves in a vector
          // (Champ_Inputs). The vector would then point to garbage...
          // Adopt the new one and destroy the old one (Cond_lim_utilisateur)
          les_conditions_limites(rang).adopt(*sa);
          //les_conditions_limites(rang)=(*sa);
          delete sa;
        }
      les_conditions_limites(rang)->associer_fr_dis_base(domaine_dis().frontiere_dis(rang));

      // Test to prevent the use of 'Raccord_distant_homogene' in sequential computation
      const Frontiere& frontiere=ledomaine.frontiere(rang);
      if ((frontiere.que_suis_je()=="Raccord_distant_homogene") && Process::is_sequential())
        {
          Cerr<<"At least one connection (adjacent boundary on two domains) is of type 'Raccord distant homogene'." << finl;
          Cerr<<"Use 'Raccord local homogene' to define the connections in sequential computing"<<finl;
          Cerr<<"And 'Raccord distant homogene' for the connections in parallel computing." << finl;
          exit();
        }

      nb_clim++;
    }
  if (nb_clim!=n)
    {
      domaine_dis().ecrire_noms_bords(Cerr);
      Cerr << "It misses " << n-nb_clim << " boundaries conditions " << finl;
      Cerr << "We read " << nb_clim << " boundaries conditions " << finl;
      Cerr << "We waited " << n << " boundary conditions " << finl;
      exit();
    }
  for (int num_cl=0; num_cl<nb_clim; num_cl++)
    {
      les_conditions_limites(num_cl)->associer_domaine_cl_dis_base(*this);
      les_conditions_limites(num_cl)->verifie_ch_init_nb_comp();
      //      const Frontiere_dis_base& la_fr_dis = les_conditions_limites(num_cl).frontiere_dis();
      //      const Frontiere& frontiere=la_fr_dis.frontiere();
      //      const Domaine& domaine2=frontiere.domaine();
    }

  return is;
}

/*! @brief Returns 1 if the object contains a boundary condition with the specified Name.
 *
 *     Returns 0 otherwise.
 *
 * @param (Nom& type) the name of the boundary condition to search for
 * @return (int) 1 if the boundary condition with the specified name was found, 0 otherwise.
 */
int Domaine_Cl_dis_base::contient_Cl(const Nom& type)
{
  for (auto& itr : les_conditions_limites_)
    if (itr.get_info()->has_base(type)) return 1;
  return 0;
}

/*! @brief Returns a reference to the discretized domain associated with the boundary conditions.
 *
 * This Domaine_dis is associated through the associated equation
 *     and not directly with the Domaine_Cl_dis_base object.
 *
 * @return (Domaine_dis_base&) the discretized domain associated with the equation associated with the boundary conditions.
 */
Domaine_dis_base& Domaine_Cl_dis_base::domaine_dis()
{
  return equation().domaine_dis();
}

/*! @brief Returns a reference to the discretized domain associated with the boundary conditions.
 *
 * This Domaine_dis is associated through the associated equation
 *     and not directly with the Domaine_Cl_dis_base object.
 *     (const version)
 *
 * @return (Domaine_dis_base&) the discretized domain associated with the equation associated with the boundary conditions.
 */
const Domaine_dis_base& Domaine_Cl_dis_base::domaine_dis() const
{
  return equation().domaine_dis();
}

/*! @brief Changes the i-th future time of all BCs.
 *
 */
void Domaine_Cl_dis_base::changer_temps_futur(double temps,int i)
{
  for (int j=0; j<nb_cond_lim(); j++)
    les_conditions_limites_[j]->changer_temps_futur(temps,i);
}

/*! @brief Changes the i-th future time of all BCs.
 *
 */
void Domaine_Cl_dis_base::set_temps_defaut(double temps)
{
  for (int j=0; j<nb_cond_lim(); j++)
    les_conditions_limites_[j]->set_temps_defaut(temps);
}

/*! @brief Rotates the wheel of all BCs up to the given time
 *
 */
int Domaine_Cl_dis_base::avancer(double temps)
{
  int ok=1;
  for (int j=0; j<nb_cond_lim(); j++)
    ok = ok && les_conditions_limites_[j]->avancer(temps);
  return ok;
}

/*! @brief Rotates the wheel of all BCs back to the given time
 *
 */
int Domaine_Cl_dis_base::reculer(double temps)
{
  int ok=1;
  for (int j=0; j<nb_cond_lim(); j++)
    ok = ok && les_conditions_limites_[j]->reculer(temps);
  return ok;
}

/*! @brief Performs a time update of all boundary conditions.
 *
 * @param (double temps) the time step for update
 */
void Domaine_Cl_dis_base::mettre_a_jour(double temps)
{
  les_conditions_limites_.mettre_a_jour(temps);
}

/* @brief See ICoCo::ProblemTrio::resetTime()
 */
void Domaine_Cl_dis_base::resetTime(double temps)
{
  les_conditions_limites_.resetTime(temps);
}

/*! @brief Performs a time update for sub-time-steps of a time scheme (e.g. in RungeKutta)
 *
 *     for all boundary conditions returning 1 through the method
 *     int Cond_Lim_base::a_mettre_a_jour_ss_pas_dt();
 *
 * @param (double temps) the time step for update
 */
void Domaine_Cl_dis_base::mettre_a_jour_ss_pas_dt(double temps)
{
  for (auto &itr : les_conditions_limites_)
    {
      Cond_lim_base& la_cl = itr.valeur();
      if (la_cl.a_mettre_a_jour_ss_pas_dt() == 1)
        la_cl.mettre_a_jour(temps);
    }
}

/*! @brief Initializes the BCs. Unlike the update methods, the
 *
 *     initialize methods of BCs cannot depend on the outside
 *     (it may not be initialized itself)
 *
 * @return (int) 1 if OK, 0 otherwise
 */
int Domaine_Cl_dis_base::initialiser(double temps)
{
  return les_conditions_limites_.initialiser(temps);
}

/*! @brief Computes the exchange coefficients for thermally coupled problems
 *
 * @return (int) always returns 1
 */
int Domaine_Cl_dis_base::calculer_coeffs_echange(double temps)
{
  les_conditions_limites_.calculer_coeffs_echange(temps);
  return 1;
}

/*! @brief Calls Cond_lim_base::completer() on each boundary condition
 *
 */
void Domaine_Cl_dis_base::completer()
{
  les_conditions_limites_.completer(domaine_dis());
  completer(domaine_dis());
}

/*! @brief Returns the boundary condition associated with a given real face.
 *
 * Stores in face_locale the face number within the boundary.
 *  Triggers an error if the face does not have a BC.
 *
 */
const Cond_lim_base& Domaine_Cl_dis_base::condition_limite_de_la_face_reelle(int face_globale, int& face_locale) const
{
  for (int i=0; i<nb_cond_lim(); i++)
    {
      const Frontiere& fr=les_conditions_limites(i)->frontiere_dis().frontiere();
      if (face_globale>=fr.num_premiere_face() && face_globale < fr.num_premiere_face()+fr.nb_faces())
        {
          face_locale=face_globale-fr.num_premiere_face();
          return les_conditions_limites(i).valeur();
        }
    }
  assert(0); // the face does not have a BC
  return les_conditions_limites(0).valeur(); // For compilation
}

/*! @brief Returns the boundary condition associated with a given virtual face.
 *
 * Stores in face_locale the face number within the boundary.
 *  Triggers an error if the face does not have a BC.
 *
 */
const Cond_lim_base& Domaine_Cl_dis_base::condition_limite_de_la_face_virtuelle(int face_globale, int& face_locale) const
{
  for (int i=0; i<nb_cond_lim(); i++)
    {
      const Frontiere& fr=les_conditions_limites(i)->frontiere_dis().frontiere();
      const ArrOfInt& faces_virt=fr.get_faces_virt();
      for (int j=0; j<faces_virt.size_array(); j++)
        if (face_globale==faces_virt[j])
          {
            face_locale=fr.nb_faces()+j;
            return les_conditions_limites(i).valeur();
          }
    }
  assert(0); // the face does not have a BC
  return les_conditions_limites(0).valeur(); // For compilation
}

/*! @brief Returns the boundary condition associated with a boundary of the given name.
 *
 * Triggers an error if no boundary has this name.
 *
 */
Cond_lim_base& Domaine_Cl_dis_base::condition_limite_de_la_frontiere(Nom frontiere)
{
  for (int i=0; i<nb_cond_lim(); i++)
    {
      const Frontiere& fr=les_conditions_limites(i)->frontiere_dis().frontiere();
      if (fr.le_nom()==frontiere)
        return les_conditions_limites(i).valeur();
    }
  assert(0); // No boundary with this name
  exit();
  return les_conditions_limites(0).valeur(); // For compilation
}

/*! @brief Returns the boundary condition associated with a boundary of the given name.
 *
 * Triggers an error if no boundary has this name.
 *
 */
const Cond_lim_base& Domaine_Cl_dis_base::condition_limite_de_la_frontiere(Nom frontiere) const
{
  for (int i=0; i<nb_cond_lim(); i++)
    {
      const Frontiere& fr=les_conditions_limites(i)->frontiere_dis().frontiere();
      if (fr.le_nom()==frontiere)
        return les_conditions_limites(i).valeur();
    }
  assert(0); // No boundary with this name
  exit();
  return les_conditions_limites(0).valeur(); // For compilation
}

/*! @brief Computes the growth rate of unsteady BCs between t1 and t2.
 *
 */
void Domaine_Cl_dis_base::calculer_derivee_en_temps(double t1, double t2)
{
  for (int i=0; i<nb_cond_lim(); i++)
    {
      Champ_front_base& champ=les_conditions_limites(i)->champ_front();
      if (champ.instationnaire()) champ.calculer_derivee_en_temps(t1,t2);
    }
}

/*! @brief Returns the i-th boundary condition.
 *
 * (const version)
 *
 * @param (int i) the rank of the i-th boundary condition
 * @return (Cond_lim&) the i-th boundary condition
 */
const Cond_lim& Domaine_Cl_dis_base::les_conditions_limites(int i) const
{
  return les_conditions_limites_[i];
}

/*! @brief Returns the i-th boundary condition.
 *
 * @param (int i) the rank of the i-th boundary condition
 * @return (Cond_lim&) the i-th boundary condition
 */
Cond_lim& Domaine_Cl_dis_base::les_conditions_limites(int i)
{
  return les_conditions_limites_[i];
}

/*! @brief Returns the array of boundary conditions.
 *
 * @return (Conds_lim&) the array of boundary conditions
 */
Conds_lim& Domaine_Cl_dis_base::les_conditions_limites()
{
  return les_conditions_limites_;
}

/*! @brief Returns the array of boundary conditions.
 *
 * (const version)
 *
 * @return (Conds_lim&) the array of boundary conditions
 */
const Conds_lim& Domaine_Cl_dis_base::les_conditions_limites() const
{
  return les_conditions_limites_;
}

/*! @brief Returns the number of boundary conditions.
 *
 * @return (int) the number of boundary conditions
 */
int Domaine_Cl_dis_base::nb_cond_lim() const
{
  return les_conditions_limites_.size();
}

int  Domaine_Cl_dis_base::nb_faces_Cl() const
{
  return domaine().nb_faces_frontiere();
}

Domaine& Domaine_Cl_dis_base::domaine()
{
  return domaine_dis().domaine();
}
const Domaine& Domaine_Cl_dis_base::domaine() const
{
  return domaine_dis().domaine();
}

/*! @brief Given a boundary face index in the Domaine_VF, returns the boundary condition to which this face
 *
 *   belongs, for 0 <= num_face < nb_faces_Cl().
 *
 */
const Cond_lim& Domaine_Cl_dis_base::la_cl_de_la_face(int num_face) const
{
  // Generic algorithm: we iterate over the boundaries until we find
  // the one that contains the face.
  // The faces of boundary 0, then those of boundary 1, etc., follow
  // each other in the Domaine_VF.
  assert(num_face >= 0);
  int i = 0;
  const int nb_cl = les_conditions_limites_.size();
  for (i = 0; i < nb_cl; i++)
    {
      const Cond_lim_base& cl = les_conditions_limites_[i].valeur();
      const Frontiere&      fr = cl.frontiere_dis().frontiere();
      int num_premiere_face = fr.num_premiere_face();
      int nb_faces          = fr.nb_faces();

      if (num_face >= num_premiere_face && num_face < num_premiere_face + nb_faces)
        break;
    }
  if (i == nb_cl)
    {
      Cerr << "Error in Domaine_Cl_dis_base::la_cl_de_la_face(num_face="
           << num_face << ")\n This face is not on a boundary of cond_lim."
           << finl;
      exit();
    }
  return les_conditions_limites_[i];
}

void Domaine_Cl_dis_base::nommer(const Nom& un_nom)
{
  nom_ = un_nom;
}

void Domaine_Cl_dis_base::associer_inconnue(const Champ_Inc_base& inco)
{
  mon_inconnue=inco;
}

const Champ_Inc_base& Domaine_Cl_dis_base::inconnue() const
{
  return mon_inconnue;
}

Champ_Inc_base& Domaine_Cl_dis_base::inconnue()
{
  return mon_inconnue;
}
