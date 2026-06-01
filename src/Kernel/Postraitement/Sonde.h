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

#ifndef Sonde_included
#define Sonde_included

#include <TRUSTArrays.h>
#include <TRUST_Ref.h>
#include <TRUSTTab.h>
#include <SFichier.h>
#include <Motcle.h>

class Operateur_Statistique_tps_base;
class Champ_Generique_base;
class Postraitement;
class Champ_base;
#include <Domaine_forward.h>

/*! @brief class Sonde. This class allows tracking the evolution of a field over time.
 *
 *      The set of points at which the field is to be probed and
 *      the observation period are chosen. Postraitement objects hold
 *      probes on the fields to observe; a probe also holds a
 *      reference to a post-processing object.
 *
 * @sa Postraitement Sondes
 */
class Sonde : public Objet_U
{
  Declare_instanciable_sans_constructeur_ni_destructeur(Sonde);

public :

  Sonde();
  Sonde(const Nom& );
  void associer_post(const Postraitement& );
  virtual void initialiser();
  virtual void mettre_a_jour(double temps, double tinit);
  virtual void postraiter();
  void ouvrir_fichier();
  virtual void completer();
  inline void fermer_fichier();
  inline const Champ_Generique_base& le_champ() const;
  inline const DoubleTab& les_positions_sondes_initiales() const; // Initial positions
  inline const DoubleTab& les_positions_sondes() const; // Positions after displacement
  inline const DoubleTab& les_positions() const; // Proc-local positions
  inline const IntVect& les_poly() const;
  inline void fixer_periode(double);
  inline double temps() const;
  inline SFichier& fichier();
  inline ~Sonde() override;
  inline const Nom& get_nom() const { return nom_; }
  inline const Nom& get_type() const { return type_; }
  inline const int& get_dim() const { return dim ; }
  inline void nommer(const Nom& n) override { nom_ = n; }

  virtual void ajouter_bords(const DoubleTab& coords_bords);
  virtual void init_bords();
  virtual void mettre_a_jour_bords();
  // Boundary processing (option "gravcl")
  void resetTime(double time) { nb_bip = time/periode; };

protected :
  /** Retrieve the domain to be used for the probe */
  virtual const Domaine& get_domaine_geom() const;
  /** Retrieve the field names - overriden in IJK */
  virtual const Noms get_noms_champ() const;
  /** Retrieve the number of component for the field - overriden in IJK */
  virtual int get_nb_compo_champ() const;
  /** Retrieve the current time for the field - overriden in IJK */
  virtual double get_temps_champ() const;
  /** Should exit if invalid type for a probe - overriden for IJK */
  virtual void validate_type(const Motcle& loc) const { }
  /** Should exit if invalid position has been set for a probe - overriden for IJK */
  virtual void validate_position() const { }
  virtual void create_champ_generique(Entree& is, const Motcle& motlu);
  /** Fix probe position - used by IJK */
  virtual void fix_probe_position() { }
  /** Fix probe position when keyword grav or gravcl is used */
  virtual void fix_probe_position_grav();
  /** Fill the array 'valeurs_locales' used by each proc to store prob local values */
  virtual void fill_local_values();
  /** Update the underlying field source if needed */
  virtual void update_source(double un_temps);

  OBS_PTR(Postraitement) mon_post;
  Nom nom_;                               ///< the probe name
  Nom nom_fichier_;                       ///< the name of the file containing the probe
  int dim;                                ///< the dimension of the probe (point:0,segment:1,plan:2,volume:3)
  OBS_PTR(Champ_Generique_base) mon_champ;
  OBS_PTR(Operateur_Statistique_tps_base) operateur_statistique_;        // Reference to an optional statistical operator
  /** Index of the component to probe. If ncomp = -1 the probe applies to all
   * components of the field */
  int ncomp;
  DoubleTab les_positions_sondes_initiales_;   ///< coordinates of the initial point probes
  DoubleTab les_positions_sondes_;             ///< coordinates of probes across the whole domain after displacement (master only)
  DoubleTab les_positions_;               ///< coordinates of probes local to each proc
  int numero_elem_;                       ///< equals -1 if undefined, otherwise the element number on the master
  IntVect elem_;                          ///< elements containing the local point probes
  double periode;                         ///< sampling period
  /** keys for probe typing (probes redefined at nodes or from values at vertices or at the center of gravity or at vertices)
   */
  bool nodes = false;
  bool chsom = false;
  bool grav = false;
  bool gravcl = false; // Values at centers of gravity (like grav) but with optional addition of boundary values via the CL domain of the post-processed field
  bool som = false;
  DoubleTab valeurs_locales,valeurs_sur_maitre;     ///< valeurs_locales: values on each proc, valeurs_sur_maitre: values gathered on the master
  double nb_bip;
  SFichier le_fichier_;
  Motcle nom_champ_lu_;
  ArrsOfInt participant ;            // on the master: participant[pe][i] -> the i-th point on pe corresponds to the participant[pe][i]-th position
  Nom type_;
  int orientation_faces_;

  // Boundary processing (option "gravcl")
  ArrOfInt faces_bords_;                  ///< array containing the indices of the boundary faces hit by the probe
  IntTab rang_cl_;                        ///< for a given face, index of the CL that this face bears
  int nbre_points1 = -1,nbre_points2 = -1,nbre_points3 = -1;        ///< used to create sonde_segment, sonde_plan, etc...
};


/*! @brief The elapsed time.
 *
 * @return (double) the elapsed time
 */
inline double Sonde::temps() const
{
  return nb_bip*periode;
}



/*! @brief Closes the file to which the probe writes.
 *
 */
inline void Sonde::fermer_fichier()
{
  if (fichier().is_open()) le_fichier_.close();
}


/*! @brief Sets the period at which the field is probed.
 *
 * @param (double pe) the probing period for the field
 */
inline void Sonde::fixer_periode(double pe)
{
  periode=pe;
}


/*! @brief Returns the associated field.
 *
 * @return (Champ_base&) the associated field
 */
inline const Champ_Generique_base& Sonde::le_champ() const
{
  return mon_champ.valeur();
}


/*! @brief Returns the array of positions of the field that are probed.
 *
 * @return (DoubleTab&) the probed positions
 */

inline const DoubleTab& Sonde::les_positions_sondes_initiales() const
{
  return les_positions_sondes_initiales_;
}

inline const DoubleTab& Sonde::les_positions_sondes() const
{
  return les_positions_sondes_;
}

inline const DoubleTab& Sonde::les_positions() const
{
  return les_positions_;
}

/*! @brief Returns the array of elements that are probed.
 *
 * @return (IntVect&) the elements that are probed
 */
inline const IntVect& Sonde::les_poly() const
{
  return elem_;
}


/*! @brief Returns an output file stream pointing to the output file used by the probe.
 *
 * @return (SFichier&) the output file used by the probe
 */
inline SFichier& Sonde::fichier()
{
  return le_fichier_;
}

/*! @brief Destructor.
 *
 * Closes the file before destroying the object.
 *
 */
inline Sonde::~Sonde()
{
  fermer_fichier();
}


#endif
