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

#ifndef Turbulence_paroi_base_included
#define Turbulence_paroi_base_included

#include <Champs_compris_interface.h>
#include <Champs_compris.h>
#include <Champ_Inc_base.h>
#include <TRUST_Ref.h>
#include <YAML_data.h>

class Modele_turbulence_hyd_base;
class Probleme_base;
class EcrFicPartage;
class Param;

/*! @brief Base class for the hierarchy of wall-law models computing turbulent quantities near walls.
 *         Contains a reference to a turbulence model.
 *
 * @sa Paroi_std_hyd_VDF Paroi_std_scal_hyd_VDF
 */
class Turbulence_paroi_base: public Champs_compris_interface, public Objet_U
{
  Declare_base(Turbulence_paroi_base);
public:

  static void typer_lire_turbulence_paroi(OWN_PTR(Turbulence_paroi_base)&, const Modele_turbulence_hyd_base&, Entree& );

  virtual void set_param(Param& param) const override { /* Do nothing */ }
  inline void associer_modele(const Modele_turbulence_hyd_base&);
  virtual void associer(const Domaine_dis_base&, const Domaine_Cl_dis_base&)=0;
  virtual void completer() { }
  virtual int init_lois_paroi() =0;
  inline int calculer_hyd(Champ_Inc_base&);
  inline int calculer_hyd(Champ_Fonc_base&, Champ_Fonc_base&);
  virtual int calculer_hyd(DoubleTab&) =0;
  virtual int calculer_hyd(DoubleTab&, DoubleTab&) =0;
  virtual int calculer_hyd_BiK(DoubleTab&, DoubleTab&) =0;
  inline virtual DoubleTab& corriger_derivee_impl(DoubleTab& d) const { return d; }
  inline virtual void imprimer_ustar(Sortie&) const { }
  virtual void imprimer_premiere_ligne_ustar(int, const LIST(Nom)&, const Nom&) const;
  virtual void imprimer_ustar_mean_only(Sortie&, int, const LIST(Nom)&, const Nom&) const;
  // added to account for Cisaillement_paroi in the base class

  inline const DoubleTab& Cisaillement_paroi() const;
  inline const DoubleVect& tab_u_star() const;
  inline double tab_u_star(int face) const;
  inline const DoubleVect& tab_d_plus() const;
  inline double tab_d_plus(int face) const;

  //OC 01/2006: added save/resume functionality: useful for TBLE for now.
  int sauvegarder(Sortie&) const override { return 0; }
  int reprendre(Entree&) override { return 0; }
  virtual std::vector<YAML_data> data_a_sauvegarder() const { return std::vector<YAML_data>(); }

  void creer_champ(const Motcle& motlu) override;
  const Champ_base& get_champ(const Motcle& nom) const override;
  void get_noms_champs_postraitables(Noms& nom, Option opt = NONE) const override;
  bool has_champ(const Motcle& nom, OBS_PTR(Champ_base) &ref_champ) const override;
  bool has_champ(const Motcle& nom) const override;

  // Writing u_star, Cisaillement_paroi, etc. to a separate file...
  void ouvrir_fichier_partage(EcrFicPartage&, const Nom&) const;
  void ouvrir_fichier_partage(EcrFicPartage&, const Nom&, const Nom&) const;
  // indicates whether shear is used or not
  virtual bool use_shear() const { return true; } // Generally true except for example for paroi_negligeable_XXX

protected:
  OBS_PTR(Modele_turbulence_hyd_base) mon_modele_turb_hyd;
  DoubleTab Cisaillement_paroi_;         // values of tangential stresses at
  // walls computed locally from u*
  DoubleVect tab_u_star_;                // values of u* computed locally
  DoubleVect tab_d_plus_;                // values of d+ computed locally
  mutable OWN_PTR(Champ_Fonc_base)  champ_u_star_;                                // field for post-processing
  mutable int nb_impr_ = 0, nb_impr0_ = 0;                        // print counter
  Champs_compris champs_compris_;
  OBS_PTR(Domaine_VF) le_dom_dis_;
  OBS_PTR(Domaine_Cl_dis_base) le_dom_Cl_dis_;
};

/*! @brief Associates a turbulence model to the object.
 *
 * @brief Associates a hydraulic turbulence model to the object.
 * @param (Modele_turbulence_hyd_base& le_modele) the hydraulic turbulence model to associate with the object
 */
inline void Turbulence_paroi_base::associer_modele(const Modele_turbulence_hyd_base& le_modele)
{
  mon_modele_turb_hyd = le_modele;
}

/*! @brief Simple call to int calculer_hyd(DoubleTab& ).
 *
 * @param (Champ_Inc_base& ch)
 * @return (int) propagated return code
 */
inline int Turbulence_paroi_base::calculer_hyd(Champ_Inc_base& ch)
{
  return calculer_hyd(ch.valeurs());
}

/*! @brief Simple call to int calculer_hyd(DoubleTab&, DoubleTab&).
 *
 * @param (Champ_Inc_base& ch1)
 * @param (Champ_Inc_base& ch2)
 * @return (int) propagated return code
 */
inline int Turbulence_paroi_base::calculer_hyd(Champ_Fonc_base& ch1, Champ_Fonc_base& ch2)
{
  return calculer_hyd(ch1.valeurs(), ch2.valeurs());
}

inline const DoubleTab& Turbulence_paroi_base::Cisaillement_paroi() const
{
  return Cisaillement_paroi_;
}

inline const DoubleVect& Turbulence_paroi_base::tab_u_star() const
{
  return tab_u_star_;
}

inline double Turbulence_paroi_base::tab_u_star(int face) const
{
  return tab_u_star_(face);
}

inline const DoubleVect& Turbulence_paroi_base::tab_d_plus() const
{
  return tab_d_plus_;
}

inline double Turbulence_paroi_base::tab_d_plus(int face) const
{
  return tab_d_plus_(face);
}

inline void erreur_non_convergence()
{
  Cerr << "TRUST stopped cause in the function calculer_u_plus of wall law." << finl;
  Cerr << "The iterative process of u* did not converge" << finl;
  if (Process::is_parallel())
    Cerr << " on the processor" << Process::me();
  else
    Cerr << ".";
  Cerr << finl << "The hydraulic calculation is may be diverging." << finl;
  Process::exit();
}

#endif
