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

#include <Milieu_composite_Euler.h>
#include <Schema_Temps_base.h>
#include <Champ_Composite.h>
#include <TRUSTTab_parts.h>
#include <Momentum_Euler.h>
#include <Discret_Thyd.h>
#include <Fluide_base.h>
#include <Pb_Euler.h>
#include <EChaine.h>
#include <Param.h>

Implemente_instanciable(Momentum_Euler,"Momentum_Euler|QDM_Euler",Navier_Stokes_std);
// XD momentum_euler eqn_base qdm_euler -1 Momentum conservation equation for a multi-phase Euler problem where the unknown is the velocity
// XD attr termes_non_conservatifs bloc_op_non_conservativtifs non_conservative_terms 1 Keyword to alter the non-conservative scheme.

Sortie& Momentum_Euler::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Momentum_Euler::readOn(Entree& is)
{
  assert(l_inco_ch_ && le_fluide);
  Equation_base::readOn(is);

  if (!sub_type(Pb_Euler, probleme()))
    {
      Cerr << "Equation " << que_suis_je() << " can only used with a pb of type Pb_Euler not " << probleme().que_suis_je() << " !!" << finl;
      Process::exit();
    }

  // si monophasique et termes_non_conservatifs non-lu ... on type sans error !
  const bool is_single_phase = (ref_cast(Pb_Euler, probleme()).nb_phases() == 1);
  if (is_single_phase && !terme_nconserv_)
    {
      EChaine str(" { negligeable } ");
      str >> terme_nconserv_;
      terme_nconserv_.associer_eqn(*this);
    }
  else if (!is_single_phase && !terme_nconserv_)
    {
      Cerr << "Error while reading " << que_suis_je() << " !!! non_conservative_terms operator is not read although it is required !!! Fix your data file ..." << finl;
      Process::exit();
    }

  terme_convectif.set_fichier("Convection_qdm");
  terme_convectif.set_description("Momentum flow rate=Integral(rho*u*u*ndS) [N] if SI units used");

  terme_nconserv_.set_fichier("Non_conservative_qdm");
  terme_nconserv_.set_description("Conribution of non_conservative operator in QDM equation");

  return is;
}

void Momentum_Euler::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("convection", (this), Param::REQUIRED);
  param.ajouter_non_std("termes_non_conservatifs|non_conservative_terms", (this));
}

int Momentum_Euler::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "termes_non_conservatifs|non_conservative_terms")
    {
      Cerr << "Reading and typing of the termes_non_conservatifs operator : " << finl;
      is >> terme_nconserv_;
      terme_nconserv_.associer_eqn(*this);
      return 1;
    }
  else
    return Navier_Stokes_std::lire_motcle_non_standard(mot, is);
}

const Operateur& Momentum_Euler::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_convectif;
    case 1:
      return terme_nconserv_;
    default :
      Cerr << "Error for Momentum_Euler::operateur(int i)" << finl;
      Cerr << "Momentum_Euler has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_convectif;
}

Operateur& Momentum_Euler::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_convectif;
    case 1:
      return terme_nconserv_;
    default :
      Cerr << "Error for Momentum_Euler::operateur(int i)" << finl;
      Cerr << "Momentum_Euler has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_convectif;
}

bool Momentum_Euler::initTimeStep(double dt)
{
  Schema_Temps_base& sch = schema_temps();
  // Mise a jour du temps dans la pression
  for (int i = 1; i <= sch.nb_valeurs_futures(); i++)
    {
      la_pression->changer_temps_futur(sch.temps_futur(i), i);
      la_pression->futur(i) = la_pression->valeurs();
    }
  return Equation_base::initTimeStep(dt);
}

void Momentum_Euler::abortTimeStep()
{
  Equation_base::abortTimeStep();
}

const Champ_Don_base& Momentum_Euler::diffusivite_pour_transport() const
{
  Process::exit("Momentum_Euler::diffusivite_pour_transport() should not be called !!! \n");
  return le_fluide->viscosite_dynamique();
}

const Champ_base& Momentum_Euler::diffusivite_pour_pas_de_temps() const
{
  Process::exit("Momentum_Euler::diffusivite_pour_transport() should not be called !!! \n");
  return le_fluide->viscosite_cinematique();
}

void Momentum_Euler::completer()
{
  Cerr << " Momentum_Euler::completer" << finl;
  Equation_base::completer();
  la_pression->associer_domaine_cl_dis(le_dom_Cl_dis);
  Cerr << "unknow field type  " << inconnue().que_suis_je() << finl;
  Cerr << "unknow field name  " << inconnue().le_nom() << finl;
  Cerr << "equation type " << inconnue().equation().que_suis_je() << finl;

  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());

  vitesse_normale_.resize(dom.nb_faces(), 2 * pb.nb_phases()); // dimension du tab à changer pour 3 pahses
  dom.creer_tableau_faces(vitesse_normale_);
  assert(vitesse_normale_.dimension(0) == dom.nb_faces());
  assert(vitesse_normale_.dimension_tot(0) == dom.nb_faces_tot());
  assert(vitesse_normale_.line_size() == 2 * pb.nb_phases());
}

void Momentum_Euler::get_noms_champs_postraitables(Noms& noms, Option opt) const
{
  Equation_base::get_noms_champs_postraitables(noms, opt);

  Noms noms_compris = champs_compris_.liste_noms_compris();

  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());
  for (int i = 0; i < pb.nb_phases(); i++)
    noms_compris.add(noms_vit_phases_[i]);

  if (opt == DESCRIPTION)
    Cerr << " Momentum_Euler : " << noms_compris << finl;
  else
    noms.add(noms_compris);
}

void Momentum_Euler::creer_champ(const Motcle& motlu)
{
  Equation_base::creer_champ(motlu);
  int i = noms_vit_phases_.rang(motlu);
  if (i >= 0 && !vit_phases_[i])
    {
      discretisation().discretiser_champ("vitesse", domaine_dis(), noms_vit_phases_[i], "m/s", dimension, 1, 0, vit_phases_[i]);
      champs_compris_.ajoute_champ(vit_phases_[i]);
    }
}

void Momentum_Euler::verifie_ch_init_nb_comp(const Champ_Inc_base& ch_ref, const int nb_comp) const
{
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());
  const Nature_du_champ nature = ch_ref.nature_du_champ();
  const Nom& nom = ch_ref.le_nom();

  if (nature == vectoriel)
    {
      if (nb_comp != pb.nb_phases() * Objet_U::dimension)
        {
          Cerr << "The nature of the field " << nom << " unknown to the equation name " << le_nom() << " is vector." << finl;
          Cerr << "The number of components readed for this field " << nb_comp << " is not compatible with its nature." << finl;
          Cerr << "It should read " << pb.nb_phases() * Objet_U::dimension << " components for this field." << finl;
          Process::exit();
        }
    }
  else
    Equation_base::verifie_ch_init_nb_comp(ch_ref, nb_comp);
}

Entree& Momentum_Euler::lire_cond_init(Entree& is)
{
  Cerr << "Reading of initial conditions\n";
  Nom nom;
  is >> nom;
  if (nom != "{")
    {
      Cerr << que_suis_je() << ": expected { instead of " << nom << finl;
      Process::exit();
    }
  int vit_lu = 0, press_lu = 0;
  for (is >> nom; nom != "}"; is >> nom)
    if (nom == "vitesse" || nom == "velocity")
      {
        OWN_PTR(Champ_Don_base) src;
        is >> src;
        if (src->que_suis_je() == "Champ_Composite")
          {
            const int nb_phases = ref_cast(Pb_Euler, probleme()).nb_phases(), nb_dim = ref_cast(Champ_Composite,src.valeur()).get_champ_composite_dim();
            if (nb_dim != nb_phases)
              {
                Cerr << que_suis_je() << ": velocity initial condition Champ_Composite should have " << nb_phases << " fields and not " << nb_dim << " !" << finl;
                Process::exit();
              }
          }

        verifie_ch_init_nb_comp(la_vitesse, src->nb_comp());
        la_vitesse->affecter(src), vit_lu = 1;
        la_vitesse->passe() = la_vitesse->valeurs();
      }
    else if (nom == "pression" || nom == "pressure")
      {
        OWN_PTR(Champ_Don_base) src;
        is >> src, verifie_ch_init_nb_comp(la_pression, src->nb_comp());
        la_pression->affecter(src);
        la_pression->passe() = la_pression->valeurs();
        press_lu = 1;
      }
    else
      {
        Cerr << que_suis_je() << ": expected vitesse|velocity|pression|pressure instead of " << nom << finl;
        Process::exit();
      }

  if (!vit_lu)
    {
      Cerr << que_suis_je() << ": velocity initial condition not found." << finl;
      Process::exit();
    }
  if (!press_lu)
    {
      Cerr << que_suis_je() << ": pressure initial condition not found." << finl;
      Process::exit();
    }

  return is;
}

int Momentum_Euler::preparer_calcul()
{
  Equation_base::preparer_calcul(); //pour eviter Navier_Stokes_std::preparer_calcul() !

  // XXX Elie Saikali : utile pour cas reprise !
  const double temps = schema_temps().temps_courant();
  la_pression->changer_temps(temps);

  return 1;
}

double Momentum_Euler::alpha_res() const
{
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());
  if (pb.nb_phases() == 1)
    return 0.;
  return -1.;
}

void Momentum_Euler::discretiser()
{
  Cerr << "Momentum_Euler discretization" << finl;
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());

  const double temps = schema_temps().temps_courant();
  const int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  const int N = pb.nb_phases();

  dis.vitesse(schema_temps(), domaine_dis(), la_vitesse, N);
  la_vitesse->fixer_nature_du_champ(vectoriel);
  la_vitesse->add_synonymous(Nom("velocity"));
  champs_compris_.ajoute_champ(la_vitesse);

  /* vitesse par phase */
  noms_vit_phases_.dimensionner(N);
  vit_phases_.resize(N);

  for (int i = 0; i < N; i++)
    {
      noms_vit_phases_[i] = Nom("vitesse_") + pb.nom_phase(i);
      if (!vit_phases_[i])
        {
          discretisation().discretiser_champ("vitesse", domaine_dis(), noms_vit_phases_[i], "m/s", dimension, 1, 0, vit_phases_[i]);
          champs_compris_.ajoute_champ(vit_phases_[i]);
        }
    }

  Cerr << "Pressure discretization" << finl;
  dis.discretiser_champ("pression", domaine_dis(), "pression", "Pa.m3/kg", N, nb_valeurs_temp, temps, la_pression);
  la_pression->fixer_nature_du_champ(N == 1 ? scalaire : multi_scalaire);
  la_pression->associer_eqn(*this);

  /* nom pression par phase */
  for (int i = 0; i < pb.nb_phases(); i++)
    la_pression->fixer_nom_compo(i, Nom("pression_") + pb.nom_phase(i));

  champs_compris_.ajoute_champ(la_pression);

  // on copie la structure //
  vitesse_son_ = la_pression->valeurs(); // nb_elem_tot * nb_phase

  Cerr << "Unknown alpha_rho_u discretization" << finl;
  dis.discretiser_champ("vitesse", domaine_dis(), "alpha_rho_u", "kg/sm3", dimension * N, nb_valeurs_temp, temps, l_inco_ch_);
  champs_compris_.ajoute_champ(l_inco_ch_);

  Equation_base::discretiser();
  Cerr << "Momentum_Euler discretization ==> ok" << finl;
}

int Momentum_Euler::sauvegarder(Sortie& os) const
{
  int bytes = 0;
  bytes += Equation_base::sauvegarder(os);
  sauver();

  return bytes;
}

void Momentum_Euler::mettre_a_jour_champs_conserves(double temps, int reset)
{
  Equation_base::mettre_a_jour_champs_conserves(temps);
  calculer_vitesse();
  calculer_vitesse_normale();
  mettre_a_jour_p_c();
}

void Momentum_Euler::mettre_a_jour_p_c()
{
  ref_cast(Milieu_composite_Euler,milieu()).calculer_pression(la_pression->valeurs());
  ref_cast(Milieu_composite_Euler,milieu()).calculer_vitesse_son(vitesse_son_);
}

double Momentum_Euler::calculer_pas_de_temps() const
{
  const Schema_Temps_base& sh = schema_temps();
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const IntTab& elem_faces = dom.elem_faces();
  const DoubleVect& surf = dom.face_surfaces();
  const DoubleVect& vol = dom.volumes();
  const int nb_phases = ref_cast(Pb_Euler,probleme()).nb_phases();
  const DoubleTab& u_n = vitesse_normale();
  double dt = sh.pas_temps_max();

  DoubleTrav dt_e(dom.nb_elem());
  for (int n = 0; n < nb_phases; n++)
    {
      for (int e = 0; e < dom.nb_elem(); e++)
        {
          double som = 0;
          for (int i = 0; i < elem_faces.line_size(); i++)
            {
              int f = elem_faces(e, i);
              som += (f >= 0) ? (fabs(u_n(f, n)) + vitesse_son_(e, n)) * surf(f) : 0;
            }
          dt_e(e) = vol(e) / som;
        }
      dt = std::min(dt, mp_min_vect(dt_e)); // @suppress("Invalid arguments") // @suppress("Function cannot be resolved")
    }
  return ((sh.facteur_securite_pas() * dt) < sh.temps_max() - sh.temps_courant()) ? dt : (sh.temps_max() - sh.temps_courant()) / sh.facteur_securite_pas();
}

void Momentum_Euler::init_alpha_rho_u()
{
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const DoubleTab& alpha_rho = ref_cast(Pb_Euler,probleme()).equation_masse().inconnue().valeurs();
  const DoubleTab& vit = vitesse().valeurs();
  const int Nb_phase = ref_cast(Pb_Euler, probleme()).nb_phases(), D = Objet_U::dimension;
  assert(vit.line_size() == Nb_phase * D);
  assert(vit.dimension(0) == dom.nb_elem());

  DoubleTab& alpha_rhoU = inconnue().valeurs();
  for (int n = 0; n < Nb_phase; n++)
    {
      assert(vit_phases_[n]);
      DoubleTab& vit_phase = vit_phases_[n]->valeurs();

      for (int j = 0; j < dom.nb_elem(); j++)
        for (int d = 0; d < D; d++)
          {
            vit_phase(j, d) = vit(j, Nb_phase * d + n);
            alpha_rhoU(j, Nb_phase * d + n) = vit(j, Nb_phase * d + n) * alpha_rho(j, n);
          }
      vit_phase.echange_espace_virtuel();
    }
  alpha_rhoU.echange_espace_virtuel();
}

void Momentum_Euler::calculer_vitesse()
{
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const DoubleTab& alpha_rho_U = inconnue().valeurs();
  const DoubleTab& alpha_rho = ref_cast(Pb_Euler,probleme()).equation_masse().inconnue().valeurs();
  const int Nb_phase = ref_cast(Pb_Euler, probleme()).nb_phases(), D = Objet_U::dimension;
  assert(alpha_rho_U.line_size() == Nb_phase * D);
  assert(alpha_rho_U.dimension(0) == dom.nb_elem());

  DoubleTab& vit = vitesse().valeurs();
  for (int n = 0; n < Nb_phase; n++)
    {
      assert(vit_phases_[n]);
      DoubleTab& vit_phase = vit_phases_[n]->valeurs();

      for (int j = 0; j < dom.nb_elem(); j++)
        for (int d = 0; d < D; d++)
          {
            vit_phase(j, d) = alpha_rho_U(j, Nb_phase * d + n) / alpha_rho(j, n);
            vit(j, Nb_phase * d + n) = alpha_rho_U(j, Nb_phase * d + n) / alpha_rho(j, n);
          }
      vit_phase.echange_espace_virtuel();
    }
  vit.echange_espace_virtuel();
}

const Champ_Inc_base& Momentum_Euler::vitesse_phase(const int i) const
{
  assert(i >= 0 && i < (int )vit_phases_.size());
  return vit_phases_[i].valeur();
}

Champ_Inc_base& Momentum_Euler::vitesse_phase(const int i)
{
  assert(i >= 0 && i < (int )vit_phases_.size());
  return vit_phases_[i].valeur();
}

void Momentum_Euler::calculer_vitesse_normale()
{
  assert(Objet_U::dimension == 2);
  const int Nb_phase = ref_cast(Pb_Euler, probleme()).nb_phases();
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const DoubleTab& vit = vitesse().valeurs();
  const IntTab& f_e = dom.face_voisins();
  for (int n = 0; n < Nb_phase; n++)
    for (int f = 0; f < dom.nb_faces(); f++)
      {
        // attention :  u_n( faces ,  0 ) est la vitesse normale left (au sens de la normale sortante) et u_n( faces ,  0 ) est vit norma right
        // CE N EST PAS	 la dimesion de l espae
        int el = f_e(f, 0), er = f_e(f, 1);
        double nx = dom.face_normales(f, 0) / dom.face_surfaces(f);
        double ny = dom.face_normales(f, 1) / dom.face_surfaces(f);
        vitesse_normale_(f, n) = (el >= 0) ? vit(el, n) * nx + vit(el, n + Nb_phase) * ny : -123.123;
        vitesse_normale_(f, n + Nb_phase) = (er >= 0) ? vit(er, n) * nx + vit(er, n + Nb_phase) * ny : -123.123;
      }

  vitesse_normale_.echange_espace_virtuel();
}
