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
#include <Param.h>
#include <vector>

Implemente_instanciable(Momentum_Euler,"Momentum_Euler|QDM_Euler",Navier_Stokes_std);
// XD momentum_euler eqn_base qdm_euler -1 Momentum conservation equation for a multi-phase Euler problem where the unknown is the velocity
// XD attr termes_non_conservatifs bloc_op_non_conservativtifs non_conservative_terms 1 Keyword to alter the non-conservative scheme.

Sortie& Momentum_Euler::printOn(Sortie& is) const { return Equation_base::printOn(is); }

Entree& Momentum_Euler::readOn(Entree& is)
{
  terme_nconserv_.associer_eqn(*this);
  Equation_base::readOn(is);

  terme_convectif.set_fichier("Convection_qdm");
  terme_convectif.set_description("Momentum flow rate=Integral(rho*u*u*ndS) [N] if SI units used");
  terme_diffusif.set_fichier("Contrainte_visqueuse");
  terme_diffusif.set_description("Friction drag exerted by the fluid=Integral(-mu*(grad(u) +grad(u)^T)*ndS) [N] if SI units used");

  assert(le_fluide.non_nul());
  if (!sub_type(Fluide_base, le_fluide.valeur()))
    {
      Cerr << "ERROR : the Momentum_Euler equation can be associated only to a fluid." << finl;
      exit();
    }

  terme_convectif->set_incompressible(1);

  return is;
}

void Momentum_Euler::set_param(Param& param) const
{
  Equation_base::set_param(param);
  param.ajouter_non_std("diffusion", (this));
  param.ajouter_non_std("convection", (this));
  param.ajouter_condition("is_read_diffusion", "The diffusion operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_condition("is_read_convection", "The convection operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_non_std("solveur_pression", (this)); // XD attr solveur_pression solveur_sys_base solveur_pression 1 Linear pressure system resolution method.
  param.ajouter_non_std("termes_non_conservatifs|non_conservative_terms", (this));
}

int Momentum_Euler::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "termes_non_conservatifs|non_conservative_terms")
    {
      Cerr << "Reading and typing of the termes_non_conservatifs operator : " << finl;
      is >> terme_nconserv_;
    }
  else
    return Navier_Stokes_std::lire_motcle_non_standard(mot, is);
  return 1;
}

int Momentum_Euler::has_interface_blocs() const
{
  int ok = Navier_Stokes_std::has_interface_blocs();
  return ok;
}

void Momentum_Euler::mettre_a_jour(double temps)
{
  Equation_base::mettre_a_jour(temps);
}

bool Momentum_Euler::initTimeStep(double dt)
{
  Schema_Temps_base& sch = schema_temps();
  ConstDoubleTab_parts ppart(pression().valeurs());
  /* si pression_pa() est plus petit que pression() (ex. : variables auxiliaires PolyMAC_P0P1NC), alors on ne copie que la 1ere partie */
  const DoubleTab& p_red = pression_pa().valeurs().dimension_tot(0) < pression().valeurs().dimension_tot(0) ? ppart[0] : pression().valeurs();
  for (int i = 1; i <= sch.nb_valeurs_futures(); i++)
    {
      // Mise a jour du temps dans la pression
      pression().changer_temps_futur(sch.temps_futur(i), i);
      pression().futur(i) = pression().valeurs();
      pression_pa().changer_temps_futur(sch.temps_futur(i), i);
      pression_pa().futur(i) = p_red;
    }
  return Equation_base::initTimeStep(dt);
}

void Momentum_Euler::abortTimeStep()
{
  Equation_base::abortTimeStep();
}

void Momentum_Euler::discretiser_vitesse()
{
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());

  dis.vitesse(schema_temps(), domaine_dis(), la_vitesse, ref_cast(Pb_Euler, probleme()).nb_phases());

  noms_vit_phases_.dimensionner(pb.nb_phases()), vit_phases_.resize(pb.nb_phases());
  for (int i = 0; i < pb.nb_phases(); i++)
    {
      noms_vit_phases_[i] = Nom("vitesse_") + pb.nom_phase(i);
      if (vit_phases_[i].est_nul())
        {
          discretisation().discretiser_champ("vitesse", domaine_dis(), noms_vit_phases_[i], "m/s", dimension, 1, 0, vit_phases_[i]);
          champs_compris_.ajoute_champ(vit_phases_[i]);
        }
    }
}

void Momentum_Euler::discretiser_grad_p()
{
// Ne fait rien ! Est appele par defaut dans Navier_Stokes_std.discretiser() mais pas requis en Pb_Multiphase
// La dicretisation par dans le Momentum_Euler.creer_champ()
}

const Champ_Don_base& Momentum_Euler::diffusivite_pour_transport() const
{
  return le_fluide->viscosite_dynamique();
}

const Champ_base& Momentum_Euler::diffusivite_pour_pas_de_temps() const
{
  return le_fluide->viscosite_cinematique();
}

const Champ_base& Momentum_Euler::vitesse_pour_transport() const
{
  return la_vitesse;
}

/*! @brief Complete l'equation base, associe la pression a l'equation,
 *
 *     complete la divergence, le gradient et le solveur pression.
 *     Ajout de 2 termes sources: l'un representant la force centrifuge
 *     dans le cas axi-symetrique,l'autre intervenant dans la resolution
 *     en 2D axisymetrique
 *
 */
void Momentum_Euler::completer()
{
  Cerr << " Navier_Stokes_std::completer_deb" << finl;
  Equation_base::completer();
  la_pression->associer_domaine_cl_dis(le_dom_Cl_dis);
  Cerr << " Navier_Stokes_std::completer_fin" << finl;
  Cerr << "unknow field type  " << inconnue().que_suis_je() << finl;
  Cerr << "unknow field name  " << inconnue().le_nom() << finl;
  Cerr << "equation type " << inconnue().equation().que_suis_je() << finl;

  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());

  vitesse_son_.resize(dom.nb_elem_tot(), pb.nb_phases());

  vitesse_normale_.resize(dom.nb_faces_tot(), 2 * pb.nb_phases()); // dimension du tab à changer pour 3 pahses
}

void Momentum_Euler::get_noms_champs_postraitables(Noms& noms, Option opt) const
{
  Navier_Stokes_std::get_noms_champs_postraitables(noms, opt);

  Noms noms_compris;
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());
  for (int i = 0; i < pb.nb_phases(); i++)
    {
      noms_compris.add(noms_vit_phases_[i]);
    }
  if (opt == DESCRIPTION)
    Cerr << " Momentum_Euler : " << noms_compris << finl;
  else
    noms.add(noms_compris);
}

void Momentum_Euler::creer_champ(const Motcle& motlu)
{
  Navier_Stokes_std::creer_champ(motlu);
  int i = noms_vit_phases_.rang(motlu);
  if (i >= 0 && vit_phases_[i].est_nul())
    {
      discretisation().discretiser_champ("vitesse", domaine_dis(), noms_vit_phases_[i], "m/s", dimension, 1, 0, vit_phases_[i]);
      champs_compris_.ajoute_champ(vit_phases_[i]);
    }
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
        la_pression_en_pa->passe() = la_pression_en_pa->valeurs() = la_pression->passe() = la_pression->valeurs();
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
  pression().changer_temps(temps);
  //pression_pa().changer_temps(temps);

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
  const Discret_Thyd& dis = ref_cast(Discret_Thyd, discretisation());
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());
  int nb_valeurs_temp = schema_temps().nb_valeurs_temporelles();
  double temps = schema_temps().temps_courant();

  Cerr << "Hydraulic equation discretization " << finl;

  Cerr << "Velocity discretization" << finl;
  discretiser_vitesse();
  la_vitesse->add_synonymous(Nom("velocity"));
  champs_compris_.ajoute_champ(la_vitesse);

  //dis.pression(schema_temps(), domaine_dis(), la_pression);
  Cerr << "Pressure discretization" << finl;
  dis.discretiser_champ("pression", domaine_dis(), "pression", "Pa.m3/kg", pb.nb_phases(), nb_valeurs_temp, temps, la_pression);
  //la_pression->add_synonymous(Nom("P_star"));

  la_pression->fixer_nature_du_champ(pb.nb_phases() == 1 ? scalaire : pb.nb_phases() == dimension ? vectoriel : multi_scalaire); //pfft
  for (int i = 0; i < pb.nb_phases(); i++)
    la_pression->fixer_nom_compo(i, Nom("pression_") + pb.nom_phase(i));
  champs_compris_.ajoute_champ(la_pression);
  la_pression->associer_eqn(*this);

  Cerr << "QDM discretization" << finl;
  dis.discretiser_champ("vitesse", domaine_dis(), "alpha_rho_u", "kg/sm3", dimension * pb.nb_phases(), nb_valeurs_temp, temps, l_inco_ch_);
  champs_compris_.ajoute_champ(l_inco_ch_);

  //////////////////////////////////

  dis.pression_en_pa(schema_temps(), domaine_dis(), la_pression_en_pa);
  la_pression_en_pa->add_synonymous(Nom("Pressure"));
  champs_compris_.ajoute_champ(la_pression_en_pa);

  Equation_base::discretiser();
}

int Momentum_Euler::sauvegarder(Sortie& os) const
{
  int bytes = 0;
  bytes += Equation_base::sauvegarder(os);
  sauver();

  return bytes;
}

DoubleTab& Momentum_Euler::corriger_derivee_expl(DoubleTab& derivee)
{
  return derivee;
}

DoubleTab& Momentum_Euler::corriger_derivee_impl(DoubleTab& derivee)
{
  return derivee;
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
  ref_cast(Milieu_composite_Euler,milieu()).calculer_pression(pression().valeurs());
  ref_cast(Milieu_composite_Euler,milieu()).calculer_vitesse_son(vitesse_son());
}

double Momentum_Euler::calculer_pas_de_temps() const
{
  const Schema_Temps_base& sh = schema_temps();
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const IntTab& elem_faces = dom.elem_faces();
  const DoubleVect& surf = dom.face_surfaces();
  const DoubleVect& vol = dom.volumes();
  const int nb_phases = ref_cast(Pb_Euler,probleme()).nb_phases();
  const DoubleTab& c = vitesse_son();
  const DoubleTab& u_n = vitesse_normale();
  double dt = sh.pas_temps_max();

  DoubleTab dt_e(dom.nb_elem());
  for (int n = 0; n < nb_phases; n++)
    {
      for (int e = 0; e < dom.nb_elem(); e++)
        {
          double som = 0;
          for (int i = 0; i < elem_faces.line_size(); i++)
            {
              int f = elem_faces(e, i);
              som += (f >= 0) ? (fabs(u_n(f, n)) + c(e, n)) * surf(f) : 0;
            }
          dt_e(e) = vol(e) / som;
        }
      dt = std::min(dt, mp_min_vect(dt_e)); // @suppress("Invalid arguments") // @suppress("Function cannot be resolved")
    }
  return ((sh.facteur_securite_pas() * dt) < sh.temps_max() - sh.temps_courant()) ? dt : (sh.temps_max() - sh.temps_courant()) / sh.facteur_securite_pas();
}

void Momentum_Euler::init_alpha_rho_u()
{
  const DoubleTab& alpha_rho = ref_cast(Pb_Euler,probleme()).equation_masse().inconnue().valeurs();
  DoubleTab& alpha_rhoU = inconnue().valeurs();

  const int Nb_phase = ref_cast(Pb_Euler, probleme()).nb_phases(), D = Objet_U::dimension;
  for (int n = 0; n < Nb_phase; n++)
    {
      assert(vit_phases_[n].non_nul());
      DoubleTab_parts psrc(vitesse().valeurs()), pdst(vit_phases_[n]->valeurs());
      for (int i = 0; i < std::min(psrc.size(), pdst.size()); i++)
        {
          DoubleTab& src = psrc[i], &dst = pdst[i];
          assert(src.line_size() == Nb_phase * D);
          for (int j = 0; j < src.dimension_tot(0); j++)
            for (int d = 0; d < D; d++)
              {
                dst(j, d) = src(j, Nb_phase * d + n);
                alpha_rhoU(j, Nb_phase * d + n) = src(j, Nb_phase * d + n) * alpha_rho(j, n);
              }
        }
    }
}

void Momentum_Euler::calculer_vitesse()
{
  const DoubleTab& alpha_rho = ref_cast(Pb_Euler,probleme()).equation_masse().inconnue().valeurs();
  DoubleTab& U = vitesse().valeurs();
  const int Nb_phase = ref_cast(Pb_Euler, probleme()).nb_phases(), D = Objet_U::dimension;

  for (int n = 0; n < Nb_phase; n++)
    {
      DoubleTab_parts psrc(inconnue().valeurs()), pdst(vit_phases_[n]->valeurs());

      for (int i = 0; i < std::min(psrc.size(), pdst.size()); i++)
        {
          DoubleTab& src = psrc[i], &dst = pdst[i];
          assert(src.line_size() == Nb_phase * D);
          for (int j = 0; j < src.dimension_tot(0); j++)
            for (int d = 0; d < D; d++)
              {
                dst(j, d) = src(j, Nb_phase * d + n) / alpha_rho(j, n);
                U(j, Nb_phase * d + n) = src(j, Nb_phase * d + n) / alpha_rho(j, n);
              }
        }
    }
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
  const int Nb_phase = ref_cast(Pb_Euler, probleme()).nb_phases();
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const DoubleTab& U = vitesse().valeurs();
  const IntTab& f_e = dom.face_voisins();
  DoubleTab& u_n = vitesse_normale();
  assert(Objet_U::dimension == 2);
  for (int n = 0; n < Nb_phase; n++)
    {
      for (int f = 0; f < dom.nb_faces_tot(); f++)
        {
          // attention :  u_n( faces ,  0 ) est la vitesse normale left (au sens de la normale sortante) et u_n( faces ,  0 ) est vit norma right
          // CE N EST PAS	 la dimesion de l espae
          int el = f_e(f, 0), er = f_e(f, 1);
          double nx = dom.face_normales(f, 0) / dom.face_surfaces(f);
          double ny = dom.face_normales(f, 1) / dom.face_surfaces(f);
          u_n(f, n) = (el >= 0) ? U(el, n) * nx + U(el, n + Nb_phase) * ny : 123.123;
          u_n(f, n + Nb_phase) = (er >= 0) ? U(er, n) * nx + U(er, n + Nb_phase) * ny : 123.123;
        }
    }
}

void Momentum_Euler::flux(const int f, const int left_or_right, DoubleTab& res) const
{
  //left_or_right = 0 : left et 1 right;
  const Pb_Euler& pb = ref_cast(Pb_Euler, probleme());
  const Domaine_VF& dom = ref_cast(Domaine_VF, domaine_dis());
  const DoubleTab& vit_normale = vitesse_normale();
  const DoubleTab& p = pression().valeurs();
  const DoubleTab& alpha_rhoU = inconnue().valeurs();
  const DoubleTab& alpha = pb.equation_fraction().inconnue().valeurs();
  const int nb_phases = pb.nb_phases();
  const int e = dom.face_voisins(f, left_or_right);
  assert(res.dimension_tot(0) == alpha_rhoU.line_size());
  assert(Objet_U::dimension == 2);

  for (int n = 0; n < nb_phases; n++)
    {
      for (int d = 0; d < Objet_U::dimension; d++)
        {
          double vect_n_compo = dom.face_normales(f, d) / dom.face_surfaces(f);
          res(n + nb_phases * d) = vect_n_compo * alpha(e, n) * p(e, n) + alpha_rhoU(e, n + nb_phases * d) * vit_normale(f, n + left_or_right * nb_phases);
        }
    }
}

const Operateur& Momentum_Euler::operateur(int i) const
{
  if (i == 2)
    return terme_nconserv_;
  else
    return Navier_Stokes_std::operateur(i);
  return terme_convectif;
}

Operateur& Momentum_Euler::operateur(int i)
{
  if (i == 2)
    return terme_nconserv_;
  else
    return Navier_Stokes_std::operateur(i);
  return terme_convectif;
}
