#include <Fluide_Incompressible.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Milieu_composite_Euler.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Pb_Multiphase.h>
#include <Interprete.h>
#include <Domaine_VF.h>
#include<Fluide_reel_base.h>
#include <Pb_Euler.h>
#include <Momentum_Euler.h>
Implemente_instanciable(Milieu_composite_Euler, "Milieu_composite_Euler", Fluide_base);
// XD Milieu_composite_Euler listobj Milieu_composite_Euler -1 milieu_base 0 Composite medium made of several sub mediums.

Sortie& Milieu_composite_Euler::printOn(Sortie& os) const { return os; }

Entree& Milieu_composite_Euler::readOn(Entree& is)
{
  int i = 0;
  std::vector<std::pair<std::string, int>> especes; // string pour phase_nom. int : 0 pour liquide et 1 pour gaz
  Cerr << "Reading Milieu_composite_Euler..." << finl;
  Nom mot;
  is >> mot;
  if (mot != "{")
    {
      Cerr << "Milieu_composite_Euler : { expected instead of " << mot << finl;
      Process::exit();
    }

  for (is >> mot; mot != "}"; is >> mot)
    {
      if (Motcle(mot) == "POROSITES_CHAMP") is >> ch_porosites_;
      else if (Motcle(mot) == "DIAMETRE_HYD_CHAMP") is >> ch_diametre_hyd_;
      else if (Motcle(mot) == "POROSITES")
        {
          Cerr << "You should use porosites_champ and not porosites ! Call the 911 !" << finl;
          Process::exit();
        }
      else if (Motcle(mot) == "GRAVITE")
        {
          Cerr << que_suis_je() << " : gravity should not be defined in Pb_Multiphase ! Use source_qdm if you want gravity in QDM equation !" << finl;
          Process::exit();
        }
      else if (!mot.debute_par("saturation") && !mot.debute_par("interface")) // on ajout les phases
        {
          noms_phases_.add(mot);
          Cerr << "Milieu_composite_Euler: add phase " << mot << " ... " << finl;
          OWN_PTR(Fluide_base) fluide;
          fluide.typer_lire_simple(is, "Typing the fluid medium ...");

          if (fluide->has_porosites())
            {
              Cerr << que_suis_je() + " : porosity should be defined only once in the Milieu_composite_Euler block, not in " + fluide->que_suis_je() << finl;
              Process::exit();
            }
          if (fluide->a_gravite())
            {
              Cerr << que_suis_je() + " : gravity should be defined only once in the Milieu_composite_Euler block, not in " + fluide->que_suis_je() << finl;
              Process::exit();
            }
          if (fluide->has_hydr_diam())
            {
              Cerr << que_suis_je() + " : diametre_hyd_champ should be defined only once in the Milieu_composite_Euler block, not in " + fluide->que_suis_je() << finl;
              Process::exit();
            }
          if (!sub_type(Fluide_base,fluide.valeur()))
            {
              Cerr << que_suis_je() + " : only real fluids should be read and not a fluid of type " + fluide->que_suis_je() << finl;
              Process::exit();
            }

          fluide->set_id_composite(i++);
          fluide->nommer(mot); // XXX
          fluides_.push_back(fluide);
          especes.push_back(check_fluid_name(fluide->le_nom()));
        }
      else if (mot.debute_par("saturation")) // on ajout la saturation
        {
          has_saturation_ = true;
          Cerr << "Milieu_composite_Euler: add saturation " << mot << " ... " << finl;
          sat_lu_.typer_lire_simple(is, "Typing the saturation ...");
        }
      else // on ajout l'interface
        {
          has_interface_ = true;
          Cerr << "Milieu_composite_Euler: add interface " << mot << " ... " << finl;
          inter_lu_.typer_lire_simple(is, "Typing the interface ...");
          //inter_BN_.typer_lire_simple(is, "Typing the interface ...");

        }
    }

  // Sais pas s'il faut tester ca ou non ... Yannick, Antoine ! help :/
  if (has_saturation() && has_interface())
    {
      Cerr << "You define both interface and saturation in Milieu_composite_Euler ???" << finl;
      Process::exit();
    }

  // Traitement pour les interfaces
  const int N = (int)fluides_.size();
  for (int n = 0; n < N; n++)
    {
      std::vector<Interface_base *> inter;
      for (int m = 0; m < N; m++)
        {
          const std::string espn = especes[n].first, espm = especes[m].first;
          const int pn = especes[n].second, pm = especes[m].second;

          if (pn != pm && (has_interface() || (espn == espm && has_saturation())))
            {
              Cerr << "Interface between fluid " << n << " : " << fluides_[n]->le_nom() << " and " << m << " : " << fluides_[m]->le_nom() << finl;
              inter.push_back(&ref_cast(Interface_base, has_saturation_ ? sat_lu_.valeur() : inter_lu_.valeur()));
              const Saturation_base *sat = sub_type(Saturation_base, *inter.back()) ? &ref_cast(Saturation_base, *inter.back()) : nullptr;
              if (sat && sat->get_Pref() > 0) // pour loi en e = e0 + cp * (T - T0)
                {
                  const double hn = pn ? sat->Hvs(sat->get_Pref()) : sat->Hls(sat->get_Pref()),
                               hm = pm ? sat->Hvs(sat->get_Pref()) : sat->Hls(sat->get_Pref()),
                               T0 = sat->Tsat(sat->get_Pref());
                  fluides_[n]->set_h0_T0(hn, T0), fluides_[m]->set_h0_T0(hm, T0);
                }
            }
          else inter.push_back(nullptr);
        }
      tab_interface_.push_back(inter);
    }
  return is;
}







std::pair<std::string, int> Milieu_composite_Euler::check_fluid_name(const Nom& name)
{
  int phase = name.debute_par("gaz");
  Nom espece = phase ? name.getSuffix("gaz_") : name.getSuffix("liquide_");

  // Suppression des suffixes "_group1" et "_group2" pour iate 2 groups
  std::string nomjdd = espece.getString();
  std::string suffix1 = "_group1";
  std::string suffix2 = "_group2";
  if (nomjdd.size() >= suffix1.size() && nomjdd.compare(nomjdd.size() - suffix1.size(), suffix1.size(), suffix1) == 0)
    {
      nomjdd.erase(nomjdd.size() - suffix1.size(), suffix1.size());
    }
  else if (nomjdd.size() >= suffix2.size() && nomjdd.compare(nomjdd.size() - suffix2.size(), suffix2.size(), suffix2) == 0)
    {
      nomjdd.erase(nomjdd.size() - suffix2.size(), suffix2.size());
    }

  return std::make_pair(nomjdd, phase);
}

bool Milieu_composite_Euler::has_interface(int k, int l) const
{
  if (tab_interface_[k][l]) return true;
  return false;
}

Interface_base& Milieu_composite_Euler::get_interface(int k, int l) const
{
  return *(tab_interface_[k][l]);
}

bool Milieu_composite_Euler::has_saturation(int k, int l) const
{
  if (tab_interface_[k][l] && sub_type(Saturation_base, *tab_interface_[k][l])) return true;
  return false;
}

Saturation_base& Milieu_composite_Euler::get_saturation(int k, int l) const
{
  return ref_cast(Saturation_base, *(tab_interface_[k][l]));
}

const Fluide_base& Milieu_composite_Euler::get_fluid(const int i) const
{
  assert(i >= 0 && i < (int )fluides_.size());
  return fluides_[i].valeur();
}


Fluide_base& Milieu_composite_Euler::get_fluid(const int i)
{
  assert(i >= 0 && i < (int )fluides_.size());
  return fluides_[i].valeur();
}

int Milieu_composite_Euler::initialiser(const double temps)
{
  for (auto &itr : fluides_) itr->initialiser(temps);

  const bool res_en_T = equation_.count("alpha_energie_tot") ? true : false;

  //const Equation_base& eqn = res_en_T ? equation("temperature") : equation("enthalpie");
  const Equation_base& eqn = equation("alpha_energie_tot");
  Champ_Inc_base& ch_rho = ref_cast(Champ_Inc_base, ch_rho_.valeur()),
                  &ch_e = ref_cast(Champ_Inc_base, ch_e_int_.valeur()),
                   &ch_h_ou_T = ref_cast(Champ_Inc_base, ch_h_ou_T_.valeur());

  ch_rho.associer_eqn(eqn);
  ch_rho.init_champ_calcule(*this, calculer_masse_volumique);

  ch_e.associer_eqn(eqn);
  ch_e.init_champ_calcule(*this, calculer_energie_interne);

  ch_h_ou_T.associer_eqn(eqn);
  res_en_T ? ch_h_ou_T.init_champ_calcule(*this, calculer_enthalpie) : ch_h_ou_T.init_champ_calcule(*this, calculer_temperature_multiphase);

  t_init_ = temps;

  // XXX Elie Saikali : utile pour cas reprise !
  ch_rho_->changer_temps(temps);
  ch_e_int_->changer_temps(temps);
  ch_h_ou_T_->changer_temps(temps);
  return 1;
  //return Milieu_base::initialiser_porosite(temps);
}

void Milieu_composite_Euler::preparer_calcul()
{
  mettre_a_jour(t_init_);
  fluid_properties_initialised_ = true;
}

void Milieu_composite_Euler::discretiser(const Probleme_base& pb, const  Discretisation_base& dis)
{

  Cerr << "Composite medium discretization." << finl;

  const int N = (int)fluides_.size(), nc = pb.equation(0).inconnue().nb_valeurs_temporelles();
  const Domaine_dis_base& domaine_dis = pb.equation(0).domaine_dis();
  const double temps = pb.schema_temps().temps_courant();

  res_en_T_ = sub_type(Pb_Euler,pb) ? ref_cast(Pb_Euler,pb).resolution_en_T() : true;

  for (auto& itr : fluides_) itr->discretiser(pb, dis);  // a voir

  /* masse volumique, energie interne, enthalpie : champ_inc */
  OWN_PTR(Champ_Inc_base) rho_inc, ei_inc, h_ou_T_inc;
  dis.discretiser_champ("champ_elem", domaine_dis, "masse_volumique", "kg/m^3", N, nc, temps, rho_inc);
  dis.discretiser_champ("champ_elem", domaine_dis, "energie_interne", "J/m^3", N, nc, temps, ei_inc);
  if (res_en_T_)
    dis.discretiser_champ("champ_elem", domaine_dis, "enthalpie", "J/m^3", N, nc, temps, h_ou_T_inc);
  else
    dis.discretiser_champ("champ_elem", domaine_dis, "temperature", "C", N, nc, temps, h_ou_T_inc);

  ch_rho_ = rho_inc, ch_e_int_ = ei_inc, ch_h_ou_T_ = h_ou_T_inc;

  /* autres champs : champ_fonc */
  dis.discretiser_champ("champ_elem", domaine_dis, "viscosite_dynamique", "kg/m/s", N, temps, ch_mu_);
  dis.discretiser_champ("champ_elem", domaine_dis, "viscosite_cinematique", "m2/s", N, temps, ch_nu_);
  dis.discretiser_champ("champ_elem", domaine_dis, "diffusivite", "m2/s", N, temps, ch_alpha_);
  dis.discretiser_champ("champ_elem", domaine_dis, "alpha_fois_rho", "kg/m/s", N, temps, ch_alpha_fois_rho_);
  dis.discretiser_champ("champ_elem", domaine_dis, "conductivite", "W/m/K", N, temps, ch_lambda_);
  dis.discretiser_champ("champ_elem", domaine_dis, "capacite_calorifique", "J/kg/K", N, temps, ch_Cp_);
  dis.discretiser_champ("champ_elem", domaine_dis, "masse_volumique_melange", "kg/m^3", 1, temps, rho_m_);
  dis.discretiser_champ("champ_elem", domaine_dis, "enthalpie_melange", "J/m^3", 1, temps, h_m_);

  champs_compris_.ajoute_champ(ch_rho_);
  champs_compris_.ajoute_champ(ch_e_int_);
  champs_compris_.ajoute_champ(ch_h_ou_T_);

  std::vector<OWN_PTR(Champ_Don_base)* > fields = {&ch_mu_, &ch_nu_, &ch_lambda_, &ch_alpha_, &ch_alpha_fois_rho_, &ch_Cp_, &rho_m_, &h_m_};
  for (auto && f: fields) champs_compris_.ajoute_champ((*f).valeur());

  // on discretise les champs sigma / Tsat si besoin ...
  if ((int) tab_interface_.size() > 0)
    {
      Cerr << "Milieu_composite_Euler_Euler::" << __func__ << " ==> Surface tension discretization ..." << finl;
      if (has_saturation_) Cerr << "Milieu_composite_Euler_Euler::" << __func__ << " ==> Saturation temperature discretization ..." << finl;

      for (int k = 0; k < N; k++)
        for (int l = k + 1; l < N; l++)
          {
            int phase = fluides_[k]->le_nom().debute_par("gaz");
            Nom espece = phase ? fluides_[k].le_nom().getSuffix("gaz_") : fluides_[k].le_nom().getSuffix("liquide_");
            if (has_interface(k, l)) // OK si interf/saturation
              {
                Interface_base& inter = get_interface(k, l);
                inter.assoscier_pb(pb);
                Nom sig_nom = Nom("surface_tension_") + espece;
                inter.discretiser_sigma(sig_nom, temps);
                champs_compris_.ajoute_champ(inter.get_sigma_champ());
              }

            if (has_saturation(k, l)) // OK si saturation seulement
              {
                Saturation_base& sat = get_saturation(k, l);
                Nom Tsat_nom = Nom("Tsat_") + espece;
                sat.discretiser_Tsat(Tsat_nom, temps);
                champs_compris_.ajoute_champ(sat.get_Tsat_champ());
              }
          }
    }

  // Finalement, on discretise la porosite + diametre_hydro
  Milieu_base::discretiser_porosite(pb,dis);
  Milieu_base::discretiser_diametre_hydro(pb,dis);
  inter_lu_->assoscier_pb(pb);
}

void Milieu_composite_Euler::mettre_a_jour(double temps)
{

//  //  pr ession().mettre_a_jour(temps);
//
//
//
//  for (auto& itr : fluides_) itr->mettre_a_jour(temps);
//
//  ch_rho_->mettre_a_jour(temps);
//  ch_e_int_->mettre_a_jour(temps);
//  ch_h_ou_T_->mettre_a_jour(temps);
//
//  std::vector<OWN_PTR(Champ_Don_base)* > fields = {&ch_mu_, &ch_nu_, &ch_lambda_, &ch_alpha_, &ch_alpha_fois_rho_, &ch_Cp_, &rho_m_, &h_m_};
//  for (auto && f: fields) (*f)->mettre_a_jour(temps);
//
//  mettre_a_jour_tabs_Euler();
//
//  // MAJ sigma & Tsat
//  if ((int) tab_interface_.size() > 0)
//    {
//      const int N = (int) fluides_.size();
//      for (int k = 0; k < N; k++)
//        for (int l = k + 1; l < N; l++)
//          if (has_interface(k, l)) // OK si interf/saturation
//            get_interface(k, l).mettre_a_jour(temps);
//    }
//
//  Milieu_base::mettre_a_jour_porosite(temps);
}

//void Milieu_composite_Euler::mettre_a_jour_tabs()
//{}

void Milieu_composite_Euler::associer_equation(const Equation_base *eqn) const
{
  Fluide_base::associer_equation(eqn);
  // on fait suivre aux milieux sous-jacents (les fluide_reel en ont besoin!)
  for (const auto& itr : fluides_) itr->associer_equation(eqn);
}

int Milieu_composite_Euler::check_unknown_range() const
{
  int ok = 1;
  for (int i = 0; ok && i < (int)fluides_.size(); i++) ok &= fluides_[i]->check_unknown_range(); //chaque fluide controle
  return ok;
}

void Milieu_composite_Euler::calculer_masse_volumique(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv)
{
  const Milieu_composite_Euler& mil = ref_cast(Milieu_composite_Euler, obj);
  const Domaine_VF& zvf = ref_cast(Domaine_VF, mil.equation_.begin()->second->domaine_dis());
  int i, Ni = val.dimension_tot(0), Nb = bval.dimension_tot(0), n, N = (int)mil.fluides_.size();
  std::vector<const DoubleTab *> split(N);
  for (n = 0; n < N; n++) split[n] = &mil.fluides_[n]->masse_volumique().valeurs();
  for (i = 0; i < Ni; i++)
    for (n = 0; n < N; n++) val(i, n) = (*split[n])(i * (split[n]->dimension(0) > 1), 0);

  std::vector<DoubleTab> bsplit(N);
  for (n = 0; n < N; n++)
    if (mil.fluides_[n]->masse_volumique().a_un_domaine_dis_base())
      bsplit[n] = mil.fluides_[n]->masse_volumique().valeur_aux_bords();
    else bsplit[n].resize(bval.dimension_tot(0), 1), mil.fluides_[n]->masse_volumique().valeur_aux(zvf.xv_bord(), bsplit[n]);
  for (i = 0; i < Nb; i++)
    for (n = 0; n < N; n++) bval(i, n) = bsplit[n](i * (split[n]->dimension(0) > 1), 0);

  /* derivees */
  std::vector<const tabs_t *> split_der(N);
  for (n = 0; n < N; n++) split_der[n] = sub_type(Champ_Inc_base, mil.fluides_[n]->masse_volumique()) ? &ref_cast(Champ_Inc_base, mil.fluides_[n]->masse_volumique()).derivees() : nullptr;
  std::set<std::string> noms_der;
  for (n = 0; n < N; n++)
    if (split_der[n])
      for (auto &&n_d : *split_der[n]) noms_der.insert(n_d.first);
  for (auto &&nom : noms_der)
    {
      for (n = 0; n < N; n++) split[n] = split_der[n] && split_der[n]->count(nom) ? &split_der[n]->at(nom) : nullptr;
      DoubleTab& der = deriv[nom];
      for (der.resize(Ni, N), i = 0; i < Ni; i++)
        for (n = 0; n < N; n++) der(i, n) = split[n] ? (*split[n])(i * (split[n]->dimension(0) > 1)) : 0;
    }
}

void Milieu_composite_Euler::calculer_energie_interne(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv)
{
  const Milieu_composite_Euler& mil = ref_cast(Milieu_composite_Euler, obj);
  int i, Ni = val.dimension_tot(0), Nb = bval.dimension_tot(0), n, N = (int)mil.fluides_.size();
  std::vector<const DoubleTab *> split(N);
  for (n = 0; n < N; n++) split[n] = &mil.fluides_[n]->energie_interne().valeurs();
  for (i = 0; i < Ni; i++)
    for (n = 0; n < N; n++) val(i, n) = (*split[n])(i * (split[n]->dimension(0) > 1), 0);

  std::vector<DoubleTab> bsplit(N);
  for (n = 0; n < N; n++) bsplit[n] = mil.fluides_[n]->energie_interne().valeur_aux_bords();
  for (i = 0; i < Nb; i++)
    for (n = 0; n < N; n++) bval(i, n) = bsplit[n](i * (split[n]->dimension(0) > 1), 0);

  /* derivees */
  std::vector<const tabs_t *> split_der(N);
  for (n = 0; n < N; n++) split_der[n] = sub_type(Champ_Inc_base, mil.fluides_[n]->energie_interne()) ? &ref_cast(Champ_Inc_base, mil.fluides_[n]->energie_interne()).derivees() : nullptr;
  std::set<std::string> noms_der;
  for (n = 0; n < N; n++)
    if (split_der[n])
      for (auto &&n_d : *split_der[n]) noms_der.insert(n_d.first);
  for (auto &&nom : noms_der)
    {
      for (n = 0; n < N; n++) split[n] = split_der[n] && split_der[n]->count(nom) ? &split_der[n]->at(nom) : nullptr;
      DoubleTab& der = deriv[nom];
      for (der.resize(Ni, N), i = 0; i < Ni; i++)
        for (n = 0; n < N; n++) der(i, n) = split[n] ? (*split[n])(i * (split[n]->dimension(0) > 1)) : 0;
    }
}

void Milieu_composite_Euler::calculer_enthalpie(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv)
{
  const Milieu_composite_Euler& mil = ref_cast(Milieu_composite_Euler, obj);
  int i, Ni = val.dimension_tot(0), Nb = bval.dimension_tot(0), n, N = (int)mil.fluides_.size();
  std::vector<const DoubleTab *> split(N);
  for (n = 0; n < N; n++) split[n] = &mil.fluides_[n]->enthalpie().valeurs();
  for (i = 0; i < Ni; i++)
    for (n = 0; n < N; n++) val(i, n) = (*split[n])(i * (split[n]->dimension(0) > 1), 0);

  std::vector<DoubleTab> bsplit(N);
  for (n = 0; n < N; n++) bsplit[n] = mil.fluides_[n]->enthalpie().valeur_aux_bords();
  for (i = 0; i < Nb; i++)
    for (n = 0; n < N; n++) bval(i, n) = bsplit[n](i * (split[n]->dimension(0) > 1), 0);

  /* derivees */
  std::vector<const tabs_t *> split_der(N);
  for (n = 0; n < N; n++) split_der[n] = sub_type(Champ_Inc_base, mil.fluides_[n]->enthalpie()) ? &ref_cast(Champ_Inc_base, mil.fluides_[n]->enthalpie()).derivees() : nullptr;
  std::set<std::string> noms_der;
  for (n = 0; n < N; n++)
    if (split_der[n])
      for (auto &&n_d : *split_der[n]) noms_der.insert(n_d.first);
  for (auto &&nom : noms_der)
    {
      for (n = 0; n < N; n++) split[n] = split_der[n] && split_der[n]->count(nom) ? &split_der[n]->at(nom) : nullptr;
      DoubleTab& der = deriv[nom];
      for (der.resize(Ni, N), i = 0; i < Ni; i++)
        for (n = 0; n < N; n++) der(i, n) = split[n] ? (*split[n])(i * (split[n]->dimension(0) > 1)) : 0;
    }
}

void Milieu_composite_Euler::calculer_temperature_multiphase(const Objet_U& obj, DoubleTab& val, DoubleTab& bval, tabs_t& deriv)
{
  const Milieu_composite_Euler& mil = ref_cast(Milieu_composite_Euler, obj);
  int i, Ni = val.dimension_tot(0), Nb = bval.dimension_tot(0), n, N = (int)mil.fluides_.size();
  std::vector<const DoubleTab *> split(N);
  for (n = 0; n < N; n++) split[n] = &mil.fluides_[n]->temperature_multiphase().valeurs();
  for (i = 0; i < Ni; i++)
    for (n = 0; n < N; n++) val(i, n) = (*split[n])(i * (split[n]->dimension(0) > 1), 0);

  std::vector<DoubleTab> bsplit(N);
  for (n = 0; n < N; n++) bsplit[n] = mil.fluides_[n]->temperature_multiphase().valeur_aux_bords();
  for (i = 0; i < Nb; i++)
    for (n = 0; n < N; n++) bval(i, n) = bsplit[n](i * (split[n]->dimension(0) > 1), 0);

  /* derivees */
  std::vector<const tabs_t *> split_der(N);
  for (n = 0; n < N; n++) split_der[n] = sub_type(Champ_Inc_base, mil.fluides_[n]->temperature_multiphase()) ? &ref_cast(Champ_Inc_base, mil.fluides_[n]->temperature_multiphase()).derivees() : nullptr;
  std::set<std::string> noms_der;
  for (n = 0; n < N; n++)
    if (split_der[n])
      for (auto &&n_d : *split_der[n]) noms_der.insert(n_d.first);
  for (auto &&nom : noms_der)
    {
      for (n = 0; n < N; n++) split[n] = split_der[n] && split_der[n]->count(nom) ? &split_der[n]->at(nom) : nullptr;
      DoubleTab& der = deriv[nom];
      for (der.resize(Ni, N), i = 0; i < Ni; i++)
        for (n = 0; n < N; n++) der(i, n) = split[n] ? (*split[n])(i * (split[n]->dimension(0) > 1)) : 0;
    }
}

void Milieu_composite_Euler::abortTimeStep()
{
  Fluide_base::abortTimeStep();
  if (ch_e_int_.non_nul()) ch_e_int_->abortTimeStep();
  if (ch_h_ou_T_.non_nul()) ch_h_ou_T_->abortTimeStep();
}

bool Milieu_composite_Euler::initTimeStep(double dt)
{
  for (auto& itr : fluides_) itr->initTimeStep(dt);

  if (!equation_.size()) return true; //pas d'equation associee -> ???
  const Schema_Temps_base& sch = equation_.begin()->second->schema_temps(); //on recupere le schema en temps par la 1ere equation

  /* champs dont on doit creer des cases */
  std::vector<Champ_Inc_base *> vch;
  if (ch_rho_.non_nul() && sub_type(Champ_Inc_base, ch_rho_.valeur()))
    vch.push_back(&ref_cast(Champ_Inc_base, ch_rho_.valeur()));
  if (ch_e_int_.non_nul() && sub_type(Champ_Inc_base, ch_e_int_.valeur()))
    vch.push_back(&ref_cast(Champ_Inc_base, ch_e_int_.valeur()));
  if (ch_h_ou_T_.non_nul() && sub_type(Champ_Inc_base, ch_h_ou_T_.valeur()))
    vch.push_back(&ref_cast(Champ_Inc_base, ch_h_ou_T_.valeur()));

  for (auto &pch : vch)
    for (int i = 1; i <= sch.nb_valeurs_futures(); i++)
      {
        pch->changer_temps_futur(sch.temps_futur(i), i);
        pch->futur(i) = pch->valeurs();
      }
  return true;
}


void Milieu_composite_Euler::mettre_a_jour_tabs_Euler()
{
//  const int N = (int)fluides_.size();
//
//  /* viscosite dynamique */
//  {
//    DoubleTab& tab = ch_mu_->valeurs();
//    const int Nl = ch_mu_->valeurs().dimension_tot(0);
//    for (int n = 0; n < N; n++)
//      {
//        const Champ_base& ch_n = fluides_[n]->viscosite_dynamique();
//        const int cch = sub_type(Champ_Uniforme, ch_n);
//        for (int i = 0; i < Nl; i++) tab(i, n) = ch_n.valeurs()(!cch * i, 0);
//      }
//  }
//
//  /* viscosite cinematique */
//  {
//    DoubleTab& tab = ch_nu_->valeurs();
//    const int Nl = ch_nu_->valeurs().dimension_tot(0);
//    for (int n = 0; n < N; n++)
//      {
//        const Champ_base& ch_n = fluides_[n]->viscosite_cinematique();
//        const int cch = sub_type(Champ_Uniforme, ch_n);
//        for (int i = 0; i < Nl; i++) tab(i, n) = ch_n.valeurs()(!cch * i, 0);
//      }
//  }
//
//  /* diffusivite */
//  {
//    DoubleTab& tab = ch_alpha_->valeurs();
//    const int Nl = ch_alpha_->valeurs().dimension_tot(0);
//    for (int n = 0; n < N; n++)
//      {
//        const Champ_base& ch_n = fluides_[n]->diffusivite();
//        const int cch = sub_type(Champ_Uniforme, ch_n);
//        for (int i = 0; i < Nl; i++) tab(i, n) = ch_n.valeurs()(!cch * i, 0);
//      }
//  }
//
//  /* alpha fois rho */
//  {
//    DoubleTab& tab = ch_alpha_fois_rho_->valeurs();
//    const int Nl = ch_alpha_fois_rho_->valeurs().dimension_tot(0);
//    for (int n = 0; n < N; n++)
//      {
//        const Champ_base& ch_n = fluides_[n]->diffusivite_fois_rho();
//        const int cch = sub_type(Champ_Uniforme, ch_n);
//        for (int i = 0; i < Nl; i++) tab(i, n) = ch_n.valeurs()(!cch * i, 0);
//      }
//  }
//
//  /* conductivite */
//  {
//    DoubleTab& tab = ch_lambda_->valeurs();
//    const int Nl = ch_lambda_->valeurs().dimension_tot(0);
//    for (int n = 0; n < N; n++)
//      {
//        const Champ_base& ch_n = fluides_[n]->conductivite();
//        const int cch = sub_type(Champ_Uniforme, ch_n);
//        for (int i = 0; i < Nl; i++) tab(i, n) = ch_n.valeurs()(!cch * i, 0);
//      }
//  }
//
//  /* capacite calorifique */
//  {
//    DoubleTab& tab = ch_Cp_->valeurs();
//    const int Nl = ch_Cp_->valeurs().dimension_tot(0);
//    for (int n = 0; n < N; n++)
//      {
//        const Champ_base& ch_n = fluides_[n]->capacite_calorifique();
//        const int cch = sub_type(Champ_Uniforme, ch_n);
//        for (int i = 0; i < Nl; i++) tab(i, n) = ch_n.valeurs()(!cch * i, 0);
//      }
//  }
//
//  /* calcul des quantites melange */
//  DoubleTab& trm = rho_m_->valeurs(), &thm = h_m_->valeurs();
//  trm = 0, thm = 0;
//
//  const DoubleTab& a = equation("densite").inconnue().valeurs(), &r = ch_rho_->valeurs(),
//                   &ent = ch_h_ou_T_->valeurs();
//
//  const int Nl = rho_m_->valeurs().dimension_tot(0);
//
//  // masse volumique
//  for (int n = 0; n < N; n++)
//    for (int i = 0; i < Nl; i++) trm(i) += a(i, n) * r(i, n);
//
//  // enthalpie
//  for (int n = 0; n < N; n++)
//    for (int i = 0; i < Nl; i++) thm(i) += a(i, n) * r(i, n) * ent(i, n);
//  for (int i = 0; i < Nl; i++) thm(i) /= trm(i);
}

//////////////////////////////////////////::
///////////////////////////////////////////
//////////////////////////////////////////

void Milieu_composite_Euler::init_energie_tot(DoubleTab& alpha_energie_tot_jdd) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler,equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& alpha = equation("alpha").inconnue().valeurs();
  const DoubleTab& p = qdm.pression().valeurs();
  const int Nb_phase = (int)fluides_.size();
  const int Nb_elem = qdm.domaine_dis().nb_elem_tot();

  for (int n = 0 ; n < Nb_phase; n++)
    {
      const DoubleTab& U = qdm.vitesse_phase(n).valeurs();
      const Fluide_reel_base& phase=ref_cast(Fluide_reel_base,get_fluid(n));
      for (int i = 0; i <Nb_elem ; i++)
        {
          double nom_u2 =0 ;
          for (int d=0 ; d < dimension ; d++ ) nom_u2+= U(i,d)*U(i,d);
          alpha_energie_tot_jdd(i,n) = alpha(i,n) * phase.init_energie_tot(rho(i,n),nom_u2,p(i,n));
        }
    }
}

void Milieu_composite_Euler::calculer_pression(DoubleTab& p) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler,equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& alpha = equation("alpha").inconnue().valeurs();
  const DoubleTab& alpha_rhoE = equation("alpha_energie_tot").inconnue().valeurs();
  DoubleTab rhoE = alpha_rhoE;
  tab_divide_any_shape(rhoE,alpha); // @suppress("Function cannot be resolved")
  const int& Nb_phase = (int)fluides_.size();
  const int& Nb_elem = qdm.domaine_dis().nb_elem_tot();

  for (int n = 0 ; n < Nb_phase; n++)
    {
      const DoubleTab& U = qdm.vitesse_phase(n).valeurs();
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base,get_fluid(n));
      for (int i = 0; i <Nb_elem ; i++)
        {
          double nom_u2 =0 ;
          for (int d=0 ; d < dimension ; d++ ) nom_u2+= U(i,d)*U(i,d);
          p(i,n) = phase.calculer_pression(rho(i,n),nom_u2,rhoE(i,n));
          assert(p(i,n) > 0);
        }
    }
}

void Milieu_composite_Euler::calculer_vitesse_son(DoubleTab& c) const
{
  const Momentum_Euler& qdm = ref_cast(Momentum_Euler,equation("alpha_rho_u"));
  const DoubleTab& rho = ref_cast(Density_Euler,equation("alpha_rho")).densite().valeurs();
  const DoubleTab& p = qdm.pression().valeurs();

  const int& Nb_phase = (int)fluides_.size();
  const int& Nb_elem = qdm.domaine_dis().nb_elem_tot();
  for (int n = 0; n < Nb_phase; n++)
    {
      const Fluide_reel_base& phase = ref_cast(Fluide_reel_base,get_fluid(n));
      for (int i = 0; i <Nb_elem ; i++)
        {
          c(i,n) = phase.calculer_vitesse_son(rho(i,n),p(i,n));
        }
    }
}





