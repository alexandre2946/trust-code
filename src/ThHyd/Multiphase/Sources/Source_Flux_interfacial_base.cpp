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

#include <Source_Flux_interfacial_base.h>
#include <Viscosite_turbulente_base.h>
#include <Flux_interfacial_base.h>
#include <Changement_phase_base.h>
#include <Operateur_Diff_base.h>
#include <Champ_Inc_P0_base.h>
#include <Aire_interfaciale.h>
#include <Milieu_composite.h>
#include <Champ_Face_base.h>
#include <Champ_Uniforme.h>
#include <Pb_Multiphase.h>
#include <Synonyme_info.h>
#include <Matrix_tools.h>
#include <EOS_to_TRUST.h>
#include <Array_tools.h>
#include <Domaine_VF.h>
#include <Domaine.h>

Implemente_base(Source_Flux_interfacial_base,"Source_Flux_interfacial_base", Sources_Multiphase_base);

Sortie& Source_Flux_interfacial_base::printOn(Sortie& os) const { return os; }

Entree& Source_Flux_interfacial_base::readOn(Entree& is)
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  if (!pbm.has_correlation("flux_interfacial"))
    Process::exit(que_suis_je() + " : a flux_interfacial correlation must be defined in the global correlations { } block!");

  const bool res_en_T = pbm.resolution_en_T();
  if (!res_en_T) Process::exit("Source_Flux_interfacial_base::readOn NOT YET PORTED TO ENTHALPY EQUATION ! TODO FIXME !!");

  for (int n = 0; n < pbm.nb_phases(); n++)
    {
      if ((( pbm.nom_phase(n).finit_par("group1"))) || (( pbm.nom_phase(n).finit_par("group2"))))  mod2grp = 1.;
    }

  correlation_ = pbm.get_correlation("flux_interfacial");

  dv_min = ref_cast(Flux_interfacial_base, correlation_.valeur()).dv_min();

  return is;
}

void Source_Flux_interfacial_base::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Champ_Inc_P0_base& ch = ref_cast(Champ_Inc_P0_base, equation().inconnue());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, equation().domaine_dis());
  const DoubleTab& inco = ch.valeurs();

  /* we must be able to add/subtract equations between components */
  int i, e, n, N = inco.line_size();
  if (N == 1) return;
  std::set<int> idx;
  for (auto &&n_m : matrices)
    if (n_m.first.find("/") == std::string::npos) /* ignore unknowns coming from other problems */
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first).valeurs();
        const int M = dep.line_size();
        Stencil sten(0, 2);

        if (n_m.first == "temperature" || n_m.first == "pression" || n_m.first == "alpha" || n_m.first == "interfacial_area" ) /* temperature/pression: dependance locale */
          for (e = 0; e < domaine.nb_elem(); e++)
            for (n = 0; n < N; n++)
              for (int m = 0; m < M; m++) sten.append_line(N * e + n, M * e + m);
        if (mat.nb_colonnes())
          for (e = 0; e < domaine.nb_elem(); e++) /* other variables: components may be mixed */
            {
              for (idx.clear(), n = 0, i = N * e; n < N; n++, i++)
                for (auto j = mat.get_tab1()(i) - 1; j < mat.get_tab1()(i + 1) - 1; j++)
                  idx.insert(mat.get_tab2()(j) - 1); // idx: set of columns on which at least one of the N component rows at e depends
              for (n = 0, i = N * e; n < N; n++, i++)
                for (auto &&x : idx) sten.append_line(i, x); // add this dependency to all rows
            }
        tableau_trier_retirer_doublons(sten);
        Matrix_tools::allocate_morse_matrix(inco.size_totale(), dep.size_totale(), sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Source_Flux_interfacial_base::completer()
{
  const Domaine_VF& domaine = ref_cast(Domaine_VF, equation().domaine_dis());
  int N = equation().inconnue().valeurs().line_size();
  if (!sub_type(Source_Flux_interfacial_base, equation().sources().dernier().valeur()))
    Process::exit(que_suis_je() + " : Source_Flux_interfacial_base must be the last source term in the source term declaration list of the " + equation().que_suis_je() + " equation ! ");

  if (sub_type(Energie_Multiphase, equation()))
    {
      qpi_.resize(0, N, N), domaine.domaine().creer_tableau_elements(qpi_);
      dT_qpi_.resize(0, N, N, N), domaine.domaine().creer_tableau_elements(dT_qpi_);
      da_qpi_.resize(0, N, N, N), domaine.domaine().creer_tableau_elements(da_qpi_);
      dp_qpi_.resize(0, N, N), domaine.domaine().creer_tableau_elements(dp_qpi_);
    }

  if(ref_cast(Operateur_Diff_base, equation().probleme().equation(0).operateur(0).l_op_base()).is_turb()) is_turb_ = 1;
}

DoubleTab& Source_Flux_interfacial_base::qpi() const
{
  if (sub_type(Energie_Multiphase, equation())) return qpi_;
  return ref_cast(Source_Flux_interfacial_base, ref_cast(Energie_Multiphase, ref_cast(Pb_Multiphase, equation().probleme()).equation_energie()).sources().dernier().valeur()).qpi();
}

DoubleTab& Source_Flux_interfacial_base::dT_qpi() const
{
  if (sub_type(Energie_Multiphase, equation())) return dT_qpi_;
  return ref_cast(Source_Flux_interfacial_base, ref_cast(Energie_Multiphase, ref_cast(Pb_Multiphase, equation().probleme()).equation_energie()).sources().dernier().valeur()).dT_qpi();
}

DoubleTab& Source_Flux_interfacial_base::da_qpi() const
{
  if (sub_type(Energie_Multiphase, equation())) return da_qpi_;
  return ref_cast(Source_Flux_interfacial_base, ref_cast(Energie_Multiphase, ref_cast(Pb_Multiphase, equation().probleme()).equation_energie()).sources().dernier().valeur()).da_qpi();
}

DoubleTab& Source_Flux_interfacial_base::dp_qpi() const
{
  if (sub_type(Energie_Multiphase, equation())) return dp_qpi_;
  return ref_cast(Source_Flux_interfacial_base, ref_cast(Energie_Multiphase, ref_cast(Pb_Multiphase, equation().probleme()).equation_energie()).sources().dernier().valeur()).dp_qpi();
}

void Source_Flux_interfacial_base::mettre_a_jour(double temps)
{
  qpi_ = 0, dT_qpi_ = 0, da_qpi_ = 0, dp_qpi_ = 0;
}

void Source_Flux_interfacial_base::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Pb_Multiphase& pbm = ref_cast(Pb_Multiphase, equation().probleme());
  const Champ_Inc_P0_base& ch = ref_cast(Champ_Inc_P0_base, equation().inconnue());
  // Matrice_Morse *mat = matrices.count(ch.le_nom().getString()) ? matrices.at(ch.le_nom().getString()) : nullptr;
  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, equation().domaine_dis());
  const DoubleVect& pe = milc.porosite_elem(), &ve = domaine.volumes();
  const tabs_t& der_h = ref_cast(Champ_Inc_base, milc.enthalpie()).derivees();
  const Champ_base& ch_rho = milc.masse_volumique();
  const Champ_Inc_base& ch_alpha = pbm.equation_masse().inconnue(), &ch_a_r = pbm.equation_masse().champ_conserve(),
                        &ch_temp = pbm.equation_energie().inconnue(), &ch_p = ref_cast(QDM_Multiphase, pbm.equation_qdm()).pression(),
                         *pch_rho = sub_type(Champ_Inc_base, ch_rho) ? &ref_cast(Champ_Inc_base, ch_rho) : nullptr;

  const DoubleTab& inco = ch.valeurs(), &alpha = ch_alpha.valeurs(), &press = ch_p.valeurs(), &temp  = ch_temp.valeurs(), &temp_p  = ch_temp.passe(),
                   &h = milc.enthalpie().valeurs(), *dP_h = der_h.count("pression") ? &der_h.at("pression") : nullptr, *dT_h = der_h.count("temperature") ? &der_h.at("temperature") : nullptr,
                    &lambda = milc.conductivite().passe(), &mu = milc.viscosite_dynamique().passe(), &rho = milc.masse_volumique().passe(), &Cp = milc.capacite_calorifique().passe(),
                     &p_ar = ch_a_r.passe(), &a_r = ch_a_r.valeurs(), &qi = qpi(), &dTqi = dT_qpi(), &daqi = da_qpi(), &dpqi = dp_qpi(),
                      *d_bulles = (equation().probleme().has_champ("diametre_bulles")) ? &equation().probleme().get_champ("diametre_bulles").valeurs() : nullptr,
                       *k_turb = (equation().probleme().has_champ("k")) ? &equation().probleme().get_champ("k").passe() : nullptr;

  Matrice_Morse *Mp = matrices.count("pression")    ? matrices.at("pression")    : nullptr,
                 *Mt = matrices.count("temperature") ? matrices.at("temperature") : nullptr,
                  *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr,
                   *Mai = matrices.count("interfacial_area") ? matrices.at("interfacial_area") : nullptr;

  int i, col, e, d, D = dimension, k, l, n, N = inco.line_size(), is_therm;
  const int cL = (lambda.dimension_tot(0) == 1), cM = (mu.dimension_tot(0) == 1), cR = (rho.dimension_tot(0) == 1), cCp = (Cp.dimension_tot(0) == 1);

  const Flux_interfacial_base& correlation_fi = ref_cast(Flux_interfacial_base, correlation_.valeur());
  const Changement_phase_base *correlation_G = pbm.has_correlation("changement_phase") ? &ref_cast(Changement_phase_base, pbm.get_correlation("changement_phase")) : nullptr;
  double dt = equation().schema_temps().pas_de_temps(), alpha_min = 1.e-6;

  // Turbulent viscosity for correlations that need it
  DoubleTrav nut;
  if (is_turb_)
    {
      nut.resize(0, N);
      MD_Vector_tools::creer_tableau_distribue(equation().inconnue().valeurs().get_md_vector(), nut); //Necessary to compare size in eddy_viscosity()
      const Viscosite_turbulente_base& corr_visc_turb = ref_cast(Viscosite_turbulente_base, *ref_cast(Operateur_Diff_base, equation().probleme().equation(0).operateur(0).l_op_base()).correlation_viscosite_turbulente());
      corr_visc_turb.eddy_viscosity(nut);
    }

  /* phase-change limiter: limit gamma to avoid alpha_k < 0 in any phase */
  /* to do so, assemble the mass equation without phase change */
  DoubleTrav sec_m(alpha); // residuals
  std::map<std::string, Matrice_Morse> mat_m_stockees;
  matrices_t mat_m; // derivatives
  for (auto &&n_m : matrices)
    if (n_m.first.find("/") == std::string::npos) /* ignore unknowns coming from other problems */
      {
        Matrice_Morse& dst = mat_m_stockees[n_m.first], &src = *n_m.second;
        dst.get_set_tab1().ref_array(src.get_set_tab1()); // same stencils as those of the current equation
        dst.get_set_tab2().ref_array(src.get_set_tab2());
        dst.set_nb_columns(src.nb_colonnes());
        dst.get_set_coeff().resize(src.get_set_coeff().size()); //zero coefficients
        mat_m[n_m.first] = &dst;
      }
  const Masse_Multiphase& eq_m = ref_cast(Masse_Multiphase, pbm.equation_masse());
  for (i = 0; i < eq_m.nombre_d_operateurs(); i++) /* all operators */
    eq_m.operateur(i).l_op_base().ajouter_blocs(mat_m, sec_m, semi_impl);
  for (i = 0; i < eq_m.sources().size(); i++)
    if (!sub_type(Source_Flux_interfacial_base, eq_m.sources()(i).valeur())) /* all sources except the interfacial flux */
      eq_m.sources()(i)->ajouter_blocs(mat_m, sec_m, semi_impl);
  std::vector<std::array<Matrice_Morse *, 2>> vec_m; //vector "source matrix, destination matrix"
  for (auto &&n_m : matrices)
    if (n_m.first.find("/") == std::string::npos && mat_m[n_m.first]->get_tab1().size_array() > 1) vec_m.push_back({{ mat_m[n_m.first], n_m.second }});

  /* elements */
  //coefficients and many derivatives...
  DoubleTrav dT_G(N), da_G(N), nv(N, N);
  Flux_interfacial_base::input_t in;
  Flux_interfacial_base::output_t out;
  DoubleTab& hi = out.hi, &dT_hi = out.dT_hi, &da_hi = out.da_hi, &dP_hi = out.dp_hi;
  hi.resize(N, N), dT_hi.resize(N, N, N), da_hi.resize(N, N, N), dP_hi.resize(N, N), in.v.resize(N, D);

  // And for the span methods of the Saturation class
  const int nbelem = domaine.nb_elem(), nb_max_sat =  N * (N-1) /2; // arithmetic series!
  DoubleTrav Ts_tab(nbelem,nb_max_sat), dPTs_tab(nbelem,nb_max_sat), Hvs_tab(nbelem,nb_max_sat), Hls_tab(nbelem,nb_max_sat), dPHvs_tab(nbelem,nb_max_sat), dPHls_tab(nbelem,nb_max_sat), Lvap_tab(nbelem,nb_max_sat), dP_Lvap_tab(nbelem,nb_max_sat), Sigma_tab(nbelem,nb_max_sat);

  // fill velocity at elem tab
  DoubleTab pvit_elem(0, N * D);
  domaine.domaine().creer_tableau_elements(pvit_elem);
  const Champ_Face_base& ch_vit = ref_cast(Champ_Face_base,ref_cast(Pb_Multiphase, equation().probleme()).equation_qdm().inconnue());
  ch_vit.get_elem_vector_field(pvit_elem);

  // fill the arrays ...
  for (k = 0; k < N; k++)
    for (l = k + 1; l < N; l++)
      if (milc.has_saturation(k, l))
        {
          Saturation_base& z_sat = milc.get_saturation(k, l);
          const int ind_trav = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // upper triangular matrix indexing!
          // XXX XXX XXX
          // Warning: this is dangerous! For now, we assume the pressure field has 1 component. However, the size of res is nb_max_sat*nbelem !!
          // Also, the Span passed for the pressure field is nbelem, not nbelem_tot ....
          assert(press.line_size() == 1);

          // retrieve Tsat and sigma ...
          const DoubleTab& sig = z_sat.get_sigma_tab(), &tsat = z_sat.get_Tsat_tab();

          // fill in the good case
          for (int ii = 0; ii < nbelem; ii++)
            {
              Ts_tab(ii, ind_trav) = tsat(ii);
              Sigma_tab(ii, ind_trav) = sig(ii);
            }

          // fill in the rest
          MSatSpanD sats_all = { };
          sats_all.insert( { SAT::T_SAT_DP, dPTs_tab.get_span() });
          sats_all.insert( { SAT::HV_SAT, Hvs_tab.get_span() });
          sats_all.insert( { SAT::HL_SAT, Hls_tab.get_span() });
          sats_all.insert( { SAT::HV_SAT_DP, dPHvs_tab.get_span() });
          sats_all.insert( { SAT::HL_SAT_DP, dPHls_tab.get_span() });
          sats_all.insert( { SAT::LV_SAT, Lvap_tab.get_span() });
          sats_all.insert( { SAT::LV_SAT_DP, dP_Lvap_tab.get_span() });

          ConstDoubleTab_parts ppart(press);
          z_sat.compute_all_flux_interfacial_pb_multiphase(ppart[0].get_span() /* elem reel */, sats_all, nb_max_sat, ind_trav);
        }

  for (e = 0; e < domaine.nb_elem(); e++)
    {
      double vol = pe(e) * ve(e), x, G = 0, dh = milc.diametre_hydraulique_elem()(e), dP_G = 0.; // E. Saikali: initialise dP_G here otherwise memory leak ...
      for (in.v = 0, d = 0; d < D; d++)
        for (n = 0; n < N; n++)
          in.v(n, d) = pvit_elem(e, N * d + n);
      for (nv = 0, d = 0; d < D; d++)
        for (n = 0; n < N; n++)
          for (k = 0 ; k<N ; k++) nv(n, k) += std::pow(pvit_elem(e, N * d + n) - ((n!=k) ? pvit_elem(e, N * d + k) : 0) , 2); // nv(n,n) = ||v(n)||, nv(n, k!=n) = ||v(n)-v(k)||
      for (n = 0; n < N; n++)
        for (k = 0 ; k<N ; k++) nv(n, k) = std::max(sqrt(nv(n, k)), dv_min);
      // interface exchange coefficients (explicit)
      in.dh = dh, in.alpha = &alpha(e, 0), in.T = &temp(e, 0),  in.T_passe = &temp_p(e, 0), in.p = press(e, 0), in.nv = &nv(0, 0), in.h = &h(e, 0), in.dT_h = dT_h ? &(*dT_h)(e, 0) : nullptr, in.dP_h = dP_h ? &(*dP_h)(e, 0) : nullptr;
      in.lambda = &lambda(!cL * e, 0), in.mu = &mu(!cM * e, 0), in.rho = &rho(!cR * e, 0), in.Cp = &Cp(!cCp * e, 0), in.e = e, in.Lvap = &Lvap_tab(e, 0), in.dP_Lvap = &dP_Lvap_tab(e, 0), in.sigma = &Sigma_tab(e,0), in.Tsat = &Ts_tab(e,0), in.dP_Tsat = &dPTs_tab(e,0);
      in.d_bulles = (d_bulles) ? &(*d_bulles)(e,0) : nullptr, in.k_turb = (k_turb) ? &(*k_turb)(e,0) : nullptr, in.nut = (is_turb_) ? &nut(e,0) : nullptr;
      correlation_fi.coeffs(in, out);

      for (k = 0; k < N; k++)
        for (l = k + 1; l < N; l++)
          if (milc.has_saturation(k, l)) //flux between phase k <-> phase l if saturation
            {
              Saturation_base& sat = milc.get_saturation(k, l);
              const int i_sat = (k*(N-1)-(k-1)*(k)/2) + (l-k-1); // upper triangular matrix indexing

              if (correlation_G) /* phase change rate from a correlation */
                G = correlation_G->calculer(k, l, dh, &alpha(e, 0), &temp(e, 0), press(e, 0), &nv(0), &lambda(!cL * e, 0), &mu(!cM * e, 0), &rho(!cM * e, 0), &Cp(!cCp * e, 0),
                                            sat, dT_G, da_G, dP_G);

              /* thermal limit */
              double Tk = temp(e, k), Tl = temp(e, l), Ts = Ts_tab(e,i_sat), dP_Ts = dPTs_tab(e,i_sat), // temperature on each side + Tsat + derivatives
                     phi = hi(k, l) * (Tk - Ts) + hi(l, k) * (Tl - Ts) + qi(e, k, l) / vol, L = (phi < 0 ? h(e, l) : Hvs_tab(e,i_sat)) - (phi > 0 ? h(e, k) : Hls_tab(e,i_sat));
              if ((is_therm = !correlation_G || std::fabs(G) > std::fabs(phi / L))) G = phi / L;
              /* phase enthalpies (in the direction corresponding to the chosen mode for G) */

              /* Is G limited by evanescence on side k or l? */
              int n_lim = G > 0 ? k : l, sgn = G > 0 ? 1 : -1; // outgoing phase
              double Glim = sec_m(e, n_lim) / vol + p_ar(e, n_lim) / dt; //maximum acceptable phase change rate for this phase
              if (std::fabs(G) < Glim) n_lim = -1;       //G does not make the phase evanescent -> no limiting (n_lim = -2)
              else G = (G > 0 ? 1 : -1) * Glim;//the phase would become evanescent due to G -> clamp to G_lim (n_lim = k / l)

              double hk = G > 0 ? h(e, k) : Hls_tab(e,i_sat), dTk_hk = G > 0 && dT_h ? (*dT_h)(e, k) : 0, dP_hk = G > 0 ? (dP_h ? (*dP_h)(e, k) : 0) : dPHls_tab(e,i_sat),
                       hl = G < 0 ? h(e, l) : Hvs_tab(e,i_sat), dTl_hl = G < 0 && dT_h ? (*dT_h)(e, l) : 0, dP_hl = G < 0 ? (dP_h ? (*dP_h)(e, l) : 0) : dPHvs_tab(e,i_sat);
              if (n_lim < 0 && is_therm) /* derivatives of G at the thermal limit */
              {
                  double dP_phi = dP_hi(k, l) * (Tk - Ts) + dP_hi(l, k) * (Tl - Ts) - (hi(k, l) + hi(l, k)) * dP_Ts + dpqi(e, k, l) / vol,
                         dTk_L = -dTk_hk, dTl_L = dTl_hl, dP_L = dP_hl - dP_hk;
                  dP_G = (dP_phi - G * dP_L) / L;
                  for (n = 0; n < N; n++) dT_G(n) = ((n == k) * (hi(k, l) - G * dTk_L) + (n == l) * (hi(l, k) - G * dTl_L) + dT_hi(k, l, n) * (Tk - Ts) + dT_hi(l, k, n) * (Tl - Ts) + dTqi(e, k, l, n)) / L;
                  for (n = 0; n < N; n++) da_G(n) = (da_hi(k, l, n) * (Tk - Ts) + da_hi(l, k, n) * (Tl - Ts) + daqi(e, k, l, n)) / L;
                }


              if (sub_type(Masse_Multiphase, equation())) //mass equation -> phase change
                {
                  for (i = 0; i < 2; i++) secmem(e, i ? l : k) -= vol * (i ? -1 : 1) * G;
                  if (n_lim < 0) /* G by thermal limit */
                    {
                      if (Ma)
                        for (i = 0; i < 2; i++)
                          for (n = 0; n < N; n++)  // derivatives w.r.t. alpha
                            (*Ma)(N * e + (i ? l : k), N * e + n) += vol * (i ? -1 : 1) * da_G(n);
                      if (Mt)
                        for (i = 0; i < 2; i++)
                          for (n = 0; n < N; n++)  // derivatives w.r.t. T
                            (*Mt)(N * e + (i ? l : k), N * e + n) += vol * (i ? -1 : 1) * dT_G(n);
                      if (Mp)
                        for (i = 0; i < 2; i++) // derivatives w.r.t. p
                          (*Mp)(N * e + (i ? l : k), e) += vol * (i ? -1 : 1) * dP_G;
                    }
                  else for (auto &s_d : vec_m) /* G by evanescence */
                      for (auto j = s_d[0]->get_tab1()(N * e + n_lim) - 1; j < s_d[0]->get_tab1()(N * e + n_lim + 1) - 1; j++)
                        for (col = s_d[0]->get_tab2()(j) - 1, x = -s_d[0]->get_coeff()(j), i = 0; i < 2; i++)
                          (*s_d[1])(N * e + (i ? l : k), col) += (i ? -1 : 1) * sgn * x;
                }
              else if (sub_type(Energie_Multiphase, equation())) // energy equation -> heat transfer
                {
                  // thermal limit applies on one side: c (=0,1) / n_c (=k,l) / sign of outgoing flux from phase k: s_c (=1,-1)
                  int c = (a_r(e, k) > a_r(e, l)), n_c = c ? l : k, n_d = c ? k : l, s_c = c ? -1 : 1;
                  double Tc = c ? Tl : Tk, hc = c ? hl : hk, dT_hc = c ? dTl_hl : dTk_hk, dP_hc = c ? dP_hl : dP_hk;
                  for (i = 0; i < 2; i++) secmem(e, i ? l : k)-= vol * (i ? -1 : 1) * (s_c *        hi(n_c, n_d) * (Tc - Ts)                              + G * hc)                                              - (i != c) * qi(e, k, l);
                  /* derivatives (including those of G, except in the limit case) */
                  if (Ma)
                    for (i = 0; i < 2; i++)
                      for (n = 0; n < N; n++) // derivatives w.r.t. alpha
                        (*Ma)(N * e + (i ? l : k), N * e + n) += vol * (i ? -1 : 1) * (s_c *  da_hi(n_c, n_d, n) * (Tc - Ts)                              + (n_lim < 0) * da_G(n) * hc)                          - (i != c) * daqi(e, k, l, n);
                  if (Mt)
                    for (i = 0; i < 2; i++)
                      for (n = 0; n < N; n++) // derivatives w.r.t. T
                        (*Mt)(N * e + (i ? l : k), N * e + n) += vol * (i ? -1 : 1) * (s_c * (dT_hi(n_c, n_d, n) * (Tc - Ts) + (n == n_c) * hi(n_c, n_d)) + (n_lim < 0) * dT_G(n) * hc + G * (n == n_c) * dT_hc) - (i != c) * dTqi(e, k, l, n);
                  if (Mp)
                    for (i = 0; i < 2; i++) // derivatives w.r.t. p
                      (*Mp)(N * e + (i ? l : k), e)           += vol * (i ? -1 : 1) * (s_c * (dP_hi(n_c, n_d)    * (Tc - Ts) - hi(n_c, n_d) * dP_Ts)      + (n_lim < 0) * dP_G * hc + G * dP_hc)                 - (i != c) * dpqi(e, k, l);
                  if (n_lim >= 0)
                    for (auto &s_d : vec_m) /* derivatives of G in the evanescent case */
                      for (auto j = s_d[0]->get_tab1()(N * e + n_lim) - 1; j < s_d[0]->get_tab1()(N * e + n_lim + 1) - 1; j++)
                        for (col = s_d[0]->get_tab2()(j) - 1, x = -s_d[0]->get_coeff()(j), i = 0; i < 2; i++)
                          (*s_d[1])(N * e + (i ? l : k), col) += (i ? -1 : 1) * hc * sgn * x;
                }
              else if (sub_type(Aire_interfaciale, equation())) //interfacial area equation; looks like the mass equation; not a conservation equation!
                if (0.> mod2grp)
                  {
                    if (k==0) // k is the carrier phase
                      if (alpha(e, l) > alpha_min) // if the phase l is present
                        {
                          secmem(e, l) += vol * 2./3. * inco(e, l) / (alpha(e, l) * (pch_rho ? (*pch_rho).valeurs()(e, l) : rho(e, l))) * G ;
                          if (n_lim < 0) /* G by thermal limit */
                            {
                              if (Ma)   // derivatives w.r.t. alpha
                                {
                                  (*Ma)(N * e + l , N * e + l) -= vol * 2./3. * inco(e, l) / (-alpha(e, l)*alpha(e, l) * (pch_rho ? (*pch_rho).valeurs()(e, l) : rho(e, l))) * G;
                                  for (n = 0; n < N; n++) (*Ma)(N * e + l , N * e + n) -=
                                      vol * 2./3. * inco(e, l) / (alpha(e, l)              * (pch_rho ? (*pch_rho).valeurs()(e, l) : rho(e, l))) * da_G(n);
                                }
                              if (Mt)   // derivatives w.r.t. T
                                {
                                  if (pch_rho) (*Mt)(N * e + l , N * e + l) -=
                                      vol * 2./3. * inco(e, l) / alpha(e, l)*-pch_rho->derivees().at("temperature")(e, l)/((*pch_rho).valeurs()(e, l)*(*pch_rho).valeurs()(e, l)) * G;
                                  for (n = 0; n < N; n++) (*Mt)(N * e + l , N * e + n) -=
                                      vol * 2./3. * inco(e, l) / (alpha(e, l) * (pch_rho ? (*pch_rho).valeurs()(e, l) : rho(e, l))) * dT_G(n);
                                }
                              if (Mp)  // derivatives w.r.t. p
                                {
                                  if (pch_rho) (*Mp)(N * e + l , e) -=
                                      vol * 2./3. * inco(e, l) / alpha(e, l)*-pch_rho->derivees().at("pression")(e, l)/((*pch_rho).valeurs()(e, l)*(*pch_rho).valeurs()(e, l)) * G;
                                  for (n = 0; n < N; n++) (*Mp)(N * e + l , N * e + n) -=
                                      vol * 2./3. * inco(e, l) / (alpha(e, l) * (pch_rho ? (*pch_rho).valeurs()(e, l) : rho(e, l))) * dP_G;
                                }
                              if (Mai) // derivatives w.r.t. ai
                                {
                                  (*Mai)(N * e + l , N * e + l) -= vol * 2./3. / (alpha(e, l)* (pch_rho ? (*pch_rho).valeurs()(e, l) : rho(e, l))) * G;
                                } // dAi_G to add
                            }
                          else for (auto &s_d : vec_m) /* G by evanescence */
                              for (auto j = s_d[0]->get_tab1()(N * e + n_lim) - 1; j < s_d[0]->get_tab1()(N * e + n_lim + 1) - 1; j++)
                                for (col = s_d[0]->get_tab2()(j) - 1, x = -s_d[0]->get_coeff()(j), i = 0; i < 2; i++)
                                  (*s_d[1])(N * e + l , col) -= 2./3. * inco(e, l) / (alpha(e, l) * (pch_rho ? (*pch_rho).valeurs()(e, l) : rho(e, l))) * sgn * x;
                        }
                  }
            }
          else if (sub_type(Energie_Multiphase, equation())) /* no saturation: energy exchanges only */
            {
              // outgoing flux from phase k: hm * (Tk - Tl)
              double hm = 1. / (1. / hi(k, l) + 1. / hi(l, k)), Tk = temp(e, k), Tl = temp(e, l);
              for (i = 0; i < 2; i++) secmem(e, i ? l : k) -= vol * (i ? -1 : 1) * hm * (Tk - Tl);
              if (Mt)
                for (i = 0; i < 2; i++) // only T derivatives
                  {
                    (*Mt)(N * e + (i ? l : k), N * e + k) += vol * (i ? -1 : 1) * hm;
                    (*Mt)(N * e + (i ? l : k), N * e + l) += vol * (i ? 1 : -1) * hm;
                  }
            }
    }
}
