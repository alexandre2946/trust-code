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

#include <Discretisation_base.h>
#include <Operateur_Diff_base.h>
#include <Schema_Temps_base.h>
#include <Aire_interfaciale.h>
#include <Masse_Multiphase.h>
#include <Neumann_val_ext.h>
#include <MD_Vector_tools.h>
#include <Assembleur_base.h>
#include <TRUSTTab_parts.h>
#include <QDM_Multiphase.h>
#include <Pb_Multiphase.h>
#include <Pb_Conduction.h>
#include <MD_Vector_std.h>
#include <Matrix_tools.h>
#include <Matrice_Bloc.h>
#include <Array_tools.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>
#include <Dirichlet.h>
#include <SFichier.h>
#include <EChaine.h>
#include <Lapack.h>
#include <Debog.h>
#include <SETS.h>

#include <TablePrinter.h>
using bprinter::TablePrinter;

#ifndef NDEBUG
template class std::map<std::string, matrices_t>;
template class std::map<std::string, Matrice_Morse>;
#endif

Implemente_instanciable_sans_constructeur(SETS, "SETS", Simpler);
// XD sets simpler sets INHERITS_BRACE Stability-Enhancing Two-Step solver which is useful for a multiphase problem. Ref
// XD_CONT : J. H. MAHAFFY, A stability-enhancing two-step method for fluid flow calculations, Journal of Computational
// XD_CONT Physics, 46, 3, 329 (1982).
// XD attr criteres_convergence bloc_criteres_convergence criteres_convergence OPT Set the convergence thresholds for
// XD_CONT each unknown (i.e: alpha, temperature, velocity and pressure). The default values are respectively 0.01, 0.1,
// XD_CONT 0.01 and 100
// XD attr iter_min entier iter_min OPT Number of minimum iterations (default value 1)
// XD attr iter_max entier iter_max OPT Number of maximum iterations (default value 10)
// XD attr seuil_convergence_implicite floattant seuil_convergence_implicite OPT Convergence criteria.
// XD attr nb_corrections_max entier nb_corrections_max OPT Maximum number of corrections performed by the PISO
// XD_CONT algorithm to achieve the projection of the velocity field. The algorithm may perform less corrections then
// XD_CONT nb_corrections_max if the accuracy of the projection is sufficient. (By default nb_corrections_max is set to
// XD_CONT 21).
// XD attr facsec_diffusion_for_sets floattant facsec_diffusion_for_sets OPT facsec to impose on the diffusion time step
// XD_CONT in sets while the total time step stays smaller than the convection time step.
// XD bloc_criteres_convergence bloc_lecture nul NO_BRACE Not set

Implemente_instanciable_sans_constructeur(ICE, "ICE", SETS);
// XD ice sets ice INHERITS_BRACE Implicit Continuous-fluid Eulerian solver which is useful for a multiphase problem.
// XD_CONT Robust pressure reduction resolution.
// XD attr pression_degeneree entier pression_degeneree OPT Set to 1 if the pressure field is degenerate (ex. :
// XD_CONT incompressible fluid with no imposed-pressure BCs). Default: autodetected
// XD attr reduction_pression|pressure_reduction entier reduction_pression OPT Set to 1 if the user wants a resolution
// XD_CONT with a pressure reduction. Otherwise, the flag is to be set to 0 so that the complete matrix is considered.
// XD_CONT The default value of this flag is 1.
// XD criteres_convergence interprete nul NO_BRACE convergence criteria
// XD attr aco chaine(into=["{"]) aco REQ Opening curly bracket.
// XD attr inco chaine inco OPT Unknown (i.e: alpha, temperature, velocity and pressure)
// XD attr val floattant val OPT Convergence threshold
// XD attr acof chaine(into=["}"]) acof REQ Closing curly bracket.

SETS::SETS()
{
  sets_ = 1;
}

ICE::ICE()
{
  sets_ = 0;
}

Sortie& SETS::printOn(Sortie& os) const
{
  return Simpler::printOn(os);
}

Entree& SETS::readOn(Entree& is)
{
  /* default values of the convergence criteria */
  crit_conv_ = {{ "alpha", 1e-2 }, { "temperature", 1e-1 },
    { "enthalpie", 1e2 }, { "vitesse", 1e-2 },
    { "pression", 100 }, { "k", 1e-2 },
    { "tau", 1e-2 }, { "omega", 1e-2 },
    { "k_WIT", 1e-2 }, { "interfacial_area", 1e2 }
  };
  return Simpler::readOn(is);
}

Entree& SETS::lire(const Motcle& mot, Entree& is)
{
  if (mot == Motcle("criteres_convergence"))
    {
      Nom nom;
      is >> nom;
      if (nom != "{")
        Process::exit(Nom("SETS::lire() : { expected instead of ") + nom);

      for (is >> nom; nom != "}"; is >> nom)
        {
          double val;
          is >> val;
          crit_conv_[nom.getString()] = val;
        }
    }
  else if (mot == "iter_min")
    is >> iter_min_;
  else if (mot == "iter_max")
    is >> iter_max_;
  else if (mot == "pression_degeneree")
    is >> p_degen_;
  else if (mot == "pressure_reduction" || mot == "reduction_pression")
    is >> pressure_reduction_;
  else
    return Simpler::lire(mot, is); // does the parent class know this keyword?
  return is;
}

Sortie& ICE::printOn(Sortie& os) const
{
  return SETS::printOn(os);
}

Entree& ICE::readOn(Entree& is)
{
  return SETS::readOn(is);
}

static inline double corriger_incr_alpha(const DoubleTab& alpha, DoubleTab& incr, double& err_a_sum)
{
  const int N = alpha.line_size();
  double a_sum, corr_max = 0;
  DoubleTrav n_a(N);

  err_a_sum = 0;

  for (int i = 0; i < alpha.dimension_tot(0); i++)
    {
      a_sum = 0;
      for (int n = 0; n < N; n++)
        {
          n_a(n) = alpha(i, n) + incr(i, n);
          a_sum += n_a(n);
        }
      err_a_sum = std::max(err_a_sum, std::abs(a_sum - 1));

      a_sum = 0;
      for (int n = 0; n < N; n++)
        {
          n_a(n) = std::max(n_a(n), 0.);
          a_sum += n_a(n);
        }
      if (static_cast<bool>(a_sum))
        for (int n = 0; n < N; n++)
          n_a(n) /= a_sum;
      else
        for (int n = 0; n < N; n++)
          n_a(n) = 1. / N;

      for (int n = 0; n < N; n++)
        {
          corr_max = std::max(corr_max,
                              std::fabs(alpha(i, n) + incr(i, n) - n_a(n)));
          incr(i, n) = n_a(n) - alpha(i, n);
        }
    }
  err_a_sum = Process::mp_max(err_a_sum);
  return Process::mp_max(corr_max);
}

double SETS::unknown_positivation(const DoubleTab& uk, DoubleTab& incr)
{
  double corr_max = 0;
  const int N = uk.line_size();
  for (int i = 0; i < uk.dimension_tot(0); i++)
    for (int n = 0; n < N; n++)
      if (uk(i, n) + incr(i, n) < 0)
        {
          if (-uk(i, n) - incr(i, n) > corr_max)
            corr_max = std::abs(uk(i, n) + incr(i, n));
          incr(i, n) = -uk(i, n);
        }
  return Process::mp_max(corr_max);
}

#ifdef PETSCKSP_H
void SETS::init_cv_ctx(const DoubleTab& secmem, const DoubleVect& norme)
{
  cv_ctx = (SETS::cv_test_t *) calloc(1, sizeof(SETS::cv_test_t));
  cv_ctx->obj = this;
  cv_ctx->eps_alpha = crit_conv_["alpha"];
  norm = norme;
  residu = secmem;
  cv_ctx->t = nullptr;
  cv_ctx->v = nullptr;
  /* numbering to retrieve the residual: same approach as in Solv_Petsc */
  ArrOfBit items_to_keep;
  const int size = secmem.size_array();
  auto idx = mppartial_sum(secmem.get_md_vector()->get_sequential_items_flags(items_to_keep, secmem.line_size()));
  ix.resize(size);
  for (int i = 0; i < size; i++)
    if (items_to_keep[i])
      ix[i] = idx, idx++;
    else
      ix[i] = -1;
  KSPConvergedDefaultCreate(&cv_ctx->defctx);
}

#if PETSC_VERSION_GE(3,24,0)
PetscErrorCode SETS::destroy_cvctx(void **mctx)
{
  SETS::cv_test_t *ctx = (SETS::cv_test_t *)*mctx;
  if (ctx->v)
    VecDestroy(&ctx->v);
  if (ctx->t)
    VecDestroy(&ctx->t);
  PetscErrorCode err = KSPConvergedDefaultDestroy(&ctx->defctx);
  free(ctx);
  return err;
}
#else
PetscErrorCode SETS::destroy_cvctx(void *mctx)
{
  SETS::cv_test_t *ctx = (SETS::cv_test_t *)mctx;
  if (ctx->v)
    VecDestroy(&ctx->v);
  if (ctx->t)
    VecDestroy(&ctx->t);
  PetscErrorCode err =  KSPConvergedDefaultDestroy(ctx->defctx);
  free(ctx);
  return err;
}
#endif

/* test de convergence */
PetscErrorCode SETS::convergence_test(KSP ksp, PetscInt it, PetscReal rnorm, KSPConvergedReason *reason,void *mctx)
{
  SETS::cv_test_t *ctx = (SETS::cv_test_t *)mctx;
  if (ctx->t == nullptr) /* ctx->t, ctx-v non initialises -> on les cree */
    {
      Mat m;
      KSPGetOperators(ksp,&m,nullptr);
      MatCreateVecs(m, &ctx->v, &ctx->t);
      VecSetOption(ctx->v, VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    }
  // PETSc 3.20 -> 3.22 change bug ? KSPBuildResidual corrupted after KSPConvergedDefault call
  // So we switch order, annoying cause KSPBuildResidual is expensive...
  Vec resi;
  KSPBuildResidual(ksp, ctx->t, ctx->v, &resi);//residu
  //Cerr << "Residual before:" << finl;
  //VecView(resi,PETSC_VIEWER_STDOUT_WORLD);
  PetscErrorCode ret = KSPConvergedDefault(ksp, it, rnorm, reason, &ctx->defctx); //call the default convergence test
  if (ret || *reason <= 0) return ret; //not yet converged -> nothing to do
  /* sinon -> on verfie que sum alpha = 1 est aussi OK */
  //KSPBuildResidual(ksp, ctx->t, ctx->v, &resi);//residu
  //Cerr << "Residual after:" << finl;
  //VecView(resi,PETSC_VIEWER_STDOUT_WORLD);
  VecGetValues(resi, ctx->obj->ix.size_array(), ctx->obj->ix.addr(), ctx->obj->residu.addr());
  bool ok = true;
  for (int i = 0; ok && i < ctx->obj->norm.size(); i++)
    if (ctx->obj->ix[i] >= 0)
      ok &= std::abs(ctx->obj->residu[i]) < ctx->obj->norm[i] * ctx->eps_alpha;
  if (!Process::mp_and(ok)) *reason = KSP_CONVERGED_ITERATING;
  return ret;
}
#endif

bool SETS::iterer_eqn(Equation_base& eqn,
                      const DoubleTab& inut,
                      DoubleTab& current,
                      double dt,
                      int numero_iteration,
                      int& ok)
{
  /* only Pb_Multiphase or Pb_conduction are handled */
  if (!sub_type(Pb_Multiphase, eqn.probleme()) && !sub_type(Pb_Conduction, eqn.probleme()))
    Process::exit(que_suis_je() + " cannot be applied to the problem "
                  + eqn.probleme().le_nom() + " of type " + eqn.probleme().que_suis_je() + "!");

  /* equations not solved directly: Masse_Multiphase and Aire_interfaciale (always), Energie_Multiphase (in ICE), Convection_Diffusion_std i.e. turbulent quantities (in ICE) */
  if (sub_type(Masse_Multiphase, eqn) || sub_type(Aire_interfaciale, eqn) || (!sets_ && sub_type(Energie_Multiphase, eqn))|| (!sets_ && sub_type(Convection_Diffusion_std, eqn)))
    {
      if (eqn.positive_unkown())
        {
          DoubleTrav incr;
          incr = current;
          incr -= eqn.inconnue().valeurs();
          unknown_positivation(eqn.inconnue().valeurs(), incr);
          incr += eqn.inconnue().valeurs();
          current = incr;
        }
      return true;
    }

  /* QDM_Multiphase: solved by iterer_NS */
  int cv {0};
  if (sub_type(QDM_Multiphase, eqn))
    {
      Matrice_Morse matrix_unused;
      DoubleTrav secmem_unused;
      iterer_NS(eqn, current, ref_cast(QDM_Multiphase, eqn).pression().valeurs(),
                dt, matrix_unused, 0, secmem_unused, numero_iteration, cv, ok);
      return cv;
    }

  /* remaining case: thermal or turbulent equation of a Pb_Multi or Pb_conduction -> set semi_impl if needed, then solve */
  const std::string& nom_inco = eqn.inconnue().le_nom().getString();
  const std::string nom_pb_inco = eqn.probleme().le_nom().getString() + "/" + nom_inco;

  tabs_t semi_impl; /* in ICE, temperatures of all problems are explicit */
  const Operateur_Diff_base& op_diff = ref_cast(Operateur_Diff_base, eqn.operateur(0).l_op_base());

  if (!sets_)
    for (auto &&op_ext : op_diff.op_ext)
      semi_impl[nom_inco + (op_ext != &op_diff ? "/" + op_ext->equation().probleme().le_nom().getString() : "")] = op_ext->equation().inconnue().passe();

  if (!mat_pred_.count(nom_pb_inco)) /* matrix: size allocated on first pass */
    eqn.dimensionner_blocs( { { nom_inco, &mat_pred_[nom_pb_inco] } }, semi_impl);

  /* assembly and resolution */
  DoubleTrav secmem(current);
  SolveurSys& solv = get_and_set_parametre_implicite(eqn).solveur();
  eqn.assembler_blocs_avec_inertie( { { nom_inco, &mat_pred_[nom_pb_inco] } }, secmem, semi_impl);
  mat_pred_[nom_pb_inco].ajouter_multvect(current, secmem); //convert increment -> variable
  solv->reinit();
  solv.resoudre_systeme(mat_pred_[nom_pb_inco], secmem, current);

  if (eqn.positive_unkown())
    {
      DoubleTrav incr;
      incr = current;
      incr -= eqn.inconnue().valeurs();
      unknown_positivation(eqn.inconnue().valeurs(), incr);
      incr += eqn.inconnue().valeurs();
      current = incr;
    }

  /* update */
  eqn.probleme().mettre_a_jour(eqn.schema_temps().temps_courant());
  eqn.inconnue().futur() = current;
  return true;
}

//Input: Un ; Pn
//Output: Un+1 = U***_k ; Pn+1 = P**_k
//n denotes a time step

void SETS::iterer_NS(Equation_base& eqn, DoubleTab& current,
                     DoubleTab& pression, double dt,
                     Matrice_Morse& matrice_unused, double seuil_resol,
                     DoubleTrav& secmem_unused, int nb_ite, int& cv, int& ok)
{
  Pb_Multiphase& pb = *(Pb_Multiphase*) &ref_cast(Pb_Multiphase, eqn.probleme());
  const bool res_en_T = pb.resolution_en_T();

  const int n_eq = pb.nombre_d_equations();
  int& it = iteration_;
  QDM_Multiphase& eq_qdm = ref_cast(QDM_Multiphase, eqn);
  double err_a_sum = 0;

  std::vector<Equation_base*> eq_list; // equation order
  for (int i = 0; i < n_eq; i++)
    eq_list.push_back(&pb.equation(i));
  std::map<std::string, Equation_base*> eqs; // eqs[unknown] = equation
  std::vector<std::string> noms; // unknown order: same as equations, then pressure
  for (int i = 0; i < n_eq; i++)
    {
      noms.push_back(eq_list[i]->inconnue().le_nom().getString());
      eqs[noms[i]] = eq_list[i];
    }
  noms.push_back("pression"); // no equation associated with pressure!

  std::map<std::string, Champ_Inc_base*> inco; // all Champ_Inc
  for (auto &i_eq : eqs)
    inco[i_eq.first] = &i_eq.second->inconnue();
  inco["pression"] = &eq_qdm.pression();
  // Newton initialization with current values (stored in passe() in implicit mode), sizing incr / sec
  double t = eqn.schema_temps().temps_courant();
  pb.mettre_a_jour(t); //unknowns -> medium -> conserved fields

  /* valeurs semi-implicites : inconnues (alpha, v, T, ..) et champs conserves (alpha_rho, alpha_rho_e, ..) */
  tabs_t semi_impl;
  for (auto &&i_eq : eqs)
    if (i_eq.second != &eq_qdm)
      {
        semi_impl[i_eq.first] = i_eq.second->inconnue().passe();
        if (i_eq.second->has_champ_conserve())
          semi_impl[i_eq.second->champ_conserve().le_nom().getString()] = i_eq.second->champ_conserve().passe();
        if (i_eq.second->has_champ_convecte())
          semi_impl[i_eq.second->champ_convecte().le_nom().getString()] = i_eq.second->champ_convecte().passe();
      }
  // in SETS, replace the past velocity value with the one from a prediction step
  if (sets_ && !first_call_)
    {
      DoubleTrav secmem(current);
      /* implicit assembly, velocities only */
      if (!mat_pred_.count("vitesse"))
        eq_qdm.dimensionner_blocs( { { "vitesse", &mat_pred_["vitesse"] } }); // first pass: sizing
      eq_qdm.assembler_blocs_avec_inertie( { { "vitesse", &mat_pred_["vitesse"] } }, secmem);

      /* solve and store the predicted velocity in current */
      SolveurSys& solv_qdm = get_and_set_parametre_implicite(eqn).solveur();
      solv_qdm->reinit();
      mat_pred_["vitesse"].ajouter_multvect(current, secmem); //convert increment -> variable to satisfy iterative solvers
      solv_qdm.resoudre_systeme(mat_pred_["vitesse"], secmem, current);
      semi_impl["vitesse"] = current;

      if (facsec_diffusion_for_sets() > 0)
        for (auto &&i_eq : eqs)
          if (i_eq.second != &eq_qdm)
            {
              Cerr << "Prediction of " << i_eq.first << finl;
              const DoubleTab inut_loc(0);
              iterer_eqn(*i_eq.second, inut_loc, i_eq.second->inconnue().valeurs(), dt, 0, ok);
              /*
               semi_impl[i_eq.first] = (*i_eq.second).inconnue().valeurs();
               if (i_eq.second->has_champ_conserve()) semi_impl[i_eq.second->champ_conserve().le_nom().getString()] = i_eq.second->champ_conserve().passe();
               if (i_eq.second->has_champ_convecte()) semi_impl[i_eq.second->champ_convecte().le_nom().getString()] = i_eq.second->champ_convecte().passe();
               */
            }
    }
  else
    semi_impl["vitesse"] = eq_qdm.inconnue().passe();
  eqn.solv_masse().corriger_solution(current, current, 0); //for PolyMAC_MPFA: vf -> ve

  // first pass: sizing mat_semi_impl and mdv_semi_impl, filling p_degen_
  if (!mat_semi_impl_.nb_lignes())
    {
      mat_semi_impl_.dimensionner(n_eq + 1, n_eq + 1); //last row -> continuity
      for (int i = 0; i <= n_eq; i++)
        for (int j = 0; j <= n_eq; j++)
          mat_semi_impl_.get_bloc(i, j).typer("Matrice_Morse"), mats_[noms[i]][noms[j]] = &ref_cast(Matrice_Morse, mat_semi_impl_.get_bloc(i, j).valeur());
      for (auto &&i_eq : eqs)
        i_eq.second->dimensionner_blocs(mats_[i_eq.first], semi_impl); //semi-implicit option
      eq_qdm.assembleur_pression()->dimensionner_continuite(mats_["pression"]);
      /* if mat_semi_impl is used directly: size the empty blocks */
      if (!pressure_reduction_)
        for (int i = 0; i <= n_eq; i++)
          for (int j = 0; j <= n_eq; j++)
            if (!mats_[noms[i]][noms[j]]->nb_lignes() && !mats_[noms[i]][noms[j]]->nb_colonnes())
              mats_[noms[i]][noms[j]]->dimensionner(inco[noms[i]]->valeurs().size_totale(), inco[noms[j]]->valeurs().size_totale(), 0);

      MD_Vector_composite mdc;
      for (int i = 0; i <= n_eq; i++)
        mdc.add_part(inco[noms[i]]->valeurs().get_md_vector(), inco[noms[i]]->valeurs().line_size());
      mdv_semi_impl_.copy(mdc);

      /* set p_degen if not read: for incompressible flow without imposed pressure BCs, pressure is degenerate */
      if (p_degen_ < 0)
        {
          p_degen_ = sub_type(Fluide_base, eq_qdm.milieu());
          for (int i = 0; i < eq_qdm.domaine_Cl_dis().nb_cond_lim(); i++)
            p_degen_ &= !sub_type(Neumann_val_ext, eq_qdm.domaine_Cl_dis().les_conditions_limites(i).valeur());
        }
    }

  /* increments and residuals for the Newton algorithm */
  DoubleTrav v_inco, v_incr, v_sec; //format "complete vector"
  MD_Vector_tools::creer_tableau_distribue(mdv_semi_impl_, v_incr);
  MD_Vector_tools::creer_tableau_distribue(mdv_semi_impl_, v_sec);

  DoubleTab_parts p_incr(v_incr);
  DoubleTab_parts p_sec(v_sec);
  ptabs_t incr, sec; //format "std::map"
  for (int i = 0; i <= n_eq; i++)
    {
      incr[noms[i]] = &p_incr[i];
      sec[noms[i]] = &p_sec[i];
    }
  if (!pressure_reduction_) /* v_inco: filled by copying unknowns when no pressure reduction is performed */
    {
      MD_Vector_tools::creer_tableau_distribue(mdv_semi_impl_, v_inco);
      DoubleTab_parts p_inco(v_inco);
      for (int i = 0; i <= n_eq; i++)
        p_inco[i] = inco[noms[i]]->valeurs();
    }
  /* Newton: assembly of mat_semi_impl -> assembly of the pressure matrix -> solve -> back-substitution */
  TablePrinter tp(&get_Cerr().get_ostream());
  if (!Process::me())
    {
      tp.AddColumn("it", 5);
      Nom fichier(nom_du_cas() + ".newton_evol");
      if (!header_written_)
        {
          Cerr << "SETS => File " << fichier << " is created ..." << finl;
          SFichier newton_evol(fichier);
          newton_evol << "# File showing the simulation time, time step, Newton iterations, status and convergence of the increments." << finl;
          newton_evol << "# Time \t Time_step \t Newton \t Status \t ";
        }

      for (auto &&n_i : inco)
        {
          tp.AddColumn(n_i.first, std::max(12, (int) n_i.first.length()));
          if (!header_written_)
            {
              SFichier newton_evol(fichier, ios::app);
              newton_evol << n_i.first << "\t ";
            }
        }
      header_written_ = true;
      tp.PrintHeader();
    }

  cv = 0;
  for (it = 0; it < iter_min_ || (!cv && it < iter_max_); it++)
    {
      /* fill via assembler_blocs */
      // energy equation first, so that q_pi can be used in other equations
      res_en_T ?
      eqs["temperature"]->assembler_blocs_avec_inertie(mats_["temperature"], *sec["temperature"], semi_impl) :
      eqs["enthalpie"]->assembler_blocs_avec_inertie(mats_["enthalpie"], *sec["enthalpie"], semi_impl);
      // the others
      for (auto &n_e : eqs)
        if (n_e.first != "temperature" && n_e.first != "enthalpie")
          n_e.second->assembler_blocs_avec_inertie(mats_[n_e.first], *sec[n_e.first], semi_impl);

      // continuity equation sum_k alpha_k = 1
      eq_qdm.assembleur_pression()->assembler_continuite(mats_["pression"], *sec["pression"]);

      SolveurSys& solv_p = eq_qdm.solveur_pression(); /* use the "solveur_pression" of the QDM equation */
      if (pressure_reduction_) /* pressure reduction */
        {
          /* express the other unknowns (x) as functions of p: velocity, then temperature / pressure */
          tabs_t b_p;
          std::vector<std::set<std::pair<std::string, int>>> ordre;
          if (eq_qdm.domaine_dis().le_nom() == "PolyMAC_MPFA")
            ordre.push_back( { { "vitesse", 1 } }); //if PolyMAC_MPFA: start with ve
          ordre.push_back( { { "vitesse", 0 } }), ordre.push_back( { }); //then vf, then all other unknowns simultaneously
          for (auto &&nom : noms)
            if (nom != "vitesse" && nom != "pression")
              ordre.back().insert( { { nom, 0 } });
          ok = mp_min(eliminer(ordre, "pression", mats_, sec, A_p_, b_p));
          if (!ok)
            {
              Cerr << "Elimination failed!";
              break; //if the elimination fails, we exit
            }

          /* assembly of the pressure system */
          assembler("pression", A_p_, b_p, mats_, sec, matrice_pression_, *sec["pression"], p_degen_);
#ifdef PETSCKSP_H
          if (!cv_ctx && sub_type(Solv_Petsc, solv_p.valeur()))
            {
              init_cv_ctx(*sec["pression"], eq_qdm.assembleur_pression()->norme_continuite());
              ref_cast(Solv_Petsc, solv_p.valeur()).set_convergence_test(convergence_test, cv_ctx, destroy_cvctx);
            }
#endif

          /* solve: only if the error on alpha (in secmem_pression) exceeds a threshold */
          if (mp_max_abs_vect(*sec["pression"]) > 1e-16)
            {
              matrice_pression_.ajouter_multvect(inco["pression"]->valeurs(), *sec["pression"]); //convert increment -> variable to satisfy iterative solvers
              solv_p->reinit(), solv_p->set_return_on_error(1); /* to avoid an exit() on failure */
              ok = (solv_p.resoudre_systeme(matrice_pression_, *sec["pression"], *incr["pression"]) >= 0);
              if (!ok)
                {
                  Cerr << "Solver failed!";
                  break; //solver failed -> exit
                }
              incr["pression"]->echange_espace_virtuel();
              *incr["pression"] -= inco["pression"]->valeurs();
            }
          else
            *incr["pression"] = 0;

          // increments of the other variables
          for (auto &&n_v : b_p)
            {
              *incr[n_v.first] = n_v.second; //constant part
              A_p_[n_v.first].ajouter_multvect(*incr["pression"], *incr[n_v.first]); //dependence on pressure increments
              incr[n_v.first]->echange_espace_virtuel();
            }
        }
      else /* no pressure reduction: go directly through mat_semi_impl */
        {
          mat_semi_impl_.ajouter_multvect(v_inco, v_sec); //convert increment -> variable to satisfy iterative solvers
          solv_p->reinit(), solv_p->set_return_on_error(1); /* to avoid an exit() on failure */
          ok = (solv_p.resoudre_systeme(mat_semi_impl_, v_sec, v_incr) >= 0);
          if (!ok)
            {
              Cerr << "Solver failed!";
              break; //solver failed -> exit
            }
          v_incr -= v_inco; //back to increments
        }

      eqn.solv_masse().corriger_solution(*incr["vitesse"], *incr["vitesse"], 1); //for PolyMAC_MPFA: used to correct ve

      if (!Process::me())
        tp << it + 1;
      incr_var_convergence_.clear(); // clear if new time-step or Newton incr !
      for (auto &&n_v : incr)
        {
          const double vm = mp_min_vect(*n_v.second);
          const double vM = mp_max_vect(*n_v.second);
          const double x = std::fabs(vM) > std::fabs(vm) ? vM : vm;
          if (!Process::me())
            {
              tp << x;
              incr_var_convergence_.push_back(x);
            }
        }

      /* convergence check */
      cv = corriger_incr_alpha(inco["alpha"]->valeurs(), *incr["alpha"], err_a_sum) < crit_conv_["alpha"];
      for (auto &&n_v : incr)
        if (crit_conv_.count(n_v.first))
          cv &= mp_max_abs_vect(*n_v.second) < crit_conv_.at(n_v.first);

      /* updates: unknowns -> medium -> conserved fields -> sources */
      for (auto &&n_i : inco)
        n_i.second->valeurs() += *incr[n_i.first];
      if (p_degen_)
        inco["pression"]->valeurs() -= mp_min_vect(inco["pression"]->valeurs()); // Use the minimum pressure as reference so that the same reference is used in sequential and parallel
      if (!(ok = err_a_sum < crit_conv_["alpha"]))
        {
          Cerr << "Error in alpha!";
          break; //if we have exceeded the medium bounds on (p, T) or lack precision, we must exit
        }
      if (!(ok = eq_qdm.milieu().check_unknown_range()))
        {
          Cerr << "Out of bounds!";
          break; //if we have exceeded the medium bounds on (p, T) or lack precision, we must exit
        }
      pb.mettre_a_jour(t); //unknowns -> medium -> conserved fields
    }

  if (!Process::me())
    {
      tp.PrintFooter();
      Nom fichier(nom_du_cas() + ".newton_evol");
      SFichier newton_evol(fichier, ios::app);
      newton_evol.setf(ios::scientific);
      newton_evol << finl << t << "\t " << eqn.schema_temps().pas_de_temps() << "\t " << it << "\t ";

      /* status ! */
      if (!cv)
        newton_evol << "*KO*\t ";
      else
        newton_evol << "OK\t ";

      for (auto &itr : incr_var_convergence_)
        newton_evol << itr << "\t ";
    }

  //ha ha ha
  if (ok && cv)
    {
      for (int i = 0; i < pb.nombre_d_equations(); i++)
        if (pb.equation(i).positive_unkown())
          {
            std::string nom_inco = pb.equation(i).inconnue().le_nom().getString();
            unknown_positivation(inco[nom_inco]->valeurs(), *incr[nom_inco]);
          }

      pb.mettre_a_jour(t); //unknowns -> medium -> conserved fields
      for (auto &&n_i : inco)
        n_i.second->futur() = n_i.second->valeurs();
      eq_qdm.pression().futur() = eq_qdm.pression().valeurs();
      ConstDoubleTab_parts ppart(inco["pression"]->valeurs());
      // in multiphase, pressure is already in Pa
      /* if pression_pa() is smaller than pression() (e.g. auxiliary variables in PolyMAC_HFV), copy only the first part */
      eq_qdm.pression_pa().valeurs() = eq_qdm.pression_pa().valeurs().dimension_tot(0) < inco["pression"]->valeurs().dimension_tot(0) ? ppart[0] : inco["pression"]->valeurs(); // in multiphase, pressure is already in Pa
      first_call_ = 0;
    }
  else
    {
      for (auto &&n_i : inco)
        n_i.second->futur() = n_i.second->valeurs() = n_i.second->passe();
      if (err_a_sum > crit_conv_["alpha"])
        Cerr << que_suis_je() + ": pressure solver inaccuracy detected (requested precision " << crit_conv_["alpha"] << " , achieved precision " << err_a_sum << " )" << finl;
      ok = 0;
    }
  return;
}

int SETS::eliminer(const std::vector<std::set<std::pair<std::string, int>>> ordre,
                   const std::string inco_p,
                   const std::map<std::string, matrices_t>& mats,
                   const ptabs_t& sec,
                   std::map<std::string, Matrice_Morse>& A_p,
                   tabs_t& b_p)
{

  /* split the unknowns from sec into parts via DoubleTab_parts */
  std::map<std::pair<std::string, int>, int> offs; //offs[{inco, bloc}]: offset of block k of unknown inco
  std::map<std::pair<std::string, int>, std::array<int, 2>> dims; //dims[{inco, bloc}]: dimension 0 / line_size() of block k of inco
  for (auto &&n_v : sec)
    {
      ConstDoubleTab_parts part(*n_v.second);
      for (int i = 0; i < part.size(); i++)
        {
          offs[ { n_v.first, i }] = offs.count( { n_v.first, i - 1 }) ? offs[ { n_v.first, i - 1 }] + dims[ { n_v.first, i - 1 }][0] * dims[ { n_v.first, i - 1 }][1] : 0;
          dims[ { n_v.first, i }] = { part[i].dimension_tot(0), part[i].line_size() };
        }
    }

  /* loop over blocks */
  std::set<std::pair<std::string, int>> e_ib; //list of { variable, block } already eliminated
  std::set<std::string> e_i; //list of variables already eliminated
  int infoo {0}; // in fortran call
  int prems = !A_p.size(); //if A_p is empty, first pass -> sizing is needed
  int oMg {0}, M {0};
  const Matrice_Morse *A;
  char trans = 'T';

  for (auto &&bloc : ordre)
    {
      std::set<std::string> i_bloc, dep; //block variables, already-eliminated variables on which the block depends
      for (auto &&i_b : bloc)
        i_bloc.insert(i_b.first);
      for (auto &&i_b : bloc)
        for (auto &&v_m : mats.at(i_b.first))
          if (v_m.second->nb_colonnes() && v_m.first != inco_p && e_i.count(v_m.first) && !i_bloc.count(v_m.first))
            dep.insert(v_m.first);

      /* rows of the block to process */
      std::pair<std::string, int> i_b0 = *bloc.begin(); //first unknown of the block
      IntTrav calc(dims[i_b0][0]); //calc[i] = 1 if item i must be processed
      A = mats.at(i_b0.first).at(i_b0.first);
      oMg = offs[i_b0];
      M = dims[i_b0][1];
      for (int i = 0; i < calc.size_array(); i++)
        if (A->get_tab1()(oMg + M * i + 1) > A->get_tab1()(oMg + M * i))
          calc(i) = 1; //process all rows whose matrix is filled

      if (prems) //first pass -> sizing of A_p
        {
          /* verification de la compatibilite des inconnues du bloc -> avec les MD_Vector renseignes dans sec */
          for (auto &&i_b : bloc)
            if (dims[i_b0][0] != dims[i_b][0])
              {
                Cerr << "SETS::eliminer() : discretization of unknowns" << i_b0.first << "/" << i_b0.second << " and " << i_b.first << "/" << i_b.second << " incompatible!" << finl;
                Process::exit();
              }

          std::vector<std::set<int>> stencil(calc.size_array()); //stencil[i] -> stencil of item i (to be multiplied by line_size() of each variable)
          for (auto &&i_b : bloc)
            for (auto &&v_m : mats.at(i_b.first))
              if (v_m.second->nb_coeff())
                {
                  oMg = offs[i_b];
                  M = dims[i_b][1];
                  if (v_m.first == inco_p) //direct dependence on inco_p
                    {
                      for (int i = 0; i < calc.size_array(); i++)
                        if (calc[i])
                          for (auto j = v_m.second->get_tab1()(oMg + M * i) - 1; j < v_m.second->get_tab1()(oMg + M * (i + 1)) - 1; j++) //dependences of all rows
                            stencil[i].insert(v_m.second->get_tab2()(j) - 1);
                    }
                  else if (e_i.count(v_m.first) || i_bloc.count(v_m.first)) //dependence on a partially / fully eliminated variable
                    {
                      A = A_p.count(v_m.first) ? &A_p.at(v_m.first) : nullptr;
                      for (int i = 0; i < calc.size_array(); i++)
                        if (calc[i])
                          for (auto j = v_m.second->get_tab1()(oMg + M * i) - 1; j < v_m.second->get_tab1()(oMg + M * (i + 1)) - 1; j++)
                            {
                              const int jb = v_m.second->get_tab2()(j) - 1;
                              int k = 0;
                              while (offs.count( { v_m.first, k + 1 }) && jb >= offs.at( { v_m.first, k + 1 }))
                                k++; //l: block of the unknown on which we depend
                              const int oNg = offs[ { v_m.first, k }];
                              const int N = dims[ { v_m.first, k }][1];
                              const int n = jb - oNg - N * i;
                              if (bloc.count( { v_m.first, k }) && n >= 0 && n < N)
                                continue; //(variable, block) currently being eliminated and block-diagonal coefficient -> ok
                              else if (e_ib.count( { v_m.first, k }))  //(variable, block) eliminated -> dependence on inco_p in A_p
                                for (auto l = A->get_tab1()(jb) - 1; l < A->get_tab1()(jb + 1) - 1; l++)
                                  stencil[i].insert(A->get_tab2()(l) - 1);
                              else
                                {
                                  Cerr << "SETS::eliminer() : dependance ( " << i_b.first << "/" << i_b.second << " , " << v_m.first << "/" << k << " ) interdite!" << finl;
                                  Process::exit();
                                }
                            }
                    }
                  else
                    {
                      Cerr << "SETS::eliminer() : dependance ( " << i_b.first << "/" << i_b.second << " , " << v_m.first << " ) interdite!" << finl;
                      Process::exit();
                    }
                }

          for (auto &&i_bl : bloc) //stencil per unknown -> expanding by line_size()
            {
              Stencil sten(0, 2);

              oMg = offs[i_bl];
              M = dims[i_bl][1];
              for (int i = 0; i < calc.size_array(); i++)
                if (calc[i])
                  for (int m = 0; m < M; m++)
                    for (auto &&col : stencil[i])
                      sten.append_line(oMg + M * i + m, col);
              Matrice_Morse mat2;
              Matrix_tools::allocate_morse_matrix(sec.at(i_bl.first)->size_totale(), sec.at(inco_p)->size_totale(), sten, mat2);
              A_p[i_bl.first].nb_colonnes() ? A_p[i_bl.first] += mat2 : A_p[i_bl.first] = mat2; //A_p may already be partially created
            }
        }

      std::vector<std::string> vbloc(i_bloc.begin(), i_bloc.end());
      std::vector<std::string> vdep(dep.begin(), dep.end()); //as a list
      const int nv = (int) vbloc.size(), nd = (int) vdep.size(); //number of block variables, dependencies, total size
      std::vector<int> size, off_l, off_g; //per (variable, block): size in local system, offset in local system, offset in global systems
      int nb {0};
      for (auto &&i_b : bloc)
        {
          off_l.push_back(nb);
          size.push_back(dims[i_b][1]);
          off_g.push_back(offs[i_b]);
          nb += size.back();
        }

      std::vector<const Matrice_Morse*> pmat(nv), dAp(nd); //per variable: direct dependence on inco_p, for block variables / already resolved
      std::vector<std::vector<const Matrice_Morse*>> mat(nv), dmat(nv); //matrices of block variables, of dependencies
      std::vector<const DoubleTab*> vsec(nv), dbp(nd); //right-hand sides of variables, b_p vectors of variables / dependencies

      std::vector<Matrice_Morse*> Ap(nv); //results: d{inco} = A_p[inco].d{inco_p} + b_p[inco] for the block variables
      std::vector<DoubleTab*> bp(nv);

      for (int i = 0; i < nv; i++)
        {
          Ap[i] = &A_p[vbloc[i]];
          vsec[i] = sec.at(vbloc[i]);
          if (!b_p.count(vbloc[i]))
            b_p[vbloc[i]] = *vsec[i]; //create the b_p
          bp[i] = &b_p[vbloc[i]];
          const matrices_t& line = mats.at(vbloc[i]);
          pmat[i] = line.count(inco_p) && line.at(inco_p)->nb_colonnes() ? line.at(inco_p) : nullptr;
          mat[i].resize(nv);
          for (int j = 0; j < nv; j++)
            mat[i][j] = line.count(vbloc[j]) && line.at(vbloc[j])->nb_colonnes() ? line.at(vbloc[j]) : nullptr;
          dmat[i].resize(nd);
          for (int j = 0; j < nd; j++)
            dmat[i][j] = line.count(vdep[j]) && line.at(vdep[j])->nb_colonnes() ? line.at(vdep[j]) : nullptr;
        }

      //b_p / A_p of the dependencies
      for (int i = 0; i < nd; i++)
        {
          dbp[i] = &b_p.at(vdep[i]);
          dAp[i] = &A_p.at(vdep[i]);
        }

      DoubleTrav D(nb, nb); // bloc diagonal
      DoubleTrav S; // seconds membres

      IntTrav piv(nb);
      for (int i = 0; i < calc.size_array(); i++)
        if (calc[i])
          {
            const auto deb = Ap[0]->get_tab1()(off_g[0] + size[0] * i) - 1;
            const auto fin = Ap[0]->get_tab1()(off_g[0] + size[0] * i + 1) - 1;
            const int ic = (int)(fin - deb);
            const int nc = ic + 1;
            S.resize(nc, nb), S = 0; //right-hand side: S(i, .) -> dependence on the i-th stencil column of the block Ap, S(ic, .) -> constant part
            //right-hand side of the equations
            for (int j = 0; j < nv; j++)
              {
                M = size[j];
                oMg = off_g[j];
                const int oMl = off_l[j];
                for (int m = 0; m < M; m++)
                  S(ic, oMl + m) = vsec[j]->addr()[oMg + M * i + m];
              }

            /* remplissage par les matrices du bloc : diagonale, second membre (si partie d'une variable deja eliminee) */
            D = 0;
            for (int j = 0; j < nv; j++)
              {
                M = size[j];
                oMg = off_g[j];
                const int oMl = off_l[j];
                for (int k = 0; k < nv; k++)
                  if (mat[j][k])
                    {
                      const int N = size[k];
                      const int oNg = off_g[k];
                      const int oNl = off_l[k];
                      for (int m = 0; m < M; m++)
                        for (auto l = mat[j][k]->get_tab1()(oMg + M * i + m) - 1; l < mat[j][k]->get_tab1()(oMg + M * i + m + 1) - 1; l++)
                          {
                            const int jb = mat[j][k]->get_tab2()(l) - 1;
                            const int n = jb - oNg - N * i;
                            if (n >= 0 && n < N) //we are in the diagonal block
                              D(oMl + m, oNl + n) = mat[j][k]->get_coeff()(l);
                            else //dependence on an already eliminated block
                              {
                                const double coeff = mat[j][k]->get_coeff()(l);
                                S(ic, oMl + m) -= coeff * bp[k]->addr()[jb];
                                auto pos = deb - 1;
                                for (auto lb = Ap[k]->get_tab1()(jb) - 1; lb < Ap[k]->get_tab1()(jb + 1) - 1; lb++)
                                  {
                                    const int col = Ap[k]->get_tab2()(lb);
                                    pos++;
                                    while(Ap[0]->get_tab2()(pos) != col && pos < fin)
                                      pos++;
                                    assert(Ap[0]->get_tab2()(pos) == col);
                                    S((int)(pos - deb), oMl + m) -= coeff * Ap[k]->get_coeff()(lb);
                                  }
                              }
                          }
                    }
              }
            //direct dependence on inco_p -> in S([0, ic[, .)
            for (int j = 0; j < nv; j++)
              if (pmat[j])
                {
                  M = size[j];
                  oMg = off_g[j];
                  const int oMl = off_l[j];
                  for (int m = 0; m < M; m++)
                    {
                      auto pos = deb - 1;
                      for (auto k = pmat[j]->get_tab1()(oMg + M * i + m) - 1; k < pmat[j]->get_tab1()(oMg + M * i + m + 1) - 1; k++)
                        {
                          int col = pmat[j]->get_tab2()(k);
                          pos++;
                          while (Ap[0]->get_tab2()(pos) != col && pos < fin)
                            pos++;
                          assert(Ap[0]->get_tab2()(pos) == col);
                          S((int)(pos - deb), oMl + m) -= pmat[j]->get_coeff()(k);
                        }
                    }
                }
            //dependence on an out-of-block eliminated variable -> b_p contributes to S(0, .), A_p contributes to S(1..nc, .)
            for (int j = 0; j < nv; j++)
              for (int k = 0; k < nd; k++)
                if (dmat[j][k])
                  {
                    M = size[j];
                    oMg = off_g[j];
                    const int oMl = off_l[j];
                    for (int m = 0; m < M; m++)
                      for (auto l = dmat[j][k]->get_tab1()(oMg + M * i + m) - 1; l < dmat[j][k]->get_tab1()(oMg + M * i + m + 1) - 1; l++)
                        {
                          const double coeff = dmat[j][k]->get_coeff()(l);
                          const int jb = dmat[j][k]->get_tab2()(l) - 1;
                          S(ic, oMl + m) -= coeff * dbp[k]->addr()[jb]; //"constant" part
                          auto pos = deb - 1;
                          for (auto lb = dAp[k]->get_tab1()(jb) - 1; lb < dAp[k]->get_tab1()(jb + 1) - 1; lb++) //"inco_p dependence" part
                            {
                              const int col = dAp[k]->get_tab2()(lb);
                              pos++;
                              while (Ap[0]->get_tab2()(pos) != col && pos < fin)
                                pos++;
                              assert(Ap[0]->get_tab2()(pos) == col);
                              S((int)(pos - deb), oMl + m) -= coeff * dAp[k]->get_coeff()(lb);
                            }
                        }
                  }
            /* factorisation et resolution */
            // DoubleTrav D_back = D;
            F77NAME(dgetrf)(&nb, &nb, &D(0, 0), &nb, &piv(0), &infoo);
            if (infoo > 0)
              return 0; //singularity encountered -> exit before dividing by 0
            F77NAME(dgetrs)(&trans, &nb, &nc, &D(0, 0), &nb, &piv(0), &S(0, 0), &nb, &infoo);

            /* stockage : S(0, .) dans b_p, S(1..nc, .) dans A_p */
            for (int j = 0; j < nv; j++)
              {
                M = size[j];
                oMg = off_g[j];
                const int oMl = off_l[j];
                for (int m = 0; m < M; m++)
                  {
                    bp[j]->addr()[oMg + M * i + m] = S(ic, oMl + m);
                    {
                      auto l = Ap[j]->get_tab1()(oMg + M * i + m) - 1;
                      for (int k = 0; k < ic; k++, l++)
                        Ap[j]->get_set_coeff()(l) = S(k, oMl + m);
                    }
                  }
              }
          }

      for (auto &&i_b : bloc)
        e_ib.insert(i_b), e_i.insert(i_b.first);
    }
  return 1;
}

void SETS::assembler(const std::string inco_p,
                     const std::map<std::string, Matrice_Morse>& A_p,
                     const tabs_t& b_p,
                     const std::map<std::string, matrices_t>& mats,
                     const ptabs_t& sec,
                     Matrice_Morse& P,
                     DoubleTab& secmem,
                     int p_degen)
{
  secmem = *sec.at(inco_p);
  const int np = secmem.dimension_tot(0);
  const int M = secmem.line_size();

  /* calc(i) = 1 if we must fill rows [N * i, (N + 1) * i[ of the matrix */
  ArrOfBit calc(np);
  if (secmem.get_md_vector())
    secmem.get_md_vector()->get_sequential_items_flags(calc);
  else
    calc = 1;

  if (!P.nb_colonnes()) //sizing on first pass
    {
      Stencil stencil(0, 2);

      for (auto &&n_m : mats.at(inco_p))
        if (n_m.second && n_m.second->nb_colonnes())
          {
            const Matrice_Morse& Mp = *n_m.second, *Ap = n_m.first != inco_p ? &A_p.at(n_m.first) : nullptr;
            for (int i = 0; i < np; i++)
              if (calc[i])
                {
                  int ib = M * i;
                  for (int m = 0; m < M; m++, ib++)
                    for (auto j = Mp.get_tab1()(ib) - 1; j < Mp.get_tab1()(ib + 1) - 1; j++)
                      {
                        const int k = Mp.get_tab2()(j) - 1;
                        for (auto l = (Ap ? Ap->get_tab1()(k) - 1 : 0); l < (Ap ? Ap->get_tab1()(k + 1) - 1 : 1); l++)
                          stencil.append_line(ib, Ap ? Ap->get_tab2()(l) - 1 : k);
                      }
                }
          }
      tableau_trier_retirer_doublons(stencil);
      Matrix_tools::allocate_morse_matrix(secmem.size_totale(), secmem.size_totale(), stencil, P);
    }

  /* remplissage de P / secmem */
  P.get_set_coeff() = 0;
  for (auto &&n_m : mats.at(inco_p))
    if (n_m.first == inco_p && n_m.second && n_m.second->nb_colonnes())
      P += *n_m.second; /* dependance directe en inco_p -> on ajoute */
    else if (n_m.second && n_m.second->nb_colonnes()) /* dependance en une autre inconnue -> produit Mp.Ap + modif second membre */
      {
        const Matrice_Morse& Mp = *n_m.second, &Ap = A_p.at(n_m.first);
        const DoubleTab& bp = b_p.at(n_m.first);
        for (int i = 0; i < np; i++)
          if (calc[i])
            {
              int ib = M * i;
              for (int m = 0; m < M; m++, ib++)
                {
                  const auto deb = P.get_tab1()(ib) - 1;
                  const auto fin = P.get_tab1()(ib + 1) - 1;
                  for (auto j = Mp.get_tab1()(ib) - 1; j < Mp.get_tab1()(ib + 1) - 1; j++)
                    {
                      const int k = Mp.get_tab2()(j) - 1;
                      secmem(i, m) -= Mp.get_coeff()(j) * bp.addr()[k];
                      auto pos = deb - 1;
                      for (auto l = Ap.get_tab1()(k) - 1; l < Ap.get_tab1()(k + 1) - 1; l++)
                        {
                          const int col = Ap.get_tab2()(l);
                          pos++;
                          while (P.get_tab2()(pos) != col && pos < fin)
                            pos++;
                          assert(P.get_tab2()(pos) == col);
                          P.get_set_coeff()(pos) += Mp.get_coeff()(j) * Ap.get_coeff()(l);
                        }
                    }
                }
            }
      }
  const double diag = P.get_coeff()(0);
  if (p_degen && !Process::me())
    for (auto i = 0; i < P.get_tab1()(1) - 1; i++)
      P.get_set_coeff()(i) += diag; //de-degenerate the matrix
}
