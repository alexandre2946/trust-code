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

#include <Simple.h>
#include <Navier_Stokes_std.h>
#include <EChaine.h>
#include <Matrice_Bloc.h>
#include <Matrice_Morse_Sym.h>
#include <Assembleur_base.h>
#include <Schema_Temps_base.h>
#include <Schema_Euler_Implicite.h>
#include <Fluide_Dilatable_base.h>
#include <Probleme_base.h>
#include <MD_Vector_composite.h>
#include <MD_Vector_tools.h>
#include <TRUSTTab_parts.h>
#include <SETS.h>
#include <TRUSTTrav.h>
#include <Discretisation_base.h>

Implemente_instanciable_sans_constructeur(Simple,"Simple",Simpler_Base);
// XD simple piso simple INHERITS_BRACE SIMPLE type algorithm
// XD attr relax_pression floattant relax_pression OPT Value between 0 and 1 (by default 1), this keyword is used only
// XD_CONT by the SIMPLE algorithm for relaxing the increment of pressure.

Simple::Simple()
{
  alpha_ = 1.;
  beta_ = 1.;
  with_d_rho_dt_ = 1;
  Ustar_old.resize(0);
}

Sortie& Simple::printOn(Sortie& os ) const
{
  return Simpler_Base::printOn(os);
}

Entree& Simple::readOn(Entree& is )
{
  return Simpler_Base::readOn(is);
}

Entree& Simple::lire(const Motcle& motlu,Entree& is)
{
  Motcles les_mots(1);
  {
    les_mots[0] = "relax_pression";
  }

  int rang=les_mots.search(motlu);
  switch(rang)
    {
    case 0:
      {
        is >> beta_;
        break;
      }

    default :
      {
        return Simpler_Base::lire(motlu, is);
      }
    }

  return is;
}

void diviser_par_rho_np1_face(Equation_base& eqn,DoubleTab& tab_array)
{
  Fluide_Dilatable_base& fluide_dil = ref_cast(Fluide_Dilatable_base,eqn.milieu());
  int nbdim = tab_array.nb_dim();
  int taille_0_tot = tab_array.dimension_tot(0);
  int dim = tab_array.dimension(1);
  CDoubleArrView rho = static_cast<const ArrOfDouble&>(fluide_dil.rho_face_np1()).view_ro();
  if (nbdim==1)
    {
      DoubleArrView array = static_cast<ArrOfDouble&>(tab_array).view_rw();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, taille_0_tot), KOKKOS_LAMBDA(const int i)
      {
        for (int j = 0; j < dim; j++)
          array(i) /= rho(i);
      });
    }
  else
    {
      DoubleTabView array = tab_array.view_rw();
      Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), range_1D(0, taille_0_tot), KOKKOS_LAMBDA(const int i)
      {
        for (int j=0; j<dim; j++)
          array(i,j) /= rho(i);
      });
    }
  end_gpu_timer(__KERNEL_NAME__);
}

void iterer_eqn_expl(Equation_base& eqn,int nb_iter,double dt,DoubleTab& current,DoubleTab& dudt,
                     int& converge,const int no_qdm_)
{
  if (nb_iter>1)
    {
      converge = 1;
      return;
    }
  if ((no_qdm_))
    {
      Cout<<eqn.que_suis_je()<<" equation is not solved."<<finl;
      eqn.valider_iteration();
      converge = 1;
      return;

    }
  double dt_stab = eqn.calculer_pas_de_temps();

  // See Schema_Temps_base::limpr for information on modf
  double n_sous_ite;
  modf((dt/dt_stab), &n_sous_ite);
  if (dt>dt_stab*n_sous_ite) n_sous_ite=n_sous_ite+1.;
  double dt_calc = dt/n_sous_ite;
  assert(dt_calc<=dt_stab);
  for (double i=0; i<n_sous_ite; i=i+1.)
    {
      eqn.derivee_en_temps_inco(dudt);
      dudt *= dt_calc;
      current += dudt;
      eqn.valider_iteration();
    }
  Cout <<"It is solved with " << n_sous_ite << " explicit sub-iterations."<<finl;
  converge=1;
}

void iterer_eqn_expl_diffusion_implicite(Equation_base& eqn,int nb_iter,double dt,DoubleTab& current,DoubleTab& dudt,
                                         int& converge,const int no_qdm_,const double seuil)
{
  if (nb_iter>1)
    {
      converge = 1;
      return;
    }
  if ((no_qdm_))
    {
      Cout<<eqn.que_suis_je()<<" equation is not solved."<<finl;
      eqn.valider_iteration();
      converge = 1;
      return;
    }

  //dudt = current;
  Schema_Temps_base& sch = eqn.probleme().schema_temps();
  double seuil_sa = sch.seuil_diffusion_implicite();
  int flag_sa = sch.diffusion_implicite();
  sch.set_seuil_diffusion_implicite() = seuil;
  sch.set_diffusion_implicite() = 1;
  eqn.derivee_en_temps_inco(dudt);
  dudt *= dt;
  current += dudt;
  eqn.valider_iteration();

  //Restore the implicit diffusion parameters to their previous values
  sch.set_seuil_diffusion_implicite() = seuil_sa;
  sch.set_diffusion_implicite() = flag_sa;
  converge = 1;
}

bool Simple::iterer_eqn(Equation_base& eqn,const DoubleTab& inut,DoubleTab& current,
                        double dt,int nb_iter, int& ok)
{
  DoubleTab dudt(current);
  Parametre_implicite& param = get_and_set_parametre_implicite(eqn);


  int converge=0;
  int nb_pas_dt = eqn.schema_temps().nb_pas_dt();

  double temps = eqn.schema_temps().temps_courant();
  int freq = param.equation_frequence_resolue(temps);

  if (freq!=0 && (nb_pas_dt%freq)!=0)
    {
      Cout<<"Time step : "<<nb_pas_dt<<" - "<<eqn.que_suis_je()<<" equation is not solved."<<finl;
      eqn.valider_iteration();
      return true;
    }
  if (eqn.equation_non_resolue())
    {
      Cout<<eqn.que_suis_je()<<" equation is not solved."<<finl;
      // compute the derivative once to obtain the boundary fluxes
      if  (eqn.schema_temps().nb_pas_dt()==0)
        {
          DoubleTab toto(eqn.inconnue().valeurs());
          eqn.derivee_en_temps_inco(toto);
        }
      eqn.valider_iteration();
      converge = 1;
      return true;
    }
  if (param.calcul_explicite())
    {
      if (param.seuil_diffusion_implicite()<0)
        iterer_eqn_expl(eqn,nb_iter,dt,current,dudt,converge,no_qdm_);
      else
        iterer_eqn_expl_diffusion_implicite(eqn,nb_iter,dt,current,dudt,converge,no_qdm_,param.seuil_diffusion_implicite());

      return (converge==1);
    }

  double& seuil_convergence_implicite = param.seuil_convergence_implicite();
  double seuil_verification_solveur = param.seuil_verification_solveur();
  double seuil_test_preliminaire_solveur = param.seuil_test_preliminaire_solveur();
  assert(seuil_verification_solveur>0);
  SolveurSys& solveur = param.solveur();
  double seuil_convg = seuil_convergence_implicite;
  if (is_seuil_convg_variable)
    {
      double residu = ref_cast(Schema_Euler_Implicite,eqn.schema_temps()).residu_old();

      if ((residu*dt*10<seuil_convg)&&(nb_iter==1))
        {
          seuil_convergence_implicite /= 1.1;
          seuil_convg = seuil_convergence_implicite;
        }
      Cout<<"seuil_convg "<<seuil_convg<<finl;
    }

  converge = 0;
  if ((no_qdm_) && (sub_type(Navier_Stokes_std,eqn)))
    {
      Matrice& matrice_en_pression_2 = ref_cast(Navier_Stokes_std, eqn).matrice_pression();
      if (!matrice_en_pression_2) matrice_en_pression_2.detach();

      converge = 1;
      return true;
    }

  //////////////////////////////////////////////////////////////////////////////////////////
  // Implicit resolution - via iterer_NS for Navier_Stokes
  //                     - ...solveur.resoudre_systeme()... for other equations
  /////////////////////////////////////////////////////////////////////////////////////////

  dudt = current; // to be able to test convergence.
  Matrice_Morse matrice;
  if (!(sub_type(Navier_Stokes_std,eqn) && sub_type(SETS, *this))) //SETS and ICE manage their own matrices
    {
      eqn.dimensionner_matrice(matrice);
      matrice.get_set_coeff() = 0;
    }

  DoubleTrav resu(current);

  if( sub_type(Navier_Stokes_std,eqn))
    {
      Navier_Stokes_std& eqnNS = ref_cast(Navier_Stokes_std,eqn);
      DoubleTab& pression = eqnNS.pression().valeurs();
      DoubleTrav secmem(pression);
      pression.echange_espace_virtuel();
      iterer_NS(eqnNS,current,pression,dt,matrice,seuil_verification_solveur,secmem,nb_iter,converge, ok);
    }
  else
    {
      solveur->reinit();
      DoubleTrav resu_temp(current); /* residu en increments */
      if (eqn.has_interface_blocs()) /* if assembler_blocs is available */
        {
          if (eqn.discretisation().is_poly_family() || eqn.que_suis_je().debute_par("Equation_flux"))
            {
              eqn.assembler_blocs_avec_inertie({{ eqn.inconnue().le_nom().getString(), &matrice }}, resu_temp, { });
              resu = resu_temp;
              matrice.ajouter_multvect(current, resu);
            }
          else
            {
              eqn.assembler_blocs_avec_inertie({{ eqn.inconnue().le_nom().getString(), &matrice }}, resu, { });
              resu_temp = 0;
              matrice.ajouter_multvect(current,resu_temp);
              resu_temp -= resu;
            }

        }
      else
        {
          eqn.assembler_avec_inertie(matrice,current,resu);
          resu_temp = 0;
          matrice.ajouter_multvect(current,resu_temp);
          resu_temp -= resu;
        }

      if (seuil_test_preliminaire_solveur>0)
        {
          double norme_b=mp_norme_vect(resu_temp);
          if (norme_b<seuil_test_preliminaire_solveur)
            {
              //  GF the following test might be sensible?
              //      if ( nb_iter>1)
              {
                converge = 1;
                Cout<<"Resolution of system is not necessary: "<< norme_b <<" < "<< seuil_test_preliminaire_solveur <<finl;
              }
            }
        }
      if (converge==0)
        {
          int con = 0;
          while (con==0)
            {
              con = 1;
              solveur.resoudre_systeme(matrice,resu,current);
              if (eqn.positive_unkown())
                for (int i = 0; i < current.dimension_tot(0); i++)
                  for (int j = 0; j < current.line_size(); j++)
                    current(i, j) = std::max(current(i, j), 0.);
              ok = eqn.milieu().check_unknown_range(); //verify that the unknown is within the medium's bounds

              if (ok)
                {
                  resu_temp = 0;
                  matrice.ajouter_multvect(current,resu_temp);
                  resu_temp -= resu;
                  double norme_resu = mp_norme_vect(resu_temp) ;

                  if (norme_resu>seuil_verification_solveur) con = 0;
                  if (con==0)
                    {
                      Cout<<"Residu of iterative solving : "<<norme_resu<<" instead of "<<seuil_verification_solveur<<finl;
                      Cout<<"Iterating on the linear system again:"<< finl;
                      seuil_verification_solveur *= 2;
                      con = 0;
                    }
                }
              else current = eqn.inconnue().passe(); //if ok == 0, restore the previous value of the unknown
            }
          converge = 0;
        }
    }

  ///////////////////////////////////////////////////////////////////////
  // Convergence test of the solution between two successive iterations
  // Not applied for the N_S unknown with the PISO and Implicit algorithm
  ///////////////////////////////////////////////////////////////////////

  if(!converge && ok)
    {
      // allows checking what happens
      // in particular the positivity of K and eps
      eqn.valider_iteration();
      dudt -= current;
      double dudt_norme = mp_norme_vect(dudt);
      converge = (dudt_norme < seuil_convg);
      if (!converge)
        {
          Cout<<eqn.que_suis_je()<<" is not converged at the implicit iteration "<<nb_iter<<" ( ||uk-uk-1|| = "<<dudt_norme<<" > implicit threshold "<<seuil_convg<<" )"<<finl;
          if (nb_iter>=10) Cout << "Consider lowering facsec_max value. Look at the reference manual for advice to set facsec_max value according to the problem type." << finl;
        }
      else
        Cout<<eqn.que_suis_je()<<" is converged at the implicit iteration "<<nb_iter<<" ( ||uk-uk-1|| = "<<dudt_norme<<" < implicit threshold "<<seuil_convg<<" )"<<finl;
    }

  if(ok && (eqn.discretisation().is_poly_family() || eqn.probleme().que_suis_je().debute_par("Pb_Multiphase"))) eqn.probleme().mettre_a_jour(eqn.schema_temps().temps_courant());
  solveur->reinit();
  return (ok && converge==1);
}

bool Simple::iterer_eqs(LIST(OBS_PTR(Equation_base)) eqs, int nb_iter, int& ok)
{
  // retrieve the linear system solver
  Parametre_implicite& param = get_and_set_parametre_implicite(eqs[0]);
  SolveurSys& solveur = param.solveur();
  double seuil_convg = param.seuil_convergence_implicite();
  int i, j, bs = 0; //bs : common line size of arrays if > 0, 0 otherwise

  /* key for memoization */
  list_of_eq_ptr_t key(eqs.size());
  for (i = 0; i < eqs.size(); i++) key[i] = (intptr_t) &eqs[i].valeur();

  int init = !mbloc.count(key); //first pass
  Matrice_Bloc& Mglob = mbloc[key];

  if (init)
    for (Mglob.dimensionner(eqs.size(), eqs.size()), i = 0; i < eqs.size(); i++)
      for (j = 0; j < eqs.size(); j++) Mglob.get_bloc(i, j).typer("Matrice_Morse");

  /* for interface_blocs: if all equations have this interface, we use it */
  int interface_blocs_ok = 1;
  for (i = 0; i < eqs.size(); i++) interface_blocs_ok &= eqs[i]->has_interface_blocs();
  std::vector<matrices_t> mats(eqs.size()); //matrix row for equation i
  for (i = 0; i < eqs.size(); i++)
    for (j = 0; j < eqs.size(); j++)
      {
        Nom nom_i = eqs[j]->inconnue().le_nom();
        // field from another problem: add a suffix
        if (eqs[i]->probleme().le_nom().getString() != eqs[j]->probleme().le_nom().getString()) nom_i += Nom("/") + eqs[j]->probleme().le_nom();
        mats[i][nom_i.getString()] = &ref_cast(Matrice_Morse, Mglob.get_bloc(i, j).valeur());
      }

  //Do the unknowns/residuals have the same shape?
  for (bs = eqs[0]->inconnue().valeurs().line_size(), i = 1; i < eqs.size(); i++)
    if (eqs[i]->inconnue().valeurs().line_size() != bs) bs = 0;

  //MD_Vector global
  MD_Vector_composite mdc; //version composite
  for (i = 0; i < eqs.size(); i++)
    mdc.add_part(eqs[i]->inconnue().valeurs().get_md_vector(), bs ? 0 : eqs[i]->inconnue().valeurs().line_size());
  MD_Vector mdv;
  mdv.copy(mdc);

  if (init) //first pass -> sizing of MD_Vector and matrices
    {
      /* dimensionnement de la matrice globale */
      if (interface_blocs_ok)
        for (i = 0; i < eqs.size(); i++)
          for (eqs[i]->dimensionner_blocs(mats[i], {}), j = 0; j < eqs.size(); j++)
            {
              Matrice_Morse& mat = ref_cast(Matrice_Morse, Mglob.get_bloc(i, j).valeur());
              if (!mat.nb_colonnes())
                mat.dimensionner(eqs[i]->inconnue().valeurs().size_totale(), eqs[j]->inconnue().valeurs().size_totale(), 0);
            }
      else for (i = 0; i < eqs.size(); i++)
          for (j = 0; j < eqs.size(); j++)
            {
              Matrice_Morse& mat = ref_cast(Matrice_Morse, Mglob.get_bloc(i, j).valeur()), mat2;
              int nl = eqs[i]->inconnue().valeurs().size_totale(), nc = eqs[j]->inconnue().valeurs().size_totale();
              if (i == j) eqs[i]->dimensionner_matrice(mat);
              eqs[i]->dimensionner_termes_croises(i == j ? mat2 : mat, eqs[j]->probleme(), nl, nc);
              if (i == j) mat += mat2;
            }
    }
  else for (i = 0; i < eqs.size(); i++)
      for (j = 0; j < eqs.size(); j++) //subsequent passes -> just reallocate the coeff() arrays
        {
          Matrice_Morse& mat = ref_cast(Matrice_Morse, Mglob.get_bloc(i, j).valeur());
          mat.get_set_coeff().resize(mat.get_set_tab2().size_array());
        }

  //work arrays
  DoubleTrav inconnues, residus, dudt;
  if (bs) inconnues.resize(0, bs), residus.resize(0, bs), dudt.resize(0, bs); //so that aggregated arrays have the correct line_size() if it exists
  MD_Vector_tools::creer_tableau_distribue(mdv, inconnues);
  MD_Vector_tools::creer_tableau_distribue(mdv, residus);
  MD_Vector_tools::creer_tableau_distribue(mdv, dudt);
  DoubleTab_parts residu_parts(residus), inconnues_parts(inconnues), dudt_parts(dudt);

  //fill unknowns
  for(i = 0; i < eqs.size(); i++) inconnues_parts[i] = eqs[i]->inconnue().valeurs();
  dudt = inconnues;

  //fill matrices
  if (interface_blocs_ok)
    {
      for (i = 0; i < eqs.size(); i++)
        {
          eqs[i]->assembler_blocs_avec_inertie(mats[i], residu_parts[i], {});
          if (!eqs[i]->discretisation().is_poly_family())
            {
              for (j = 0; j < eqs.size(); j++)
                {
                  Nom nom_j = eqs[j]->inconnue().le_nom();
                  if (eqs[i]->probleme().le_nom().getString() != eqs[j]->probleme().le_nom().getString())
                    {
                      nom_j += Nom("/") + eqs[j]->probleme().le_nom();
                      mats[i][nom_j.getString()]->ajouter_multvect(inconnues_parts[j], residu_parts[i]);
                    }
                }
            }
        }
      if (eqs[0]->discretisation().is_poly_family()) Mglob.ajouter_multvect(inconnues, residus); //to avoid solving in increments
    }
  else for(i = 0; i < eqs.size(); i++)
      for (j = 0; j < eqs.size(); j++)
        {
          Matrice_Morse& mat = ref_cast(Matrice_Morse, Mglob.get_bloc(i, j).valeur());
          eqs[i]->ajouter_termes_croises(inconnues_parts[i], eqs[j]->probleme(), inconnues_parts[j], residu_parts[i]);
          eqs[i]->contribuer_termes_croises(inconnues_parts[i], eqs[j]->probleme(), inconnues_parts[j], mat);
          /* if i == j, then assembler_avec_inertie() handles the matrix/vector product: otherwise it must be done manually */
          if (i == j) eqs[i]->assembler_avec_inertie(mat, inconnues_parts[i], residu_parts[i]);
          else mat.ajouter_multvect(inconnues_parts[j], residu_parts[i]);
        }

  // resolution
  solveur->reinit();
  solveur.resoudre_systeme(Mglob, residus, inconnues);
  inconnues.echange_espace_virtuel();

  // update
  // Optimization: combine N mp_norme_vect into 1 collective call
  // First pass: compute local squared norms
  ArrOfDouble dudt_carres((int)eqs.size());
  for(i = 0; i < eqs.size(); i++)
    {
      dudt_parts[i] -= inconnues_parts[i];
      dudt_carres[(int)i] = local_carre_norme_vect(dudt_parts[i]);
    }
  // Single MPI reduction for all norms
  Process::mp_sum_for_each_item(dudt_carres);
  // Second pass: use the norms and do updates
  bool converge = true;
  for(i = 0; i < eqs.size(); i++)
    {
      double dudt_norme = sqrt(dudt_carres[(int)i]);
      eqs[i]->inconnue().valeurs() = inconnues_parts[i];

      converge &= (dudt_norme < seuil_convg);
      if (!converge)
        {
          Cout<<eqs[i]->que_suis_je()<<" is not converged at the implicit iteration "<<nb_iter<<" ( ||uk-uk-1|| = "<<dudt_norme<<" > implicit threshold "<<seuil_convg<<" )"<<finl;
          if (nb_iter>=10) Cout << "Consider lowering facsec_max value. Look at the reference manual for advice to set facsec_max value according to the problem type." << finl;
        }
      else
        Cout<<eqs[i]->que_suis_je()<<" is converged at the implicit iteration "<<nb_iter<<" ( ||uk-uk-1|| = "<<dudt_norme<<" < implicit threshold "<<seuil_convg<<" )"<<finl;
      eqs[i]->inconnue().futur() = eqs[i]->inconnue().valeurs();
      const double t = eqs[i]->schema_temps().temps_courant() + eqs[i]->schema_temps().pas_de_temps();
      eqs[i]->domaine_Cl_dis().imposer_cond_lim(eqs[i]->inconnue(), t);
      eqs[i]->inconnue().valeurs() = eqs[i]->inconnue().futur();
      eqs[i]->inconnue().Champ_base::changer_temps(t);
    }
  for(i = 0; i < eqs.size(); i++) eqs[i]->probleme().mettre_a_jour(eqs[i]->schema_temps().temps_courant());

  //deallocate the coefficient arrays
  for (i = 0; i < eqs.size(); i++)
    for (j = 0; j < eqs.size(); j++) ref_cast(Matrice_Morse, Mglob.get_bloc(i, j).valeur()).get_set_coeff().reset();

  return converge;
}

void Simple::calculer_correction_en_vitesse(const DoubleTrav& correction_en_pression,DoubleTrav& gradP,DoubleTrav& correction_en_vitesse,const Matrice_Morse& matrice,const Operateur_Grad& gradient)
{
  int deux_entrees = 0;
  if (correction_en_vitesse.nb_dim()==2) deux_entrees = 1;
  gradient->multvect(correction_en_pression,gradP);
  int nb_comp = 1;
  if(deux_entrees)
    nb_comp = correction_en_vitesse.dimension(1);

  ConstDoubleTab_parts part(correction_en_vitesse);
  int nb_ligne_reel = part[0].dimension(0);
  if (deux_entrees==0)
    {
      // D(Uk-1)^-1 resu
      int i,j;
      for(i=0; i<nb_ligne_reel; i++)

        for (j=0; j<nb_comp; j++)
          {
            //k=tab1(i*nb_comp+j)-1;
            correction_en_vitesse(i) = -gradP(i)/matrice(i*nb_comp+j,i*nb_comp+j);
          }
    }
  else
    {
      int i,j;
      for(i=0; i<nb_ligne_reel; i++)
        for (j=0; j<nb_comp; j++)
          {
            //k = tab1(i*nb_comp+j)-1;
            correction_en_vitesse(i,j) = -gradP(i,j)/matrice(i*nb_comp+j,i*nb_comp+j);
          }
    }
  correction_en_vitesse.echange_espace_virtuel();
}


//Input: Uk-1 ; Pk-1
//Output: Uk ; Pk
//k denotes an iteration

void Simple::iterer_NS(Equation_base& eqn,DoubleTab& current,DoubleTab& pression,
                       double dt,Matrice_Morse& matrice,double seuil_resol,DoubleTrav& secmem,int nb_ite,int& converge, int& ok)
{
  Parametre_implicite& param = get_and_set_parametre_implicite(eqn);
  SolveurSys& solveur = param.solveur();

  Navier_Stokes_std& eqnNS = ref_cast(Navier_Stokes_std,eqn);
  eqnNS.reassembler_pression_si_necessaire();
  DoubleTrav gradP(current);
  DoubleTrav correction_en_pression(pression);
  DoubleTrav correction_en_vitesse(current);
  DoubleTrav resu(current);
  int is_dilat = eqn.probleme().is_dilatable();

  //int deux_entrees = 0;
  //if (current.nb_dim()==2) deux_entrees = 1;
  Operateur_Grad& gradient = eqnNS.operateur_gradient();
  Operateur_Div& divergence = eqnNS.operateur_divergence();

  /* int nb_comp = 1;
     int nb_dim = current.nb_dim();
     if (nb_dim==2)
     nb_comp = current.dimension(1);
  */

  gradient.calculer(pression,gradP);
  //Build matrix and residual
  //matrice = A[Uk-1] = M/dt + CONV + DIFF
  //resu = A[Uk-1]Uk-1 -(A[Uk-1]Uk-1-Ss) + Sv + (M/dt)Uk-1 -BtPk-1
  if (eqnNS.has_interface_blocs()) //if the interface_blocs is available, use it
    eqnNS.assembler_blocs_avec_inertie({{ "vitesse", &matrice }}, resu);
  else //otherwise, go through ajouter/contribuer
    {
      resu -= gradP;
      eqnNS.assembler_avec_inertie(matrice,current,resu);
    }

  solveur->reinit();

  //Solve the system A[Uk-1]U* = -BtP* + Sv + Ss + (M/dt)Uk-1
  //current = U*
  solveur.resoudre_systeme(matrice,resu,current);

  //Velocity field relaxation U*
  //U* = alpha U*_new + (1-alpha)*U*_old
  if (nb_ite==1)
    Ustar_old = current;
  current *= alpha_ ;
  current.ajoute(1.-alpha_,Ustar_old);
  current.echange_espace_virtuel();
  Ustar_old = current;

  //Build matrice_en_pression_2 = BD-1Bt[Uk-1]
  Matrice& matrice_en_pression_2 = eqnNS.matrice_pression();
  assembler_matrice_pression_implicite(eqnNS,matrice,matrice_en_pression_2);
  SolveurSys& solveur_pression_ = eqnNS.solveur_pression();
  solveur_pression_->reinit();

  //Compute secmem = BU* (incompressible) BU* -drho/dt (quasi-compressible)
  if (is_dilat)
    {
      if (with_d_rho_dt_)
        {
          Fluide_Dilatable_base& fluide_dil = ref_cast(Fluide_Dilatable_base,eqn.milieu());
          fluide_dil.secmembre_divU_Z(secmem);
          secmem *= -1;
        }
      else secmem = 0;
      divergence.ajouter(current,secmem);
    }
  else
    divergence.calculer(current,secmem);
  secmem *= -1;
  secmem.echange_espace_virtuel();


  //Solve the system (BD-1Bt)P' = BU* (incompressible)
  //                 (BD-1Bt)P' = BU* -drho/dt (quasi-compressible)
  //correction_en_pression = P'
  solveur_pression_.resoudre_systeme(matrice_en_pression_2.valeur(),
                                     secmem,correction_en_pression);

  //Solve DU' = BP'
  //correction_en_vitesse = U'
  calculer_correction_en_vitesse(correction_en_pression,gradP,correction_en_vitesse,matrice,gradient);

  //Pressure correction P = P* + beta_*P'
  //Velocity correction U = U* + beta_u*U' (beta_u=1)

  pression.ajoute(beta_,correction_en_pression);
  eqnNS.assembleur_pression()->modifier_solution(pression);

  current += correction_en_vitesse;

  if (is_dilat)
    diviser_par_rho_np1_face(eqn,current);
}
