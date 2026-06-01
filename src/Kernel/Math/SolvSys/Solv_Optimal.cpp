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
#include <Solv_Optimal.h>
#include <Param.h>
#include <LecFicDistribueBin.h>
#include <LecFicDiffuse.h>
#include <Matrice_Bloc.h>
#include <Matrice_Morse_Sym.h>
#include <SFichier.h>
#include <petsc_for_kernel.h>
#undef setbit // Otherwise conflict with PETSc
#include <MD_Vector_tools.h>
#include <Perf_counters.h>

Implemente_instanciable(Test_solveur,"Test_solveur",Interprete);
// XD test_solveur interprete test_solveur BRACE To test several solvers
// XD attr fichier_secmem chaine fichier_secmem OPT Filename containing the second member B
// XD attr fichier_matrice chaine fichier_matrice OPT Filename containing the matrix A
// XD attr fichier_solution chaine fichier_solution OPT Filename containing the solution x
// XD attr nb_test entier nb_test OPT Number of tests to measure the time resolution (one preconditionnement)
// XD attr impr rien impr OPT To print the convergence solver
// XD attr solveur solveur_sys_base solveur OPT To specify a solver
// XD attr fichier_solveur chaine fichier_solveur OPT To specify a file containing a list of solvers
// XD attr genere_fichier_solveur floattant genere_fichier_solveur OPT To create a file of the solver with a threshold
// XD_CONT convergence
// XD attr seuil_verification floattant seuil_verification OPT Check if the solution satisfy ||Ax-B||<precision
// XD attr pas_de_solution_initiale rien pas_de_solution_initiale OPT Resolution isn\'t initialized with the solution x
// XD attr ascii rien ascii OPT Ascii files

Implemente_instanciable_sans_constructeur_ni_destructeur(Solv_Optimal,"Solv_Optimal",solv_iteratif);
// XD optimal solveur_sys_base optimal BRACE Optimal is a solver which tests several solvers of the previous list to
// XD_CONT choose the fastest one for the considered linear system.
// XD attr seuil floattant seuil REQ Convergence threshold
// XD attr impr rien impr OPT To print the convergency of the fastest solver
// XD attr quiet rien quiet OPT To disable printing of information
// XD attr save_matrice|save_matrix rien save_matrice OPT To save the linear system (A, x, B) into a file
// XD attr frequence_recalc entier frequence_recalc OPT To set a time step period (by default, 100) for re-checking the
// XD_CONT fatest solver
// XD attr nom_fichier_solveur chaine nom_fichier_solveur OPT To specify the file containing the list of the tested
// XD_CONT solvers
// XD attr fichier_solveur_non_recree rien fichier_solveur_non_recree OPT To avoid the creation of the file containing
// XD_CONT the list

Sortie& Test_solveur::printOn(Sortie& s ) const
{
  return s;
}

Entree& Test_solveur::readOn(Entree& is )
{
  return is;
}

void test_un_solveur(SolveurSys& solveur, const Matrice_Base& matrice, const DoubleVect& secmem, DoubleVect& solution, int nmax, ArrOfDouble& temps, double seuil_verification=DMAXFLOAT)
{
  DoubleVect solution_ref(solution);
  int n=temps.size_array();
  statistics().create_custom_counter("Custom solver",1);
  for (int i=0; i<n; i++)
    {
      solution=solution_ref;

      // etape de resolution
      statistics().begin_count("Custom solver",statistics().get_last_opened_counter_level()+1);
      double t_0,t;
      t_0 = statistics().get_total_time("Custom solver");
      Cout<<"------------------------------------"<<finl;
      Cout<<"Try " << i << " of solver " << solveur <<finl;
      //solveur->fixer_limpr(0);
      solveur.nommer("test_solver");

      solveur.resoudre_systeme(matrice,secmem,solution);
      t = statistics().get_time_since_last_open("Custom solver");
      statistics().end_count("Custom solver");
      // on recupere un delta time et non un time absolu !!
      DoubleVect test(secmem);
      test*=-1;
      matrice.ajouter_multvect(solution,test);
      //test-=secmem;
      double norme=mp_norme_vect(test);
      double time_resol = t-t_0;
      Cout<<"CPU= " <<time_resol<<" s , ||Ax -b||= " << norme << finl;
      Process::imprimer_ram_totale();
      temps[i]=Process::mp_max(time_resol);

      if (norme>seuil_verification)
        {
          Cerr<<"residue calculated greater than the threshold value indicated ("<<seuil_verification<<")"<<finl;
          Cerr<<" one will not use the solver "<<solveur<<finl;
          temps=1.e37;
        }
    }
}
int test_solveur(SolveurSys& solveur,  const Matrice_Base& matrice , const DoubleVect& secmem , DoubleVect& solution , int nmax, ArrOfDouble& temps, Entree& list_solveur, double seuil_verification=DMAXFLOAT)
{
  DoubleVect solution_ref(solution);
  int numero=0,numero_best=-1;
  double best_time=1e36;
  Motcle motsolveur("solveur"),motlu;
  list_solveur >> motlu;
  int dernier=temps.size_array()-1;
  while (list_solveur.good())
    {
      if (motlu!=motsolveur)
        {
          Cerr<<" One expected "<<motsolveur <<" and not "<<motlu<<finl;
          Process::exit();
        }
      numero++;
      SolveurSys newsolveur;
      list_solveur>>newsolveur;
      solution=solution_ref;
      test_un_solveur(newsolveur, matrice, secmem ,solution, nmax, temps, seuil_verification);
      if (temps[dernier]<best_time)
        {
          solveur=newsolveur;
          best_time=temps[dernier];
          numero_best=numero;
        }
      list_solveur>>motlu;
    }
  if (numero_best==-1)
    {
      Cerr<<" None of the tested solvers give the expected residue" <<finl;
      Cerr<< " does the same solver must be kept ? "<<finl;
      // exit();
    }
  temps=best_time;
  return numero_best;
}
/*! @brief Generates a solver test file that differs based on whether the matrix can be solved with or without GCP.
 *
 */
void generate_defaut(const Matrice_Base& matrice, const double seuil, Sortie& sortie, int limpr=0)
{
  Nom impr(" impr " );
  if (limpr==0) impr=" ";
  if (limpr==-1) impr=" quiet ";
  if (Process::je_suis_maitre())
    {
      if((!sub_type(Matrice_Morse_Sym,matrice))&&(!sub_type(Matrice_Bloc,matrice)))
        {
          sortie <<" solveur gmres { diag seuil "<<seuil <<" "<<impr<<"}"<<finl;
#ifdef __PETSCKSP_H
          sortie <<" solveur petsc bicgstab { precond diag { }                  seuil "<<seuil <<" "<<impr<<"}"<<finl;
#endif
        }
      else
        {
          // pressure-type matrix (?)
          // Update of tested solvers 24/05/2012
          sortie <<" solveur gcp       { precond ssor       { omega 1.6 }       seuil "<<seuil <<" "<<impr<<"}"<<finl;
#ifdef __PETSCKSP_H
          sortie <<" solveur petsc gcp { precond ssor       { omega 1.6 }       seuil "<<seuil <<" "<<impr<<"}"<<finl;
          if (Process::nproc()<512)
            sortie <<" solveur petsc cholesky { impr }"<< finl;
          else
            {
              // For very large runs, switch from Cholesky to BICGSTAB ILU_SP(1) block preconditioner
              sortie <<" solveur petsc bicgstab { precond block_jacobi_icc { level 1 } seuil "<<seuil <<" "<<impr<<"}"<<finl;
              // Consider CG ILU_SP(1) by block because BICGSTAB may struggle to converge during the initial projection...
              sortie <<" solveur petsc gcp { precond block_jacobi_icc { level 1 } seuil "<<seuil <<" "<<impr<<"}"<<finl;
            }
          // SPAI has never proven effective (like all Hypre preconditioners), so it is removed
          //sortie <<" solveur petsc gcp { precond spai       { level 2 epsilon 0.2 } seuil "<<seuil <<" "<<impr<<"}"<<finl;
#endif
        }
    }
}
Entree& Test_solveur::interpreter(Entree& is)
{
  Matrice matrice;
  DoubleVect secmem,solution;
  Nom fichier_secmem("Secmem.sa");
  Nom fichier_solution("Solution.sa");
  Nom fichier_matrice("Matrice.sa");
  Nom fichier_solveur;
  bool pas_de_solution_init=false;
  bool ascii = false;
  bool limpr_ = false;
  double seuil_verification=DMAXFLOAT;
  SolveurSys solveur;
  double seuil_list=0;
  int nb_test=2; // Each solver is tested twice because the first run may be penalised by preconditioning cost
  Param  param((*this).que_suis_je());
  param.ajouter("fichier_secmem",&fichier_secmem);  // filename containing the right-hand side (Secmem.sa by default)
  param.ajouter("fichier_matrice",&fichier_matrice);  // filename containing the matrix (Matrice.sa by default)
  param.ajouter("fichier_solution",&fichier_solution);  // filename containing the solution (Solution.sa by default)
  param.ajouter("nb_test",&nb_test);  // number of solves to measure time (single preconditioning)
  param.ajouter_flag("impr",&limpr_); // enable solver output
  param.ajouter("solveur",&solveur); // specify a solver
  param.ajouter("fichier_solveur",&fichier_solveur); // specify a file containing solvers
  param.ajouter("genere_fichier_solveur",&seuil_list); // generate the solver file with a given threshold
  param.ajouter("seuil_verification",&seuil_verification); // check if the solution satisfies ||Ax-b|| < seuil_verification
  param.ajouter_flag("pas_de_solution_initiale",&pas_de_solution_init); // pas_de_solution_initiale: do not initialize the solve with the current solution
  param.ajouter_flag("ascii",&ascii); // when files are in ASCII format
  param.lire_avec_accolades_depuis(is);
  int binaire=1;
  if (ascii)
    binaire=0;
  // Re-read the matrix and the right-hand side
  {
    LecFicDistribue entree;
    entree.set_bin(binaire);
    entree.ouvrir(fichier_matrice);
    entree>>matrice;
    Cout<<" size of system "<<matrice.valeur( ).nb_colonnes()<<finl;
    //matrice->imprimer_formatte(Cout);
  }
  {
    LecFicDistribue entree;
    entree.set_bin(binaire);
    entree.ouvrir(fichier_secmem);
    MD_Vector_tools::restore_vector_with_md(secmem,entree);
    Cout<<" size of system "<<secmem.size_totale()<<finl;
  }


  if (pas_de_solution_init)
    {
      solution=secmem;
      solution=0.;
    }
  else
    {
      LecFicDistribue entree;
      entree.set_bin(binaire);
      entree.ouvrir(fichier_solution);

      MD_Vector_tools::restore_vector_with_md(solution,entree);
      Cout<<" size of system "<<solution.size_totale()<<finl;
      solution.set_md_vector(secmem.get_md_vector());
    }

  if (seuil_list!=0)
    {
      if (fichier_solveur==Nom())
        fichier_solveur="list_solveur";
      if (je_suis_maitre())
        {
          SFichier list_solveur(fichier_solveur);
          generate_defaut(matrice,seuil_list,list_solveur,limpr_);
        }
    }

  secmem.echange_espace_virtuel();
  solution.echange_espace_virtuel();
  ArrOfDouble temps(nb_test);
  if (fichier_solveur==Nom())
    test_un_solveur(solveur,  matrice , secmem , solution , -10, temps,seuil_verification);
  else
    {
      LecFicDiffuse list_solveur(fichier_solveur);
      int numero_best=test_solveur(solveur,  matrice , secmem , solution  , -10, temps,list_solveur,seuil_verification);
      Cout <<"------------------------------------------------"<<finl;
      Cout <<"Best solver : number "<<numero_best<<" "<<solveur<<finl;
      Cout <<"Best CPU time = "<<temps[0]<<finl;
      Cout <<"------------------------------------------------"<<finl;
    }
  return is;
}

static int numero_solv_optimal=0; // solver number, used to give distinct default solver file names to each solver
Solv_Optimal::Solv_Optimal():n_resol_(0),n_reinit_(0)
{
  freq_recalc_ = (int)(pow(2.0,(double)((sizeof(int)*8)-1))-1);
  freq_recalc_ = 100;
  fichier_solveur_="solveurs_";
  fichier_solveur_+=Nom(numero_solv_optimal);
  numero_solv_optimal++;

}
Solv_Optimal::~Solv_Optimal()
{
  if (le_solveur_) Cerr<<" The solver used by Solv_Optimal was "<<le_solveur_<<finl;
}
Sortie& Solv_Optimal::printOn(Sortie& s ) const
{
  return s;
}

Entree& Solv_Optimal::readOn(Entree& is )
{
  bool impr = false;
  bool quiet = false;
  Param param((*this).que_suis_je());
  param.ajouter("seuil",&seuil_,Param::REQUIRED); // convergence threshold
  param.ajouter_flag("impr",&impr); // enable solver output
  param.ajouter_flag("quiet",&quiet);
  param.ajouter("save_matrice|save_matrix",&save_matrice_); // save the linear system A, x, b
  param.ajouter("frequence_recalc",&freq_recalc_); // frequency for re-evaluating the optimal solver
  param.ajouter("nom_fichier_solveur",&fichier_solveur_); // filename containing the tested solvers
  param.ajouter_flag("fichier_solveur_non_recree",&fichier_solveur_non_recree_); // if set, the file is not created at the start of the computation
  param.lire_avec_accolades_depuis(is);
  fixer_limpr(impr);
  if (quiet)
    fixer_limpr(-1);

  return is;
}

/*! @brief Key method of Solv_Optimal: at the first iteration,
 *
 *     generates the file fichier_solveur_ containing the list of solvers to test,
 *     picks the first solver from the file.
 *     At iteration 3 and every freq_recalc_ iterations thereafter, finds the fastest solver,
 *     taking into account whether the matrix has changed.
 *     Calls test_solveur.
 *
 */
void Solv_Optimal::prepare_resol(const Matrice_Base& matrice, const DoubleVect& secmem, DoubleVect& solution, int nmax)
{
  if (n_resol_==0)
    {
      // First generate the default file, then create the solver from the first entry in the file

      if ((!fichier_solveur_non_recree_)&&(je_suis_maitre()))
        {
          SFichier list_solveur(fichier_solveur_);
          generate_defaut(matrice,seuil_,list_solveur,limpr());
        }
      LecFicDiffuse list_solveur2(fichier_solveur_);
      Nom mot;
      list_solveur2>>mot;
      list_solveur2>>le_solveur_;
    }
  // OK
  n_resol_++;
  if ((n_resol_>=3)&&((n_resol_-3)%freq_recalc_==0))
    {
      LecFicDiffuse list_solveur(fichier_solveur_);
      int nb_ite;
      if (n_reinit_<2)
        nb_ite=2; // Constant matrix (solved twice; only the second time is kept, free of any preconditioning overhead)
      else
        nb_ite=1; // Non-constant matrix (solved only once)
      ArrOfDouble temps(nb_ite);
      statistics().end_count(STD_COUNTERS::system_solver,0,0);
      int numero_best=test_solveur(le_solveur_,  matrice , secmem , solution  , nmax, temps,list_solveur,seuil_);
      statistics().begin_count(STD_COUNTERS::system_solver,statistics().get_last_opened_counter_level()+1);
      Cout <<"------------------------------------------------"<<finl;
      Cout <<"Best solver : number "<<numero_best<<" "<<le_solveur_<<finl;
      Cout <<"Best CPU time = "<<temps[0]<<finl;
      Cout <<"------------------------------------------------"<<finl;
      if (je_suis_maitre())
        {
          Nom best("best_");
          best+=fichier_solveur_;
          SFichier sf(best);
          sf << le_solveur_ <<finl;
        }
    }
}
int Solv_Optimal::resoudre_systeme(const Matrice_Base& matrice, const DoubleVect& secmem, DoubleVect& solution)
{
  prepare_resol(matrice,secmem,solution,-20);
  return le_solveur_->resoudre_systeme( matrice, secmem,  solution);
}
int Solv_Optimal::resoudre_systeme(const Matrice_Base& matrice, const DoubleVect& secmem, DoubleVect& solution, int nmax)
{
  prepare_resol(matrice,secmem,solution, nmax);
  return le_solveur_->resoudre_systeme( matrice, secmem,  solution, nmax );
}


void Solv_Optimal::reinit()
{
  n_reinit_++;
  if (le_solveur_)
    le_solveur_->reinit();
}




