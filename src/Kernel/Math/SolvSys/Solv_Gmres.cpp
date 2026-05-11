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

#include <Solv_Gmres.h>
#include <Matrice_Morse_Sym.h>
#include <Matrice_Bloc.h>
#include <Motcle.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Solv_Gmres,"Solv_Gmres",solv_iteratif);
// XD solv_gmres solveur_sys_base gmres BRACE Preconditioned GMRES.
// XD attr diag rien diag OPT Keyword to precondition with the diagonal
// XD attr seuil floattant seuil OPT Value of the final residue. The solver ceases iterations when the Euclidean residue standard ||Ax-B|| is less than this value. default value 1e-12.
// XD attr impr rien impr OPT Keyword which is used to request display of the Euclidean residue standard each time this iterates
// XD attr save_matrice|save_matrix rien save_matrice OPT To save the matrix in a file.
// XD attr quiet rien quiet OPT To not displaying any outputs of the solver.
// XD attr nb_it_max entier nb_it_max OPT Keyword to set the maximum iterations number for the solver.
// XD attr controle_residu entier controle_residu OPT Keyword of Boolean type (by default 0). If set to 1, check the convergence after solve

Solv_Gmres::Solv_Gmres()
{
  seuil_ = 1.e-12;
  nb_it_max_ = 1000000;
  controle_residu_ =0;
  dim_espace_Krilov_=10;
}

// printOn et readOn
Sortie& Solv_Gmres::printOn(Sortie& s ) const
{
  s<<" { seuil "<<seuil_;
  if (precond_diag_)
    s <<" diag ";
  else if (is_local_gmres) s<<" sans_precond ";

  if (controle_residu_) s<< " controle_residu "<<controle_residu_;
  if (nb_it_max_!=1000000) s<<" nb_it_max "<<nb_it_max_;
  if (limpr()==1) s<<" impr ";
  if (limpr()==-1) s<<" quiet ";
  if (save_matrice_) s<< " save_matrice ";
  s<<" dim_espace_krilov "<<dim_espace_Krilov_;
  s<<" } ";
  return s;
}

Entree& Solv_Gmres::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Solv_Gmres::set_param(Param& param) const
{
  param.ajouter_non_std("impr",(this));
  param.ajouter("seuil",&seuil_);
  param.ajouter_non_std("diag",(this));
  param.ajouter_non_std("sans_precond",(this));
  param.ajouter("nb_it_max",&nb_it_max_);
  param.ajouter("controle_residu",&controle_residu_);
  param.ajouter("save_matrice|save_matrix",&save_matrice_);
  param.ajouter("dim_espace_krilov",&dim_espace_Krilov_);
  param.ajouter_non_std("quiet",(this));
}

int Solv_Gmres::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  int retval = 1;

  if (mot=="impr") fixer_limpr(1);
  else if (mot=="quiet") fixer_limpr(-1);
  else if (mot=="diag")
    {
      is_local_gmres=true;
      precond_diag_=1;
    }
  else if (mot=="sans_precond")
    {
      is_local_gmres=true;
      precond_diag_=false;
    }
  else retval = -1;

  return retval;
}



int Solv_Gmres::resoudre_systeme(const Matrice_Base& la_matrice,
                                 const DoubleVect& secmem,
                                 DoubleVect& solution)
{
  if(sub_type(Matrice_Morse,la_matrice))
    {
      const Matrice_Morse& matrice = ref_cast(Matrice_Morse, la_matrice);
      return Gmres(matrice,secmem,solution);
    }
  else
    {
      if(sub_type(Matrice_Bloc,la_matrice))
        {
          const Matrice_Bloc& matrice = ref_cast(Matrice_Bloc,la_matrice);
          if(matrice.nb_bloc_lignes()>1)
            {
              Cerr<<"Solv_Gmres : WARNING : one is not able to carry out Gmres by blocks"<<finl;
              exit();
              return(-1);
            }

          if (Process::is_parallel())
            {
              Cerr<<"Solv_Gmres : WARNING : one is not able to carry out parallel calculation with Gmres"<<finl;
              exit();
              return(-1);
            }

          const Matrice_Morse& MB00 = ref_cast(Matrice_Morse,matrice.get_bloc(0,0).valeur());
          int retour= Gmres(MB00,secmem,solution);
          return retour;
        }
      else
        {
          Cerr<<"Solv_Gmres : WARNING : only linear systems based on Matrice_Morse_Sym or Matrice_Bloc type matrixes can be solved"<<finl;
          exit();
          return(-1);
        }
    }
}

int Solv_Gmres::gmres_local(const Matrice_Morse& A, const DoubleVect& b, DoubleVect& tab_x)
{
  // PL c'est pas joli et c'est de moi mais je ne comprends pas pourquoi
  // l'utilisation de b.size_reelle() plante sur la matrice en pression depuis la 1.6.0 (non teste)
  // Qu'est ce qui est fait en implicite pour b.size_reelle() soit >=0 ???
  const int ns=(b.size_reelle_ok()?b.size_reelle():b.size_array());
  int nb_ligne_tot=(int)Process::mp_sum((double) ns);

  // A present dans le jdd
  double epsGMRES=1.e-10*0;
  //int nkr_min = 10;
  //int nkr=std::max(nkr_min,nb_ligne_tot/2);                         // dimension de l'espace de Krylov
  int nkr = dim_espace_Krilov_;
  int nit1_min = 20;
  int nit1=std::max(nit1_min,nb_ligne_tot);
  int nit=std::min(nb_it_max_,nit1);
  double rec_min = seuil_;
  double rec_max = 0.1  ;
  double res2_old=-1;
  if (v.size()==0)
    {
      v.dimensionner(nkr);                         // Krilov vectors
      h.resize(nkr + 1, nkr);                // Heisenberg maatrix of coefficients
      r.resize(nkr + 1);
      h_loc.resize(nkr);
      dh_loc.resize(nkr+1);
    }

  if (tab_Diag.size_array()!=ns)
    {
      tab_v0 = tab_x;
      tab_v1 = tab_x;
      tab_Diag.resize(ns);
    }

  // Initialisation
  tab_v0 = 0.;
  tab_v1 = 0.;
  const bool precond_diag = precond_diag_;
  {
    Matrice_Morse_View matrice;
    matrice.set(const_cast<Matrice_Morse&>(A));
    DoubleArrView Diag = tab_Diag.view_wo();
    Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), ns, KOKKOS_LAMBDA(
                           const int i)
    {
      Diag[i] = precond_diag ? 1. / matrice(i, i) : 1.;
    });
    end_gpu_timer(__KERNEL_NAME__);
  }

  A.multvect_(tab_x,tab_v0);
  tab_v0 *= -1.;
  tab_v0 += b;

  // Reduce 2 mp_sum calls to 1 by computing local norms before and after GPU kernel
  double res0 = local_carre_norme_vect(tab_v0);
  {
    CDoubleArrView Diag = tab_Diag.view_ro();
    DoubleArrView v0 = tab_v0.view_rw();
    Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), ns, KOKKOS_LAMBDA(
                           const int i)
    {
      v0(i) *= Diag(i);
    });
    end_gpu_timer(__KERNEL_NAME__);
  }
  double res = local_carre_norme_vect(tab_v0);
  // Single collective operation
  Process::mp_sum_for_each(res0, res);
  res0 = sqrt(res0);
  res = sqrt(res);

  if (limpr()==1)
    Cout<<"Gmres : initial residual = "<<res0<<finl;
  // See http://stackoverflow.com/questions/3437085/check-nan-number
  // May be could be interesting to implement isnan function somewhere
  if (res0!=res0)
    {
      Cerr << "Nan detected in Solv_Gmres::gmres_local()" << finl;
      Cerr << "Contact TRUST support." << finl;
      Process::exit();
    }
  rec_min = (rec_min<res*epsGMRES) ? res*epsGMRES : rec_min;
  rec_min = (rec_min<rec_max) ? rec_min : rec_max ;
  bool legacy = getenv("TRUST_GMRES_REDUCE_COLLECTIVES") == nullptr;

  // iterations
  for(int it=0; it<nit; it++)
    {
      if (res==0) return 0; // nothing to do
      int nk = nkr;

      //...Orthogonalisation of Arnoldi
      tab_v0 /= res;
      r = 0. ;
      r[0] = res;
      h = 0.;
      for(int j=0; j<nkr; j++)
        {
          tab_v0.echange_espace_virtuel();
          A.multvect_(tab_v0,tab_v1);
          {
            CDoubleArrView Diag = tab_Diag.view_ro();
            DoubleArrView v1 = tab_v1.view_rw();
            Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), ns, KOKKOS_LAMBDA(
                                   const int i)
            {
              v1(i) *= Diag(i);
            });
            end_gpu_timer(__KERNEL_NAME__);
          }
          v[j] = tab_v0;
          tab_v0 = tab_v1 ;
          // Modifie par DJ
          //---------------
          double tem;
          if (legacy)
            {
              // GMRES using classical Gram–Schmidt
              for (int i = 0; i <= j; i++)
                {
                  h(i, j) += mp_prodscal(tab_v0, v[i]);
                  tab_v0.ajoute(-h(i, j), v[i], VECT_REAL_ITEMS);
                }
              tem=mp_norme_vect(tab_v0);
            }
          else
            {
              // Communication-reduced GMRES using classical Gram–Schmidt with reorthogonalization (CGS-2), reducing global MPI collectives but at higher computational cost
              // Compute local dot products
              for (int i = 0; i <= j; i++)
                h_loc[i] = local_prodscal(tab_v0, v[i]);
              Process::mp_sum_for_each_item(h_loc, j+1); // One collective
              // Orthoganalization
              for (int i = 0; i <= j; i++)
                {
                  h(i, j) = h_loc[i]; // Store in Hessenberg
                  tab_v0.ajoute(-h(i, j), v[i], VECT_REAL_ITEMS);
                }
              // Compute correction terms
              for (int i = 0; i <= j; i++)
                dh_loc[i] = local_prodscal(tab_v0, v[i]);
              dh_loc[j + 1] = local_carre_norme_vect(tab_v0);
              Process::mp_sum_for_each_item(dh_loc, j+2); // One collective
              // Accumulate + correct
              for (int i = 0; i <= j; i++)
                {
                  h(i, j) += dh_loc[i];
                  tab_v0.ajoute(-dh_loc[i], v[i], VECT_REAL_ITEMS);
                }
              tem=std::sqrt(dh_loc[j + 1]); // Save one more collective
            }

          h(j+1,j) = tem;
          if(tem<rec_min)
            {
              nk = j+1;
              goto l5;
            }
          tab_v0 /= tem;
        }
      //...Triangularisation
l5:
      for(int i=0; i<nk; i++)
        {
          int im = i+1;
          double tem = 1./sqrt(h(i,i)*h(i,i) + h(im,i)*h(im,i));
          double ccos = h(i,i) * tem;
          double ssin = - h(im,i) * tem;
          for(int j=i; j<nk; j++)
            {
              tem = h(i,j);
              h(i,j) = ccos * tem - ssin * h(im,j);
              h(im,j) =  ssin * tem + ccos * h(im,j);
            }
          r[im] = ssin * r[i];
          r[i] *= ccos;
        }

      //...Solution of linear system
      for(int i=nk-1; i>=0; i--)
        {
          r[i] /= h(i,i);
          for(int i0=i-1; i0>=0; i0--)
            r[i0] -= h(i0,i)* r[i];
        }
      for(int i=0; i<nk; i++)
        tab_x.ajoute(r[i], v[i], VECT_REAL_ITEMS);

      tab_x.echange_espace_virtuel();
      A.multvect_(tab_x,tab_v0);
      tab_v0 *= -1. ;
      tab_v0 += b;

      // calcul du residu sans le precond....
      double res2=mp_norme_vect(tab_v0);
      if ((it>0) && (controle_residu_==1) && (sup_strict(res2,res2_old)))
        {
          Cout << "The Gmres iterative system is stopped after : " << it+1 <<" iterations "<<finl;
          Cout << "since an increase of the residue is detected."<< finl;
          return it;
        }

      res2_old = res2;
      if (limpr()==1)
        Cout<<" - At it = "<< it+1 <<", residu scalar = "<< res2 << finl;

      // Test d'arret sur le residu
      if(res2<rec_min)
        {
          // Ajoute par DJ
          //--------------
          if (limpr()>-1)
            {
              Cout << "Gmres : Number of iterations to reach convergence : " << it+1 << finl;
              double residu_relatif = (res0>0?res2/res0:res2);
              Cout << "Final residue: " << res2 << " ( " << residu_relatif << " )" << finl;
            }
          return it+1;
        }
      // Test d'arret sur le nombre d'iterations max
      else if (it==nit-1)
        {
          if (limpr()>-1)
            {
              Cout << "Gmres : Stopped after "<< it+1 <<" iterations (=nb_it_max)"<< finl;
              double residu_relatif = (res0>0?res2/res0:res2);
              Cout << "Final residue: " << res2 << " ( " << residu_relatif << " )" << finl;
            }
          if (it == (nb_ligne_tot-1))
            {
              Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< finl;
              Cerr << "!!! Gmres stopped after a number of iterations equal to the matrix size. "<< finl;
              Cerr << "!!! Either your matrix is ill-conditioned (try cholesky instead), or your convergence threshold is too low. "<< finl;
              Cerr << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<< finl;
              Process::exit(-1);
            }
          return it+1;
        }

      // Calcul du residu avec preconditionnement
      {
        CDoubleArrView Diag = tab_Diag.view_ro();
        DoubleArrView v0 = tab_v0.view_rw();
        Kokkos::parallel_for(start_gpu_timer(__KERNEL_NAME__), ns, KOKKOS_LAMBDA(
                               const int i)
        {
          v0(i) *= Diag(i);
        });
        end_gpu_timer(__KERNEL_NAME__);
      }
      res = mp_norme_vect(tab_v0);
    }
  return -1;
}


int Solv_Gmres::Gmres(const Matrice_Morse& matrice,
                      const DoubleVect& secmem,
                      DoubleVect& solution)
{
  if (!is_local_gmres)
    return matrice.inverse(secmem, solution, seuil_);
  else
    return gmres_local(matrice,secmem,solution);
}


