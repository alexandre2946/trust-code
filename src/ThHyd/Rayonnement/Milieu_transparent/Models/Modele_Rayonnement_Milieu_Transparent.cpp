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

#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Discretisation_base.h>
#include <Domaine_Cl_dis_base.h>
#include <Frontiere_dis_base.h>
#include <communications.h>
#include <Probleme_base.h>
#include <SFichierBin.h>
#include <EFichierBin.h>
#include <Fluide_base.h>
#include <Interprete.h>
#include <sys/stat.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Modele_Rayonnement_Milieu_Transparent, "Modele_Rayonnement_Milieu_Transparent", Objet_U);

extern "C" {
  void F77DECLARE(dgemv)(char *trans, integer *M, integer *N, double *alpha, double *const A, integer *lda, const double *dx, integer *incx, double *beta, double *dy, integer *incy);
}

Sortie& Modele_Rayonnement_Milieu_Transparent::printOn(Sortie& os) const
{
  return os;
}

Entree& Modele_Rayonnement_Milieu_Transparent::readOn(Entree& is)
{
  Cerr << "Reading params of " << que_suis_je() << finl;
  Nom fichier_face_rayo, fichier_fij;

  Param param(que_suis_je());
  param.ajouter("fichier_face_rayo", &fichier_face_rayo, Param::REQUIRED);
  param.ajouter("fichier_fij", &fichier_fij, Param::REQUIRED);
  param.ajouter("fichier_matrice", &nom_fic_mat_ray_inv_);
  param.ajouter("relaxation", &relaxation_);
  param.ajouter("nom_pb_rayonnant", &nom_pb_rayonnant_);
  param.ajouter_flag("format_binaire", &fic_mat_ray_inv_bin_);
  param.lire_avec_accolades_depuis(is);

  if (fic_mat_ray_inv_bin_)
    lire_fichiers(fichier_face_rayo, fichier_fij, nom_fic_mat_ray_inv_);
  else
    lire_fichiers(fichier_face_rayo, fichier_fij);

  // si on a lu le nom du pb rayonnant on indique au pb qu'il est rayonnant
  if (nom_pb_rayonnant_ != "??")
    {
      Probleme_base& pb = ref_cast(Probleme_base, interprete().objet(nom_pb_rayonnant_));
      if (sub_type(Fluide_base, pb.milieu()))
        {
          Fluide_base& fluide = ref_cast(Fluide_base, pb.milieu());
          fluide.fixer_type_rayo();
        }
      else
        Process::exit("Le nom du pb rayonnant n'est pas un pb fluide !! \n");
    }
  return is;
}

void Modele_Rayonnement_Milieu_Transparent::lire_fichiers(Nom& fich_faces_rayo, Nom& fich_fij, Nom& fich_mat)
{
  Cerr << "Modele_Rayonnement_Milieu_Transparent::lire_fichiers" << finl;
  Cerr << "fichier_face_rayo = " << fich_faces_rayo << finl;
  Cerr << "fichier_fij = " << fich_fij << finl;
  Cerr << "fichier_matrice = " << fich_mat << finl;

  struct stat f, e, m;

  // 3 fichiers a tester :
  const char *facteur_file = fich_fij;
  const char *emissivite_file = fich_faces_rayo;
  const char *matrice_inverse_file = fich_mat;

  if (stat(facteur_file, &f))
    {
      Cerr << facteur_file << " doesn't exist." << finl;
      Process::exit();
    }

  if (stat(emissivite_file, &e))
    {
      Cerr << emissivite_file << " doesn't exist." << finl;
      Process::exit();
    }

  if (stat(matrice_inverse_file, &m))
    lire_matrice_inv_ = false;
  else
    lire_matrice_inv_ = true;

  // Teste si facteur_file ou emissivite_file plus recent que matrice_inverse_file
  if (lire_matrice_inv_ && (f.st_mtime > m.st_mtime || e.st_mtime > m.st_mtime))
    lire_matrice_inv_ = false;

  // Pour l'instant, on se contente de relire les deux premiers fichiers. Le dernier nom (nom3)
  // sert de test pour indiquer qu'il faut inverser la matrice
  lire_fichiers(fich_faces_rayo, fich_fij);
}

void Modele_Rayonnement_Milieu_Transparent::lire_fichiers(Nom& fich_faces_rayo, Nom& fich_fij)
{
  Cerr << "Modele_Rayonnement_Milieu_Transparent::lire_fichiers" << finl;
  Cerr << "fichier_face_rayo = " << fich_faces_rayo << finl;
  Cerr << "fichier_fij = " << fich_fij << finl;

  EFichier fic1(fich_faces_rayo);
  EFichier fic2(fich_fij);

  // lecture du nombre de face rayonnante et des faces rayonnantes dans fic1
  Cerr << "Lecture du fichier : " << fich_faces_rayo << finl;

  fic1 >> nb_faces_totales_ >> nb_faces_rayonnantes_;
  Cerr << "Vous avez defini " << nb_faces_rayonnantes_ << " faces rayonnantes sur " << nb_faces_totales_ << " faces en tout" << finl;
  les_faces_rayonnantes_.dimensionner(nb_faces_totales_);

  les_flux_radiatifs_.resize(nb_faces_rayonnantes());

  int irayo = 0;
  int jrayo = 0;
  for (int i = 0; i < nb_faces_totales_; i++)
    {
      fic1 >> les_faces_rayonnantes_[i];
      if (les_faces_rayonnantes_[i].emissivite() != -1)
        {
          if (les_faces_rayonnantes_[i].emissivite() < 0.)
            {
              Cerr << "erreur dans " << fich_faces_rayo << finl;
              Cerr << "emissivite <0 et != -1 :" << les_faces_rayonnantes_[i].emissivite() << finl;
              Process::exit();
            }
          if (irayo >= nb_faces_rayonnantes_)
            {
              Cerr << "erreur dans " << fich_faces_rayo << finl;
              Cerr << "nb_faces_rayonnantes utilisateur " << nb_faces_rayonnantes_ << finl;
              Cerr << "nb_faces_rayonnantes deduites des emissivites " << irayo + 1 << finl;
              // G.F.
              Cerr << "on a corrige emissivite=0 ne veut pas dire pas de rayonnement" << finl;
              Cerr << "Il faut mettre maintenant l'emissivite a -1" << finl;
              Cerr << "pour ne pas tenir compte de certains bords " << finl;
              Process::exit();
            }
          irayo++;
          jrayo++;

        }
      else
        jrayo++;
    }

  if (irayo != nb_faces_rayonnantes_)
    {
      Cerr << "erreur dans " << fich_faces_rayo << finl;
      Cerr << "nb_faces_rayonnantes utilisateur " << nb_faces_rayonnantes_ << finl;
      Cerr << "nb_faces_rayonnantes deduites des emissivites " << irayo << finl;
      Process::exit();
    }

  if (!lire_matrice_inv_) // lecture fichier fij dans ce cas !!
    {
      Cerr << "Lecture du fichier : " << fich_fij << finl;

      fic2 >> ordre_mat_forme_;
      if (ordre_mat_forme_ != nb_faces_totales_)
        {
          Cerr << "Le nombre de faces dans la matrice des facteurs de formes est different de celui du fichier des faces rayonnantes" << finl;
          Cerr << "Verifier votre fichier " << fich_fij << finl;
          Process::exit();
        }

      Cerr << "l'ordre de la matrice est " << ordre_mat_forme_ << finl;
      // Dimensionnement de la "matrice" des facteurs de formes.
      les_facteurs_de_forme_.resize(nb_faces_rayonnantes(), nb_faces_rayonnantes());
      irayo = 0;
      jrayo = 0;
      double poub;
      for (int ii = 0; ii < nb_faces_totales_; ii++)
        {
          if (les_faces_rayonnantes_[ii].emissivite() != -1)
            {
              for (int j = 0; j < nb_faces_totales_; j++)
                if (les_faces_rayonnantes_[j].emissivite() != -1)
                  {
                    fic2 >> les_facteurs_de_forme_(irayo, jrayo);
                    jrayo++;
                  }
                else
                  fic2 >> poub;
              irayo++;
            }
          else
            {
              for (int j = 0; j < nb_faces_totales_; j++)
                fic2 >> poub;
            }
          jrayo = 0;
        }
    }
  else
    {
      Cerr << "On ne lit pas la matrice des facteurs de forme puisque l'on va " << finl;
      Cerr << "lire directement la matrice inverse qui a deja ete calculee" << finl;
      Cerr << "lors d'un calcul precedent." << finl;
    }
  Cerr << "La lecture des fichiers du prepro est terminee. Fichiers corrects." << finl;
}

void Modele_Rayonnement_Milieu_Transparent::calculer_temperatures()
{
  for (int i = 0; i < nb_faces_totales(); i++)
    if (les_faces_rayonnantes_[i].emissivite() != -1)
      les_faces_rayonnantes_[i].calculer_temperature();
}

void Modele_Rayonnement_Milieu_Transparent::calculer_flux_radiatifs()
{
  int jrayo = 0;

  if (processeur_rayonnant() != -1)
    {
      assert(me() == 0);
      ArrOfDouble secmem(nb_faces_rayonnantes());

      // Remplissage du second membre.
      int irayo = 0;

      for (int i = 0; i < nb_faces_totales(); i++)
        {
          const Face_Rayonnante& Facei = les_faces_rayonnantes_[i];
          if (Facei.emissivite() != -1)
            {
              //                    secmem(irayo) = Facei.emissivite()*SIGMA*(pow(Facei.T_face_rayo(),4));
              secmem[irayo] = (pow(Facei.T_face_rayo(), 4));
              irayo++;
            }
        }

      int ii, jj;
      // La matrice_rayo a ete inversee au debut du calcule dans Modele_Rayonnement_Milieu_Transparent::preparer_calcul.
      // Il ne reste ici qu'a calculer le produit du second membre avec la matrice inverse.
      if (1 == 1)
        {
          for (ii = 0; ii < nb_faces_rayonnantes(); ii++)
            {
              les_flux_radiatifs_(ii) = 0.;
              for (jj = 0; jj < nb_faces_rayonnantes(); jj++)
                les_flux_radiatifs_(ii) += matrice_rayo_(ii, jj) * secmem[jj];
            }
        }
      else
        {
          double alpha = 1, beta = 0;
          integer n = nb_faces_rayonnantes(), inc = 1;

          F77NAME(dgemv)((char*) "T", &n, &n, &alpha, matrice_rayo_.addr(), &n, secmem.addr(), &inc, &beta, les_flux_radiatifs_.addr(), &inc);
        }
    }
  envoyer_broadcast(les_flux_radiatifs_, 0);

  for (int i = 0; i < nb_faces_totales_; i++)
    {
      Face_Rayonnante& Facei = les_faces_rayonnantes_[i];
      if (Facei.emissivite() != -1)
        {
          Facei.mettre_a_jour_flux_radiatif(les_flux_radiatifs_(jrayo));
          jrayo++;
        }
    }
}

void Modele_Rayonnement_Milieu_Transparent::imprimer_flux_radiatifs(Sortie& os) const
{
  if (Process::me()) return; /* seulement maitre qui imprime ! */

  Nom fichier1(nom_du_cas());
  fichier1 += "_Flux_radiatif.out";
  Nom fichier2(nom_du_cas());

  Nom espace = "\t\t";

  fichier2 += "_Temperature_rayonnante.out";
  if (deja_imprime_ == 0)
    {
      SFichier os1(fichier1);
      SFichier os2(fichier2);
      deja_imprime_ = 1;
      os1 << "# Impression sur les bords rayonnants de l'equation du transfert radiatif" << finl << "# Flux radiatif [W]" << finl << "# Bord:";
      os2 << "# Impression sur les bords rayonnants de l'equation du transfert radiatif" << finl << "# Temperature de bord moyenne en K" << finl << "# Bord:";
      for (int i = 0; i < nb_faces_totales(); i++)
        if (les_faces_rayonnantes_[i].emissivite() != -1)
          {

            os1 << espace << les_faces_rayonnantes_[i].nom_bord_rayo_lu();
            os2 << espace << les_faces_rayonnantes_[i].nom_bord_rayo_lu();
          }
      os1 << espace << "Total" << finl << "# Temps" << finl;
      os2 << finl << "# Temps" << finl;
    }
  SFichier os1(fichier1, ios::app);
  SFichier os2(fichier2, ios::app);
  //const int& precision=sch.precision_impr();
  int precision = 3;
  os1.precision(precision);
  os1.setf(ios::scientific);
  os2.precision(precision);
  os2.setf(ios::scientific);

  double flux_tot = 0;

  os1 << temps_;
  os2 << temps_;
  for (int i = 0; i < nb_faces_totales(); i++)
    if (les_faces_rayonnantes_[i].emissivite() != -1)
      flux_tot += les_faces_rayonnantes_[i].imprimer_flux_radiatif(os, os1, os2);

  os1 << "\t" << flux_tot << finl;
  os2 << finl;
  os << "Bilan flux radiatifs : " << flux_tot << " W" << finl;
}

double Modele_Rayonnement_Milieu_Transparent::flux_radiatif(int num_face) const
{
  if (corres_.size() == 0)
    {
      // on recupere le domaine
      int i0 = 0;
      // on cherche la premiere cond_lim rayo
      while (((les_faces_rayonnantes_[i0].ensembles_faces_bord(0).nb_faces_bord() == 0) || (les_faces_rayonnantes_[i0].emissivite() == -1)) && (i0 < nb_faces_totales()))
        i0++;
      int nbre_face_de_bord;

      if (i0 == nb_faces_totales())
        nbre_face_de_bord = 0;
      else
        nbre_face_de_bord = les_faces_rayonnantes_[i0].ensembles_faces_bord(0).la_cl_base().domaine_Cl_dis().nb_faces_Cl();

      corres_.resize(nbre_face_de_bord);
      corres_ = -1;
      for (int i = 0; i < nb_faces_totales(); i++)
        if (les_faces_rayonnantes_[i].emissivite() != -1)
          {
            const Ensemble_Faces_base& ensemble = les_faces_rayonnantes_[i].ensembles_faces_bord(0);
            if (ensemble.nb_faces_bord() != 0)
              {
                //const IntVect&  num_face_ens= ensemble.Table_faces ();
                const Frontiere& la_front = ensemble.la_cl_base().frontiere_dis().frontiere();
                int ndeb = la_front.num_premiere_face();
                int nbfaces = la_front.nb_faces();
                for (int face = 0; face < nbfaces; face++)
                  {
                    if (ensemble.contient(face))
                      {
                        int nglob = ndeb + face;
                        if (corres_[nglob] != -1)
                          {
                            Cerr << me() << "la face " << nglob << " semble etre contenu par les faces_rayonnantes " << i << " et " << corres_[nglob] << finl;
                            Cerr << me() << les_faces_rayonnantes_[i].nom_bord_rayo_lu() << " " << les_faces_rayonnantes_[corres_[nglob]].nom_bord_rayo_lu() << finl;
                            Cerr << ensemble.nb_faces_bord() << finl;
                            Process::exit();
                          }
                        else
                          corres_[nglob] = i;
                      }
                  }
              }
          }

      for (int i = 0; i < nbre_face_de_bord; i++)
        if (corres_[i] == -1)
          Cout << "Face " << i << " sans groupes " << finl;

    }
  return les_faces_rayonnantes_[corres_[num_face]].flux_radiatif();
}

void Modele_Rayonnement_Milieu_Transparent::preparer_calcul()
{
  if (je_suis_maitre())
    {
      // dimensionnement de la matrice de rayonnement.
      int irayo = 0;
      int jrayo = 0;

      matrice_rayo_.resize(nb_faces_rayonnantes(), nb_faces_rayonnantes());

      if (!lire_matrice_inv_)
        {
          for (int i = 0; i < nb_faces_totales_; i++)
            {
              for (int j = 0; j < nb_faces_totales_; j++)
                {
                  if (les_faces_rayonnantes_[i].emissivite() != -1)
                    {
                      if ((jrayo < nb_faces_rayonnantes()))
                        {
                          if (irayo == jrayo)
                            matrice_rayo_(irayo, jrayo) = 1 - (1 - les_faces_rayonnantes_[i].emissivite()) * les_facteurs_de_forme_(irayo, jrayo);
                          else
                            matrice_rayo_(irayo, jrayo) = (les_faces_rayonnantes_[i].emissivite() - 1) * les_facteurs_de_forme_(irayo, jrayo);
                        }
                      jrayo++;
                    }
                }

              if (les_faces_rayonnantes_[i].emissivite() != -1)
                irayo++;

              jrayo = 0;
            }
        }
      else
        {
          // Puisque l'on va lire la matrice de rayonnement inverse dans un fichier, il n'est pas necessaire de remplire matrice_rayo
        }

      Nom version("version_2beta3");

      // INVERSION DE LA MATRICE UNE FOIS POUR TOUT
      if (!lire_matrice_inv_)
        {
          Cerr << "Inversion de la matrice de rayonnement au debut du calcul" << finl;
          // On commence par decomposer la matrice de rayonnement en une decomposition LU.
          DoubleVect sol_tmp(nb_faces_rayonnantes());
          {
            IntVect index(nb_faces_rayonnantes());
            DoubleTab lu_dec(nb_faces_rayonnantes(), nb_faces_rayonnantes());
            DoubleVect secmem_tmp(nb_faces_rayonnantes());

            int cvg = matrice_rayo_.decomp_LU(nb_faces_rayonnantes(), index, lu_dec);

            // Puis on inverse la matrice_rayo
            if (cvg == 1)
              {
                for (jrayo = 0; jrayo < nb_faces_rayonnantes(); jrayo++)
                  {
                    for (irayo = 0; irayo < nb_faces_rayonnantes(); irayo++)
                      secmem_tmp(irayo) = 0.;
                    secmem_tmp(jrayo) = 1.;

                    lu_dec.resoud_LU(nb_faces_rayonnantes(), index, secmem_tmp, sol_tmp);
                    // On recopie le resultat de la resolution ci-dessus dans la colonne jrayo de matrice_rayo
                    for (irayo = 0; irayo < nb_faces_rayonnantes(); irayo++)
                      matrice_rayo_(irayo, jrayo) = sol_tmp(irayo);
                  }
                // matrice_rayo contient maintenant l'inverse de la matrice_rayo initiale
              }
          }

          {
            // on change la matrice;
            les_facteurs_de_forme_ *= -1;
            for (jrayo = 0; jrayo < nb_faces_rayonnantes(); jrayo++)
              les_facteurs_de_forme_(jrayo, jrayo) += 1;

            for (jrayo = 0; jrayo < nb_faces_rayonnantes(); jrayo++)
              {
                // on sauve la colonne
                for (irayo = 0; irayo < nb_faces_rayonnantes(); irayo++)
                  sol_tmp(irayo) = matrice_rayo_(irayo, jrayo);
                //calcul de (I-Fij)*M-1
                for (irayo = 0; irayo < nb_faces_rayonnantes(); irayo++)
                  {
                    double res = 0;
                    for (int krayo = 0; krayo < nb_faces_rayonnantes(); krayo++)
                      res += les_facteurs_de_forme_(irayo, krayo) * sol_tmp(krayo);
                    matrice_rayo_(irayo, jrayo) = res;
                  }
              }

            // on continue en multipliant par sigma delta_i_j emissivite(j)
            jrayo = 0;
            for (int j = 0; j < nb_faces_totales(); j++)
              {
                const Face_Rayonnante& Facej = les_faces_rayonnantes_[j];
                if (Facej.emissivite() != -1)
                  {
                    for (irayo = 0; irayo < nb_faces_rayonnantes(); irayo++)
                      matrice_rayo_(irayo, jrayo) *= Facej.emissivite() * SIGMA_;
                    jrayo++;
                  }
              }
          }

          // On a fini on peut vider Fij
          les_facteurs_de_forme_.resize(0, 0);

          // Impression "jolie" de la matrice dans un fichier nom_fic_mat_ray_inv_

          if (nom_fic_mat_ray_inv_ != "??")
            {
              if (!fic_mat_ray_inv_bin_)
                {
                  SFichier fic(nom_fic_mat_ray_inv_);
                  fic.setf(ios::scientific);
                  fic.precision(10);
                  fic << version << finl << nb_faces_rayonnantes() << finl;
                  for (irayo = 0; irayo < nb_faces_rayonnantes(); irayo++)
                    {
                      for (jrayo = 0; jrayo < nb_faces_rayonnantes(); jrayo++)
                        fic << "   " << matrice_rayo_(irayo, jrayo);
                      fic << "  " << finl;
                    }
                }
              else
                // Ecriture le la matrice inverse en format binaire
                {
                  SFichierBin fic_bin(nom_fic_mat_ray_inv_);
                  fic_bin << version << finl << matrice_rayo_;
                  fic_bin.flush();
                }
            }
        }
      else
        {
          if (!fic_mat_ray_inv_bin_)
            {
              Cerr << "Lecture du fichier ASCI " << nom_fic_mat_ray_inv_ << finl;
              EFichier fic(nom_fic_mat_ray_inv_);

              int i, j;
              Cerr << "Lecture du fichier : " << nom_fic_mat_ray_inv_ << finl;
              Nom toto;
              fic >> toto;
              if (toto != version)
                {
                  Cerr << "Il faut detruire votre fichier " << nom_fic_mat_ray_inv_ << " la matrice stockee a changee." << finl;
                  Process::exit();
                }
              fic >> ordre_mat_forme_;
              if (ordre_mat_forme_ != nb_faces_rayonnantes())
                {
                  Cerr << "L'ordre de la matrice inverse de rayonnement est different" << finl;
                  Cerr << "du nombre de faces rayonnantes" << finl;
                  Cerr << "Verifiez votre fichier " << nom_fic_mat_ray_inv_ << finl;
                  Process::exit();
                }
              Cerr << "l'ordre de la matrice est " << ordre_mat_forme_ << finl;

              for (i = 0; i < ordre_mat_forme_; i++)
                for (j = 0; j < ordre_mat_forme_; j++)
                  fic >> matrice_rayo_(i, j);
            }
          else
            {
              Cerr << "Lecture du fichier BINAIRE " << nom_fic_mat_ray_inv_ << finl;
              EFichierBin fic_bin(nom_fic_mat_ray_inv_);

              Cerr << "Lecture du fichier : " << nom_fic_mat_ray_inv_ << finl;
              Nom toto;
              fic_bin >> toto;
              if (toto != version)
                {
                  Cerr << "Il faut detruire votre fichier " << nom_fic_mat_ray_inv_ << " la matrice stockee a changee." << finl;
                  Process::exit();
                }
              fic_bin >> matrice_rayo_;
              ordre_mat_forme_ = matrice_rayo_.dimension(0);
              if (ordre_mat_forme_ != nb_faces_rayonnantes())
                {
                  Cerr << "L'ordre de la matrice inverse de rayonnement est different" << finl;
                  Cerr << "du nombre de faces rayonnantes" << finl;
                  Cerr << "Verifiez votre fichier " << nom_fic_mat_ray_inv_ << finl;
                  Process::exit();
                }
              Cerr << "l'ordre de la matrice est " << ordre_mat_forme_ << finl;
            }
          Cerr << "La lecture de la matrice de rayonnement inverse dans le fichier " << finl;
          Cerr << nom_fic_mat_ray_inv_ << " est terminee. Fichier correct " << finl;
        }
    }
}

void Modele_Rayonnement_Milieu_Transparent::mettre_a_jour(double temps)
{
  temps_ = temps;
  calculer_temperatures();
  calculer_flux_radiatifs();
}

void Modele_Rayonnement_Milieu_Transparent::discretiser(const Discretisation_base& dis, const Domaine& domaine)
{
  for (int i = 0; i < nb_faces_totales(); i++)
    {
      Face_Rayonnante& face_rayo = face_rayonnante(i);
      Nom type("Ensemble_Faces_");
      Nom discr = dis.que_suis_je();
      if (discr == "VEFPreP1B")
        discr = "VEF";

      type += discr;
      for (int j = 0; j < face_rayo.nb_ensembles_faces(); j++)
        {
          if (face_rayo.nom_bord_rayo() != face_rayo.nom_bord_rayo_lu())
            {
              Ensemble_Faces_base& faces_j = face_rayo.ensembles_faces_bord(j);
              faces_j.lire(face_rayo.nom_bord_rayo_lu(), face_rayo.nom_bord_rayo(), domaine);
            }
        }
    }
}
