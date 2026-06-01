/****************************************************************************
* Copyright (c) 2025, CEA
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

#include <Traitement_particulier_NS_THI_new.h>
#include <MD_Vector_base.h>
#include <Domaine.h>

Implemente_base(Traitement_particulier_NS_THI_new,"Traitement_particulier_NS_THI_new",Traitement_particulier_NS_base);

Sortie& Traitement_particulier_NS_THI_new::printOn(Sortie& is) const
{
  return is;
}

Entree& Traitement_particulier_NS_THI_new::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_THI_new::lire(Entree& is)
{
  Motcle accouverte = "{" , accfermee = "}" ;
  Motcle valec="val_Ec";
  Motcle facon="facon_init";
  Motcle motbidon, motlu;
  is >> motbidon ;
  if (motbidon == accouverte)
    {
      Motcles les_mots(2);
      {
        les_mots[0] = "init_Ec";
        les_mots[1] = "calc_spectre";
      }
      {
        is >> motlu;
        while(motlu != accfermee)
          {
            int rang=les_mots.search(motlu);
            switch(rang)
              {
              case 0 :
                {
                  is >> init;  // init_Ec = 1
                  if (init!=0)
                    {
                      is >> motlu;
                      if (motlu==valec)
                        {
                          is >> Ec_init; // read the value of Ec_init
                          Cerr << "With initialization of kinetic energy Ec_init= " << Ec_init << finl;
                          is >> motlu;
                          if (motlu==facon)
                            {
                              is >> fac_init; // Which value is the initialization based on??
                              Cerr << "With initialization of kinetic energy on Ecspat (init_fac==0) or Ecspec (init_fac==1) : fac_init=" << fac_init << finl;
                            }
                        }
                      else
                        {
                          Cerr << "Error while reading Traitement_particulier_NS_THI_new_VDF" << finl;
                          Cerr << "The only possible keyword here is: val_Ec" << finl;
                          Cerr << "You read:" << motlu << finl;
                          exit();
                        }
                    }
                  break;
                }
              case 1 :
                {
                  is >> oui_calc_spectre;
                  if (oui_calc_spectre!=0) Cerr << "Computing spectra." << finl;
                  break;
                }

              default :
                {
                  Cerr << "Error while reading Traitement_particulier_NS_THI_new";
                  Cerr << "Possible keywords are: init_Ec, calc_spectre, { and }" << finl;
                  Cerr << "You read:" << motlu << finl;
                  exit();
                  break;
                }
              }
            is >> motlu;
          }
        is >> motlu;
        if (motlu != accfermee)
          {
            Cerr << "Error while reading Traitement_particulier_NS_THI_new";
            Cerr << "We expected a }" << finl;
            exit();
          }
      }
    }
  else
    {
      Cerr << "Error while reading Traitement_particulier_NS_THI_new";
      Cerr << "We expected a {" << finl;
      exit();
    }
  return is;
}


void Traitement_particulier_NS_THI_new::preparer_calcul_particulier()
{
  if ((oui_calc_spectre != 0)||(fac_init!=0))
    init_calc_spectre();
  if (init==1) renorm_Ec();
  if (oui_calc_spectre != 0)
    calcul_spectre();
}

void Traitement_particulier_NS_THI_new::post_traitement_particulier()
{
  if (oui_calc_spectre != 0)
    //    renorm_Ec();
    calcul_spectre();
}

void Traitement_particulier_NS_THI_new::en_cours_de_resolution(int nb_op, DoubleTab& u, DoubleTab& u_av, double dt)
{
  // it is actually never called, cf. below (Caroline)
  if (oui_transf == 1)
    {
      //      calcul_spectre_operateur( nb_op, u, u_av,dt);

      // Patrick

      Cerr << " finally Traitement_particulier_NS_THI_new::en_cours_de_resolution  is indeed called (?) " << finl;
      exit();
    }
  return;
}

int& Traitement_particulier_NS_THI_new::calcul_nb_som_dir(const Domaine& domaine)
{
  const char* methode_actuelle="Traitement_particulier_NS_THI_new::calcul_nb_som_dir";

  // Used to compute the number of common vertices in parallel
  double nb_som = static_cast<double>(domaine.md_vector_sommets()->nb_items_seq_tot());
  // Sum over all processors

  double nb=pow(nb_som*1.,1./3.);
  nb_som_dir = (int)(nb);
  if (nb_som_dir*nb_som_dir*nb_som_dir != nb_som)
    {
      nb_som_dir=nb_som_dir+1;
      if (nb_som_dir*nb_som_dir*nb_som_dir != nb_som)
        msg_erreur_maillage(methode_actuelle);
    }
  nb_som_dir=nb_som_dir-1;
  return nb_som_dir;
}

void Traitement_particulier_NS_THI_new::msg_erreur_maillage(const char* methode_actuelle)
{
  if (je_suis_maitre())
    {
      Cerr << finl;
      Cerr << "Problem in " << methode_actuelle << " :" << finl;
      Cerr << "Your mesh does not seem to have the same number of nodes" << finl;
      Cerr << "in all directions. Check your data file... " << finl << finl;
      // Abort for parallel cases!
      abort();
    }
}
