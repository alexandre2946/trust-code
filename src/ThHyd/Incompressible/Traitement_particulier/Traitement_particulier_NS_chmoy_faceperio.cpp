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

#include <Traitement_particulier_NS_chmoy_faceperio.h>
#include <LecFicDistribue.h>
#include <Navier_Stokes_std.h>
#include <Schema_Temps_base.h>

Implemente_base_sans_constructeur_ni_destructeur(Traitement_particulier_NS_chmoy_faceperio,"Traitement_particulier_NS_chmoy_faceperio",Traitement_particulier_NS_base);
// XD chmoy_faceperio traitement_particulier_base chmoy_faceperio NO_BRACE non documente
// XD attr bloc bloc_lecture bloc REQ not_set

/*! @brief Prints the object to an output stream.
 *
 * @param is an output stream
 * @return the modified output stream
 */
Sortie& Traitement_particulier_NS_chmoy_faceperio::printOn(Sortie& is) const
{
  return is;
}


/*! @brief Reads the object from an input stream.
 *
 * @param is an input stream
 * @return the modified input stream
 */
Entree& Traitement_particulier_NS_chmoy_faceperio::readOn(Entree& is)
{
  return is;
}

Entree& Traitement_particulier_NS_chmoy_faceperio::lire(Entree& is)
{
  Motcle accouverte = "{" , accfermee = "}" ;
  Motcle motbidon, motlu;
  is >> motbidon ;
  if (motbidon == accouverte)
    {
      Motcles les_mots(1);
      les_mots[0] = "stats";

      is >> motlu;
      while(motlu != accfermee)
        {
          int rang=les_mots.search(motlu);
          switch(rang)
            {
            case 0 :
              {
                // For stats
                is >> temps_deb;     // start time for temporal stats computation
                is >> temps_fin;     // end time for temporal stats computation
                // initialisation for temporal stats computation
                oui_stat=1;          // =1 : computing temporal stats
                Cerr << "Reading time statitics parameters..." << finl;
                Cerr << "Initial time : " << temps_deb << " End time : " << temps_fin << finl;
                // Check whether values specific to spatial stats computations have been provided
                break;
              }
            default :
              {
                Cerr << "Default case..." << finl;
                Cerr << "Possible keywords are "<< les_mots <<" { and }" << finl;
                Cerr << "You read:" << motlu << finl;
                break;
              }
            }
          is >> motlu;
        }
      is >> motlu;
      if (motlu != accfermee)
        {
          Cerr << "Error while reading in Traitement_particulier_NS_canal" << finl;;
          Cerr << "We were expecting a }" << finl;
          exit();
        }
    }
  else
    {
      Cerr << "Error while reading in Traitement_particulier_NS_canal" << finl;
      Cerr << "We were expecting a {" << finl;
      exit();
    }

  return is;
}


void Traitement_particulier_NS_chmoy_faceperio::preparer_calcul_particulier()
{

  if(Objet_U::dimension!=3)
    {
      Cerr << " Traitement_particulier_NS_chmoy_faceperio : not designed for calculations other than 3D " << finl;
      exit();
    }

  if (oui_stat != 0)
    init_calcul_stats();

  double temps = mon_equation->inconnue().temps();
  int Nbfaces,num_face;

  if (temps>temps_deb)
    {
      Nom fichier = "chmoy_face_perio";
      ifstream fic(fichier);
      if (!fic)
        {
          Cerr << " no file : chmoy_face_perio  - resuming calculations not possible" << finl;
          exit();
        }
      else
        {
          LecFicDistribue fic2(fichier);
          fic2 >> Nbfaces;

          for (int i=0; i<Nbfaces; i++)
            fic2 >> num_face >> chmoy_faceperio(num_face,0) >> chmoy_faceperio(num_face,1) >> chmoy_faceperio(num_face,2);
        }
    }
}

void Traitement_particulier_NS_chmoy_faceperio::post_traitement_particulier()
{
  double temps = mon_equation->inconnue().temps();
  double dt = mon_equation->schema_temps().pas_de_temps();
  if (temps>temps_deb && temps<temps_fin)
    calcul_chmoy_faceperio(temps_deb,temps,dt);
}


