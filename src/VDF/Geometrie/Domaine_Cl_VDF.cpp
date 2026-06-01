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

#include <Dirichlet_entree_fluide_leaves.h>
#include <Sortie_libre_pression_imposee.h>
#include <Dirichlet_paroi_defilante.h>
#include <Dirichlet_paroi_fixe.h>
#include <Champ_Face_VDF.h>
#include <Domaine_Cl_VDF.h>
#include <Champ_P0_VDF.h>
#include <Domaine_VDF.h>
#include <Periodique.h>
#include <Option_VDF.h>
#include <Navier.h>
#include <Debog.h>

Implemente_instanciable(Domaine_Cl_VDF,"Domaine_Cl_VDF",Domaine_Cl_dis_base);

Sortie& Domaine_Cl_VDF::printOn(Sortie& os) const
{
  Domaine_Cl_dis_base::printOn(os);
  return os << "type_arete_bord_ : " << type_arete_bord_ << finl << "num_Cl_face_ : " << num_Cl_face_ << finl;
}

Entree& Domaine_Cl_VDF::readOn(Entree& is) { return Domaine_Cl_dis_base::readOn(is); }

void Domaine_Cl_VDF::associer(const Domaine_dis_base& dom_dis)
{
  const Domaine_VDF& domaine_vdf = ref_cast(Domaine_VDF, dom_dis);
  type_arete_bord_.resize(domaine_vdf.nb_aretes_bord());
  type_arete_coin_.resize(domaine_vdf.nb_aretes_coin());
  num_Cl_face_.resize(domaine_vdf.nb_faces_bord());
}

void Domaine_Cl_VDF::completer(const Domaine_dis_base& un_domaine_dis)
{
  if (sub_type(Domaine_VDF,un_domaine_dis))
    {
      const Domaine_VDF& le_dom_VDF = ref_cast(Domaine_VDF,un_domaine_dis);

      //  Fill the integer array type_arete_bord_ and the intermediate array les_faces_Cl.
      //  les_faces_Cl gives the boundary condition for each boundary face using the following conventions:
      //    0 for a wall condition
      //    1 for a fluid inlet or outlet condition ("fluid" face)
      //    2 for a Navier condition (symmetry/slip wall)
      //    3 for a periodicity condition
      //    0 for any other BC.
      //  Only boundary conditions that affect the momentum equation are considered.

      int nb_aretes_bord = le_dom_VDF.nb_aretes_bord();
      IntVect les_faces_Cl;
      le_dom_VDF.creer_tableau_faces_bord(les_faces_Cl);

      // Loop over boundary conditions to fill les_faces_Cl:
      for (int n_bord = 0; n_bord < le_dom_VDF.nb_front_Cl(); n_bord++)
        {
          // for each boundary condition, determine its type
          const Cond_lim_base& la_cl = les_conditions_limites_[n_bord].valeur();

          int numero_cl = 0;

          if ((sub_type(Dirichlet_paroi_fixe, la_cl)) || (sub_type(Dirichlet_paroi_defilante, la_cl)))
            numero_cl = 0;
          else if ((sub_type(Dirichlet_entree_fluide, la_cl)) || (sub_type(Neumann_sortie_libre, la_cl)))
            numero_cl = 1;
          else if (sub_type(Navier, la_cl)) // (symmetry/slip wall)
            numero_cl = 2;
          else if (sub_type(Periodique, la_cl))
            numero_cl = 3;

          const Frontiere& fr = la_cl.frontiere_dis().frontiere();
          const int ndeb = fr.num_premiere_face(), nfin = ndeb + fr.nb_faces();
          for (int i = ndeb; i < nfin; i++)
            {
              les_faces_Cl[i] = numero_cl;
              num_Cl_face_[i] = n_bord;
            }
        }

      les_faces_Cl.echange_espace_virtuel();

      // Loop over boundary edges to fill type_arete_bord_

      type_arete_bord_ = TypeAreteBordVDF::VIDE; // initialize

      int face1, face2, rang1, rang2;
      int ndeb = le_dom_VDF.premiere_arete_bord(), nfin = ndeb + nb_aretes_bord;
      int num_arete_;

      int decal_virt = le_dom_VDF.nb_faces();
      int sz_faces_Cl = les_faces_Cl.size_reelle();
      const ArrOfInt& ind_faces_virt_bord = le_dom_VDF.domaine().ind_faces_virt_bord();

      for (int num_arete = ndeb; num_arete < nfin; num_arete++)
        {
          num_arete_ = num_arete - ndeb;
          face1 = le_dom_VDF.Qdm(num_arete, 0);
          face2 = le_dom_VDF.Qdm(num_arete, 1);
          rang1 = face1;
          rang2 = face2;

          if (rang1 >= sz_faces_Cl)
            rang1 = ind_faces_virt_bord[rang1 - decal_virt];

          if (rang2 >= sz_faces_Cl)
            rang2 = ind_faces_virt_bord[rang2 - decal_virt];

          if (les_faces_Cl[rang1] == 0) // wall
            {
              if (les_faces_Cl[rang2] == 0) // wall
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::PAROI_PAROI;
              else if (les_faces_Cl[rang2] == 1) // fluid inlet/outlet
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::PAROI_FLUIDE;
              else if (les_faces_Cl[rang2] == 2) // navier
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::PAROI_NAVIER;
              else
                {
                  Cerr << "Processing an edge separating two faces of the same orientation: " << finl;
                  Cerr << "one of these faces has a wall-type boundary condition " << finl;
                  Cerr << "and the other has a Periodicity condition" << finl;
                  Cerr << "TRUST cannot handle this situation" << finl;
                  Process::Journal() << "ERREUR Faces_CL" << finl;
                  Process::Journal() << "face1:" << face1 << " , face2:" << face2 << finl;
                  Process::Journal() << "rang1:" << rang1 << " , rang2:" << rang2 << finl;
                  Process::Journal() << "les_faces_Cl[rang1]:" << les_faces_Cl[rang1] << " , les_faces_Cl[rang2]:" << les_faces_Cl[rang2] << finl;
                  Process::exit();
                }
            }
          else if (les_faces_Cl[rang1] == 1) // fluid inlet/outlet
            {
              if (les_faces_Cl[rang2] == 0) // wall
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::PAROI_FLUIDE;
              else if (les_faces_Cl[rang2] == 1) // fluid
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::FLUIDE_FLUIDE;
              else if (les_faces_Cl[rang2] == 2) // navier
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::FLUIDE_NAVIER;
              else
                {
                  Cerr << "Processing an edge separating two faces of the same orientation: " << finl;
                  Cerr << "one of these faces has a fluid boundary condition " << finl;
                  Cerr << "and the other has a Periodicity condition" << finl;
                  Cerr << "TRUST cannot handle this situation" << finl;
                  Process::exit();
                }
            }
          else if (les_faces_Cl[rang1] == 2) // navier
            {
              if (les_faces_Cl[rang2] == 0) // wall
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::PAROI_NAVIER;
              else if (les_faces_Cl[rang2] == 1) // inlet
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::FLUIDE_NAVIER;
              else if (les_faces_Cl[rang2] == 2) // navier
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::NAVIER_NAVIER;
              else
                {
                  Cerr << "Processing an edge separating two faces of the same orientation: " << finl;
                  Cerr << "one of these faces has a Symmetry-type boundary condition " << finl;
                  Cerr << "and the other has a Periodicity condition" << finl;
                  Cerr << "TRUST cannot handle this situation" << finl;
                  Process::exit();
                }
            }
          else if (les_faces_Cl[rang1] == 3)
            {
              if (les_faces_Cl[rang2] == 3)
                type_arete_bord_[num_arete_] = TypeAreteBordVDF::PERIO_PERIO;
              else
                {
                  Cerr << "On traite une arete qui separe deux faces de meme orientation : " << finl;
                  Cerr << "l'une de ces faces porte une condition limite de type Periodicite " << finl;
                  Cerr << "et l'autre porte une condition d'un autre type" << finl;
                  Cerr << "TRUST ne sait pas traiter cette situation" << finl;
                  Process::exit();
                }
            }
        }

      // MODIFS CA : 21/09/99
      // Loop over corner edges to fill type_arete_coin_

      type_arete_coin_= TypeAreteCoinVDF::VIDE; // initialize

      ndeb = le_dom_VDF.premiere_arete_coin();
      nfin = ndeb + le_dom_VDF.nb_aretes_coin();
      int fac1, fac2, fac3, fac4;
      ArrOfInt f(2);
      for (int num_arete = ndeb; num_arete < nfin; num_arete++)
        {
          num_arete_ = num_arete - ndeb;

          // Find the 2 values different from -1
          fac1 = le_dom_VDF.Qdm(num_arete_, 0);
          fac2 = le_dom_VDF.Qdm(num_arete_, 1);
          fac3 = le_dom_VDF.Qdm(num_arete_, 2);
          fac4 = le_dom_VDF.Qdm(num_arete_, 3);

          f = -2;
          int i = 0;
          if (fac1 != -1)
            {
              f(i) = fac1;
              i++;
            }
          if (fac2 != -1)
            {
              f(i) = fac2;
              i++;
            }
          if (fac3 != -1)
            {
              f(i) = fac3;
              i++;
            }
          if (fac4 != -1)
            {
              f(i) = fac4;
              i++;
            }

          rang1 = f(0);
          rang2 = f(1);

          if (rang1 >= sz_faces_Cl)
            rang1 -= decal_virt;

          if (rang2 >= sz_faces_Cl)
            rang2 -= decal_virt;

          if (les_faces_Cl[rang1] == 3) // periodic face
            {
              if (les_faces_Cl[rang2] == 3)  // periodic face
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PERIO_PERIO;
              else if (les_faces_Cl[rang2] == 0) // wall face
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PERIO_PAROI;
              else if (les_faces_Cl[rang2] == 1) // free outlet or fluid inlet
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PERIO_FLUIDE;
              else
                {
                  Cerr << "Processing an edge separating two faces in a corner: " << finl;
                  Cerr << "one of these faces has a periodicity boundary condition " << finl;
                  Cerr << "and the other has a condition other than periodicity or wall" << finl;
                  Cerr << "This modification has not yet been implemented!!!" << finl;
                  // exit();
                }
            }
          else if (les_faces_Cl[rang1] == 0)  // wall face
            {
              if (les_faces_Cl[rang2] == 3)   // periodic face
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PERIO_PAROI;
              else if (les_faces_Cl[rang2] == 0) // wall face
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PAROI_PAROI;
              else if (les_faces_Cl[rang2] == 1)   // fluid face
                {
                  if (Option_VDF::traitement_coins)
                    type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PAROI_FLUIDE;
                  else
                    /* {
                     Cerr << "On traite une arete qui separe deux faces dans un coin : " << finl;
                     Cerr << "l'une de ces faces porte une condition limite de type paroi " << finl;
                     Cerr << "et l'autre porte une condition autre que periodicite" << finl;
                     Cerr << "On n a pas encore fait la modif!!!" << finl;
                     } */
                    // Do not change anything. This edge type does not require special treatment.
                    type_arete_coin_[num_arete_] = TypeAreteCoinVDF::VIDE;
                }
            }
          else if (les_faces_Cl[rang1] == 1)  // free outlet or fluid inlet face
            {
              if (les_faces_Cl[rang2] == 2)  // symmetry face
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::FLUIDE_NAVIER;
              else if (les_faces_Cl[rang2] == 0) // wall face
                {
                  if (Option_VDF::traitement_coins)
                    type_arete_coin_[num_arete_] = TypeAreteCoinVDF::FLUIDE_PAROI;
                  else
                    // Do not change anything. This edge type does not require special treatment.
                    type_arete_coin_[num_arete_] = TypeAreteCoinVDF::VIDE;
                }
              else if (les_faces_Cl[rang2] == 1) // free outlet or fluid inlet face
                {
                  if (Option_VDF::traitement_coins)
                    type_arete_coin_[num_arete_] = TypeAreteCoinVDF::FLUIDE_FLUIDE;
                  else
                    // Do not change anything. This edge type does not require special treatment.
                    type_arete_coin_[num_arete_] = TypeAreteCoinVDF::VIDE;
                }
              // Modif AC : 27/02/03
              else if (les_faces_Cl[rang2] == 3) // periodic face
                {
                  type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PERIO_FLUIDE;
                }
              // End Modif AC : 27/02/03
              else
                {
                  Cerr << "This case is not handled for les_faces_Cl[rang1] = 1 " << finl;
                  Cerr << "les_faces_Cl[rang2] = " << les_faces_Cl[rang2] << finl;
                  exit();
                }
            }
          else if (les_faces_Cl[rang1] == 2)   // symmetry face
            {
              if (les_faces_Cl[rang2] == 0)  // wall face
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::PAROI_NAVIER;
              else if (les_faces_Cl[rang2] == 1)  // free outlet or fluid inlet face
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::FLUIDE_NAVIER;
              else if (les_faces_Cl[rang2] == 2) // symmetry
                type_arete_coin_[num_arete_] = TypeAreteCoinVDF::NAVIER_NAVIER;
              else if (les_faces_Cl[rang2] == 3) // periodicity
                {
                  // Do not change anything. This edge type does not require special treatment.
                  type_arete_coin_[num_arete_] = TypeAreteCoinVDF::VIDE;
                }
              else
                {
                  Cerr << "This case is not handled for les_faces_Cl[rang1] = 2 " << finl;
                  Cerr << "les_faces_Cl[rang2] = " << les_faces_Cl[rang2] << finl;
                  Process::exit();
                }
            }
        }
    }

  else
    {
      Cerr << "Domaine_Cl_VDF::completer() expects a Domaine_VDF argument\n";
      Process::exit();
    }
}

/*! @brief Enforces the boundary conditions at time "temps" on the Champ_Inc.
 *
 * @param ch The field on which boundary conditions are imposed.
 * @param temps The current time.
 */
void Domaine_Cl_VDF::imposer_cond_lim(Champ_Inc_base& ch, double temps)
{
  static int init=0;
  DoubleTab& ch_tab = ch.valeurs(temps);
  const int N = ch_tab.line_size();
  if (sub_type(Champ_P0_VDF,ch)) { /* Do nothing */}
  else if(ch.nature_du_champ()==scalaire) { /* Do nothing */}
  else if (sub_type(Champ_Face_VDF,ch))
    {
      Champ_Face_VDF& ch_face = ref_cast(Champ_Face_VDF, ch);
      const Domaine_VDF& mon_dom_VDF = ch_face.domaine_vdf();
      int ndeb,nfin, num_face;

      for(int i=0; i<nb_cond_lim(); i++)
        {
          const Cond_lim_base& la_cl = les_conditions_limites(i).valeur();
          if (sub_type(Periodique,la_cl))
            {
              if (init == 0)
                {
                  // Ensure that the field has the same value
                  // on two periodic faces that are facing each other
                  const Periodique& la_cl_perio = ref_cast(Periodique,la_cl);
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
                  ndeb = le_bord.num_premiere_face();
                  nfin = ndeb + le_bord.nb_faces();
                  int voisine;
                  double moy;
                  for (num_face=ndeb; num_face<nfin; num_face++)
                    {
                      voisine = la_cl_perio.face_associee(num_face-ndeb) + ndeb;
                      if ( ch_tab[num_face] != ch_tab[voisine] )
                        {
                          //                           Cerr << "dans Domaine_Cl_VDF::imposer_cond_lim : on reajuste les vitesses!! pour la face num=" << num_face << finl;
                          //                           Cerr << "difference = ch_tab[num_face]-ch_tab[voisine]=" << ch_tab[num_face]-ch_tab[voisine] << finl;
                          moy = 0.5*(ch_tab[num_face] + ch_tab[voisine]);
                          ch_tab[num_face] = moy;
                          ch_tab[voisine] = moy;
                        }
                    }
                  // This should not be done on the first BC but once all BCs have been processed once, for multi-perio case with non-periodic IC
                  // init = 1;
                }
            }
          else if( sub_type(Navier,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();
              for (num_face=ndeb; num_face<nfin; num_face++)
                for (int n = 0; n < N; n++)
                  ch_tab(num_face, n) = 0;
            }
          else if ( sub_type(Dirichlet_entree_fluide,la_cl) )
            {
              const Dirichlet_entree_fluide& la_cl_diri = ref_cast(Dirichlet_entree_fluide,la_cl);
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();

              for (num_face = ndeb; num_face < nfin; num_face++)
                for (int n = 0; n < N; n++)
                  {
                    // WEC : optimizable (for each face searches for the right time!)
                    ch_tab(num_face, n) = la_cl_diri.val_imp_au_temps(temps, num_face - ndeb, N * mon_dom_VDF.orientation(num_face) + n);
                  }
            }
          else if ( sub_type(Dirichlet_paroi_fixe,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();
              for (num_face=ndeb; num_face<nfin; num_face++)
                for (int n = 0; n < N; n++)
                  ch_tab(num_face, n) = 0;
            }
          else if ( sub_type(Dirichlet_paroi_defilante,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();
              for (num_face=ndeb; num_face<nfin; num_face++)
                for (int n = 0; n < N; n++)
                  ch_tab(num_face, n) = 0;
            }
        }
      init = 1;
    }
  else
    {
      Cerr << "The type OWN_PTR(Champ_Inc_base) " <<  ch.que_suis_je() << " is not supported in VDF\n";
      exit();
    }
  ch_tab.echange_espace_virtuel();
  Debog::verifier("Domaine_Cl_VDF::imposer_cond_lim ch_tab",ch_tab);
}


int Domaine_Cl_VDF::nb_faces_sortie_libre() const
{
  int compteur = 0;
  for (const auto &itr : les_conditions_limites_)
    {
      if (sub_type(Sortie_libre_pression_imposee, itr.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF, itr->frontiere_dis());
          compteur += le_bord.nb_faces();
        }
    }
  return compteur;
}

Domaine_VDF& Domaine_Cl_VDF::domaine_VDF()
{
  return ref_cast(Domaine_VDF, domaine_dis());
}

const Domaine_VDF& Domaine_Cl_VDF::domaine_VDF() const
{
  return ref_cast(Domaine_VDF, domaine_dis());
}

int Domaine_Cl_VDF::nb_faces_bord() const
{
  return domaine_VDF().nb_faces_bord();
}
