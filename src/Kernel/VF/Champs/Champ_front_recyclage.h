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

#ifndef Champ_front_recyclage_included
#define Champ_front_recyclage_included

#include <Ch_front_var_instationnaire_dep.h>
#include <TRUSTArrays.h>
#include <TRUSTTabs.h>
#include <TRUST_Ref.h>
#include <Parser_U.h>
#include <Noms.h>

class Domaine_Cl_dis_base;
class Domaine_dis_base;
class Champ_Inc_base;
class Front_dis_base;
class Equation_base;
class Milieu_base;
class Param;

/*! @brief Champ_front_recyclage
 *
 *                         delt_dist                                             delt_dist
 *         pb1           <---------->             pb2                            <-------->    pb
 *                                  ch_fr2                                     ch_fr
 *      _____________________         ____________________________               ____________________________
 *     |                 |   |       |                            |              |        |                  |
 *     |         dom1    |   |       |            dom2            |              |        |     dom          |
 *     |         ch1     |   |       |bord2                       |              |bord    |     ch           |
 *     |_________________|___|       |____________________________|              |________|__________________|
 *                     plan1                                                             plan
 *
 *                              Fig. 1                                                        Fig. 2
 *
 *      The goal of this class is to evaluate the values of a boundary field (ch_fr2) on the boundary of a domain (bord2)
 *      by exploiting the values of a field 1 (ch1, called the evaluator field) evaluated in a plane (plan1) at a distance
 *      delt_dist from bord2 (Fig. 1).
 *
 */

//     Problems pb2 and pb1, domains dom2 and dom1, and fields 2 and 1 can be identical
//     (pb2=pb1=pb dom2=dom1=dom and ch2=ch1=ch)
//     in which case the values of the boundary field (ch_fr) on the boundary (bord) will be built
//     from the values of field ch (which becomes the evaluator field) computed in the plane (plan)
//     at a distance delt_dist from the boundary (Fig. 2).

//     The expression for the values assigned to the boundary field on bord2 (or bord) is:
//     val_ch_fr2(dir) = ampli_moy_imposee(dir)*moyenne_imposee(dir)
//                       + ampli_fluct(dir)*(val_evaluateur(dir)-ampli_moy_recyclee(dir)*moyenne_recyclee(dir))
//
//     val_ch_fr2         :  values taken by the boundary field ch_fr2 (or ch_fr)
//     moyenne_imposee    :  mean of the boundary field (can be imposed analytically or read from a file)
//     val_evaluateur     :  values of the evaluator field in plane1 (or plane) evaluated by interpolation
//     moyenne_recyclee   :  mean of the evaluator field (can be computed by surface averaging or from a special treatment)
//     ampli_moy_imposee  :  amplification factor for the imposed mean
//     ampli_moy_recyclee :  amplification factor for the recycled mean
//     ampli_fluct        :  amplification factor for the recycled fluctuation
//     dir                :  direction

//     User syntax:
//     Champ_front_recyclage
//     {
//      pb_champ_evaluateur nom_pb1 nom_inco1 nb_compo1
//      [ moyenne_imposee methode_moy [fichier] nom_fich1 (nom_fich2) ]
//      [ moyenne_recyclee methode_recyc [fichier] nom_fich1 (nom_fich2) ]
//      [ direction_anisotrope direction ]
//      [ distance_plan dist0 dist1 (dist2) ]
//      [ ampli_fluctuation nb_comp ampli_fluc0 ampli_fluc1 (ampli_fluc2) ]
//      [ ampli_moyenne_imposee nb_comp ampli_moy0 ampli_moy1 (ampli_moy2) ]
//      [ ampli_moyenne_recyclee nb_comp ampli_recy0 ampli_recy1 (ampli_recy2) ]
//     }
//
//     methode_moy = 1 (keyword profil)
//                     to impose an analytical profile
//     methode_moy = 2 (keyword interpolation):
//                       reads from a file and builds a mean field
//                     by interpolating the data read.
//                     The mean is built for a preferred direction
//                     (direction_anisotrope) and is 0 for other directions
//     methode_moy = 3 (keyword connexion_approchee)
//                     reads from a file and retains the value of the
//                     variable read by connection with the closest point
//                     to the considered boundary face
//     methode_moy = 4 (keyword connexion_exacte)
//                     reads a geometry file containing the coordinates of points
//                     located in the evaluation plane and reads from a separate file
//                     the mean values. The mean values read
//                     are stored when the exact correspondence between
//                     facing points is verified.
//     methode_moy = 5 (keyword logarithmique)
//                     builds the mean using a wall law (logarithmic)
//
//
//    methode_moy = 2 and methode_moy = 3:
//                     a single file to read containing positions and values of the variable
//    methode_moy = 4: two files to read: the first containing the values of the variable
//                                         and the second containing the positions
//
//
//    methode_recyc = 1 (keyword surfacique)
//                      surface average of the recycled values
//                     (the average is computed on bord2 where the values are retrieved)
//    methode_recyc = 2 (keyword interpolation):
//                       see methode_moy = 2
//    methode_recyc = 3 (keyword connexion_approchee)
//                     see methode_moy = 3
//    methode_recyc = 4 (keyword connexion_exacte)
//                     see methode_moy = 4
//
//    methode_recyc = 2 and methode_moy = 3:
//                     a single file to read containing positions and values of the variable
//    methode_recyc = 4: two files to read: the first containing the values of the variable
//                                           and the second containing the positions
//
//////////////////////////////////////////////////////////////////////////////

class Champ_front_recyclage : public Ch_front_var_instationnaire_dep
{

  Declare_instanciable_sans_constructeur(Champ_front_recyclage);

public:

  Champ_front_recyclage();
  int lire_info_moyenne_imposee(Entree& is);
  int lire_info_moyenne_recyclee(Entree& is);
  void calcul_moyenne_imposee(const DoubleTab& tab,double temps);
  void calcul_moyenne_recyclee(const DoubleTab& tab,double temps);
  void initialiser_moyenne_imposee(DoubleTab& moyenne);
  void initialiser_moyenne_recyclee(DoubleTab& moyenne);
  void associer_champ_evaluateur(const Nom&, const Motcle&);
  int initialiser(double temps, const Champ_Inc_base& inco) override;
  void mettre_a_jour(double temps) override;

  static void get_coord_faces(const Frontiere_dis_base& fr_vf,
                              DoubleTab& coords,
                              const DoubleVect& delt_dist);

  void lire_fichier_format1(DoubleTab& moyenne,
                            const Frontiere_dis_base& fr_vf,
                            const Nom& nom_fich);
  void lire_fichier_format2(DoubleTab& moyenne,
                            const Frontiere_dis_base& fr_vf,
                            const Nom& nom_fich);
  void lire_fichier_format3(DoubleTab& moyenne,
                            const Frontiere_dis_base& fr_vf,
                            const Nom& nom_fich1, const Nom& nom_fich2);
  double UPb(double y,Nom nom_fich);

protected :
  void set_param(Param& param) const override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;

  OBS_PTR(Champ_Inc_base) l_inconnue1;  //Reference to the unknown field (ch1) used as evaluator
  //in the plane where the values are retrieved

  DoubleVect delt_dist;             //distance vector between bord2 and the plane (plan1)
  //for computing the values of ch1

  Nom nom_pb1;                      //name of the evaluator problem (pb1)

  Motcle nom_inco1;                 //name of the evaluator unknown field (ch1)

  DoubleTab moyenne_imposee_;            //See description above
  DoubleTab moyenne_recyclee_;
  DoubleVect ampli_fluct_;
  DoubleVect ampli_moy_imposee_;
  DoubleVect ampli_moy_recyclee_;

  int methode_moy_impos_;            //method to evaluate moyenne_imposee_
  int methode_moy_recycl_;       //method to evaluate moyenne_recyclee_

  Nom fich_impos_,fich_recycl_;     //Names of files optionally used
  Nom fich_maillage_;                    //to evaluate moyenne_imposee_ and moyenne_recyclee_

  int ndir;                            //direction of anisotropy

  VECT(Parser_U) profil_2;             //Parser and expressions for imposing an analytical mean
  Noms fcts_xyz;

  double u_tau,diametre,visco_cin;  //parameters for imposing a mean using a logarithmic profile

  // Set of points where inconnue1 must be evaluated (only coordinates
  //  included in the local domain1), sorted by the processor to which
  //  the evaluation result must be sent. Note that there is not necessarily
  //  equality between the local unknown on the face and the remote unknown on the face if
  //  the surface meshes of the local and remote boundaries are not identical...
  DoubleTabs inconnues1_coords_to_eval_;

  // For each point where inconnue1 must be evaluated, index of the element in which
  //  this point lies (always by destination processor)
  ArrsOfInt inconnues1_elems_;

  // Upon reception of values, indices of the boundary faces where the
  //  result received from each processor must be stored
  ArrsOfInt inconnues2_faces_;
};

#endif
