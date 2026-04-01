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
#include <Pb_Couple_Rayonnement.h>
#include <Paroi_Rayo_transp.h>
#include <Probleme_base.h>
#include <Fluide_base.h>

Implemente_instanciable(Pb_Couple_Rayonnement, "Pb_Couple_Rayonnement", Probleme_Couple);

Entree& Pb_Couple_Rayonnement::readOn(Entree& is) { return is; }

Sortie& Pb_Couple_Rayonnement::printOn(Sortie& os) const { return Probleme_Couple::printOn(os); }

void Pb_Couple_Rayonnement::initialize()
{
  completer();
  Probleme_Couple::initialize();
  le_modele_rayo().preparer_calcul();
}

int Pb_Couple_Rayonnement::associer_(Objet_U& ob)
{
  int set_type_rayo = 0;
  if (Probleme_Couple::associer_(ob))
    {
      set_type_rayo=1;

    }
  else
    {
      if( sub_type(Modele_Rayonnement_Milieu_Transparent, ob))
        {
          set_type_rayo=1;
          Cerr << "association du modele au pbc" << finl;
          le_modele_rayo_associe(ref_cast(Modele_Rayonnement_Milieu_Transparent, ob));
        }
      else
        return 0;
    }
  if (set_type_rayo==1)
    return 1;
  else
    return 0;
}

void Pb_Couple_Rayonnement::le_modele_rayo_associe(const Modele_Rayonnement_Milieu_Transparent& un_modele_de_rayonnement)
{
  le_modele_de_rayo_ = un_modele_de_rayonnement;
}

void Pb_Couple_Rayonnement::associer_cl_base(const Cond_lim_base& les_cl)
{
  les_cl_ = les_cl;
}

int Pb_Couple_Rayonnement::postraiter(int force)
{
  int ok = Probleme_Couple::postraiter(force);
  if (!ok)
    return 0;

  // Impression en plus du modele de rayonnement
  const Modele_Rayonnement_Milieu_Transparent& mod_rayo = le_modele_de_rayo_.valeur();
  if (mod_rayo.processeur_rayonnant() != -1)
    if (schema_temps().limpr())
      {
        Cout << "Impression des flux radiatifs sur les bords de rayonnement" << finl;
        Cout << "----------------------------------------------------------------" << finl;
        le_modele_de_rayo_->imprimer_flux_radiatifs(Cout);
      }
  return 1;
}

void Pb_Couple_Rayonnement::validateTimeStep()
{
  Probleme_Couple::validateTimeStep();
  le_modele_rayo().mettre_a_jour(presentTime());
}

void Pb_Couple_Rayonnement::completer()
{
  le_modele_de_rayo_->discretiser(ref_cast(Probleme_base,probleme(0)).discretisation(), ref_cast(Probleme_base,probleme(0)).domaine());
  Modele_Rayonnement_Milieu_Transparent& mod_rayo = le_modele_de_rayo_.valeur();

  int nb_pb_ray = 0;
  int compte_nb_bords_rayo = 0;

  for (int l = 0; l < nb_problemes(); l++)
    {
      Probleme_base& le_pb = ref_cast(Probleme_base, probleme(l));
      if (sub_type(Fluide_base, le_pb.milieu()))
        {
          if (ref_cast(Fluide_base, le_pb.milieu()).is_rayo_transp())
            nb_pb_ray++;
        }
    }
  if (nb_pb_ray > 1)
    {
      Cerr << "Pb_Couple_Rayonnement::completer - We can only treat 1 transparent medium at present. You defined " << nb_pb_ray << " !!!" << finl;
      Process::exit();
    }
  else if (nb_pb_ray == 0)
    Process::exit("Pb_Couple_Rayonnement::completer - You should define the transparent medium using the flag transparent_medium_radiation !!!\n");

  for (int l = 0; l < nb_problemes(); l++)
    {
      Probleme_base& le_pb = ref_cast(Probleme_base, probleme(l));

      bool is_pb_fluide = false;

      if (sub_type(Fluide_base, le_pb.milieu()))
        if (ref_cast(Fluide_base, le_pb.milieu()).is_rayo_transp())
          is_pb_fluide = true;

      if (is_pb_fluide)
        {
          le_modele_de_rayo_->set_nom_pb_rayonnant(le_pb.le_nom());
          Cerr << "Le probleme rayonnant trouve est : " << le_pb.le_nom() << finl;
        }

      for (int j = 0; j < le_pb.nombre_d_equations(); j++)
        {
          Domaine_Cl_dis_base& la_zcl = le_pb.equation(j).domaine_Cl_dis();
          for (int num_cl = 0; num_cl < la_zcl.nb_cond_lim(); num_cl++)
            {
              Cond_lim_base& la_cl = la_zcl.les_conditions_limites(num_cl).valeur();

              Cond_lim_rayo_milieu_transp *la_cl_rayo;
              if (la_cl.is_bc_rayo_milieu_transp(la_cl_rayo))
                {
                  ((*la_cl_rayo)).associer_modele_rayo(mod_rayo);

                  // on associe la cl liee au pb fluide
                  if (is_pb_fluide)
                    {
                      int ok = 0;
                      for (int i = 0; i < mod_rayo.nb_faces_totales(); i++)
                        {
                          if (mod_rayo.face_rayonnante(i).nom_bord_rayo() == la_zcl.les_conditions_limites(num_cl)->frontiere_dis().le_nom())
                            //if (la_cl.frontiere_dis().frontiere().nb_faces()!=0)
                            {
                              if (mod_rayo.face_rayonnante(i).emissivite() != -1)
                                ok = 1;
                              //Cerr<< mod_rayo.face_rayonnante(i).nom_bord_rayo()<<" associe a "<<la_zcl.les_conditions_limites(num_cl).frontiere_dis().le_nom()<<finl;
                              mod_rayo.face_rayonnante(i).ensembles_faces_bord(0).associer_les_cl(la_cl);
                              compte_nb_bords_rayo += 1;
                            }
                        }
                      if (ok == 0)
                        {
                          Cerr << "La condition limite de nom " << la_zcl.les_conditions_limites(num_cl)->frontiere_dis().le_nom()
                               << " est definie comme rayonnante, mais n'est pas dans la liste des faces rayonnantes ou son emissivite vaut -1" << finl;
                          Process::exit();
                        }
                    }
                }
            }
        }
      if (is_pb_fluide)
        {
          for (int i = 0; i < mod_rayo.nb_faces_totales(); i++)
            {
              if (!mod_rayo.face_rayonnante(i).ensembles_faces_bord(0).is_ok() && (mod_rayo.face_rayonnante(i).emissivite() != -1))
                {
                  Cerr << "Le bord " << mod_rayo.face_rayonnante(i).nom_bord_rayo_lu() << " n'a pas ete asssocie a une condition limite rayonnante." << finl;
                  Cerr << "Soit vous mettez une condition limite rayonnante pour " << mod_rayo.face_rayonnante(i).nom_bord_rayo() << finl;
                  Cerr << "Soit vous affectez une emissivite de -1 a ce bord." << finl;
                  Cerr << finl;
                }
            }
        }
    }

  if (compte_nb_bords_rayo != mod_rayo.nb_faces_rayonnantes())
    abort();

  if (nproc() == 1)
    mod_rayo.associer_processeur_rayonnant(me());
  else
    {
      if (compte_nb_bords_rayo != 0)
        {
          LIST(Nom) collectnoms;
          for (int i = 0; i < mod_rayo.nb_faces_rayonnantes(); i++)
            {
              if (mod_rayo.face_rayonnante(i).ensembles_faces_bord(0).nb_faces_bord() != 0)
                collectnoms.add(mod_rayo.face_rayonnante(i).nom_bord_rayo_lu());
            }
          Cerr << me() << collectnoms << finl;
          // on verifie que l'on a bien tous les noms
          if (me() == 0)
            mod_rayo.associer_processeur_rayonnant(me());
          else
            mod_rayo.associer_processeur_rayonnant(-1);
        }
      else
        {
          //tout est ok
          Cerr << "On redimenssionne le tableau de faces de bord" << finl;
          Cerr << "compte_nb_bords_rayo = " << compte_nb_bords_rayo << finl;
          Cerr << "mod_rayo.nb_faces_rayonnantes() = " << mod_rayo.nb_faces_rayonnantes() << finl;
          mod_rayo.associer_processeur_rayonnant(-1);
        }
    }
}
